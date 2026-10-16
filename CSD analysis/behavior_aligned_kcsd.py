"""
Paired onset/offset 1D kCSD for behavior and play-bout LFP exports.

Baseline correction belongs to this loader, not to the numerical kCSD
estimator.  For every event and channel, one baseline is measured from the
pre-onset LFP and subtracted from both the onset-aligned and offset-aligned
epochs.  The two corrected evoked voltages are then passed independently to
the general 1D kCSD implementation in ``csd_1d_pipeline.py``.

Each alignment leaf must contain:
  - lfp.npy        (n_events, n_channels, n_time)
  - x_coords.npy   (n_channels,)
  - y_coords.npy   (n_channels,)
  - T0.npy         scalar alignment sample
  - fs.npy         optional scalar sampling rate (default: 2500 Hz)
  - event_ids.npy  strongly recommended; must match between onset and offset

Alignment folders may occur anywhere below a root, for example:

    behavior/Pounce/onset/session/probe2/
    behavior/Pounce/offset/session/probe2/
    play_bout/onset/session/probe2/
    play_bout/offset/session/probe2/

Discovery replaces the ``onset`` path component with ``offset``, so behavior
and play-bout datasets use exactly the same processing code.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from csd_1d_pipeline import (
    DELTA_CSD_DEFAULTS,
    DatasetConfig,
    apply_normalize,
    column_y,
    compute_fd_csd_columns,
    compute_kcsd_1d_columns,
    load_raw,
    plot_column,
    prepare_channels,
    probe_geometry_kind,
    save_fd_outputs,
    save_kcsd_outputs,
)


PathLike = Union[str, os.PathLike]

BEHAVIOR_ALIGNED_CSD_DEFAULTS: Dict[str, Any] = {
    **DELTA_CSD_DEFAULTS,
    "t_pre": 500.0,
    "t_post": 500.0,
    "baseline_correct": True,
    # This window is relative to behavior/play-bout onset, not offset.
    "onset_baseline_window_ms": (-500.0, 0.0),
    "output_subdir": "csd_1d_onset_baseline",
}


def _load_event_ids(folder: PathLike) -> Optional[np.ndarray]:
    path = os.path.join(os.fspath(folder), "event_ids.npy")
    if not os.path.isfile(path):
        return None
    return np.asarray(np.load(path, allow_pickle=False)).reshape(-1)


def load_paired_raw(
    onset_folder: PathLike,
    offset_folder: PathLike,
) -> Dict[str, Any]:
    """Load and validate paired raw onset/offset LFP stacks."""
    onset, x_on, y_on, t0_on, fs_on = load_raw(onset_folder)
    offset, x_off, y_off, t0_off, fs_off = load_raw(offset_folder)

    if onset.shape[:2] != offset.shape[:2]:
        raise ValueError(
            "Onset/offset event or channel mismatch: "
            f"onset={onset.shape}, offset={offset.shape}"
        )
    if not np.array_equal(x_on, x_off, equal_nan=True) or not np.array_equal(
        y_on, y_off, equal_nan=True
    ):
        raise ValueError("Onset and offset channel coordinates do not match")
    if fs_on is not None and fs_off is not None and not np.isclose(fs_on, fs_off):
        raise ValueError(f"Onset/offset sampling rates differ: {fs_on} vs {fs_off}")

    onset_ids = _load_event_ids(onset_folder)
    offset_ids = _load_event_ids(offset_folder)
    if (onset_ids is None) != (offset_ids is None):
        raise ValueError(
            "event_ids.npy must be present in both onset and offset folders, or neither"
        )
    if onset_ids is not None:
        if len(onset_ids) != onset.shape[0]:
            raise ValueError(
                f"event_ids length {len(onset_ids)} != n_events {onset.shape[0]}"
            )
        if not np.array_equal(onset_ids, offset_ids):
            raise ValueError(
                "Onset and offset event_ids.npy differ; refusing to apply mismatched baselines"
            )
        pairing = "event_ids"
    else:
        pairing = "row_order"

    return {
        "onset": onset,
        "offset": offset,
        "x": x_on,
        "y": y_on,
        "t0_onset": t0_on,
        "t0_offset": t0_off,
        "fs_loaded": fs_on if fs_on is not None else fs_off,
        "event_ids": onset_ids,
        "pairing": pairing,
    }


def _window_bounds(
    t0: int,
    fs: float,
    start_ms: float,
    end_ms: float,
    n_time: int,
    label: str,
) -> Tuple[int, int]:
    start = t0 + int(round(start_ms * fs / 1000.0))
    end = t0 + int(round(end_ms * fs / 1000.0))
    if start < 0 or end > n_time or start >= end:
        raise RuntimeError(
            f"{label} window [{start}:{end}) exceeds signal length {n_time} "
            f"(T0={t0}, relative window={start_ms:g}…{end_ms:g} ms)"
        )
    return start, end


def epoch_paired_with_onset_baseline(
    onset_lfp: np.ndarray,
    offset_lfp: np.ndarray,
    t0_onset: int,
    t0_offset: int,
    fs: float,
    t_pre: float,
    t_post: float,
    step_ms: float,
    *,
    onset_baseline_window_ms: Tuple[float, float] = (-500.0, 0.0),
    baseline_correct: bool = True,
    normalize: str = "none",
) -> Dict[str, np.ndarray]:
    """
    Epoch paired alignments and apply the onset baseline to both.

    Baseline is event- and channel-specific with shape
    ``(n_events, n_channels, 1)``.  It is never estimated from the offset
    epoch.  Rows must represent the same events in the same order.
    """
    if onset_lfp.ndim != 3 or offset_lfp.ndim != 3:
        raise ValueError("onset_lfp and offset_lfp must be 3-D event × channel × time")
    if onset_lfp.shape[:2] != offset_lfp.shape[:2]:
        raise ValueError(
            f"Onset/offset event or channel mismatch: {onset_lfp.shape} vs {offset_lfp.shape}"
        )
    if onset_lfp.shape[0] == 0:
        raise ValueError("Cannot estimate paired kCSD with zero events")

    onset_start, onset_end = _window_bounds(
        t0_onset,
        fs,
        -float(t_pre),
        float(t_post),
        onset_lfp.shape[2],
        "onset epoch",
    )
    offset_start, offset_end = _window_bounds(
        t0_offset,
        fs,
        -float(t_pre),
        float(t_post),
        offset_lfp.shape[2],
        "offset epoch",
    )
    baseline_start, baseline_end = _window_bounds(
        t0_onset,
        fs,
        float(onset_baseline_window_ms[0]),
        float(onset_baseline_window_ms[1]),
        onset_lfp.shape[2],
        "onset baseline",
    )

    onset_epochs = onset_lfp[:, :, onset_start:onset_end].astype(float, copy=True)
    offset_epochs = offset_lfp[:, :, offset_start:offset_end].astype(float, copy=True)
    if onset_epochs.shape[2] != offset_epochs.shape[2]:
        raise ValueError(
            f"Onset/offset epoch lengths differ: {onset_epochs.shape[2]} vs "
            f"{offset_epochs.shape[2]}"
        )

    onset_baseline = onset_lfp[:, :, baseline_start:baseline_end].mean(
        axis=2, keepdims=True
    )
    if baseline_correct:
        onset_epochs -= onset_baseline
        offset_epochs -= onset_baseline

    onset_evoked = apply_normalize(onset_epochs.mean(axis=0), normalize)
    offset_evoked = apply_normalize(offset_epochs.mean(axis=0), normalize)
    pre_samples = int(round(float(t_pre) * fs / 1000.0))
    time_ms = (np.arange(onset_epochs.shape[2]) - pre_samples) / fs * 1000.0
    step = max(1, int(round(step_ms * fs / 1000.0)))
    time_indices = np.arange(0, onset_epochs.shape[2], step)

    return {
        "onset_evoked": onset_evoked,
        "offset_evoked": offset_evoked,
        "onset_epochs": onset_epochs,
        "offset_epochs": offset_epochs,
        "onset_baseline": onset_baseline,
        "time_ms": time_ms,
        "time_indices": time_indices,
    }


def _save_pair_metadata(
    folder: PathLike,
    output_subdir: str,
    metadata: Dict[str, Any],
    onset_baseline: np.ndarray,
) -> None:
    output = os.path.join(os.fspath(folder), output_subdir)
    os.makedirs(output, exist_ok=True)
    np.save(os.path.join(output, "onset_baseline.npy"), onset_baseline)
    with open(
        os.path.join(output, "paired_baseline_metadata.json"),
        "w",
        encoding="utf-8",
    ) as stream:
        json.dump(metadata, stream, indent=2, sort_keys=True)


def _plot_alignment(
    label: str,
    name: str,
    folder: PathLike,
    output_subdir: str,
    csd: Dict[float, np.ndarray],
    x: np.ndarray,
    y: np.ndarray,
    time_ms: np.ndarray,
    time_indices: np.ndarray,
    *,
    sort_by_depth: bool,
    show: bool,
    save_plots: bool,
) -> None:
    figure_dir = os.path.join(os.fspath(folder), output_subdir, "figures")
    if save_plots:
        os.makedirs(figure_dir, exist_ok=True)
    first_x = float(sorted(csd)[0])
    for ux in sorted(csd):
        ux = float(ux)
        save_path = None
        if save_plots:
            save_path = os.path.join(figure_dir, f"kcsd_x{ux}.png")
        plot_column(
            csd[ux],
            time_ms,
            time_indices,
            column_y(x, y, ux, sort_by_depth=sort_by_depth),
            title=f"{name} — {label} (onset baseline) — 1D kCSD x={ux}",
            cbar_label="kCSD (a.u.)",
            save_path=save_path,
            show=show and ux == first_x,
        )


def run_behavior_aligned_kcsd(
    onset_folder: PathLike,
    offset_folder: PathLike,
    *,
    name: Optional[str] = None,
    stop_on_missing_event_ids: bool = True,
    **overrides: Any,
) -> Dict[str, Any]:
    """
    Run paired onset/offset kCSD using one shared pre-onset baseline.

    ``overrides`` accepts the general ``DatasetConfig`` fields plus
    ``onset_baseline_window_ms``.  The default is -500…0 ms relative to onset.
    kCSD itself remains unchanged and receives the two prepared evoked LFPs.
    """
    valid_overrides = set(DatasetConfig.__dataclass_fields__) | {
        "onset_baseline_window_ms"
    }
    unknown = sorted(set(overrides) - valid_overrides)
    if unknown:
        raise TypeError(f"Unknown paired kCSD option(s): {', '.join(unknown)}")

    settings = dict(BEHAVIOR_ALIGNED_CSD_DEFAULTS)
    settings.update(overrides)
    baseline_window = tuple(settings.pop("onset_baseline_window_ms"))
    if (
        len(baseline_window) != 2
        or baseline_window[0] >= baseline_window[1]
        or baseline_window[1] > 0
    ):
        raise ValueError(
            "onset_baseline_window_ms must be (start, end) with start < end <= 0"
        )

    paired = load_paired_raw(onset_folder, offset_folder)
    if stop_on_missing_event_ids and paired["pairing"] != "event_ids":
        raise ValueError(
            "event_ids.npy is required when stop_on_missing_event_ids=True"
        )

    cfg_values = {
        key: value
        for key, value in settings.items()
        if key in DatasetConfig.__dataclass_fields__ and key not in {"folder", "name"}
    }
    cfg = DatasetConfig(folder=os.fspath(onset_folder), name=name, **cfg_values)
    if paired["fs_loaded"] is not None:
        cfg.fs = float(paired["fs_loaded"])

    n_events = paired["onset"].shape[0]
    combined = np.concatenate((paired["onset"], paired["offset"]), axis=0)
    combined, x, y, channel_meta = prepare_channels(
        combined, paired["x"], paired["y"], cfg
    )
    onset = combined[:n_events]
    offset = combined[n_events:]

    prepared = epoch_paired_with_onset_baseline(
        onset,
        offset,
        paired["t0_onset"],
        paired["t0_offset"],
        cfg.fs,
        cfg.t_pre,
        cfg.t_post,
        cfg.step_ms,
        onset_baseline_window_ms=baseline_window,
        baseline_correct=cfg.baseline_correct,
        normalize=cfg.normalize,
    )

    resolved_name = name or os.path.relpath(
        os.fspath(onset_folder), os.path.commonpath((onset_folder, offset_folder))
    )
    metadata: Dict[str, Any] = {
        "name": resolved_name,
        "onset_folder": os.path.abspath(os.fspath(onset_folder)),
        "offset_folder": os.path.abspath(os.fspath(offset_folder)),
        "baseline_source": "onset",
        "onset_baseline_window_ms": list(baseline_window),
        "baseline_correct": bool(cfg.baseline_correct),
        "pairing": paired["pairing"],
        "n_events": int(n_events),
        "geometry": probe_geometry_kind(x),
        "channel_processing": channel_meta,
        "config": asdict(cfg),
    }
    result: Dict[str, Any] = {
        "name": resolved_name,
        "meta": metadata,
        "event_ids": paired["event_ids"],
        "x": x,
        "y": y,
        **prepared,
        "onset_csd_kcsd": None,
        "offset_csd_kcsd": None,
        "onset_csd_fd": None,
        "offset_csd_fd": None,
    }

    if cfg.compute_kcsd:
        csd_kwargs = {
            "y_scale": cfg.y_scale,
            "h": cfg.h,
            "n_src_factor": cfg.n_src_factor,
            "smooth_sigma_time": cfg.smooth_sigma_time,
            "sort_by_depth": cfg.sort_by_depth,
        }
        onset_csd = compute_kcsd_1d_columns(
            prepared["onset_evoked"], x, y, prepared["time_indices"], **csd_kwargs
        )
        offset_csd = compute_kcsd_1d_columns(
            prepared["offset_evoked"], x, y, prepared["time_indices"], **csd_kwargs
        )
        result["onset_csd_kcsd"] = onset_csd
        result["offset_csd_kcsd"] = offset_csd

        if cfg.save:
            for folder, csd, evoked in (
                (onset_folder, onset_csd, prepared["onset_evoked"]),
                (offset_folder, offset_csd, prepared["offset_evoked"]),
            ):
                save_kcsd_outputs(
                    folder,
                    csd,
                    x,
                    y,
                    evoked,
                    prepared["time_ms"],
                    prepared["time_indices"],
                    output_subdir=cfg.output_subdir,
                    smooth_sigma_time=cfg.smooth_sigma_time,
                    sort_by_depth=cfg.sort_by_depth,
                )

        if cfg.plot or cfg.save_plots:
            _plot_alignment(
                "onset",
                resolved_name,
                onset_folder,
                cfg.output_subdir,
                onset_csd,
                x,
                y,
                prepared["time_ms"],
                prepared["time_indices"],
                sort_by_depth=cfg.sort_by_depth,
                show=cfg.plot,
                save_plots=cfg.save_plots,
            )
            _plot_alignment(
                "offset",
                resolved_name,
                offset_folder,
                cfg.output_subdir,
                offset_csd,
                x,
                y,
                prepared["time_ms"],
                prepared["time_indices"],
                sort_by_depth=cfg.sort_by_depth,
                show=cfg.plot,
                save_plots=cfg.save_plots,
            )

    if cfg.compute_fd:
        fd_kwargs = {
            "tissue_sigma": cfg.tissue_sigma,
            "lfp_in_uV": cfg.lfp_in_uV,
            "sigma_spatial_ch": cfg.sigma_spatial_ch,
            "sort_by_depth": cfg.sort_by_depth,
        }
        onset_fd = compute_fd_csd_columns(prepared["onset_evoked"], x, y, **fd_kwargs)
        offset_fd = compute_fd_csd_columns(prepared["offset_evoked"], x, y, **fd_kwargs)
        result["onset_csd_fd"] = onset_fd
        result["offset_csd_fd"] = offset_fd
        if cfg.save:
            save_fd_outputs(onset_folder, onset_fd, output_subdir=cfg.output_subdir)
            save_fd_outputs(offset_folder, offset_fd, output_subdir=cfg.output_subdir)

    if cfg.save:
        for folder in (onset_folder, offset_folder):
            _save_pair_metadata(
                folder,
                cfg.output_subdir,
                metadata,
                prepared["onset_baseline"],
            )

    return result


def discover_behavior_aligned_pairs(root: PathLike) -> List[Dict[str, str]]:
    """Find onset leaves and their mirrored offset leaves below ``root``."""
    root = os.path.abspath(os.fspath(root))
    required = {"lfp.npy", "x_coords.npy", "y_coords.npy", "T0.npy"}
    found: List[Dict[str, str]] = []
    for dirpath, _dirnames, filenames in os.walk(root):
        if not required.issubset(filenames):
            continue
        rel_parts = list(os.path.relpath(dirpath, root).split(os.sep))
        onset_positions = [i for i, part in enumerate(rel_parts) if part == "onset"]
        if not onset_positions:
            continue
        pair_parts = rel_parts.copy()
        pair_parts[onset_positions[-1]] = "offset"
        offset_folder = os.path.join(root, *pair_parts)
        if not os.path.isdir(offset_folder):
            continue
        if not required.issubset(set(os.listdir(offset_folder))):
            continue
        name_parts = rel_parts.copy()
        del name_parts[onset_positions[-1]]
        found.append(
            {
                "name": "/".join(name_parts),
                "onset_folder": dirpath,
                "offset_folder": offset_folder,
            }
        )
    found.sort(key=lambda item: item["name"])
    return found


def run_behavior_aligned_kcsd_from_root(
    root: PathLike,
    *,
    stop_on_error: bool = False,
    **overrides: Any,
) -> List[Dict[str, Any]]:
    """Discover and process all behavior and play-bout onset/offset pairs."""
    pairs = discover_behavior_aligned_pairs(root)
    if not pairs:
        raise FileNotFoundError(
            f"No paired onset/offset CSD-ready folders found under {root}"
        )

    results: List[Dict[str, Any]] = []
    for pair in pairs:
        try:
            results.append(
                run_behavior_aligned_kcsd(
                    pair["onset_folder"],
                    pair["offset_folder"],
                    name=pair["name"],
                    **overrides,
                )
            )
        except Exception as exc:  # noqa: BLE001 - batch processing should continue
            if stop_on_error:
                raise
            results.append({"name": pair["name"], "error": str(exc), **pair})
    return results


__all__: Sequence[str] = (
    "BEHAVIOR_ALIGNED_CSD_DEFAULTS",
    "discover_behavior_aligned_pairs",
    "epoch_paired_with_onset_baseline",
    "load_paired_raw",
    "run_behavior_aligned_kcsd",
    "run_behavior_aligned_kcsd_from_root",
)
