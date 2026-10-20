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

``estimate_animal_pooled_kcsd`` is separate from that per-session run. It
concatenates trials from sessions of the same animal for play bouts and the
play, aggression, and exploration categories, forms the trial-weighted PSTH,
and only then estimates kCSD.
"""

from __future__ import annotations

import csv
import json
import os
import re
import shutil
import time
from dataclasses import asdict
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.interpolate import interp1d

from area_csd_consistency import animal_key_from_session
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
    "t_pre": 1000.0,
    "t_post": 2000.0,
    "baseline_correct": True,
    # This window is relative to behavior/play-bout onset, not offset.
    "onset_baseline_window_ms": (-1000.0, 0.0),
    "output_subdir": "csd_1d_onset_baseline",
    # Batch root runs should not retain every trial matrix in Python RAM.
    "retain_epochs": False,
    # Behavior exploration uses one averaged-column image per physical probe.
    "save_column_plots": False,
    "save_mean_plot": True,
}

ANIMAL_POOLED_FOLDER = "animal_pooled"
ANIMAL_POOLED_DATASETS = (
    ("category", "play"),
    ("category", "aggression"),
    ("category", "exploration"),
    ("play_bout", "all"),
)


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
    # Memory-map trial stacks: channel preparation/epoching reads from disk
    # without first materializing both complete arrays in RAM.
    onset, x_on, y_on, t0_on, fs_on = load_raw(onset_folder, mmap_mode="r")
    offset, x_off, y_off, t0_off, fs_off = load_raw(offset_folder, mmap_mode="r")

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
    onset_baseline_window_ms: Tuple[float, float] = (-1000.0, 0.0),
    baseline_correct: bool = True,
    normalize: str = "none",
    retain_epochs: bool = True,
    event_chunk_size: int = 8,
) -> Dict[str, Any]:
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

    onset_view = onset_lfp[:, :, onset_start:onset_end]
    offset_view = offset_lfp[:, :, offset_start:offset_end]
    if onset_view.shape[2] != offset_view.shape[2]:
        raise ValueError(
            f"Onset/offset epoch lengths differ: {onset_view.shape[2]} vs "
            f"{offset_view.shape[2]}"
        )

    work_dtype = np.result_type(onset_lfp.dtype, offset_lfp.dtype, np.float32)
    onset_baseline = onset_lfp[:, :, baseline_start:baseline_end].mean(
        axis=2, keepdims=True
    ).astype(work_dtype, copy=False)
    if retain_epochs:
        onset_epochs = onset_view.astype(work_dtype, copy=True)
        offset_epochs = offset_view.astype(work_dtype, copy=True)
        if baseline_correct:
            onset_epochs -= onset_baseline
            offset_epochs -= onset_baseline
        onset_evoked = onset_epochs.mean(axis=0)
        offset_evoked = offset_epochs.mean(axis=0)
    else:
        if event_chunk_size < 1:
            raise ValueError("event_chunk_size must be at least 1")
        epoch_shape = onset_view.shape[1:]
        onset_sum = np.zeros(epoch_shape, dtype=np.float64)
        offset_sum = np.zeros(epoch_shape, dtype=np.float64)
        for first in range(0, onset_view.shape[0], event_chunk_size):
            last = min(first + event_chunk_size, onset_view.shape[0])
            onset_chunk = onset_view[first:last].astype(work_dtype, copy=True)
            offset_chunk = offset_view[first:last].astype(work_dtype, copy=True)
            if baseline_correct:
                baseline_chunk = onset_baseline[first:last]
                onset_chunk -= baseline_chunk
                offset_chunk -= baseline_chunk
            onset_sum += onset_chunk.sum(axis=0, dtype=np.float64)
            offset_sum += offset_chunk.sum(axis=0, dtype=np.float64)
        onset_evoked = onset_sum / onset_view.shape[0]
        offset_evoked = offset_sum / offset_view.shape[0]
        onset_epochs = None
        offset_epochs = None

    onset_evoked = apply_normalize(onset_evoked, normalize)
    offset_evoked = apply_normalize(offset_evoked, normalize)
    pre_samples = int(round(float(t_pre) * fs / 1000.0))
    time_ms = (np.arange(onset_view.shape[2]) - pre_samples) / fs * 1000.0
    step = max(1, int(round(step_ms * fs / 1000.0)))
    time_indices = np.arange(0, onset_view.shape[2], step)

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


def average_kcsd_columns(
    csd: Dict[float, np.ndarray],
    x: np.ndarray,
    y: np.ndarray,
    evoked: np.ndarray,
) -> Dict[str, np.ndarray]:
    """
    Interpolate staggered x columns onto their overlapping depths and average.

    Columns receive equal weight. The common depth axis is the union of
    observed electrode depths inside the range covered by every column, so no
    extrapolation is used.
    """
    if not csd:
        raise ValueError("Cannot average an empty kCSD column dictionary")

    column_data = []
    for ux in sorted(csd):
        indices = np.where(x == ux)[0]
        y_column = np.asarray(y[indices], dtype=float)
        order = np.argsort(y_column)
        y_column = y_column[order]
        csd_column = np.asarray(csd[ux], dtype=float)[order]
        voltage_column = np.asarray(evoked[indices], dtype=float)[order]
        if len(np.unique(y_column)) != len(y_column):
            raise ValueError(f"Duplicate y coordinates within x={ux}")
        column_data.append((float(ux), y_column, csd_column, voltage_column))

    overlap_min = max(data[1].min() for data in column_data)
    overlap_max = min(data[1].max() for data in column_data)
    if overlap_min >= overlap_max:
        raise ValueError(
            "The probe columns have no overlapping y range and cannot be averaged"
        )

    common_y = np.unique(
        np.concatenate(
            [
                data[1][(data[1] >= overlap_min) & (data[1] <= overlap_max)]
                for data in column_data
            ]
        )
    )
    if common_y.size < 2:
        raise ValueError("Fewer than two shared depth positions for column averaging")

    csd_interpolated = []
    voltage_interpolated = []
    for _ux, y_column, csd_column, voltage_column in column_data:
        csd_interpolated.append(
            interp1d(y_column, csd_column, axis=0, bounds_error=True)(common_y)
        )
        voltage_interpolated.append(
            interp1d(y_column, voltage_column, axis=0, bounds_error=True)(common_y)
        )

    csd_stack = np.stack(csd_interpolated, axis=0)
    voltage_stack = np.stack(voltage_interpolated, axis=0)
    n_columns = len(column_data)
    if n_columns > 1:
        csd_std = csd_stack.std(axis=0, ddof=1)
    else:
        csd_std = np.zeros_like(csd_stack[0])

    return {
        "csd_mean": csd_stack.mean(axis=0),
        "csd_std": csd_std,
        "voltage_mean": voltage_stack.mean(axis=0),
        "y": common_y,
        "x_center": np.asarray(
            np.mean([data[0] for data in column_data]), dtype=float
        ),
        "n_columns_contributing": np.full(
            common_y.shape, n_columns, dtype=np.int16
        ),
    }


def _save_averaged_column_outputs(
    folder: PathLike,
    averaged: Dict[str, np.ndarray],
    time_ms: np.ndarray,
    time_indices: np.ndarray,
    *,
    output_subdir: str,
) -> str:
    output = os.path.join(os.fspath(folder), output_subdir)
    os.makedirs(output, exist_ok=True)
    np.save(os.path.join(output, "csd_mean_columns.npy"), averaged["csd_mean"])
    np.save(os.path.join(output, "csd_std_columns.npy"), averaged["csd_std"])
    np.save(os.path.join(output, "y_mean_columns.npy"), averaged["y"])
    np.save(
        os.path.join(output, "n_columns_contributing.npy"),
        averaged["n_columns_contributing"],
    )
    np.save(os.path.join(output, "V_mean_columns.npy"), averaged["voltage_mean"])
    np.save(os.path.join(output, "x_center.npy"), averaged["x_center"])
    np.save(os.path.join(output, "time_mean_columns.npy"), time_ms[time_indices])
    np.save(os.path.join(output, "time_voltage_mean_columns.npy"), time_ms)
    return output


def _plot_averaged_alignment(
    label: str,
    name: str,
    folder: PathLike,
    output_subdir: str,
    averaged: Dict[str, np.ndarray],
    time_ms: np.ndarray,
    time_indices: np.ndarray,
    *,
    clim: Tuple[float, float],
    show: bool,
    save_plot: bool,
) -> None:
    save_path = None
    if save_plot:
        figure_dir = os.path.join(os.fspath(folder), output_subdir, "figures")
        os.makedirs(figure_dir, exist_ok=True)
        save_path = os.path.join(figure_dir, "kcsd_xmean.png")
    plot_column(
        averaged["csd_mean"],
        time_ms,
        time_indices,
        averaged["y"],
        title=f"{name} — {label} — mean across probe columns",
        cbar_label="mean kCSD (a.u.)",
        clim=clim,
        save_path=save_path,
        show=show,
    )
    if save_path is not None:
        print(f"    Saved averaged figure: {save_path}", flush=True)


def _remove_individual_column_figures(
    folder: PathLike,
    output_subdir: str,
) -> int:
    """Remove legacy per-x PNGs while retaining the averaged xmean image."""
    figure_dir = os.path.join(os.fspath(folder), output_subdir, "figures")
    if not os.path.isdir(figure_dir):
        return 0
    removed = 0
    for name in os.listdir(figure_dir):
        if (
            name.startswith("kcsd_x")
            and name.lower().endswith(".png")
            and name != "kcsd_xmean.png"
        ):
            os.remove(os.path.join(figure_dir, name))
            removed += 1
    return removed


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
        if save_path is not None:
            print(f"    Saved figure: {save_path}", flush=True)


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
    ``onset_baseline_window_ms``. The default is -1000…0 ms relative to onset;
    response epochs run from -1000 to +2000 ms.
    kCSD itself remains unchanged and receives the two prepared evoked LFPs.
    """
    valid_overrides = set(DatasetConfig.__dataclass_fields__) | {
        "onset_baseline_window_ms",
        "retain_epochs",
        "event_chunk_size",
        "save_column_plots",
        "save_mean_plot",
    }
    unknown = sorted(set(overrides) - valid_overrides)
    if unknown:
        raise TypeError(f"Unknown paired kCSD option(s): {', '.join(unknown)}")

    settings = dict(BEHAVIOR_ALIGNED_CSD_DEFAULTS)
    settings.update(overrides)
    baseline_window = tuple(settings.pop("onset_baseline_window_ms"))
    retain_epochs = bool(settings.pop("retain_epochs"))
    event_chunk_size = int(settings.pop("event_chunk_size", 8))
    save_column_plots = bool(settings.pop("save_column_plots"))
    save_mean_plot = bool(settings.pop("save_mean_plot"))
    if (
        len(baseline_window) != 2
        or baseline_window[0] >= baseline_window[1]
        or baseline_window[1] > 0
    ):
        raise ValueError(
            "onset_baseline_window_ms must be (start, end) with start < end <= 0"
        )

    resolved_name = name or os.path.relpath(
        os.fspath(onset_folder), os.path.commonpath((onset_folder, offset_folder))
    )
    started = time.perf_counter()
    print(f"  Loading paired raw LFP: {resolved_name}", flush=True)
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
    print(
        f"  Preparing {n_events} paired events "
        f"({cfg.t_pre:g} ms before, {cfg.t_post:g} ms after)",
        flush=True,
    )
    onset, x, y, channel_meta = prepare_channels(
        paired["onset"], paired["x"], paired["y"], cfg
    )
    offset, x_offset, y_offset, offset_channel_meta = prepare_channels(
        paired["offset"], paired["x"], paired["y"], cfg
    )
    if not np.array_equal(x, x_offset, equal_nan=True) or not np.array_equal(
        y, y_offset, equal_nan=True
    ):
        raise RuntimeError("Onset/offset channel preparation produced different geometry")
    if channel_meta != offset_channel_meta:
        raise RuntimeError("Onset/offset channel preparation settings differ")

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
        retain_epochs=retain_epochs,
        event_chunk_size=event_chunk_size,
    )
    print(
        f"  Evoked LFP ready: channels={len(y)}, "
        f"CSD frames={len(prepared['time_indices'])}",
        flush=True,
    )
    metadata: Dict[str, Any] = {
        "name": resolved_name,
        "onset_folder": os.path.abspath(os.fspath(onset_folder)),
        "offset_folder": os.path.abspath(os.fspath(offset_folder)),
        "baseline_source": "onset",
        "onset_baseline_window_ms": list(baseline_window),
        "baseline_correct": bool(cfg.baseline_correct),
        "retain_epochs": retain_epochs,
        "event_chunk_size": event_chunk_size,
        "save_column_plots": save_column_plots,
        "save_mean_plot": save_mean_plot,
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
        "onset_evoked": prepared["onset_evoked"],
        "offset_evoked": prepared["offset_evoked"],
        "onset_baseline": prepared["onset_baseline"],
        "time_ms": prepared["time_ms"],
        "time_indices": prepared["time_indices"],
        "onset_csd_kcsd": None,
        "offset_csd_kcsd": None,
        "onset_csd_mean_columns": None,
        "offset_csd_mean_columns": None,
        "onset_csd_fd": None,
        "offset_csd_fd": None,
    }
    if retain_epochs:
        result["onset_epochs"] = prepared["onset_epochs"]
        result["offset_epochs"] = prepared["offset_epochs"]

    if cfg.compute_kcsd:
        csd_kwargs = {
            "y_scale": cfg.y_scale,
            "h": cfg.h,
            "n_src_factor": cfg.n_src_factor,
            "smooth_sigma_time": cfg.smooth_sigma_time,
            "sort_by_depth": cfg.sort_by_depth,
        }
        print(f"  Computing onset kCSD: {resolved_name}", flush=True)
        onset_csd = compute_kcsd_1d_columns(
            prepared["onset_evoked"],
            x,
            y,
            prepared["time_indices"],
            progress_label=f"{resolved_name} | onset",
            **csd_kwargs,
        )
        print(
            f"  Computing offset kCSD: {resolved_name} (same onset baseline)",
            flush=True,
        )
        offset_csd = compute_kcsd_1d_columns(
            prepared["offset_evoked"],
            x,
            y,
            prepared["time_indices"],
            progress_label=f"{resolved_name} | offset",
            **csd_kwargs,
        )
        result["onset_csd_kcsd"] = onset_csd
        result["offset_csd_kcsd"] = offset_csd
        print("  Interpolating and averaging probe columns", flush=True)
        onset_averaged = average_kcsd_columns(
            onset_csd, x, y, prepared["onset_evoked"]
        )
        offset_averaged = average_kcsd_columns(
            offset_csd, x, y, prepared["offset_evoked"]
        )
        result["onset_csd_mean_columns"] = onset_averaged
        result["offset_csd_mean_columns"] = offset_averaged
        metadata["n_columns_averaged"] = int(
            onset_averaged["n_columns_contributing"][0]
        )
        metadata["x_center"] = float(onset_averaged["x_center"])

        if cfg.save:
            for alignment, folder, csd, evoked in (
                ("onset", onset_folder, onset_csd, prepared["onset_evoked"]),
                ("offset", offset_folder, offset_csd, prepared["offset_evoked"]),
            ):
                print(f"  Saving {alignment} matrices: {folder}", flush=True)
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
            print("  Saving averaged-column matrices", flush=True)
            _save_averaged_column_outputs(
                onset_folder,
                onset_averaged,
                prepared["time_ms"],
                prepared["time_indices"],
                output_subdir=cfg.output_subdir,
            )
            _save_averaged_column_outputs(
                offset_folder,
                offset_averaged,
                prepared["time_ms"],
                prepared["time_indices"],
                output_subdir=cfg.output_subdir,
            )

        if (cfg.plot or cfg.save_plots) and save_column_plots:
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

        if (cfg.plot or cfg.save_plots) and save_mean_plot:
            if cfg.save_plots and not save_column_plots:
                removed = sum(
                    _remove_individual_column_figures(folder, cfg.output_subdir)
                    for folder in (onset_folder, offset_folder)
                )
                if removed:
                    print(
                        f"  Removed {removed} legacy per-column PNGs; "
                        "averaged xmean figures are the primary view",
                        flush=True,
                    )
            shared_abs_max = max(
                float(np.nanmax(np.abs(onset_averaged["csd_mean"]))),
                float(np.nanmax(np.abs(offset_averaged["csd_mean"]))),
            )
            if not np.isfinite(shared_abs_max) or shared_abs_max <= 0:
                shared_abs_max = 1.0
            shared_clim = (-shared_abs_max, shared_abs_max)
            print(
                "  Saving averaged onset/offset figures with shared color limits",
                flush=True,
            )
            _plot_averaged_alignment(
                "onset",
                resolved_name,
                onset_folder,
                cfg.output_subdir,
                onset_averaged,
                prepared["time_ms"],
                prepared["time_indices"],
                clim=shared_clim,
                show=cfg.plot,
                save_plot=cfg.save_plots,
            )
            _plot_averaged_alignment(
                "offset",
                resolved_name,
                offset_folder,
                cfg.output_subdir,
                offset_averaged,
                prepared["time_ms"],
                prepared["time_indices"],
                clim=shared_clim,
                show=cfg.plot,
                save_plot=cfg.save_plots,
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

    print(
        f"  Finished {resolved_name} in {time.perf_counter() - started:.1f} s",
        flush=True,
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
        # Animal-level PSTHs are written back under this root. They are not
        # another session and must not be estimated or averaged again.
        if ANIMAL_POOLED_FOLDER in rel_parts:
            continue
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


def _safe_figure_name_part(value: str) -> str:
    """Make one metadata value safe and sortable in Windows filenames."""
    value = re.sub(r'[<>:"/\\|?*]+', "_", str(value))
    value = re.sub(r"\s+", "_", value.strip())
    return value or "unknown"


def organize_behavior_aligned_figures(
    root: PathLike,
    *,
    output_folder_name: str = "sorted_figures",
    csd_output_subdir: str = "csd_1d_onset_baseline",
    refresh: bool = True,
) -> Dict[str, Any]:
    """
    Copy every behavior kCSD PNG into session-first and behavior-first views.

    Original figures are never changed. With ``refresh=True``, only the
    generated ``sorted_figures`` folder is rebuilt, preventing stale copies.
    """
    root = os.path.abspath(os.fspath(root))
    pairs = discover_behavior_aligned_pairs(root)
    if not pairs:
        raise FileNotFoundError(
            f"No paired onset/offset CSD-ready folders found under {root}"
        )

    output_root = os.path.join(root, output_folder_name)
    if refresh and os.path.isdir(output_root):
        shutil.rmtree(output_root)
    os.makedirs(output_root, exist_ok=True)

    rows: List[Dict[str, str]] = []
    destinations: Dict[str, str] = {}
    missing_figure_folders = 0
    print(
        f"Organizing PNGs from {len(pairs)} paired datasets into {output_root}",
        flush=True,
    )

    for pair_number, pair in enumerate(pairs, start=1):
        parts = pair["name"].replace("\\", "/").split("/")
        if len(parts) < 3:
            print(
                f"[{pair_number}/{len(pairs)}] Skipping unrecognized path: "
                f"{pair['name']}",
                flush=True,
            )
            continue
        scope, behavior, session = parts[:3]
        probe = parts[3] if len(parts) > 3 else "NPX1"
        print(
            f"[{pair_number}/{len(pairs)}] Sorting figures: "
            f"{scope}/{behavior} | {session} | {probe}",
            flush=True,
        )

        safe_scope = _safe_figure_name_part(scope)
        safe_behavior = _safe_figure_name_part(behavior)
        safe_session = _safe_figure_name_part(session)
        safe_probe = _safe_figure_name_part(probe)

        for alignment, leaf in (
            ("onset", pair["onset_folder"]),
            ("offset", pair["offset_folder"]),
        ):
            figure_dir = os.path.join(leaf, csd_output_subdir, "figures")
            if not os.path.isdir(figure_dir):
                missing_figure_folders += 1
                continue

            all_figure_names = sorted(
                name
                for name in os.listdir(figure_dir)
                if name.lower().endswith(".png") and name.startswith("kcsd_x")
            )
            # The averaged physical-probe view is the primary exploration
            # image. Fall back to legacy column PNGs only for older runs.
            if "kcsd_xmean.png" in all_figure_names:
                figure_names = ["kcsd_xmean.png"]
            else:
                figure_names = all_figure_names
            for figure_name in figure_names:
                source = os.path.join(figure_dir, figure_name)
                xpos = _safe_figure_name_part(
                    os.path.splitext(figure_name)[0].removeprefix("kcsd_")
                )
                session_name = (
                    f"{safe_session}__{safe_scope}__{safe_behavior}__"
                    f"{safe_probe}__{xpos}.png"
                )
                behavior_name = (
                    f"{safe_scope}__{safe_behavior}__{safe_session}__"
                    f"{safe_probe}__{xpos}.png"
                )
                by_session = os.path.join(
                    output_root, "by_session", alignment, session_name
                )
                by_behavior = os.path.join(
                    output_root, "by_behavior", alignment, behavior_name
                )

                for destination in (by_session, by_behavior):
                    prior_source = destinations.get(destination)
                    if prior_source is not None and prior_source != source:
                        raise RuntimeError(
                            "Figure-name collision: "
                            f"{prior_source} and {source} both map to {destination}"
                        )
                    destinations[destination] = source
                    os.makedirs(os.path.dirname(destination), exist_ok=True)
                    shutil.copy2(source, destination)
                    print(
                        f"    Saved sorted copy: "
                        f"{os.path.relpath(destination, output_root)}",
                        flush=True,
                    )

                rows.append(
                    {
                        "alignment": alignment,
                        "scope": scope,
                        "behavior": behavior,
                        "session": session,
                        "probe": probe,
                        "xpos": xpos,
                        "original_file": os.path.relpath(source, root),
                        "by_session_file": os.path.relpath(by_session, root),
                        "by_behavior_file": os.path.relpath(by_behavior, root),
                    }
                )

    index_path = os.path.join(output_root, "figure_index.csv")
    fieldnames = [
        "alignment",
        "scope",
        "behavior",
        "session",
        "probe",
        "xpos",
        "original_file",
        "by_session_file",
        "by_behavior_file",
    ]
    with open(index_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"Figure sorting complete: {len(rows)} source PNGs, "
        f"{2 * len(rows)} copies, index={index_path}",
        flush=True,
    )
    if missing_figure_folders:
        print(
            f"Note: {missing_figure_folders} alignment folders had no PNG "
            "(for example, plotting may have been disabled).",
            flush=True,
        )
    return {
        "output_root": output_root,
        "index_path": index_path,
        "n_source_figures": len(rows),
        "n_copies": 2 * len(rows),
        "missing_figure_folders": missing_figure_folders,
    }


def organize_animal_pooled_figures(
    root: PathLike,
    *,
    pooled_folder_name: str = ANIMAL_POOLED_FOLDER,
    output_folder_name: str = "sorted_figures",
    csd_output_subdir: str = "csd_1d_onset_baseline",
    refresh: bool = True,
    missing_ok: bool = False,
) -> Dict[str, Any]:
    """
    Copy animal-pooled kCSD PNGs into one onset folder and one offset folder.

    Original figures stay beside each pooled PSTH. Names follow the session
    organizer, with the animal key in the session position:

    ``by_animal/onset/Batch1RatA2__category__play__NPX1__xmean.png``
    ``by_behavior/onset/category__play__Batch1RatA2__NPX1__xmean.png``
    """
    root = os.path.abspath(os.fspath(root))
    pooled_root = os.path.join(root, pooled_folder_name)
    if not os.path.isdir(pooled_root):
        raise FileNotFoundError(
            f"No animal-pooled results at {pooled_root}. Run "
            "estimate_animal_pooled_kcsd first."
        )

    figure_dirs: List[str] = []
    for dirpath, _dirnames, filenames in os.walk(pooled_root):
        relative = os.path.relpath(dirpath, pooled_root).split(os.sep)
        if output_folder_name in relative:
            continue
        if (
            relative[-1] == "figures"
            and len(relative) >= 2
            and relative[-2] == csd_output_subdir
            and any(name.lower().endswith(".png") and name.startswith("kcsd_x") for name in filenames)
        ):
            figure_dirs.append(dirpath)
    figure_dirs.sort()

    output_root = os.path.join(pooled_root, output_folder_name)
    if not figure_dirs:
        message = (
            f"No animal-pooled kCSD PNGs under {pooled_root}. "
            "Estimation may have been run with save_plots=False."
        )
        if missing_ok:
            print(message, flush=True)
            return {
                "output_root": output_root,
                "index_path": None,
                "n_source_figures": 0,
                "n_copies": 0,
            }
        raise FileNotFoundError(message)

    if refresh and os.path.isdir(output_root):
        try:
            shutil.rmtree(output_root)
        except PermissionError as exc:
            print(
                f"Could not delete {output_root}: {exc}. "
                "Existing figure copies will be overwritten.",
                flush=True,
            )
    os.makedirs(output_root, exist_ok=True)

    rows: List[Dict[str, str]] = []
    destinations: Dict[str, str] = {}
    print(
        f"Copying animal-pooled kCSD figures from {len(figure_dirs)} "
        f"folders into {output_root}",
        flush=True,
    )
    for figure_dir in figure_dirs:
        relative = os.path.relpath(figure_dir, pooled_root).split(os.sep)
        leaf = relative[:-2]
        if len(leaf) not in (4, 5) or leaf[2] not in {"onset", "offset"}:
            print(f"  Skipping unrecognized figure folder: {figure_dir}", flush=True)
            continue
        scope, behavior, alignment, animal = leaf[:4]
        probe = leaf[4] if len(leaf) == 5 else "NPX1"
        all_figure_names = sorted(
            name
            for name in os.listdir(figure_dir)
            if name.lower().endswith(".png") and name.startswith("kcsd_x")
        )
        if "kcsd_xmean.png" in all_figure_names:
            figure_names = ["kcsd_xmean.png"]
        else:
            figure_names = all_figure_names
        safe_scope = _safe_figure_name_part(scope)
        safe_behavior = _safe_figure_name_part(behavior)
        safe_animal = _safe_figure_name_part(animal)
        safe_probe = _safe_figure_name_part(probe)
        for figure_name in figure_names:
            source = os.path.join(figure_dir, figure_name)
            xpos = _safe_figure_name_part(
                os.path.splitext(figure_name)[0].removeprefix("kcsd_")
            )
            animal_name = (
                f"{safe_animal}__{safe_scope}__{safe_behavior}__"
                f"{safe_probe}__{xpos}.png"
            )
            behavior_name = (
                f"{safe_scope}__{safe_behavior}__{safe_animal}__"
                f"{safe_probe}__{xpos}.png"
            )
            by_animal = os.path.join(output_root, "by_animal", alignment, animal_name)
            by_behavior = os.path.join(
                output_root, "by_behavior", alignment, behavior_name
            )
            for destination in (by_animal, by_behavior):
                prior_source = destinations.get(destination)
                if prior_source is not None and prior_source != source:
                    raise RuntimeError(
                        "Figure-name collision: "
                        f"{prior_source} and {source} both map to {destination}"
                    )
                destinations[destination] = source
                os.makedirs(os.path.dirname(destination), exist_ok=True)
                shutil.copy2(source, destination)
                print(
                    f"    Saved sorted copy: "
                    f"{os.path.relpath(destination, output_root)}",
                    flush=True,
                )
            rows.append(
                {
                    "alignment": alignment,
                    "scope": scope,
                    "behavior": behavior,
                    "animal": animal,
                    "probe": probe,
                    "xpos": xpos,
                    "original_file": os.path.relpath(source, root),
                    "by_animal_file": os.path.relpath(by_animal, root),
                    "by_behavior_file": os.path.relpath(by_behavior, root),
                }
            )

    index_path = os.path.join(output_root, "figure_index.csv")
    fieldnames = [
        "alignment",
        "scope",
        "behavior",
        "animal",
        "probe",
        "xpos",
        "original_file",
        "by_animal_file",
        "by_behavior_file",
    ]
    with open(index_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"Animal-pooled figure sorting complete: {len(rows)} source PNGs, "
        f"{2 * len(rows)} copies, index={index_path}",
        flush=True,
    )
    return {
        "output_root": output_root,
        "index_path": index_path,
        "n_source_figures": len(rows),
        "n_copies": 2 * len(rows),
    }


def _load_and_average_saved_columns(
    folder: PathLike,
    *,
    output_subdir: str,
) -> Tuple[Dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """Load existing per-column kCSD matrices; do not estimate kCSD again."""
    output = os.path.join(os.fspath(folder), output_subdir)
    if not os.path.isdir(output):
        raise FileNotFoundError(f"Missing saved kCSD folder: {output}")

    matches = []
    for name in os.listdir(output):
        match = re.fullmatch(r"csd_x(.+)\.npy", name)
        if match:
            try:
                matches.append((float(match.group(1)), match.group(1)))
            except ValueError:
                continue
    matches.sort()
    if not matches:
        raise FileNotFoundError(f"No saved csd_x*.npy matrices under {output}")

    csd: Dict[float, np.ndarray] = {}
    x_parts = []
    y_parts = []
    voltage_parts = []
    voltage_time = None
    for ux, tag in matches:
        csd[ux] = np.load(os.path.join(output, f"csd_x{tag}.npy"))
        y_column = np.load(os.path.join(output, f"y_x{tag}.npy")).reshape(-1)
        voltage = np.load(os.path.join(output, f"V_x{tag}.npy"))
        this_voltage_time = np.load(
            os.path.join(output, f"time_voltage_x{tag}.npy")
        ).reshape(-1)
        if voltage_time is None:
            voltage_time = this_voltage_time
        elif not np.array_equal(voltage_time, this_voltage_time):
            raise ValueError(f"Column voltage time axes differ under {output}")
        if csd[ux].shape[0] != len(y_column) or voltage.shape[0] != len(y_column):
            raise ValueError(f"Channel/depth mismatch for saved x={ux} under {output}")
        x_parts.append(np.full(len(y_column), ux, dtype=float))
        y_parts.append(y_column)
        voltage_parts.append(voltage)

    csd_time = np.load(os.path.join(output, "time.npy")).reshape(-1)
    if any(matrix.shape[1] != len(csd_time) for matrix in csd.values()):
        raise ValueError(f"Saved CSD time-axis mismatch under {output}")
    assert voltage_time is not None
    averaged = average_kcsd_columns(
        csd,
        np.concatenate(x_parts),
        np.concatenate(y_parts),
        np.concatenate(voltage_parts, axis=0),
    )
    return averaged, csd_time, voltage_time


def _save_postprocessed_average(
    folder: PathLike,
    averaged: Dict[str, np.ndarray],
    csd_time: np.ndarray,
    voltage_time: np.ndarray,
    *,
    output_subdir: str,
) -> None:
    output = os.path.join(os.fspath(folder), output_subdir)
    np.save(os.path.join(output, "csd_mean_columns.npy"), averaged["csd_mean"])
    np.save(os.path.join(output, "csd_std_columns.npy"), averaged["csd_std"])
    np.save(os.path.join(output, "y_mean_columns.npy"), averaged["y"])
    np.save(
        os.path.join(output, "n_columns_contributing.npy"),
        averaged["n_columns_contributing"],
    )
    np.save(os.path.join(output, "V_mean_columns.npy"), averaged["voltage_mean"])
    np.save(os.path.join(output, "x_center.npy"), averaged["x_center"])
    np.save(os.path.join(output, "time_mean_columns.npy"), csd_time)
    np.save(os.path.join(output, "time_voltage_mean_columns.npy"), voltage_time)


def average_saved_behavior_kcsd_from_root(
    root: PathLike,
    *,
    output_subdir: str = "csd_1d_onset_baseline",
    save_plots: bool = True,
    remove_column_plots: bool = True,
    stop_on_error: bool = False,
) -> List[Dict[str, Any]]:
    """
    Average already-saved per-column matrices without rerunning kCSD.

    This is the intended migration path for completed notebook runs.
    """
    root = os.path.abspath(os.fspath(root))
    pairs = discover_behavior_aligned_pairs(root)
    if not pairs:
        raise FileNotFoundError(
            f"No paired onset/offset CSD-ready folders found under {root}"
        )

    print(
        f"Post-processing {len(pairs)} saved datasets; kCSD will NOT be rerun",
        flush=True,
    )
    results: List[Dict[str, Any]] = []
    for pair_number, pair in enumerate(pairs, start=1):
        print(
            f"\n[{pair_number}/{len(pairs)}] Averaging saved matrices: "
            f"{pair['name']}",
            flush=True,
        )
        try:
            onset_average, onset_time, onset_voltage_time = (
                _load_and_average_saved_columns(
                    pair["onset_folder"], output_subdir=output_subdir
                )
            )
            offset_average, offset_time, offset_voltage_time = (
                _load_and_average_saved_columns(
                    pair["offset_folder"], output_subdir=output_subdir
                )
            )
            if not np.array_equal(onset_time, offset_time):
                raise ValueError("Onset and offset saved CSD time axes differ")

            _save_postprocessed_average(
                pair["onset_folder"],
                onset_average,
                onset_time,
                onset_voltage_time,
                output_subdir=output_subdir,
            )
            _save_postprocessed_average(
                pair["offset_folder"],
                offset_average,
                offset_time,
                offset_voltage_time,
                output_subdir=output_subdir,
            )

            shared_abs_max = max(
                float(np.nanmax(np.abs(onset_average["csd_mean"]))),
                float(np.nanmax(np.abs(offset_average["csd_mean"]))),
            )
            if not np.isfinite(shared_abs_max) or shared_abs_max <= 0:
                shared_abs_max = 1.0
            shared_clim = (-shared_abs_max, shared_abs_max)

            removed = 0
            if remove_column_plots:
                removed = sum(
                    _remove_individual_column_figures(folder, output_subdir)
                    for folder in (pair["onset_folder"], pair["offset_folder"])
                )
            if save_plots:
                selected_indices = np.arange(len(onset_time))
                _plot_averaged_alignment(
                    "onset",
                    pair["name"],
                    pair["onset_folder"],
                    output_subdir,
                    onset_average,
                    onset_time,
                    selected_indices,
                    clim=shared_clim,
                    show=False,
                    save_plot=True,
                )
                _plot_averaged_alignment(
                    "offset",
                    pair["name"],
                    pair["offset_folder"],
                    output_subdir,
                    offset_average,
                    offset_time,
                    selected_indices,
                    clim=shared_clim,
                    show=False,
                    save_plot=True,
                )
            n_columns = int(onset_average["n_columns_contributing"][0])
            print(
                f"  Averaged {n_columns} columns; removed {removed} "
                "individual-column PNGs",
                flush=True,
            )
            results.append(
                {
                    "name": pair["name"],
                    "n_columns": n_columns,
                    "removed_column_plots": removed,
                }
            )
        except Exception as exc:  # noqa: BLE001 - batch post-processing continues
            print(f"  ERROR: {exc}", flush=True)
            if stop_on_error:
                raise
            results.append({"name": pair["name"], "error": str(exc)})

    n_failed = sum("error" in result for result in results)
    print(
        f"\nSaved-matrix post-processing finished: "
        f"{len(results) - n_failed} successful, {n_failed} failed",
        flush=True,
    )
    return results


def _order_channels(
    x: np.ndarray,
    y: np.ndarray,
    onset: np.ndarray,
    offset: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    order = np.lexsort((np.asarray(y, dtype=float), np.asarray(x, dtype=float)))
    return x[order], y[order], onset[order], offset[order]


def _same_electrode_grid(records: Sequence[Dict[str, Any]]) -> bool:
    reference_x = records[0]["x"]
    reference_y = records[0]["y"]
    return all(
        record["x"].shape == reference_x.shape
        and np.allclose(record["x"], reference_x)
        and np.allclose(record["y"], reference_y)
        for record in records[1:]
    )


def _column_profiles(
    x: np.ndarray,
    y: np.ndarray,
    onset: np.ndarray,
    offset: np.ndarray,
) -> Dict[float, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    profiles = {}
    for ux in np.unique(np.asarray(x, dtype=float)):
        indices = np.where(x == ux)[0]
        order = np.argsort(y[indices])
        yy = np.asarray(y[indices][order], dtype=float)
        onset_column = np.asarray(onset[indices][order], dtype=float)
        offset_column = np.asarray(offset[indices][order], dtype=float)
        unique_y, inverse = np.unique(yy, return_inverse=True)
        if len(unique_y) < len(yy):
            onset_sum = np.zeros((len(unique_y), onset_column.shape[1]), dtype=float)
            offset_sum = np.zeros_like(onset_sum)
            counts = np.zeros(len(unique_y), dtype=float)
            for index, group in enumerate(inverse):
                onset_sum[group] += onset_column[index]
                offset_sum[group] += offset_column[index]
                counts[group] += 1.0
            onset_column = onset_sum / counts[:, None]
            offset_column = offset_sum / counts[:, None]
            yy = unique_y
        profiles[float(ux)] = (yy, onset_column, offset_column)
    return profiles


def _weighted_mean(records: Sequence[Dict[str, Any]], key: str) -> np.ndarray:
    total = float(sum(int(record["n_events"]) for record in records))
    accumulator = np.zeros_like(records[0][key], dtype=float)
    for record in records:
        accumulator += int(record["n_events"]) * np.asarray(record[key], dtype=float)
    return accumulator / total


def trial_weighted_psth(records: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Average session PSTHs with weight equal to the number of trials.

    This is the mean that concatenation of the trials would produce. Sessions
    with the same electrode coordinates are averaged directly. When depths
    differ between days, each probe column is interpolated onto depths inside
    the range covered by every session, and those aligned PSTHs are averaged.
    Values are not extrapolated past a session's own electrodes.
    """
    if not records:
        raise ValueError("No session PSTHs to combine")
    ordered = []
    for record in records:
        x, y, onset, offset = _order_channels(
            np.asarray(record["x"], dtype=float),
            np.asarray(record["y"], dtype=float),
            np.asarray(record["onset_evoked"], dtype=float),
            np.asarray(record["offset_evoked"], dtype=float),
        )
        ordered.append({**record, "x": x, "y": y, "onset_evoked": onset, "offset_evoked": offset})
    reference_time = np.asarray(ordered[0]["time_ms"], dtype=float)
    for record in ordered[1:]:
        if not np.allclose(record["time_ms"], reference_time):
            raise ValueError(
                f"Time axis for {record['session']} does not match "
                f"{ordered[0]['session']}"
            )
    total_events = int(sum(int(record["n_events"]) for record in ordered))
    if _same_electrode_grid(ordered):
        return {
            "x": ordered[0]["x"],
            "y": ordered[0]["y"],
            "onset_evoked": _weighted_mean(ordered, "onset_evoked"),
            "offset_evoked": _weighted_mean(ordered, "offset_evoked"),
            "time_ms": reference_time,
            "n_events": total_events,
            "aligned_by": "identical_electrodes",
        }

    column_sets = [_column_profiles(r["x"], r["y"], r["onset_evoked"], r["offset_evoked"]) for r in ordered]
    shared_x = set(column_sets[0])
    for columns in column_sets[1:]:
        shared_x &= set(columns)
    if not shared_x:
        raise ValueError("Sessions have no probe column in common")
    x_values = []
    y_values = []
    onset_rows = []
    offset_rows = []
    for ux in sorted(shared_x):
        profiles = [columns[ux] for columns in column_sets]
        overlap_min = max(float(profile[0][0]) for profile in profiles)
        overlap_max = min(float(profile[0][-1]) for profile in profiles)
        if overlap_min >= overlap_max:
            raise ValueError(f"Column x={ux:g} has no shared depth range")
        observed = [
            profile[0][(profile[0] >= overlap_min) & (profile[0] <= overlap_max)]
            for profile in profiles
        ]
        common_y = np.unique(np.concatenate(observed))
        if common_y.size < 2:
            raise ValueError(f"Column x={ux:g} has fewer than two shared depths")
        onset_acc = np.zeros((common_y.size, reference_time.size), dtype=float)
        offset_acc = np.zeros_like(onset_acc)
        for record, (yy, onset_column, offset_column) in zip(ordered, profiles):
            onset_acc += record["n_events"] * interp1d(
                yy, onset_column, axis=0, bounds_error=True
            )(common_y)
            offset_acc += record["n_events"] * interp1d(
                yy, offset_column, axis=0, bounds_error=True
            )(common_y)
        onset_acc /= total_events
        offset_acc /= total_events
        x_values.append(np.full(common_y.size, ux, dtype=float))
        y_values.append(common_y)
        onset_rows.append(onset_acc)
        offset_rows.append(offset_acc)
    return {
        "x": np.concatenate(x_values),
        "y": np.concatenate(y_values),
        "onset_evoked": np.vstack(onset_rows),
        "offset_evoked": np.vstack(offset_rows),
        "time_ms": reference_time,
        "n_events": total_events,
        "aligned_by": "shared_column_depths",
    }


def _session_psth(
    onset_folder: str,
    offset_folder: str,
    *,
    epoch_settings: Dict[str, Any],
) -> Dict[str, Any]:
    paired = load_paired_raw(onset_folder, offset_folder)
    fs = float(paired["fs_loaded"]) if paired["fs_loaded"] is not None else float(
        epoch_settings["fs"]
    )
    prepared = epoch_paired_with_onset_baseline(
        paired["onset"],
        paired["offset"],
        paired["t0_onset"],
        paired["t0_offset"],
        fs,
        float(epoch_settings["t_pre"]),
        float(epoch_settings["t_post"]),
        float(epoch_settings["step_ms"]),
        onset_baseline_window_ms=tuple(epoch_settings["onset_baseline_window_ms"]),
        baseline_correct=bool(epoch_settings["baseline_correct"]),
        normalize=str(epoch_settings["normalize"]),
        retain_epochs=False,
    )
    return {
        "x": np.asarray(paired["x"], dtype=float),
        "y": np.asarray(paired["y"], dtype=float),
        "onset_evoked": prepared["onset_evoked"],
        "offset_evoked": prepared["offset_evoked"],
        "time_ms": prepared["time_ms"],
        "n_events": int(paired["onset"].shape[0]),
        "fs": fs,
    }


def _write_pooled_leaf(
    folder: str,
    evoked: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    time_ms: np.ndarray,
    fs: float,
) -> None:
    os.makedirs(folder, exist_ok=True)
    pre_samples = int(np.argmin(np.abs(time_ms)))
    if not np.isclose(time_ms[pre_samples], 0.0):
        raise ValueError("Pooled PSTH time axis does not contain 0 ms")
    np.save(os.path.join(folder, "lfp.npy"), evoked.astype(np.float32)[None, :, :])
    np.save(os.path.join(folder, "x_coords.npy"), np.asarray(x, dtype=float))
    np.save(os.path.join(folder, "y_coords.npy"), np.asarray(y, dtype=float))
    np.save(os.path.join(folder, "T0.npy"), np.array(pre_samples))
    np.save(os.path.join(folder, "fs.npy"), np.array(fs))
    np.save(os.path.join(folder, "time_ms.npy"), np.asarray(time_ms, dtype=float))
    np.save(os.path.join(folder, "event_ids.npy"), np.array([1], dtype=np.int64))


def estimate_animal_pooled_kcsd(
    root: PathLike,
    *,
    output_folder_name: str = ANIMAL_POOLED_FOLDER,
    datasets: Sequence[Tuple[str, str]] = ANIMAL_POOLED_DATASETS,
    stop_on_error: bool = False,
    refresh: bool = True,
    **overrides: Any,
) -> List[Dict[str, Any]]:
    """
    Pool trials from sessions of one animal, then estimate paired kCSD.

    Only play bouts and the play, aggression, and exploration categories are
    included. Individual behaviors are left as session-level exports. The
    pooled PSTH weights each session by its trial count, which is the same
    mean as concatenating the trials before averaging. kCSD is fit to that
    pooled PSTH, not to each session and then averaged.
    """
    root = os.path.abspath(os.fspath(root))
    selected = {tuple(item) for item in datasets}
    grouped: Dict[Tuple[str, str, str, str], List[Dict[str, str]]] = {}
    for pair in discover_behavior_aligned_pairs(root):
        parts = pair["name"].replace("\\", "/").split("/")
        if output_folder_name in parts or len(parts) < 3:
            continue
        scope, behavior, session = parts[:3]
        if (scope, behavior) not in selected:
            continue
        probe = parts[3] if len(parts) > 3 else ""
        try:
            animal_key = animal_key_from_session(session)
        except ValueError as exc:
            print(f"Skipping {pair['name']}: {exc}", flush=True)
            continue
        grouped.setdefault((scope, behavior, animal_key, probe), []).append(
            {**pair, "session": session, "probe": probe}
        )
    if not grouped:
        raise FileNotFoundError(
            "No play-bout or category onset/offset exports were found under "
            f"{root}. Individual behavior folders are not pooled."
        )

    settings = dict(BEHAVIOR_ALIGNED_CSD_DEFAULTS)
    settings.update(overrides)
    epoch_settings = {
        "fs": float(settings.get("fs", 2500.0)),
        "t_pre": float(settings["t_pre"]),
        "t_post": float(settings["t_post"]),
        "step_ms": float(settings["step_ms"]),
        "onset_baseline_window_ms": tuple(settings["onset_baseline_window_ms"]),
        "baseline_correct": bool(settings["baseline_correct"]),
        "normalize": str(settings["normalize"]),
    }
    output_root = os.path.join(root, output_folder_name)
    if refresh and os.path.isdir(output_root):
        try:
            shutil.rmtree(output_root)
        except PermissionError as exc:
            print(
                f"Could not delete {output_root}: {exc}. "
                "Existing pooled files will be overwritten.",
                flush=True,
            )
    os.makedirs(output_root, exist_ok=True)

    print(
        f"Pooling trials within animal for {len(grouped)} "
        "play-bout/category probe groups, then estimating kCSD",
        flush=True,
    )
    results: List[Dict[str, Any]] = []
    index_rows: List[Dict[str, Any]] = []
    for group_number, key in enumerate(sorted(grouped), start=1):
        scope, behavior, animal_key, probe = key
        members = sorted(grouped[key], key=lambda item: item["session"])
        label = f"{scope}/{behavior} {animal_key} {probe or 'NPX1'}"
        print(
            f"\n[{group_number}/{len(grouped)}] {label}: "
            f"{len(members)} sessions",
            flush=True,
        )
        try:
            psths = []
            for member in members:
                print(f"  Reading trials: {member['session']}", flush=True)
                psth = _session_psth(
                    member["onset_folder"],
                    member["offset_folder"],
                    epoch_settings=epoch_settings,
                )
                psth["session"] = member["session"]
                psths.append(psth)
                print(f"    {psth['n_events']} trials", flush=True)
            pooled = trial_weighted_psth(psths)
            tail = [animal_key] if not probe else [animal_key, probe]
            onset_folder = os.path.join(output_root, scope, behavior, "onset", *tail)
            offset_folder = os.path.join(output_root, scope, behavior, "offset", *tail)
            _write_pooled_leaf(
                onset_folder,
                pooled["onset_evoked"],
                pooled["x"],
                pooled["y"],
                pooled["time_ms"],
                psths[0]["fs"],
            )
            _write_pooled_leaf(
                offset_folder,
                pooled["offset_evoked"],
                pooled["x"],
                pooled["y"],
                pooled["time_ms"],
                psths[0]["fs"],
            )
            source_rows = [
                {"session": item["session"], "n_events": int(item["n_events"])}
                for item in psths
            ]
            pooling_record = {
                "animal_key": animal_key,
                "scope": scope,
                "behavior": behavior,
                "probe": probe or "NPX1",
                "n_sessions": len(psths),
                "n_events": int(pooled["n_events"]),
                "aligned_by": pooled["aligned_by"],
                "sessions": source_rows,
            }
            for folder in (onset_folder, offset_folder):
                with open(os.path.join(folder, "animal_pooling.json"), "w", encoding="utf-8") as stream:
                    json.dump(pooling_record, stream, indent=2)
            print(
                f"  Pooled {pooled['n_events']} trials by {pooled['aligned_by']}; "
                "estimating kCSD",
                flush=True,
            )
            result = run_behavior_aligned_kcsd(
                onset_folder,
                offset_folder,
                name=f"{output_folder_name}/{scope}/{behavior}/{animal_key}"
                + (f"/{probe}" if probe else ""),
                **overrides,
            )
            result["animal_pooling"] = pooling_record
            results.append(result)
            index_rows.append(
                {
                    "scope": scope,
                    "behavior": behavior,
                    "animal_key": animal_key,
                    "probe": probe or "NPX1",
                    "n_sessions": len(psths),
                    "n_events": int(pooled["n_events"]),
                    "sessions": "; ".join(
                        f"{item['session']} ({item['n_events']})" for item in source_rows
                    ),
                    "aligned_by": pooled["aligned_by"],
                    "onset_folder": os.path.relpath(onset_folder, root),
                }
            )
        except Exception as exc:  # noqa: BLE001 - one animal should not stop the rest
            print(f"  ERROR {label}: {exc}", flush=True)
            if stop_on_error:
                raise
            results.append({"name": label, "error": str(exc)})
    index_path = os.path.join(output_root, "index.csv")
    with open(index_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "scope",
                "behavior",
                "animal_key",
                "probe",
                "n_sessions",
                "n_events",
                "sessions",
                "aligned_by",
                "onset_folder",
            ],
        )
        writer.writeheader()
        writer.writerows(index_rows)
    n_failed = sum("error" in result for result in results)
    print(
        f"\nAnimal-pooled kCSD finished: {len(results) - n_failed} successful, "
        f"{n_failed} failed. Index: {index_path}",
        flush=True,
    )
    organize_animal_pooled_figures(
        root,
        pooled_folder_name=output_folder_name,
        csd_output_subdir=str(settings.get("output_subdir", "csd_1d_onset_baseline")),
        refresh=refresh,
        missing_ok=True,
    )
    return results


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

    print(f"Found {len(pairs)} paired datasets under {root}", flush=True)
    batch_started = time.perf_counter()
    results: List[Dict[str, Any]] = []
    for pair_number, pair in enumerate(pairs, start=1):
        pair_started = time.perf_counter()
        print(
            f"\n[{pair_number}/{len(pairs)}] Processing {pair['name']}",
            flush=True,
        )
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
            print(f"  ERROR after {time.perf_counter() - pair_started:.1f} s: {exc}", flush=True)
            if stop_on_error:
                raise
            results.append({"name": pair["name"], "error": str(exc), **pair})
    n_failed = sum("error" in result for result in results)
    print(
        f"\nBatch finished in {time.perf_counter() - batch_started:.1f} s: "
        f"{len(results) - n_failed} successful, {n_failed} failed",
        flush=True,
    )
    return results


__all__: Sequence[str] = (
    "BEHAVIOR_ALIGNED_CSD_DEFAULTS",
    "average_kcsd_columns",
    "average_saved_behavior_kcsd_from_root",
    "ANIMAL_POOLED_FOLDER",
    "ANIMAL_POOLED_DATASETS",
    "discover_behavior_aligned_pairs",
    "epoch_paired_with_onset_baseline",
    "estimate_animal_pooled_kcsd",
    "load_paired_raw",
    "organize_animal_pooled_figures",
    "organize_behavior_aligned_figures",
    "run_behavior_aligned_kcsd",
    "run_behavior_aligned_kcsd_from_root",
    "trial_weighted_psth",
)
