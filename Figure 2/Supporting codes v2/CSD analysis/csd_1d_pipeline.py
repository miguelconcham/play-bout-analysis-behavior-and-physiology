"""
General 1D CSD pipeline for Neuropixels LFP stacks.

Each dataset folder is expected to contain:
  - lfp.npy        (n_stim, n_chan, n_time)
  - x_coords.npy   (n_chan,)
  - y_coords.npy   (n_chan,)
  - T0.npy         scalar sample index of stimulus onset

Typical per-session settings (channel range, reference, bad channel) are
passed via a DatasetConfig dataclass / dict. Use run_dataset() for one
folder or run_many() to process a list of sessions.
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

try:
    from scipy.interpolate import interp1d
    from scipy.ndimage import gaussian_filter1d
except ImportError as exc:
    raise ImportError(
        "SciPy failed to import — usually a NumPy/SciPy version mismatch.\n"
        "If you also use `icsd` (needs numpy<1.23), install the pinned stack:\n"
        "  python -m pip install -U -r requirements.txt\n"
        "  # numpy>=1.21,<1.23  and  scipy>=1.9,<1.10  and  kcsd\n"
        "Do NOT upgrade to scipy 1.10+ while keeping numpy 1.22 / icsd.\n"
        "Then restart the kernel.\n"
        f"Original error: {exc}"
    ) from exc


def _patch_scipy_simps() -> None:
    """kcsd validation imports ``simps``; SciPy≥1.14 renamed it to ``simpson``."""
    import scipy.integrate as _integrate

    if not hasattr(_integrate, "simps") and hasattr(_integrate, "simpson"):
        _integrate.simps = _integrate.simpson  # type: ignore[attr-defined]


def _import_kcsd1d():
    """
    Import KCSD1D using only the ``kcsd`` package (kCSD_env).

    Prefer ``kcsd.KCSD.KCSD1D`` so we never import ``kcsd.validation``
    (that path breaks on SciPy≥1.14 because of ``simps``). Fall back to
    ``from kcsd import KCSD1D`` after patching ``simps``.
    """
    # 1) Direct class import — works even if scipy.simps is missing
    try:
        from kcsd.KCSD import KCSD1D as cls
        return cls
    except Exception:
        pass

    # 2) Package import with simps shim (matches older notebooks)
    try:
        _patch_scipy_simps()
        from kcsd import KCSD1D as cls
        return cls
    except Exception:
        return None


KCSD1D = _import_kcsd1d()


PathLike = Union[str, os.PathLike]


@dataclass
class DatasetConfig:
    """Per-session settings. Only ``folder`` is required."""

    folder: str
    name: Optional[str] = None

    # Timing — these are GENERIC defaults (Lombard/stim style -200/+800 ms).
    # For delta peak/trough runs, run_delta_csd_from_root() overlays
    # DELTA_CSD_DEFAULTS which sets t_pre=t_post=500 (-500…+500 ms).
    fs: float = 2500.0
    t_pre: float = 200.0  # ms (overridden to 500 for delta)
    t_post: float = 800.0  # ms (overridden to 500 for delta)
    step_ms: float = 5.0  # ms between CSD frames

    # Channel selection (applied after bad-channel removal)
    chan_start: Optional[int] = None  # inclusive; None -> 0
    chan_end: Optional[int] = None  # exclusive; None -> all

    # Channels to drop (0-based indices into the *original* probe)
    channels_remove: Tuple[int, ...] = ()

    # Reference: int index, "deepest", "shallowest", or None (no re-ref)
    ref_chan: Union[int, str, None] = "deepest"
    ref_rng_seed: Optional[int] = 0  # used when several deepest channels exist

    # Epoch handling
    baseline_correct: bool = True
    # Reserved for later (delta authenticity CSD currently uses "none")
    # Options: "none" | "zscore_channel" | "zscore_trial"
    normalize: str = "none"

    # 1D kCSD — defaults match the working "CSD to share" notebook
    compute_kcsd: bool = True
    y_scale: str = "unit"  # (y-ymin)/(ymax-ymin), as in the share notebook
    h: float = 0.1
    n_src_factor: float = 2.0  # n_src = factor * n_ch  (share notebook: 2*n_ch)
    smooth_sigma_time: float = 1.0  # gaussian_filter1d over time
    sort_by_depth: bool = True

    # Finite-difference CSD (optional; share notebook uses kCSD only)
    compute_fd: bool = False
    tissue_sigma: float = 0.3  # S/m
    lfp_in_uV: bool = True
    sigma_spatial_ch: float = 1.0

    # I/O
    save: bool = True
    save_plots: bool = True
    output_subdir: str = "csd_1d"
    plot: bool = False

    def resolved_name(self) -> str:
        if self.name:
            return self.name
        return os.path.basename(os.path.normpath(self.folder))


def _as_config(cfg: Union[DatasetConfig, Dict[str, Any]]) -> DatasetConfig:
    if isinstance(cfg, DatasetConfig):
        return cfg
    return DatasetConfig(**cfg)


def load_raw(
    folder: PathLike,
    *,
    mmap_mode: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, Optional[float]]:
    """
    Load a CSD-ready folder.

    Supports both layouts produced by Estimate_delta_peak_trough_csd:
      - NPX1: one shank, full channel count in lfp.npy
      - NPX2: one probe/shank folder with reduced n_chan
    Optional ``fs.npy`` overrides the config sampling rate when present.
    """
    folder = os.fspath(folder)
    lfp = np.load(os.path.join(folder, "lfp.npy"), mmap_mode=mmap_mode)
    x = np.load(os.path.join(folder, "x_coords.npy")).flatten()
    y = np.load(os.path.join(folder, "y_coords.npy")).flatten()
    t0 = int(np.asarray(np.load(os.path.join(folder, "T0.npy"))).reshape(-1)[0])
    fs_file = os.path.join(folder, "fs.npy")
    fs_loaded: Optional[float] = None
    if os.path.isfile(fs_file):
        fs_loaded = float(np.asarray(np.load(fs_file)).reshape(-1)[0])
    if lfp.ndim != 3:
        raise ValueError(f"lfp.npy must be (n_stim, n_chan, n_time); got {lfp.shape}")
    if not (lfp.shape[1] == len(x) == len(y)):
        raise ValueError(
            f"Channel mismatch: lfp n_chan={lfp.shape[1]}, x={len(x)}, y={len(y)}"
        )
    return lfp, x, y, t0, fs_loaded


def probe_geometry_kind(x: np.ndarray) -> str:
    """Return 'NPX1_single_shank', 'NPX2_single_shank', or 'multi_column'."""
    ux = np.unique(x[~np.isnan(x)])
    if ux.size <= 2:
        # NPX1 has two x columns on one shank tip; NPX2 shank also has two sites/row
        return "single_shank"
    return "multi_column"


def _resolve_ref_chan(
    y: np.ndarray,
    ref_chan: Union[int, str, None],
    rng_seed: Optional[int],
) -> Optional[int]:
    if ref_chan is None:
        return None
    if isinstance(ref_chan, (int, np.integer)):
        idx = int(ref_chan)
        if idx < 0 or idx >= len(y):
            raise IndexError(f"ref_chan={idx} out of range for n_chan={len(y)}")
        return idx
    key = str(ref_chan).lower()
    if key == "deepest":
        target = y.min()
    elif key == "shallowest":
        target = y.max()
    else:
        raise ValueError(
            "ref_chan must be int, None, 'deepest', or 'shallowest'; "
            f"got {ref_chan!r}"
        )
    candidates = np.where(y == target)[0]
    if len(candidates) == 1:
        return int(candidates[0])
    rng = np.random.default_rng(rng_seed)
    return int(rng.choice(candidates))


def prepare_channels(
    lfp: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    cfg: DatasetConfig,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    Apply bad-channel removal, referencing, then channel-range selection.

    Order matters: reference is chosen on the probe *after* bad channels are
    removed, then the analysis range is sliced. This avoids the indexing bugs
    in the session-specific notebooks where ``lfp_refed`` was kept full-size
    while ``x``/``y`` were sliced.
    """
    meta: Dict[str, Any] = {"n_chan_raw": int(lfp.shape[1])}

    keep = np.ones(lfp.shape[1], dtype=bool)
    for ch in cfg.channels_remove:
        if 0 <= ch < keep.size:
            keep[ch] = False
    removed = np.where(~keep)[0].tolist()
    # Preserve memory-mapped LFP when no bad channels are removed. Boolean
    # indexing would eagerly copy the complete event stack into RAM.
    if removed:
        lfp = lfp[:, keep, :]
        x = x[keep]
        y = y[keep]
    meta["channels_removed"] = removed

    ref_idx = _resolve_ref_chan(y, cfg.ref_chan, cfg.ref_rng_seed)
    meta["ref_chan_after_removal"] = ref_idx
    if ref_idx is not None:
        lfp = lfp - lfp[:, ref_idx, :][:, None, :]

    # Range is applied after removal (same as the session notebooks).
    # Overshooting end is clamped like NumPy slicing (e.g. end=384 on 383 ch).
    start = 0 if cfg.chan_start is None else int(cfg.chan_start)
    end = lfp.shape[1] if cfg.chan_end is None else int(cfg.chan_end)
    start = max(0, start)
    end = min(lfp.shape[1], end)
    if start >= end:
        raise ValueError(
            f"Invalid channel range [{start}:{end}) for n_chan={lfp.shape[1]} "
            f"(after removing {removed})"
        )
    lfp = lfp[:, start:end, :]
    x = x[start:end]
    y = y[start:end]
    meta["chan_start"] = start
    meta["chan_end"] = end
    meta["n_chan"] = int(lfp.shape[1])
    return lfp, x, y, meta


def apply_normalize(evoked: np.ndarray, mode: str) -> np.ndarray:
    """Optional post-average normalization (hook for later delta CSD work)."""
    mode = (mode or "none").lower()
    if mode in ("none", "", "off"):
        return evoked
    if mode == "zscore_channel":
        mu = evoked.mean(axis=1, keepdims=True)
        sd = evoked.std(axis=1, keepdims=True)
        sd = np.where(sd < 1e-12, 1.0, sd)
        return (evoked - mu) / sd
    raise ValueError(
        f"Unknown normalize={mode!r}; use 'none' or 'zscore_channel' "
        "(more modes can be added later)."
    )


def epoch_and_average(
    lfp: np.ndarray,
    t0: int,
    fs: float,
    t_pre: float,
    t_post: float,
    step_ms: float,
    *,
    baseline_correct: bool = True,
    normalize: str = "none",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns
    -------
    evoked : (n_chan, n_time_win)
    time_ms : (n_time_win,)
    time_indices : indices into the window used for CSD frames
    epochs : (n_stim, n_chan, n_time_win) (baseline-corrected if requested)
    """
    pre_samp = int(round(t_pre * fs / 1000.0))
    post_samp = int(round(t_post * fs / 1000.0))
    start = t0 - pre_samp
    end = t0 + post_samp
    n_time = lfp.shape[2]
    if start < 0 or end > n_time:
        raise RuntimeError(
            f"Window [{start}:{end}) exceeds signal length {n_time} "
            f"(T0={t0}, t_pre={t_pre}, t_post={t_post})"
        )

    epochs = lfp[:, :, start:end]
    if baseline_correct and pre_samp > 0:
        baseline = epochs[:, :, :pre_samp].mean(axis=2, keepdims=True)
        epochs = epochs - baseline
    evoked = epochs.mean(axis=0)
    evoked = apply_normalize(evoked, normalize)
    time_ms = (np.arange(evoked.shape[1]) - pre_samp) / fs * 1000.0
    step = max(1, int(round(step_ms * fs / 1000.0)))
    time_indices = np.arange(0, evoked.shape[1], step)
    return evoked, time_ms, time_indices, epochs


def _normalize_y(y_ch: np.ndarray, y_scale: str) -> np.ndarray:
    y_ch = np.asarray(y_ch, dtype=float)
    if y_scale == "unit":
        span = y_ch.max() - y_ch.min()
        if span <= 0:
            y_norm = np.zeros_like(y_ch)
        else:
            y_norm = (y_ch - y_ch.min()) / span
    elif y_scale == "mm":
        y_norm = (y_ch - y_ch.min()) / 1000.0
    else:
        raise ValueError(f"y_scale must be 'mm' or 'unit'; got {y_scale!r}")
    return y_norm[:, None]


def compute_kcsd_1d_columns(
    evoked: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    time_indices: np.ndarray,
    *,
    y_scale: str = "mm",
    h: float = 0.25,
    n_src_factor: float = 4.0,
    smooth_sigma_time: float = 0.0,
    sort_by_depth: bool = True,
    progress_label: str = "",
) -> Dict[float, np.ndarray]:
    """
    1D kCSD per shank / unique x.

    All selected time frames for one column are solved together. KCSD1D
    accepts ``pots`` as channels × time, so constructing one model per frame
    is unnecessary and is hundreds of times slower for multi-second epochs.
    """
    if KCSD1D is None:
        raise ImportError(
            "kcsd is not installed. Install with `pip install kcsd` "
            "or set compute_kcsd=False and use finite-difference CSD."
        )

    unique_x = np.unique(x)
    csd_dict: Dict[float, np.ndarray] = {}

    for column_number, ux in enumerate(unique_x, start=1):
        idx = np.where(x == ux)[0]
        y_ch = y[idx].astype(float)
        v_ch = evoked[idx, :].astype(float)

        if sort_by_depth:
            order = np.argsort(y_ch)
            y_ch = y_ch[order]
            v_ch = v_ch[order, :]

        n_ch = len(y_ch)
        n_time_sel = len(time_indices)
        label = f"{progress_label}: " if progress_label else ""
        print(
            f"    {label}kCSD column {column_number}/{len(unique_x)} "
            f"(x={float(ux):g}, channels={n_ch}, frames={n_time_sel})",
            flush=True,
        )

        # Same as share notebook: unit-normalize depth, smooth V, KCSD1D, interp back
        y_norm = _normalize_y(y_ch, y_scale)
        if smooth_sigma_time and smooth_sigma_time > 0:
            v_s = gaussian_filter1d(v_ch, sigma=float(smooth_sigma_time), axis=1)
        else:
            v_s = v_ch

        n_src = max(2, int(round(n_src_factor * n_ch)))  # share notebook: 2 * n_ch
        k = KCSD1D(
            ele_pos=y_norm,
            pots=v_s[:, time_indices],
            h=h,
            n_src_init=n_src,
        )
        csd_grid = np.asarray(k.values(estimate="CSD"))
        if csd_grid.ndim == 1:
            csd_grid = csd_grid[:, None]
        grid_y = np.linspace(y_ch.min(), y_ch.max(), csd_grid.shape[0])
        interpolate = interp1d(
            grid_y,
            csd_grid,
            axis=0,
            kind="linear",
            fill_value="extrapolate",
        )
        csd_stack = np.asarray(interpolate(y_ch))
        if csd_stack.shape != (n_ch, n_time_sel):
            raise RuntimeError(
                f"Unexpected vectorized kCSD shape {csd_stack.shape}; "
                f"expected {(n_ch, n_time_sel)}"
            )

        csd_dict[float(ux)] = csd_stack
        print(
            f"    {label}completed column {column_number}/{len(unique_x)}",
            flush=True,
        )

    return csd_dict


def compute_fd_csd_columns(
    evoked: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    *,
    tissue_sigma: float = 0.3,
    lfp_in_uV: bool = True,
    sigma_spatial_ch: float = 1.0,
    sort_by_depth: bool = True,
) -> Dict[float, np.ndarray]:
    """Finite-difference CSD = -σ ∂²V/∂y² per unique x column."""
    unique_x = np.unique(x)
    out: Dict[float, np.ndarray] = {}

    for ux in unique_x:
        idx = np.where(x == ux)[0]
        y_col = y[idx].astype(float)
        v_col = evoked[idx, :].astype(float)
        if sort_by_depth:
            order = np.argsort(y_col)
            y_col = y_col[order]
            v_col = v_col[order, :]

        y_m = y_col * 1e-6
        v = v_col * 1e-6 if lfp_in_uV else v_col.copy()
        if sigma_spatial_ch and sigma_spatial_ch > 0:
            v = gaussian_filter1d(v, sigma=float(sigma_spatial_ch), axis=0, mode="nearest")
        dvy = np.gradient(v, y_m, axis=0)
        d2 = np.gradient(dvy, y_m, axis=0)
        out[float(ux)] = -tissue_sigma * d2

    return out


def stack_columns_for_matlab(
    csd_dict: Dict[float, np.ndarray],
) -> Tuple[np.ndarray, np.ndarray]:
    """Pad columns to common n_y and return (n_y, n_time, n_x), sorted x."""
    all_x = np.array(sorted(csd_dict.keys()), dtype=float)
    n_y = max(v.shape[0] for v in csd_dict.values())
    n_time = next(iter(csd_dict.values())).shape[1]
    arr = np.zeros((n_y, n_time, len(all_x)))
    for i, xv in enumerate(all_x):
        col = csd_dict[xv]
        if col.shape[0] < n_y:
            padded = np.zeros((n_y, n_time))
            padded[-col.shape[0] :, :] = col
            arr[:, :, i] = padded
        else:
            arr[:, :, i] = col
    return arr, all_x


def column_y(x: np.ndarray, y: np.ndarray, ux: float, sort_by_depth: bool = True) -> np.ndarray:
    y_ch = y[np.where(x == ux)[0]].astype(float)
    if sort_by_depth:
        y_ch = np.sort(y_ch)
    return y_ch


def save_kcsd_outputs(
    folder: PathLike,
    csd_dict: Dict[float, np.ndarray],
    x: np.ndarray,
    y: np.ndarray,
    evoked: np.ndarray,
    time_ms: np.ndarray,
    time_indices: np.ndarray,
    *,
    output_subdir: str = "csd_1d",
    smooth_sigma_time: float = 0.0,
    sort_by_depth: bool = True,
) -> str:
    folder = os.fspath(folder)
    save_folder = os.path.join(folder, output_subdir)
    os.makedirs(save_folder, exist_ok=True)

    for ux, csd_stack in csd_dict.items():
        tag = f"x{ux}"
        np.save(os.path.join(save_folder, f"csd_{tag}.npy"), csd_stack)
        y_ch = column_y(x, y, ux, sort_by_depth=sort_by_depth)
        np.save(os.path.join(save_folder, f"y_{tag}.npy"), y_ch)

        idx = np.where(x == ux)[0]
        v_ch = evoked[idx, :].astype(float)
        if sort_by_depth:
            order = np.argsort(y[idx])
            v_ch = v_ch[order, :]
        if smooth_sigma_time and smooth_sigma_time > 0:
            v_ch = gaussian_filter1d(v_ch, sigma=float(smooth_sigma_time), axis=1)
        np.save(os.path.join(save_folder, f"V_{tag}.npy"), v_ch)
        np.save(os.path.join(save_folder, f"time_voltage_{tag}.npy"), time_ms)

    np.save(os.path.join(save_folder, "time.npy"), time_ms[time_indices])

    csd_array, all_x = stack_columns_for_matlab(csd_dict)
    np.save(os.path.join(folder, "csd_1d_array.npy"), csd_array)
    np.save(os.path.join(folder, "csd_1d_x.npy"), all_x)
    np.save(os.path.join(folder, "csd_1d_time.npy"), time_ms[time_indices])
    return save_folder


def save_fd_outputs(
    folder: PathLike,
    csd_fd_dict: Dict[float, np.ndarray],
    *,
    output_subdir: str = "csd_1d",
) -> str:
    folder = os.fspath(folder)
    save_folder = os.path.join(folder, output_subdir)
    os.makedirs(save_folder, exist_ok=True)
    for ux, csd in csd_fd_dict.items():
        np.save(os.path.join(save_folder, f"csd_fd_x{ux}.npy"), csd)
    arr, all_x = stack_columns_for_matlab(csd_fd_dict)
    np.save(os.path.join(folder, "csd_fd_array.npy"), arr)
    np.save(os.path.join(folder, "csd_fd_x.npy"), all_x)
    return save_folder


def plot_column(
    csd_stack: np.ndarray,
    time_ms: np.ndarray,
    time_indices: np.ndarray,
    y_ch: np.ndarray,
    *,
    title: str = "",
    cbar_label: str = "CSD",
    clim: Optional[Tuple[float, float]] = None,
    horizontal_lines: Optional[Sequence[float]] = None,
    save_path: Optional[PathLike] = None,
    show: bool = True,
):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(
        csd_stack,
        origin="lower",
        aspect="auto",
        extent=[time_ms[time_indices[0]], time_ms[time_indices[-1]], y_ch.min(), y_ch.max()],
        cmap="RdBu_r",
    )
    if clim is not None:
        im.set_clim(*clim)
    if horizontal_lines:
        for depth in horizontal_lines:
            ax.axhline(
                float(depth),
                color="white",
                linestyle="--",
                linewidth=1.0,
                zorder=3,
            )
    fig.colorbar(im, ax=ax, label=cbar_label)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Depth (µm)")
    ax.set_title(title)
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def run_dataset(cfg: Union[DatasetConfig, Dict[str, Any]]) -> Dict[str, Any]:
    """Full load → preprocess → CSD → optional save/plot for one folder."""
    cfg = _as_config(cfg)
    name = cfg.resolved_name()
    print(f"\n=== {name} ===")
    print(f"folder: {cfg.folder}")

    lfp, x, y, t0, fs_loaded = load_raw(cfg.folder)
    if fs_loaded is not None:
        cfg.fs = fs_loaded
    lfp, x, y, meta = prepare_channels(lfp, x, y, cfg)
    meta["geometry"] = probe_geometry_kind(x)
    meta["n_events"] = int(lfp.shape[0])
    meta["baseline_correct"] = bool(cfg.baseline_correct)
    meta["normalize"] = cfg.normalize
    print(
        f"channels: {meta['n_chan']} ({meta['geometry']}, n_events={meta['n_events']}) "
        f"(removed {meta['channels_removed']}, "
        f"range [{meta['chan_start']}:{meta['chan_end']}), "
        f"ref={meta['ref_chan_after_removal']}, "
        f"baseline={cfg.baseline_correct}, normalize={cfg.normalize})"
    )
    print(
        f"time window: {-cfg.t_pre:.0f} … +{cfg.t_post:.0f} ms around T0 "
        f"(t_pre={cfg.t_pre}, t_post={cfg.t_post})"
    )

    evoked, time_ms, time_indices, epochs_bc = epoch_and_average(
        lfp,
        t0,
        cfg.fs,
        cfg.t_pre,
        cfg.t_post,
        cfg.step_ms,
        baseline_correct=cfg.baseline_correct,
        normalize=cfg.normalize,
    )
    # Spacing within each shank (unique x); global sort is misleading for NPX.
    spacings = []
    for ux in np.unique(x):
        ys = np.sort(y[x == ux])
        if len(ys) > 1:
            spacings.append(np.median(np.diff(ys)))
    median_spacing = float(np.median(spacings)) if spacings else float("nan")
    print(f"evoked {evoked.shape}, median within-shank spacing {median_spacing:.1f} µm")

    result: Dict[str, Any] = {
        "name": name,
        "config": asdict(cfg),
        "meta": meta,
        "x": x,
        "y": y,
        "evoked": evoked,
        "time_ms": time_ms,
        "time_indices": time_indices,
        "epochs_bc": epochs_bc,
        "median_spacing_um": median_spacing,
        "csd_kcsd": None,
        "csd_fd": None,
    }

    if cfg.compute_kcsd:
        print("Computing 1D kCSD …")
        csd_dict = compute_kcsd_1d_columns(
            evoked,
            x,
            y,
            time_indices,
            y_scale=cfg.y_scale,
            h=cfg.h,
            n_src_factor=cfg.n_src_factor,
            smooth_sigma_time=cfg.smooth_sigma_time,
            sort_by_depth=cfg.sort_by_depth,
            progress_label=name,
        )
        result["csd_kcsd"] = csd_dict
        if cfg.save:
            save_folder = save_kcsd_outputs(
                cfg.folder,
                csd_dict,
                x,
                y,
                evoked,
                time_ms,
                time_indices,
                output_subdir=cfg.output_subdir,
                smooth_sigma_time=cfg.smooth_sigma_time,
                sort_by_depth=cfg.sort_by_depth,
            )
            print(f"Saved kCSD matrices → {save_folder}")
            print(
                f"  (per-column csd_x*.npy, stacked csd_1d_array.npy "
                f"shape n_y × n_time × n_x, plus csd_1d_time.npy)"
            )
        if cfg.plot or cfg.save_plots:
            fig_dir = os.path.join(cfg.folder, cfg.output_subdir, "figures")
            if cfg.save_plots:
                os.makedirs(fig_dir, exist_ok=True)
            # Save/show one figure per shank column (unique x)
            for ux in sorted(csd_dict.keys()):
                ux = float(ux)
                y_ch = column_y(x, y, ux, sort_by_depth=cfg.sort_by_depth)
                plot_path = None
                if cfg.save_plots:
                    plot_path = os.path.join(fig_dir, f"kcsd_x{ux}.png")
                plot_column(
                    csd_dict[ux],
                    time_ms,
                    time_indices,
                    y_ch,
                    title=f"{name} — 1D kCSD x={ux}",
                    cbar_label="kCSD (a.u.)",
                    save_path=plot_path,
                    # only show the first column interactively to avoid many popups
                    show=cfg.plot and ux == float(sorted(csd_dict.keys())[0]),
                )
            if cfg.save_plots:
                print(f"Saved kCSD figures → {fig_dir}")

    if cfg.compute_fd:
        print("Computing finite-difference CSD …")
        csd_fd = compute_fd_csd_columns(
            evoked,
            x,
            y,
            tissue_sigma=cfg.tissue_sigma,
            lfp_in_uV=cfg.lfp_in_uV,
            sigma_spatial_ch=cfg.sigma_spatial_ch,
            sort_by_depth=cfg.sort_by_depth,
        )
        result["csd_fd"] = csd_fd
        if cfg.save:
            save_fd_outputs(cfg.folder, csd_fd, output_subdir=cfg.output_subdir)
            print(f"Saved FD CSD → {os.path.join(cfg.folder, cfg.output_subdir)}")
        if cfg.plot or cfg.save_plots:
            fig_dir = os.path.join(cfg.folder, cfg.output_subdir, "figures")
            if cfg.save_plots:
                os.makedirs(fig_dir, exist_ok=True)
            for ux in sorted(csd_fd.keys()):
                ux = float(ux)
                y_ch = column_y(x, y, ux, sort_by_depth=cfg.sort_by_depth)
                fd_indices = np.arange(csd_fd[ux].shape[1])
                plot_path = None
                if cfg.save_plots:
                    plot_path = os.path.join(fig_dir, f"fd_x{ux}.png")
                plot_column(
                    csd_fd[ux],
                    time_ms,
                    fd_indices,
                    y_ch,
                    title=f"{name} — FD CSD x={ux}",
                    cbar_label="CSD (A/m^3)",
                    save_path=plot_path,
                    show=cfg.plot and ux == float(sorted(csd_fd.keys())[0]),
                )
            if cfg.save_plots:
                print(f"Saved FD figures → {fig_dir}")

    print(f"Done: {name}")
    return result


def merge_config(
    defaults: Union[DatasetConfig, Dict[str, Any]],
    overrides: Union[DatasetConfig, Dict[str, Any]],
) -> DatasetConfig:
    """Merge ``overrides`` on top of ``defaults`` (dict keys only overwrite)."""
    if isinstance(overrides, DatasetConfig):
        return overrides
    if isinstance(defaults, DatasetConfig):
        base = asdict(defaults)
    else:
        base = dict(defaults)
    base.update(overrides)
    if "folder" not in base or not base["folder"]:
        raise ValueError("merged config needs a 'folder'")
    return DatasetConfig(**base)


def run_many(
    configs: Sequence[Union[DatasetConfig, Dict[str, Any]]],
    *,
    defaults: Optional[Union[DatasetConfig, Dict[str, Any]]] = None,
    stop_on_error: bool = False,
) -> List[Dict[str, Any]]:
    """
    Process many sessions.

    If ``defaults`` is given, each *dict* config is merged on top of those
    defaults (handy when most sessions share fs / window / kCSD params).
    """
    results: List[Dict[str, Any]] = []
    for raw in configs:
        if defaults is not None and isinstance(raw, dict):
            cfg = merge_config(defaults, raw)
        else:
            cfg = _as_config(raw)
        try:
            results.append(run_dataset(cfg))
        except Exception as exc:  # noqa: BLE001 - batch runner should continue
            print(f"ERROR on {cfg.resolved_name()}: {exc}")
            if stop_on_error:
                raise
            results.append(
                {"name": cfg.resolved_name(), "error": str(exc), "config": asdict(cfg)}
            )
    return results


# ---------------------------------------------------------------------------
# Delta peak/trough CSD folders (from Estimate_delta_peak_trough_csd.m)
# ---------------------------------------------------------------------------

# kCSD settings copied from the working "CSD to share to cursor.ipynb"
SHARE_NOTEBOOK_KCSD: Dict[str, Any] = {
    "fs": 2500.0,
    "step_ms": 5.0,
    "compute_kcsd": True,
    "compute_fd": False,       # share notebook only runs 1D kCSD
    "y_scale": "unit",
    "h": 0.1,
    "n_src_factor": 2.0,
    "smooth_sigma_time": 1.0,
    "sort_by_depth": True,
    "save": True,
    "save_plots": True,
    "plot": True,
}

# Delta peak/trough exports: same kCSD maths as the share notebook.
# Time window is -500 … +500 ms around T0 (matches MATLAB Estimate cuts),
# with baseline correction from the pre-T0 half of that window.
DELTA_CSD_DEFAULTS: Dict[str, Any] = {
    **SHARE_NOTEBOOK_KCSD,
    "t_pre": 500.0,   # ms before T0  →  -500 ms   << look here for delta window
    "t_post": 500.0,  # ms after T0   →  +500 ms
    "ref_chan": None,  # still no re-reference for authenticity analysis
    "baseline_correct": True,
    "normalize": "none",
    "channels_remove": (),
    "chan_start": None,
    "chan_end": None,
    "save": True,
    "save_plots": True,  # PNG under <session>/csd_1d/figures/
    "plot": False,
}


def discover_delta_csd_datasets(
    root: PathLike,
    *,
    alignments: Sequence[str] = ("peak", "trough"),
) -> List[Dict[str, Any]]:
    """
    Walk the MATLAB export tree and build DATASETS entries.

    Expected layout (top-level peak/ and trough/ folders)::

        <root>/peak/<animal>/                 # NPX1
        <root>/peak/<animal>/probe2/          # NPX2 shank
        <root>/trough/<animal>/
        <root>/trough/<animal>/probe3/

    Also accepts the older leaf layout ``.../<animal>/peak/`` if present.
    """
    root = os.fspath(root)
    found: List[Dict[str, Any]] = []
    if not os.path.isdir(root):
        return found

    required = ("lfp.npy", "x_coords.npy", "y_coords.npy", "T0.npy")
    for dirpath, _dirnames, filenames in os.walk(root):
        if not all(f in filenames for f in required):
            continue
        rel = os.path.relpath(dirpath, root).replace("\\", "/")
        parts = rel.split("/")
        # Pooled leaves live under animal_pooled/ and are not session exports.
        if "animal_pooled" in parts:
            continue
        # New layout: peak|trough is the first path component
        # Old layout: peak|trough is the last path component
        if not parts:
            continue
        if parts[0] not in alignments and parts[-1] not in alignments:
            continue
        name = rel.replace("/", "_")
        found.append({"name": name, "folder": dirpath})
    found.sort(key=lambda d: d["name"])
    return found


def run_delta_csd_from_root(
    root: PathLike,
    *,
    stop_on_error: bool = False,
    plot: Optional[bool] = None,
    **overrides: Any,
) -> List[Dict[str, Any]]:
    """
    One-call entry: pass the Estimate output folder and run CSD on every
    peak/ and trough/ leaf dataset found underneath.

    Example
    -------
    run_delta_csd_from_root(
        r"\\\\experimentfs.home.bccn-berlin.de\\\\...\\\\delta peak trough CSD"
    )
    """
    root = os.fspath(root)
    if not os.path.isdir(root):
        raise FileNotFoundError(f"Delta CSD root not found: {root}")

    datasets = discover_delta_csd_datasets(root)
    if not datasets:
        raise FileNotFoundError(
            f"No peak/trough folders with lfp.npy under {root}. "
            "Expected <root>/peak/<session>/ and <root>/trough/<session>/."
        )

    defaults = dict(DELTA_CSD_DEFAULTS)
    defaults.update(overrides)
    if plot is not None:
        defaults["plot"] = plot

    print(f"Delta CSD root: {root}")
    print(
        f"Effective window: -{defaults['t_pre']:.0f} … +{defaults['t_post']:.0f} ms | "
        f"baseline_correct={defaults['baseline_correct']} | "
        f"ref_chan={defaults.get('ref_chan')} | "
        f"save_plots={defaults.get('save_plots')}"
    )
    print(
        "(Note: DatasetConfig class defaults t_pre=200/t_post=800 are NOT used here; "
        "DELTA_CSD_DEFAULTS overrides them.)"
    )
    print(f"Discovered {len(datasets)} dataset folder(s)")
    for d in datasets:
        print(f"  - {d['name']}")

    results = run_many(datasets, defaults=defaults, stop_on_error=stop_on_error)

    # Also collect a flat copy of all PNGs under <root>/csd_figures/
    if defaults.get("save_plots", True):
        fig_root = os.path.join(root, "csd_figures")
        os.makedirs(fig_root, exist_ok=True)
        n_copied = 0
        for r in results:
            if "error" in r:
                continue
            src_dir = os.path.join(r["config"]["folder"], r["config"]["output_subdir"], "figures")
            if not os.path.isdir(src_dir):
                continue
            dest_dir = os.path.join(fig_root, r["name"])
            os.makedirs(dest_dir, exist_ok=True)
            for fn in os.listdir(src_dir):
                if not fn.lower().endswith(".png"):
                    continue
                src = os.path.join(src_dir, fn)
                dst = os.path.join(dest_dir, fn)
                try:
                    import shutil

                    shutil.copy2(src, dst)
                    n_copied += 1
                except OSError as exc:
                    print(f"Could not copy {src} → {dst}: {exc}")
        print(f"Copied {n_copied} figure(s) → {fig_root}")

    return results


# ---------------------------------------------------------------------------
# Example session table matching the uploaded notebooks
# ---------------------------------------------------------------------------

EXAMPLE_DATASETS: List[Dict[str, Any]] = [
    {
        "name": "230516_USV",
        "folder": (
            "//experimentfs.home.bccn-berlin.de/experiment/share/Miguel/"
            "LOMBARD EFFECT NPX files and analysis/CSD analysis/230516/USV CSD"
        ),
        "chan_start": 0,
        "chan_end": 383,
        "channels_remove": (191,),
        "ref_chan": 0,
        "y_scale": "unit",
        "h": 0.1,
        "n_src_factor": 2.0,
        "smooth_sigma_time": 1.0,
    },
    {
        "name": "230516",
        "folder": (
            "//experimentfs.home.bccn-berlin.de/experiment/share/Miguel/"
            "LOMBARD EFFECT NPX files and analysis/CSD analysis/230516"
        ),
        "chan_start": 100,
        "chan_end": 260,
        "channels_remove": (191,),
        "ref_chan": 0,
        "y_scale": "unit",
        "h": 0.1,
        "n_src_factor": 2.0,
        "smooth_sigma_time": 1.0,
    },
    {
        "name": "230606_TONE",
        "folder": (
            "//experimentfs.home.bccn-berlin.de/experiment/share/Miguel/"
            "LOMBARD EFFECT NPX files and analysis/CSD analysis/230606/TONE responses"
        ),
        "chan_start": 75,
        "chan_end": 275,
        "channels_remove": (191,),
        "ref_chan": 0,
        "y_scale": "unit",
        "h": 0.1,
        "n_src_factor": 2.0,
        "smooth_sigma_time": 1.0,
    },
    {
        "name": "230608_USV",
        "folder": (
            "//experimentfs.home.bccn-berlin.de/experiment/share/Miguel/"
            "LOMBARD EFFECT NPX files and analysis/CSD analysis/230608/USV sitm"
        ),
        "chan_start": 100,
        "chan_end": 384,
        "channels_remove": (191,),
        "ref_chan": 0,
        "y_scale": "unit",
        "h": 0.1,
        "n_src_factor": 2.0,
        "smooth_sigma_time": 1.0,
    },
    {
        "name": "250320_USV",
        "folder": (
            "//experimentfs.home.bccn-berlin.de/experiment/share/Miguel/"
            "LOMBARD EFFECT NPX files and analysis/CSD analysis/250320/USV stim"
        ),
        "chan_start": 0,
        "chan_end": 384,
        "channels_remove": (),
        "ref_chan": "deepest",
        "y_scale": "mm",
        "h": 0.25,
        "n_src_factor": 4.0,
        "smooth_sigma_time": 0.0,
    },
]


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI: python csd_1d_pipeline.py <delta-root>   OR   --delta-root PATH"""
    import argparse

    p = argparse.ArgumentParser(
        description=(
            "Run 1D CSD. Preferred: pass the Estimate output folder "
            "(.../delta peak trough CSD) as a positional path or --delta-root."
        )
    )
    p.add_argument(
        "delta_root_positional",
        nargs="?",
        default=None,
        help="Estimate output folder containing peak/ and trough/",
    )
    p.add_argument("--folder", action="append", default=[], help="Single dataset folder (repeatable)")
    p.add_argument(
        "--delta-root",
        default=None,
        help="Same as positional: root from Estimate_delta_peak_trough_csd",
    )
    p.add_argument("--chan-start", type=int, default=None)
    p.add_argument("--chan-end", type=int, default=None)
    p.add_argument("--remove", type=int, action="append", default=[], help="Bad channel index")
    p.add_argument("--ref", default="deepest", help="int | deepest | shallowest | none")
    p.add_argument("--no-baseline", action="store_true", help="Disable pre-event baseline")
    p.add_argument("--normalize", default="none", help="none | zscore_channel")
    p.add_argument("--no-kcsd", action="store_true")
    p.add_argument("--no-fd", action="store_true")
    p.add_argument("--no-save", action="store_true")
    p.add_argument("--plot", action="store_true")
    p.add_argument("--examples", action="store_true", help="Run EXAMPLE_DATASETS table")
    args = p.parse_args(argv)

    if args.examples:
        run_many(EXAMPLE_DATASETS, stop_on_error=False)
        return 0

    delta_root = args.delta_root or args.delta_root_positional
    if delta_root:
        overrides: Dict[str, Any] = {}
        if args.no_kcsd:
            overrides["compute_kcsd"] = False
        if args.no_fd:
            overrides["compute_fd"] = False
        if args.no_save:
            overrides["save"] = False
        if args.normalize:
            overrides["normalize"] = args.normalize
        run_delta_csd_from_root(
            delta_root,
            stop_on_error=False,
            plot=args.plot if args.plot else None,
            **overrides,
        )
        return 0

    if not args.folder:
        p.error(
            "Provide the delta CSD root folder (positional or --delta-root), "
            "or --folder / --examples"
        )

    ref: Union[int, str, None]
    if args.ref.lower() == "none":
        ref = None
    else:
        try:
            ref = int(args.ref)
        except ValueError:
            ref = args.ref

    configs = [
        DatasetConfig(
            folder=f,
            chan_start=args.chan_start,
            chan_end=args.chan_end,
            channels_remove=tuple(args.remove),
            ref_chan=ref,
            baseline_correct=not args.no_baseline,
            normalize=args.normalize,
            compute_kcsd=not args.no_kcsd,
            compute_fd=not args.no_fd,
            save=not args.no_save,
            plot=args.plot,
            save_plots=True,
        )
        for f in args.folder
    ]
    run_many(configs, stop_on_error=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
