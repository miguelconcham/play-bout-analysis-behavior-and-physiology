"""
Cross-session anatomical consistency of saved behavior-aligned kCSD.

This module does not estimate kCSD. It reads ``csd_mean_columns.npy`` maps
that were already averaged across staggered probe columns, assigns those
depths to anatomical areas with the same rules as
``GENERATE_PSTH_AREA_MAPS.m``, and asks whether the saved kCSD sign repeats
across sessions.

Polarity is reported exactly as stored (``RdBu_r`` in the existing figures).
Positive and negative are not relabeled as source or sink.

Each probe map is also smoothed with a Gaussian of 250 ms by 100 µm FWHM
before a one-sample t-test across sessions. Clusters that pass a pointwise
threshold are kept only when they span at least 25 ms and 10% of the
relative-depth axis. A probe segment whose unique depths span less than
200 µm is left out of that session.

Each accepted session map is written to ``area_consistency/session_cache``
before the group figures are saved. ``finalize_area_consistency_from_cache``
rebuilds the means, clusters, and copies from that cache. A refresh does not
delete the cache, and a failed delete on experimentFS is overwritten instead
of stopping the run.

Area rules, matching the MATLAB generator:

- Animal key: session folder ``<batch> <date> <repeated> ...``.
  ``Single2`` -> ``Batch`` + second batch character + ``Single2``.
  Otherwise -> ``Batch`` + second batch character + repeated animal
  + fourth batch character. MATLAB indexes those characters from 1.
- A depth ``y`` belongs to a row when
  ``y >= depth_start`` and ``y < depth_end + 1``.
- NPX1 also requires ``Probe_Area == PAG``. NPX2 also requires
  ``ProbeNum`` equal to the ``probeN`` folder. NPX2 does not filter
  ``Probe_Area``.
- If several rows match one depth, the narrowest depth span wins. If that
  is still a tie, the depth is left unassigned.
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
import shutil
import stat
import warnings
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np


PathLike = Union[str, os.PathLike]

# FWHM = 2 * sqrt(2 * ln 2) * sigma.
_FWHM_PER_SIGMA = 2.0 * float(np.sqrt(2.0 * np.log(2.0)))
DEFAULT_SMOOTH_FWHM_MS = 250.0
DEFAULT_SMOOTH_FWHM_UM = 100.0
DEFAULT_MIN_DEPTH_SPAN_UM = 200.0
DEFAULT_MIN_SESSIONS_FOR_TTEST = 5
DEFAULT_TTEST_ALPHA = 0.05
DEFAULT_MIN_CLUSTER_DURATION_MS = 25.0
DEFAULT_MIN_CLUSTER_DEPTH_FRACTION = 0.10
DEFAULT_DEPTH_GRID_STEP_UM = 20.0
SESSION_CACHE_NAME = "session_cache"
_MIN_KERNEL_SUPPORT = 0.5
_FOUR_CONNECTED = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=int)

_REQUIRED_COLUMNS = (
    "animalname",
    "area",
    "probenum",
    "depth_start",
    "depth_end",
    "probe_area",
)


@dataclass(frozen=True)
class AreaLimit:
    """One row of ``Area_limits_GoodLooking.xlsx``."""

    animal_name: str
    area: str
    probe_num: Optional[int]
    depth_start: float
    depth_end: float
    probe_area: str


def animal_key_from_session(session_name: str) -> str:
    """
    Map a session folder name to ``AnimalName`` in the area table.

    Examples
    --------
    ``B1A2 2024-01-01 RatA notes`` -> ``Batch1RatA2``
    ``B1XX 2024-01-01 Single2`` -> ``Batch1Single2``
    """
    tokens = str(session_name).strip().split()
    if len(tokens) < 3:
        raise ValueError(
            "Cannot parse implanted animal from session folder "
            f"{session_name!r}; expected '<batch> <date> <repeated animal> ...'"
        )
    animal_batch = tokens[0]
    repeated_animal = tokens[2]
    if len(animal_batch) < 2:
        raise ValueError(
            f"Batch token {animal_batch!r} has no second character for the "
            "MATLAB AnimalName key"
        )
    if repeated_animal == "Single2":
        return "Batch" + animal_batch[1] + repeated_animal
    if len(animal_batch) < 4:
        raise ValueError(
            f"Batch token {animal_batch!r} has no fourth character for the "
            "MATLAB AnimalName key"
        )
    return "Batch" + animal_batch[1] + repeated_animal + animal_batch[3]


def load_area_limits(path: PathLike) -> List[AreaLimit]:
    """
    Read area rows from the workbook used by the MATLAB generators.

    ``Area_limits_GoodLooking.xlsx`` is source data. This function opens it
    with a read-only byte copy and never saves, overwrites, or truncates it.
    """
    path = os.fspath(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Area table not found: {path}. Expected "
            "Data/Area_limits_GoodLooking.xlsx on experimentFS."
        )
    try:
        import openpyxl
    except ImportError as exc:
        raise ImportError(
            "Reading Area_limits_GoodLooking.xlsx requires openpyxl. "
            "Install it in the kCSD_env kernel with: python -m pip install openpyxl"
        ) from exc

    with open(path, "rb") as stream:
        payload = io.BytesIO(stream.read())
    workbook = openpyxl.load_workbook(payload, read_only=True, data_only=True)
    try:
        sheet = workbook.active
        iterator = sheet.iter_rows(values_only=True)
        header = next(iterator, None)
        if not header:
            raise ValueError(f"Area table {path} has no header row")
        columns = {
            str(name).strip().lower(): index
            for index, name in enumerate(header)
            if name is not None and str(name).strip()
        }
        missing = [name for name in _REQUIRED_COLUMNS if name not in columns]
        if missing:
            raise ValueError(
                f"Area table {path} is missing columns {missing}. "
                f"Found {sorted(columns)}."
            )
        rows: List[AreaLimit] = []
        for values in iterator:
            if values is None:
                continue
            animal = _cell_text(values, columns["animalname"])
            area = _cell_text(values, columns["area"])
            if not animal or not area:
                continue
            try:
                depth_start = float(values[columns["depth_start"]])
                depth_end = float(values[columns["depth_end"]])
            except (TypeError, ValueError):
                continue
            if not np.isfinite(depth_start) or not np.isfinite(depth_end):
                continue
            rows.append(
                AreaLimit(
                    animal_name=animal,
                    area=area,
                    probe_num=_as_probe_num(values[columns["probenum"]]),
                    depth_start=depth_start,
                    depth_end=depth_end,
                    probe_area=_cell_text(values, columns["probe_area"]),
                )
            )
    finally:
        workbook.close()
    if not rows:
        raise ValueError(f"Area table {path} has no usable rows")
    return rows


def select_area_for_depth(
    y: float,
    rows: Sequence[AreaLimit],
    *,
    npx_type: int,
    probe_num: int,
) -> Optional[AreaLimit]:
    """
    Return the single area row that owns depth ``y``.

    ``rows`` must already be limited to one animal. NPX1 keeps ``PAG``
    probe areas. NPX2 keeps the matching probe number.
    """
    matches: List[AreaLimit] = []
    for row in rows:
        if row.depth_end < row.depth_start:
            continue
        if not (y >= row.depth_start and y < row.depth_end + 1):
            continue
        if npx_type == 1:
            if row.probe_area != "PAG":
                continue
        elif npx_type == 2:
            if row.probe_num != probe_num:
                continue
        else:
            raise ValueError(f"Unsupported NPX type: {npx_type}")
        matches.append(row)
    if not matches:
        return None
    narrowest = min(row.depth_end - row.depth_start for row in matches)
    tied = [
        row
        for row in matches
        if np.isclose(row.depth_end - row.depth_start, narrowest)
    ]
    identities = {(row.area, row.depth_start, row.depth_end) for row in tied}
    if len(identities) != 1:
        return None
    return tied[0]


def resample_area_segment(
    y: np.ndarray,
    csd: np.ndarray,
    depth_start: float,
    depth_end: float,
    relative_depth: np.ndarray,
) -> Optional[np.ndarray]:
    """
    Resample one area segment onto relative depth without extrapolation.

    Relative depth is ``(y - depth_start) / (depth_end - depth_start)``.
    Query positions outside the observed relative range are NaN.
    """
    y = np.asarray(y, dtype=float).reshape(-1)
    csd = np.asarray(csd, dtype=float)
    if csd.ndim != 2 or csd.shape[0] != len(y):
        raise ValueError(
            f"CSD shape {csd.shape} does not match {len(y)} depths"
        )
    span = float(depth_end) - float(depth_start)
    if not np.isfinite(span) or span <= 0:
        return None
    relative = (y - float(depth_start)) / span
    resampled = np.full((len(relative_depth), csd.shape[1]), np.nan, dtype=float)
    any_finite = False
    for time_index in range(csd.shape[1]):
        good = np.isfinite(relative) & np.isfinite(csd[:, time_index])
        if int(np.count_nonzero(good)) < 2:
            continue
        column = _interpolate_rows(
            relative[good],
            csd[good, time_index : time_index + 1],
            relative_depth,
        )
        if column is None:
            continue
        resampled[:, time_index] = column[:, 0]
        any_finite = True
    if not any_finite:
        return None
    return resampled


def smooth_csd_in_depth_and_time(
    y: np.ndarray,
    csd: np.ndarray,
    time_ms: np.ndarray,
    *,
    fwhm_um: float = DEFAULT_SMOOTH_FWHM_UM,
    fwhm_ms: float = DEFAULT_SMOOTH_FWHM_MS,
    grid_step_um: float = DEFAULT_DEPTH_GRID_STEP_UM,
    min_support: float = _MIN_KERNEL_SUPPORT,
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """
    Gaussian-smooth one probe segment on a uniform micrometer grid.

    Widths are full widths at half maximum. Samples are removed where half
    or more of the kernel lies outside the observed depths or time range.
    """
    from scipy import ndimage

    y = np.asarray(y, dtype=float).reshape(-1)
    csd = np.asarray(csd, dtype=float)
    time_ms = np.asarray(time_ms, dtype=float).reshape(-1)
    if csd.shape != (len(y), len(time_ms)):
        raise ValueError(
            f"CSD shape {csd.shape} does not match y={len(y)} and time={len(time_ms)}"
        )
    order = np.argsort(y)
    y = y[order]
    csd = csd[order]
    unique_y, inverse = np.unique(y, return_inverse=True)
    if len(unique_y) < len(y):
        totals = np.zeros((len(unique_y), csd.shape[1]), dtype=float)
        counts = np.zeros(len(unique_y), dtype=float)
        for index, group in enumerate(inverse):
            totals[group] += csd[index]
            counts[group] += 1.0
        csd = totals / counts[:, None]
        y = unique_y
    if len(y) < 2 or y[-1] <= y[0]:
        return None
    step = float(grid_step_um)
    if step <= 0:
        raise ValueError("grid_step_um must be positive")
    n_grid = max(2, int(round((y[-1] - y[0]) / step)) + 1)
    grid = np.linspace(y[0], y[-1], n_grid)
    on_grid = _interpolate_rows(y, csd, grid)
    if on_grid is None:
        return None
    dy = float(grid[1] - grid[0])
    dt = _axis_step(time_ms)
    sigma_y = 0.0 if fwhm_um <= 0 or dy <= 0 else (fwhm_um / _FWHM_PER_SIGMA) / dy
    sigma_t = 0.0 if fwhm_ms <= 0 or dt <= 0 else (fwhm_ms / _FWHM_PER_SIGMA) / dt
    smoothed = _nan_aware_gaussian(
        on_grid,
        (sigma_y, sigma_t),
        min_support=min_support,
        ndimage=ndimage,
    )
    return grid, smoothed


def infer_saved_kcsd_clusters(
    stack: np.ndarray,
    time_ms: np.ndarray,
    relative_depth: np.ndarray,
    *,
    alpha: float = DEFAULT_TTEST_ALPHA,
    min_sessions: int = DEFAULT_MIN_SESSIONS_FOR_TTEST,
    min_duration_ms: float = DEFAULT_MIN_CLUSTER_DURATION_MS,
    min_depth_fraction: float = DEFAULT_MIN_CLUSTER_DEPTH_FRACTION,
) -> Dict[str, object]:
    """One-sample t-test of smoothed session maps, then a cluster-size filter."""
    t_stat, p_value, n_sessions, mean = one_sample_ttest_map(stack)
    mask, clusters = label_significant_clusters(
        mean,
        p_value,
        n_sessions,
        time_ms,
        relative_depth,
        t_stat=t_stat,
        alpha=alpha,
        min_sessions=min_sessions,
        min_duration_ms=min_duration_ms,
        min_depth_fraction=min_depth_fraction,
    )
    return {
        "mean": mean,
        "t": t_stat,
        "p": p_value,
        "n_sessions": n_sessions,
        "mask": mask,
        "clusters": clusters,
    }


def one_sample_ttest_map(
    stack: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """T-test each depth/time bin of session maps against zero."""
    from scipy import stats

    stack = np.asarray(stack, dtype=float)
    finite = np.isfinite(stack)
    n_sessions = np.sum(finite, axis=0).astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean = np.nanmean(stack, axis=0)
        std = np.nanstd(stack, axis=0, ddof=1)
    standard_error = np.full(mean.shape, np.nan, dtype=float)
    comparable = n_sessions >= 2
    standard_error[comparable] = std[comparable] / np.sqrt(n_sessions[comparable])
    t_stat = np.full(mean.shape, np.nan, dtype=float)
    finite_error = comparable & np.isfinite(standard_error) & (standard_error > 0)
    t_stat[finite_error] = mean[finite_error] / standard_error[finite_error]
    perfect = comparable & np.isfinite(standard_error) & (standard_error == 0) & (mean != 0)
    t_stat[perfect] = np.sign(mean[perfect]) * np.inf
    p_value = np.full(mean.shape, np.nan, dtype=float)
    df = n_sessions - 1
    testable = np.isfinite(t_stat) & (df >= 1)
    if np.any(testable):
        p_value[testable] = 2.0 * stats.t.sf(np.abs(t_stat[testable]), df[testable])
    p_value[perfect] = 0.0
    return t_stat, p_value, n_sessions, mean


def label_significant_clusters(
    mean: np.ndarray,
    p_value: np.ndarray,
    n_sessions: np.ndarray,
    time_ms: np.ndarray,
    relative_depth: np.ndarray,
    *,
    t_stat: Optional[np.ndarray] = None,
    alpha: float = DEFAULT_TTEST_ALPHA,
    min_sessions: int = DEFAULT_MIN_SESSIONS_FOR_TTEST,
    min_duration_ms: float = DEFAULT_MIN_CLUSTER_DURATION_MS,
    min_depth_fraction: float = DEFAULT_MIN_CLUSTER_DEPTH_FRACTION,
) -> Tuple[np.ndarray, List[Dict[str, object]]]:
    """
    Keep 4-connected same-sign clusters whose bounding box is large enough.

    Positive and negative bins are clustered separately. The duration and
    depth tests use the inclusive bin width, so five 5 ms frames count as
    25 ms and four relative-depth bins of 0.025 count as 0.10.
    """
    from scipy import ndimage

    mean = np.asarray(mean, dtype=float)
    p_value = np.asarray(p_value, dtype=float)
    n_sessions = np.asarray(n_sessions, dtype=float)
    time_ms = np.asarray(time_ms, dtype=float).reshape(-1)
    relative_depth = np.asarray(relative_depth, dtype=float).reshape(-1)
    significant = (p_value < alpha) & (n_sessions >= min_sessions) & np.isfinite(mean)
    mask = np.zeros(mean.shape, dtype=np.int8)
    clusters: List[Dict[str, object]] = []
    dt = _axis_step(time_ms)
    drel = _axis_step(relative_depth)
    next_id = 1
    for sign, selected in (
        (1, significant & (mean > 0)),
        (-1, significant & (mean < 0)),
    ):
        labeled, n_labels = ndimage.label(selected, structure=_FOUR_CONNECTED)
        for label in range(1, int(n_labels) + 1):
            depth_index, time_index = np.where(labeled == label)
            duration_ms = (
                float(time_ms[time_index.max()] - time_ms[time_index.min()]) + dt
            )
            depth_fraction = (
                float(
                    relative_depth[depth_index.max()]
                    - relative_depth[depth_index.min()]
                )
                + drel
            )
            if duration_ms < min_duration_ms or depth_fraction < min_depth_fraction:
                continue
            mask[depth_index, time_index] = sign
            values = mean[depth_index, time_index]
            peak_offset = int(np.nanargmax(np.abs(values)))
            peak_depth = int(depth_index[peak_offset])
            peak_time = int(time_index[peak_offset])
            clusters.append(
                {
                    "cluster_id": next_id,
                    "saved_kcsd_sign": _sign_name(sign),
                    "time_start_ms": float(time_ms[time_index.min()]),
                    "time_end_ms": float(time_ms[time_index.max()]),
                    "duration_ms": duration_ms,
                    "relative_depth_start": float(relative_depth[depth_index.min()]),
                    "relative_depth_end": float(relative_depth[depth_index.max()]),
                    "depth_fraction": depth_fraction,
                    "n_bins": int(depth_index.size),
                    "n_sessions_min": int(np.min(n_sessions[depth_index, time_index])),
                    "peak_mean": float(mean[peak_depth, peak_time]),
                    "peak_t": (
                        float(t_stat[peak_depth, peak_time])
                        if t_stat is not None
                        else np.nan
                    ),
                    "peak_p": float(p_value[peak_depth, peak_time]),
                }
            )
            next_id += 1
    return mask, clusters


def _nan_aware_gaussian(values, sigma, *, min_support, ndimage):
    """Smooth finite samples and drop bins whose kernel is mostly unsupported."""
    valid = np.isfinite(values).astype(float)
    filled = np.where(np.isfinite(values), values, 0.0)
    sigma_values = np.atleast_1d(np.asarray(sigma, dtype=float))
    if np.all(sigma_values == 0):
        smoothed = np.array(filled, dtype=float, copy=True)
        support = valid
    else:
        smoothed = ndimage.gaussian_filter(
            filled, sigma=sigma, mode="constant", cval=0.0
        )
        support = ndimage.gaussian_filter(
            valid, sigma=sigma, mode="constant", cval=0.0
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        smoothed = np.divide(
            smoothed,
            support,
            out=np.full(smoothed.shape, np.nan, dtype=float),
            where=support > min_support,
        )
    return smoothed


def _axis_step(values: np.ndarray) -> float:
    diffs = np.diff(np.asarray(values, dtype=float))
    positive = diffs[np.isfinite(diffs) & (diffs > 0)]
    if positive.size == 0:
        return 0.0
    return float(np.median(positive))


def _ensure_writable(path: str) -> None:
    if os.path.exists(path):
        try:
            os.chmod(path, stat.S_IWRITE | stat.S_IREAD)
        except OSError:
            pass


def _remove_tree(path: str) -> bool:
    """
    Delete a folder tree.

    experimentFS often raises PermissionError for ``shutil.rmtree`` when a
    file is open or marked read-only. Returning False lets the caller
    overwrite files instead of stopping.
    """
    if not os.path.exists(path):
        return True

    def _retry(function, failed_path, _exc_info):
        _ensure_writable(failed_path)
        function(failed_path)

    try:
        if os.path.isdir(path):
            shutil.rmtree(path, onerror=_retry)
        else:
            _ensure_writable(path)
            os.remove(path)
        return True
    except PermissionError as exc:
        print(
            f"Could not delete {path}: {exc}. "
            "The existing files will be overwritten. On experimentFS this "
            "usually means a file in that folder is open or read-only.",
            flush=True,
        )
        return False


def _prepare_area_output(output_root: str, *, refresh: bool) -> None:
    """Create the area folder without deleting the session cache."""
    try:
        os.makedirs(output_root, exist_ok=True)
    except PermissionError as exc:
        raise PermissionError(
            f"Cannot create {output_root}. The behavior-aligned CSD folder "
            "on experimentFS is not writable from this kernel."
        ) from exc
    if not refresh:
        return
    for name in os.listdir(output_root):
        if name == SESSION_CACHE_NAME:
            continue
        _remove_tree(os.path.join(output_root, name))


def _processing_metadata(
    *,
    n_relative_depths: int,
    smooth_fwhm_ms: float,
    smooth_fwhm_um: float,
    min_depth_span_um: float,
    depth_grid_step_um: float,
) -> Dict[str, object]:
    return {
        "n_relative_depths": int(n_relative_depths),
        "smooth_fwhm_ms": float(smooth_fwhm_ms),
        "smooth_fwhm_um": float(smooth_fwhm_um),
        "min_depth_span_um": float(min_depth_span_um),
        "depth_grid_step_um": float(depth_grid_step_um),
    }


def _metadata_matches(saved: Dict[str, object], expected: Dict[str, object]) -> bool:
    for key, value in expected.items():
        other = saved.get(key)
        if isinstance(value, float):
            try:
                if other is None or not np.isclose(float(other), value):
                    return False
            except (TypeError, ValueError):
                return False
        elif other != value:
            return False
    return True


def _leaf_cache_dir(cache_root: str, parsed: Dict[str, str]) -> str:
    probe = parsed["probe"] or "NPX1"
    return os.path.join(
        cache_root,
        _safe_path_part(parsed["scope"]),
        _safe_path_part(parsed["behavior"]),
        parsed["alignment"],
        _safe_path_part(parsed["session"]),
        _safe_path_part(probe),
    )


def _save_leaf_cache(
    cache_root: str,
    parsed: Dict[str, str],
    assignment: Dict[str, object],
    processing: Dict[str, object],
) -> None:
    folder = _leaf_cache_dir(cache_root, parsed)
    os.makedirs(folder, exist_ok=True)
    areas = assignment["areas"]
    smoothed = assignment.get("areas_smoothed") or {}
    stored = {}
    for area, matrix in areas.items():
        safe = _safe_path_part(str(area))
        stored[safe] = str(area)
        _write_npy(os.path.join(folder, f"{safe}__csd.npy"), matrix)
        if area in smoothed:
            _write_npy(os.path.join(folder, f"{safe}__smoothed.npy"), smoothed[area])
    _write_npy(os.path.join(folder, "time_ms.npy"), assignment["time_ms"])
    payload = {
        **processing,
        "scope": parsed["scope"],
        "behavior": parsed["behavior"],
        "alignment": parsed["alignment"],
        "session": parsed["session"],
        "probe": parsed["probe"] or "NPX1",
        "areas": stored,
    }
    path = os.path.join(folder, "cache.json")
    _ensure_writable(path)
    with open(path, "w", encoding="utf-8") as stream:
        json.dump(payload, stream)


def _load_leaf_cache(
    folder: str,
    processing: Dict[str, object],
) -> Optional[Dict[str, object]]:
    meta_path = os.path.join(folder, "cache.json")
    if not os.path.isfile(meta_path):
        return None
    with open(meta_path, encoding="utf-8") as stream:
        meta = json.load(stream)
    if not _metadata_matches(meta, processing):
        return None
    areas = {}
    smoothed = {}
    for safe, area in meta.get("areas", {}).items():
        csd_path = os.path.join(folder, f"{safe}__csd.npy")
        if not os.path.isfile(csd_path):
            return None
        areas[area] = np.load(csd_path)
        smoothed_path = os.path.join(folder, f"{safe}__smoothed.npy")
        if os.path.isfile(smoothed_path):
            smoothed[area] = np.load(smoothed_path)
    return {
        "status": "ok",
        "areas": areas,
        "areas_smoothed": smoothed,
        "time_ms": np.load(os.path.join(folder, "time_ms.npy")),
        "parsed": {
            "scope": meta["scope"],
            "behavior": meta["behavior"],
            "alignment": meta["alignment"],
            "session": meta["session"],
            "probe": "" if meta.get("probe") == "NPX1" else meta.get("probe", ""),
        },
    }


def _iter_cache_folders(cache_root: str) -> List[str]:
    found = []
    if not os.path.isdir(cache_root):
        return found
    for dirpath, _dirnames, filenames in os.walk(cache_root):
        if "cache.json" in filenames:
            found.append(dirpath)
    found.sort()
    return found


def _write_npy(path: str, array: np.ndarray) -> None:
    _ensure_writable(path)
    np.save(path, array)


def _copy_overwrite(source: str, destination: str) -> None:
    _ensure_writable(destination)
    shutil.copy2(source, destination)


def _absorb_leaf_assignment(
    parsed: Dict[str, str],
    assignment: Dict[str, object],
    source: str,
    *,
    session_parts: Dict[Tuple[str, str, str, str, str], List[np.ndarray]],
    session_smoothed_parts: Dict[Tuple[str, str, str, str, str], List[np.ndarray]],
    session_time: Dict[Tuple[str, str, str, str, str], np.ndarray],
    skipped: List[Dict[str, str]],
    stop_on_error: bool,
) -> None:
    if assignment["status"] != "ok":
        print(f"  Skipped: {assignment['reason']}", flush=True)
        skipped.append({"leaf": source, "reason": str(assignment["reason"])})
        return
    for area, matrix in assignment["areas"].items():
        key = (
            parsed["scope"],
            parsed["behavior"],
            parsed["alignment"],
            parsed["session"],
            str(area),
        )
        label = (
            f"{parsed['scope']}/{parsed['behavior']} {parsed['alignment']} "
            f"{parsed['session']} {parsed['probe'] or 'NPX1'}"
        )
        if key in session_time and not np.allclose(
            session_time[key], assignment["time_ms"]
        ):
            reason = f"Time axis mismatch for {label} area {area}"
            if stop_on_error:
                raise ValueError(reason)
            print(f"  Skipped area {area}: {reason}", flush=True)
            skipped.append({"leaf": source, "reason": reason})
            continue
        session_parts.setdefault(key, []).append(np.asarray(matrix, dtype=float))
        session_time[key] = np.asarray(assignment["time_ms"], dtype=float)
        smoothed = assignment["areas_smoothed"].get(area)
        if smoothed is not None:
            session_smoothed_parts.setdefault(key, []).append(
                np.asarray(smoothed, dtype=float)
            )


def summarize_area_csd_across_sessions(
    root: PathLike,
    area_table_path: Optional[PathLike] = None,
    *,
    area_rows: Optional[Sequence[Union[AreaLimit, Dict[str, object]]]] = None,
    output_subdir: str = "csd_1d_onset_baseline",
    output_folder_name: str = "area_consistency",
    n_relative_depths: int = 41,
    min_depths: int = 2,
    min_sessions_for_consistency: int = 2,
    consistency_threshold: float = 0.8,
    magnitude_percentile: float = 75.0,
    smooth_fwhm_ms: float = DEFAULT_SMOOTH_FWHM_MS,
    smooth_fwhm_um: float = DEFAULT_SMOOTH_FWHM_UM,
    min_depth_span_um: float = DEFAULT_MIN_DEPTH_SPAN_UM,
    depth_grid_step_um: float = DEFAULT_DEPTH_GRID_STEP_UM,
    min_sessions_for_ttest: int = DEFAULT_MIN_SESSIONS_FOR_TTEST,
    ttest_alpha: float = DEFAULT_TTEST_ALPHA,
    min_cluster_duration_ms: float = DEFAULT_MIN_CLUSTER_DURATION_MS,
    min_cluster_depth_fraction: float = DEFAULT_MIN_CLUSTER_DEPTH_FRACTION,
    save_plots: bool = True,
    refresh: bool = True,
    stop_on_error: bool = False,
    reuse_session_cache: bool = True,
    from_cache_only: bool = False,
) -> Dict[str, object]:
    """
    Pool saved averaged kCSD maps by anatomical area across sessions.

    Requires ``csd_mean_columns.npy`` unless ``from_cache_only`` is set.
    Per-column matrices are not refit and ``KCSD1D`` is never called. Run
    ``average_saved_behavior_kcsd_from_root`` first when a completed kCSD
    run has not yet written the column average.

    The unsmoothed mean and sign-consistency maps are always saved. A second
    map smooths each probe in micrometers and milliseconds, tests session
    means against zero, and outlines clusters that are long enough in time
    and relative depth.

    Accepted session maps are stored in ``session_cache`` before any group
    folder is written. ``refresh`` replaces group outputs and leaves that
    cache in place. If experimentFS refuses to delete an old folder, the
    new files are written over the old ones.
    """
    if n_relative_depths < 2:
        raise ValueError("n_relative_depths must be at least 2")
    if min_depths < 2:
        raise ValueError("min_depths must be at least 2")
    if min_depth_span_um < 0:
        raise ValueError("min_depth_span_um must be non-negative")
    if not 0 < ttest_alpha < 1:
        raise ValueError("ttest_alpha must be between 0 and 1")
    if min_sessions_for_ttest < 2:
        raise ValueError("min_sessions_for_ttest must be at least 2")
    root = os.path.abspath(os.fspath(root))
    processing = _processing_metadata(
        n_relative_depths=n_relative_depths,
        smooth_fwhm_ms=smooth_fwhm_ms,
        smooth_fwhm_um=smooth_fwhm_um,
        min_depth_span_um=min_depth_span_um,
        depth_grid_step_um=depth_grid_step_um,
    )
    cache_root = os.path.join(root, output_folder_name, SESSION_CACHE_NAME)
    rows_by_animal: Dict[str, List[AreaLimit]] = {}
    if not from_cache_only:
        rows = _coerce_area_rows(area_rows, area_table_path)
        for row in rows:
            rows_by_animal.setdefault(row.animal_name, []).append(row)

    discovered = {"means": [], "missing_means": []}
    if not from_cache_only:
        discovered = _discover_averaged_leaves(
            root,
            output_subdir=output_subdir,
            output_folder_name=output_folder_name,
        )
        if not discovered["means"]:
            extra = ""
            if discovered["missing_means"]:
                extra = (
                    f" Found {len(discovered['missing_means'])} folders with "
                    "per-column csd_x*.npy but no csd_mean_columns.npy."
                )
            raise FileNotFoundError(
                "No saved csd_mean_columns.npy maps under "
                f"{root}.{extra} Run average_saved_behavior_kcsd_from_root "
                "before the area summary. This step does not estimate kCSD."
            )

    if from_cache_only:
        print(
            "Resuming area consistency from session_cache. "
            "csd_mean_columns.npy will not be read and kCSD will NOT be rerun.",
            flush=True,
        )
    else:
        print(
            "Area consistency reads saved csd_mean_columns.npy only; "
            "kCSD will NOT be rerun",
            flush=True,
        )
    print(
        "Sign is the saved kCSD polarity (RdBu_r). It is not relabeled "
        "as sink or source.",
        flush=True,
    )
    if discovered["missing_means"]:
        print(
            f"Skipping {len(discovered['missing_means'])} folders that still "
            "need column averaging (no csd_mean_columns.npy).",
            flush=True,
        )

    relative_depth = np.linspace(0.0, 1.0, int(n_relative_depths))
    session_parts: Dict[Tuple[str, str, str, str, str], List[np.ndarray]] = {}
    session_smoothed_parts: Dict[Tuple[str, str, str, str, str], List[np.ndarray]] = {}
    session_time: Dict[Tuple[str, str, str, str, str], np.ndarray] = {}
    skipped: List[Dict[str, str]] = []

    if from_cache_only:
        folders = _iter_cache_folders(cache_root)
        if not folders:
            raise FileNotFoundError(
                f"No session cache at {cache_root}. The previous run stopped "
                "before any cached session map was written, so it cannot be "
                "resumed. Rerun summarize_area_csd_across_sessions once. That "
                "cell saves each session while it reads csd_mean_columns.npy. "
                "After the cache exists, this step rebuilds the means, "
                "clusters, and figure copies without reading those matrices."
            )
        loaded_any = False
        for folder_number, folder in enumerate(folders, start=1):
            label = os.path.relpath(folder, cache_root)
            print(f"[{folder_number}/{len(folders)}] Cache {label}", flush=True)
            try:
                loaded = _load_leaf_cache(folder, processing)
            except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
                print(f"  Skipped: {exc}", flush=True)
                skipped.append({"leaf": folder, "reason": str(exc)})
                continue
            if loaded is None:
                reason = (
                    "Session cache does not match the current smoothing, "
                    "depth-span, or relative-depth settings"
                )
                print(f"  Skipped: {reason}", flush=True)
                skipped.append({"leaf": folder, "reason": reason})
                continue
            loaded_any = True
            _absorb_leaf_assignment(
                loaded["parsed"],
                loaded,
                folder,
                session_parts=session_parts,
                session_smoothed_parts=session_smoothed_parts,
                session_time=session_time,
                skipped=skipped,
                stop_on_error=stop_on_error,
            )
        if not loaded_any:
            raise FileNotFoundError(
                f"No usable session cache at {cache_root}. The saved maps "
                "do not match the current smoothing, depth-span, or "
                "relative-depth settings, or the cache files could not be "
                "read. Rerun summarize_area_csd_across_sessions so the cache "
                "is rebuilt."
            )
    else:
        for leaf_number, leaf in enumerate(discovered["means"], start=1):
            parsed = _parse_mean_leaf(root, leaf, output_subdir)
            label = (
                f"{parsed['scope']}/{parsed['behavior']} {parsed['alignment']} "
                f"{parsed['session']} {parsed['probe'] or 'NPX1'}"
            )
            print(
                f"[{leaf_number}/{len(discovered['means'])}] Reading {label}",
                flush=True,
            )
            cache_folder = _leaf_cache_dir(cache_root, parsed)
            assignment = None
            if reuse_session_cache:
                assignment = _load_leaf_cache(cache_folder, processing)
                if assignment is not None:
                    print(f"  Using saved session cache: {label}", flush=True)
            if assignment is None:
                try:
                    assignment = _segment_leaf(
                        leaf,
                        parsed,
                        rows_by_animal,
                        relative_depth,
                        min_depths=min_depths,
                        min_depth_span_um=min_depth_span_um,
                        smooth_fwhm_ms=smooth_fwhm_ms,
                        smooth_fwhm_um=smooth_fwhm_um,
                        depth_grid_step_um=depth_grid_step_um,
                    )
                except Exception as exc:
                    if stop_on_error:
                        raise
                    print(f"  Skipped: {exc}", flush=True)
                    skipped.append({"leaf": leaf, "reason": str(exc)})
                    continue
                if assignment["status"] == "ok":
                    _save_leaf_cache(cache_root, parsed, assignment, processing)
            _absorb_leaf_assignment(
                parsed,
                assignment,
                leaf,
                session_parts=session_parts,
                session_smoothed_parts=session_smoothed_parts,
                session_time=session_time,
                skipped=skipped,
                stop_on_error=stop_on_error,
            )

    grouped: Dict[Tuple[str, str, str, str], List[Tuple[str, np.ndarray, np.ndarray]]] = {}
    grouped_smoothed: Dict[Tuple[str, str, str, str], Dict[str, np.ndarray]] = {}
    for key, matrices in session_parts.items():
        scope, behavior, alignment, session, area = key
        session_mean = _nanmean_stack(np.stack(matrices, axis=0))
        if not np.any(np.isfinite(session_mean)):
            continue
        group_key = (scope, behavior, alignment, area)
        grouped.setdefault(group_key, []).append(
            (session, session_mean, session_time[key])
        )
        if key in session_smoothed_parts:
            smoothed_mean = _nanmean_stack(
                np.stack(session_smoothed_parts[key], axis=0)
            )
            if np.any(np.isfinite(smoothed_mean)):
                grouped_smoothed.setdefault(group_key, {})[session] = smoothed_mean

    output_root = os.path.join(root, output_folder_name)
    _prepare_area_output(output_root, refresh=refresh)

    summaries = []
    stable_rows: List[Dict[str, object]] = []
    cluster_rows: List[Dict[str, object]] = []
    plot_groups: Dict[Tuple[str, str, str], Dict[str, Dict[str, object]]] = {}
    for group_key in sorted(grouped):
        scope, behavior, alignment, area = group_key
        observations = sorted(grouped[group_key], key=lambda item: item[0])
        reference_time = observations[0][2]
        usable = []
        for session, matrix, time_ms in observations:
            if np.allclose(time_ms, reference_time):
                usable.append((session, matrix))
            else:
                reason = (
                    f"{scope}/{behavior} {alignment} {area}: session {session} "
                    "time axis differs from the first session"
                )
                print(f"  Skipped: {reason}", flush=True)
                skipped.append({"leaf": session, "reason": reason})
        if not usable:
            continue
        smoothed_by_session = grouped_smoothed.get(group_key, {})
        smoothed_stack = np.stack(
            [
                smoothed_by_session.get(
                    session, np.full(matrix.shape, np.nan, dtype=float)
                )
                for session, matrix in usable
            ],
            axis=0,
        )
        inference = infer_saved_kcsd_clusters(
            smoothed_stack,
            reference_time,
            relative_depth,
            alpha=ttest_alpha,
            min_sessions=min_sessions_for_ttest,
            min_duration_ms=min_cluster_duration_ms,
            min_depth_fraction=min_cluster_depth_fraction,
        )
        print(
            f"Processing {scope}/{behavior} {alignment} {area}: "
            f"{len(usable)} sessions, {len(inference['clusters'])} clusters",
            flush=True,
        )
        stack = np.stack([matrix for _session, matrix in usable], axis=0)
        stats = _across_session_stats(stack)
        folder = os.path.join(
            output_root,
            _safe_path_part(scope),
            _safe_path_part(behavior),
            alignment,
            _safe_path_part(area),
        )
        os.makedirs(folder, exist_ok=True)
        _save_group_arrays(folder, stats, relative_depth, reference_time)
        _save_inference_arrays(folder, inference)
        for cluster in inference["clusters"]:
            cluster_rows.append(
                {
                    "scope": scope,
                    "behavior": behavior,
                    "alignment": alignment,
                    "area": area,
                    **cluster,
                }
            )
        sessions_path = os.path.join(folder, "sessions.txt")
        _ensure_writable(sessions_path)
        with open(sessions_path, "w", encoding="utf-8") as stream:
            stream.write("".join(f"{session}\n" for session, _matrix in usable))
        summary, feature_rows = _group_tables(
            scope,
            behavior,
            alignment,
            area,
            [session for session, _matrix in usable],
            stats,
            relative_depth,
            reference_time,
            min_sessions_for_consistency=min_sessions_for_consistency,
            consistency_threshold=consistency_threshold,
            magnitude_percentile=magnitude_percentile,
        )
        summary["folder"] = folder
        summary["n_clusters"] = len(inference["clusters"])
        summaries.append(summary)
        stable_rows.extend(feature_rows)
        plot_groups.setdefault((scope, behavior, area), {})[alignment] = {
            "folder": folder,
            "mean": stats["mean"],
            "consistency": stats["sign_consistency"],
            "n_sessions": stats["n_sessions"],
            "n_sessions_total": len(usable),
            "smoothed_mean": inference["mean"],
            "cluster_mask": inference["mask"],
            "smooth_fwhm_ms": smooth_fwhm_ms,
            "smooth_fwhm_um": smooth_fwhm_um,
            "min_sessions_for_ttest": min_sessions_for_ttest,
            "ttest_alpha": ttest_alpha,
            "min_cluster_duration_ms": min_cluster_duration_ms,
            "min_cluster_depth_fraction": min_cluster_depth_fraction,
        }

    color_limits = _shared_alignment_color_limits(plot_groups)
    if save_plots:
        _write_area_figures(
            plot_groups,
            color_limits,
            min_sessions_for_consistency=min_sessions_for_consistency,
        )

    summary_csv = os.path.join(output_root, "summary.csv")
    stable_csv = os.path.join(output_root, "stable_features.csv")
    _write_csv(
        summary_csv,
        [
            "scope",
            "behavior",
            "alignment",
            "area",
            "n_sessions",
            "sessions",
            "max_abs_mean",
            "relative_depth_at_max_abs",
            "time_ms_at_max_abs",
            "saved_kcsd_sign_at_max_abs",
            "sign_consistency_at_max_abs",
            "n_stable_bins",
            "n_clusters",
        ],
        summaries,
    )
    clusters_csv = os.path.join(output_root, "clusters.csv")
    _write_csv(
        clusters_csv,
        [
            "scope",
            "behavior",
            "alignment",
            "area",
            "cluster_id",
            "saved_kcsd_sign",
            "time_start_ms",
            "time_end_ms",
            "duration_ms",
            "relative_depth_start",
            "relative_depth_end",
            "depth_fraction",
            "n_bins",
            "n_sessions_min",
            "peak_mean",
            "peak_t",
            "peak_p",
        ],
        cluster_rows,
    )
    _write_csv(
        stable_csv,
        [
            "scope",
            "behavior",
            "alignment",
            "area",
            "relative_depth",
            "time_ms",
            "mean_csd",
            "saved_kcsd_sign",
            "sign_consistency",
            "n_sessions",
            "fraction_positive",
            "fraction_negative",
        ],
        stable_rows,
    )
    print(
        f"Area summary complete: {len(summaries)} area maps, "
        f"{len(stable_rows)} stable bins, {len(cluster_rows)} clusters, "
        f"index={summary_csv}",
        flush=True,
    )
    sorted_figures = organize_area_consistency_figures(
        root,
        area_folder_name=output_folder_name,
        refresh=True,
        missing_ok=True,
    )
    return {
        "output_root": output_root,
        "summary_csv": summary_csv,
        "stable_features_csv": stable_csv,
        "clusters_csv": clusters_csv,
        "sorted_figures": sorted_figures,
        "n_groups": len(summaries),
        "n_missing_mean": len(discovered["missing_means"]),
        "missing_mean_folders": discovered["missing_means"],
        "skipped": skipped,
        "groups": summaries,
        "color_limits": color_limits,
    }


def finalize_area_consistency_from_cache(
    root: PathLike,
    *,
    output_folder_name: str = "area_consistency",
    n_relative_depths: int = 41,
    min_sessions_for_consistency: int = 2,
    consistency_threshold: float = 0.8,
    magnitude_percentile: float = 75.0,
    smooth_fwhm_ms: float = DEFAULT_SMOOTH_FWHM_MS,
    smooth_fwhm_um: float = DEFAULT_SMOOTH_FWHM_UM,
    min_depth_span_um: float = DEFAULT_MIN_DEPTH_SPAN_UM,
    depth_grid_step_um: float = DEFAULT_DEPTH_GRID_STEP_UM,
    min_sessions_for_ttest: int = DEFAULT_MIN_SESSIONS_FOR_TTEST,
    ttest_alpha: float = DEFAULT_TTEST_ALPHA,
    min_cluster_duration_ms: float = DEFAULT_MIN_CLUSTER_DURATION_MS,
    min_cluster_depth_fraction: float = DEFAULT_MIN_CLUSTER_DEPTH_FRACTION,
    save_plots: bool = True,
    refresh: bool = True,
) -> Dict[str, object]:
    """
    Write area means, clusters, and sorted figures from ``session_cache``.

    Does not read ``csd_mean_columns.npy``, the area spreadsheet, or refit
    kCSD. The cache is filled by ``summarize_area_csd_across_sessions`` as
    each saved map is read. A run that stopped before the first cache file
    was written has nothing to resume; rerun that function once.
    """
    return summarize_area_csd_across_sessions(
        root,
        output_folder_name=output_folder_name,
        n_relative_depths=n_relative_depths,
        min_sessions_for_consistency=min_sessions_for_consistency,
        consistency_threshold=consistency_threshold,
        magnitude_percentile=magnitude_percentile,
        smooth_fwhm_ms=smooth_fwhm_ms,
        smooth_fwhm_um=smooth_fwhm_um,
        min_depth_span_um=min_depth_span_um,
        depth_grid_step_um=depth_grid_step_um,
        min_sessions_for_ttest=min_sessions_for_ttest,
        ttest_alpha=ttest_alpha,
        min_cluster_duration_ms=min_cluster_duration_ms,
        min_cluster_depth_fraction=min_cluster_depth_fraction,
        save_plots=save_plots,
        refresh=refresh,
        from_cache_only=True,
        reuse_session_cache=True,
    )


def organize_area_consistency_figures(
    root: PathLike,
    *,
    area_folder_name: str = "area_consistency",
    output_folder_name: str = "by_behavior",
    refresh: bool = True,
    missing_ok: bool = False,
) -> Dict[str, object]:
    """
    Copy merged area PNGs into one onset folder and one offset folder.

    Original area folders stay in place. Filenames start with the behavior,
    so a file browser lists the same behavior together:

    ``by_behavior/onset/play__category__LPAG__mean_csd.png``
    """
    root = os.path.abspath(os.fspath(root))
    area_root = os.path.join(root, area_folder_name)
    if not os.path.isdir(area_root):
        raise FileNotFoundError(
            f"No merged area maps at {area_root}. Run "
            "summarize_area_csd_across_sessions first."
        )

    sources = _discover_area_figures(area_root, output_folder_name)
    if not sources:
        message = (
            "No mean_csd.png, sign_consistency.png, or "
            f"smoothed_mean_clusters.png files under {area_root}"
        )
        if missing_ok:
            print(message, flush=True)
            return {
                "output_root": os.path.join(area_root, output_folder_name),
                "index_path": None,
                "n_copies": 0,
            }
        raise FileNotFoundError(message)

    output_root = os.path.join(area_root, output_folder_name)
    if refresh:
        _remove_tree(output_root)
    os.makedirs(output_root, exist_ok=True)

    rows: List[Dict[str, str]] = []
    destinations: Dict[str, str] = {}
    print(
        f"Copying {len(sources)} merged area figures into {output_root}",
        flush=True,
    )
    for source in sources:
        relative = os.path.relpath(source["folder"], area_root).split(os.sep)
        scope, behavior, alignment, area = relative
        kind = source["kind"]
        filename = (
            f"{_safe_path_part(behavior)}__{_safe_path_part(scope)}__"
            f"{_safe_path_part(area)}__{kind}.png"
        )
        destination = os.path.join(output_root, alignment, filename)
        prior = destinations.get(destination)
        if prior is not None and prior != source["path"]:
            raise RuntimeError(
                "Area-figure name collision: "
                f"{prior} and {source['path']} both map to {destination}"
            )
        destinations[destination] = source["path"]
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        _copy_overwrite(source["path"], destination)
        print(
            f"  {alignment}: {filename}",
            flush=True,
        )
        rows.append(
            {
                "alignment": alignment,
                "scope": scope,
                "behavior": behavior,
                "area": area,
                "kind": kind,
                "original_file": os.path.relpath(source["path"], root),
                "sorted_file": os.path.relpath(destination, root),
            }
        )

    index_path = os.path.join(output_root, "figure_index.csv")
    _write_csv(
        index_path,
        [
            "alignment",
            "scope",
            "behavior",
            "area",
            "kind",
            "original_file",
            "sorted_file",
        ],
        rows,
    )
    print(
        f"Merged area figures ready: {len(rows)} copies in onset/offset folders, "
        f"index={index_path}",
        flush=True,
    )
    return {
        "output_root": output_root,
        "index_path": index_path,
        "n_copies": len(rows),
    }


def _discover_area_figures(
    area_root: str,
    output_folder_name: str,
) -> List[Dict[str, str]]:
    figure_names = {
        "mean_csd.png": "mean_csd",
        "sign_consistency.png": "sign_consistency",
        "smoothed_mean_clusters.png": "smoothed_clusters",
    }
    found: List[Dict[str, str]] = []
    for dirpath, _dirnames, filenames in os.walk(area_root):
        relative_parts = os.path.relpath(dirpath, area_root).split(os.sep)
        if relative_parts[0] in {output_folder_name, SESSION_CACHE_NAME}:
            continue
        if len(relative_parts) != 4 or relative_parts[2] not in {"onset", "offset"}:
            continue
        for filename, kind in figure_names.items():
            if filename in filenames:
                found.append(
                    {
                        "folder": dirpath,
                        "path": os.path.join(dirpath, filename),
                        "kind": kind,
                    }
                )
    found.sort(key=lambda item: item["path"])
    return found


def _coerce_area_rows(
    area_rows: Optional[Sequence[Union[AreaLimit, Dict[str, object]]]],
    area_table_path: Optional[PathLike],
) -> List[AreaLimit]:
    if area_rows is not None:
        coerced = []
        for row in area_rows:
            if isinstance(row, AreaLimit):
                coerced.append(row)
            else:
                coerced.append(AreaLimit(**row))
        return coerced
    if area_table_path is None:
        raise ValueError("Pass area_table_path or area_rows")
    return load_area_limits(area_table_path)


def _cell_text(values: Sequence[object], index: int) -> str:
    if index >= len(values) or values[index] is None:
        return ""
    return str(values[index]).strip()


def _as_probe_num(value: object) -> Optional[int]:
    if value is None or value == "":
        return None
    number = float(value)
    if not np.isfinite(number):
        return None
    rounded = round(number)
    if not np.isclose(number, rounded):
        raise ValueError(f"ProbeNum is not an integer: {value}")
    return int(rounded)


def _discover_averaged_leaves(
    root: str,
    *,
    output_subdir: str,
    output_folder_name: str,
) -> Dict[str, List[str]]:
    means: List[str] = []
    missing_means: List[str] = []
    for dirpath, _dirnames, filenames in os.walk(root):
        if os.path.basename(dirpath) != output_subdir:
            continue
        relative = os.path.relpath(dirpath, root).split(os.sep)
        if output_folder_name in relative:
            continue
        names = set(filenames)
        if "csd_mean_columns.npy" in names:
            means.append(dirpath)
        elif any(
            name.startswith("csd_x") and name.endswith(".npy") for name in names
        ):
            missing_means.append(dirpath)
    means.sort()
    missing_means.sort()
    return {"means": means, "missing_means": missing_means}


def _parse_mean_leaf(root: str, leaf: str, output_subdir: str) -> Dict[str, str]:
    if os.path.basename(leaf) != output_subdir:
        raise ValueError(f"Expected maps inside {output_subdir}: {leaf}")
    alignment_leaf = os.path.dirname(leaf)
    parts = os.path.relpath(alignment_leaf, root).split(os.sep)
    if len(parts) not in (4, 5) or parts[2] not in {"onset", "offset"}:
        raise ValueError(
            "Unrecognized behavior-aligned CSD path "
            f"{os.path.relpath(alignment_leaf, root)}"
        )
    return {
        "scope": parts[0],
        "behavior": parts[1],
        "alignment": parts[2],
        "session": parts[3],
        "probe": parts[4] if len(parts) == 5 else "",
    }


def _probe_identity(probe_label: str) -> Tuple[int, int]:
    if probe_label == "" or probe_label.upper() == "NPX1":
        return 1, 1
    match = re.fullmatch(r"probe(\d+)", probe_label, flags=re.IGNORECASE)
    if match is None:
        raise ValueError(f"Unrecognized probe folder: {probe_label}")
    return 2, int(match.group(1))


def _segment_leaf(
    leaf: str,
    parsed: Dict[str, str],
    rows_by_animal: Dict[str, List[AreaLimit]],
    relative_depth: np.ndarray,
    *,
    min_depths: int,
    min_depth_span_um: float,
    smooth_fwhm_ms: float,
    smooth_fwhm_um: float,
    depth_grid_step_um: float,
) -> Dict[str, object]:
    csd = np.load(os.path.join(leaf, "csd_mean_columns.npy"))
    y = np.load(os.path.join(leaf, "y_mean_columns.npy")).reshape(-1)
    time_ms = np.load(os.path.join(leaf, "time_mean_columns.npy")).reshape(-1)
    if csd.ndim != 2 or csd.shape[0] != len(y) or csd.shape[1] != len(time_ms):
        raise ValueError(
            f"Saved mean CSD shape {csd.shape} does not match y={len(y)} "
            f"and time={len(time_ms)}"
        )
    try:
        animal_key = animal_key_from_session(parsed["session"])
    except ValueError as exc:
        return {"status": "skip", "reason": str(exc)}
    animal_rows = rows_by_animal.get(animal_key)
    if not animal_rows:
        return {
            "status": "skip",
            "reason": (
                f"No area-table rows for {animal_key} "
                f"(session {parsed['session']})"
            ),
        }
    npx_type, probe_num = _probe_identity(parsed["probe"])
    grouped: Dict[Tuple[str, float, float], List[int]] = {}
    spans: Dict[Tuple[str, float, float], AreaLimit] = {}
    for index, depth in enumerate(y):
        chosen = select_area_for_depth(
            float(depth),
            animal_rows,
            npx_type=npx_type,
            probe_num=probe_num,
        )
        if chosen is None:
            continue
        identity = (chosen.area, chosen.depth_start, chosen.depth_end)
        grouped.setdefault(identity, []).append(index)
        spans[identity] = chosen

    areas: Dict[str, np.ndarray] = {}
    area_pieces: Dict[str, List[np.ndarray]] = {}
    smoothed_pieces: Dict[str, List[np.ndarray]] = {}
    n_short = 0
    for identity, indices in grouped.items():
        if len(indices) < min_depths:
            continue
        chosen = spans[identity]
        y_segment = y[indices]
        observed_span = float(np.nanmax(y_segment) - np.nanmin(y_segment))
        if observed_span < min_depth_span_um:
            n_short += 1
            continue
        csd_segment = csd[indices]
        resampled = resample_area_segment(
            y_segment,
            csd_segment,
            chosen.depth_start,
            chosen.depth_end,
            relative_depth,
        )
        if resampled is None:
            continue
        area_pieces.setdefault(chosen.area, []).append(resampled)
        smoothed = smooth_csd_in_depth_and_time(
            y_segment,
            csd_segment,
            time_ms,
            fwhm_um=smooth_fwhm_um,
            fwhm_ms=smooth_fwhm_ms,
            grid_step_um=depth_grid_step_um,
        )
        if smoothed is None:
            continue
        smoothed_y, smoothed_csd = smoothed
        smoothed_relative = resample_area_segment(
            smoothed_y,
            smoothed_csd,
            chosen.depth_start,
            chosen.depth_end,
            relative_depth,
        )
        if smoothed_relative is not None:
            smoothed_pieces.setdefault(chosen.area, []).append(smoothed_relative)
    for area, pieces in area_pieces.items():
        areas[area] = _nanmean_stack(np.stack(pieces, axis=0))
    areas_smoothed = {
        area: _nanmean_stack(np.stack(pieces, axis=0))
        for area, pieces in smoothed_pieces.items()
    }
    if not areas:
        if n_short:
            reason = (
                f"No area spanned at least {min_depth_span_um:g} µm "
                f"in {parsed['session']}"
            )
        else:
            reason = (
                f"No area had at least {min_depths} depths in {parsed['session']}"
            )
        return {"status": "skip", "reason": reason}
    return {
        "status": "ok",
        "areas": areas,
        "areas_smoothed": areas_smoothed,
        "time_ms": time_ms,
    }


def _interpolate_rows(
    src_x: np.ndarray,
    values: np.ndarray,
    query_x: np.ndarray,
) -> Optional[np.ndarray]:
    order = np.argsort(src_x)
    src_x = np.asarray(src_x, dtype=float)[order]
    values = np.asarray(values, dtype=float)[order]
    unique_x, inverse = np.unique(src_x, return_inverse=True)
    if len(unique_x) < len(src_x):
        totals = np.zeros((len(unique_x), values.shape[1]), dtype=float)
        counts = np.zeros(len(unique_x), dtype=float)
        for index, group in enumerate(inverse):
            totals[group] += values[index]
            counts[group] += 1.0
        values = totals / counts[:, None]
        src_x = unique_x
    if len(src_x) < 2:
        return None
    out = np.empty((len(query_x), values.shape[1]), dtype=float)
    for time_index in range(values.shape[1]):
        column = np.interp(query_x, src_x, values[:, time_index])
        column[(query_x < src_x[0]) | (query_x > src_x[-1])] = np.nan
        out[:, time_index] = column
    return out


def _nanmean_stack(stack: np.ndarray) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(stack, axis=0)


def _across_session_stats(stack: np.ndarray) -> Dict[str, np.ndarray]:
    finite = np.isfinite(stack)
    n_sessions = np.sum(finite, axis=0).astype(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean = np.nanmean(stack, axis=0)
        std = np.nanstd(stack, axis=0, ddof=1)
    sem = np.full(mean.shape, np.nan, dtype=float)
    enough = n_sessions >= 2
    sem[enough] = std[enough] / np.sqrt(n_sessions[enough])
    positive = np.sum(stack > 0, axis=0).astype(float)
    negative = np.sum(stack < 0, axis=0).astype(float)
    fraction_positive = np.full(mean.shape, np.nan, dtype=float)
    fraction_negative = np.full(mean.shape, np.nan, dtype=float)
    present = n_sessions > 0
    fraction_positive[present] = positive[present] / n_sessions[present]
    fraction_negative[present] = negative[present] / n_sessions[present]
    sign_consistency = np.fmax(fraction_positive, fraction_negative)
    covered = np.any(finite, axis=2)
    return {
        "mean": mean,
        "sem": sem,
        "fraction_positive": fraction_positive,
        "fraction_negative": fraction_negative,
        "sign_consistency": sign_consistency,
        "n_sessions": n_sessions,
        "n_sessions_by_depth": np.sum(covered, axis=0).astype(float),
    }


def _save_group_arrays(
    folder: str,
    stats: Dict[str, np.ndarray],
    relative_depth: np.ndarray,
    time_ms: np.ndarray,
) -> None:
    _write_npy(os.path.join(folder, "mean_csd.npy"), stats["mean"])
    _write_npy(os.path.join(folder, "sem_csd.npy"), stats["sem"])
    _write_npy(os.path.join(folder, "fraction_positive.npy"), stats["fraction_positive"])
    _write_npy(os.path.join(folder, "fraction_negative.npy"), stats["fraction_negative"])
    _write_npy(os.path.join(folder, "sign_consistency.npy"), stats["sign_consistency"])
    _write_npy(
        os.path.join(folder, "n_sessions_by_depth.npy"),
        stats["n_sessions_by_depth"],
    )
    _write_npy(os.path.join(folder, "n_sessions.npy"), stats["n_sessions"])
    _write_npy(os.path.join(folder, "relative_depth.npy"), relative_depth)
    _write_npy(os.path.join(folder, "time_ms.npy"), time_ms)


def _save_inference_arrays(folder: str, inference: Dict[str, object]) -> None:
    _write_npy(os.path.join(folder, "smoothed_mean_csd.npy"), inference["mean"])
    _write_npy(os.path.join(folder, "t_statistic.npy"), inference["t"])
    _write_npy(os.path.join(folder, "p_value.npy"), inference["p"])
    _write_npy(os.path.join(folder, "cluster_mask.npy"), inference["mask"])
    _write_csv(
        os.path.join(folder, "clusters.csv"),
        [
            "cluster_id",
            "saved_kcsd_sign",
            "time_start_ms",
            "time_end_ms",
            "duration_ms",
            "relative_depth_start",
            "relative_depth_end",
            "depth_fraction",
            "n_bins",
            "n_sessions_min",
            "peak_mean",
            "peak_t",
            "peak_p",
        ],
        inference["clusters"],
    )


def _group_tables(
    scope: str,
    behavior: str,
    alignment: str,
    area: str,
    sessions: Sequence[str],
    stats: Dict[str, np.ndarray],
    relative_depth: np.ndarray,
    time_ms: np.ndarray,
    *,
    min_sessions_for_consistency: int,
    consistency_threshold: float,
    magnitude_percentile: float,
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    mean = stats["mean"]
    consistency = stats["sign_consistency"]
    n_sessions = stats["n_sessions"]
    finite = np.isfinite(mean)
    if np.any(finite):
        flat_index = int(np.nanargmax(np.abs(mean)))
        depth_index, time_index = np.unravel_index(flat_index, mean.shape)
        peak_mean = float(mean[depth_index, time_index])
        peak = {
            "max_abs_mean": abs(peak_mean),
            "relative_depth_at_max_abs": float(relative_depth[depth_index]),
            "time_ms_at_max_abs": float(time_ms[time_index]),
            "saved_kcsd_sign_at_max_abs": _sign_name(peak_mean),
            "sign_consistency_at_max_abs": float(consistency[depth_index, time_index]),
        }
    else:
        peak = {
            "max_abs_mean": np.nan,
            "relative_depth_at_max_abs": np.nan,
            "time_ms_at_max_abs": np.nan,
            "saved_kcsd_sign_at_max_abs": "",
            "sign_consistency_at_max_abs": np.nan,
        }
    eligible = (
        finite
        & np.isfinite(consistency)
        & (n_sessions >= min_sessions_for_consistency)
    )
    feature_rows: List[Dict[str, object]] = []
    if np.any(eligible):
        threshold = float(np.percentile(np.abs(mean[eligible]), magnitude_percentile))
        selected = (
            eligible
            & (consistency >= consistency_threshold)
            & (np.abs(mean) >= threshold)
        )
        depth_indices, time_indices = np.where(selected)
        order = np.lexsort((time_indices, depth_indices))
        for depth_index, time_index in zip(depth_indices[order], time_indices[order]):
            value = float(mean[depth_index, time_index])
            feature_rows.append(
                {
                    "scope": scope,
                    "behavior": behavior,
                    "alignment": alignment,
                    "area": area,
                    "relative_depth": float(relative_depth[depth_index]),
                    "time_ms": float(time_ms[time_index]),
                    "mean_csd": value,
                    "saved_kcsd_sign": _sign_name(value),
                    "sign_consistency": float(consistency[depth_index, time_index]),
                    "n_sessions": int(n_sessions[depth_index, time_index]),
                    "fraction_positive": float(
                        stats["fraction_positive"][depth_index, time_index]
                    ),
                    "fraction_negative": float(
                        stats["fraction_negative"][depth_index, time_index]
                    ),
                }
            )
    summary = {
        "scope": scope,
        "behavior": behavior,
        "alignment": alignment,
        "area": area,
        "n_sessions": len(sessions),
        "sessions": "; ".join(sessions),
        "n_stable_bins": len(feature_rows),
        **peak,
    }
    return summary, feature_rows


def _sign_name(value: float) -> str:
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "zero"


def _symmetric_limit(images: Sequence[np.ndarray]) -> Tuple[float, float]:
    peak = 0.0
    for image in images:
        finite = np.abs(np.asarray(image, dtype=float))
        finite = finite[np.isfinite(finite)]
        if finite.size:
            peak = max(peak, float(np.max(finite)))
    if not np.isfinite(peak) or peak <= 0:
        peak = 1.0
    return (-peak, peak)


def _shared_alignment_color_limits(
    plot_groups: Dict[Tuple[str, str, str], Dict[str, Dict[str, object]]],
) -> Dict[str, Tuple[float, float]]:
    """Use one symmetric color scale for onset and offset of the same area."""
    color_limits: Dict[str, Tuple[float, float]] = {}
    for (scope, behavior, area), alignments in plot_groups.items():
        peak = 0.0
        for payload in alignments.values():
            finite = np.abs(payload["mean"][np.isfinite(payload["mean"])])
            if finite.size:
                peak = max(peak, float(np.max(finite)))
        if not np.isfinite(peak) or peak <= 0:
            peak = 1.0
        color_limits[f"{scope}/{behavior}/{area}"] = (-peak, peak)
    return color_limits


def _write_area_figures(
    plot_groups: Dict[Tuple[str, str, str], Dict[str, Dict[str, object]]],
    color_limits: Dict[str, Tuple[float, float]],
    *,
    min_sessions_for_consistency: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    for (scope, behavior, area), alignments in sorted(plot_groups.items()):
        clim = color_limits[f"{scope}/{behavior}/{area}"]
        for alignment, payload in sorted(alignments.items()):
            folder = str(payload["folder"])
            time_ms = np.load(os.path.join(folder, "time_ms.npy"))
            relative_depth = np.load(os.path.join(folder, "relative_depth.npy"))
            _save_heatmap(
                plt,
                payload["mean"],
                time_ms,
                relative_depth,
                title=(
                    f"{scope}/{behavior} — {alignment} — {area}\n"
                    f"{payload['n_sessions_total']} sessions; saved kCSD sign"
                ),
                cbar_label="kCSD (a.u.)",
                cmap="RdBu_r",
                clim=clim,
                save_path=os.path.join(folder, "mean_csd.png"),
            )
            shown_consistency = np.array(payload["consistency"], dtype=float, copy=True)
            shown_consistency[payload["n_sessions"] < min_sessions_for_consistency] = (
                np.nan
            )
            _save_heatmap(
                plt,
                shown_consistency,
                time_ms,
                relative_depth,
                title=(
                    f"{scope}/{behavior} — {alignment} — {area}\n"
                    "sign consistency; blank where fewer than "
                    f"{min_sessions_for_consistency} sessions"
                ),
                cbar_label="fraction of sessions with the majority saved sign",
                cmap="viridis",
                clim=(0.0, 1.0),
                save_path=os.path.join(folder, "sign_consistency.png"),
            )
            smoothed_clim = _symmetric_limit(
                [
                    payload["smoothed_mean"]
                    for payload in alignments.values()
                    if payload.get("smoothed_mean") is not None
                ]
            )
            _save_heatmap(
                plt,
                payload["smoothed_mean"],
                time_ms,
                relative_depth,
                title=(
                    f"{scope}/{behavior} — {alignment} — {area}\n"
                    f"smoothed FWHM {payload['smooth_fwhm_ms']:g} ms x "
                    f"{payload['smooth_fwhm_um']:g} µm; "
                    f"clusters require n≥{payload['min_sessions_for_ttest']}, "
                    f"p<{payload['ttest_alpha']:g}, "
                    f"≥{payload['min_cluster_duration_ms']:g} ms and "
                    f"≥{payload['min_cluster_depth_fraction']:.0%} of depth"
                ),
                cbar_label="smoothed kCSD (a.u.)",
                cmap="RdBu_r",
                clim=smoothed_clim,
                save_path=os.path.join(folder, "smoothed_mean_clusters.png"),
                contour=payload["cluster_mask"],
            )
            print(f"  Saved area figures: {folder}", flush=True)


def _save_heatmap(
    plt,
    image: np.ndarray,
    time_ms: np.ndarray,
    relative_depth: np.ndarray,
    *,
    title: str,
    cbar_label: str,
    cmap: str,
    clim: Tuple[float, float],
    save_path: str,
    contour: Optional[np.ndarray] = None,
) -> None:
    if not np.any(np.isfinite(image)):
        print(f"  No finite values to plot: {save_path}", flush=True)
        return
    fig, ax = plt.subplots(figsize=(10, 6))
    extent = [
        float(time_ms[0]),
        float(time_ms[-1]),
        float(relative_depth[0]),
        float(relative_depth[-1]),
    ]
    displayed = np.ma.masked_invalid(image)
    im = ax.imshow(
        displayed,
        origin="lower",
        aspect="auto",
        extent=extent,
        cmap=cmap,
    )
    im.set_clim(*clim)
    fig.colorbar(im, ax=ax, label=cbar_label)
    if contour is not None and np.any(np.asarray(contour) != 0):
        ax.contour(
            np.abs(np.asarray(contour)),
            levels=[0.5],
            origin="lower",
            extent=extent,
            colors="k",
            linewidths=0.8,
        )
    if time_ms[0] <= 0 <= time_ms[-1]:
        ax.axvline(0.0, color="k", linewidth=0.6)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Relative depth (0 = area start, 1 = area end)")
    ax.set_title(title)
    fig.tight_layout()
    try:
        _ensure_writable(save_path)
        fig.savefig(save_path, dpi=150)
    except PermissionError as exc:
        print(
            f"  Could not write {save_path}: {exc}. Close that image if it "
            "is open, then rerun the session-cache cell.",
            flush=True,
        )
    finally:
        plt.close(fig)


def _write_csv(
    path: str,
    fieldnames: Sequence[str],
    rows: Iterable[Dict[str, object]],
) -> None:
    _ensure_writable(path)
    with open(path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def _safe_path_part(value: str) -> str:
    value = re.sub(r'[<>:"/\\|?*]+', "_", str(value))
    value = re.sub(r"\s+", "_", value.strip())
    return value or "unknown"
