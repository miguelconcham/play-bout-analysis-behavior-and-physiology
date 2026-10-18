"""
Cross-session anatomical consistency of saved behavior-aligned kCSD.

This module does not estimate kCSD. It reads ``csd_mean_columns.npy`` maps
that were already averaged across staggered probe columns, assigns those
depths to anatomical areas with the same rules as
``GENERATE_PSTH_AREA_MAPS.m``, and asks whether the saved kCSD sign repeats
across sessions.

Polarity is reported exactly as stored (``RdBu_r`` in the existing figures).
Positive and negative are not relabeled as source or sink.

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
import os
import re
import shutil
import warnings
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np


PathLike = Union[str, os.PathLike]

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
    """Load the area workbook used by the MATLAB generators."""
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

    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
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
    finite = np.isfinite(relative) & np.all(np.isfinite(csd), axis=1)
    if int(np.count_nonzero(finite)) < 2:
        return None
    resampled = _interpolate_rows(relative[finite], csd[finite], relative_depth)
    if resampled is None or not np.any(np.isfinite(resampled)):
        return None
    return resampled


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
    save_plots: bool = True,
    refresh: bool = True,
    stop_on_error: bool = False,
) -> Dict[str, object]:
    """
    Pool saved averaged kCSD maps by anatomical area across sessions.

    Requires ``csd_mean_columns.npy``. Per-column matrices are not refit and
    ``KCSD1D`` is never called. Run ``average_saved_behavior_kcsd_from_root``
    first when a completed kCSD run has not yet written the column average.
    """
    if n_relative_depths < 2:
        raise ValueError("n_relative_depths must be at least 2")
    if min_depths < 2:
        raise ValueError("min_depths must be at least 2")
    root = os.path.abspath(os.fspath(root))
    rows = _coerce_area_rows(area_rows, area_table_path)
    rows_by_animal: Dict[str, List[AreaLimit]] = {}
    for row in rows:
        rows_by_animal.setdefault(row.animal_name, []).append(row)

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
    session_time: Dict[Tuple[str, str, str, str, str], np.ndarray] = {}
    skipped: List[Dict[str, str]] = []

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
        try:
            assignment = _segment_leaf(
                leaf,
                parsed,
                rows_by_animal,
                relative_depth,
                min_depths=min_depths,
            )
        except Exception as exc:
            if stop_on_error:
                raise
            print(f"  Skipped: {exc}", flush=True)
            skipped.append({"leaf": leaf, "reason": str(exc)})
            continue
        if assignment["status"] != "ok":
            print(f"  Skipped: {assignment['reason']}", flush=True)
            skipped.append({"leaf": leaf, "reason": assignment["reason"]})
            continue
        for area, matrix in assignment["areas"].items():
            key = (
                parsed["scope"],
                parsed["behavior"],
                parsed["alignment"],
                parsed["session"],
                area,
            )
            if key in session_time and not np.allclose(
                session_time[key], assignment["time_ms"]
            ):
                reason = f"Time axis mismatch for {label} area {area}"
                if stop_on_error:
                    raise ValueError(reason)
                print(f"  Skipped area {area}: {reason}", flush=True)
                skipped.append({"leaf": leaf, "reason": reason})
                continue
            session_parts.setdefault(key, []).append(matrix)
            session_time[key] = assignment["time_ms"]

    grouped: Dict[Tuple[str, str, str, str], List[Tuple[str, np.ndarray, np.ndarray]]] = {}
    for key, matrices in session_parts.items():
        scope, behavior, alignment, session, area = key
        session_mean = _nanmean_stack(np.stack(matrices, axis=0))
        if not np.any(np.isfinite(session_mean)):
            continue
        group_key = (scope, behavior, alignment, area)
        grouped.setdefault(group_key, []).append(
            (session, session_mean, session_time[key])
        )

    output_root = os.path.join(root, output_folder_name)
    if refresh and os.path.isdir(output_root):
        shutil.rmtree(output_root)
    os.makedirs(output_root, exist_ok=True)

    summaries = []
    stable_rows: List[Dict[str, object]] = []
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
        print(
            f"Processing {scope}/{behavior} {alignment} {area}: "
            f"{len(usable)} sessions",
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
        with open(os.path.join(folder, "sessions.txt"), "w", encoding="utf-8") as stream:
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
        summaries.append(summary)
        stable_rows.extend(feature_rows)
        plot_groups.setdefault((scope, behavior, area), {})[alignment] = {
            "folder": folder,
            "mean": stats["mean"],
            "consistency": stats["sign_consistency"],
            "n_sessions": stats["n_sessions"],
            "n_sessions_total": len(usable),
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
        ],
        summaries,
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
        f"{len(stable_rows)} stable bins, index={summary_csv}",
        flush=True,
    )
    return {
        "output_root": output_root,
        "summary_csv": summary_csv,
        "stable_features_csv": stable_csv,
        "n_groups": len(summaries),
        "n_missing_mean": len(discovered["missing_means"]),
        "missing_mean_folders": discovered["missing_means"],
        "skipped": skipped,
        "groups": summaries,
        "color_limits": color_limits,
    }


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
    for identity, indices in grouped.items():
        if len(indices) < min_depths:
            continue
        chosen = spans[identity]
        resampled = resample_area_segment(
            y[indices],
            csd[indices],
            chosen.depth_start,
            chosen.depth_end,
            relative_depth,
        )
        if resampled is None:
            continue
        area_pieces.setdefault(chosen.area, []).append(resampled)
    for area, pieces in area_pieces.items():
        areas[area] = _nanmean_stack(np.stack(pieces, axis=0))
    if not areas:
        return {
            "status": "skip",
            "reason": f"No area had at least {min_depths} depths in {parsed['session']}",
        }
    return {"status": "ok", "areas": areas, "time_ms": time_ms}


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
    np.save(os.path.join(folder, "mean_csd.npy"), stats["mean"])
    np.save(os.path.join(folder, "sem_csd.npy"), stats["sem"])
    np.save(os.path.join(folder, "fraction_positive.npy"), stats["fraction_positive"])
    np.save(os.path.join(folder, "fraction_negative.npy"), stats["fraction_negative"])
    np.save(os.path.join(folder, "sign_consistency.npy"), stats["sign_consistency"])
    np.save(
        os.path.join(folder, "n_sessions_by_depth.npy"),
        stats["n_sessions_by_depth"],
    )
    np.save(os.path.join(folder, "n_sessions.npy"), stats["n_sessions"])
    np.save(os.path.join(folder, "relative_depth.npy"), relative_depth)
    np.save(os.path.join(folder, "time_ms.npy"), time_ms)


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
) -> None:
    if not np.any(np.isfinite(image)):
        print(f"  No finite values to plot: {save_path}", flush=True)
        return
    fig, ax = plt.subplots(figsize=(10, 6))
    displayed = np.ma.masked_invalid(image)
    im = ax.imshow(
        displayed,
        origin="lower",
        aspect="auto",
        extent=[
            float(time_ms[0]),
            float(time_ms[-1]),
            float(relative_depth[0]),
            float(relative_depth[-1]),
        ],
        cmap=cmap,
    )
    im.set_clim(*clim)
    fig.colorbar(im, ax=ax, label=cbar_label)
    if time_ms[0] <= 0 <= time_ms[-1]:
        ax.axvline(0.0, color="k", linewidth=0.6)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Relative depth (0 = area start, 1 = area end)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    plt.close(fig)


def _write_csv(
    path: str,
    fieldnames: Sequence[str],
    rows: Iterable[Dict[str, object]],
) -> None:
    with open(path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def _safe_path_part(value: str) -> str:
    value = re.sub(r'[<>:"/\\|?*]+', "_", str(value))
    value = re.sub(r"\s+", "_", value.strip())
    return value or "unknown"
