"""Delta phase difference versus distance from the LPAG-center channel.

This is the geometry used by ``extract_delta_phase_vs_lpag_distance.m``:
for every probe that has LPAG, take the saved delta (freq1) phase-lag
column against the LPAG-center site and plot it against signed probe ``y``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np


NPX2_X_COORDS = np.array(
    [[8.0, 40.0], [258.0, 290.0], [508.0, 540.0], [758.0, 790.0]],
    dtype=float,
)


@dataclass
class PhaseDistanceTrace:
    session_name: str
    animal_id: str
    line_label: str
    probe_num: int
    lpag_center_index: int
    lpag_center_y: float
    channel_index: np.ndarray
    distance_um: np.ndarray
    phase_rad: np.ndarray
    phase_unwrapped: np.ndarray
    area: List[str]
    border_distance_um: np.ndarray
    border_from_area: List[str]
    border_to_area: List[str] = field(default_factory=list)


def session_animal_id(session_name: str) -> str:
    tokens = str(session_name).strip().split()
    if not tokens:
        return str(session_name)
    return tokens[0]


def channel_area(areas: Sequence, index: int) -> str:
    if index >= len(areas):
        return ""
    value = areas[index]
    while isinstance(value, (list, tuple)) and value:
        value = value[0]
    if value is None:
        return ""
    label = str(value).strip()
    if label in {"None", "nan", "NaN"}:
        return ""
    return label


def probe_numbers(x: np.ndarray) -> np.ndarray:
    """Assign NPX2 shank numbers, or probe 1 for NPX1."""
    x = np.asarray(x, dtype=float).reshape(-1)
    probe_num = np.full(x.shape, np.nan)
    is_npx2 = False
    for probe_n, pair in enumerate(NPX2_X_COORDS, start=1):
        on_probe = np.isin(x, pair)
        if np.any(on_probe):
            is_npx2 = True
            probe_num[on_probe] = probe_n
    if not is_npx2:
        probe_num[np.isfinite(x)] = 1
    return probe_num


def region_borders(
    distance_um: np.ndarray, areas: Sequence[str]
) -> Tuple[np.ndarray, List[str], List[str]]:
    border_distance = []
    border_from = []
    border_to = []
    distance_um = np.asarray(distance_um, dtype=float).reshape(-1)
    for index in range(1, len(distance_um)):
        from_area = areas[index - 1]
        to_area = areas[index]
        if not from_area or not to_area or from_area == to_area:
            continue
        border_distance.append(0.5 * (distance_um[index - 1] + distance_um[index]))
        border_from.append(from_area)
        border_to.append(to_area)
    return (
        np.asarray(border_distance, dtype=float),
        border_from,
        border_to,
    )


def pooled_region_borders(
    traces: Sequence[PhaseDistanceTrace], bin_um: float = 10.0
) -> Tuple[np.ndarray, List[str]]:
    """Collapse nearby session borders so one animal plot can draw them once."""
    buckets: dict = {}
    for trace in traces:
        for distance, from_area, to_area in zip(
            np.asarray(trace.border_distance_um, dtype=float).reshape(-1),
            trace.border_from_area,
            trace.border_to_area,
        ):
            if not np.isfinite(distance):
                continue
            key = int(np.round(distance / bin_um) * bin_um)
            buckets.setdefault(key, []).append((from_area, to_area))
    distances = []
    labels = []
    for key in sorted(buckets):
        pairs = buckets[key]
        counts = {}
        for pair in pairs:
            counts[pair] = counts.get(pair, 0) + 1
        from_area, to_area = max(counts, key=counts.get)
        distances.append(float(key))
        labels.append(f"{from_area} | {to_area}")
    return np.asarray(distances, dtype=float), labels


def sort_channels_by_probe_depth(
    probe_idx: np.ndarray, x: np.ndarray, y: np.ndarray
) -> np.ndarray:
    """Anatomical order along the shank: increasing y, then x at a given depth.

    LFP rows follow ``chanMap+1``, not probe depth. The Estimate heatmaps
    sort ``channel_map(:,2)`` for the same reason.
    """
    probe_idx = np.asarray(probe_idx, dtype=int).reshape(-1)
    if probe_idx.size == 0:
        return probe_idx
    coords = np.column_stack((y[probe_idx], x[probe_idx]))
    order = np.lexsort((coords[:, 1], coords[:, 0]))
    return probe_idx[order]


def extract_session_traces(
    phase_matrix: np.ndarray,
    channel_map: np.ndarray,
    areas: Sequence,
    session_name: str,
) -> List[PhaseDistanceTrace]:
    """Return one trace per probe that has an LPAG center."""
    phase_matrix = np.asarray(phase_matrix, dtype=float)
    channel_map = np.asarray(channel_map, dtype=float)
    n_chan = min(
        phase_matrix.shape[0],
        phase_matrix.shape[1],
        channel_map.shape[0],
        len(areas),
    )
    if n_chan < 2:
        return []

    x = channel_map[:n_chan, 0]
    y = channel_map[:n_chan, 1]
    valid = np.isfinite(x) & np.isfinite(y)
    probe_num = probe_numbers(x)
    area_names = [channel_area(areas, index) for index in range(n_chan)]
    is_lpag = np.array(
        [
            area_names[index] == "LPAG" and bool(valid[index])
            for index in range(n_chan)
        ],
        dtype=bool,
    )
    if not np.any(is_lpag):
        return []

    lpag_probes = np.unique(probe_num[is_lpag & np.isfinite(probe_num)])
    traces: List[PhaseDistanceTrace] = []
    animal_id = session_animal_id(session_name)
    for probe in lpag_probes:
        probe_mask = valid & (probe_num == probe)
        lpag_mask = probe_mask & is_lpag
        if not np.any(lpag_mask):
            continue
        mean_x = float(np.mean(x[lpag_mask]))
        mean_y = float(np.mean(y[lpag_mask]))
        distance_to_mean = np.abs(x - mean_x) + np.abs(y - mean_y)
        distance_to_mean = distance_to_mean.astype(float)
        distance_to_mean[~lpag_mask] = np.inf
        ref = int(np.argmin(distance_to_mean))
        if not np.isfinite(distance_to_mean[ref]):
            continue

        probe_idx = sort_channels_by_probe_depth(np.where(probe_mask)[0], x, y)
        keep = np.isfinite(phase_matrix[probe_idx, ref])
        probe_idx = sort_channels_by_probe_depth(probe_idx[keep], x, y)
        if probe_idx.size < 2:
            continue

        phase_rad = np.angle(np.exp(1j * phase_matrix[probe_idx, ref]))
        distance_um = y[probe_idx] - y[ref]
        area = [area_names[index] for index in probe_idx]
        border_distance, border_from, border_to = region_borders(distance_um, area)
        probe_int = int(probe)
        if lpag_probes.size > 1:
            line_label = f"{session_name} probe {probe_int}"
        else:
            line_label = session_name
        traces.append(
            PhaseDistanceTrace(
                session_name=session_name,
                animal_id=animal_id,
                line_label=line_label,
                probe_num=probe_int,
                lpag_center_index=ref,
                lpag_center_y=float(y[ref]),
                channel_index=probe_idx.astype(int),
                distance_um=np.asarray(distance_um, dtype=float),
                phase_rad=np.asarray(phase_rad, dtype=float),
                phase_unwrapped=np.unwrap(phase_rad),
                area=area,
                border_distance_um=border_distance,
                border_from_area=border_from,
                border_to_area=border_to,
            )
        )
    return traces


def extract_delta_phase_vs_lpag_distance(
    sessions: Iterable,
    animal_names: Optional[Sequence[str]] = None,
) -> List[PhaseDistanceTrace]:
    """Loop saved session dicts / objects that expose the GENERATE fields."""
    sessions = list(sessions)
    traces: List[PhaseDistanceTrace] = []
    for session_n, session in enumerate(sessions):
        if animal_names is not None and session_n < len(animal_names):
            session_name = str(animal_names[session_n])
        else:
            session_name = _field(session, "session_name", f"session {session_n + 1}")
        phase = _field(session, "mean_phaseLag_freq1")
        channel_map = _field(session, "channel_map")
        areas = _field(session, "areas_by_channel")
        if phase is None or channel_map is None or areas is None:
            continue
        traces.extend(
            extract_session_traces(phase, channel_map, areas, session_name)
        )
    return traces


def lpag_center_map_index(x: np.ndarray, y: np.ndarray, in_lpag: np.ndarray) -> int:
    """Closest LPAG site to the mean (x, y) of LPAG sites (Fig 2 PSTH geometry)."""
    x = np.asarray(x, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    in_lpag = np.asarray(in_lpag, dtype=bool).reshape(-1)
    if not np.any(in_lpag):
        raise ValueError("No LPAG sites")
    mean_loc = np.array([np.mean(x[in_lpag]), np.mean(y[in_lpag])])
    distance = np.abs(x - mean_loc[0]) + np.abs(y - mean_loc[1])
    distance = distance.astype(float)
    distance[~in_lpag] = np.inf
    return int(np.argmin(distance))


def lpag_ref_jobs(
    npx_type: int,
    lpag_rows: Sequence[Tuple[int, float, float]],
    xcoords: np.ndarray,
    ycoords: np.ndarray,
    chan_map: np.ndarray,
    n_lfp: int = 384,
) -> List[dict]:
    """Mirror GENERATE_COHERENCE_MAPS_LPAG_REF local_lpag_jobs (0-based Python).

    ``lpag_rows`` are (probe_num, depth_start, depth_end) with 1-based probe
    numbers. ``chan_map`` holds 0-based channel ids; LFP row = chan_map + 1.
    ``ref_in_shank`` is a 0-based index into the returned shank vectors.
    """
    xcoords = np.asarray(xcoords, dtype=float).reshape(-1)
    ycoords = np.asarray(ycoords, dtype=float).reshape(-1)
    chan_map = np.asarray(chan_map, dtype=float).reshape(-1)
    jobs: List[dict] = []
    if npx_type == 1:
        depth_start = min(row[1] for row in lpag_rows)
        depth_end = max(row[2] for row in lpag_rows)
        in_lpag = (ycoords >= depth_start) & (ycoords <= depth_end)
        if not np.any(in_lpag):
            return []
        closest = lpag_center_map_index(xcoords, ycoords, in_lpag)
        valid = np.isfinite(xcoords) & np.isfinite(ycoords) & np.isfinite(chan_map)
        lfp_rows = chan_map[valid] + 1
        keep = (lfp_rows >= 1) & (lfp_rows <= n_lfp)
        shank_lfp = lfp_rows[keep]
        det_row = chan_map[closest] + 1
        matches = np.where(shank_lfp == det_row)[0]
        if matches.size == 0:
            return []
        jobs.append(
            {
                "probe_num": 1,
                "lfp_rows": shank_lfp,
                "x_coords": xcoords[valid][keep],
                "y_coords": ycoords[valid][keep],
                "detection_lfp_row": det_row,
                "ref_in_shank": int(matches[0]),
            }
        )
        return jobs

    probes = sorted({int(row[0]) for row in lpag_rows})
    for probe_n in probes:
        if probe_n < 1 or probe_n > NPX2_X_COORDS.shape[0]:
            continue
        rows = [row for row in lpag_rows if int(row[0]) == probe_n]
        depth_start = min(row[1] for row in rows)
        depth_end = max(row[2] for row in rows)
        shank_mask = np.isin(xcoords, NPX2_X_COORDS[probe_n - 1])
        lpag_mask = shank_mask & (ycoords >= depth_start) & (ycoords <= depth_end)
        if not np.any(lpag_mask):
            continue
        closest = lpag_center_map_index(xcoords, ycoords, lpag_mask)
        if not shank_mask[closest]:
            shank_idx = np.where(shank_mask)[0]
            closest = int(
                shank_idx[
                    np.argmin(
                        np.abs(xcoords[shank_idx] - np.mean(xcoords[lpag_mask]))
                        + np.abs(ycoords[shank_idx] - np.mean(ycoords[lpag_mask]))
                    )
                ]
            )
        lfp_rows = chan_map[shank_mask] + 1
        x_sh = xcoords[shank_mask]
        y_sh = ycoords[shank_mask]
        keep = np.isfinite(lfp_rows) & (lfp_rows >= 1) & (lfp_rows <= n_lfp)
        lfp_rows = lfp_rows[keep]
        x_sh = x_sh[keep]
        y_sh = y_sh[keep]
        det_row = chan_map[closest] + 1
        matches = np.where(lfp_rows == det_row)[0]
        if matches.size == 0:
            continue
        jobs.append(
            {
                "probe_num": probe_n,
                "lfp_rows": lfp_rows,
                "x_coords": x_sh,
                "y_coords": y_sh,
                "detection_lfp_row": det_row,
                "ref_in_shank": int(matches[0]),
            }
        )
    return jobs


def pli_wpli_phase_vs_ref(
    analytic: np.ndarray, ref_ch: int
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """One-segment PLI / wPLI / phase lag versus a reference channel.

    ``analytic`` is nCh x nTime. This is the inner loop of
    ``Compute_pli_wpli_segments`` when ``ReferenceChannel`` is set: nCh
    comparisons instead of nCh x nCh pairwise maps.
    """
    analytic = np.asarray(analytic, dtype=complex)
    n_ch = analytic.shape[0]
    if ref_ch < 0 or ref_ch >= n_ch:
        raise ValueError(f"ReferenceChannel {ref_ch} is outside 0:{n_ch - 1}")
    xref = analytic[ref_ch]
    pli = np.empty(n_ch)
    wpli = np.empty(n_ch)
    phase_lag = np.empty(n_ch)
    for index in range(n_ch):
        cs = analytic[index] * np.conj(xref)
        im_part = np.imag(cs)
        pli[index] = np.abs(np.mean(np.sign(im_part)))
        wpli[index] = np.abs(np.mean(im_part)) / (np.mean(np.abs(im_part)) + np.finfo(float).eps)
        phase_lag[index] = np.angle(np.mean(cs))
    return pli, wpli, phase_lag


def extract_lpag_ref_session_trace(
    phase_vs_ref: np.ndarray,
    channel_map: np.ndarray,
    areas: Sequence,
    session_name: str,
    lpag_center_index: int,
    probe_num: int = 1,
    npx_type: int = 1,
) -> Optional[PhaseDistanceTrace]:
    """One probe of GENERATE_COHERENCE_MAPS_LPAG_REF vector output.

    ``lpag_center_index`` is 0-based into the shank vectors. Pairwise
    matrices are rejected.
    """
    phase_vs_ref = np.asarray(phase_vs_ref, dtype=float).reshape(-1)
    channel_map = np.asarray(channel_map, dtype=float)
    if channel_map.ndim != 2 or channel_map.shape[1] < 2:
        return None
    n_chan = min(phase_vs_ref.size, channel_map.shape[0], len(areas))
    if n_chan < 2:
        return None
    x = channel_map[:n_chan, 0]
    y = channel_map[:n_chan, 1]
    phase = phase_vs_ref[:n_chan]
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(phase)
    if lpag_center_index < 0 or lpag_center_index >= n_chan or not valid[lpag_center_index]:
        return None
    probe_idx = sort_channels_by_probe_depth(np.where(valid)[0], x, y)
    if probe_idx.size < 2:
        return None
    area_names = [channel_area(areas, index) for index in range(n_chan)]
    area = [area_names[index] for index in probe_idx]
    phase_rad = np.angle(np.exp(1j * phase[probe_idx]))
    distance_um = y[probe_idx] - y[lpag_center_index]
    border_distance, border_from, border_to = region_borders(distance_um, area)
    if npx_type == 2:
        line_label = f"{session_name} probe {int(probe_num)}"
    else:
        line_label = session_name
    return PhaseDistanceTrace(
        session_name=session_name,
        animal_id=session_animal_id(session_name),
        line_label=line_label,
        probe_num=int(probe_num),
        lpag_center_index=int(lpag_center_index),
        lpag_center_y=float(y[lpag_center_index]),
        channel_index=probe_idx.astype(int),
        distance_um=np.asarray(distance_um, dtype=float),
        phase_rad=np.asarray(phase_rad, dtype=float),
        phase_unwrapped=np.unwrap(phase_rad),
        area=area,
        border_distance_um=border_distance,
        border_from_area=border_from,
        border_to_area=border_to,
    )


def session_trend_outlier_mask(
    distance_um: np.ndarray,
    phase_rad: np.ndarray,
    threshold_factor: float = 3.0,
) -> np.ndarray:
    """True where a point is far from a per-session moving-median trend.

    Matches Analyze_coherence_maps: sort by distance, estimate a local trend,
    then flag residuals beyond ``threshold_factor`` scaled MAD (MATLAB
    ``isoutlier(..., 'median')``).
    """
    distance_um = np.asarray(distance_um, dtype=float).reshape(-1)
    phase_rad = np.asarray(phase_rad, dtype=float).reshape(-1)
    n = distance_um.size
    outlier = np.zeros(n, dtype=bool)
    valid = np.isfinite(distance_um) & np.isfinite(phase_rad)
    if np.count_nonzero(valid) < 8:
        return outlier
    order = np.argsort(distance_um[valid], kind="mergesort")
    valid_idx = np.flatnonzero(valid)[order]
    phase_sorted = phase_rad[valid_idx]
    win = max(7, 2 * (phase_sorted.size // 8) + 1)
    if win % 2 == 0:
        win += 1
    pad = win // 2
    padded = np.pad(phase_sorted, (pad, pad), mode="edge")
    trend = np.empty_like(phase_sorted)
    for i in range(phase_sorted.size):
        trend[i] = np.median(padded[i : i + win])
    residual = phase_sorted - trend
    med = np.median(residual)
    mad = np.median(np.abs(residual - med))
    scaled = 1.4826 * mad
    threshold = max(threshold_factor * scaled, 0.5)
    outlier[valid_idx] = np.abs(residual - med) > threshold
    return outlier


def interpolate_session_trend_outliers(
    distance_um: np.ndarray,
    phase_rad: np.ndarray,
    threshold_factor: float = 3.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Replace trend outliers by linear interpolation along distance."""
    distance_um = np.asarray(distance_um, dtype=float).reshape(-1)
    phase_clean = np.asarray(phase_rad, dtype=float).reshape(-1).copy()
    outlier = session_trend_outlier_mask(distance_um, phase_clean, threshold_factor)
    inlier = ~outlier & np.isfinite(distance_um) & np.isfinite(phase_clean)
    query = outlier & np.isfinite(distance_um)
    if np.count_nonzero(inlier) < 2 or not np.any(query):
        return phase_clean, outlier
    order = np.argsort(distance_um[inlier], kind="mergesort")
    filled = np.interp(
        distance_um[query],
        distance_um[inlier][order],
        phase_clean[inlier][order],
        left=np.nan,
        right=np.nan,
    )
    good = np.isfinite(filled)
    query_idx = np.flatnonzero(query)
    phase_clean[query_idx[good]] = filled[good]
    return phase_clean, outlier


def _field(session, name: str, default=None):
    if isinstance(session, dict):
        return session.get(name, default)
    return getattr(session, name, default)
