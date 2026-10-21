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

        probe_idx = np.where(probe_mask)[0]
        probe_idx = probe_idx[np.argsort(y[probe_idx])]
        keep = np.isfinite(phase_matrix[probe_idx, ref])
        probe_idx = probe_idx[keep]
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


def _field(session, name: str, default=None):
    if isinstance(session, dict):
        return session.get(name, default)
    return getattr(session, name, default)
