"""
Within-animal 1D kCSD for delta peak and trough LFP exports.

Peak and trough stay separate: they are different events and are not
subtracted from each other. Sessions of one implanted animal and one probe
are combined the same way as play onset. Each session's evoked LFP is
weighted by its event count, which is the mean you get by stacking the
trials. kCSD is fit only after that pooled profile exists.

The mean-across-columns figure draws a white dashed line at each brain-region
boundary from ``Area_limits_GoodLooking.xlsx``. The workbook is only read.
"""

from __future__ import annotations

import csv
import json
import os
import re
import shutil
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from area_csd_consistency import (
    AreaLimit,
    animal_key_from_session,
    load_area_limits,
)
from behavior_aligned_kcsd import (
    _save_averaged_column_outputs,
    _write_pooled_leaf,
    average_kcsd_columns,
    trial_weighted_psth,
)
from csd_1d_pipeline import (
    DELTA_CSD_DEFAULTS,
    DatasetConfig,
    discover_delta_csd_datasets,
    epoch_and_average,
    load_raw,
    plot_column,
    prepare_channels,
    run_dataset,
)


ANIMAL_POOLED_FOLDER = "animal_pooled"
_PROBE_FOLDER = re.compile(r"probe(\d+)\Z", re.IGNORECASE)


def region_border_depths(
    rows: Sequence[AreaLimit],
    *,
    y_min: float,
    y_max: float,
    npx_type: int,
    probe_num: int,
) -> np.ndarray:
    """
    Unique region edges, in micrometers, that fall inside the plotted depths.

    NPX1 keeps ``PAG`` probe areas. NPX2 keeps the matching probe number.
    A shared edge between two regions is returned once.
    """
    if y_max < y_min:
        raise ValueError(f"Depth range is reversed: {y_min} .. {y_max}")
    edges: List[float] = []
    for row in rows:
        if row.depth_end < row.depth_start:
            continue
        if npx_type == 1:
            if row.probe_area != "PAG":
                continue
        elif npx_type == 2:
            if row.probe_num != probe_num:
                continue
        else:
            raise ValueError(f"Unsupported NPX type: {npx_type}")
        for edge in (float(row.depth_start), float(row.depth_end)):
            if y_min <= edge <= y_max:
                edges.append(edge)
    if not edges:
        return np.array([], dtype=float)
    return np.unique(np.asarray(edges, dtype=float))


def probe_identity(probe_label: str) -> Tuple[int, int]:
    """Return ``(npx_type, probe_num)`` for an empty NPX1 label or ``probeN``."""
    if probe_label == "" or probe_label.upper() == "NPX1":
        return 1, 1
    match = _PROBE_FOLDER.fullmatch(probe_label)
    if match is None:
        raise ValueError(f"Unrecognized probe folder: {probe_label}")
    return 2, int(match.group(1))


def parse_delta_leaf(root: str, folder: str) -> Optional[Dict[str, str]]:
    """
    Read a peak or trough session leaf.

    Accepted layouts are ``peak/<session>``, ``peak/<session>/probeN``, and
    the older ``<session>/peak`` or ``<session>/probeN/peak`` trees.
    """
    relative = os.path.relpath(folder, root).replace("\\", "/")
    parts = [part for part in relative.split("/") if part and part != "."]
    if not parts or ANIMAL_POOLED_FOLDER in parts:
        return None
    if parts[0] in {"peak", "trough"}:
        alignment = parts[0]
        rest = parts[1:]
    elif parts[-1] in {"peak", "trough"}:
        alignment = parts[-1]
        rest = parts[:-1]
    else:
        return None
    probe = ""
    if len(rest) == 2 and _PROBE_FOLDER.fullmatch(rest[-1]):
        session, probe = rest
    elif len(rest) == 1:
        session = rest[0]
    else:
        return None
    return {
        "alignment": alignment,
        "session": session,
        "probe": probe,
        "folder": folder,
    }


def pool_delta_evoked(records: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Trial-weighted mean of one evoked profile per session.

    This uses the play-onset average: identical electrodes are averaged
    directly, and differing depths are interpolated only inside the overlap.
    """
    adapted = []
    for record in records:
        evoked = np.asarray(record["evoked"], dtype=float)
        adapted.append(
            {
                "session": record["session"],
                "n_events": int(record["n_events"]),
                "x": np.asarray(record["x"], dtype=float),
                "y": np.asarray(record["y"], dtype=float),
                "time_ms": np.asarray(record["time_ms"], dtype=float),
                "onset_evoked": evoked,
                "offset_evoked": evoked,
            }
        )
    pooled = trial_weighted_psth(adapted)
    return {
        "x": pooled["x"],
        "y": pooled["y"],
        "evoked": pooled["onset_evoked"],
        "time_ms": pooled["time_ms"],
        "n_events": int(pooled["n_events"]),
        "aligned_by": pooled["aligned_by"],
    }


def _delta_session_evoked(
    folder: str,
    *,
    epoch_settings: Dict[str, Any],
) -> Dict[str, Any]:
    lfp, x, y, t0, fs_loaded = load_raw(folder)
    fs = float(fs_loaded) if fs_loaded is not None else float(epoch_settings["fs"])
    cfg = DatasetConfig(
        folder=folder,
        fs=fs,
        ref_chan=epoch_settings["ref_chan"],
        channels_remove=tuple(epoch_settings["channels_remove"]),
        chan_start=epoch_settings["chan_start"],
        chan_end=epoch_settings["chan_end"],
        baseline_correct=bool(epoch_settings["baseline_correct"]),
        normalize=str(epoch_settings["normalize"]),
        t_pre=float(epoch_settings["t_pre"]),
        t_post=float(epoch_settings["t_post"]),
        compute_kcsd=False,
        save=False,
        save_plots=False,
        plot=False,
    )
    lfp, x, y, _meta = prepare_channels(lfp, x, y, cfg)
    evoked, time_ms, _time_indices, epochs = epoch_and_average(
        lfp,
        t0,
        fs,
        float(epoch_settings["t_pre"]),
        float(epoch_settings["t_post"]),
        float(epoch_settings["step_ms"]),
        baseline_correct=bool(epoch_settings["baseline_correct"]),
        normalize=str(epoch_settings["normalize"]),
    )
    del epochs
    return {
        "x": np.asarray(x, dtype=float),
        "y": np.asarray(y, dtype=float),
        "evoked": evoked,
        "time_ms": time_ms,
        "n_events": int(lfp.shape[0]),
        "fs": fs,
    }


def _rows_by_animal(rows: Sequence[AreaLimit]) -> Dict[str, List[AreaLimit]]:
    grouped: Dict[str, List[AreaLimit]] = {}
    for row in rows:
        grouped.setdefault(row.animal_name, []).append(row)
    return grouped


def _plot_pooled_mean(
    *,
    name: str,
    alignment: str,
    folder: str,
    output_subdir: str,
    averaged: Dict[str, np.ndarray],
    time_ms: np.ndarray,
    time_indices: np.ndarray,
    borders: np.ndarray,
    show: bool,
    save_plot: bool,
) -> None:
    limit = float(np.nanmax(np.abs(averaged["csd_mean"])))
    clim = (-limit, limit) if np.isfinite(limit) and limit > 0 else None
    save_path = None
    if save_plot:
        figure_dir = os.path.join(folder, output_subdir, "figures")
        os.makedirs(figure_dir, exist_ok=True)
        save_path = os.path.join(figure_dir, "kcsd_xmean.png")
    plot_column(
        averaged["csd_mean"],
        time_ms,
        time_indices,
        averaged["y"],
        title=f"{name} — {alignment} — mean across probe columns",
        cbar_label="mean kCSD (a.u.)",
        clim=clim,
        horizontal_lines=[float(depth) for depth in borders],
        save_path=save_path,
        show=show,
    )
    if save_path is not None:
        print(f"    Saved averaged figure: {save_path}", flush=True)


def estimate_animal_pooled_delta_kcsd(
    root: str,
    area_table: str,
    *,
    output_folder_name: str = ANIMAL_POOLED_FOLDER,
    alignments: Sequence[str] = ("peak", "trough"),
    stop_on_error: bool = False,
    refresh: bool = True,
    **overrides: Any,
) -> List[Dict[str, Any]]:
    """
    Pool peak and trough trials within each animal, then estimate kCSD.

    ``area_table`` is ``Area_limits_GoodLooking.xlsx``. It is read so the
    mean figure can mark region borders. It is never saved or overwritten.
    """
    root = os.path.abspath(os.fspath(root))
    if not os.path.isdir(root):
        raise FileNotFoundError(f"Delta CSD root not found: {root}")
    area_rows = load_area_limits(area_table)
    rows_by_animal = _rows_by_animal(area_rows)

    grouped: Dict[Tuple[str, str, str], List[Dict[str, str]]] = {}
    for dataset in discover_delta_csd_datasets(root, alignments=alignments):
        parsed = parse_delta_leaf(root, dataset["folder"])
        if parsed is None:
            print(f"Skipping unrecognized delta folder: {dataset['folder']}", flush=True)
            continue
        try:
            animal_key = animal_key_from_session(parsed["session"])
        except ValueError as exc:
            print(f"Skipping {dataset['folder']}: {exc}", flush=True)
            continue
        grouped.setdefault(
            (parsed["alignment"], animal_key, parsed["probe"]), []
        ).append(parsed)
    if not grouped:
        raise FileNotFoundError(
            f"No peak/trough session folders with lfp.npy under {root}."
        )

    settings = dict(DELTA_CSD_DEFAULTS)
    settings.update(overrides)
    save_mean_plot = bool(settings.get("save_plots", True))
    show_plot = bool(settings.get("plot", False))
    compute_kcsd = bool(settings.get("compute_kcsd", True))
    output_subdir = str(settings.get("output_subdir", "csd_1d"))
    epoch_settings = {
        "fs": float(settings.get("fs", 2500.0)),
        "t_pre": float(settings["t_pre"]),
        "t_post": float(settings["t_post"]),
        "step_ms": float(settings["step_ms"]),
        "baseline_correct": bool(settings["baseline_correct"]),
        "normalize": str(settings["normalize"]),
        "ref_chan": settings.get("ref_chan"),
        "channels_remove": tuple(settings.get("channels_remove") or ()),
        "chan_start": settings.get("chan_start"),
        "chan_end": settings.get("chan_end"),
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
        f"Pooling trials within animal for {len(grouped)} peak/trough probe groups, "
        "then estimating kCSD",
        flush=True,
    )
    results: List[Dict[str, Any]] = []
    index_rows: List[Dict[str, Any]] = []
    for group_number, key in enumerate(sorted(grouped), start=1):
        alignment, animal_key, probe = key
        members = sorted(grouped[key], key=lambda item: item["session"])
        label = f"{alignment} {animal_key} {probe or 'NPX1'}"
        print(
            f"\n[{group_number}/{len(grouped)}] {label}: {len(members)} sessions",
            flush=True,
        )
        try:
            psths = []
            for member in members:
                print(f"  Reading trials: {member['session']}", flush=True)
                psth = _delta_session_evoked(
                    member["folder"],
                    epoch_settings=epoch_settings,
                )
                psth["session"] = member["session"]
                psths.append(psth)
                print(f"    {psth['n_events']} events", flush=True)
            reference_fs = psths[0]["fs"]
            if any(not np.isclose(item["fs"], reference_fs) for item in psths[1:]):
                raise ValueError(f"Sampling rate differs across sessions in {label}")
            pooled = pool_delta_evoked(psths)
            tail = [animal_key] if not probe else [animal_key, probe]
            leaf = os.path.join(output_root, alignment, *tail)
            _write_pooled_leaf(
                leaf,
                pooled["evoked"],
                pooled["x"],
                pooled["y"],
                pooled["time_ms"],
                reference_fs,
            )
            source_rows = [
                {"session": item["session"], "n_events": int(item["n_events"])}
                for item in psths
            ]
            pooling_record: Dict[str, Any] = {
                "animal_key": animal_key,
                "alignment": alignment,
                "probe": probe or "NPX1",
                "n_sessions": len(psths),
                "n_events": int(pooled["n_events"]),
                "aligned_by": pooled["aligned_by"],
                "sessions": source_rows,
                "region_borders_um": [],
            }
            result: Dict[str, Any] = {
                "name": f"{output_folder_name}/{alignment}/{animal_key}"
                + (f"/{probe}" if probe else ""),
                "animal_pooling": pooling_record,
                "folder": leaf,
            }
            if compute_kcsd:
                print(
                    f"  Pooled {pooled['n_events']} events by {pooled['aligned_by']}; "
                    "estimating kCSD",
                    flush=True,
                )
                kcsd_result = run_dataset(
                    {
                        "folder": leaf,
                        "name": result["name"],
                        "fs": reference_fs,
                        "t_pre": epoch_settings["t_pre"],
                        "t_post": epoch_settings["t_post"],
                        "step_ms": float(settings["step_ms"]),
                        "baseline_correct": epoch_settings["baseline_correct"],
                        "normalize": epoch_settings["normalize"],
                        "ref_chan": epoch_settings["ref_chan"],
                        "channels_remove": epoch_settings["channels_remove"],
                        "chan_start": epoch_settings["chan_start"],
                        "chan_end": epoch_settings["chan_end"],
                        "compute_kcsd": True,
                        "compute_fd": bool(settings.get("compute_fd", False)),
                        "y_scale": settings["y_scale"],
                        "h": settings["h"],
                        "n_src_factor": settings["n_src_factor"],
                        "smooth_sigma_time": settings["smooth_sigma_time"],
                        "sort_by_depth": settings["sort_by_depth"],
                        "save": bool(settings.get("save", True)),
                        "save_plots": False,
                        "plot": False,
                        "output_subdir": output_subdir,
                    }
                )
                averaged = average_kcsd_columns(
                    kcsd_result["csd_kcsd"],
                    kcsd_result["x"],
                    kcsd_result["y"],
                    kcsd_result["evoked"],
                )
                if settings.get("save", True):
                    _save_averaged_column_outputs(
                        leaf,
                        averaged,
                        kcsd_result["time_ms"],
                        kcsd_result["time_indices"],
                        output_subdir=output_subdir,
                    )
                animal_rows = rows_by_animal.get(animal_key, [])
                if not animal_rows:
                    print(
                        f"  No area-table rows for {animal_key}; "
                        "the mean figure has no region borders.",
                        flush=True,
                    )
                npx_type, probe_num = probe_identity(probe)
                borders = region_border_depths(
                    animal_rows,
                    y_min=float(np.min(averaged["y"])),
                    y_max=float(np.max(averaged["y"])),
                    npx_type=npx_type,
                    probe_num=probe_num,
                )
                pooling_record["region_borders_um"] = [float(depth) for depth in borders]
                _plot_pooled_mean(
                    name=result["name"],
                    alignment=alignment,
                    folder=leaf,
                    output_subdir=output_subdir,
                    averaged=averaged,
                    time_ms=kcsd_result["time_ms"],
                    time_indices=kcsd_result["time_indices"],
                    borders=borders,
                    show=show_plot,
                    save_plot=save_mean_plot and bool(settings.get("save", True)),
                )
                result["kcsd"] = kcsd_result
                result["region_borders_um"] = pooling_record["region_borders_um"]
            with open(os.path.join(leaf, "animal_pooling.json"), "w", encoding="utf-8") as stream:
                json.dump(pooling_record, stream, indent=2)
            results.append(result)
            index_rows.append(
                {
                    "alignment": alignment,
                    "animal_key": animal_key,
                    "probe": probe or "NPX1",
                    "n_sessions": len(psths),
                    "n_events": int(pooled["n_events"]),
                    "sessions": "; ".join(
                        f"{item['session']} ({item['n_events']})" for item in source_rows
                    ),
                    "aligned_by": pooled["aligned_by"],
                    "folder": os.path.relpath(leaf, root),
                    "region_borders_um": "; ".join(
                        f"{depth:g}" for depth in pooling_record["region_borders_um"]
                    ),
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
                "alignment",
                "animal_key",
                "probe",
                "n_sessions",
                "n_events",
                "sessions",
                "aligned_by",
                "folder",
                "region_borders_um",
            ],
        )
        writer.writeheader()
        writer.writerows(index_rows)
    n_failed = sum("error" in result for result in results)
    print(
        f"\nAnimal-pooled delta kCSD finished: {len(results) - n_failed} successful, "
        f"{n_failed} failed. Index: {index_path}",
        flush=True,
    )
    return results
