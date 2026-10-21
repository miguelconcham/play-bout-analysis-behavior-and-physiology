"""Pick the newest merged / all-together coherence .mat, ignoring satellite pieces."""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

MERGED_TOKENS = (
    "merged",
    "all_together",
    "alltogeth",
    "all together",
    "combined",
)
SATELLITE_TOKENS = ("mpfc", "b1d11007", "npx2")


def is_merged_structure_name(filename: str) -> bool:
    name = str(filename).lower()
    if not name.endswith(".mat"):
        return False
    if "animal_names" in name:
        return False
    if any(token in name for token in SATELLITE_TOKENS):
        return False
    if "coherence" not in name and "pli_struct" not in name:
        return False
    return any(token in name for token in MERGED_TOKENS)


def partner_names_file(structure_name: str) -> str:
    if "coherence_structure" in structure_name:
        return structure_name.replace(
            "coherence_structure", "coherence_structure_animal_names", 1
        )
    if structure_name.lower().endswith(".mat"):
        return structure_name[:-4] + "_animal_names.mat"
    return structure_name + "_animal_names"


def find_latest_merged_coherence_files(
    entries: Sequence[Tuple[str, float]],
) -> Tuple[str, Optional[str]]:
    """Return (structure_name, names_name) from (name, mtime) entries."""
    candidates = [
        (name, mtime)
        for name, mtime in entries
        if is_merged_structure_name(name)
    ]
    if not candidates:
        raise FileNotFoundError(
            "No merged / all-together coherence file among "
            + ", ".join(name for name, _ in entries)
        )
    structure_name, _ = max(candidates, key=lambda item: item[1])
    names_name = partner_names_file(structure_name)
    present = {name for name, _ in entries}
    if names_name in present:
        return structure_name, names_name
    name_candidates = [
        (name, mtime)
        for name, mtime in entries
        if "animal_names" in name.lower()
        and any(token in name.lower() for token in MERGED_TOKENS)
    ]
    if not name_candidates:
        return structure_name, None
    fallback, _ = max(name_candidates, key=lambda item: item[1])
    return structure_name, fallback
