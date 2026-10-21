"""Pick the newest combined coherence .mat, ignoring satellite pieces.

pairwise (default)
    Newest *merged* / *all_together* / *combined* file, or the all-animal
    2026 save ``coherence_pli_struct_2026.mat``. Ignores satellite pieces
    (mpfc, B1D11007, NPX2) and ``*lpag_ref*`` files.

lpag_ref
    ``GENERATE_COHERENCE_MAPS_LPAG_REF`` vector saves such as
    ``coherence_pli_struct_lpag_ref_2026.mat``.
"""

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
SPECIAL_PARTNERS = {
    "coherence_pli_struct_2026.mat": "coherence_pli_structure_animal_names2026.mat",
    "coherence_pli_struct_lpag_ref_2026.mat": (
        "coherence_pli_struct_lpag_ref_animal_names_2026.mat"
    ),
}


def is_merged_structure_name(filename: str) -> bool:
    """True for a pairwise combined save (not LPAG-ref, not satellite)."""
    return is_structure_name(filename, kind="pairwise")


def is_lpag_ref_structure_name(filename: str) -> bool:
    return is_structure_name(filename, kind="lpag_ref")


def is_structure_name(filename: str, kind: str = "pairwise") -> bool:
    name = str(filename).lower()
    if not name.endswith(".mat"):
        return False
    if "animal_names" in name:
        return False
    if "coherence" not in name and "pli_struct" not in name:
        return False
    is_lpag_ref = "lpag_ref" in name
    if kind == "lpag_ref":
        return is_lpag_ref
    if is_lpag_ref:
        return False
    if any(token in name for token in SATELLITE_TOKENS):
        return False
    if any(token in name for token in MERGED_TOKENS):
        return True
    return "2026" in name and "pli_struct" in name


def partner_names_file(structure_name: str) -> str:
    special = SPECIAL_PARTNERS.get(structure_name)
    if special is not None:
        return special
    special = SPECIAL_PARTNERS.get(str(structure_name).lower())
    if special is not None:
        return special
    if "coherence_structure" in structure_name:
        return structure_name.replace(
            "coherence_structure", "coherence_structure_animal_names", 1
        )
    if structure_name.lower().endswith(".mat"):
        return structure_name[:-4] + "_animal_names.mat"
    return structure_name + "_animal_names"


def find_latest_merged_coherence_files(
    entries: Sequence[Tuple[str, float]],
    kind: str = "pairwise",
) -> Tuple[str, Optional[str]]:
    """Return (structure_name, names_name) from (name, mtime) entries."""
    if kind not in {"pairwise", "lpag_ref"}:
        raise ValueError("kind must be 'pairwise' or 'lpag_ref'")
    candidates = [
        (name, mtime)
        for name, mtime in entries
        if is_structure_name(name, kind=kind)
    ]
    if not candidates:
        label = "LPAG-ref" if kind == "lpag_ref" else "merged / all-together / 2026"
        raise FileNotFoundError(
            f"No {label} coherence file among "
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
        if _is_names_fallback(name, kind)
    ]
    if not name_candidates:
        return structure_name, None
    fallback, _ = max(name_candidates, key=lambda item: item[1])
    return structure_name, fallback


def _is_names_fallback(filename: str, kind: str) -> bool:
    name = str(filename).lower()
    if not name.endswith(".mat") or "animal_names" not in name:
        return False
    is_lpag_ref = "lpag_ref" in name
    if kind == "lpag_ref":
        return is_lpag_ref
    if is_lpag_ref:
        return False
    return any(token in name for token in MERGED_TOKENS) or "2026" in name
