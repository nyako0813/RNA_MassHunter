from pathlib import Path
from typing import Any

import yaml

from rna_masshunter.warnings_manager import add_warning

PROTON_MASS = 1.007276466812


def load_base_masses(path: str | Path, warnings: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def calculate_unmodified_rna_mass(sequence: str, base_masses: dict[str, Any], warnings: list[dict[str, Any]] | None = None) -> float | None:
    residues = base_masses.get("rna_residue_masses", {})
    constants = base_masses.get("constants", {})
    water = float(constants.get("water", 18.010564684))
    total = water
    for index, base in enumerate(sequence.upper(), start=1):
        if base not in residues:
            if warnings is not None:
                add_warning(warnings, "ERROR", "masses", "Unknown RNA base in sequence.", {"position": index, "base": base})
            return None
        total += float(residues[base])
    return total


def neutral_mass_from_mz(mz: float, charge: int, polarity: str = "negative", proton_mass: float = PROTON_MASS) -> float:
    z = abs(int(charge))
    if str(polarity).lower() == "negative":
        return float(mz) * z + z * proton_mass
    return float(mz) * z - z * proton_mass
