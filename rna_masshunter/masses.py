from pathlib import Path
from typing import Any

import yaml

from rna_masshunter.warnings_manager import add_warning

PROTON_MASS = 1.007276466812
DEFAULT_WATER_MASS = 18.010564684
DEFAULT_PHOSPHATE_MASS = 79.966331


def load_base_masses(path: str | Path, warnings: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if warnings is not None:
        if "rna_residue_masses" not in data:
            add_warning(warnings, "ERROR", "masses", "base_masses.yaml is missing rna_residue_masses.")
        if "constants" not in data:
            add_warning(warnings, "WARNING", "masses", "base_masses.yaml is missing constants; defaults will be used where possible.")
        elif "phosphate" not in data.get("constants", {}):
            add_warning(warnings, "WARNING", "masses", "base_masses.yaml constants.phosphate is missing; default phosphate mass was used.", DEFAULT_PHOSPHATE_MASS)
    return data


def _constant(
    base_masses: dict[str, Any],
    name: str,
    default: float,
    warnings: list[dict[str, Any]] | None = None,
) -> float:
    constants = base_masses.get("constants", {})
    if name not in constants and warnings is not None:
        add_warning(warnings, "WARNING", "masses", f"base_masses.yaml constants.{name} is missing; default value was used.", default)
    try:
        return float(constants.get(name, default))
    except (TypeError, ValueError):
        if warnings is not None:
            add_warning(warnings, "WARNING", "masses", f"base_masses.yaml constants.{name} is not numeric; default value was used.", default)
        return default


def calculate_unmodified_rna_mass(
    sequence: str,
    base_masses: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
    terminal_form: str = "default",
) -> float | None:
    residues = base_masses.get("rna_residue_masses", {})
    if not residues:
        if warnings is not None:
            add_warning(warnings, "ERROR", "masses", "RNA residue masses are missing; mass could not be calculated.")
        return None

    water = _constant(base_masses, "water", DEFAULT_WATER_MASS, warnings)
    phosphate = _constant(base_masses, "phosphate", DEFAULT_PHOSPHATE_MASS, warnings)
    total = water
    for index, base in enumerate(sequence.upper(), start=1):
        if base not in residues:
            if warnings is not None:
                add_warning(warnings, "ERROR", "masses", "Unknown RNA base in sequence.", {"position": index, "base": base})
            return None
        total += float(residues[base])

    terminal_adjustments = {
        "default": 0.0,
        "dephosphorylated": 0.0,
        "residual_phosphate": phosphate,
        "cyclic_phosphate": phosphate - water,
    }
    if terminal_form not in terminal_adjustments:
        if warnings is not None:
            add_warning(warnings, "WARNING", "masses", "Unknown terminal form; no terminal mass adjustment was applied.", terminal_form)
        terminal_form = "default"
    return total + terminal_adjustments[terminal_form]


def neutral_mass_from_mz(mz: float, charge: int, polarity: str = "negative", proton_mass: float = PROTON_MASS) -> float:
    z = abs(int(charge))
    if str(polarity).lower() == "negative":
        return float(mz) * z + z * proton_mass
    return float(mz) * z - z * proton_mass


def mz_from_neutral_mass(neutral_mass: float, charge: int, polarity: str = "negative", proton_mass: float = PROTON_MASS) -> float:
    z = abs(int(charge))
    if str(polarity).lower() == "negative":
        return (float(neutral_mass) - z * proton_mass) / z
    return (float(neutral_mass) + z * proton_mass) / z
