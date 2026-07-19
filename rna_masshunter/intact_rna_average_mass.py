"""Deterministic average-mass references for intact RNA.

The fixed conventional atomic weights are C=12.0107, H=1.00794,
N=14.0067, O=15.9994, P=30.973761998, and S=32.065 Da.  They are
kept in-process so production calculations never depend on a web service.
Average neutral M is calculated only from elemental composition.  Proton,
electron, charge, and m/z adjustments are separate display diagnostics.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType
from typing import Mapping

from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.intact_rna_mass import (
    IntactRnaMassParameters,
    IntactRnaMassResult,
    calculate_intact_rna_mass,
)

AVERAGE_ATOMIC_MASS_CONSTANT_SET = "IUPAC_CONVENTIONAL_ATOMIC_WEIGHTS_FIXED"
AVERAGE_ATOMIC_MASS_CONSTANT_VERSION = "RNA_MASSHUNTER_2026_01"
AVERAGE_ATOMIC_MASS_CONSTANT_SOURCE = (
    "Fixed conventional standard atomic weights; values embedded in module docstring"
)
AVERAGE_ATOMIC_MASSES: Mapping[str, float] = MappingProxyType({
    "C": 12.0107,
    "H": 1.00794,
    "N": 14.0067,
    "O": 15.9994,
    "P": 30.973761998,
    "S": 32.065,
})
PROTON_MASS_DA = 1.007276466621


class TheoreticalMassDefinition(str, Enum):
    MONOISOTOPIC_NEUTRAL_M = "MONOISOTOPIC_NEUTRAL_M"
    AVERAGE_NEUTRAL_M = "AVERAGE_NEUTRAL_M"
    AVERAGE_M_PLUS_H = "AVERAGE_M_PLUS_H"
    AVERAGE_M_MINUS_H = "AVERAGE_M_MINUS_H"


class ComparisonReferenceRole(str, Enum):
    PRIMARY_CANDIDATE_NEUTRAL_M = "PRIMARY_CANDIDATE_NEUTRAL_M"
    OUTPUT_SPECIES_DIAGNOSTIC_M_PLUS_H = "OUTPUT_SPECIES_DIAGNOSTIC_M_PLUS_H"
    OUTPUT_SPECIES_DIAGNOSTIC_M_MINUS_H = "OUTPUT_SPECIES_DIAGNOSTIC_M_MINUS_H"
    MONOISOTOPIC_DIAGNOSTIC_ONLY = "MONOISOTOPIC_DIAGNOSTIC_ONLY"


class MassDisplaySpecies(str, Enum):
    M = "M"
    M_PLUS_H = "M_PLUS_H"
    M_MINUS_H = "M_MINUS_H"


@dataclass(frozen=True)
class MassReference:
    theoretical_mass_definition: TheoreticalMassDefinition
    comparison_role: ComparisonReferenceRole
    mass_display_species: MassDisplaySpecies
    mass: float
    protonation_adjustment_applied: bool
    deprotonation_adjustment_applied: bool
    charge_adjustment_applied: bool = False
    electron_mass_adjustment_applied: bool = False
    mz_conversion_applied: bool = False


@dataclass(frozen=True)
class IntactRnaAverageMassResult:
    monoisotopic_result: IntactRnaMassResult
    average_neutral_molecular_mass_m: float
    average_m_plus_h: float
    average_m_minus_h: float
    average_atomic_mass_constant_set: str
    average_atomic_mass_constant_version: str
    average_atomic_mass_constant_source: str
    proton_mass_da: float
    ion_mode: str | None
    charge_state: int | None
    references: tuple[MassReference, ...]

    @property
    def elemental_composition(self) -> ElementalComposition:
        return self.monoisotopic_result.elemental_composition

    @property
    def monoisotopic_neutral_mass(self) -> float:
        return self.monoisotopic_result.monoisotopic_neutral_mass


def _composition_counts(composition: ElementalComposition | Mapping[str, int]) -> dict[str, int]:
    if isinstance(composition, ElementalComposition):
        return composition.to_dict()
    if not isinstance(composition, Mapping):
        raise TypeError("composition must be ElementalComposition or a mapping")
    counts: dict[str, int] = {}
    for raw_element, raw_count in composition.items():
        element = str(raw_element)
        if element not in AVERAGE_ATOMIC_MASSES:
            raise ValueError(f"unsupported average-mass element: {element}")
        if isinstance(raw_count, bool) or not isinstance(raw_count, int):
            raise ValueError(f"element count must be an integer: {element}={raw_count!r}")
        if raw_count < 0:
            raise ValueError(f"negative element count: {element}={raw_count}")
        counts[element] = raw_count
    return counts


def calculate_average_neutral_mass_from_composition(
    composition: ElementalComposition | Mapping[str, int],
) -> float:
    """Return average neutral molecular mass M; an empty composition returns 0.0."""
    counts = _composition_counts(composition)
    return sum(AVERAGE_ATOMIC_MASSES[element] * count for element, count in counts.items())


def calculate_intact_rna_average_mass(
    sequence: str,
    *,
    parameters: IntactRnaMassParameters,
    ion_mode: str | None = None,
    charge_state: int | None = None,
) -> IntactRnaAverageMassResult:
    """Calculate neutral M; optional acquisition descriptors never alter it."""
    normalized_mode = str(ion_mode).upper() if ion_mode is not None else None
    if normalized_mode not in {None, "POSITIVE", "NEGATIVE"}:
        raise ValueError("ion_mode must be POSITIVE, NEGATIVE, or None")
    if charge_state is not None and (
        isinstance(charge_state, bool) or not isinstance(charge_state, int)
    ):
        raise ValueError("charge_state must be an integer or None")
    mono = calculate_intact_rna_mass(sequence, parameters=parameters)
    neutral = calculate_average_neutral_mass_from_composition(mono.elemental_composition)
    plus_h = neutral + PROTON_MASS_DA
    minus_h = neutral - PROTON_MASS_DA
    references = (
        MassReference(
            TheoreticalMassDefinition.AVERAGE_NEUTRAL_M,
            ComparisonReferenceRole.PRIMARY_CANDIDATE_NEUTRAL_M,
            MassDisplaySpecies.M,
            neutral,
            False,
            False,
        ),
        MassReference(
            TheoreticalMassDefinition.AVERAGE_M_PLUS_H,
            ComparisonReferenceRole.OUTPUT_SPECIES_DIAGNOSTIC_M_PLUS_H,
            MassDisplaySpecies.M_PLUS_H,
            plus_h,
            True,
            False,
        ),
        MassReference(
            TheoreticalMassDefinition.AVERAGE_M_MINUS_H,
            ComparisonReferenceRole.OUTPUT_SPECIES_DIAGNOSTIC_M_MINUS_H,
            MassDisplaySpecies.M_MINUS_H,
            minus_h,
            False,
            True,
        ),
        MassReference(
            TheoreticalMassDefinition.MONOISOTOPIC_NEUTRAL_M,
            ComparisonReferenceRole.MONOISOTOPIC_DIAGNOSTIC_ONLY,
            MassDisplaySpecies.M,
            mono.monoisotopic_neutral_mass,
            False,
            False,
        ),
    )
    return IntactRnaAverageMassResult(
        mono,
        neutral,
        plus_h,
        minus_h,
        AVERAGE_ATOMIC_MASS_CONSTANT_SET,
        AVERAGE_ATOMIC_MASS_CONSTANT_VERSION,
        AVERAGE_ATOMIC_MASS_CONSTANT_SOURCE,
        PROTON_MASS_DA,
        normalized_mode,
        charge_state,
        references,
    )
