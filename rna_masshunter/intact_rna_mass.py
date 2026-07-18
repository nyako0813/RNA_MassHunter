"""Explicit elemental-composition model for intact canonical RNA molecules."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from hashlib import sha256
from typing import TypeVar

from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.structure_fragment import RNA_RESIDUE_COMPOSITIONS

THEORETICAL_MASS_TYPE = "MONOISOTOPIC_NEUTRAL"
ASCII_WHITESPACE = frozenset(" \t\n\r\v\f")


class FivePrimeState(str, Enum):
    OH = "OH"
    MONOPHOSPHATE = "MONOPHOSPHATE"
    DIPHOSPHATE = "DIPHOSPHATE"
    TRIPHOSPHATE = "TRIPHOSPHATE"


class ThreePrimeState(str, Enum):
    OH = "OH"
    MONOPHOSPHATE = "MONOPHOSPHATE"
    CYCLIC_PHOSPHATE = "CYCLIC_PHOSPHATE"


class RnaTopology(str, Enum):
    LINEAR = "LINEAR"
    CIRCULAR = "CIRCULAR"


class CcaPolicy(str, Enum):
    AS_PROVIDED = "AS_PROVIDED"


@dataclass(frozen=True)
class IntactRnaMassParameters:
    five_prime_state: FivePrimeState = FivePrimeState.OH
    three_prime_state: ThreePrimeState = ThreePrimeState.OH
    topology: RnaTopology = RnaTopology.LINEAR
    cca_policy: CcaPolicy = CcaPolicy.AS_PROVIDED
    convert_t_to_u: bool = False
    terminal_state_confirmed: bool = False


@dataclass(frozen=True)
class IntactRnaMassResult:
    normalized_sequence: str
    sequence_length: int
    sequence_sha256: str
    ends_with_cca: bool
    phosphodiester_bond_count: int
    five_prime_state: str
    three_prime_state: str
    topology: str
    cca_policy: str
    elemental_composition: ElementalComposition
    formula: str
    monoisotopic_neutral_mass: float
    theoretical_mass_type: str
    terminal_state_confirmed: bool
    t_to_u_conversion_applied: bool
    warnings: tuple[str, ...]


_EnumT = TypeVar("_EnumT", bound=Enum)


def _coerce_enum(value: object, enum_type: type[_EnumT], field_name: str) -> _EnumT:
    if isinstance(value, enum_type):
        return value
    try:
        return enum_type(str(value).upper())
    except (TypeError, ValueError) as exc:
        allowed = ", ".join(item.value for item in enum_type)
        raise ValueError(f"unknown {field_name}: {value!r}; expected one of {allowed}") from exc


def _normalize_sequence(sequence: str, convert_t_to_u: bool) -> tuple[str, bool]:
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    if not isinstance(convert_t_to_u, bool):
        raise ValueError("convert_t_to_u must be boolean")
    normalized = "".join(character for character in sequence if character not in ASCII_WHITESPACE).upper()
    if not normalized:
        raise ValueError("RNA sequence must not be empty")
    conversion_applied = "T" in normalized and convert_t_to_u
    if conversion_applied:
        normalized = normalized.replace("T", "U")
    for position, base in enumerate(normalized, 1):
        if base not in RNA_RESIDUE_COMPOSITIONS:
            raise ValueError(f"invalid canonical RNA base at position {position}: {base!r}")
    return normalized, conversion_applied


def _apply_delta(counts: dict[str, int], delta: dict[str, int]) -> None:
    for element, value in delta.items():
        counts[element] = counts.get(element, 0) + value
        if counts[element] == 0:
            counts.pop(element)


_FIVE_PRIME_DELTAS = {
    FivePrimeState.OH: {},
    FivePrimeState.MONOPHOSPHATE: {"H": 1, "O": 3, "P": 1},
    FivePrimeState.DIPHOSPHATE: {"H": 2, "O": 6, "P": 2},
    FivePrimeState.TRIPHOSPHATE: {"H": 3, "O": 9, "P": 3},
}
_THREE_PRIME_DELTAS = {
    ThreePrimeState.OH: {},
    ThreePrimeState.MONOPHOSPHATE: {"H": 1, "O": 3, "P": 1},
    ThreePrimeState.CYCLIC_PHOSPHATE: {"H": -1, "O": 2, "P": 1},
}
_OH_OH_FROM_RESIDUES = {"H": 1, "O": -2, "P": -1}  # +H2O - HPO3
_CIRCULAR_DELTA = {"H": -2, "O": -1}  # -H2O from the linear OH/OH reference


def calculate_intact_rna_mass(
    sequence: str,
    *,
    parameters: IntactRnaMassParameters,
) -> IntactRnaMassResult:
    """Calculate an exact neutral monoisotopic mass without charge or m/z conversion."""
    if not isinstance(parameters, IntactRnaMassParameters):
        raise TypeError("parameters must be IntactRnaMassParameters")

    five_prime = _coerce_enum(parameters.five_prime_state, FivePrimeState, "five_prime_state")
    three_prime = _coerce_enum(parameters.three_prime_state, ThreePrimeState, "three_prime_state")
    topology = _coerce_enum(parameters.topology, RnaTopology, "topology")
    cca_policy = _coerce_enum(parameters.cca_policy, CcaPolicy, "cca_policy")
    if not isinstance(parameters.terminal_state_confirmed, bool):
        raise ValueError("terminal_state_confirmed must be boolean")

    normalized, conversion_applied = _normalize_sequence(sequence, parameters.convert_t_to_u)
    if topology is RnaTopology.CIRCULAR:
        if five_prime is not FivePrimeState.OH or three_prime is not ThreePrimeState.OH:
            raise ValueError("circular RNA does not have independent 5-prime or 3-prime terminal phosphate states")
        if len(normalized) < 2:
            raise ValueError("circular RNA requires at least two residues")

    counts: dict[str, int] = {}
    for base in normalized:
        _apply_delta(counts, RNA_RESIDUE_COMPOSITIONS[base])
    _apply_delta(counts, _OH_OH_FROM_RESIDUES)
    _apply_delta(counts, _FIVE_PRIME_DELTAS[five_prime])
    _apply_delta(counts, _THREE_PRIME_DELTAS[three_prime])
    if topology is RnaTopology.CIRCULAR:
        _apply_delta(counts, _CIRCULAR_DELTA)
    negative = {element: count for element, count in counts.items() if count < 0}
    if negative:
        raise ValueError(f"invalid final elemental composition: {negative}")

    composition = ElementalComposition(counts)
    warnings: list[str] = []
    if conversion_applied:
        warnings.append("T_TO_U_CONVERSION_APPLIED")
    if not parameters.terminal_state_confirmed:
        warnings.append("TERMINAL_STATE_NOT_CONFIRMED")
    return IntactRnaMassResult(
        normalized_sequence=normalized,
        sequence_length=len(normalized),
        sequence_sha256=sha256(normalized.encode("ascii")).hexdigest(),
        ends_with_cca=normalized.endswith("CCA"),
        phosphodiester_bond_count=(len(normalized) if topology is RnaTopology.CIRCULAR else len(normalized) - 1),
        five_prime_state=five_prime.value,
        three_prime_state=three_prime.value,
        topology=topology.value,
        cca_policy=cca_policy.value,
        elemental_composition=composition,
        formula=composition.canonical_string(),
        monoisotopic_neutral_mass=composition.exact_mass,
        theoretical_mass_type=THEORETICAL_MASS_TYPE,
        terminal_state_confirmed=parameters.terminal_state_confirmed,
        t_to_u_conversion_applied=conversion_applied,
        warnings=tuple(warnings),
    )
