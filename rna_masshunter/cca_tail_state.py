"""Pure, non-formal candidate model for explicit tRNA 3-prime CCA tail states."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from hashlib import sha256
from typing import Iterable, TypeVar

from rna_masshunter.intact_rna_mass import (
    ASCII_WHITESPACE,
    IntactRnaMassParameters,
    IntactRnaMassResult,
    calculate_intact_rna_mass,
)
from rna_masshunter.structure_fragment import RNA_RESIDUE_COMPOSITIONS

CCA_TAIL_STATE_APPLIED_TO_FORMAL_SCORE = False
CCA_TAIL_STATE_APPLIED_TO_RANKING = False
CCA_TAIL_STATE_APPLIED_TO_CANDIDATE_FILTERING = False
CCA_TAIL_STATE_APPLIED_TO_FINAL_CONSENSUS = False


class RegisteredSequenceCCAMode(str, Enum):
    EXCLUDES_CCA = "EXCLUDES_CCA"
    INCLUDES_COMPLETE_CCA = "INCLUDES_COMPLETE_CCA"
    UNKNOWN = "UNKNOWN"


class CCATailState(str, Enum):
    NONE = "NONE"
    C = "C"
    CC = "CC"
    CCA = "CCA"


class CCATailStatus(str, Enum):
    CONFIRMED = "CONFIRMED"
    ASSUMED = "ASSUMED"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class CCATailVariant:
    registered_sequence: str
    registered_sequence_length: int
    registered_sequence_sha256: str
    registered_sequence_cca_mode: RegisteredSequenceCCAMode
    registered_cca_tail_state: CCATailState
    core_sequence: str
    core_sequence_length: int
    core_sequence_sha256: str
    candidate_cca_tail_state: CCATailState
    candidate_cca_tail_sequence: str
    candidate_cca_tail_length: int
    candidate_cca_tail_status: CCATailStatus
    complete_candidate_sequence: str
    complete_candidate_sequence_length: int
    complete_candidate_sequence_sha256: str
    missing_cca_suffix: str
    required_added_nucleotides: tuple[str, ...]
    cca_completion_required: bool
    cca_completion_step_count: int
    cca_state_confirmed: bool
    repair_intermediate_assigned: bool
    biological_cause_assigned: bool
    rnase_t_assigned: bool
    structure_identity_assigned: bool
    mass_match_only: bool
    formal_evidence: bool
    interpretation_warnings: tuple[str, ...]


@dataclass(frozen=True)
class CCATailVariantMass:
    variant: CCATailVariant
    intact_mass_parameters: IntactRnaMassParameters
    intact_mass_result: IntactRnaMassResult
    monoisotopic_neutral_mass: float
    charge_adjustment_applied: bool
    mz_conversion_applied: bool
    formal_evidence: bool


_EnumT = TypeVar("_EnumT", bound=Enum)
_TAIL_SEQUENCE = {
    CCATailState.NONE: "",
    CCATailState.C: "C",
    CCATailState.CC: "CC",
    CCATailState.CCA: "CCA",
}
_MISSING_SUFFIX = {
    CCATailState.NONE: "CCA",
    CCATailState.C: "CA",
    CCATailState.CC: "A",
    CCATailState.CCA: "",
}
_STATE_ORDER = (
    CCATailState.NONE,
    CCATailState.C,
    CCATailState.CC,
    CCATailState.CCA,
)
_DEFAULT_EXCLUDES_STATES = _STATE_ORDER
_DEFAULT_INCLUDES_STATES = (CCATailState.CCA, CCATailState.CC)


def _coerce_enum(value: object, enum_type: type[_EnumT], field_name: str) -> _EnumT:
    if isinstance(value, enum_type):
        return value
    try:
        return enum_type(value)
    except (TypeError, ValueError) as exc:
        allowed = ", ".join(item.value for item in enum_type)
        raise ValueError(f"unknown {field_name}: {value!r}; expected one of {allowed}") from exc


def _normalize_rna_sequence(sequence: str) -> str:
    if not isinstance(sequence, str):
        raise TypeError("registered_sequence must be a string")
    normalized = "".join(base for base in sequence if base not in ASCII_WHITESPACE).upper()
    if not normalized:
        raise ValueError("registered_sequence must not be empty")
    for position, base in enumerate(normalized, 1):
        if base not in RNA_RESIDUE_COMPOSITIONS:
            raise ValueError(f"invalid canonical RNA base at position {position}: {base!r}")
    return normalized


def derive_core_sequence(
    registered_sequence: str,
    registered_sequence_cca_mode: RegisteredSequenceCCAMode | str,
) -> str:
    """Normalize and derive an explicit core without inferring mode from suffix."""
    normalized = _normalize_rna_sequence(registered_sequence)
    mode = _coerce_enum(
        registered_sequence_cca_mode,
        RegisteredSequenceCCAMode,
        "registered_sequence_cca_mode",
    )
    if mode is RegisteredSequenceCCAMode.UNKNOWN:
        raise ValueError("UNKNOWN registered CCA mode has no explicit core boundary")
    if mode is RegisteredSequenceCCAMode.EXCLUDES_CCA:
        return normalized
    if not normalized.endswith("CCA"):
        raise ValueError("INCLUDES_COMPLETE_CCA registered sequence must end with CCA")
    core = normalized[:-3]
    if not core:
        raise ValueError("INCLUDES_COMPLETE_CCA registered sequence must contain a non-empty core")
    return core


def build_cca_tail_variant(
    registered_sequence: str,
    registered_sequence_cca_mode: RegisteredSequenceCCAMode | str,
    candidate_cca_tail_state: CCATailState | str,
    *,
    cca_tail_status: CCATailStatus | str = CCATailStatus.ASSUMED,
) -> CCATailVariant:
    """Build one mass-only CCA-tail hypothesis without biological assignment."""
    normalized = _normalize_rna_sequence(registered_sequence)
    mode = _coerce_enum(
        registered_sequence_cca_mode,
        RegisteredSequenceCCAMode,
        "registered_sequence_cca_mode",
    )
    state = _coerce_enum(candidate_cca_tail_state, CCATailState, "candidate_cca_tail_state")
    status = _coerce_enum(cca_tail_status, CCATailStatus, "cca_tail_status")
    core = derive_core_sequence(normalized, mode)
    registered_state = (
        CCATailState.NONE
        if mode is RegisteredSequenceCCAMode.EXCLUDES_CCA
        else CCATailState.CCA
    )
    tail = _TAIL_SEQUENCE[state]
    complete = core + tail
    missing = _MISSING_SUFFIX[state]
    additions = tuple(missing)
    return CCATailVariant(
        registered_sequence=normalized,
        registered_sequence_length=len(normalized),
        registered_sequence_sha256=sha256(normalized.encode("ascii")).hexdigest(),
        registered_sequence_cca_mode=mode,
        registered_cca_tail_state=registered_state,
        core_sequence=core,
        core_sequence_length=len(core),
        core_sequence_sha256=sha256(core.encode("ascii")).hexdigest(),
        candidate_cca_tail_state=state,
        candidate_cca_tail_sequence=tail,
        candidate_cca_tail_length=len(tail),
        candidate_cca_tail_status=status,
        complete_candidate_sequence=complete,
        complete_candidate_sequence_length=len(complete),
        complete_candidate_sequence_sha256=sha256(complete.encode("ascii")).hexdigest(),
        missing_cca_suffix=missing,
        required_added_nucleotides=additions,
        cca_completion_required=bool(missing),
        cca_completion_step_count=len(additions),
        cca_state_confirmed=False,
        repair_intermediate_assigned=False,
        biological_cause_assigned=False,
        rnase_t_assigned=False,
        structure_identity_assigned=False,
        mass_match_only=True,
        formal_evidence=False,
        interpretation_warnings=(
            "MASS_MATCH_ONLY",
            "CCA_STATE_NOT_CONFIRMED",
            "NO_BIOLOGICAL_CAUSE_ASSIGNED",
        ),
    )


def generate_cca_tail_variants(
    registered_sequence: str,
    registered_sequence_cca_mode: RegisteredSequenceCCAMode | str,
    candidate_states: Iterable[CCATailState | str] | None = None,
) -> tuple[CCATailVariant, ...]:
    """Generate a deterministic mode-specific set of CCA-tail hypotheses."""
    mode = _coerce_enum(
        registered_sequence_cca_mode,
        RegisteredSequenceCCAMode,
        "registered_sequence_cca_mode",
    )
    if mode is RegisteredSequenceCCAMode.UNKNOWN:
        raise ValueError("UNKNOWN registered CCA mode cannot generate candidates")
    if candidate_states is None:
        states = (
            _DEFAULT_EXCLUDES_STATES
            if mode is RegisteredSequenceCCAMode.EXCLUDES_CCA
            else _DEFAULT_INCLUDES_STATES
        )
    else:
        if isinstance(candidate_states, (str, bytes)):
            raise TypeError("candidate_states must be an iterable of CCATailState values")
        requested = {
            _coerce_enum(value, CCATailState, "candidate_cca_tail_state")
            for value in candidate_states
        }
        states = tuple(state for state in _STATE_ORDER if state in requested)
        if not states:
            raise ValueError("candidate_states must not be empty")
    return tuple(
        build_cca_tail_variant(registered_sequence, mode, state)
        for state in states
    )


def calculate_cca_tail_variant_mass(
    variant: CCATailVariant,
    intact_mass_parameters: IntactRnaMassParameters,
) -> CCATailVariantMass:
    """Calculate neutral intact mass while keeping terminal chemistry explicit."""
    if not isinstance(variant, CCATailVariant):
        raise TypeError("variant must be CCATailVariant")
    if not isinstance(intact_mass_parameters, IntactRnaMassParameters):
        raise TypeError("intact_mass_parameters must be IntactRnaMassParameters")
    result = calculate_intact_rna_mass(
        variant.complete_candidate_sequence,
        parameters=intact_mass_parameters,
    )
    return CCATailVariantMass(
        variant=variant,
        intact_mass_parameters=intact_mass_parameters,
        intact_mass_result=result,
        monoisotopic_neutral_mass=result.monoisotopic_neutral_mass,
        charge_adjustment_applied=False,
        mz_conversion_applied=False,
        formal_evidence=False,
    )
