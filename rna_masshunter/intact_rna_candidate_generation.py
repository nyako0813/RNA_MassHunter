"""Pure generation of unmodified intact-RNA CCA/terminal-state candidates.

Assumption count is deterministic: one point each for unknown sample CCA state,
non-confirmed 5-prime state, non-confirmed 3-prime state, and non-confirmed
topology. Observed data are not used and therefore never contribute to the count.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from hashlib import sha256

from rna_masshunter.cca_tail_state import (
    CCATailState,
    CCATailStatus,
    CCATailVariant,
    RegisteredSequenceCCAMode,
    generate_cca_tail_variants,
)
from rna_masshunter.intact_rna_average_mass import (
    TheoreticalMassDefinition,
    calculate_intact_rna_average_mass,
)
from rna_masshunter.intact_rna_mass import (
    CcaPolicy,
    FivePrimeState,
    IntactRnaMassParameters,
    RnaTopology,
    ThreePrimeState,
)
from rna_masshunter.sciex_sample_manifest import (
    AnalyteLevel,
    ExperimentType,
    SCIEXSampleManifest,
    SequenceStatus,
    get_measurement,
    get_rna_identity,
    get_sample,
)
from rna_masshunter.cca_tail_state import (
    CCATailState,
    CCATailStatus,   # ← 追加
    CCATailVariant,
    RegisteredSequenceCCAMode,
    generate_cca_tail_variants,
)

CANDIDATE_GENERATION_APPLIED_TO_FORMAL_SCORE = False
CANDIDATE_GENERATION_APPLIED_TO_RANKING = False
CANDIDATE_GENERATION_APPLIED_TO_CANDIDATE_FILTERING = False
CANDIDATE_GENERATION_APPLIED_TO_FINAL_CONSENSUS = False

_ALLOWED_CONFIG_STATUSES = (CCATailStatus.ASSUMED, CCATailStatus.UNKNOWN)


def _cca_candidate_states_from_config(
    cca_tail_config: dict | None,
    mode: RegisteredSequenceCCAMode,
) -> tuple[CCATailState, ...] | None:
    """Resolve candidate CCA states from config, or None to use built-in defaults."""
    settings = cca_tail_config or {}
    if not settings.get("enabled", True):
        return None
    key = (
        "excludes_cca_candidate_states"
        if mode is RegisteredSequenceCCAMode.EXCLUDES_CCA
        else "includes_cca_candidate_states"
    )
    raw = settings.get(key)
    if not raw:
        return None
    return tuple(CCATailState(value) for value in raw)


class CandidateGenerationError(ValueError):
    """Raised when manifest metadata cannot safely route candidate generation."""


class CandidateCategory(str, Enum):
    UNMODIFIED_INTACT_RNA_CCA_TERMINAL_STATE = "UNMODIFIED_INTACT_RNA_CCA_TERMINAL_STATE"


class CandidateRole(str, Enum):
    UNMODIFIED_REFERENCE_BASELINE = "UNMODIFIED_REFERENCE_BASELINE"


class CandidateSetName(str, Enum):
    PRIMARY_DEFAULT = "PRIMARY_DEFAULT"
    SECONDARY_TERMINAL_DIAGNOSTIC = "SECONDARY_TERMINAL_DIAGNOSTIC"


class AssumptionStatus(str, Enum):
    CONFIRMED = "CONFIRMED"
    ASSUMED = "ASSUMED"
    ALTERNATIVE = "ALTERNATIVE"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class TerminalCandidateSet:
    candidate_set_name: CandidateSetName
    candidate_priority: int
    five_prime_state: FivePrimeState
    five_prime_state_status: AssumptionStatus
    three_prime_state: ThreePrimeState
    three_prime_state_status: AssumptionStatus
    topology: RnaTopology
    topology_status: AssumptionStatus


PRIMARY_DEFAULT = TerminalCandidateSet(
    candidate_set_name=CandidateSetName.PRIMARY_DEFAULT,
    candidate_priority=1,
    five_prime_state=FivePrimeState.MONOPHOSPHATE,
    five_prime_state_status=AssumptionStatus.ASSUMED,
    three_prime_state=ThreePrimeState.OH,
    three_prime_state_status=AssumptionStatus.ASSUMED,
    topology=RnaTopology.LINEAR,
    topology_status=AssumptionStatus.ASSUMED,
)
SECONDARY_TERMINAL_DIAGNOSTIC = TerminalCandidateSet(
    candidate_set_name=CandidateSetName.SECONDARY_TERMINAL_DIAGNOSTIC,
    candidate_priority=2,
    five_prime_state=FivePrimeState.OH,
    five_prime_state_status=AssumptionStatus.ALTERNATIVE,
    three_prime_state=ThreePrimeState.OH,
    three_prime_state_status=AssumptionStatus.ASSUMED,
    topology=RnaTopology.LINEAR,
    topology_status=AssumptionStatus.ASSUMED,
)


@dataclass(frozen=True)
class IntactRnaTheoreticalCandidate:
    candidate_id: str
    candidate_category: CandidateCategory
    rna_identity_id: str
    sample_id: str
    measurement_id: str
    registered_sequence: str
    registered_sequence_sha256: str
    registered_sequence_cca_mode: RegisteredSequenceCCAMode
    core_sequence: str
    core_sequence_sha256: str
    cca_tail_state: CCATailState
    cca_tail_status: CCATailStatus
    cca_tail_sequence: str
    cca_completion_step_count: int
    complete_sequence: str
    complete_sequence_length: int
    complete_sequence_sha256: str
    five_prime_state: FivePrimeState
    five_prime_state_status: AssumptionStatus
    three_prime_state: ThreePrimeState
    three_prime_state_status: AssumptionStatus
    topology: RnaTopology
    topology_status: AssumptionStatus
    terminal_state_confirmed: bool
    theoretical_formula: str
    theoretical_mass: float
    theoretical_mass_type: str
    theoretical_monoisotopic_neutral_mass: float
    theoretical_average_neutral_molecular_mass_m: float
    theoretical_average_m_plus_h: float
    theoretical_average_m_minus_h: float
    primary_reference_candidate: TheoreticalMassDefinition
    observed_output_species: str
    observed_output_species_confirmed: bool
    candidate_role: CandidateRole
    native_modifications_expected: bool
    modification_mass_not_yet_applied: bool
    biological_unmodified_state_assigned: bool
    target_rna_identity_confirmed_by_mass: bool
    co_captured_rna_excluded: bool
    candidate_priority: int
    candidate_set_name: CandidateSetName
    candidate_assumption_count: int
    mass_equivalent_group_id: str
    candidate_ambiguity_count: int
    mass_equivalent_candidate_ids: tuple[str, ...]
    mass_match_only: bool
    unmodified_candidate: bool
    cca_state_confirmed: bool
    structure_identity_assigned: bool
    position_assigned: bool
    modification_assigned: bool
    biological_cause_assigned: bool
    rnase_t_assigned: bool
    applied_to_formal_score: bool
    applied_to_ranking: bool
    applied_to_candidate_filtering: bool
    applied_to_final_consensus: bool
    charge_correction_applied: bool
    proton_correction_applied: bool
    mz_conversion_applied: bool
    average_mass_used: bool
    interpretation_warnings: tuple[str, ...]


def _terminal_label_five(state: FivePrimeState) -> str:
    return {
        FivePrimeState.OH: "5OH",
        FivePrimeState.MONOPHOSPHATE: "5P",
        FivePrimeState.DIPHOSPHATE: "5PP",
        FivePrimeState.TRIPHOSPHATE: "5PPP",
    }[state]


def _terminal_label_three(state: ThreePrimeState) -> str:
    return {
        ThreePrimeState.OH: "3OH",
        ThreePrimeState.MONOPHOSPHATE: "3P",
        ThreePrimeState.CYCLIC_PHOSPHATE: "3CYCLIC_P",
    }[state]


def _candidate_id(
    rna_identity_id: str,
    variant: CCATailVariant,
    terminal_set: TerminalCandidateSet,
) -> str:
    return "__".join((
        rna_identity_id,
        variant.candidate_cca_tail_state.value,
        _terminal_label_five(terminal_set.five_prime_state),
        _terminal_label_three(terminal_set.three_prime_state),
        terminal_set.topology.value,
    ))


def _assumption_count(terminal_set: TerminalCandidateSet) -> int:
    statuses = (
        terminal_set.five_prime_state_status,
        terminal_set.three_prime_state_status,
        terminal_set.topology_status,
    )
    return 1 + sum(status is not AssumptionStatus.CONFIRMED for status in statuses)


def build_unmodified_intact_candidate(
    *,
    rna_identity_id: str,
    sample_id: str,
    measurement_id: str,
    cca_variant: CCATailVariant,
    terminal_candidate_set: TerminalCandidateSet,
) -> IntactRnaTheoreticalCandidate:
    """Combine immutable CCA and terminal assumptions into one mass-only candidate."""
    if not isinstance(cca_variant, CCATailVariant):
        raise TypeError("cca_variant must be CCATailVariant")
    if not isinstance(terminal_candidate_set, TerminalCandidateSet):
        raise TypeError("terminal_candidate_set must be TerminalCandidateSet")
    for field_name, value in (
        ("rna_identity_id", rna_identity_id),
        ("sample_id", sample_id),
        ("measurement_id", measurement_id),
    ):
        if not isinstance(value, str) or not value:
            raise ValueError(f"{field_name} must be non-empty text")

    parameters = IntactRnaMassParameters(
        five_prime_state=terminal_candidate_set.five_prime_state,
        three_prime_state=terminal_candidate_set.three_prime_state,
        topology=terminal_candidate_set.topology,
        cca_policy=CcaPolicy.AS_PROVIDED,
        convert_t_to_u=False,
        terminal_state_confirmed=False,
    )
    average_result = calculate_intact_rna_average_mass(
        cca_variant.complete_candidate_sequence,
        parameters=parameters,
    )
    mass_result = average_result.monoisotopic_result
    candidate_id = _candidate_id(rna_identity_id, cca_variant, terminal_candidate_set)
    return IntactRnaTheoreticalCandidate(
        candidate_id=candidate_id,
        candidate_category=CandidateCategory.UNMODIFIED_INTACT_RNA_CCA_TERMINAL_STATE,
        rna_identity_id=rna_identity_id,
        sample_id=sample_id,
        measurement_id=measurement_id,
        registered_sequence=cca_variant.registered_sequence,
        registered_sequence_sha256=cca_variant.registered_sequence_sha256,
        registered_sequence_cca_mode=cca_variant.registered_sequence_cca_mode,
        core_sequence=cca_variant.core_sequence,
        core_sequence_sha256=cca_variant.core_sequence_sha256,
        cca_tail_state=cca_variant.candidate_cca_tail_state,
        cca_tail_status=cca_variant.candidate_cca_tail_status,
        cca_tail_sequence=cca_variant.candidate_cca_tail_sequence,
        cca_completion_step_count=cca_variant.cca_completion_step_count,
        complete_sequence=cca_variant.complete_candidate_sequence,
        complete_sequence_length=cca_variant.complete_candidate_sequence_length,
        complete_sequence_sha256=cca_variant.complete_candidate_sequence_sha256,
        five_prime_state=terminal_candidate_set.five_prime_state,
        five_prime_state_status=terminal_candidate_set.five_prime_state_status,
        three_prime_state=terminal_candidate_set.three_prime_state,
        three_prime_state_status=terminal_candidate_set.three_prime_state_status,
        topology=terminal_candidate_set.topology,
        topology_status=terminal_candidate_set.topology_status,
        terminal_state_confirmed=False,
        theoretical_formula=mass_result.formula,
        theoretical_mass=mass_result.monoisotopic_neutral_mass,
        theoretical_mass_type=mass_result.theoretical_mass_type,
        theoretical_monoisotopic_neutral_mass=mass_result.monoisotopic_neutral_mass,
        theoretical_average_neutral_molecular_mass_m=(
            average_result.average_neutral_molecular_mass_m
        ),
        theoretical_average_m_plus_h=average_result.average_m_plus_h,
        theoretical_average_m_minus_h=average_result.average_m_minus_h,
        primary_reference_candidate=TheoreticalMassDefinition.AVERAGE_NEUTRAL_M,
        observed_output_species="UNKNOWN",
        observed_output_species_confirmed=False,
        candidate_role=CandidateRole.UNMODIFIED_REFERENCE_BASELINE,
        native_modifications_expected=True,
        modification_mass_not_yet_applied=True,
        biological_unmodified_state_assigned=False,
        target_rna_identity_confirmed_by_mass=False,
        co_captured_rna_excluded=False,
        candidate_priority=terminal_candidate_set.candidate_priority,
        candidate_set_name=terminal_candidate_set.candidate_set_name,
        candidate_assumption_count=_assumption_count(terminal_candidate_set),
        mass_equivalent_group_id="",
        candidate_ambiguity_count=1,
        mass_equivalent_candidate_ids=(candidate_id,),
        mass_match_only=True,
        unmodified_candidate=True,
        cca_state_confirmed=False,
        structure_identity_assigned=False,
        position_assigned=False,
        modification_assigned=False,
        biological_cause_assigned=False,
        rnase_t_assigned=False,
        applied_to_formal_score=False,
        applied_to_ranking=False,
        applied_to_candidate_filtering=False,
        applied_to_final_consensus=False,
        charge_correction_applied=False,
        proton_correction_applied=False,
        mz_conversion_applied=False,
        average_mass_used=False,
        interpretation_warnings=(
            "MASS_MATCH_ONLY",
            "UNMODIFIED_REFERENCE_BASELINE",
            "NATIVE_MODIFICATIONS_EXPECTED",
            "OBSERVED_OUTPUT_SPECIES_UNKNOWN",
            "CCA_STATE_NOT_CONFIRMED",
            "TERMINAL_STATE_NOT_CONFIRMED",
        ),
    )


def _assign_mass_equivalent_groups(
    candidates: tuple[IntactRnaTheoreticalCandidate, ...],
) -> tuple[IntactRnaTheoreticalCandidate, ...]:
    by_formula: dict[str, list[IntactRnaTheoreticalCandidate]] = {}
    for candidate in candidates:
        by_formula.setdefault(candidate.theoretical_formula, []).append(candidate)
    result: list[IntactRnaTheoreticalCandidate] = []
    for candidate in candidates:
        group = by_formula[candidate.theoretical_formula]
        group_id = "MEG__" + sha256(
            (candidate.theoretical_mass_type + "|" + candidate.theoretical_formula).encode("ascii")
        ).hexdigest()[:16].upper()
        ids = tuple(item.candidate_id for item in group)
        result.append(replace(
            candidate,
            mass_equivalent_group_id=group_id,
            candidate_ambiguity_count=len(group),
            mass_equivalent_candidate_ids=ids,
        ))
    return tuple(result)


def generate_candidates_for_measurement(
    manifest: SCIEXSampleManifest,
    measurement_id: str,
    *,
    include_secondary_terminal_state: bool = True,
    cca_tail_config: dict | None = None,   # ← 追加
) -> tuple[IntactRnaTheoreticalCandidate, ...]:
    """Route a full-length measurement through manifest references only."""
    if not isinstance(manifest, SCIEXSampleManifest):
        raise TypeError("manifest must be SCIEXSampleManifest")
    if not isinstance(include_secondary_terminal_state, bool):
        raise TypeError("include_secondary_terminal_state must be boolean")
    measurement = get_measurement(manifest, measurement_id)
    if measurement.experiment_type is not ExperimentType.FULL_LENGTH:
        raise CandidateGenerationError(
            f"measurement {measurement_id} is not FULL_LENGTH"
        )
    if measurement.expected_analyte_level is not AnalyteLevel.INTACT_RNA:
        raise CandidateGenerationError(
            f"measurement {measurement_id} does not expect INTACT_RNA"
        )
    sample = get_sample(manifest, measurement.sample_id)
    identity = get_rna_identity(manifest, sample.rna_identity_id)
    if identity.sequence_status is SequenceStatus.UNKNOWN or identity.sequence is None:
        raise CandidateGenerationError(
            f"RNA identity {identity.rna_identity_id} has no candidate-ready sequence"
        )
    if identity.registered_sequence_cca_mode is RegisteredSequenceCCAMode.UNKNOWN:
        raise CandidateGenerationError(
            f"RNA identity {identity.rna_identity_id} has UNKNOWN registered CCA mode"
        )

    candidate_states = _cca_candidate_states_from_config(
        cca_tail_config, identity.registered_sequence_cca_mode
    )
    cca_variants = generate_cca_tail_variants(
        identity.sequence,
        identity.registered_sequence_cca_mode,
        candidate_states,   # ← None なら既存のデフォルト挙動と完全に同じ
    )
    terminal_sets = (
        (PRIMARY_DEFAULT, SECONDARY_TERMINAL_DIAGNOSTIC)
        if include_secondary_terminal_state
        else (PRIMARY_DEFAULT,)
    )
    candidates = tuple(
        build_unmodified_intact_candidate(
            rna_identity_id=identity.rna_identity_id,
            sample_id=sample.sample_id,
            measurement_id=measurement.measurement_id,
            cca_variant=variant,
            terminal_candidate_set=terminal_set,
        )
        for terminal_set in terminal_sets
        for variant in cca_variants
    )
    return _assign_mass_equivalent_groups(candidates)