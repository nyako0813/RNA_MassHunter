"""MS1-only shadow audit of RNase-T1 fragment neutral-delta state series.

The audit deliberately reports mass-pattern candidates, not chemical identities.
It reuses the project T1 digest, ion generator, mzML profile builder, and chemical
mass registry.  Nothing returned here is eligible for formal propagation.
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import asdict, dataclass, field
from enum import Enum
from hashlib import sha256
from math import isfinite
from pathlib import Path
from statistics import median
from typing import Any, Mapping, Sequence

from rna_masshunter.sciex_intact_peak_family import DeltaMassDefinition
from rna_masshunter.sciex_mzml_source_metadata_audit import (
    MzMLSourceMetadataRecord, PolarityStatus, RepresentationStatus,
)
from rna_masshunter.sciex_t1_fragment_delta_audit import (
    ChemicalReferenceCategory, build_chemical_delta_reference_registry,
)
from rna_masshunter.sciex_t1_fragment_shadow_match import (
    T1FragmentMatchClass, T1FragmentMatchParameters, T1IonCandidate, T1IonMode,
    TheoreticalT1Fragment, generate_t1_ion_candidates,
    generate_theoretical_t1_fragments,
)
from rna_masshunter.sciex_t1_profile_peak_audit import T1PeakQualityClass
from rna_masshunter.sciex_t1_replicate_consistency_audit import (
    ReplicateRunPeak, ReplicateRunPeakProfile, build_ms1_peak_profile_for_run,
)

OPTIONAL_RESULT_KEY = "sciex_t1_fragment_state_series_audit"
ALGORITHM_VERSION = "sciex-t1-fragment-state-series-audit-v1"

_BLOCK_ORDER = (
    "INPUT_FILE_NOT_FOUND", "INPUT_FILE_UNREADABLE", "SOURCE_METADATA_RECORD_MISSING",
    "USER_MANIFEST_CONTEXT_MISSING", "SEQUENCE_MISSING", "T1_FRAGMENT_GENERATION_FAILED",
    "NO_THEORETICAL_FRAGMENTS", "NO_MS1_SPECTRA", "PROFILE_EXTRACTION_FAILED",
    "NO_DETECTED_PEAKS", "SOURCE_POLARITY_NOT_NEGATIVE", "MIXED_POLARITY_INPUT",
    "MISSING_POLARITY_METADATA", "REPRESENTATION_NOT_PROFILE",
    "MISSING_REPRESENTATION_METADATA", "NO_FRAGMENT_MATCHES", "NO_STATE_SERIES",
    "LOW_SCAN_RECURRENCE", "LOW_PROMINENCE", "INVALID_FWHM",
    "APEX_CENTROID_DISAGREEMENT", "FRAGMENT_AMBIGUITY", "CHARGE_AMBIGUITY",
    "PEAK_MULTIPLICITY", "STATE_ASSIGNMENT_AMBIGUITY", "INCOMPLETE_STATE_SERIES",
    "INSUFFICIENT_STATE_SPACING_ACCURACY", "FULL_LENGTH_SERIES_RESULT_MISSING",
    "AMBIGUOUS_FRAGMENT_LOCALIZATION", "CHEMICAL_IDENTITY_UNSUPPORTED",
)


def _blocks(values: Sequence[str]) -> tuple[str, ...]:
    found = set(values)
    return tuple(x for x in _BLOCK_ORDER if x in found) + tuple(sorted(found - set(_BLOCK_ORDER)))


def _id(prefix: str, value: str) -> str:
    return prefix + "__" + sha256(value.encode()).hexdigest()[:20].upper()


@dataclass(frozen=True, kw_only=True)
class StateSeriesSafeguards:
    shadow_analysis_only: bool = True
    mass_evidence_only: bool = True
    formal_propagation: bool = False
    chemical_identity_assigned: bool = False
    modification_assigned: bool = False
    exact_nucleotide_localization: bool = False
    exact_atom_localization: bool = False
    reaction_order_assigned: bool = False
    ms2_used_for_state_assignment: bool = False
    applied_to_formal_score: bool = False
    applied_to_ranking: bool = False
    applied_to_candidate_filtering: bool = False
    applied_to_final_consensus: bool = False


class StateLabel(str, Enum):
    BASE_STATE = "BASE_STATE"
    PLUS_16_EQUIVALENT = "PLUS_16_EQUIVALENT"
    PLUS_18_EQUIVALENT = "PLUS_18_EQUIVALENT"
    PLUS_32_EQUIVALENT = "PLUS_32_EQUIVALENT"
    PLUS_34_EQUIVALENT = "PLUS_34_EQUIVALENT"


class StateMatchStatus(str, Enum):
    STRICT = "STRICT"
    SUPPORTIVE = "SUPPORTIVE"


class SeriesAmbiguityStatus(str, Enum):
    UNAMBIGUOUS = "UNAMBIGUOUS"
    FRAGMENT_AMBIGUOUS = "FRAGMENT_AMBIGUOUS"
    CHARGE_AMBIGUOUS = "CHARGE_AMBIGUOUS"
    PEAK_MULTIPLICITY = "PEAK_MULTIPLICITY"
    STATE_ASSIGNMENT_AMBIGUOUS = "STATE_ASSIGNMENT_AMBIGUOUS"
    MULTI_AXIS_AMBIGUOUS = "MULTI_AXIS_AMBIGUOUS"


class SeriesQualityStatus(str, Enum):
    HIGH_QUALITY_STATE_FAMILY = "HIGH_QUALITY_STATE_FAMILY"
    SUPPORTIVE_STATE_FAMILY = "SUPPORTIVE_STATE_FAMILY"
    LOW_RECURRENCE_STATE_FAMILY = "LOW_RECURRENCE_STATE_FAMILY"
    AMBIGUOUS_STATE_FAMILY = "AMBIGUOUS_STATE_FAMILY"
    INSUFFICIENT_PEAK_QUALITY = "INSUFFICIENT_PEAK_QUALITY"
    BLOCKED = "BLOCKED"


class ReconciliationStatus(str, Enum):
    FULL_LENGTH_PATTERN_COMPATIBLE = "FULL_LENGTH_PATTERN_COMPATIBLE"
    PARTIALLY_COMPATIBLE_WITH_FULL_LENGTH_PATTERN = "PARTIALLY_COMPATIBLE_WITH_FULL_LENGTH_PATTERN"
    T1_PLUS16_SERIES_ONLY = "T1_PLUS16_SERIES_ONLY"
    T1_PLUS18_SERIES_ONLY = "T1_PLUS18_SERIES_ONLY"
    T1_SERIES_NOT_OBSERVED = "T1_SERIES_NOT_OBSERVED"
    INSUFFICIENT_T1_EVIDENCE = "INSUFFICIENT_T1_EVIDENCE"
    AMBIGUOUS_FRAGMENT_LOCALIZATION = "AMBIGUOUS_FRAGMENT_LOCALIZATION"


@dataclass(frozen=True)
class StateDeltaDefinition:
    state_label: StateLabel
    target_neutral_delta: float
    mass_definition: DeltaMassDefinition
    provenance: str


def build_default_state_delta_definitions() -> tuple[StateDeltaDefinition, ...]:
    refs = build_chemical_delta_reference_registry()
    oxygen = next(x for x in refs if x.reference_category is ChemicalReferenceCategory.OXYGEN_ADDITION_EQUIVALENT and x.reference_mass_definition is DeltaMassDefinition.MONOISOTOPIC_DELTA)
    water = next(x for x in refs if x.reference_category is ChemicalReferenceCategory.WATER_ADDITION_EQUIVALENT and x.reference_mass_definition is DeltaMassDefinition.MONOISOTOPIC_DELTA)
    o, w = oxygen.signed_delta_da, water.signed_delta_da
    return (
        StateDeltaDefinition(StateLabel.BASE_STATE, 0.0, DeltaMassDefinition.MONOISOTOPIC_DELTA, "UNMODIFIED_FRAGMENT_REFERENCE"),
        StateDeltaDefinition(StateLabel.PLUS_16_EQUIVALENT, o, DeltaMassDefinition.MONOISOTOPIC_DELTA, oxygen.reference_id),
        StateDeltaDefinition(StateLabel.PLUS_18_EQUIVALENT, w, DeltaMassDefinition.MONOISOTOPIC_DELTA, water.reference_id),
        StateDeltaDefinition(StateLabel.PLUS_32_EQUIVALENT, 2 * o, DeltaMassDefinition.MONOISOTOPIC_DELTA, f"2*{oxygen.reference_id}"),
        StateDeltaDefinition(StateLabel.PLUS_34_EQUIVALENT, w + o, DeltaMassDefinition.MONOISOTOPIC_DELTA, f"{water.reference_id}+{oxygen.reference_id}"),
    )


@dataclass(frozen=True)
class StateSeriesAuditParameters:
    fragment_matching: T1FragmentMatchParameters = field(default_factory=T1FragmentMatchParameters)
    minimum_recurrence_fraction: float = 0.005
    high_recurrence_fraction: float = 0.02
    minimum_relative_prominence: float = 0.001
    apex_centroid_neutral_disagreement_da: float = 0.05

    def validate(self) -> None:
        self.fragment_matching.validate()
        if not 0 <= self.minimum_recurrence_fraction <= self.high_recurrence_fraction <= 1:
            raise ValueError("invalid recurrence thresholds")
        if self.apex_centroid_neutral_disagreement_da <= 0:
            raise ValueError("positive apex/centroid threshold required")


@dataclass(frozen=True, kw_only=True)
class T1FragmentIonHypothesis(StateSeriesSafeguards):
    ion_hypothesis_id: str
    fragment_id: str
    fragment_sequence: str
    start_position: int
    end_position: int
    theoretical_neutral_mass: float
    fragment_length: int
    base_composition: tuple[tuple[str, int], ...]
    contains_g: bool
    cleavage_start: str
    cleavage_end: str
    terminal_chemistry: str
    generation_status: str
    generation_block_reasons: tuple[str, ...]
    ion_mode: T1IonMode
    charge: int
    adduct_hypothesis: str
    theoretical_mz: float
    ion_hypothesis_status: str
    ion_hypothesis_block_reasons: tuple[str, ...]
    source_candidate: T1IonCandidate = field(repr=False, compare=False)


@dataclass(frozen=True, kw_only=True)
class T1FragmentStateCandidate(StateSeriesSafeguards):
    state_candidate_id: str
    fragment_id: str
    fragment_sequence: str
    start_position: int
    end_position: int
    charge: int
    ion_mode: T1IonMode
    state_label: StateLabel
    target_neutral_delta: float
    expected_mz_delta: float
    expected_mz: float
    observed_peak_id: str
    observed_mz: float
    observed_centroid_mz: float | None
    neutral_delta_from_base: float
    centroid_neutral_delta_from_base: float | None
    delta_error_da: float
    centroid_delta_error_da: float | None
    delta_error_ppm: float | None
    scan_recurrence_fraction: float
    intensity_rank: int
    prominence: float | None
    relative_prominence: float | None
    fwhm: float | None
    state_match_status: StateMatchStatus
    candidate_count_for_peak: int
    candidate_count_for_fragment_ion_state: int
    distinct_fragment_count_for_peak: int
    distinct_charge_count_for_peak: int
    distinct_state_count_for_peak: int
    state_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class T1FragmentPeakMatch(StateSeriesSafeguards):
    match_id: str
    fragment_id: str
    peak_id: str
    fragment_sequence: str
    start_position: int
    end_position: int
    theoretical_neutral_mass: float
    ion_mode: T1IonMode
    charge: int
    theoretical_mz: float
    observed_apex_mz: float
    observed_centroid_mz: float | None
    delta_mz: float
    absolute_delta_mz: float
    ppm_error: float
    apex_match_status: StateMatchStatus
    centroid_match_status: str
    scan_recurrence_fraction: float
    prominence: float | None
    fwhm: float | None
    candidate_count_for_peak: int
    candidate_count_for_fragment_ion: int
    fragment_ambiguity_status: str
    charge_ambiguity_status: str
    match_quality_status: str
    match_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class T1FragmentStateFamily(StateSeriesSafeguards):
    state_family_id: str
    fragment_id: str
    fragment_sequence: str
    start_position: int
    end_position: int
    localization_level: str
    charge: int
    ion_mode: T1IonMode
    base_peak_id: str
    base_observed_mz: float
    detected_state_labels: tuple[StateLabel, ...]
    detected_state_count: int
    state_series_pattern: str
    observed_neutral_deltas: tuple[float, ...]
    state_mass_errors: tuple[float, ...]
    state_peak_ids: tuple[str, ...]
    state_recurrence_fractions: tuple[float, ...]
    state_intensity_ranks: tuple[int, ...]
    missing_expected_states: tuple[StateLabel, ...]
    extra_unresolved_states: tuple[str, ...]
    series_ambiguity_status: SeriesAmbiguityStatus
    series_quality_status: SeriesQualityStatus
    series_confidence: str
    series_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class T1FullLengthSeriesReconciliation(StateSeriesSafeguards):
    reconciliation_id: str
    state_family_id: str
    full_length_normalized_deltas: tuple[float, ...]
    t1_detected_state_labels: tuple[StateLabel, ...]
    reconciliation_status: ReconciliationStatus
    localization_level: str
    reconciliation_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class T1StateRunSummary(StateSeriesSafeguards):
    source_id: str
    rna_identity: str
    digest_type: str
    context_source: str
    context_confidence: str
    input_path: str
    status: str
    ms1_spectra_used: int
    ms2_spectra_present: int
    ms2_spectra_excluded: int
    unknown_ms_level_excluded: int
    aggregation_method: str
    grid_method: str
    per_scan_normalization: str
    smoothing_method: str
    baseline_method: str
    peak_detection_method: str
    detected_peak_count: int
    selected_peak_count: int
    theoretical_fragment_count: int
    negative_ion_hypothesis_count: int
    positive_ion_hypothesis_count: int
    positive_ion_hypothesis_status: str
    block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class T1FragmentStateSeriesSummary(StateSeriesSafeguards):
    fragment_match_count: int
    strict_match_count: int
    supportive_match_count: int
    unambiguous_match_count: int
    fragment_ambiguous_match_count: int
    charge_ambiguous_match_count: int
    median_absolute_match_error_mz: float | None
    state_family_count: int
    plus_16_family_count: int
    plus_18_family_count: int
    plus_32_family_count: int
    plus_34_family_count: int
    high_quality_family_count: int
    supportive_family_count: int
    low_quality_or_ambiguous_family_count: int
    full_length_compatible_family_count: int
    full_length_partial_family_count: int
    overall_block_reasons: tuple[str, ...]


@dataclass(frozen=True)
class T1FragmentStateSeriesAuditResult:
    parameters: StateSeriesAuditParameters
    run_profile: ReplicateRunPeakProfile
    run_summary: T1StateRunSummary
    theoretical_fragments: tuple[TheoreticalT1Fragment, ...]
    ion_hypotheses: tuple[T1FragmentIonHypothesis, ...]
    fragment_matches: tuple[T1FragmentPeakMatch, ...]
    state_candidates: tuple[T1FragmentStateCandidate, ...]
    state_families: tuple[T1FragmentStateFamily, ...]
    reconciliations: tuple[T1FullLengthSeriesReconciliation, ...]
    summary: T1FragmentStateSeriesSummary
    algorithm_version: str = ALGORITHM_VERSION
    formal_propagation: bool = False


def build_t1_run_peak_profile(mzml_path: Path, *, source_metadata_record: MzMLSourceMetadataRecord | None = None, detection_config: Mapping[str, Any] | None = None) -> ReplicateRunPeakProfile:
    return build_ms1_peak_profile_for_run(Path(mzml_path), metadata_record=source_metadata_record, detection_config=detection_config)


def generate_t1_fragment_ion_hypotheses(fragments: Sequence[TheoreticalT1Fragment], *, ion_mode: str = "negative", charge_config: Mapping[str, Any] | None = None, observed_mz_range: tuple[float, float] | None = None) -> tuple[T1FragmentIonHypothesis, ...]:
    if str(ion_mode).lower() not in {"negative", T1IonMode.NEGATIVE_DEPROTONATED.value.lower()}:
        return ()
    config = dict(charge_config or {})
    params = config.pop("parameters", None) or T1FragmentMatchParameters(**config)
    by_id = {x.theoretical_t1_fragment_id: x for x in fragments}
    ions = generate_t1_ion_candidates(fragments, parameters=params)
    output = []
    for ion in ions:
        if ion.ion_mode is not T1IonMode.NEGATIVE_DEPROTONATED:
            continue
        if observed_mz_range is not None and not observed_mz_range[0] <= ion.theoretical_mz <= observed_mz_range[1]:
            continue
        fragment = by_id[ion.theoretical_t1_fragment_id]
        output.append(T1FragmentIonHypothesis(
            ion_hypothesis_id=ion.ion_candidate_id, fragment_id=ion.theoretical_t1_fragment_id,
            fragment_sequence=ion.fragment_sequence, start_position=ion.start_position,
            end_position=ion.end_position, theoretical_neutral_mass=fragment.neutral_monoisotopic_mass,
            fragment_length=len(fragment.fragment_sequence),
            base_composition=tuple((base, fragment.fragment_sequence.count(base)) for base in "ACGU"),
            contains_g="G" in fragment.fragment_sequence, cleavage_start=fragment.five_prime_state,
            cleavage_end=fragment.three_prime_state, terminal_chemistry=fragment.cleavage_context,
            generation_status="GENERATED_BY_EXISTING_RNASE_T1_DIGEST", generation_block_reasons=(),
            ion_mode=ion.ion_mode, charge=ion.charge, adduct_hypothesis="[M-zH]z-",
            theoretical_mz=ion.theoretical_mz, ion_hypothesis_status="ELIGIBLE_NEGATIVE_SOURCE",
            ion_hypothesis_block_reasons=("SOURCE_POLARITY_NEGATIVE",), source_candidate=ion,
        ))
    return tuple(sorted(output, key=lambda x: (x.theoretical_mz, x.fragment_id, x.charge, x.ion_hypothesis_id)))


def _state_candidates(fragment_ions: Sequence[T1FragmentIonHypothesis], peaks: Sequence[ReplicateRunPeak], definitions: Sequence[StateDeltaDefinition], params: StateSeriesAuditParameters) -> tuple[T1FragmentStateCandidate, ...]:
    ordered_peaks = tuple(sorted(peaks, key=lambda x: (x.apex_mz, x.peak_id)))
    masses = [x.apex_mz for x in ordered_peaks]
    raw: list[dict[str, Any]] = []
    for ion in sorted(fragment_ions, key=lambda x: (x.theoretical_mz, x.fragment_id, x.charge)):
        for definition in definitions:
            expected = ion.theoretical_mz + definition.target_neutral_delta / ion.charge
            lo = bisect_left(masses, expected - params.fragment_matching.exploratory_tolerance_mz)
            hi = bisect_right(masses, expected + params.fragment_matching.exploratory_tolerance_mz)
            for peak in ordered_peaks[lo:hi]:
                error_mz = peak.apex_mz - expected
                error_da = error_mz * ion.charge
                centroid_error = (peak.centroid_mz - expected) * ion.charge if peak.centroid_mz is not None else None
                raw.append(dict(ion=ion, definition=definition, peak=peak, error_da=error_da, centroid_error=centroid_error,
                                status=StateMatchStatus.STRICT if abs(error_mz) <= params.fragment_matching.strict_tolerance_mz else StateMatchStatus.SUPPORTIVE))
    by_peak: dict[str, list[dict[str, Any]]] = {}
    by_ion_state: dict[tuple[str, StateLabel], list[dict[str, Any]]] = {}
    for row in raw:
        by_peak.setdefault(row["peak"].peak_id, []).append(row)
        by_ion_state.setdefault((row["ion"].ion_hypothesis_id, row["definition"].state_label), []).append(row)
    output = []
    for row in raw:
        ion, definition, peak = row["ion"], row["definition"], row["peak"]
        same_peak, same_state = by_peak[peak.peak_id], by_ion_state[(ion.ion_hypothesis_id, definition.state_label)]
        blocks = list(peak.detection_block_reasons)
        if peak.scan_recurrence_fraction < params.minimum_recurrence_fraction: blocks.append("LOW_SCAN_RECURRENCE")
        if peak.relative_prominence is None or peak.relative_prominence < params.minimum_relative_prominence: blocks.append("LOW_PROMINENCE")
        if peak.fwhm is None or peak.fwhm <= 0: blocks.append("INVALID_FWHM")
        if row["centroid_error"] is not None and abs(row["centroid_error"] - row["error_da"]) > params.apex_centroid_neutral_disagreement_da: blocks.append("APEX_CENTROID_DISAGREEMENT")
        fragments = {x["ion"].fragment_id for x in same_peak}; charges = {x["ion"].charge for x in same_peak}; states = {x["definition"].state_label for x in same_peak}
        if len(fragments) > 1: blocks.append("FRAGMENT_AMBIGUITY")
        if len(charges) > 1: blocks.append("CHARGE_AMBIGUITY")
        if len(same_state) > 1: blocks.append("PEAK_MULTIPLICITY")
        if len(states) > 1: blocks.append("STATE_ASSIGNMENT_AMBIGUITY")
        observed_delta = (peak.apex_mz - ion.theoretical_mz) * ion.charge
        centroid_delta = (peak.centroid_mz - ion.theoretical_mz) * ion.charge if peak.centroid_mz is not None else None
        output.append(T1FragmentStateCandidate(
            state_candidate_id=_id("T1STATE", f"{ion.ion_hypothesis_id}|{definition.state_label.value}|{peak.peak_id}"),
            fragment_id=ion.fragment_id, fragment_sequence=ion.fragment_sequence,
            start_position=ion.start_position, end_position=ion.end_position, charge=ion.charge,
            ion_mode=ion.ion_mode, state_label=definition.state_label,
            target_neutral_delta=definition.target_neutral_delta,
            expected_mz_delta=definition.target_neutral_delta / ion.charge, expected_mz=ion.theoretical_mz + definition.target_neutral_delta / ion.charge,
            observed_peak_id=peak.peak_id, observed_mz=peak.apex_mz, observed_centroid_mz=peak.centroid_mz,
            neutral_delta_from_base=observed_delta, centroid_neutral_delta_from_base=centroid_delta,
            delta_error_da=row["error_da"], centroid_delta_error_da=row["centroid_error"],
            delta_error_ppm=(row["error_da"] / definition.target_neutral_delta * 1e6 if definition.target_neutral_delta else None),
            scan_recurrence_fraction=peak.scan_recurrence_fraction, intensity_rank=peak.intensity_rank,
            prominence=peak.prominence, relative_prominence=peak.relative_prominence, fwhm=peak.fwhm,
            state_match_status=row["status"], candidate_count_for_peak=len(same_peak),
            candidate_count_for_fragment_ion_state=len(same_state), distinct_fragment_count_for_peak=len(fragments),
            distinct_charge_count_for_peak=len(charges), distinct_state_count_for_peak=len(states),
            state_block_reasons=_blocks(blocks),
        ))
    return tuple(sorted(output, key=lambda x: (x.fragment_id, x.charge, x.state_label.value, abs(x.delta_error_da), x.observed_mz, x.observed_peak_id)))


def match_t1_fragment_ions_to_peaks(fragment_ions: Sequence[T1FragmentIonHypothesis], peaks: Sequence[ReplicateRunPeak], *, matching_config: Mapping[str, Any] | None = None) -> tuple[T1FragmentPeakMatch, ...]:
    config = dict(matching_config or {})
    params = config.pop("parameters", None) or StateSeriesAuditParameters(**config)
    params.validate()
    candidates = _state_candidates(fragment_ions, peaks, (build_default_state_delta_definitions()[0],), params)
    output = []
    counts_by_ion: dict[str, int] = {}
    for candidate in candidates:
        counts_by_ion[candidate.fragment_id + f"|{candidate.charge}"] = counts_by_ion.get(candidate.fragment_id + f"|{candidate.charge}", 0) + 1
    neutral_by_fragment = {x.fragment_id: x.theoretical_neutral_mass for x in fragment_ions}
    for x in candidates:
        blocks = list(x.state_block_reasons)
        fragment_status = "FRAGMENT_AMBIGUOUS" if x.distinct_fragment_count_for_peak > 1 else "UNAMBIGUOUS"
        charge_status = "CHARGE_AMBIGUOUS" if x.distinct_charge_count_for_peak > 1 else "UNAMBIGUOUS"
        quality = "AMBIGUOUS" if fragment_status != "UNAMBIGUOUS" or charge_status != "UNAMBIGUOUS" else x.state_match_status.value
        centroid_status = "NOT_RECORDED" if x.observed_centroid_mz is None else ("AGREES" if "APEX_CENTROID_DISAGREEMENT" not in blocks else "DISAGREES")
        output.append(T1FragmentPeakMatch(
            match_id=_id("T1MATCH", x.state_candidate_id), fragment_id=x.fragment_id,
            peak_id=x.observed_peak_id, fragment_sequence=x.fragment_sequence,
            start_position=x.start_position, end_position=x.end_position,
            theoretical_neutral_mass=neutral_by_fragment[x.fragment_id], ion_mode=x.ion_mode,
            charge=x.charge, theoretical_mz=x.expected_mz, observed_apex_mz=x.observed_mz,
            observed_centroid_mz=x.observed_centroid_mz, delta_mz=x.observed_mz - x.expected_mz,
            absolute_delta_mz=abs(x.observed_mz - x.expected_mz),
            ppm_error=(x.observed_mz - x.expected_mz) / x.expected_mz * 1e6,
            apex_match_status=x.state_match_status, centroid_match_status=centroid_status,
            scan_recurrence_fraction=x.scan_recurrence_fraction, prominence=x.prominence, fwhm=x.fwhm,
            candidate_count_for_peak=x.candidate_count_for_peak,
            candidate_count_for_fragment_ion=counts_by_ion[x.fragment_id + f"|{x.charge}"],
            fragment_ambiguity_status=fragment_status, charge_ambiguity_status=charge_status,
            match_quality_status=quality, match_block_reasons=_blocks(blocks),
        ))
    return tuple(sorted(output, key=lambda x: (x.observed_apex_mz, abs(x.delta_mz), x.fragment_id, x.charge, x.match_id)))


def _ambiguity(rows: Sequence[T1FragmentStateCandidate]) -> SeriesAmbiguityStatus:
    axes = []
    if any(x.distinct_fragment_count_for_peak > 1 for x in rows): axes.append(SeriesAmbiguityStatus.FRAGMENT_AMBIGUOUS)
    if any(x.distinct_charge_count_for_peak > 1 for x in rows): axes.append(SeriesAmbiguityStatus.CHARGE_AMBIGUOUS)
    if any(x.candidate_count_for_fragment_ion_state > 1 for x in rows): axes.append(SeriesAmbiguityStatus.PEAK_MULTIPLICITY)
    if any(x.distinct_state_count_for_peak > 1 for x in rows): axes.append(SeriesAmbiguityStatus.STATE_ASSIGNMENT_AMBIGUOUS)
    return SeriesAmbiguityStatus.UNAMBIGUOUS if not axes else axes[0] if len(axes) == 1 else SeriesAmbiguityStatus.MULTI_AXIS_AMBIGUOUS


def _pattern(labels: set[StateLabel]) -> str:
    ordered = [StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT, StateLabel.PLUS_18_EQUIVALENT, StateLabel.PLUS_32_EQUIVALENT, StateLabel.PLUS_34_EQUIVALENT]
    if StateLabel.BASE_STATE not in labels or len(labels) < 2: return "UNRESOLVED_SERIES"
    known = {
        frozenset((StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT)): "BASE__PLUS16",
        frozenset((StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT)): "BASE__PLUS18",
        frozenset((StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT, StateLabel.PLUS_32_EQUIVALENT)): "BASE__PLUS16__PLUS32",
        frozenset((StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT, StateLabel.PLUS_34_EQUIVALENT)): "BASE__PLUS18__PLUS34",
    }
    return known.get(frozenset(labels), "PARTIAL_SERIES__" + "__".join(x.value.replace("_EQUIVALENT", "").replace("_STATE", "") for x in ordered if x in labels))


def build_t1_fragment_state_families(fragment_matches: Sequence[T1FragmentPeakMatch] | None = None, *, state_delta_definitions: Sequence[StateDeltaDefinition] | None = None, state_candidates: Sequence[T1FragmentStateCandidate] | None = None, parameters: StateSeriesAuditParameters | None = None) -> tuple[T1FragmentStateFamily, ...]:
    del fragment_matches  # state candidates preserve the necessary state-axis information
    params = parameters or StateSeriesAuditParameters(); params.validate()
    definitions = tuple(state_delta_definitions or build_default_state_delta_definitions())
    groups: dict[tuple[str, int, T1IonMode], list[T1FragmentStateCandidate]] = {}
    for row in tuple(state_candidates or ()):
        groups.setdefault((row.fragment_id, row.charge, row.ion_mode), []).append(row)
    expected = {x.state_label for x in definitions}
    families = []
    for key, rows in sorted(groups.items(), key=lambda x: (x[0][0], x[0][1], x[0][2].value)):
        labels = {x.state_label for x in rows}
        if StateLabel.BASE_STATE not in labels or len(labels) < 2: continue
        chosen = []
        for label in sorted(labels, key=lambda x: x.value):
            candidates = [x for x in rows if x.state_label is label]
            chosen.append(min(candidates, key=lambda x: (abs(x.delta_error_da), -x.scan_recurrence_fraction, x.intensity_rank, x.observed_peak_id)))
        ambiguity = _ambiguity(rows)
        blocks = [b for x in chosen for b in x.state_block_reasons]
        missing = expected - labels
        if StateLabel.PLUS_32_EQUIVALENT in labels and StateLabel.PLUS_16_EQUIVALENT not in labels: blocks.append("INCOMPLETE_STATE_SERIES")
        if StateLabel.PLUS_34_EQUIVALENT in labels and StateLabel.PLUS_18_EQUIVALENT not in labels: blocks.append("INCOMPLETE_STATE_SERIES")
        if ambiguity is not SeriesAmbiguityStatus.UNAMBIGUOUS: blocks.append("STATE_ASSIGNMENT_AMBIGUITY")
        low_recurrence = any(x.scan_recurrence_fraction < params.minimum_recurrence_fraction for x in chosen)
        poor_shape = any(x.fwhm is None or x.fwhm <= 0 or x.relative_prominence is None for x in chosen)
        if ambiguity is not SeriesAmbiguityStatus.UNAMBIGUOUS: quality, confidence = SeriesQualityStatus.AMBIGUOUS_STATE_FAMILY, "LOW"
        elif low_recurrence: quality, confidence = SeriesQualityStatus.LOW_RECURRENCE_STATE_FAMILY, "LOW"
        elif poor_shape: quality, confidence = SeriesQualityStatus.INSUFFICIENT_PEAK_QUALITY, "LOW"
        elif "APEX_CENTROID_DISAGREEMENT" in blocks: quality, confidence = SeriesQualityStatus.SUPPORTIVE_STATE_FAMILY, "MEDIUM"
        elif all(x.state_match_status is StateMatchStatus.STRICT and x.scan_recurrence_fraction >= params.high_recurrence_fraction for x in chosen): quality, confidence = SeriesQualityStatus.HIGH_QUALITY_STATE_FAMILY, "HIGH"
        else: quality, confidence = SeriesQualityStatus.SUPPORTIVE_STATE_FAMILY, "MEDIUM"
        base = next(x for x in chosen if x.state_label is StateLabel.BASE_STATE)
        ordered = sorted(chosen, key=lambda x: (x.target_neutral_delta, x.observed_peak_id))
        families.append(T1FragmentStateFamily(
            state_family_id=_id("T1FAMILY", f"{key[0]}|{key[1]}|{key[2].value}|{'|'.join(x.observed_peak_id for x in ordered)}"),
            fragment_id=base.fragment_id, fragment_sequence=base.fragment_sequence,
            start_position=base.start_position, end_position=base.end_position,
            localization_level="FRAGMENT_RANGE_ONLY", charge=base.charge, ion_mode=base.ion_mode,
            base_peak_id=base.observed_peak_id, base_observed_mz=base.observed_mz,
            detected_state_labels=tuple(x.state_label for x in ordered), detected_state_count=len(ordered),
            state_series_pattern=_pattern(labels), observed_neutral_deltas=tuple(x.neutral_delta_from_base for x in ordered),
            state_mass_errors=tuple(x.delta_error_da for x in ordered), state_peak_ids=tuple(x.observed_peak_id for x in ordered),
            state_recurrence_fractions=tuple(x.scan_recurrence_fraction for x in ordered),
            state_intensity_ranks=tuple(x.intensity_rank for x in ordered),
            missing_expected_states=tuple(sorted(missing, key=lambda x: x.value)), extra_unresolved_states=(),
            series_ambiguity_status=ambiguity, series_quality_status=quality, series_confidence=confidence,
            series_block_reasons=_blocks(blocks),
        ))
    return tuple(sorted(families, key=lambda x: (-x.detected_state_count, x.fragment_id, x.charge, x.base_observed_mz, x.state_family_id)))


def reconcile_t1_state_families_with_full_length_series(state_families: Sequence[T1FragmentStateFamily], full_length_series: Sequence[float] | None) -> tuple[T1FullLengthSeriesReconciliation, ...]:
    normalized = () if not full_length_series else tuple(float(x) - float(full_length_series[0]) for x in full_length_series)
    if not state_families:
        status = ReconciliationStatus.T1_SERIES_NOT_OBSERVED if normalized else ReconciliationStatus.INSUFFICIENT_T1_EVIDENCE
        blocks = ("NO_STATE_SERIES",) if normalized else ("NO_STATE_SERIES", "FULL_LENGTH_SERIES_RESULT_MISSING")
        return (T1FullLengthSeriesReconciliation(
            reconciliation_id=_id("T1RECON", "NO_STATE_FAMILY"), state_family_id="NO_STATE_FAMILY",
            full_length_normalized_deltas=normalized, t1_detected_state_labels=(),
            reconciliation_status=status, localization_level="NOT_APPLICABLE",
            reconciliation_block_reasons=_blocks(blocks),
        ),)
    output = []
    for family in sorted(state_families, key=lambda x: x.state_family_id):
        labels = set(family.detected_state_labels); blocks = []
        if not normalized:
            status = ReconciliationStatus.INSUFFICIENT_T1_EVIDENCE; blocks.append("FULL_LENGTH_SERIES_RESULT_MISSING")
        elif family.series_ambiguity_status is not SeriesAmbiguityStatus.UNAMBIGUOUS:
            status = ReconciliationStatus.AMBIGUOUS_FRAGMENT_LOCALIZATION; blocks.append("AMBIGUOUS_FRAGMENT_LOCALIZATION")
        elif {StateLabel.PLUS_18_EQUIVALENT, StateLabel.PLUS_34_EQUIVALENT} <= labels:
            status = ReconciliationStatus.FULL_LENGTH_PATTERN_COMPATIBLE
        elif StateLabel.PLUS_18_EQUIVALENT in labels and labels & {StateLabel.PLUS_16_EQUIVALENT, StateLabel.PLUS_34_EQUIVALENT}:
            status = ReconciliationStatus.PARTIALLY_COMPATIBLE_WITH_FULL_LENGTH_PATTERN
        elif StateLabel.PLUS_16_EQUIVALENT in labels:
            status = ReconciliationStatus.T1_PLUS16_SERIES_ONLY
        elif StateLabel.PLUS_18_EQUIVALENT in labels:
            status = ReconciliationStatus.T1_PLUS18_SERIES_ONLY
        else:
            status = ReconciliationStatus.T1_SERIES_NOT_OBSERVED
        output.append(T1FullLengthSeriesReconciliation(
            reconciliation_id=_id("T1RECON", family.state_family_id), state_family_id=family.state_family_id,
            full_length_normalized_deltas=normalized, t1_detected_state_labels=family.detected_state_labels,
            reconciliation_status=status, localization_level="FRAGMENT_RANGE_ONLY",
            reconciliation_block_reasons=_blocks(blocks),
        ))
    return tuple(output)


def _source_blocks(metadata: MzMLSourceMetadataRecord | None, profile: ReplicateRunPeakProfile, sequence: str) -> list[str]:
    blocks = list(profile.block_reasons)
    if metadata is None: blocks.append("SOURCE_METADATA_RECORD_MISSING")
    else:
        if metadata.context_source != "USER_PROVIDED_RUNTIME_MANIFEST": blocks.append("USER_MANIFEST_CONTEXT_MISSING")
        if metadata.polarity_status is PolarityStatus.MIXED_POLARITY: blocks.append("MIXED_POLARITY_INPUT")
        elif metadata.polarity_status is PolarityStatus.NOT_RECORDED: blocks.append("MISSING_POLARITY_METADATA")
        elif metadata.polarity_status is not PolarityStatus.NEGATIVE_ONLY: blocks.append("SOURCE_POLARITY_NOT_NEGATIVE")
        if metadata.representation_status is RepresentationStatus.NOT_RECORDED: blocks.append("MISSING_REPRESENTATION_METADATA")
        elif metadata.representation_status is not RepresentationStatus.PROFILE_ONLY: blocks.append("REPRESENTATION_NOT_PROFILE")
    if not sequence: blocks.append("SEQUENCE_MISSING")
    if profile.ms1_spectra_used == 0: blocks.append("NO_MS1_SPECTRA")
    if not profile.peaks: blocks.append("NO_DETECTED_PEAKS")
    return blocks


def audit_t1_fragment_state_series(mzml_path: Path, sequence: str, *, manifest: Any, rna_identity_id: str, base_masses: Mapping[str, Any], source_metadata_record: MzMLSourceMetadataRecord | None = None, runtime_context: Mapping[str, Any] | None = None, full_length_series: Sequence[float] | None = None, candidate_states: Sequence[str] | None = None, detection_config: Mapping[str, Any] | None = None, parameters: StateSeriesAuditParameters | None = None, state_delta_definitions: Sequence[StateDeltaDefinition] | None = None) -> T1FragmentStateSeriesAuditResult:
    params = parameters or StateSeriesAuditParameters(); params.validate()
    runtime = dict(runtime_context or {})
    profile = build_t1_run_peak_profile(Path(mzml_path), source_metadata_record=source_metadata_record, detection_config=detection_config)
    blocks = _source_blocks(source_metadata_record, profile, sequence)
    fragments: tuple[TheoreticalT1Fragment, ...] = ()
    if sequence and not any(x in blocks for x in ("SOURCE_POLARITY_NOT_NEGATIVE", "MIXED_POLARITY_INPUT", "MISSING_POLARITY_METADATA", "REPRESENTATION_NOT_PROFILE", "MISSING_REPRESENTATION_METADATA")):
        identity = __import__("rna_masshunter.sciex_sample_manifest", fromlist=["get_rna_identity"]).get_rna_identity(manifest, rna_identity_id)
        if identity.sequence != sequence: raise ValueError("runtime sequence does not match manifest RNA identity")
        try: fragments = generate_theoretical_t1_fragments(manifest, rna_identity_id, base_masses, candidate_states=candidate_states)
        except Exception:
            blocks.append("T1_FRAGMENT_GENERATION_FAILED")
    if not fragments: blocks.append("NO_THEORETICAL_FRAGMENTS")
    selected_ids = {p.peak_id for p in profile.peaks if p.detection_status in {T1PeakQualityClass.MAJOR_SHARP.value, T1PeakQualityClass.MAJOR_BROAD.value, T1PeakQualityClass.MINOR_SHARP.value} and "LOW_SCAN_RECURRENCE" not in p.detection_block_reasons}
    selected = tuple(p for p in profile.peaks if p.peak_id in selected_ids)
    observed_range = (min((p.apex_mz for p in profile.peaks), default=0.0), max((p.apex_mz for p in profile.peaks), default=0.0))
    ions = generate_t1_fragment_ion_hypotheses(fragments, charge_config={"parameters": params.fragment_matching}, observed_mz_range=observed_range if profile.peaks else None)
    definitions = tuple(state_delta_definitions or build_default_state_delta_definitions())
    states = _state_candidates(ions, selected, definitions, params)
    matches = match_t1_fragment_ions_to_peaks(ions, selected, matching_config={"parameters": params})
    families = build_t1_fragment_state_families(state_candidates=states, state_delta_definitions=definitions, parameters=params)
    reconciliations = reconcile_t1_state_families_with_full_length_series(families, full_length_series)
    if not matches: blocks.append("NO_FRAGMENT_MATCHES")
    if not families: blocks.append("NO_STATE_SERIES")
    recon_by_id = {x.state_family_id: x for x in reconciliations}
    summary = T1FragmentStateSeriesSummary(
        fragment_match_count=len(matches), strict_match_count=sum(x.apex_match_status is StateMatchStatus.STRICT for x in matches),
        supportive_match_count=sum(x.apex_match_status is StateMatchStatus.SUPPORTIVE for x in matches),
        unambiguous_match_count=sum(x.fragment_ambiguity_status == x.charge_ambiguity_status == "UNAMBIGUOUS" for x in matches),
        fragment_ambiguous_match_count=sum(x.fragment_ambiguity_status != "UNAMBIGUOUS" for x in matches),
        charge_ambiguous_match_count=sum(x.charge_ambiguity_status != "UNAMBIGUOUS" for x in matches),
        median_absolute_match_error_mz=(median(x.absolute_delta_mz for x in matches) if matches else None),
        state_family_count=len(families), plus_16_family_count=sum(StateLabel.PLUS_16_EQUIVALENT in x.detected_state_labels for x in families),
        plus_18_family_count=sum(StateLabel.PLUS_18_EQUIVALENT in x.detected_state_labels for x in families),
        plus_32_family_count=sum(StateLabel.PLUS_32_EQUIVALENT in x.detected_state_labels for x in families),
        plus_34_family_count=sum(StateLabel.PLUS_34_EQUIVALENT in x.detected_state_labels for x in families),
        high_quality_family_count=sum(x.series_quality_status is SeriesQualityStatus.HIGH_QUALITY_STATE_FAMILY for x in families),
        supportive_family_count=sum(x.series_quality_status is SeriesQualityStatus.SUPPORTIVE_STATE_FAMILY for x in families),
        low_quality_or_ambiguous_family_count=sum(x.series_quality_status not in {SeriesQualityStatus.HIGH_QUALITY_STATE_FAMILY, SeriesQualityStatus.SUPPORTIVE_STATE_FAMILY} for x in families),
        full_length_compatible_family_count=sum(recon_by_id[x.state_family_id].reconciliation_status is ReconciliationStatus.FULL_LENGTH_PATTERN_COMPATIBLE for x in families),
        full_length_partial_family_count=sum(recon_by_id[x.state_family_id].reconciliation_status is ReconciliationStatus.PARTIALLY_COMPATIBLE_WITH_FULL_LENGTH_PATTERN for x in families),
        overall_block_reasons=_blocks(blocks + ["CHEMICAL_IDENTITY_UNSUPPORTED"]),
    )
    metadata = source_metadata_record
    run_summary = T1StateRunSummary(
        source_id=profile.run_label, rna_identity=(metadata.rna_identity if metadata else rna_identity_id),
        digest_type=(metadata.digest_type if metadata else runtime.get("Digest_Type", "UNKNOWN")),
        context_source=(metadata.context_source if metadata else runtime.get("Context_Source", "UNKNOWN")),
        context_confidence=(metadata.context_confidence if metadata else runtime.get("Context_Confidence", "UNKNOWN")),
        input_path=str(mzml_path), status="BLOCKED" if any(x in blocks for x in ("NO_MS1_SPECTRA", "SOURCE_POLARITY_NOT_NEGATIVE", "MIXED_POLARITY_INPUT", "MISSING_POLARITY_METADATA", "REPRESENTATION_NOT_PROFILE", "MISSING_REPRESENTATION_METADATA")) else "COMPLETED",
        ms1_spectra_used=profile.ms1_spectra_used, ms2_spectra_present=profile.ms2_spectra_excluded,
        ms2_spectra_excluded=profile.ms2_spectra_excluded, unknown_ms_level_excluded=profile.missing_ms_level_spectra,
        aggregation_method=profile.aggregation_method, grid_method=profile.mz_grid_method,
        per_scan_normalization=profile.intensity_normalization_method, smoothing_method=profile.smoothing_method,
        baseline_method=profile.baseline_method, peak_detection_method=profile.peak_detection_method,
        detected_peak_count=profile.detected_peak_count, selected_peak_count=len(selected),
        theoretical_fragment_count=len(fragments), negative_ion_hypothesis_count=len(ions), positive_ion_hypothesis_count=0,
        positive_ion_hypothesis_status="POSITIVE_ION_HYPOTHESIS_BLOCKED_BY_SOURCE_METADATA",
        block_reasons=_blocks(blocks),
    )
    return T1FragmentStateSeriesAuditResult(params, profile, run_summary, fragments, ions, matches, states, families, reconciliations, summary)


def _record(value: Any) -> dict[str, Any]:
    row = asdict(value)
    row.pop("source_candidate", None)
    def normalize(item: Any) -> Any:
        if isinstance(item, Enum): return item.value
        if isinstance(item, dict): return {k: normalize(v) for k, v in item.items()}
        if isinstance(item, (tuple, list)): return [normalize(v) for v in item]
        return item
    return normalize(row)


def audit_optional_result(result: T1FragmentStateSeriesAuditResult) -> dict[str, Any]:
    peak_records = [_record(x) for x in result.run_profile.peaks]
    safeguards = _record(StateSeriesSafeguards())
    peak_records = [{**x, **safeguards} for x in peak_records]
    theoretical = [_record(x) for x in result.ion_hypotheses]
    return {
        "run_summary_records": [_record(result.run_summary)],
        "peak_records": peak_records,
        "theoretical_fragment_records": theoretical,
        "fragment_match_records": [_record(x) for x in result.fragment_matches],
        "state_family_records": [_record(x) for x in result.state_families],
        "reconciliation_records": [_record(x) for x in result.reconciliations],
        "summary_records": [_record(result.summary)],
    }
