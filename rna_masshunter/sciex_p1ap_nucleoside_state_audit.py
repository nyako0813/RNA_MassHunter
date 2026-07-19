"""Positive-mode P1+AP nucleoside neutral-delta state-series shadow audit.

This optional audit reuses the established P1/SAP chemical candidate generator and
SCIEX MS1 profile builder.  Matches are mass-compatible hypotheses only: chemical
identity, reaction order, and nucleotide/atom localization are never assigned.
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from enum import Enum
from hashlib import sha256
from math import isfinite
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from pyteomics import mzml

from rna_masshunter.masses import PROTON_MASS
from rna_masshunter.p1_sap_chemical_state_audit import (
    MODEL_NOT_DEFINED, NUCLEOSIDE_FORMULAS, _precursor,
    generate_chemical_state_candidates,
)
from rna_masshunter.sciex_intact_peak_family import (
    POTASSIUM_MONOISOTOPIC_MASS_DA, SODIUM_MONOISOTOPIC_MASS_DA,
)
from rna_masshunter.sciex_mzml_source_metadata_audit import (
    MzMLSourceMetadataRecord, PolarityStatus, RepresentationStatus,
)
from rna_masshunter.sciex_t1_fragment_state_series_audit import (
    StateDeltaDefinition, StateLabel, build_default_state_delta_definitions,
)
from rna_masshunter.sciex_t1_profile_peak_audit import T1PeakQualityClass
from rna_masshunter.sciex_t1_replicate_consistency_audit import (
    ReplicateAuditParameters, ReplicateRunPeak, ReplicateRunPeakProfile,
    _decode_array, _rt_minutes, build_ms1_peak_profile_from_spectra,
)

OPTIONAL_RESULT_KEY = "sciex_p1ap_nucleoside_state_audit"
ALGORITHM_VERSION = "sciex-p1ap-nucleoside-state-audit-v1"
ISOTOPE_SPACING_DA = 1.00335483507

_BLOCK_ORDER = (
    "INPUT_FILE_NOT_FOUND", "INPUT_FILE_UNREADABLE", "SOURCE_METADATA_RECORD_MISSING",
    "USER_MANIFEST_CONTEXT_MISSING", "NO_MS1_SPECTRA", "PROFILE_EXTRACTION_FAILED",
    "NO_DETECTED_PEAKS", "SOURCE_POLARITY_NOT_POSITIVE", "MIXED_POLARITY_INPUT",
    "MISSING_POLARITY_METADATA", "REPRESENTATION_NOT_PROFILE",
    "MISSING_REPRESENTATION_METADATA", "NUCLEOSIDE_REGISTRY_MISSING",
    "MODIFICATION_REGISTRY_MISSING", "NO_THEORETICAL_NUCLEOSIDES",
    "NO_ION_HYPOTHESES", "NO_CANDIDATE_MATCHES", "LOW_SCAN_RECURRENCE",
    "LOW_PROMINENCE", "INVALID_FWHM", "APEX_CENTROID_DISAGREEMENT",
    "IDENTITY_AMBIGUITY", "ADDUCT_AMBIGUITY", "POSSIBLE_ISOTOPE_EXPLANATION",
    "POSSIBLE_ADDUCT_EXPLANATION", "POSSIBLE_BACKGROUND_EXPLANATION",
    "STATE_ASSIGNMENT_AMBIGUITY", "NO_STATE_SERIES", "INCOMPLETE_STATE_SERIES",
    "INSUFFICIENT_STATE_SPACING_ACCURACY", "INSUFFICIENT_RT_EVIDENCE",
    "FULL_LENGTH_SERIES_RESULT_MISSING", "T1_STATE_RESULT_MISSING",
    "AMBIGUOUS_NUCLEOSIDE_IDENTITY", "STRUCTURALLY_INCOMPATIBLE_COMBINATION",
    "CHEMICAL_IDENTITY_UNSUPPORTED",
)


def _blocks(values: Iterable[str]) -> tuple[str, ...]:
    found = set(values)
    return tuple(x for x in _BLOCK_ORDER if x in found) + tuple(sorted(found - set(_BLOCK_ORDER)))


def _id(prefix: str, text: str) -> str:
    return prefix + "__" + sha256(text.encode()).hexdigest()[:20].upper()


@dataclass(frozen=True, kw_only=True)
class P1APSafeguards:
    shadow_analysis_only: bool = True
    mass_evidence_only: bool = True
    formal_propagation: bool = False
    chemical_identity_assigned: bool = False
    modification_assigned: bool = False
    exact_nucleotide_localization: bool = False
    exact_atom_localization: bool = False
    reaction_order_assigned: bool = False
    ms2_used_for_formal_identity: bool = False
    ms2_used_for_state_assignment: bool = False
    applied_to_formal_score: bool = False
    applied_to_ranking: bool = False
    applied_to_candidate_filtering: bool = False
    applied_to_final_consensus: bool = False


class NucleosideCandidateClass(str, Enum):
    NEUTRAL_NUCLEOSIDE = "NEUTRAL_NUCLEOSIDE"
    MODIFIED_NEUTRAL_NUCLEOSIDE = "MODIFIED_NEUTRAL_NUCLEOSIDE"
    MASS_ONLY_MODIFIED_NUCLEOSIDE = "MASS_ONLY_MODIFIED_NUCLEOSIDE"
    MONOPHOSPHATE_RESIDUAL = "MONOPHOSPHATE_RESIDUAL"
    CYCLIC_PHOSPHATE_RESIDUAL = "CYCLIC_PHOSPHATE_RESIDUAL"
    BASE_ONLY_OR_IN_SOURCE_FRAGMENT = "BASE_ONLY_OR_IN_SOURCE_FRAGMENT"
    ADDUCTED_NUCLEOSIDE = "ADDUCTED_NUCLEOSIDE"
    UNRESOLVED_LOW_MASS_SPECIES = "UNRESOLVED_LOW_MASS_SPECIES"
    STRUCTURALLY_BLOCKED = "STRUCTURALLY_BLOCKED"


class MatchStatus(str, Enum):
    STRICT = "STRICT"
    SUPPORTIVE = "SUPPORTIVE"


class RTCoelutionStatus(str, Enum):
    COELUTING = "COELUTING"
    PARTIALLY_OVERLAPPING = "PARTIALLY_OVERLAPPING"
    DISTINCT_RETENTION = "DISTINCT_RETENTION"
    INSUFFICIENT_RT_EVIDENCE = "INSUFFICIENT_RT_EVIDENCE"


class StateFamilyQuality(str, Enum):
    HIGH_QUALITY_STATE_FAMILY = "HIGH_QUALITY_STATE_FAMILY"
    SUPPORTIVE_STATE_FAMILY = "SUPPORTIVE_STATE_FAMILY"
    AMBIGUOUS_STATE_FAMILY = "AMBIGUOUS_STATE_FAMILY"
    LOW_QUALITY_STATE_FAMILY = "LOW_QUALITY_STATE_FAMILY"
    BLOCKED = "BLOCKED"


@dataclass(frozen=True)
class P1APAuditParameters:
    strict_tolerance_ppm: float = 5.0
    supportive_tolerance_ppm: float = 10.0
    ms2_precursor_tolerance_ppm: float = 20.0
    scan_relative_intensity_threshold: float = 0.0005
    minimum_recurrence_fraction: float = 0.005
    high_recurrence_fraction: float = 0.02
    minimum_relative_prominence: float = 0.001
    apex_centroid_disagreement_da: float = 0.01
    rt_coelution_apex_tolerance_min: float = 0.08

    def validate(self) -> None:
        if self.strict_tolerance_ppm <= 0 or self.supportive_tolerance_ppm < self.strict_tolerance_ppm:
            raise ValueError("invalid P1/AP matching tolerances")
        if not 0 <= self.scan_relative_intensity_threshold <= 1:
            raise ValueError("invalid scan intensity threshold")
        if not 0 <= self.minimum_recurrence_fraction <= self.high_recurrence_fraction <= 1:
            raise ValueError("invalid recurrence thresholds")


@dataclass(frozen=True, kw_only=True)
class NucleosideCandidate(P1APSafeguards):
    candidate_id: str
    candidate_name: str
    parent_base: str
    candidate_class: NucleosideCandidateClass
    modification_components: tuple[str, ...]
    molecular_formula: str
    theoretical_neutral_mass: float | None
    mass_provenance: str
    source_registry: str
    structure_constraint_status: str
    candidate_block_reasons: tuple[str, ...]
    eligible_for_mass_matching: bool
    source_candidate: Mapping[str, Any] = field(repr=False, compare=False)


@dataclass(frozen=True, kw_only=True)
class NucleosideIonHypothesis(P1APSafeguards):
    ion_hypothesis_id: str
    candidate_id: str
    candidate_name: str
    parent_base: str
    candidate_class: NucleosideCandidateClass
    ion_mode: str
    charge: int
    adduct_type: str
    adduct_mass: float
    theoretical_mz: float
    ion_hypothesis_status: str
    ion_hypothesis_block_reasons: tuple[str, ...]


@dataclass(frozen=True)
class TargetRTSummary:
    target_id: str
    first_supporting_rt: float | None
    last_supporting_rt: float | None
    rt_apex: float | None
    rt_centroid: float | None
    supporting_scan_count: int
    total_ms1_scan_count: int
    scan_recurrence_fraction: float
    rt_span: float | None
    rt_profile_status: str


@dataclass(frozen=True)
class P1APRunPeakProfile:
    profile: ReplicateRunPeakProfile
    rt_evidence: Mapping[str, TargetRTSummary]
    compatible_ms2_counts: Mapping[str, int]
    nearest_ms2_precursor: Mapping[str, tuple[float, float]]
    target_count: int


@dataclass(frozen=True, kw_only=True)
class NucleosideStateCandidate(P1APSafeguards):
    state_candidate_id: str
    candidate_id: str
    candidate_name: str
    parent_base: str
    candidate_class: NucleosideCandidateClass
    ion_hypothesis_id: str
    ion_mode: str
    charge: int
    adduct_type: str
    state_label: StateLabel
    target_neutral_delta: float
    expected_mz: float
    observed_peak_id: str
    observed_apex_mz: float
    observed_centroid_mz: float | None
    observed_neutral_delta: float
    delta_error_da: float
    delta_error_ppm: float | None
    match_status: MatchStatus
    intensity_rank: int
    scan_recurrence_fraction: float
    prominence: float | None
    fwhm: float | None
    rt_apex: float | None
    rt_centroid: float | None
    first_supporting_rt: float | None
    last_supporting_rt: float | None
    rt_span: float | None
    candidate_count_for_peak: int
    peak_count_for_candidate_state: int
    distinct_identity_count_for_peak: int
    distinct_adduct_count_for_peak: int
    isotope_compatibility: str
    known_adduct_compatibility: str
    alternative_explanation_count: int
    alternative_explanation_status: str
    state_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class NucleosidePeakMatch(P1APSafeguards):
    match_id: str
    candidate_id: str
    candidate_name: str
    parent_base: str
    candidate_class: NucleosideCandidateClass
    observed_peak_id: str
    theoretical_neutral_mass: float
    ion_mode: str
    charge: int
    adduct_type: str
    theoretical_mz: float
    observed_apex_mz: float
    observed_centroid_mz: float | None
    delta_mz: float
    absolute_delta_mz: float
    ppm_error: float
    apex_match_status: MatchStatus
    centroid_match_status: str
    intensity_rank: int
    scan_recurrence_fraction: float
    rt_apex: float | None
    rt_centroid: float | None
    candidate_count_for_peak: int
    peak_count_for_candidate_ion: int
    identity_ambiguity_status: str
    adduct_ambiguity_status: str
    match_quality_status: str
    compatible_ms2_spectrum_count: int
    nearest_ms2_precursor_mz: float | None
    ms2_precursor_delta: float | None
    ms2_availability_status: str
    match_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class NucleosideStateFamily(P1APSafeguards):
    state_family_id: str
    base_candidate_id: str
    base_candidate_name: str
    parent_base: str
    ion_mode: str
    charge: int
    adduct_type: str
    base_peak_id: str
    base_observed_mz: float
    detected_state_labels: tuple[StateLabel, ...]
    detected_state_count: int
    observed_neutral_deltas: tuple[float, ...]
    expected_neutral_deltas: tuple[float, ...]
    delta_errors: tuple[float, ...]
    base_rt_apex: float | None
    state_rt_apices: tuple[float | None, ...]
    delta_rts: tuple[float | None, ...]
    rt_overlap_fraction: float | None
    rt_coelution_status: RTCoelutionStatus
    distinct_chromatographic_feature_status: str
    missing_expected_states: tuple[StateLabel, ...]
    extra_unresolved_states: tuple[str, ...]
    identity_ambiguity_status: str
    adduct_ambiguity_status: str
    isotope_compatibility: str
    known_adduct_compatibility: str
    charge_scaled_spacing_compatibility: bool
    alternative_explanation_count: int
    alternative_explanation_status: str
    series_quality_status: StateFamilyQuality
    series_confidence: str
    series_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class CrossDigestReconciliation(P1APSafeguards):
    reconciliation_id: str
    state_family_id: str
    p1ap_state_status: str
    t1_reconciliation_status: str
    full_length_reconciliation_status: str
    full_length_normalized_deltas: tuple[float, ...]
    localization_status: str
    reconciliation_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class P1APRunSummary(P1APSafeguards):
    source_id: str
    rna_identity: str
    digest_type: str
    digest_detail: str
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
    rt_evidence_method: str
    block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class P1APNucleosideStateSummary(P1APSafeguards):
    canonical_nucleoside_hypothesis_count: int
    modified_nucleoside_hypothesis_count: int
    structure_blocked_candidate_count: int
    ion_hypothesis_count: int
    candidate_match_count: int
    strict_match_count: int
    supportive_match_count: int
    unique_matched_peak_count: int
    unambiguous_match_count: int
    identity_ambiguous_match_count: int
    adduct_ambiguous_match_count: int
    state_family_count: int
    plus16_family_count: int
    plus18_family_count: int
    plus32_family_count: int
    plus34_family_count: int
    high_quality_state_family_count: int
    supportive_state_family_count: int
    ambiguous_state_family_count: int
    ms2_available_candidate_count: int
    median_absolute_mass_error_mz: float | None
    full_length_reconciliation_status: str
    t1_reconciliation_status: str
    overall_evidence_status: str
    overall_confidence: str
    overall_block_reasons: tuple[str, ...]


@dataclass(frozen=True)
class P1APNucleosideStateAuditResult:
    parameters: P1APAuditParameters
    run_profile: P1APRunPeakProfile
    run_summary: P1APRunSummary
    candidates: tuple[NucleosideCandidate, ...]
    ion_hypotheses: tuple[NucleosideIonHypothesis, ...]
    matches: tuple[NucleosidePeakMatch, ...]
    state_candidates: tuple[NucleosideStateCandidate, ...]
    state_families: tuple[NucleosideStateFamily, ...]
    reconciliations: tuple[CrossDigestReconciliation, ...]
    summary: P1APNucleosideStateSummary
    algorithm_version: str = ALGORITHM_VERSION
    formal_propagation: bool = False


def generate_nucleoside_candidates(*, project_root: Path, sequence: str, modification_registry: Sequence[Any] | None = None, structure_constraints: Sequence[Mapping[str, Any]] | None = None) -> tuple[NucleosideCandidate, ...]:
    if not NUCLEOSIDE_FORMULAS:
        return ()
    rows, _ = generate_chemical_state_candidates(sequence, list(modification_registry or ()), Path(project_root), charges=(1,))
    neutral = [x for x in rows if x.get("Product_Type") == "monomer" and x.get("Chemical_Family") in {"DEPHOSPHORYLATED", "SULFUR_CONTAINING_NON_PT_ALTERNATIVE"}]
    output = []
    for row in neutral:
        state = str(row["Nucleoside_Modification_State"]); canonical = state == "unmodified"
        formula = str(row.get("Elemental_Composition") or MODEL_NOT_DEFINED)
        formula_known = formula != MODEL_NOT_DEFINED
        candidate_class = NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE if canonical else NucleosideCandidateClass.MODIFIED_NEUTRAL_NUCLEOSIDE if formula_known else NucleosideCandidateClass.MASS_ONLY_MODIFIED_NUCLEOSIDE
        provenance = "p1_sap_chemical_state_audit.NUCLEOSIDE_FORMULAS" if canonical else "modification_transforms_v2.yaml:composition_delta" if formula_known else "modifications.yaml:modified_nucleoside_mass_mono"
        output.append(NucleosideCandidate(
            candidate_id=str(row["Chemical_State_ID"]), candidate_name=(str(row["Base_or_Oligomer_Composition"]) if canonical else state),
            parent_base=str(row["Base_or_Oligomer_Composition"]).split(":", 1)[0], candidate_class=candidate_class,
            modification_components=() if canonical else tuple(x for x in state.split("+") if x), molecular_formula=formula,
            theoretical_neutral_mass=float(row["Neutral_Mass"]) if row.get("Neutral_Mass") is not None else None,
            mass_provenance=provenance, source_registry="p1_sap_chemical_state_audit.generate_chemical_state_candidates",
            structure_constraint_status="CANONICAL_FORMULA" if canonical else "TRANSFORM_COMPOSITION_VALIDATED" if formula_known else "MASS_ONLY_STRUCTURE_UNRESOLVED",
            candidate_block_reasons=(), eligible_for_mass_matching=row.get("Neutral_Mass") is not None,
            source_candidate=row,
        ))
    for row in rows:
        if row.get("Product_Type") != "monomer" or row.get("Chemical_Family") != "RESIDUAL_NORMAL_PHOSPHATE":
            continue
        state = str(row["Nucleoside_Modification_State"]); base_key = str(row["Base_or_Oligomer_Composition"]); formula = str(row.get("Elemental_Composition") or MODEL_NOT_DEFINED)
        output.append(NucleosideCandidate(
            candidate_id=str(row["Chemical_State_ID"]), candidate_name=f"{base_key}:{state}:residual_phosphate",
            parent_base=base_key.split(":", 1)[0], candidate_class=NucleosideCandidateClass.MONOPHOSPHATE_RESIDUAL,
            modification_components=tuple(x for x in state.split("+") if x and x != "unmodified"), molecular_formula=formula,
            theoretical_neutral_mass=float(row["Neutral_Mass"]) if row.get("Neutral_Mass") is not None else None,
            mass_provenance="p1_sap_chemical_state_audit.PHOSPHATE_MONOESTER",
            source_registry="p1_sap_chemical_state_audit.generate_chemical_state_candidates",
            structure_constraint_status="DIAGNOSTIC_RESIDUAL_PHOSPHATE_HYPOTHESIS",
            candidate_block_reasons=("CHEMICAL_IDENTITY_UNSUPPORTED",), eligible_for_mass_matching=False,
            source_candidate=row,
        ))
    for serial, item in enumerate(structure_constraints or (), 1):
        if str(item.get("status")) != "STRUCTURALLY_INCOMPATIBLE_COMBINATION":
            continue
        output.append(NucleosideCandidate(
            candidate_id=str(item.get("candidate_id") or f"BLOCKED_{serial}"), candidate_name=str(item.get("candidate_name") or "BLOCKED_COMBINATION"),
            parent_base=str(item.get("parent_base") or "UNKNOWN"), candidate_class=NucleosideCandidateClass.STRUCTURALLY_BLOCKED,
            modification_components=tuple(map(str, item.get("components") or ())), molecular_formula=MODEL_NOT_DEFINED,
            theoretical_neutral_mass=None, mass_provenance="NONE", source_registry="runtime_structure_constraints",
            structure_constraint_status="STRUCTURALLY_INCOMPATIBLE_COMBINATION",
            candidate_block_reasons=("STRUCTURALLY_INCOMPATIBLE_COMBINATION",), eligible_for_mass_matching=False,
            source_candidate=dict(item),
        ))
    return tuple(sorted(output, key=lambda x: (x.candidate_class.value, x.parent_base, x.theoretical_neutral_mass if x.theoretical_neutral_mass is not None else float("inf"), x.candidate_id)))


_ADDUCTS = {"H": (PROTON_MASS, "[M+H]+"), "Na": (SODIUM_MONOISOTOPIC_MASS_DA, "[M+Na]+"), "K": (POTASSIUM_MONOISOTOPIC_MASS_DA, "[M+K]+")}


def generate_positive_nucleoside_ion_hypotheses(candidates: Sequence[NucleosideCandidate], *, ion_config: Mapping[str, Any] | None = None) -> tuple[NucleosideIonHypothesis, ...]:
    config = dict(ion_config or {})
    adducts = tuple(config.pop("adducts", ("H",)))
    if config:
        raise ValueError(f"unsupported ion_config keys: {sorted(config)}")
    unsupported = set(adducts) - set(_ADDUCTS)
    if unsupported:
        raise ValueError(f"unsupported adduct hypotheses: {sorted(unsupported)}")
    output = []
    for candidate in candidates:
        if not candidate.eligible_for_mass_matching or candidate.theoretical_neutral_mass is None:
            continue
        for adduct in adducts:
            mass, label = _ADDUCTS[adduct]
            output.append(NucleosideIonHypothesis(
                ion_hypothesis_id=f"{candidate.candidate_id}__POSITIVE_{adduct}_Z1", candidate_id=candidate.candidate_id,
                candidate_name=candidate.candidate_name, parent_base=candidate.parent_base, candidate_class=candidate.candidate_class,
                ion_mode="POSITIVE", charge=1, adduct_type=label, adduct_mass=mass,
                theoretical_mz=candidate.theoretical_neutral_mass + mass, ion_hypothesis_status="ELIGIBLE_POSITIVE_SOURCE",
                ion_hypothesis_block_reasons=("SOURCE_POLARITY_POSITIVE",),
            ))
    return tuple(sorted(output, key=lambda x: (x.theoretical_mz, x.candidate_id, x.adduct_type)))


def generate_negative_nucleoside_ion_hypotheses(candidates: Sequence[NucleosideCandidate]) -> tuple[NucleosideIonHypothesis, ...]:
    del candidates
    return ()


def _target_id(ion_id: str, label: StateLabel) -> str:
    return f"{ion_id}|{label.value}"


def _targets(ions: Sequence[NucleosideIonHypothesis], definitions: Sequence[StateDeltaDefinition]) -> tuple[tuple[str, float], ...]:
    return tuple(sorted(((_target_id(ion.ion_hypothesis_id, definition.state_label), ion.theoretical_mz + definition.target_neutral_delta / ion.charge) for ion in ions for definition in definitions), key=lambda x: (x[1], x[0])))


def build_p1ap_ms1_peak_profile(mzml_path: Path, *, source_metadata_record: MzMLSourceMetadataRecord | None = None, detection_config: Mapping[str, Any] | None = None, rt_targets: Sequence[tuple[str, float]] = (), ms2_ion_hypotheses: Sequence[NucleosideIonHypothesis] = (), parameters: P1APAuditParameters | None = None) -> P1APRunPeakProfile:
    params = parameters or P1APAuditParameters(); params.validate()
    path = Path(mzml_path)
    if not path.is_file():
        empty = ReplicateRunPeakProfile(path.stem, str(path), "BLOCKED", "NOT_APPLICABLE", 0, 0, 0, 0, 0, "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", 0, (), "UNKNOWN", "UNKNOWN", ("INPUT_FILE_NOT_FOUND",))
        return P1APRunPeakProfile(empty, {}, {}, {}, len(rt_targets))
    target_rows = tuple(sorted(rt_targets, key=lambda x: (x[1], x[0])))
    observations: dict[str, list[tuple[float | None, float, float]]] = defaultdict(list)
    ion_rows = tuple(sorted(ms2_ion_hypotheses, key=lambda x: (x.theoretical_mz, x.ion_hypothesis_id)))
    ion_mzs = [x.theoretical_mz for x in ion_rows]
    ms2_counts: Counter[str] = Counter(); nearest: dict[str, tuple[float, float]] = {}

    def spectra() -> Iterable[Mapping[str, Any]]:
        with mzml.MzML(str(path), decode_binary=False) as reader:
            for spectrum in reader:
                try: level = int(spectrum.get("ms level"))
                except (TypeError, ValueError):
                    yield spectrum; continue
                if level == 1:
                    mz_values = _decode_array(spectrum.get("m/z array")); intensities = _decode_array(spectrum.get("intensity array"))
                    copied = dict(spectrum); copied["m/z array"] = mz_values; copied["intensity array"] = intensities
                    if len(mz_values) and len(mz_values) == len(intensities):
                        base = float(max(np.max(intensities), 0.0)); rt = _rt_minutes(spectrum)
                        if base > 0:
                            for target_id, target_mz in target_rows:
                                tolerance = target_mz * params.supportive_tolerance_ppm * 1e-6
                                lo = int(np.searchsorted(mz_values, target_mz - tolerance, side="left")); hi = int(np.searchsorted(mz_values, target_mz + tolerance, side="right"))
                                if hi <= lo: continue
                                local = intensities[lo:hi]; rel = int(np.argmax(local)); intensity = float(local[rel])
                                if intensity / base >= params.scan_relative_intensity_threshold:
                                    observations[target_id].append((rt, float(mz_values[lo + rel]), intensity))
                    yield copied
                else:
                    if level == 2 and ion_rows:
                        precursor, charge, _, _ = _precursor(spectrum)
                        if precursor is not None:
                            tolerance = precursor * params.ms2_precursor_tolerance_ppm * 1e-6
                            lo = bisect_left(ion_mzs, precursor - tolerance); hi = bisect_right(ion_mzs, precursor + tolerance)
                            for ion in ion_rows[lo:hi]:
                                if charge not in (None, 0, ion.charge): continue
                                error = precursor - ion.theoretical_mz; ms2_counts[ion.ion_hypothesis_id] += 1
                                if ion.ion_hypothesis_id not in nearest or abs(error) < abs(nearest[ion.ion_hypothesis_id][1]): nearest[ion.ion_hypothesis_id] = (precursor, error)
                    yield spectrum
    config = dict(detection_config or {})
    replicate_params = config.pop("parameters", None) or ReplicateAuditParameters(**config.pop("replicate_parameters", {}))
    peak_params = config.pop("peak_detection_parameters", None)
    if config: raise ValueError(f"unsupported detection_config keys: {sorted(config)}")
    run_label = source_metadata_record.technical_run_label if source_metadata_record and source_metadata_record.technical_run_label != "UNKNOWN" else path.stem
    profile = build_ms1_peak_profile_from_spectra(spectra(), run_label=run_label, input_path=str(path), metadata_record=source_metadata_record, parameters=replicate_params, peak_detection_parameters=peak_params)
    rt_evidence = {}
    for target_id, _ in target_rows:
        rows = observations.get(target_id, []); usable = [(rt, mz_value, intensity) for rt, mz_value, intensity in rows if rt is not None]
        if usable:
            apex = max(usable, key=lambda x: (x[2], -x[0])); total = sum(x[2] for x in usable)
            first, last = min(x[0] for x in usable), max(x[0] for x in usable)
            centroid = sum(x[0] * x[2] for x in usable) / total if total else None
            status = "TARGETED_SCAN_LEVEL_RT_EVIDENCE"
        else:
            first = last = centroid = None; apex = (None, None, None); status = "INSUFFICIENT_RT_EVIDENCE"
        rt_evidence[target_id] = TargetRTSummary(target_id, first, last, apex[0], centroid, len(rows), profile.ms1_spectra_used, len(rows) / profile.ms1_spectra_used if profile.ms1_spectra_used else 0.0, last - first if first is not None and last is not None else None, status)
    return P1APRunPeakProfile(profile, rt_evidence, dict(ms2_counts), nearest, len(target_rows))


def _selected_peaks(profile: ReplicateRunPeakProfile) -> tuple[ReplicateRunPeak, ...]:
    allowed = {T1PeakQualityClass.MAJOR_SHARP.value, T1PeakQualityClass.MAJOR_BROAD.value, T1PeakQualityClass.MINOR_SHARP.value}
    return tuple(p for p in profile.peaks if p.detection_status in allowed and "LOW_SCAN_RECURRENCE" not in p.detection_block_reasons)


def _state_search_peaks(profile: ReplicateRunPeakProfile) -> tuple[ReplicateRunPeak, ...]:
    selected = set(_selected_peaks(profile))
    isotope = {p for p in profile.peaks if p.detection_status == T1PeakQualityClass.ISOTOPE_OR_ENVELOPE_COMPONENT.value and "LOW_SCAN_RECURRENCE" not in p.detection_block_reasons}
    return tuple(sorted(selected | isotope, key=lambda x: (x.apex_mz, x.peak_id)))


def _state_candidates(ions: Sequence[NucleosideIonHypothesis], peaks: Sequence[ReplicateRunPeak], definitions: Sequence[StateDeltaDefinition], rt_evidence: Mapping[str, TargetRTSummary], params: P1APAuditParameters) -> tuple[NucleosideStateCandidate, ...]:
    ordered = tuple(sorted(peaks, key=lambda x: (x.apex_mz, x.peak_id))); masses = [x.apex_mz for x in ordered]
    raw = []
    for ion in sorted(ions, key=lambda x: (x.theoretical_mz, x.ion_hypothesis_id)):
        for definition in definitions:
            expected = ion.theoretical_mz + definition.target_neutral_delta / ion.charge
            tolerance = expected * params.supportive_tolerance_ppm * 1e-6
            for peak in ordered[bisect_left(masses, expected - tolerance):bisect_right(masses, expected + tolerance)]:
                error = (peak.apex_mz - expected) * ion.charge
                raw.append((ion, definition, peak, expected, error, MatchStatus.STRICT if abs(peak.apex_mz - expected) <= expected * params.strict_tolerance_ppm * 1e-6 else MatchStatus.SUPPORTIVE))
    by_peak: dict[str, list[Any]] = defaultdict(list); by_state: dict[tuple[str, StateLabel], list[Any]] = defaultdict(list)
    for row in raw: by_peak[row[2].peak_id].append(row); by_state[(row[0].ion_hypothesis_id, row[1].state_label)].append(row)
    output = []
    for ion, definition, peak, expected, error, status in raw:
        same_peak = by_peak[peak.peak_id]; same_state = by_state[(ion.ion_hypothesis_id, definition.state_label)]
        identities = {x[0].candidate_id for x in same_peak}; adducts = {x[0].adduct_type for x in same_peak}; blocks = list(peak.detection_block_reasons)
        centroid_error = (peak.centroid_mz - expected) * ion.charge if peak.centroid_mz is not None else None
        if peak.scan_recurrence_fraction < params.minimum_recurrence_fraction: blocks.append("LOW_SCAN_RECURRENCE")
        if peak.relative_prominence is None or peak.relative_prominence < params.minimum_relative_prominence: blocks.append("LOW_PROMINENCE")
        if peak.fwhm is None or peak.fwhm <= 0: blocks.append("INVALID_FWHM")
        if centroid_error is not None and abs(centroid_error - error) > params.apex_centroid_disagreement_da: blocks.append("APEX_CENTROID_DISAGREEMENT")
        if len(identities) > 1: blocks.append("IDENTITY_AMBIGUITY")
        if len(adducts) > 1: blocks.append("ADDUCT_AMBIGUITY")
        isotope = peak.detection_status == T1PeakQualityClass.ISOTOPE_OR_ENVELOPE_COMPONENT.value
        if isotope: blocks.append("POSSIBLE_ISOTOPE_EXPLANATION")
        if len(same_state) > 1: blocks.append("STATE_ASSIGNMENT_AMBIGUITY")
        trace = rt_evidence.get(_target_id(ion.ion_hypothesis_id, definition.state_label), TargetRTSummary("", None, None, None, None, 0, 0, 0, None, "INSUFFICIENT_RT_EVIDENCE"))
        if trace.rt_apex is None: blocks.append("INSUFFICIENT_RT_EVIDENCE")
        background = trace.rt_span is not None and trace.rt_span > 5.0 and trace.scan_recurrence_fraction > 0.1
        if background: blocks.append("POSSIBLE_BACKGROUND_EXPLANATION")
        alternatives = int(isotope) + int(len(adducts) > 1) + int(len(identities) > 1) + int(background)
        output.append(NucleosideStateCandidate(
            state_candidate_id=_id("P1APSTATE", f"{ion.ion_hypothesis_id}|{definition.state_label.value}|{peak.peak_id}"),
            candidate_id=ion.candidate_id, candidate_name=ion.candidate_name, parent_base=ion.parent_base,
            candidate_class=ion.candidate_class, ion_hypothesis_id=ion.ion_hypothesis_id, ion_mode=ion.ion_mode,
            charge=ion.charge, adduct_type=ion.adduct_type, state_label=definition.state_label,
            target_neutral_delta=definition.target_neutral_delta, expected_mz=expected,
            observed_peak_id=peak.peak_id, observed_apex_mz=peak.apex_mz, observed_centroid_mz=peak.centroid_mz,
            observed_neutral_delta=(peak.apex_mz - ion.theoretical_mz) * ion.charge, delta_error_da=error,
            delta_error_ppm=(error / definition.target_neutral_delta * 1e6 if definition.target_neutral_delta else None),
            match_status=status, intensity_rank=peak.intensity_rank, scan_recurrence_fraction=peak.scan_recurrence_fraction,
            prominence=peak.prominence, fwhm=peak.fwhm, rt_apex=trace.rt_apex, rt_centroid=trace.rt_centroid,
            first_supporting_rt=trace.first_supporting_rt, last_supporting_rt=trace.last_supporting_rt, rt_span=trace.rt_span,
            candidate_count_for_peak=len(same_peak), peak_count_for_candidate_state=len(same_state),
            distinct_identity_count_for_peak=len(identities), distinct_adduct_count_for_peak=len(adducts),
            isotope_compatibility="POSSIBLE" if isotope else "NOT_SUPPORTED_BY_PROFILE_ANNOTATION",
            known_adduct_compatibility="AMBIGUOUS" if len(adducts) > 1 else "NO_ALTERNATIVE_CONFIGURED_ADDUCT_MATCH",
            alternative_explanation_count=alternatives, alternative_explanation_status="ALTERNATIVES_PRESENT" if alternatives else "NO_ENUMERATED_ALTERNATIVE",
            state_block_reasons=_blocks(blocks),
        ))
    return tuple(sorted(output, key=lambda x: (x.candidate_id, x.adduct_type, x.state_label.value, abs(x.delta_error_da), x.observed_apex_mz, x.observed_peak_id)))


def match_nucleoside_ions_to_peaks(ion_hypotheses: Sequence[NucleosideIonHypothesis], peaks: Sequence[ReplicateRunPeak], *, matching_config: Mapping[str, Any] | None = None, rt_evidence: Mapping[str, TargetRTSummary] | None = None, ms2_counts: Mapping[str, int] | None = None, nearest_ms2: Mapping[str, tuple[float, float]] | None = None, neutral_masses: Mapping[str, float] | None = None) -> tuple[NucleosidePeakMatch, ...]:
    config = dict(matching_config or {}); params = config.pop("parameters", None) or P1APAuditParameters(**config); params.validate()
    states = _state_candidates(ion_hypotheses, peaks, (build_default_state_delta_definitions()[0],), rt_evidence or {}, params)
    count_by_ion = Counter(x.ion_hypothesis_id for x in states); masses = dict(neutral_masses or {})
    output = []
    for x in states:
        identity = "IDENTITY_AMBIGUOUS" if x.distinct_identity_count_for_peak > 1 else "UNAMBIGUOUS"
        adduct = "ADDUCT_AMBIGUOUS" if x.distinct_adduct_count_for_peak > 1 else "UNAMBIGUOUS"
        if identity != "UNAMBIGUOUS": quality = "IDENTITY_AMBIGUOUS_MATCH"
        elif adduct != "UNAMBIGUOUS": quality = "ADDUCT_AMBIGUOUS_MATCH"
        elif "LOW_SCAN_RECURRENCE" in x.state_block_reasons: quality = "LOW_RECURRENCE_MATCH"
        elif "POSSIBLE_BACKGROUND_EXPLANATION" in x.state_block_reasons: quality = "POSSIBLE_BACKGROUND_MATCH"
        elif x.match_status is MatchStatus.STRICT and "APEX_CENTROID_DISAGREEMENT" not in x.state_block_reasons: quality = "HIGH_QUALITY_NUCLEOSIDE_MATCH"
        else: quality = "SUPPORTIVE_NUCLEOSIDE_MATCH"
        centroid_status = "NOT_RECORDED" if x.observed_centroid_mz is None else "DISAGREES" if "APEX_CENTROID_DISAGREEMENT" in x.state_block_reasons else "AGREES"
        nearest = (nearest_ms2 or {}).get(x.ion_hypothesis_id); ms2 = int((ms2_counts or {}).get(x.ion_hypothesis_id, 0))
        output.append(NucleosidePeakMatch(
            match_id=_id("P1APMATCH", x.state_candidate_id), candidate_id=x.candidate_id,
            candidate_name=x.candidate_name, parent_base=x.parent_base, candidate_class=x.candidate_class,
            observed_peak_id=x.observed_peak_id, theoretical_neutral_mass=float(masses.get(x.candidate_id, x.expected_mz - PROTON_MASS)),
            ion_mode=x.ion_mode, charge=x.charge, adduct_type=x.adduct_type, theoretical_mz=x.expected_mz,
            observed_apex_mz=x.observed_apex_mz, observed_centroid_mz=x.observed_centroid_mz,
            delta_mz=x.observed_apex_mz - x.expected_mz, absolute_delta_mz=abs(x.observed_apex_mz - x.expected_mz),
            ppm_error=(x.observed_apex_mz - x.expected_mz) / x.expected_mz * 1e6,
            apex_match_status=x.match_status, centroid_match_status=centroid_status,
            intensity_rank=x.intensity_rank, scan_recurrence_fraction=x.scan_recurrence_fraction,
            rt_apex=x.rt_apex, rt_centroid=x.rt_centroid, candidate_count_for_peak=x.candidate_count_for_peak,
            peak_count_for_candidate_ion=count_by_ion[x.ion_hypothesis_id], identity_ambiguity_status=identity,
            adduct_ambiguity_status=adduct, match_quality_status=quality, compatible_ms2_spectrum_count=ms2,
            nearest_ms2_precursor_mz=nearest[0] if nearest else None, ms2_precursor_delta=nearest[1] if nearest else None,
            ms2_availability_status="PRECURSOR_COMPATIBLE_MS2_PRESENT_NOT_INTERPRETED" if ms2 else "NO_PRECURSOR_COMPATIBLE_MS2",
            match_block_reasons=x.state_block_reasons,
        ))
    return tuple(sorted(output, key=lambda x: (x.observed_apex_mz, x.absolute_delta_mz, x.candidate_id, x.adduct_type, x.match_id)))


def _rt_relationship(base: NucleosideStateCandidate, states: Sequence[NucleosideStateCandidate], tolerance: float) -> tuple[RTCoelutionStatus, float | None, tuple[float | None, ...]]:
    deltas = tuple((x.rt_apex - base.rt_apex) if x.rt_apex is not None and base.rt_apex is not None else None for x in states)
    if base.first_supporting_rt is None or base.last_supporting_rt is None or any(x.first_supporting_rt is None or x.last_supporting_rt is None for x in states):
        return RTCoelutionStatus.INSUFFICIENT_RT_EVIDENCE, None, deltas
    overlaps = []
    for x in states:
        overlap = max(0.0, min(base.last_supporting_rt, x.last_supporting_rt) - max(base.first_supporting_rt, x.first_supporting_rt))
        union = max(base.last_supporting_rt, x.last_supporting_rt) - min(base.first_supporting_rt, x.first_supporting_rt)
        overlaps.append(overlap / union if union > 0 else 1.0)
    fraction = min(overlaps, default=1.0)
    if all(delta is not None and abs(delta) <= tolerance for delta in deltas): status = RTCoelutionStatus.COELUTING
    elif fraction > 0: status = RTCoelutionStatus.PARTIALLY_OVERLAPPING
    else: status = RTCoelutionStatus.DISTINCT_RETENTION
    return status, fraction, deltas


def build_nucleoside_state_families(matches: Sequence[NucleosidePeakMatch] | None = None, *, state_delta_definitions: Sequence[StateDeltaDefinition] | None = None, state_candidates: Sequence[NucleosideStateCandidate] | None = None, parameters: P1APAuditParameters | None = None) -> tuple[NucleosideStateFamily, ...]:
    del matches
    params = parameters or P1APAuditParameters(); definitions = tuple(state_delta_definitions or build_default_state_delta_definitions())
    expected = {x.state_label for x in definitions}; groups: dict[tuple[str, str, int], list[NucleosideStateCandidate]] = defaultdict(list)
    for row in state_candidates or (): groups[(row.candidate_id, row.adduct_type, row.charge)].append(row)
    output = []
    for key, rows in sorted(groups.items()):
        labels = {x.state_label for x in rows}
        if StateLabel.BASE_STATE not in labels or len(labels) < 2: continue
        chosen = [min((x for x in rows if x.state_label is label), key=lambda x: (abs(x.delta_error_da), -x.scan_recurrence_fraction, x.intensity_rank, x.observed_peak_id)) for label in sorted(labels, key=lambda x: x.value)]
        base = next(x for x in chosen if x.state_label is StateLabel.BASE_STATE); ordered = sorted(chosen, key=lambda x: (x.target_neutral_delta, x.observed_peak_id))
        identity_ambiguous = any(x.distinct_identity_count_for_peak > 1 for x in rows); adduct_ambiguous = any(x.distinct_adduct_count_for_peak > 1 for x in rows)
        isotope = any(x.isotope_compatibility == "POSSIBLE" for x in ordered); alternatives = sum(x.alternative_explanation_count for x in ordered)
        rt_status, overlap, delta_rts = _rt_relationship(base, ordered, params.rt_coelution_apex_tolerance_min)
        blocks = [b for x in ordered for b in x.state_block_reasons]
        if StateLabel.PLUS_32_EQUIVALENT in labels and StateLabel.PLUS_16_EQUIVALENT not in labels: blocks.append("INCOMPLETE_STATE_SERIES")
        if StateLabel.PLUS_34_EQUIVALENT in labels and StateLabel.PLUS_18_EQUIVALENT not in labels: blocks.append("INCOMPLETE_STATE_SERIES")
        if identity_ambiguous or adduct_ambiguous or isotope: blocks.append("STATE_ASSIGNMENT_AMBIGUITY")
        if rt_status is RTCoelutionStatus.INSUFFICIENT_RT_EVIDENCE: blocks.append("INSUFFICIENT_RT_EVIDENCE")
        low = any(x.scan_recurrence_fraction < params.minimum_recurrence_fraction or x.fwhm is None or x.fwhm <= 0 for x in ordered)
        if identity_ambiguous or adduct_ambiguous or isotope: quality, confidence = StateFamilyQuality.AMBIGUOUS_STATE_FAMILY, "LOW"
        elif low: quality, confidence = StateFamilyQuality.LOW_QUALITY_STATE_FAMILY, "LOW"
        elif "APEX_CENTROID_DISAGREEMENT" in blocks: quality, confidence = StateFamilyQuality.SUPPORTIVE_STATE_FAMILY, "MEDIUM"
        elif all(x.match_status is MatchStatus.STRICT and x.scan_recurrence_fraction >= params.high_recurrence_fraction for x in ordered): quality, confidence = StateFamilyQuality.HIGH_QUALITY_STATE_FAMILY, "HIGH"
        else: quality, confidence = StateFamilyQuality.SUPPORTIVE_STATE_FAMILY, "MEDIUM"
        output.append(NucleosideStateFamily(
            state_family_id=_id("P1APFAMILY", f"{key}|{'|'.join(x.observed_peak_id for x in ordered)}"),
            base_candidate_id=base.candidate_id, base_candidate_name=base.candidate_name, parent_base=base.parent_base,
            ion_mode=base.ion_mode, charge=base.charge, adduct_type=base.adduct_type,
            base_peak_id=base.observed_peak_id, base_observed_mz=base.observed_apex_mz,
            detected_state_labels=tuple(x.state_label for x in ordered), detected_state_count=len(ordered),
            observed_neutral_deltas=tuple(x.observed_neutral_delta for x in ordered), expected_neutral_deltas=tuple(x.target_neutral_delta for x in ordered), delta_errors=tuple(x.delta_error_da for x in ordered),
            base_rt_apex=base.rt_apex, state_rt_apices=tuple(x.rt_apex for x in ordered), delta_rts=delta_rts,
            rt_overlap_fraction=overlap, rt_coelution_status=rt_status,
            distinct_chromatographic_feature_status="NOT_CHEMICALLY_INTERPRETED_" + rt_status.value,
            missing_expected_states=tuple(sorted(expected - labels, key=lambda x: x.value)), extra_unresolved_states=(),
            identity_ambiguity_status="IDENTITY_AMBIGUOUS" if identity_ambiguous else "UNAMBIGUOUS",
            adduct_ambiguity_status="ADDUCT_AMBIGUOUS" if adduct_ambiguous else "UNAMBIGUOUS",
            isotope_compatibility="POSSIBLE" if isotope else "NOT_SUPPORTED_BY_PROFILE_ANNOTATION",
            known_adduct_compatibility="POSSIBLE" if adduct_ambiguous else "NO_ALTERNATIVE_CONFIGURED_ADDUCT_MATCH",
            charge_scaled_spacing_compatibility=True, alternative_explanation_count=alternatives,
            alternative_explanation_status="ALTERNATIVES_PRESENT" if alternatives else "NO_ENUMERATED_ALTERNATIVE",
            series_quality_status=quality, series_confidence=confidence, series_block_reasons=_blocks(blocks),
        ))
    return tuple(sorted(output, key=lambda x: (-x.detected_state_count, x.series_quality_status.value, x.base_candidate_id, x.adduct_type, x.base_observed_mz)))


def _t1_status(t1_result: Any | None) -> str:
    if t1_result is None: return "INSUFFICIENT_CROSS_DIGEST_EVIDENCE"
    if isinstance(t1_result, str): return t1_result
    summary = getattr(t1_result, "summary", t1_result)
    count = getattr(summary, "state_family_count", None)
    if count is None and isinstance(summary, Mapping): count = summary.get("state_family_count", summary.get("State_Family_Count"))
    return "T1_SERIES_NOT_OBSERVED" if count == 0 else "T1_STATE_EVIDENCE_PRESENT" if count else "INSUFFICIENT_CROSS_DIGEST_EVIDENCE"


def reconcile_p1ap_with_t1_and_full_length(state_families: Sequence[NucleosideStateFamily], *, t1_result: Any | None = None, full_length_series: Sequence[float] | None = None) -> tuple[CrossDigestReconciliation, ...]:
    normalized = () if not full_length_series else tuple(float(x) - float(full_length_series[0]) for x in full_length_series); t1 = _t1_status(t1_result)
    families: Sequence[NucleosideStateFamily | None] = state_families or (None,)
    output = []
    for family in families:
        blocks = []; labels = set(family.detected_state_labels) if family else set()
        if family is None: p1 = "NO_P1AP_STATE_EVIDENCE"
        elif t1 == "T1_SERIES_NOT_OBSERVED": p1 = "P1AP_STATE_EVIDENCE_WITHOUT_T1_SERIES"
        elif t1 == "T1_STATE_EVIDENCE_PRESENT": p1 = "P1AP_AND_T1_PARTIALLY_COMPATIBLE"
        else: p1 = "P1AP_STATE_EVIDENCE_T1_LOCALIZATION_UNSUPPORTED"; blocks.append("T1_STATE_RESULT_MISSING")
        if not normalized: full = "INSUFFICIENT_P1AP_EVIDENCE"; blocks.append("FULL_LENGTH_SERIES_RESULT_MISSING")
        elif family is None: full = "P1AP_STATE_PATTERN_NOT_OBSERVED"
        elif family.identity_ambiguity_status != "UNAMBIGUOUS": full = "AMBIGUOUS_NUCLEOSIDE_IDENTITY"; blocks.append("AMBIGUOUS_NUCLEOSIDE_IDENTITY")
        elif {StateLabel.PLUS_18_EQUIVALENT, StateLabel.PLUS_34_EQUIVALENT} <= labels: full = "FULL_LENGTH_DELTA_PATTERN_COMPATIBLE"
        elif StateLabel.PLUS_18_EQUIVALENT in labels and labels & {StateLabel.PLUS_16_EQUIVALENT, StateLabel.PLUS_34_EQUIVALENT}: full = "PARTIALLY_COMPATIBLE_WITH_FULL_LENGTH_PATTERN"
        elif StateLabel.PLUS_16_EQUIVALENT in labels: full = "P1AP_PLUS16_ONLY"
        elif StateLabel.PLUS_18_EQUIVALENT in labels: full = "P1AP_PLUS18_ONLY"
        else: full = "P1AP_STATE_PATTERN_NOT_OBSERVED"
        output.append(CrossDigestReconciliation(
            reconciliation_id=_id("P1APRECON", family.state_family_id if family else "NO_STATE_FAMILY"),
            state_family_id=family.state_family_id if family else "NO_STATE_FAMILY", p1ap_state_status=p1,
            t1_reconciliation_status=t1, full_length_reconciliation_status=full,
            full_length_normalized_deltas=normalized, localization_status="NUCLEOSIDE_CLASS_ONLY_T1_FRAGMENT_LOCALIZATION_UNSUPPORTED",
            reconciliation_block_reasons=_blocks(blocks),
        ))
    return tuple(output)


def _source_blocks(metadata: MzMLSourceMetadataRecord | None, profile: ReplicateRunPeakProfile) -> list[str]:
    blocks = list(profile.block_reasons)
    if metadata is None: blocks.append("SOURCE_METADATA_RECORD_MISSING")
    else:
        if metadata.context_source != "USER_PROVIDED_RUNTIME_MANIFEST": blocks.append("USER_MANIFEST_CONTEXT_MISSING")
        if metadata.polarity_status is PolarityStatus.MIXED_POLARITY: blocks.append("MIXED_POLARITY_INPUT")
        elif metadata.polarity_status is PolarityStatus.NOT_RECORDED: blocks.append("MISSING_POLARITY_METADATA")
        elif metadata.polarity_status is not PolarityStatus.POSITIVE_ONLY: blocks.append("SOURCE_POLARITY_NOT_POSITIVE")
        if metadata.representation_status is RepresentationStatus.NOT_RECORDED: blocks.append("MISSING_REPRESENTATION_METADATA")
        elif metadata.representation_status is not RepresentationStatus.PROFILE_ONLY: blocks.append("REPRESENTATION_NOT_PROFILE")
    if not profile.ms1_spectra_used: blocks.append("NO_MS1_SPECTRA")
    if not profile.peaks: blocks.append("NO_DETECTED_PEAKS")
    return blocks


def audit_p1ap_nucleoside_state_series(mzml_path: Path, *, project_root: Path, sequence: str, modification_registry: Sequence[Any] | None = None, source_metadata_record: MzMLSourceMetadataRecord | None = None, runtime_context: Mapping[str, Any] | None = None, t1_result: Any | None = None, full_length_series: Sequence[float] | None = None, ion_config: Mapping[str, Any] | None = None, structure_constraints: Sequence[Mapping[str, Any]] | None = None, detection_config: Mapping[str, Any] | None = None, parameters: P1APAuditParameters | None = None) -> P1APNucleosideStateAuditResult:
    params = parameters or P1APAuditParameters(); params.validate(); runtime = dict(runtime_context or {})
    candidates = generate_nucleoside_candidates(project_root=project_root, sequence=sequence, modification_registry=modification_registry, structure_constraints=structure_constraints)
    ions = generate_positive_nucleoside_ion_hypotheses(candidates, ion_config=ion_config)
    definitions = build_default_state_delta_definitions(); targets = _targets(ions, definitions)
    run_profile = build_p1ap_ms1_peak_profile(mzml_path, source_metadata_record=source_metadata_record, detection_config=detection_config, rt_targets=targets, ms2_ion_hypotheses=ions, parameters=params)
    profile = run_profile.profile; blocks = _source_blocks(source_metadata_record, profile); selected = _selected_peaks(profile)
    if not candidates: blocks.append("NO_THEORETICAL_NUCLEOSIDES")
    if not ions: blocks.append("NO_ION_HYPOTHESES")
    states = _state_candidates(ions, _state_search_peaks(profile), definitions, run_profile.rt_evidence, params)
    neutral_masses = {x.candidate_id: x.theoretical_neutral_mass for x in candidates if x.theoretical_neutral_mass is not None}
    matches = match_nucleoside_ions_to_peaks(ions, selected, matching_config={"parameters": params}, rt_evidence=run_profile.rt_evidence, ms2_counts=run_profile.compatible_ms2_counts, nearest_ms2=run_profile.nearest_ms2_precursor, neutral_masses=neutral_masses)
    families = build_nucleoside_state_families(state_candidates=states, state_delta_definitions=definitions, parameters=params)
    reconciliations = reconcile_p1ap_with_t1_and_full_length(families, t1_result=t1_result, full_length_series=full_length_series)
    if not matches: blocks.append("NO_CANDIDATE_MATCHES")
    if not families: blocks.append("NO_STATE_SERIES")
    metadata = source_metadata_record; summary_recon = reconciliations[0]
    run_summary = P1APRunSummary(
        source_id=profile.run_label, rna_identity=metadata.rna_identity if metadata else runtime.get("RNA_Identity", "UNKNOWN"),
        digest_type=metadata.digest_type if metadata else runtime.get("Digest_Type", "UNKNOWN"), digest_detail=runtime.get("Digest_Detail", "UNKNOWN"),
        context_source=metadata.context_source if metadata else runtime.get("Context_Source", "UNKNOWN"), context_confidence=metadata.context_confidence if metadata else runtime.get("Context_Confidence", "UNKNOWN"),
        input_path=str(mzml_path), status="BLOCKED" if any(x in blocks for x in ("NO_MS1_SPECTRA", "SOURCE_POLARITY_NOT_POSITIVE", "MIXED_POLARITY_INPUT", "MISSING_POLARITY_METADATA", "REPRESENTATION_NOT_PROFILE", "MISSING_REPRESENTATION_METADATA")) else "COMPLETED",
        ms1_spectra_used=profile.ms1_spectra_used, ms2_spectra_present=profile.ms2_spectra_excluded, ms2_spectra_excluded=profile.ms2_spectra_excluded,
        unknown_ms_level_excluded=profile.missing_ms_level_spectra, aggregation_method=profile.aggregation_method, grid_method=profile.mz_grid_method,
        per_scan_normalization=profile.intensity_normalization_method, smoothing_method=profile.smoothing_method, baseline_method=profile.baseline_method,
        peak_detection_method=profile.peak_detection_method, detected_peak_count=profile.detected_peak_count, selected_peak_count=len(selected),
        rt_evidence_method="TARGETED_SCAN_LEVEL_XIC_SUPPORT_WITHOUT_LC_IDENTITY_ASSIGNMENT", block_reasons=_blocks(blocks),
    )
    canonical = sum(x.candidate_class is NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE for x in candidates)
    modified = sum(x.candidate_class in {NucleosideCandidateClass.MODIFIED_NEUTRAL_NUCLEOSIDE, NucleosideCandidateClass.MASS_ONLY_MODIFIED_NUCLEOSIDE} for x in candidates)
    ambiguous_family = sum(x.series_quality_status is StateFamilyQuality.AMBIGUOUS_STATE_FAMILY for x in families)
    summary = P1APNucleosideStateSummary(
        canonical_nucleoside_hypothesis_count=canonical, modified_nucleoside_hypothesis_count=modified,
        structure_blocked_candidate_count=sum(x.candidate_class is NucleosideCandidateClass.STRUCTURALLY_BLOCKED for x in candidates), ion_hypothesis_count=len(ions),
        candidate_match_count=len(matches), strict_match_count=sum(x.apex_match_status is MatchStatus.STRICT for x in matches), supportive_match_count=sum(x.apex_match_status is MatchStatus.SUPPORTIVE for x in matches),
        unique_matched_peak_count=len({x.observed_peak_id for x in matches}), unambiguous_match_count=sum(x.identity_ambiguity_status == x.adduct_ambiguity_status == "UNAMBIGUOUS" for x in matches),
        identity_ambiguous_match_count=sum(x.identity_ambiguity_status != "UNAMBIGUOUS" for x in matches), adduct_ambiguous_match_count=sum(x.adduct_ambiguity_status != "UNAMBIGUOUS" for x in matches),
        state_family_count=len(families), plus16_family_count=sum(StateLabel.PLUS_16_EQUIVALENT in x.detected_state_labels for x in families), plus18_family_count=sum(StateLabel.PLUS_18_EQUIVALENT in x.detected_state_labels for x in families),
        plus32_family_count=sum(StateLabel.PLUS_32_EQUIVALENT in x.detected_state_labels for x in families), plus34_family_count=sum(StateLabel.PLUS_34_EQUIVALENT in x.detected_state_labels for x in families),
        high_quality_state_family_count=sum(x.series_quality_status is StateFamilyQuality.HIGH_QUALITY_STATE_FAMILY for x in families), supportive_state_family_count=sum(x.series_quality_status is StateFamilyQuality.SUPPORTIVE_STATE_FAMILY for x in families), ambiguous_state_family_count=ambiguous_family,
        ms2_available_candidate_count=len({x.candidate_id for x in matches if x.compatible_ms2_spectrum_count > 0}), median_absolute_mass_error_mz=median(x.absolute_delta_mz for x in matches) if matches else None,
        full_length_reconciliation_status=summary_recon.full_length_reconciliation_status, t1_reconciliation_status=summary_recon.t1_reconciliation_status,
        overall_evidence_status="P1AP_STATE_EVIDENCE_DETECTED" if families else "NO_P1AP_STATE_EVIDENCE",
        overall_confidence="LOW" if ambiguous_family or not families else "HIGH" if all(x.series_quality_status is StateFamilyQuality.HIGH_QUALITY_STATE_FAMILY for x in families) else "MEDIUM",
        overall_block_reasons=_blocks(blocks + ["CHEMICAL_IDENTITY_UNSUPPORTED"]),
    )
    return P1APNucleosideStateAuditResult(params, run_profile, run_summary, candidates, ions, matches, states, families, reconciliations, summary)


def _record(value: Any) -> dict[str, Any]:
    row = asdict(value); row.pop("source_candidate", None)
    def normalize(item: Any) -> Any:
        if isinstance(item, Enum): return item.value
        if isinstance(item, dict): return {k: normalize(v) for k, v in item.items()}
        if isinstance(item, (tuple, list)): return [normalize(v) for v in item]
        return item
    return normalize(row)


def audit_optional_result(result: P1APNucleosideStateAuditResult) -> dict[str, Any]:
    safeguards = _record(P1APSafeguards())
    by_peak: dict[str, list[NucleosideStateCandidate]] = defaultdict(list)
    for state in result.state_candidates: by_peak[state.observed_peak_id].append(state)
    peaks = []
    for peak in result.run_profile.profile.peaks:
        candidates = sorted(by_peak.get(peak.peak_id, ()), key=lambda x: (x.state_label is not StateLabel.BASE_STATE, abs(x.delta_error_da), x.candidate_id))
        state = candidates[0] if candidates else None
        peaks.append({**_record(peak), **safeguards,
            "rt_apex": state.rt_apex if state else None, "rt_centroid": state.rt_centroid if state else None,
            "rt_span": state.rt_span if state else None,
            "rt_profile_status": "TARGETED_SCAN_LEVEL_RT_EVIDENCE" if state and state.rt_apex is not None else "NOT_RECORDED_UNTARGETED_PROFILE_PEAK"})
    return {
        "run_summary_records": [_record(result.run_summary)], "peak_records": peaks,
        "candidate_records": [_record(x) for x in result.candidates], "match_records": [_record(x) for x in result.matches],
        "state_family_records": [_record(x) for x in result.state_families], "reconciliation_records": [_record(x) for x in result.reconciliations],
        "summary_records": [_record(result.summary)],
    }
