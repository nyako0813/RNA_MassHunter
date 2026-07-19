"""Shadow-only MS2 identity audit for positive-mode P1/AP nucleosides.

The audit deliberately separates precursor compatibility, selective binary decode,
formula-safe product generation, and product matching.  It never assigns an exact
chemical/isomer identity and never propagates into formal scoring or localization.
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field, replace
from enum import Enum
from math import isfinite
from pathlib import Path
from statistics import median
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from pyteomics import mzml

from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.masses import PROTON_MASS
from rna_masshunter.p1_sap_chemical_state_audit import MODEL_NOT_DEFINED
from rna_masshunter.p1_sap_feature_quality import _parse_composition
from rna_masshunter.sciex_p1ap_nucleoside_state_audit import (
    MatchStatus, NucleosideCandidateClass, P1APNucleosideStateAuditResult,
)
from rna_masshunter.sciex_t1_profile_peak_audit import (
    T1PeakDetectionParameters, detect_t1_profile_peaks,
)
from rna_masshunter.sciex_t1_replicate_consistency_audit import _decode_array, _rt_minutes

OPTIONAL_RESULT_KEY = "sciex_p1ap_nucleoside_ms2_identity_audit"
ALGORITHM_VERSION = "sciex-p1ap-nucleoside-ms2-identity-audit-v1"

_BLOCK_ORDER = (
    "INPUT_FILE_NOT_FOUND", "INPUT_FILE_UNREADABLE", "SOURCE_METADATA_RECORD_MISSING",
    "P1AP_MS1_AUDIT_RESULT_MISSING", "NO_MS2_SPECTRA", "NO_PRECURSOR_COMPATIBLE_MS2",
    "MISSING_SELECTED_ION_MZ", "MISSING_ISOLATION_WINDOW", "MISSING_COLLISION_ENERGY",
    "PRECURSOR_OUTSIDE_ISOLATION_WINDOW", "MS2_BINARY_DECODE_FAILED",
    "NO_POSITIVE_INTENSITY_MS2_PEAKS", "MS2_PEAK_DETECTION_FAILED",
    "PRODUCT_ION_REGISTRY_MISSING", "NO_THEORETICAL_PRODUCT_IONS",
    "MASS_ONLY_CANDIDATE_NO_STRUCTURE_RULES", "NO_MATCHED_PRODUCT_IONS",
    "LOW_EXPLAINED_INTENSITY", "LOW_TOP_PEAK_COVERAGE", "PRODUCT_ION_AMBIGUITY",
    "SHARED_PRODUCT_IONS_ONLY", "NO_CANDIDATE_UNIQUE_PRODUCT_IONS",
    "INSUFFICIENT_MS2_RECURRENCE", "COLLISION_ENERGY_HETEROGENEITY",
    "MS2_NONDISCRIMINATING", "MS2_CONFLICTING",
    "EXACT_ISOMER_DISCRIMINATION_UNSUPPORTED", "CHEMICAL_IDENTITY_UNSUPPORTED",
)


def _blocks(values: Iterable[str]) -> tuple[str, ...]:
    found = {str(x) for x in values if x}
    order = {value: index for index, value in enumerate(_BLOCK_ORDER)}
    return tuple(sorted(found, key=lambda x: (order.get(x, len(order)), x)))


@dataclass(frozen=True, kw_only=True)
class MS2Safeguards:
    shadow_analysis_only: bool = True
    ms2_evidence_only: bool = True
    formal_propagation: bool = False
    chemical_identity_assigned: bool = False
    modification_assigned: bool = False
    exact_candidate_identity_confirmed: bool = False
    exact_isomer_identity_confirmed: bool = False
    exact_nucleotide_localization: bool = False
    exact_atom_localization: bool = False
    reaction_order_assigned: bool = False
    applied_to_formal_score: bool = False
    applied_to_ranking: bool = False
    applied_to_candidate_filtering: bool = False
    applied_to_final_consensus: bool = False


@dataclass(frozen=True, kw_only=True)
class MS2CandidateRecord(MS2Safeguards):
    candidate_id: str
    candidate_name: str
    parent_base: str
    candidate_class: str
    molecular_formula: str
    theoretical_neutral_mass: float | None
    theoretical_precursor_mz: float
    observed_ms1_peak_id: str
    observed_ms1_mz: float
    ms1_mass_error: float
    ms1_identity_ambiguity_status: str
    ms1_state_family_id: str
    ion_mode: str = "POSITIVE"
    charge: int = 1
    adduct_type: str = "[M+H]+"


@dataclass(frozen=True, kw_only=True)
class PrecursorCompatibleMS2Record(MS2Safeguards):
    candidate_id: str
    candidate_name: str
    ms2_spectrum_id: str
    ms2_scan_time: float | None
    selected_ion_mz: float | None
    isolation_target_mz: float | None
    isolation_lower_offset: float | None
    isolation_upper_offset: float | None
    isolation_lower_bound: float | None
    isolation_upper_bound: float | None
    candidate_theoretical_mz: float
    candidate_observed_ms1_mz: float
    selected_ion_delta: float | None
    isolation_contains_theoretical_mz: bool | None
    isolation_contains_observed_ms1_mz: bool | None
    precursor_charge: int | None
    collision_energy: float | None
    collision_energy_unit: str
    default_array_length: int | None
    precursor_compatibility_status: str
    precursor_compatibility_confidence: str
    precursor_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class ProcessedMS2Spectrum(MS2Safeguards):
    ms2_spectrum_id: str
    scan_time: float | None
    raw_peak_count: int
    positive_intensity_peak_count: int
    zero_intensity_peak_count: int
    negative_intensity_peak_count: int
    filtered_peak_count: int
    base_peak_mz: float | None
    base_peak_intensity: float
    tic: float
    mz_min: float | None
    mz_max: float | None
    profile_or_centroid_metadata: str
    ms2_preprocessing_status: str
    ms2_preprocessing_block_reasons: tuple[str, ...]
    peaks: tuple[tuple[float, float], ...] = field(repr=False, compare=False)


@dataclass(frozen=True, kw_only=True)
class NucleosideProductIonHypothesis(MS2Safeguards):
    product_ion_id: str
    candidate_id: str
    product_ion_label: str
    product_ion_class: str
    theoretical_product_mz: float
    product_formula: str
    neutral_loss_formula: str
    neutral_loss_mass: float
    rule_provenance: str
    product_ion_status: str
    product_ion_block_reasons: tuple[str, ...]
    candidate_specific: bool
    shared_candidate_ids: tuple[str, ...] = ()


@dataclass(frozen=True, kw_only=True)
class CandidatePairDiscrimination(MS2Safeguards):
    candidate_1_id: str
    candidate_2_id: str
    shared_product_ion_count: int
    candidate_1_unique_product_ion_count: int
    candidate_2_unique_product_ion_count: int
    formula_identical_status: str
    structure_rules_available: bool
    theoretical_discrimination_possible: bool
    discrimination_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class NucleosideMS2ProductMatch(MS2Safeguards):
    candidate_id: str
    ms2_spectrum_id: str
    product_ion_id: str
    product_ion_label: str
    product_ion_class: str
    theoretical_product_mz: float
    observed_product_mz: float
    delta_mz: float
    absolute_delta_mz: float
    ppm_error: float
    observed_intensity: float
    relative_intensity: float
    intensity_rank: int
    candidate_count_for_observed_peak: int
    observed_peak_count_for_product_ion: int
    match_ambiguity_status: str
    match_quality_status: str
    match_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class SpectrumMS2Summary(MS2Safeguards):
    candidate_id: str
    ms2_spectrum_id: str
    collision_energy: float | None
    precursor_compatibility_status: str
    theoretical_product_ion_count: int
    matched_product_ion_count: int
    unique_matched_product_ion_count: int
    shared_matched_product_ion_count: int
    diagnostic_product_ion_count: int
    base_related_ion_matched: bool
    ribose_loss_ion_matched: bool
    modification_specific_ion_matched: bool
    explained_intensity_fraction: float
    top_10_peak_explained_fraction: float
    median_product_mass_error: float | None
    product_ion_ambiguity_count: int
    spectrum_evidence_status: str
    spectrum_evidence_confidence: str
    spectrum_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class CollisionEnergySummary(MS2Safeguards):
    candidate_id: str
    collision_energy_value: float | None
    collision_energy_unit: str
    collision_energy_status: str
    spectrum_count_per_energy: int
    matched_product_ion_count_per_energy: int
    diagnostic_ion_recurrence_per_energy: float | None


@dataclass(frozen=True, kw_only=True)
class ProductIonRecurrence(MS2Safeguards):
    candidate_id: str
    product_ion_id: str
    supporting_ms2_spectrum_count: int
    compatible_ms2_spectrum_count: int
    ms2_recurrence_fraction: float
    collision_energy_set: tuple[str, ...]
    first_supporting_rt: float | None
    last_supporting_rt: float | None
    product_ion_recurrence_status: str


@dataclass(frozen=True, kw_only=True)
class NucleosideCandidateMS2Summary(MS2Safeguards):
    candidate_id: str
    candidate_name: str
    compatible_ms2_spectrum_count: int
    usable_ms2_spectrum_count: int
    collision_energy_count: int
    theoretical_product_ion_count: int
    recurrent_matched_product_ion_count: int
    candidate_unique_recurrent_ion_count: int
    shared_recurrent_ion_count: int
    best_spectrum_id: str
    best_spectrum_evidence_status: str
    median_explained_intensity_fraction: float | None
    median_top10_explained_fraction: float | None
    median_mass_error: float | None
    ms2_identity_evidence_status: str
    ms2_identity_confidence: str
    identity_ambiguity_before_ms2: str
    identity_ambiguity_after_ms2: str
    candidate_specific_ms2_rules_available: bool
    ms2_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class P1APMS1MS2Reconciliation(MS2Safeguards):
    reconciliation_id: str
    p1ap_ms1_state_family_id: str
    p1ap_ms2_status: str
    t1_reconciliation_status: str
    full_length_reconciliation_status: str
    state_interpretation_resolved: bool
    reconciliation_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class P1APNucleosideMS2Summary(MS2Safeguards):
    source_id: str
    status: str
    total_spectra_seen: int
    total_ms2_metadata_records: int
    precursor_compatible_record_count: int
    decoded_unique_spectrum_count: int
    candidate_count: int
    theoretical_product_ion_count: int
    product_match_count: int
    overall_evidence_status: str
    overall_confidence: str
    overall_block_reasons: tuple[str, ...]
    runtime_seconds: float


@dataclass(frozen=True)
class P1APNucleosideMS2AuditResult:
    candidates: tuple[MS2CandidateRecord, ...]
    precursor_records: tuple[PrecursorCompatibleMS2Record, ...]
    spectrum_records: tuple[ProcessedMS2Spectrum, ...]
    theoretical_product_records: tuple[NucleosideProductIonHypothesis, ...]
    discrimination_records: tuple[CandidatePairDiscrimination, ...]
    product_match_records: tuple[NucleosideMS2ProductMatch, ...]
    spectrum_summary_records: tuple[SpectrumMS2Summary, ...]
    collision_energy_records: tuple[CollisionEnergySummary, ...]
    recurrence_records: tuple[ProductIonRecurrence, ...]
    candidate_summary_records: tuple[NucleosideCandidateMS2Summary, ...]
    reconciliation_records: tuple[P1APMS1MS2Reconciliation, ...]
    summary: P1APNucleosideMS2Summary
    algorithm_version: str = ALGORITHM_VERSION
    formal_propagation: bool = False


def _number(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return result if isfinite(result) else None


def _integer(value: Any) -> int | None:
    number = _number(value)
    return int(number) if number is not None else None


def _polarity(spectrum: Mapping[str, Any]) -> str:
    if "positive scan" in spectrum:
        return "POSITIVE"
    if "negative scan" in spectrum:
        return "NEGATIVE"
    return "UNKNOWN"


def _representation(spectrum: Mapping[str, Any]) -> str:
    if "profile spectrum" in spectrum:
        return "PROFILE"
    if "centroid spectrum" in spectrum:
        return "CENTROID"
    return "UNKNOWN"


def _precursor_metadata(spectrum: Mapping[str, Any]) -> dict[str, Any]:
    precursors = ((spectrum.get("precursorList") or {}).get("precursor") or [])
    precursor = precursors[0] if precursors else {}
    ions = ((precursor.get("selectedIonList") or {}).get("selectedIon") or [])
    ion = ions[0] if ions else {}
    isolation = precursor.get("isolationWindow") or {}
    activation = precursor.get("activation") or {}
    target = _number(isolation.get("isolation window target m/z"))
    lower = _number(isolation.get("isolation window lower offset"))
    upper = _number(isolation.get("isolation window upper offset"))
    energy = _number(activation.get("collision energy"))
    energy_unit = str(activation.get("unitName") or activation.get("unit symbol") or "NOT_RECORDED")
    return {
        "selected": _number(ion.get("selected ion m/z")),
        "charge": _integer(ion.get("charge state")),
        "target": target, "lower": lower, "upper": upper,
        "lower_bound": target - lower if target is not None and lower is not None else None,
        "upper_bound": target + upper if target is not None and upper is not None else None,
        "energy": energy, "energy_unit": energy_unit,
    }


def candidate_records_from_ms1_result(result: P1APNucleosideStateAuditResult) -> tuple[MS2CandidateRecord, ...]:
    """Adapt matches plus unmatched canonical BASE-state alternatives from the MS1 audit."""
    candidates = {x.candidate_id: x for x in result.candidates}
    family_by_candidate: dict[str, str] = {}
    for family in result.state_families:
        family_by_candidate.setdefault(family.base_candidate_id, family.state_family_id)
    rows: dict[tuple[str, str], MS2CandidateRecord] = {}
    for match in result.matches:
        candidate = candidates[match.candidate_id]
        rows[(match.candidate_id, match.observed_peak_id)] = MS2CandidateRecord(
            candidate_id=match.candidate_id, candidate_name=match.candidate_name,
            parent_base=match.parent_base, candidate_class=match.candidate_class.value,
            molecular_formula=candidate.molecular_formula,
            theoretical_neutral_mass=candidate.theoretical_neutral_mass,
            theoretical_precursor_mz=match.theoretical_mz,
            observed_ms1_peak_id=match.observed_peak_id, observed_ms1_mz=match.observed_apex_mz,
            ms1_mass_error=match.delta_mz,
            ms1_identity_ambiguity_status=match.identity_ambiguity_status,
            ms1_state_family_id=family_by_candidate.get(match.candidate_id, ""),
            ion_mode=match.ion_mode, charge=match.charge, adduct_type=match.adduct_type,
        )
    # A canonical BASE-state alternative can be diagnostically present as an isotope/state peak
    # without appearing in the primary-match table (e.g. canonical G in an ambiguous +16 region).
    for state in result.state_candidates:
        if state.state_label.value != "BASE_STATE" or state.candidate_class is not NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE:
            continue
        key = (state.candidate_id, state.observed_peak_id)
        if key in rows:
            continue
        candidate = candidates[state.candidate_id]
        rows[key] = MS2CandidateRecord(
            candidate_id=state.candidate_id, candidate_name=state.candidate_name,
            parent_base=state.parent_base, candidate_class=state.candidate_class.value,
            molecular_formula=candidate.molecular_formula,
            theoretical_neutral_mass=candidate.theoretical_neutral_mass,
            theoretical_precursor_mz=state.expected_mz,
            observed_ms1_peak_id=state.observed_peak_id, observed_ms1_mz=state.observed_apex_mz,
            ms1_mass_error=state.observed_apex_mz - state.expected_mz,
            ms1_identity_ambiguity_status=("IDENTITY_AMBIGUOUS" if "IDENTITY_AMBIGUITY" in state.state_block_reasons else "UNAMBIGUOUS"),
            ms1_state_family_id=family_by_candidate.get(state.candidate_id, ""),
            ion_mode=state.ion_mode, charge=state.charge, adduct_type=state.adduct_type,
        )
    return tuple(sorted(rows.values(), key=lambda x: (x.theoretical_precursor_mz, x.candidate_id, x.observed_ms1_peak_id)))


def _iter_metadata(path: Path, metadata_spectra: Iterable[Mapping[str, Any]] | None):
    if metadata_spectra is not None:
        yield from metadata_spectra
        return
    with mzml.MzML(str(path), decode_binary=False) as reader:
        yield from reader


def select_precursor_compatible_ms2_spectra(
    mzml_path: Path, candidate_records: Sequence[MS2CandidateRecord], *,
    selection_config: Mapping[str, Any] | None = None,
    metadata_spectra: Iterable[Mapping[str, Any]] | None = None,
) -> list[PrecursorCompatibleMS2Record]:
    config = dict(selection_config or {})
    tolerance_ppm = float(config.pop("tolerance_ppm", 20.0))
    if config:
        raise ValueError(f"unsupported selection_config keys: {sorted(config)}")
    output: list[PrecursorCompatibleMS2Record] = []
    for spectrum in _iter_metadata(Path(mzml_path), metadata_spectra):
        if _integer(spectrum.get("ms level")) != 2 or _polarity(spectrum) == "NEGATIVE":
            continue
        meta = _precursor_metadata(spectrum)
        for candidate in sorted(candidate_records, key=lambda x: (x.theoretical_precursor_mz, x.candidate_id)):
            if meta["charge"] not in (None, 0, candidate.charge):
                continue
            selected_delta = meta["selected"] - candidate.theoretical_precursor_mz if meta["selected"] is not None else None
            selected_support = selected_delta is not None and abs(selected_delta) <= candidate.theoretical_precursor_mz * tolerance_ppm / 1e6
            bounds_known = meta["lower_bound"] is not None and meta["upper_bound"] is not None
            theoretical_in = meta["lower_bound"] <= candidate.theoretical_precursor_mz <= meta["upper_bound"] if bounds_known else None
            observed_in = meta["lower_bound"] <= candidate.observed_ms1_mz <= meta["upper_bound"] if bounds_known else None
            isolation_support = bool(theoretical_in or observed_in)
            if not selected_support and not isolation_support:
                continue
            blocks: list[str] = []
            if meta["selected"] is None: blocks.append("MISSING_SELECTED_ION_MZ")
            if not bounds_known: blocks.append("MISSING_ISOLATION_WINDOW")
            elif not isolation_support: blocks.append("PRECURSOR_OUTSIDE_ISOLATION_WINDOW")
            if meta["energy"] is None: blocks.append("MISSING_COLLISION_ENERGY")
            if selected_support and isolation_support:
                status, confidence = "BOTH_SELECTED_ION_AND_ISOLATION_SUPPORT", "HIGH"
            elif selected_support:
                status, confidence = "SELECTED_ION_MZ_WITHIN_TOLERANCE", "LOW" if not bounds_known else "MEDIUM"
            else:
                status, confidence = "ISOLATION_WINDOW_CONTAINS_CANDIDATE", "MEDIUM"
            output.append(PrecursorCompatibleMS2Record(
                candidate_id=candidate.candidate_id, candidate_name=candidate.candidate_name,
                ms2_spectrum_id=str(spectrum.get("id") or ""), ms2_scan_time=_rt_minutes(spectrum),
                selected_ion_mz=meta["selected"], isolation_target_mz=meta["target"],
                isolation_lower_offset=meta["lower"], isolation_upper_offset=meta["upper"],
                isolation_lower_bound=meta["lower_bound"], isolation_upper_bound=meta["upper_bound"],
                candidate_theoretical_mz=candidate.theoretical_precursor_mz,
                candidate_observed_ms1_mz=candidate.observed_ms1_mz,
                selected_ion_delta=selected_delta,
                isolation_contains_theoretical_mz=theoretical_in,
                isolation_contains_observed_ms1_mz=observed_in,
                precursor_charge=meta["charge"], collision_energy=meta["energy"],
                collision_energy_unit=meta["energy_unit"],
                default_array_length=_integer(spectrum.get("defaultArrayLength")),
                precursor_compatibility_status=status,
                precursor_compatibility_confidence=confidence,
                precursor_block_reasons=_blocks(blocks),
            ))
    return sorted(output, key=lambda x: (x.ms2_spectrum_id, x.candidate_id))


def _iter_decode(path: Path, spectrum_source: Iterable[Mapping[str, Any]] | None):
    if spectrum_source is not None:
        yield from spectrum_source
        return
    with mzml.MzML(str(path), decode_binary=False) as reader:
        yield from reader


def decode_selected_ms2_spectra(
    mzml_path: Path, selected_spectrum_ids: Sequence[str], *,
    preprocessing_config: Mapping[str, Any] | None = None,
    spectrum_source: Iterable[Mapping[str, Any]] | None = None,
) -> list[ProcessedMS2Spectrum]:
    config = dict(preprocessing_config or {})
    peak_params = config.pop("peak_detection_parameters", None)
    if isinstance(peak_params, Mapping): peak_params = T1PeakDetectionParameters(**peak_params)
    if config: raise ValueError(f"unsupported preprocessing_config keys: {sorted(config)}")
    wanted = set(map(str, selected_spectrum_ids)); output = []
    for spectrum in _iter_decode(Path(mzml_path), spectrum_source):
        spectrum_id = str(spectrum.get("id") or "")
        if spectrum_id not in wanted or _integer(spectrum.get("ms level")) != 2:
            continue
        blocks: list[str] = []
        try:
            mz_values = _decode_array(spectrum.get("m/z array"))
            intensities = _decode_array(spectrum.get("intensity array"))
        except Exception:
            output.append(ProcessedMS2Spectrum(
                ms2_spectrum_id=spectrum_id, scan_time=_rt_minutes(spectrum), raw_peak_count=0,
                positive_intensity_peak_count=0, zero_intensity_peak_count=0, negative_intensity_peak_count=0,
                filtered_peak_count=0, base_peak_mz=None, base_peak_intensity=0.0, tic=0.0,
                mz_min=None, mz_max=None, profile_or_centroid_metadata=_representation(spectrum),
                ms2_preprocessing_status="BLOCKED", ms2_preprocessing_block_reasons=("MS2_BINARY_DECODE_FAILED",), peaks=(),
            )); continue
        count = min(len(mz_values), len(intensities)); mz_values = mz_values[:count]; intensities = intensities[:count]
        finite = np.isfinite(mz_values) & np.isfinite(intensities)
        mz_values = mz_values[finite]; intensities = intensities[finite]
        positive = intensities > 0; zeros = intensities == 0; negative = intensities < 0
        positive_mz = mz_values[positive]; positive_i = intensities[positive]
        mode = _representation(spectrum); peaks: list[tuple[float, float]] = []
        if not len(positive_mz): blocks.append("NO_POSITIVE_INTENSITY_MS2_PEAKS")
        elif mode == "PROFILE" and len(positive_mz) >= 3 and np.all(np.diff(positive_mz) > 0):
            try:
                audit = detect_t1_profile_peaks(
                    positive_mz, positive_i, source_id=spectrum_id, measurement_id=spectrum_id,
                    rna_identity="UNKNOWN", parameters=peak_params,
                )
                peaks = [(p.centroid_mz if p.centroid_mz is not None else p.apex_mz, p.apex_intensity) for p in audit.selected_peaks]
                if not peaks: blocks.append("MS2_PEAK_DETECTION_FAILED")
            except ValueError:
                blocks.append("MS2_PEAK_DETECTION_FAILED")
        else:
            peaks = [(float(m), float(i)) for m, i in zip(positive_mz, positive_i, strict=False)]
        peaks.sort(key=lambda x: (x[0], -x[1]))
        base_index = int(np.argmax(positive_i)) if len(positive_i) else None
        output.append(ProcessedMS2Spectrum(
            ms2_spectrum_id=spectrum_id, scan_time=_rt_minutes(spectrum), raw_peak_count=count,
            positive_intensity_peak_count=int(np.sum(positive)), zero_intensity_peak_count=int(np.sum(zeros)),
            negative_intensity_peak_count=int(np.sum(negative)), filtered_peak_count=len(peaks),
            base_peak_mz=float(positive_mz[base_index]) if base_index is not None else None,
            base_peak_intensity=float(positive_i[base_index]) if base_index is not None else 0.0,
            tic=float(np.sum(positive_i)), mz_min=float(np.min(mz_values)) if len(mz_values) else None,
            mz_max=float(np.max(mz_values)) if len(mz_values) else None,
            profile_or_centroid_metadata=mode,
            ms2_preprocessing_status="COMPLETED" if peaks else "BLOCKED",
            ms2_preprocessing_block_reasons=_blocks(blocks), peaks=tuple(peaks),
        ))
    return sorted(output, key=lambda x: x.ms2_spectrum_id)


def _composition(text: str) -> ElementalComposition | None:
    counts, status = _parse_composition(text)
    return ElementalComposition(counts) if status == "ASSESSED" and counts else None


def _make_product(candidate: MS2CandidateRecord, label: str, ion_class: str,
                  neutral_product: ElementalComposition, neutral_loss: ElementalComposition,
                  provenance: str, candidate_specific: bool, *, add_proton: bool = True) -> NucleosideProductIonHypothesis:
    return NucleosideProductIonHypothesis(
        product_ion_id=f"{candidate.candidate_id}__{label}", candidate_id=candidate.candidate_id,
        product_ion_label=label, product_ion_class=ion_class,
        theoretical_product_mz=neutral_product.exact_mass + (PROTON_MASS if add_proton else 0.0),
        product_formula=neutral_product.canonical_string(),
        neutral_loss_formula=neutral_loss.canonical_string(), neutral_loss_mass=neutral_loss.exact_mass,
        rule_provenance=provenance, product_ion_status="ELIGIBLE",
        product_ion_block_reasons=(), candidate_specific=candidate_specific,
    )


def generate_nucleoside_product_ion_hypotheses(
    candidate_records: Sequence[MS2CandidateRecord], *,
    product_ion_registry: Mapping[str, Sequence[Mapping[str, Any]]] | None = None,
) -> list[NucleosideProductIonHypothesis]:
    output: list[NucleosideProductIonHypothesis] = []
    water = ElementalComposition({"H": 2, "O": 1})
    ribose_residue = ElementalComposition({"C": 5, "H": 8, "O": 4})
    for candidate in sorted(candidate_records, key=lambda x: (x.candidate_id, x.theoretical_precursor_mz)):
        composition = _composition(candidate.molecular_formula)
        candidate_class = candidate.candidate_class.value if isinstance(candidate.candidate_class, Enum) else str(candidate.candidate_class)
        if composition is not None and candidate_class == NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE.value:
            try:
                base = composition - ribose_residue
                output.append(_make_product(candidate, "BASE_MOLECULAR_ION", "BASE_RELATED_ION", base, ribose_residue, "FORMULA_DERIVED_RADICAL_CATION_DIAGNOSTIC", False, add_proton=False))
                output.append(_make_product(candidate, "PROTONATED_BASE", "RIBOSE_LOSS", base, ribose_residue, "FORMULA_DERIVED_CANONICAL_NUCLEOSIDE", False))
            except ValueError:
                pass
        # H2O loss is composition-conserving but non-specific; it is safe for any formula-defined candidate.
        if composition is not None:
            try:
                dehydrated = composition - water
                output.append(_make_product(candidate, "DEHYDRATED_PRECURSOR", "NEUTRAL_LOSS_H2O", dehydrated, water, "FORMULA_DERIVED", False))
            except ValueError:
                pass
        for serial, row in enumerate((product_ion_registry or {}).get(candidate.candidate_id, ()), 1):
            mz_value = _number(row.get("theoretical_product_mz"))
            if mz_value is None: continue
            output.append(NucleosideProductIonHypothesis(
                product_ion_id=str(row.get("product_ion_id") or f"{candidate.candidate_id}__CURATED_{serial}"),
                candidate_id=candidate.candidate_id, product_ion_label=str(row.get("product_ion_label") or "CURATED_PRODUCT"),
                product_ion_class=str(row.get("product_ion_class") or "FORMULA_DERIVED_PRODUCT"),
                theoretical_product_mz=mz_value, product_formula=str(row.get("product_formula") or "UNKNOWN"),
                neutral_loss_formula=str(row.get("neutral_loss_formula") or "UNKNOWN"),
                neutral_loss_mass=float(row.get("neutral_loss_mass") or 0.0),
                rule_provenance=str(row.get("rule_provenance") or "CURATED_REGISTRY"),
                product_ion_status="ELIGIBLE", product_ion_block_reasons=(),
                candidate_specific=bool(row.get("candidate_specific", True)),
            ))
    # Products at the same theoretical m/z are explicitly shared across candidates.
    final = []
    for product in output:
        shared = tuple(sorted({x.candidate_id for x in output if abs(x.theoretical_product_mz - product.theoretical_product_mz) <= max(product.theoretical_product_mz, x.theoretical_product_mz) * 5e-6}))
        final.append(replace(product, shared_candidate_ids=shared))
    return sorted(final, key=lambda x: (x.theoretical_product_mz, x.candidate_id, x.product_ion_id))


def compare_candidate_product_ions(candidate_records: Sequence[MS2CandidateRecord], products: Sequence[NucleosideProductIonHypothesis]) -> list[CandidatePairDiscrimination]:
    by_candidate = defaultdict(list)
    for product in products: by_candidate[product.candidate_id].append(product)
    rows = []
    ordered = sorted(candidate_records, key=lambda x: x.candidate_id)
    for index, left in enumerate(ordered):
        for right in ordered[index + 1:]:
            lp, rp = by_candidate[left.candidate_id], by_candidate[right.candidate_id]
            shared_left = {x.product_ion_id for x in lp if right.candidate_id in x.shared_candidate_ids}
            shared_right = {x.product_ion_id for x in rp if left.candidate_id in x.shared_candidate_ids}
            shared = min(len(shared_left), len(shared_right))
            l_unique = sum(1 for x in lp if x.candidate_specific and right.candidate_id not in x.shared_candidate_ids)
            r_unique = sum(1 for x in rp if x.candidate_specific and left.candidate_id not in x.shared_candidate_ids)
            possible = bool(l_unique or r_unique)
            blocks = [] if possible else ["NO_CANDIDATE_UNIQUE_PRODUCT_IONS", "EXACT_ISOMER_DISCRIMINATION_UNSUPPORTED"]
            rows.append(CandidatePairDiscrimination(
                candidate_1_id=left.candidate_id, candidate_2_id=right.candidate_id,
                shared_product_ion_count=shared, candidate_1_unique_product_ion_count=l_unique,
                candidate_2_unique_product_ion_count=r_unique,
                formula_identical_status="FORMULA_IDENTICAL" if left.molecular_formula == right.molecular_formula and left.molecular_formula != MODEL_NOT_DEFINED else "FORMULA_DIFFERENT_OR_UNKNOWN",
                structure_rules_available=possible, theoretical_discrimination_possible=possible,
                discrimination_block_reasons=_blocks(blocks),
            ))
    return rows


def match_product_ions_to_ms2_peaks(
    product_ions: Sequence[NucleosideProductIonHypothesis], spectra: Sequence[ProcessedMS2Spectrum], *,
    matching_config: Mapping[str, Any] | None = None,
    precursor_records: Sequence[PrecursorCompatibleMS2Record] | None = None,
) -> list[NucleosideMS2ProductMatch]:
    config = dict(matching_config or {}); tolerance_ppm = float(config.pop("tolerance_ppm", 20.0))
    if config: raise ValueError(f"unsupported matching_config keys: {sorted(config)}")
    compatible = {(x.candidate_id, x.ms2_spectrum_id) for x in precursor_records or ()}
    restrict = precursor_records is not None
    draft: list[dict[str, Any]] = []
    for spectrum in sorted(spectra, key=lambda x: x.ms2_spectrum_id):
        mzs = [x[0] for x in spectrum.peaks]
        ranks = {index: rank for rank, index in enumerate(sorted(range(len(spectrum.peaks)), key=lambda i: (-spectrum.peaks[i][1], spectrum.peaks[i][0])), 1)}
        for product in sorted(product_ions, key=lambda x: (x.theoretical_product_mz, x.candidate_id, x.product_ion_id)):
            if restrict and (product.candidate_id, spectrum.ms2_spectrum_id) not in compatible: continue
            tolerance = product.theoretical_product_mz * tolerance_ppm / 1e6
            lo = bisect_left(mzs, product.theoretical_product_mz - tolerance); hi = bisect_right(mzs, product.theoretical_product_mz + tolerance)
            if hi <= lo: continue
            peak_index = min(range(lo, hi), key=lambda i: (abs(mzs[i] - product.theoretical_product_mz), -spectrum.peaks[i][1], i))
            observed, intensity = spectrum.peaks[peak_index]; error = observed - product.theoretical_product_mz
            draft.append({"product": product, "spectrum": spectrum, "observed": observed, "intensity": intensity,
                          "peak_index": peak_index, "rank": ranks[peak_index], "error": error})
    candidates_for_peak = Counter((x["spectrum"].ms2_spectrum_id, round(x["observed"], 8), x["product"].candidate_id) for x in draft)
    distinct_candidates = defaultdict(set)
    for x in draft: distinct_candidates[(x["spectrum"].ms2_spectrum_id, round(x["observed"], 8))].add(x["product"].candidate_id)
    product_spectra = Counter((x["product"].product_ion_id, x["spectrum"].ms2_spectrum_id) for x in draft)
    product_counts = Counter(x["product"].product_ion_id for x in draft)
    output = []
    for x in draft:
        count = len(distinct_candidates[(x["spectrum"].ms2_spectrum_id, round(x["observed"], 8))])
        blocks = ["PRODUCT_ION_AMBIGUITY"] if count > 1 else []
        output.append(NucleosideMS2ProductMatch(
            candidate_id=x["product"].candidate_id, ms2_spectrum_id=x["spectrum"].ms2_spectrum_id,
            product_ion_id=x["product"].product_ion_id, product_ion_label=x["product"].product_ion_label,
            product_ion_class=x["product"].product_ion_class,
            theoretical_product_mz=x["product"].theoretical_product_mz, observed_product_mz=x["observed"],
            delta_mz=x["error"], absolute_delta_mz=abs(x["error"]), ppm_error=x["error"] / x["product"].theoretical_product_mz * 1e6,
            observed_intensity=x["intensity"], relative_intensity=x["intensity"] / x["spectrum"].base_peak_intensity if x["spectrum"].base_peak_intensity else 0.0,
            intensity_rank=x["rank"], candidate_count_for_observed_peak=count,
            observed_peak_count_for_product_ion=product_counts[x["product"].product_ion_id],
            match_ambiguity_status="SHARED_OR_AMBIGUOUS" if count > 1 or len(x["product"].shared_candidate_ids) > 1 else "CANDIDATE_SCOPED",
            match_quality_status="WITHIN_TOLERANCE", match_block_reasons=_blocks(blocks),
        ))
    return sorted(output, key=lambda x: (x.ms2_spectrum_id, x.observed_product_mz, x.candidate_id, x.product_ion_id))


def _spectrum_summaries(candidates, precursors, spectra, products, matches):
    product_by_candidate = defaultdict(list); match_by_pair = defaultdict(list)
    for x in products: product_by_candidate[x.candidate_id].append(x)
    for x in matches: match_by_pair[(x.candidate_id, x.ms2_spectrum_id)].append(x)
    spectrum_map = {x.ms2_spectrum_id: x for x in spectra}
    output = []
    for precursor in precursors:
        spectrum = spectrum_map.get(precursor.ms2_spectrum_id); scoped = match_by_pair[(precursor.candidate_id, precursor.ms2_spectrum_id)]
        theoretical = product_by_candidate[precursor.candidate_id]
        matched_mzs = {round(x.observed_product_mz, 8) for x in scoped}
        explained = sum(i for m, i in (spectrum.peaks if spectrum else ()) if round(m, 8) in matched_mzs)
        total = sum(i for _, i in (spectrum.peaks if spectrum else ()))
        top = sorted((spectrum.peaks if spectrum else ()), key=lambda x: (-x[1], x[0]))[:10]
        top_total = sum(i for _, i in top); top_explained = sum(i for m, i in top if round(m, 8) in matched_mzs)
        unique = sum(1 for x in scoped if len(next(p for p in theoretical if p.product_ion_id == x.product_ion_id).shared_candidate_ids) == 1)
        shared = len(scoped) - unique
        blocks = []
        if not theoretical:
            blocks += ["NO_THEORETICAL_PRODUCT_IONS", "MASS_ONLY_CANDIDATE_NO_STRUCTURE_RULES"]
        if theoretical and not scoped: blocks.append("NO_MATCHED_PRODUCT_IONS")
        if spectrum is None or not spectrum.peaks: blocks.append("NO_POSITIVE_INTENSITY_MS2_PEAKS")
        if scoped:
            status, confidence = "PRODUCT_ION_SUPPORT", "MEDIUM" if len(scoped) > 1 else "LOW"
        elif theoretical:
            status, confidence = "NO_PRODUCT_ION_SUPPORT", "LOW"
        else:
            status, confidence = "PRECURSOR_COMPATIBLE_ONLY", "LOW"
        output.append(SpectrumMS2Summary(
            candidate_id=precursor.candidate_id, ms2_spectrum_id=precursor.ms2_spectrum_id,
            collision_energy=precursor.collision_energy,
            precursor_compatibility_status=precursor.precursor_compatibility_status,
            theoretical_product_ion_count=len(theoretical), matched_product_ion_count=len(scoped),
            unique_matched_product_ion_count=unique, shared_matched_product_ion_count=shared,
            diagnostic_product_ion_count=len(scoped),
            base_related_ion_matched=any(x.product_ion_label in {"BASE_MOLECULAR_ION", "PROTONATED_BASE"} for x in scoped),
            ribose_loss_ion_matched=any(x.product_ion_class == "RIBOSE_LOSS" for x in scoped),
            modification_specific_ion_matched=any(next(p for p in theoretical if p.product_ion_id == x.product_ion_id).candidate_specific for x in scoped),
            explained_intensity_fraction=explained / total if total else 0.0,
            top_10_peak_explained_fraction=top_explained / top_total if top_total else 0.0,
            median_product_mass_error=median(abs(x.ppm_error) for x in scoped) if scoped else None,
            product_ion_ambiguity_count=sum(x.match_ambiguity_status != "CANDIDATE_SCOPED" for x in scoped),
            spectrum_evidence_status=status, spectrum_evidence_confidence=confidence,
            spectrum_block_reasons=_blocks(blocks),
        ))
    return tuple(sorted(output, key=lambda x: (x.candidate_id, x.ms2_spectrum_id)))


def _recurrence_and_energy(candidates, precursors, products, matches, spectra):
    precursor_by_candidate = defaultdict(list); match_by_candidate_product = defaultdict(list)
    product_by_candidate = defaultdict(list); spectrum_map = {x.ms2_spectrum_id: x for x in spectra}
    for x in precursors: precursor_by_candidate[x.candidate_id].append(x)
    for x in products: product_by_candidate[x.candidate_id].append(x)
    for x in matches: match_by_candidate_product[(x.candidate_id, x.product_ion_id)].append(x)
    recurrence = []
    for candidate in candidates:
        compatible = precursor_by_candidate[candidate.candidate_id]; total = len({x.ms2_spectrum_id for x in compatible})
        precursor_map = {x.ms2_spectrum_id: x for x in compatible}
        for product in product_by_candidate[candidate.candidate_id]:
            supported = match_by_candidate_product[(candidate.candidate_id, product.product_ion_id)]
            ids = {x.ms2_spectrum_id for x in supported}; fraction = len(ids) / total if total else 0.0
            energies = tuple(sorted({str(precursor_map[x].collision_energy) if precursor_map[x].collision_energy is not None else "NOT_RECORDED" for x in ids}))
            rts = [spectrum_map[x].scan_time for x in ids if x in spectrum_map and spectrum_map[x].scan_time is not None]
            if total < 2: status = "SINGLE_SPECTRUM_ONLY" if ids else "INSUFFICIENT_SPECTRA"
            elif fraction >= 0.5: status = "RECURRENT_HIGH"
            elif len(ids) >= 2: status = "RECURRENT_SUPPORTIVE"
            elif len(ids) == 1: status = "SINGLE_SPECTRUM_ONLY"
            else: status = "INCONSISTENT"
            recurrence.append(ProductIonRecurrence(
                candidate_id=candidate.candidate_id, product_ion_id=product.product_ion_id,
                supporting_ms2_spectrum_count=len(ids), compatible_ms2_spectrum_count=total,
                ms2_recurrence_fraction=fraction, collision_energy_set=energies,
                first_supporting_rt=min(rts) if rts else None, last_supporting_rt=max(rts) if rts else None,
                product_ion_recurrence_status=status,
            ))
    energy_rows = []
    for candidate in candidates:
        groups = defaultdict(list)
        for row in precursor_by_candidate[candidate.candidate_id]: groups[(row.collision_energy, row.collision_energy_unit)].append(row)
        for (energy, unit), rows in groups.items():
            ids = {x.ms2_spectrum_id for x in rows}; scoped_matches = [x for x in matches if x.candidate_id == candidate.candidate_id and x.ms2_spectrum_id in ids]
            unique_products = {x.product_ion_id for x in scoped_matches}
            energy_rows.append(CollisionEnergySummary(
                candidate_id=candidate.candidate_id, collision_energy_value=energy,
                collision_energy_unit=unit, collision_energy_status="RECORDED" if energy is not None else "NOT_RECORDED",
                spectrum_count_per_energy=len(ids), matched_product_ion_count_per_energy=len(unique_products),
                diagnostic_ion_recurrence_per_energy=(sum(1 for x in unique_products if sum(m.product_ion_id == x for m in scoped_matches) >= 2) / len(unique_products) if unique_products else None),
            ))
    return tuple(sorted(recurrence, key=lambda x: (x.candidate_id, x.product_ion_id))), tuple(sorted(energy_rows, key=lambda x: (x.candidate_id, x.collision_energy_value is None, x.collision_energy_value or 0.0)))


def summarize_candidate_ms2_identity_evidence(candidate_records, precursor_records, product_matches,
                                               *, product_ions=(), spectrum_summaries=(), recurrence_records=(), spectra=()):
    precursors = defaultdict(list); products = defaultdict(list); matches = defaultdict(list); summaries = defaultdict(list); recurrences = defaultdict(list)
    for x in precursor_records: precursors[x.candidate_id].append(x)
    for x in product_ions: products[x.candidate_id].append(x)
    for x in product_matches: matches[x.candidate_id].append(x)
    for x in spectrum_summaries: summaries[x.candidate_id].append(x)
    for x in recurrence_records: recurrences[x.candidate_id].append(x)
    usable_ids = {x.ms2_spectrum_id for x in spectra if x.ms2_preprocessing_status == "COMPLETED"}
    output = []
    for candidate in sorted(candidate_records, key=lambda x: x.candidate_id):
        cp = precursors[candidate.candidate_id]; ps = products[candidate.candidate_id]; ms = matches[candidate.candidate_id]; ss = summaries[candidate.candidate_id]; rr = recurrences[candidate.candidate_id]
        recurrent = [x for x in rr if x.supporting_ms2_spectrum_count >= 2]
        specific = [x for x in ps if x.candidate_specific]
        shared_recurrent = sum(len(next(p for p in ps if p.product_ion_id == x.product_ion_id).shared_candidate_ids) > 1 for x in recurrent)
        unique_recurrent = sum(len(next(p for p in ps if p.product_ion_id == x.product_ion_id).shared_candidate_ids) == 1 and next(p for p in ps if p.product_ion_id == x.product_ion_id).candidate_specific for x in recurrent)
        best = max(ss, key=lambda x: (x.matched_product_ion_count, x.explained_intensity_fraction, x.ms2_spectrum_id), default=None)
        blocks = []
        candidate_class = candidate.candidate_class.value if isinstance(candidate.candidate_class, Enum) else str(candidate.candidate_class)
        if not cp: blocks.append("NO_PRECURSOR_COMPATIBLE_MS2")
        if not ps:
            blocks += ["NO_THEORETICAL_PRODUCT_IONS"]
            if candidate_class == NucleosideCandidateClass.MASS_ONLY_MODIFIED_NUCLEOSIDE.value: blocks.append("MASS_ONLY_CANDIDATE_NO_STRUCTURE_RULES")
        if ps and not ms: blocks.append("NO_MATCHED_PRODUCT_IONS")
        if len(cp) > 1 and not recurrent: blocks.append("INSUFFICIENT_MS2_RECURRENCE")
        if len({x.collision_energy for x in cp}) > 1: blocks.append("COLLISION_ENERGY_HETEROGENEITY")
        if specific and any(x.candidate_specific and x.product_ion_id in {m.product_ion_id for m in ms} for x in ps):
            status, confidence, after = "MS2_SUPPORTS_CANDIDATE_CLASS", "MEDIUM", "PARTIALLY_REDUCED_NOT_RESOLVED"
        elif candidate_class == NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE.value and ms:
            status, confidence, after = "MS2_SUPPORTS_CANONICAL_NUCLEOSIDE_CLASS", "MEDIUM" if recurrent else "LOW", "CLASS_SUPPORTED_EXACT_IDENTITY_UNCONFIRMED"
        elif ps and ms and all(len(x.shared_candidate_ids) > 1 for x in ps if x.product_ion_id in {m.product_ion_id for m in ms}):
            status, confidence, after = "MS2_SUPPORTS_SHARED_ISOBARIC_CLASS", "LOW", "ISOBARIC_IDENTITY_UNRESOLVED"
        elif cp and not ps:
            status, confidence, after = "MS2_PRECURSOR_COMPATIBLE_ONLY", "LOW", "UNCHANGED"
        elif cp:
            status, confidence, after = "MS2_INSUFFICIENT", "LOW", "UNCHANGED"
        else:
            status, confidence, after = "MS2_BLOCKED", "NONE", "UNCHANGED"
        if status in {"MS2_PRECURSOR_COMPATIBLE_ONLY", "MS2_INSUFFICIENT"}: blocks.append("MS2_NONDISCRIMINATING")
        output.append(NucleosideCandidateMS2Summary(
            candidate_id=candidate.candidate_id, candidate_name=candidate.candidate_name,
            compatible_ms2_spectrum_count=len({x.ms2_spectrum_id for x in cp}),
            usable_ms2_spectrum_count=len({x.ms2_spectrum_id for x in cp} & usable_ids),
            collision_energy_count=len({x.collision_energy for x in cp if x.collision_energy is not None}),
            theoretical_product_ion_count=len(ps), recurrent_matched_product_ion_count=len(recurrent),
            candidate_unique_recurrent_ion_count=unique_recurrent, shared_recurrent_ion_count=shared_recurrent,
            best_spectrum_id=best.ms2_spectrum_id if best else "",
            best_spectrum_evidence_status=best.spectrum_evidence_status if best else "NOT_APPLICABLE",
            median_explained_intensity_fraction=median(x.explained_intensity_fraction for x in ss) if ss else None,
            median_top10_explained_fraction=median(x.top_10_peak_explained_fraction for x in ss) if ss else None,
            median_mass_error=median(abs(x.ppm_error) for x in ms) if ms else None,
            ms2_identity_evidence_status=status, ms2_identity_confidence=confidence,
            identity_ambiguity_before_ms2=candidate.ms1_identity_ambiguity_status,
            identity_ambiguity_after_ms2=after, candidate_specific_ms2_rules_available=bool(specific),
            ms2_block_reasons=_blocks(blocks + ["CHEMICAL_IDENTITY_UNSUPPORTED", "EXACT_ISOMER_DISCRIMINATION_UNSUPPORTED"]),
        ))
    return tuple(output)


def reconcile_p1ap_ms1_and_ms2_evidence(ms1_result, ms2_candidate_summaries, *, t1_result=None, full_length_series=None):
    if ms1_result is None: return ()
    by_id = {x.candidate_id: x for x in ms2_candidate_summaries}; rows = []
    for family in ms1_result.state_families:
        base = by_id.get(family.base_candidate_id)
        canonical_alternative = any(x.candidate_name == "G" and x.ms2_identity_evidence_status == "MS2_SUPPORTS_CANONICAL_NUCLEOSIDE_CLASS" for x in ms2_candidate_summaries)
        if base and canonical_alternative:
            status = "MS2_SUPPORTS_BASE_CANDIDATE_CLASS_AND_PLUS16_CANONICAL_G_ALTERNATIVE"
        elif base and base.ms2_identity_evidence_status == "MS2_SUPPORTS_CANONICAL_NUCLEOSIDE_CLASS":
            status = "MS2_SUPPORTS_BASE_CANDIDATE_CLASS"
        elif canonical_alternative:
            status = "MS2_SUPPORTS_PLUS16_CANONICAL_G_ALTERNATIVE"
        else: status = "MS2_INSUFFICIENT_FOR_STATE_RECONCILIATION"
        rows.append(P1APMS1MS2Reconciliation(
            reconciliation_id=f"MS1MS2__{family.state_family_id}", p1ap_ms1_state_family_id=family.state_family_id,
            p1ap_ms2_status=status, t1_reconciliation_status="T1_STATE_FAMILY_ABSENT" if t1_result is not None else "T1_RESULT_NOT_PROVIDED",
            full_length_reconciliation_status="P1AP_MS2_NONDISCRIMINATING" if full_length_series is not None else "FULL_LENGTH_RESULT_NOT_PROVIDED",
            state_interpretation_resolved=False,
            reconciliation_block_reasons=("MS2_NONDISCRIMINATING", "CHEMICAL_IDENTITY_UNSUPPORTED"),
        ))
    return tuple(sorted(rows, key=lambda x: x.reconciliation_id))


def audit_p1ap_nucleoside_ms2_identity(
    mzml_path: Path, *, p1ap_ms1_result: P1APNucleosideStateAuditResult | None,
    source_metadata_record: Any | None = None, runtime_context: Mapping[str, Any] | None = None,
    t1_result: Any | None = None, full_length_series: Any | None = None,
    selection_config: Mapping[str, Any] | None = None,
    preprocessing_config: Mapping[str, Any] | None = None,
    matching_config: Mapping[str, Any] | None = None,
    product_ion_registry: Mapping[str, Sequence[Mapping[str, Any]]] | None = None,
) -> P1APNucleosideMS2AuditResult:
    del runtime_context
    started = perf_counter(); path = Path(mzml_path)
    if p1ap_ms1_result is None:
        blocks = ["P1AP_MS1_AUDIT_RESULT_MISSING"]
        if not path.is_file(): blocks.append("INPUT_FILE_NOT_FOUND")
        if source_metadata_record is None: blocks.append("SOURCE_METADATA_RECORD_MISSING")
        summary = P1APNucleosideMS2Summary(
            source_id=path.stem, status="BLOCKED", total_spectra_seen=0, total_ms2_metadata_records=0,
            precursor_compatible_record_count=0, decoded_unique_spectrum_count=0, candidate_count=0,
            theoretical_product_ion_count=0, product_match_count=0,
            overall_evidence_status="P1AP_MS2_BLOCKED", overall_confidence="NONE",
            overall_block_reasons=_blocks(blocks), runtime_seconds=perf_counter() - started,
        )
        return P1APNucleosideMS2AuditResult((), (), (), (), (), (), (), (), (), (), (), summary)
    candidates = candidate_records_from_ms1_result(p1ap_ms1_result)
    base_blocks = []
    if not path.is_file(): base_blocks.append("INPUT_FILE_NOT_FOUND")
    if source_metadata_record is None: base_blocks.append("SOURCE_METADATA_RECORD_MISSING")
    if path.is_file():
        precursors = tuple(select_precursor_compatible_ms2_spectra(path, candidates, selection_config=selection_config))
        selected_ids = sorted({x.ms2_spectrum_id for x in precursors})
        spectra = tuple(decode_selected_ms2_spectra(path, selected_ids, preprocessing_config=preprocessing_config))
    else: precursors = (); spectra = ()
    products = tuple(generate_nucleoside_product_ion_hypotheses(candidates, product_ion_registry=product_ion_registry))
    discrimination = tuple(compare_candidate_product_ions(candidates, products))
    matches = tuple(match_product_ions_to_ms2_peaks(products, spectra, matching_config=matching_config, precursor_records=precursors))
    spectrum_summaries = _spectrum_summaries(candidates, precursors, spectra, products, matches)
    recurrence, collision = _recurrence_and_energy(candidates, precursors, products, matches, spectra)
    candidate_summaries = summarize_candidate_ms2_identity_evidence(
        candidates, precursors, matches, product_ions=products, spectrum_summaries=spectrum_summaries,
        recurrence_records=recurrence, spectra=spectra,
    )
    reconciliation = reconcile_p1ap_ms1_and_ms2_evidence(p1ap_ms1_result, candidate_summaries, t1_result=t1_result, full_length_series=full_length_series)
    if not precursors: base_blocks.append("NO_PRECURSOR_COMPATIBLE_MS2")
    if not matches: base_blocks.append("NO_MATCHED_PRODUCT_IONS")
    base_blocks += ["MS2_NONDISCRIMINATING", "CHEMICAL_IDENTITY_UNSUPPORTED", "EXACT_ISOMER_DISCRIMINATION_UNSUPPORTED"]
    metadata_count = int(getattr(getattr(p1ap_ms1_result, "run_summary", None), "ms2_spectra_present", 0) or 0)
    status = "COMPLETED" if path.is_file() and precursors else "BLOCKED"
    summary = P1APNucleosideMS2Summary(
        source_id=str(getattr(getattr(p1ap_ms1_result, "run_summary", None), "source_id", path.stem)),
        status=status, total_spectra_seen=metadata_count, total_ms2_metadata_records=metadata_count,
        precursor_compatible_record_count=len(precursors), decoded_unique_spectrum_count=len(spectra),
        candidate_count=len(candidates), theoretical_product_ion_count=len(products), product_match_count=len(matches),
        overall_evidence_status="P1AP_MS2_NUCLEOSIDE_CLASS_SUPPORT_NONDISCRIMINATING" if matches else "P1AP_MS2_PRECURSOR_ONLY_OR_INSUFFICIENT",
        overall_confidence="LOW", overall_block_reasons=_blocks(base_blocks), runtime_seconds=perf_counter() - started,
    )
    return P1APNucleosideMS2AuditResult(candidates, precursors, spectra, products, discrimination, matches,
        spectrum_summaries, collision, recurrence, candidate_summaries, reconciliation, summary)


def _record(value: Any) -> dict[str, Any]:
    row = asdict(value)
    row.pop("peaks", None)
    def normalize(item):
        if isinstance(item, Enum): return item.value
        if isinstance(item, dict): return {key: normalize(val) for key, val in item.items()}
        if isinstance(item, (tuple, list)): return [normalize(x) for x in item]
        return item
    return normalize(row)


def audit_optional_result(result: P1APNucleosideMS2AuditResult) -> dict[str, Any]:
    return {
        "precursor_records": [_record(x) for x in result.precursor_records],
        "spectrum_records": [_record(x) for x in result.spectrum_records],
        "theoretical_product_records": [_record(x) for x in result.theoretical_product_records],
        "discrimination_records": [_record(x) for x in result.discrimination_records],
        "product_match_records": [_record(x) for x in result.product_match_records],
        "spectrum_summary_records": [_record(x) for x in result.spectrum_summary_records],
        "collision_energy_records": [_record(x) for x in result.collision_energy_records],
        "recurrence_records": [_record(x) for x in result.recurrence_records],
        "candidate_summary_records": [_record(x) for x in result.candidate_summary_records],
        "reconciliation_records": [_record(x) for x in result.reconciliation_records],
        "summary_records": [_record(result.summary)],
    }
