"""Numerical, shadow-only linkage audit between a T1 txt profile and mzML MS1 profiles."""
from __future__ import annotations

from bisect import bisect_left
from dataclasses import asdict, dataclass, field, replace
from enum import Enum
from hashlib import sha256
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Mapping, Sequence
import csv
import re

import numpy as np
from pyteomics import mzml
from scipy.stats import spearmanr

from rna_masshunter.sciex_mzml_source_metadata_audit import MzMLSourceMetadataRecord, PolarityStatus
from rna_masshunter.sciex_t1_profile_peak_audit import (
    T1PeakDetectionParameters,
    T1PeakQualityClass,
    detect_t1_profile_peaks,
)
from rna_masshunter.sciex_t1_replicate_consistency_audit import (
    MatchAmbiguityStatus,
    ReplicateAuditParameters,
    ReplicateConsistencyAuditResult,
    ReplicateConsistencyStatus,
    ReplicateRunPeak,
    ReplicateRunPeakProfile,
    _rt_minutes,
    build_ms1_peak_profile_from_spectra,
    match_replicate_peaks,
)

OPTIONAL_RESULT_KEY = "sciex_t1_txt_mzml_source_linkage_audit"
ALGORITHM_VERSION = "sciex-t1-txt-mzml-source-linkage-audit-v1"

_BLOCK_ORDER = (
    "TXT_FILE_NOT_FOUND", "TXT_FILE_UNREADABLE", "TXT_PARSE_FAILED",
    "TXT_COLUMN_ASSIGNMENT_UNRESOLVED", "TXT_PROFILE_TYPE_UNRESOLVED", "TXT_EMPTY",
    "TXT_NO_NUMERIC_ROWS", "TXT_MZ_NOT_MONOTONIC", "TXT_DUPLICATE_MZ_VALUES",
    "TXT_NEGATIVE_INTENSITIES", "MZML_RUN_PROFILE_MISSING", "SOURCE_METADATA_RECORD_MISSING",
    "REPLICATE_AUDIT_RESULT_MISSING", "NO_OVERLAPPING_MZ_RANGE", "INSUFFICIENT_MZ_COVERAGE",
    "NO_TXT_PEAKS", "NO_REFERENCE_PEAKS", "NO_MATCHED_PEAKS", "AMBIGUOUS_PEAK_MATCH",
    "LOW_PROFILE_CORRELATION", "LOW_PEAK_OVERLAP", "INSUFFICIENT_DISCRIMINATING_PEAKS",
    "CONFLICTING_DISCRIMINATING_PEAKS", "MULTIPLE_HYPOTHESES_SIMILAR",
    "MEAN_SUM_HYPOTHESES_NOT_IDENTIFIABLE_AFTER_NORMALIZATION",
    "AGGREGATE_HYPOTHESIS_NON_UNIQUE", "PARTIAL_SCAN_HYPOTHESIS_NOT_TESTED",
    "INSUFFICIENT_SOURCE_LINKAGE", "SOURCE_POLARITY_CONFLICT", "USER_MANIFEST_METADATA_CONFLICT",
)
_LOW_QUALITY = {T1PeakQualityClass.LOW_SUPPORT.value, T1PeakQualityClass.SHOULDER_OR_OVERLAP.value}


def _ordered_blocks(values: Iterable[str]) -> tuple[str, ...]:
    found = set(values)
    return tuple(x for x in _BLOCK_ORDER if x in found) + tuple(sorted(found - set(_BLOCK_ORDER)))


def _safe_float(value: Any) -> float | None:
    try:
        result = float(value)
        return result if np.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def _median(values: Iterable[float | None]) -> float | None:
    usable = [float(value) for value in values if value is not None and np.isfinite(value)]
    return float(median(usable)) if usable else None


def _correlation(left: np.ndarray, right: np.ndarray) -> float | None:
    if len(left) < 3 or np.std(left) == 0 or np.std(right) == 0:
        return None
    value = float(np.corrcoef(left, right)[0, 1])
    return value if np.isfinite(value) else None


def _spearman(left: np.ndarray, right: np.ndarray) -> float | None:
    if len(left) < 3 or len(np.unique(left)) < 2 or len(np.unique(right)) < 2:
        return None
    value = float(spearmanr(left, right).statistic)
    return value if np.isfinite(value) else None


def _cosine(left: np.ndarray, right: np.ndarray) -> float | None:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.dot(left, right) / denominator) if denominator > 0 else None


def _normalize(values: np.ndarray, method: str) -> np.ndarray | None:
    positive = np.maximum(np.asarray(values, dtype=float), 0)
    denominator = float(positive.max()) if method == "BASE_PEAK" else float(positive.sum())
    return positive / denominator if denominator > 0 else None


class TxtProfileType(str, Enum):
    DENSE_PROFILE = "DENSE_PROFILE"
    SPARSE_PEAK_LIST = "SPARSE_PEAK_LIST"
    IRREGULAR_PROFILE = "IRREGULAR_PROFILE"
    UNRESOLVED = "UNRESOLVED"


class HypothesisType(str, Enum):
    RUN_1_ONLY_EXPORT = "RUN_1_ONLY_EXPORT"
    RUN_2_ONLY_EXPORT = "RUN_2_ONLY_EXPORT"
    MEAN_OF_RUNS_EXPORT = "MEAN_OF_RUNS_EXPORT"
    SUM_OF_RUNS_EXPORT = "SUM_OF_RUNS_EXPORT"
    MEDIAN_OF_RUNS_EXPORT = "MEDIAN_OF_RUNS_EXPORT"
    PARTIAL_SCAN_EXPORT_POSSIBLE = "PARTIAL_SCAN_EXPORT_POSSIBLE"
    DIFFERENT_PROCESSING_EXPORT_POSSIBLE = "DIFFERENT_PROCESSING_EXPORT_POSSIBLE"
    NO_SUPPORTED_LINKAGE = "NO_SUPPORTED_LINKAGE"
    INSUFFICIENT_EVIDENCE = "INSUFFICIENT_EVIDENCE"


class LinkageStatus(str, Enum):
    STRONG_LINK_TO_RUN_1 = "STRONG_LINK_TO_RUN_1"
    STRONG_LINK_TO_RUN_2 = "STRONG_LINK_TO_RUN_2"
    SUPPORTIVE_LINK_TO_RUN_1 = "SUPPORTIVE_LINK_TO_RUN_1"
    SUPPORTIVE_LINK_TO_RUN_2 = "SUPPORTIVE_LINK_TO_RUN_2"
    BEST_MATCH_AGGREGATE_PROFILE = "BEST_MATCH_AGGREGATE_PROFILE"
    POSSIBLE_PARTIAL_SCAN_EXPORT_RUN_1 = "POSSIBLE_PARTIAL_SCAN_EXPORT_RUN_1"
    POSSIBLE_PARTIAL_SCAN_EXPORT_RUN_2 = "POSSIBLE_PARTIAL_SCAN_EXPORT_RUN_2"
    MULTIPLE_HYPOTHESES_SIMILAR = "MULTIPLE_HYPOTHESES_SIMILAR"
    DIFFERENT_PROCESSING_EXPORT_POSSIBLE = "DIFFERENT_PROCESSING_EXPORT_POSSIBLE"
    NO_SUPPORTED_LINKAGE = "NO_SUPPORTED_LINKAGE"
    INSUFFICIENT_EVIDENCE = "INSUFFICIENT_EVIDENCE"
    BLOCKED_BY_INPUT_FORMAT = "BLOCKED_BY_INPUT_FORMAT"


@dataclass(frozen=True)
class SourceLinkageParameters:
    replicate_parameters: ReplicateAuditParameters = field(default_factory=ReplicateAuditParameters)
    strong_score: float = 0.75
    supportive_score: float = 0.55
    minimum_polarity_support_score: float = 0.55
    no_support_score: float = 0.30
    similar_score_margin: float = 0.05
    unique_score_margin: float = 0.10
    minimum_coverage_fraction: float = 0.25
    minimum_profile_correlation: float = 0.30
    minimum_peak_overlap: float = 0.20
    minimum_discriminating_support: int = 2
    weight_profile_correlation: float = 0.25
    weight_spearman: float = 0.10
    weight_cosine: float = 0.10
    weight_txt_peak_overlap: float = 0.25
    weight_top25_match: float = 0.15
    weight_strict_match: float = 0.15
    partial_window_count: int = 3
    minimum_ms1_scans_per_window: int = 3
    maximum_windows: int = 6

    def validate(self) -> None:
        scores = (self.strong_score, self.supportive_score, self.minimum_polarity_support_score,
                  self.no_support_score, self.similar_score_margin, self.unique_score_margin,
                  self.minimum_coverage_fraction, self.minimum_profile_correlation,
                  self.minimum_peak_overlap)
        if any(not 0 <= value <= 1 for value in scores):
            raise ValueError("source-linkage score parameters must be in [0,1]")
        if not self.no_support_score <= self.supportive_score <= self.strong_score:
            raise ValueError("invalid linkage score ordering")
        weights = (self.weight_profile_correlation, self.weight_spearman, self.weight_cosine,
                   self.weight_txt_peak_overlap, self.weight_top25_match, self.weight_strict_match)
        if any(value < 0 for value in weights) or abs(sum(weights) - 1.0) > 1e-12:
            raise ValueError("composite linkage weights must be nonnegative and sum to one")
        if not 1 <= self.partial_window_count <= self.maximum_windows:
            raise ValueError("invalid partial window count")


@dataclass(frozen=True)
class T1TxtProfile:
    input_path: str
    file_name: str
    file_size_bytes: int | str
    sha256: str
    encoding_status: str
    delimiter: str
    header_status: str
    column_count: int
    row_count: int
    numeric_row_count: int
    non_numeric_row_count: int
    mz_column_index: int | str
    intensity_column_index: int | str
    mz_min: float | str
    mz_max: float | str
    intensity_min: float | str
    intensity_max: float | str
    zero_intensity_count: int
    negative_intensity_count: int
    duplicate_mz_count: int
    mz_sorted_status: str
    mz_spacing_median: float | str
    mz_spacing_mad: float | str
    local_maximum_count: int
    profile_or_peaklist_status: TxtProfileType
    parse_status: str
    parse_block_reasons: tuple[str, ...]
    coordinates: np.ndarray = field(repr=False, compare=False)
    intensities: np.ndarray = field(repr=False, compare=False)
    shadow_analysis_only: bool = True
    source_linkage_audit_only: bool = True
    formal_propagation: bool = False
    polarity_propagation_applied: bool = False
    chemical_identity_assigned: bool = False
    fragment_identity_assigned: bool = False
    charge_state_confirmed: bool = False


@dataclass(frozen=True)
class T1TxtPeakProfile:
    txt_profile: T1TxtProfile
    detected_peak_count: int
    comparison_peak_count: int
    comparison_profile: ReplicateRunPeakProfile
    block_reasons: tuple[str, ...]
    formal_propagation: bool = False


@dataclass(frozen=True)
class SourceLinkageReferenceProfile:
    hypothesis_id: str
    hypothesis_type: HypothesisType
    reference_run_labels: tuple[str, ...]
    run_profile: ReplicateRunPeakProfile
    window_start_time: float | None = None
    window_end_time: float | None = None
    window_ms1_scan_count: int | None = None
    block_reasons: tuple[str, ...] = ()
    formal_propagation: bool = False


@dataclass(frozen=True, kw_only=True)
class TxtMzMLLinkageHypothesisResult:
    hypothesis_id: str
    hypothesis_type: HypothesisType
    reference_run_labels: tuple[str, ...]
    raw_profile_correlation: float | None
    base_peak_normalized_correlation: float | None
    sum_normalized_correlation: float | None
    log1p_correlation: float | None
    spearman_rank_correlation: float | None
    cosine_similarity: float | None
    comparison_grid_method: str
    comparison_grid_spacing: float | None
    overlap_mz_min: float | None
    overlap_mz_max: float | None
    txt_coverage_fraction: float
    mzml_coverage_fraction: float
    interpolation_method: str
    extrapolation_applied: bool
    txt_peak_count: int
    reference_peak_count: int
    matched_peak_count: int
    txt_only_peak_count: int
    reference_only_peak_count: int
    ambiguous_match_count: int
    strict_match_count: int
    supportive_match_count: int
    txt_overlap_fraction: float
    reference_overlap_fraction: float
    peak_jaccard: float
    median_absolute_delta_mz: float | None
    median_ppm_error: float | None
    mad_delta_mz: float | None
    top_10_txt_peak_match_fraction: float
    top_25_txt_peak_match_fraction: float
    discriminating_peak_support_count: int
    discriminating_peak_conflict_count: int
    low_quality_match_count: int
    coverage_fraction: float
    composite_linkage_score: float
    linkage_status: LinkageStatus
    linkage_confidence: str
    block_reasons: tuple[str, ...]
    window_start_time: float | None = None
    window_end_time: float | None = None
    window_ms1_scan_count: int | None = None
    shadow_analysis_only: bool = True
    source_linkage_audit_only: bool = True
    formal_propagation: bool = False
    polarity_propagation_applied: bool = False
    chemical_identity_assigned: bool = False
    fragment_identity_assigned: bool = False
    charge_state_confirmed: bool = False


@dataclass(frozen=True)
class TxtPeakLinkageEvidence:
    txt_peak_id: str
    txt_apex_mz: float
    txt_intensity: float
    txt_normalized_intensity: float
    txt_intensity_rank: int
    run_1_matched: bool
    run_1_peak_id: str | None
    run_1_delta_mz: float | None
    run_1_ppm_error: float | None
    run_1_intensity_rank: int | None
    run_1_scan_recurrence: float | None
    run_1_replicate_status: str
    run_2_matched: bool
    run_2_peak_id: str | None
    run_2_delta_mz: float | None
    run_2_ppm_error: float | None
    run_2_intensity_rank: int | None
    run_2_scan_recurrence: float | None
    run_2_replicate_status: str
    aggregate_hypothesis_matches: tuple[str, ...]
    peak_linkage_status: str
    peak_linkage_confidence: str
    block_reasons: tuple[str, ...]
    shadow_analysis_only: bool = True
    source_linkage_audit_only: bool = True
    formal_propagation: bool = False
    polarity_propagation_applied: bool = False
    chemical_identity_assigned: bool = False
    fragment_identity_assigned: bool = False
    charge_state_confirmed: bool = False


@dataclass(frozen=True)
class DiscriminatingPeakEvidence:
    observed_mz: float
    discrimination_type: str
    run_1_status: str
    run_2_status: str
    run_1_normalized_intensity: float | None
    run_2_normalized_intensity: float | None
    run_1_rank: int | None
    run_2_rank: int | None
    run_1_recurrence: float | None
    run_2_recurrence: float | None
    txt_matched: bool
    txt_intensity: float | None
    txt_rank: int | None
    supports_hypothesis: str
    discrimination_confidence: str
    shadow_analysis_only: bool = True
    source_linkage_audit_only: bool = True
    formal_propagation: bool = False
    polarity_propagation_applied: bool = False
    chemical_identity_assigned: bool = False
    fragment_identity_assigned: bool = False
    charge_state_confirmed: bool = False


@dataclass(frozen=True, kw_only=True)
class TxtMzMLSourceLinkageSummary:
    txt_file: str
    rna_identity: str
    digest_type: str
    context_source: str
    candidate_run_count: int
    best_hypothesis: str
    second_best_hypothesis: str
    best_linkage_status: LinkageStatus
    best_linkage_confidence: str
    best_composite_score: float
    second_best_composite_score: float
    score_margin: float
    profile_evidence: str
    peak_evidence: str
    discriminating_peak_evidence: str
    partial_scan_evidence: str
    aggregate_evidence: str
    source_linkage_confirmed: bool
    exact_run_linkage_confirmed: bool
    source_polarity: str
    source_polarity_evidence: str
    common_source_polarity_supported: bool
    polarity_propagation_eligible: bool
    polarity_propagation_applied: bool
    polarity_propagation_block_reasons: tuple[str, ...]
    overall_block_reasons: tuple[str, ...]
    shadow_analysis_only: bool = True
    source_linkage_audit_only: bool = True
    formal_propagation: bool = False
    chemical_identity_assigned: bool = False
    fragment_identity_assigned: bool = False
    charge_state_confirmed: bool = False


@dataclass(frozen=True)
class TxtMzMLSourceLinkageAuditResult:
    parameters: SourceLinkageParameters
    txt_profile: T1TxtProfile
    txt_peak_profile: T1TxtPeakProfile
    reference_profiles: tuple[SourceLinkageReferenceProfile, ...]
    hypothesis_results: tuple[TxtMzMLLinkageHypothesisResult, ...]
    peak_evidence: tuple[TxtPeakLinkageEvidence, ...]
    discriminating_evidence: tuple[DiscriminatingPeakEvidence, ...]
    summary: TxtMzMLSourceLinkageSummary
    algorithm_version: str = ALGORITHM_VERSION
    formal_propagation: bool = False


def _empty_txt(path: Path, status: str, blocks: Sequence[str]) -> T1TxtProfile:
    return T1TxtProfile(str(path), path.name, "NOT_AVAILABLE", "NOT_AVAILABLE", "NOT_AVAILABLE", "NOT_AVAILABLE",
        "NOT_AVAILABLE", 0, 0, 0, 0, "NOT_AVAILABLE", "NOT_AVAILABLE", "NOT_AVAILABLE", "NOT_AVAILABLE",
        "NOT_AVAILABLE", "NOT_AVAILABLE", 0, 0, 0, "UNRESOLVED", "NOT_AVAILABLE", "NOT_AVAILABLE", 0,
        TxtProfileType.UNRESOLVED, status, _ordered_blocks(blocks), np.array([], dtype=float), np.array([], dtype=float))


def _split_line(line: str, delimiter: str) -> list[str]:
    if delimiter == "COMMA":
        return [item.strip() for item in next(csv.reader([line]))]
    if delimiter == "TAB":
        return [item.strip() for item in line.split("\t")]
    return [item for item in re.split(r"\s+", line.strip()) if item]


def _header_indices(fields: Sequence[str]) -> tuple[int | None, int | None]:
    normalized = [re.sub(r"[^a-z0-9]+", "", value.casefold()) for value in fields]
    mz_names = {"mz", "masscharge", "masstocharge", "massovercharge"}
    intensity_names = {"intensity", "abundance", "signal"}
    mz = next((index for index, value in enumerate(normalized) if value in mz_names), None)
    intensity = next((index for index, value in enumerate(normalized) if value in intensity_names), None)
    return mz, intensity


def _classify_profile(coordinates: np.ndarray, intensities: np.ndarray, duplicates: int, sorted_ok: bool) -> TxtProfileType:
    if len(coordinates) < 2 or not sorted_ok:
        return TxtProfileType.UNRESOLVED
    diff = np.diff(coordinates)
    positive = diff[diff > 0]
    if not len(positive):
        return TxtProfileType.UNRESOLVED
    spacing = float(np.median(positive))
    gap_fraction = float(np.mean(positive > spacing * 5)) if spacing > 0 else 1.0
    if len(coordinates) >= 1000 and duplicates == 0 and gap_fraction <= 0.01:
        return TxtProfileType.DENSE_PROFILE
    if len(coordinates) <= 1000 and (spacing >= 0.02 or gap_fraction > 0.05):
        return TxtProfileType.SPARSE_PEAK_LIST
    return TxtProfileType.IRREGULAR_PROFILE


def parse_t1_txt_profile(path: Path) -> T1TxtProfile:
    source = Path(path)
    if not source.exists():
        return _empty_txt(source, "BLOCKED", ("TXT_FILE_NOT_FOUND",))
    if not source.is_file():
        return _empty_txt(source, "BLOCKED", ("TXT_FILE_UNREADABLE",))
    try:
        data = source.read_bytes()
    except OSError:
        return _empty_txt(source, "BLOCKED", ("TXT_FILE_UNREADABLE",))
    digest = sha256(data).hexdigest()
    try:
        text = data.decode("utf-8-sig")
        encoding = "UTF_8"
    except UnicodeDecodeError:
        try:
            text = data.decode("utf-8-sig", errors="replace")
            encoding = "UTF_8_WITH_REPLACEMENT"
        except Exception:
            return _empty_txt(source, "BLOCKED", ("TXT_PARSE_FAILED",))
    logical = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith(("#", ";"))]
    if not logical:
        return replace(_empty_txt(source, "BLOCKED", ("TXT_EMPTY",)), file_size_bytes=len(data), sha256=digest, encoding_status=encoding)
    first = logical[0]
    delimiter = "COMMA" if "," in first else "TAB" if "\t" in first else "WHITESPACE"
    first_fields = _split_line(first, delimiter)
    first_numeric = len(first_fields) >= 2 and all(_safe_float(value) is not None for value in first_fields[:2])
    header = not first_numeric
    mz_index, intensity_index = _header_indices(first_fields) if header else (0, 1) if len(first_fields) == 2 else (None, None)
    blocks: list[str] = []
    if mz_index is None or intensity_index is None or mz_index == intensity_index:
        blocks.append("TXT_COLUMN_ASSIGNMENT_UNRESOLVED")
    values: list[tuple[float, float]] = []
    nonnumeric = 0
    rows = logical[1:] if header else logical
    expected_columns = len(first_fields)
    if not blocks:
        for line in rows:
            fields = _split_line(line, delimiter)
            if len(fields) != expected_columns or max(mz_index, intensity_index) >= len(fields):
                nonnumeric += 1
                continue
            mz_value = _safe_float(fields[mz_index])
            intensity = _safe_float(fields[intensity_index])
            if mz_value is None or intensity is None:
                nonnumeric += 1
                continue
            values.append((mz_value, intensity))
    if not values:
        blocks.append("TXT_NO_NUMERIC_ROWS")
    coordinates = np.asarray([value[0] for value in values], dtype=float)
    intensities = np.asarray([value[1] for value in values], dtype=float)
    duplicates = len(coordinates) - len(np.unique(coordinates)) if len(coordinates) else 0
    sorted_ok = bool(len(coordinates) < 2 or np.all(np.diff(coordinates) > 0))
    negative = int(np.sum(intensities < 0)) if len(intensities) else 0
    if not sorted_ok:
        blocks.append("TXT_MZ_NOT_MONOTONIC")
    if duplicates:
        blocks.append("TXT_DUPLICATE_MZ_VALUES")
    if negative:
        blocks.append("TXT_NEGATIVE_INTENSITIES")
    diff = np.diff(coordinates) if len(coordinates) > 1 else np.array([], dtype=float)
    positive = diff[diff > 0]
    spacing = float(np.median(positive)) if len(positive) else "NOT_AVAILABLE"
    spacing_mad = float(np.median(np.abs(positive - np.median(positive)))) if len(positive) else "NOT_AVAILABLE"
    local_maxima = int(np.sum((intensities[1:-1] > intensities[:-2]) & (intensities[1:-1] > intensities[2:]))) if len(intensities) >= 3 else 0
    profile_type = _classify_profile(coordinates, intensities, duplicates, sorted_ok)
    if profile_type is TxtProfileType.UNRESOLVED:
        blocks.append("TXT_PROFILE_TYPE_UNRESOLVED")
    status = "COMPLETED" if values and not set(blocks) & {"TXT_COLUMN_ASSIGNMENT_UNRESOLVED", "TXT_NO_NUMERIC_ROWS", "TXT_MZ_NOT_MONOTONIC"} else "BLOCKED"
    return T1TxtProfile(
        str(source), source.name, len(data), digest, encoding, delimiter,
        "PRESENT" if header else "ABSENT_TWO_COLUMN_ASSUMED", expected_columns,
        len(rows), len(values), nonnumeric, mz_index if mz_index is not None else "NOT_AVAILABLE",
        intensity_index if intensity_index is not None else "NOT_AVAILABLE",
        float(coordinates.min()) if len(coordinates) else "NOT_AVAILABLE",
        float(coordinates.max()) if len(coordinates) else "NOT_AVAILABLE",
        float(intensities.min()) if len(intensities) else "NOT_AVAILABLE",
        float(intensities.max()) if len(intensities) else "NOT_AVAILABLE",
        int(np.sum(intensities == 0)), negative, duplicates,
        "STRICTLY_INCREASING" if sorted_ok else "NOT_MONOTONIC", spacing, spacing_mad,
        local_maxima, profile_type, status, _ordered_blocks(blocks), coordinates, intensities,
    )


def _peaks_to_replicate_profile(label: str, detected: Sequence[Any], *, comparison_grid=None,
                                raw_profile=None, normalized_profile=None) -> ReplicateRunPeakProfile:
    interim: list[ReplicateRunPeak] = []
    for serial, peak in enumerate(sorted(detected, key=lambda item: (item.apex_mz, item.t1_peak_id)), 1):
        interim.append(ReplicateRunPeak(
            run_label=label, peak_id=f"TXTLINKPEAK__{label}__{serial:05d}", apex_mz=peak.apex_mz,
            centroid_mz=peak.centroid_mz, raw_apex_intensity=peak.apex_intensity,
            normalized_apex_intensity=peak.relative_apex_intensity,
            raw_integrated_intensity=peak.integrated_intensity,
            normalized_integrated_intensity=peak.relative_integrated_intensity,
            relative_intensity=peak.relative_apex_intensity, intensity_rank=0,
            prominence=peak.prominence, relative_prominence=peak.relative_prominence,
            fwhm=peak.fwhm_mz, left_bound_mz=peak.left_boundary_mz, right_bound_mz=peak.right_boundary_mz,
            supporting_ms1_scan_count=0, total_ms1_scan_count=0, scan_recurrence_fraction=0,
            first_supporting_scan_time=None, last_supporting_scan_time=None,
            detection_status=peak.peak_quality_class.value, detection_block_reasons=(),
        ))
    ranked = sorted(interim, key=lambda item: (-item.normalized_apex_intensity, item.apex_mz, item.peak_id))
    ranks = {peak.peak_id: index for index, peak in enumerate(ranked, 1)}
    peaks = tuple(replace(peak, intensity_rank=ranks[peak.peak_id]) for peak in interim)
    return ReplicateRunPeakProfile(
        run_label=label, input_path="SHADOW_DERIVED", status="COMPLETED", aggregation_method="PROVIDED_PROFILE",
        ms1_spectra_total=0, ms1_spectra_used=0, ms1_spectra_excluded=0, ms2_spectra_excluded=0,
        missing_ms_level_spectra=0, mz_grid_method="PROVIDED_GRID", intensity_normalization_method="BASE_PEAK",
        baseline_method="INHERITED", smoothing_method="INHERITED", peak_detection_method="detect_t1_profile_peaks",
        detected_peak_count=len(peaks), peaks=peaks, polarity_status="UNKNOWN", representation_status="PROFILE",
        block_reasons=(), comparison_mz_grid=np.asarray(comparison_grid, dtype=float) if comparison_grid is not None else None,
        comparison_raw_profile=np.asarray(raw_profile, dtype=float) if raw_profile is not None else None,
        comparison_normalized_profile=np.asarray(normalized_profile, dtype=float) if normalized_profile is not None else None,
    )


def build_txt_profile_peaks(profile: T1TxtProfile, *, detection_config: Mapping[str, Any] | None = None) -> T1TxtPeakProfile:
    config = dict(detection_config or {})
    detection_parameters = config.pop("peak_detection_parameters", None)
    if config:
        raise ValueError(f"unsupported detection_config keys: {sorted(config)}")
    blocks = list(profile.parse_block_reasons)
    detected_count = 0
    comparison: ReplicateRunPeakProfile
    if profile.parse_status != "COMPLETED" or not len(profile.coordinates):
        blocks.append("NO_TXT_PEAKS")
        comparison = _peaks_to_replicate_profile("000_TXT_PROFILE", ())
    elif profile.profile_or_peaklist_status in {TxtProfileType.DENSE_PROFILE, TxtProfileType.IRREGULAR_PROFILE}:
        result = detect_t1_profile_peaks(
            profile.coordinates, np.maximum(profile.intensities, 0), source_id="TXT_PROFILE",
            measurement_id="TXT_PROFILE", rna_identity="USER_MANIFEST_TRNA",
            parameters=detection_parameters,
        )
        detected_count = len(result.peaks)
        comparison = _peaks_to_replicate_profile(
            "000_TXT_PROFILE", result.selected_peaks, comparison_grid=profile.coordinates,
            raw_profile=profile.intensities,
            normalized_profile=_normalize(profile.intensities, "BASE_PEAK"),
        )
    elif profile.profile_or_peaklist_status is TxtProfileType.SPARSE_PEAK_LIST:
        base = float(np.max(np.maximum(profile.intensities, 0))) if len(profile.intensities) else 0
        order = sorted(range(len(profile.coordinates)), key=lambda index: (-profile.intensities[index], profile.coordinates[index]))
        rank = {index: serial for serial, index in enumerate(order, 1)}
        peaks = tuple(ReplicateRunPeak(
            run_label="000_TXT_PROFILE", peak_id=f"TXTLINKPEAK__000_TXT_PROFILE__{serial:05d}",
            apex_mz=float(profile.coordinates[index]), centroid_mz=float(profile.coordinates[index]),
            raw_apex_intensity=float(profile.intensities[index]), normalized_apex_intensity=float(profile.intensities[index] / base) if base else 0,
            raw_integrated_intensity=float(profile.intensities[index]), normalized_integrated_intensity=float(profile.intensities[index] / np.sum(profile.intensities)) if np.sum(profile.intensities) else 0,
            relative_intensity=float(profile.intensities[index] / base) if base else 0, intensity_rank=rank[index],
            prominence=None, relative_prominence=None, fwhm=None, left_bound_mz=float(profile.coordinates[index]),
            right_bound_mz=float(profile.coordinates[index]), supporting_ms1_scan_count=0, total_ms1_scan_count=0,
            scan_recurrence_fraction=0, first_supporting_scan_time=None, last_supporting_scan_time=None,
            detection_status="SPARSE_INPUT_POINT", detection_block_reasons=("MISSING_PEAK_SHAPE_METRICS",),
        ) for serial, index in enumerate(np.argsort(profile.coordinates), 1) if profile.intensities[index] > 0)
        detected_count = len(peaks)
        comparison = replace(_peaks_to_replicate_profile("000_TXT_PROFILE", ()), detected_peak_count=len(peaks), peaks=peaks)
    else:
        blocks.extend(("TXT_PROFILE_TYPE_UNRESOLVED", "NO_TXT_PEAKS"))
        comparison = _peaks_to_replicate_profile("000_TXT_PROFILE", ())
    return T1TxtPeakProfile(profile, detected_count, len(comparison.peaks), comparison, _ordered_blocks(blocks))


def _reference_from_run(profile: ReplicateRunPeakProfile, index: int) -> SourceLinkageReferenceProfile:
    kind = HypothesisType.RUN_1_ONLY_EXPORT if index == 0 else HypothesisType.RUN_2_ONLY_EXPORT
    return SourceLinkageReferenceProfile(kind.value, kind, (profile.run_label,), profile)


def _derived_profile(label: str, grid: np.ndarray, raw: np.ndarray, normalized: np.ndarray) -> ReplicateRunPeakProfile:
    result = detect_t1_profile_peaks(grid, normalized, source_id=label, measurement_id=label,
                                     rna_identity="USER_MANIFEST_TRNA")
    return _peaks_to_replicate_profile(label, result.peaks, comparison_grid=grid,
                                       raw_profile=raw, normalized_profile=normalized)


def build_replicate_aggregate_profiles(run_profiles: Sequence[ReplicateRunPeakProfile]) -> tuple[SourceLinkageReferenceProfile, ...]:
    ordered = tuple(sorted(run_profiles, key=lambda profile: (profile.run_label, profile.input_path)))
    if len(ordered) < 2 or any(profile.comparison_mz_grid is None or profile.comparison_normalized_profile is None for profile in ordered):
        return ()
    low = max(float(profile.comparison_mz_grid[0]) for profile in ordered)
    high = min(float(profile.comparison_mz_grid[-1]) for profile in ordered)
    if high <= low:
        return ()
    steps = [float(np.median(np.diff(profile.comparison_mz_grid))) for profile in ordered]
    step = max(steps)
    grid = np.arange(low, high + step / 2, step)
    normalized_rows = np.vstack([np.interp(grid, profile.comparison_mz_grid, profile.comparison_normalized_profile) for profile in ordered])
    raw_rows = np.vstack([np.interp(grid, profile.comparison_mz_grid, profile.comparison_raw_profile if profile.comparison_raw_profile is not None else profile.comparison_normalized_profile) for profile in ordered])
    methods = (
        (HypothesisType.MEAN_OF_RUNS_EXPORT, np.mean(raw_rows, axis=0), np.mean(normalized_rows, axis=0)),
        (HypothesisType.SUM_OF_RUNS_EXPORT, np.sum(raw_rows, axis=0), np.sum(normalized_rows, axis=0)),
        (HypothesisType.MEDIAN_OF_RUNS_EXPORT, np.median(raw_rows, axis=0), np.median(normalized_rows, axis=0)),
    )
    output = []
    labels = tuple(profile.run_label for profile in ordered)
    for kind, raw_values, values in methods:
        normalized = _normalize(values, "BASE_PEAK")
        profile = _derived_profile(f"AGGREGATE_{kind.value}", grid, raw_values, normalized if normalized is not None else values)
        output.append(SourceLinkageReferenceProfile(kind.value, kind, labels, profile,
            block_reasons=("MEAN_SUM_HYPOTHESES_NOT_IDENTIFIABLE_AFTER_NORMALIZATION", "AGGREGATE_HYPOTHESIS_NON_UNIQUE")))
    return tuple(output)


def _match_orientation(match, txt_label: str):
    txt_left = match.run_1_label == txt_label
    txt_id = match.run_1_peak_id if txt_left else match.run_2_peak_id
    ref_id = match.run_2_peak_id if txt_left else match.run_1_peak_id
    delta = match.delta_mz if txt_left else -match.delta_mz if match.delta_mz is not None else None
    ppm = delta / (match.run_1_apex_mz if txt_left else match.run_2_apex_mz) * 1e6 if delta is not None else None
    return txt_id, ref_id, delta, ppm


def compare_txt_to_reference_profile(
    txt_profile: T1TxtPeakProfile, reference_profile: SourceLinkageReferenceProfile, *,
    matching_config: Mapping[str, Any] | None = None,
) -> TxtMzMLLinkageHypothesisResult:
    config = dict(matching_config or {})
    parameters = config.pop("parameters", None) or SourceLinkageParameters()
    parameters.validate()
    if config:
        raise ValueError(f"unsupported matching_config keys: {sorted(config)}")
    txt = txt_profile.comparison_profile
    ref = reference_profile.run_profile
    blocks = list(reference_profile.block_reasons)
    raw_corr = base_corr = sum_corr = log_corr = rank_corr = cosine = None
    overlap_min = overlap_max = spacing = None
    txt_coverage = ref_coverage = 0.0
    if txt.comparison_mz_grid is not None and ref.comparison_mz_grid is not None:
        overlap_min = max(float(txt.comparison_mz_grid[0]), float(ref.comparison_mz_grid[0]))
        overlap_max = min(float(txt.comparison_mz_grid[-1]), float(ref.comparison_mz_grid[-1]))
        if overlap_max > overlap_min:
            mask = (ref.comparison_mz_grid >= overlap_min) & (ref.comparison_mz_grid <= overlap_max)
            grid = ref.comparison_mz_grid[mask]
            spacing = float(np.median(np.diff(grid))) if len(grid) > 1 else None
            txt_raw = np.interp(grid, txt.comparison_mz_grid, txt.comparison_raw_profile)
            ref_raw_source = ref.comparison_raw_profile if ref.comparison_raw_profile is not None else ref.comparison_normalized_profile
            ref_raw = np.asarray(ref_raw_source[mask], dtype=float)
            txt_bp = _normalize(txt_raw, "BASE_PEAK")
            ref_bp = _normalize(ref_raw, "BASE_PEAK")
            txt_sum = _normalize(txt_raw, "SUM")
            ref_sum = _normalize(ref_raw, "SUM")
            raw_corr = _correlation(txt_raw, ref_raw)
            if txt_bp is not None and ref_bp is not None:
                base_corr = _correlation(txt_bp, ref_bp)
                log_corr = _correlation(np.log1p(txt_bp), np.log1p(ref_bp))
                rank_corr = _spearman(txt_bp, ref_bp)
                cosine = _cosine(txt_bp, ref_bp)
            if txt_sum is not None and ref_sum is not None:
                sum_corr = _correlation(txt_sum, ref_sum)
            txt_span = float(txt.comparison_mz_grid[-1] - txt.comparison_mz_grid[0])
            ref_span = float(ref.comparison_mz_grid[-1] - ref.comparison_mz_grid[0])
            overlap_span = overlap_max - overlap_min
            txt_coverage = overlap_span / txt_span if txt_span > 0 else 0
            ref_coverage = overlap_span / ref_span if ref_span > 0 else 0
        else:
            blocks.append("NO_OVERLAPPING_MZ_RANGE")
    else:
        blocks.append("NO_OVERLAPPING_MZ_RANGE")
    coverage = min(txt_coverage, ref_coverage)
    if coverage < parameters.minimum_coverage_fraction:
        blocks.append("INSUFFICIENT_MZ_COVERAGE")
    matches = match_replicate_peaks(txt, ref, parameters=parameters.replicate_parameters)
    paired = []
    txt_matched_ids: set[str] = set()
    ref_matched_ids: set[str] = set()
    for match in matches:
        txt_id, ref_id, delta, ppm = _match_orientation(match, txt.run_label)
        if txt_id and ref_id:
            paired.append((match, txt_id, ref_id, delta, ppm))
            txt_matched_ids.add(txt_id)
            ref_matched_ids.add(ref_id)
    if not txt.peaks:
        blocks.append("NO_TXT_PEAKS")
    if not ref.peaks:
        blocks.append("NO_REFERENCE_PEAKS")
    if not paired:
        blocks.append("NO_MATCHED_PEAKS")
    ambiguous = sum(match.match_ambiguity_status is not MatchAmbiguityStatus.UNAMBIGUOUS_ONE_TO_ONE for match, *_ in paired)
    if ambiguous:
        blocks.append("AMBIGUOUS_PEAK_MATCH")
    strict = sum(abs(delta) <= max(parameters.replicate_parameters.strict_absolute_tolerance_da,
        next(peak.apex_mz for peak in txt.peaks if peak.peak_id == txt_id) * parameters.replicate_parameters.strict_ppm_tolerance * 1e-6) + 1e-12
        for _, txt_id, _, delta, _ in paired)
    txt_overlap = len(txt_matched_ids) / len(txt.peaks) if txt.peaks else 0
    ref_overlap = len(ref_matched_ids) / len(ref.peaks) if ref.peaks else 0
    union = len(txt.peaks) + len(ref.peaks) - len(paired)
    jaccard = len(paired) / union if union else 0
    if base_corr is None or base_corr < parameters.minimum_profile_correlation:
        blocks.append("LOW_PROFILE_CORRELATION")
    if txt_overlap < parameters.minimum_peak_overlap:
        blocks.append("LOW_PEAK_OVERLAP")
    ranked_txt = sorted(txt.peaks, key=lambda peak: (peak.intensity_rank, peak.apex_mz, peak.peak_id))
    top10 = sum(peak.peak_id in txt_matched_ids for peak in ranked_txt[:10]) / min(10, len(ranked_txt)) if ranked_txt else 0
    top25 = sum(peak.peak_id in txt_matched_ids for peak in ranked_txt[:25]) / min(25, len(ranked_txt)) if ranked_txt else 0
    deltas = np.asarray([delta for _, _, _, delta, _ in paired], dtype=float)
    med_delta = float(np.median(deltas)) if len(deltas) else None
    mad_delta = float(np.median(np.abs(deltas - med_delta))) if len(deltas) else None
    strict_fraction = strict / len(paired) if paired else 0
    score = (
        parameters.weight_profile_correlation * max(base_corr or 0, 0)
        + parameters.weight_spearman * max(rank_corr or 0, 0)
        + parameters.weight_cosine * max(cosine or 0, 0)
        + parameters.weight_txt_peak_overlap * txt_overlap
        + parameters.weight_top25_match * top25
        + parameters.weight_strict_match * strict_fraction
    )
    if reference_profile.hypothesis_type is HypothesisType.RUN_1_ONLY_EXPORT:
        status = LinkageStatus.STRONG_LINK_TO_RUN_1 if score >= parameters.strong_score else LinkageStatus.SUPPORTIVE_LINK_TO_RUN_1 if score >= parameters.supportive_score else LinkageStatus.DIFFERENT_PROCESSING_EXPORT_POSSIBLE if score >= parameters.no_support_score else LinkageStatus.NO_SUPPORTED_LINKAGE
    elif reference_profile.hypothesis_type is HypothesisType.RUN_2_ONLY_EXPORT:
        status = LinkageStatus.STRONG_LINK_TO_RUN_2 if score >= parameters.strong_score else LinkageStatus.SUPPORTIVE_LINK_TO_RUN_2 if score >= parameters.supportive_score else LinkageStatus.DIFFERENT_PROCESSING_EXPORT_POSSIBLE if score >= parameters.no_support_score else LinkageStatus.NO_SUPPORTED_LINKAGE
    elif reference_profile.hypothesis_type is HypothesisType.PARTIAL_SCAN_EXPORT_POSSIBLE:
        run2 = "RUN_2" in reference_profile.reference_run_labels[0]
        status = (LinkageStatus.POSSIBLE_PARTIAL_SCAN_EXPORT_RUN_2 if run2 else LinkageStatus.POSSIBLE_PARTIAL_SCAN_EXPORT_RUN_1) if score >= parameters.supportive_score else LinkageStatus.NO_SUPPORTED_LINKAGE
    else:
        status = LinkageStatus.BEST_MATCH_AGGREGATE_PROFILE if score >= parameters.supportive_score else LinkageStatus.NO_SUPPORTED_LINKAGE
    confidence = "HIGH" if score >= parameters.strong_score else "MEDIUM" if score >= parameters.supportive_score else "LOW"
    low_quality = sum(match.replicate_consistency_status is ReplicateConsistencyStatus.INSUFFICIENT_PEAK_QUALITY for match, *_ in paired)
    return TxtMzMLLinkageHypothesisResult(
        hypothesis_id=reference_profile.hypothesis_id, hypothesis_type=reference_profile.hypothesis_type,
        reference_run_labels=reference_profile.reference_run_labels,
        raw_profile_correlation=raw_corr, base_peak_normalized_correlation=base_corr,
        sum_normalized_correlation=sum_corr, log1p_correlation=log_corr,
        spearman_rank_correlation=rank_corr, cosine_similarity=cosine,
        comparison_grid_method="REFERENCE_MZML_GRID_OVERLAP_ONLY", comparison_grid_spacing=spacing,
        overlap_mz_min=overlap_min, overlap_mz_max=overlap_max,
        txt_coverage_fraction=txt_coverage, mzml_coverage_fraction=ref_coverage,
        interpolation_method="LINEAR_INTERPOLATION_TXT_TO_REFERENCE_GRID",
        extrapolation_applied=False, txt_peak_count=len(txt.peaks), reference_peak_count=len(ref.peaks),
        matched_peak_count=len(paired), txt_only_peak_count=len(txt.peaks) - len(txt_matched_ids),
        reference_only_peak_count=len(ref.peaks) - len(ref_matched_ids), ambiguous_match_count=ambiguous,
        strict_match_count=strict, supportive_match_count=len(paired) - strict,
        txt_overlap_fraction=txt_overlap, reference_overlap_fraction=ref_overlap, peak_jaccard=jaccard,
        median_absolute_delta_mz=_median(abs(delta) for _, _, _, delta, _ in paired),
        median_ppm_error=_median(ppm for _, _, _, _, ppm in paired), mad_delta_mz=mad_delta,
        top_10_txt_peak_match_fraction=top10, top_25_txt_peak_match_fraction=top25,
        discriminating_peak_support_count=0, discriminating_peak_conflict_count=0,
        low_quality_match_count=low_quality, coverage_fraction=coverage,
        composite_linkage_score=float(score), linkage_status=status, linkage_confidence=confidence,
        block_reasons=_ordered_blocks(blocks), window_start_time=reference_profile.window_start_time,
        window_end_time=reference_profile.window_end_time, window_ms1_scan_count=reference_profile.window_ms1_scan_count,
    )


def evaluate_partial_scan_hypotheses(
    txt_profile: T1TxtPeakProfile, mzml_path: Path, *,
    window_config: Mapping[str, Any],
) -> tuple[TxtMzMLLinkageHypothesisResult, ...]:
    config = dict(window_config)
    parameters = config.pop("parameters", None) or SourceLinkageParameters()
    metadata_record = config.pop("metadata_record", None)
    start = config.pop("start_time", None)
    end = config.pop("end_time", None)
    count = int(config.pop("window_count", parameters.partial_window_count))
    if config:
        raise ValueError(f"unsupported window_config keys: {sorted(config)}")
    if metadata_record is not None:
        start = start if start is not None else _safe_float(metadata_record.minimum_scan_start_time)
        end = end if end is not None else _safe_float(metadata_record.maximum_scan_start_time)
    if start is None or end is None or end <= start or count < 1 or count > parameters.maximum_windows:
        return ()
    boundaries = np.linspace(float(start), float(end), count + 1)
    label = metadata_record.technical_run_label if metadata_record is not None else Path(mzml_path).stem
    output: list[TxtMzMLLinkageHypothesisResult] = []
    buffer: list[Mapping[str, Any]] = []
    current = 0

    def flush(index: int) -> None:
        nonlocal buffer
        if len(buffer) < parameters.minimum_ms1_scans_per_window:
            buffer = []
            return
        window_label = f"{label}__WINDOW_{index + 1}"
        profile = build_ms1_peak_profile_from_spectra(
            buffer, run_label=window_label, input_path=str(mzml_path), metadata_record=metadata_record,
            parameters=parameters.replicate_parameters,
        )
        reference = SourceLinkageReferenceProfile(
            f"PARTIAL_SCAN__{label}__{index + 1}", HypothesisType.PARTIAL_SCAN_EXPORT_POSSIBLE,
            (label,), profile, float(boundaries[index]), float(boundaries[index + 1]),
            profile.ms1_spectra_used,
        )
        output.append(compare_txt_to_reference_profile(txt_profile, reference,
            matching_config={"parameters": parameters}))
        buffer = []

    try:
        with mzml.MzML(str(mzml_path), decode_binary=False) as reader:
            for spectrum in reader:
                try:
                    level = int(spectrum.get("ms level", 0))
                except (TypeError, ValueError):
                    continue
                if level != 1:
                    continue
                rt = _rt_minutes(spectrum)
                if rt is None or rt < boundaries[0] or rt > boundaries[-1]:
                    continue
                index = min(count - 1, max(0, int(np.searchsorted(boundaries, rt, side="right") - 1)))
                while current < index:
                    flush(current)
                    current += 1
                buffer.append(spectrum)
            flush(current)
    except Exception:
        return ()
    return tuple(sorted(output, key=lambda result: result.hypothesis_id))


def _nearest_txt_peak(txt_peaks: Sequence[ReplicateRunPeak], mz_value: float,
                      parameters: ReplicateAuditParameters) -> ReplicateRunPeak | None:
    ordered = sorted(txt_peaks, key=lambda peak: (peak.apex_mz, peak.peak_id))
    masses = [peak.apex_mz for peak in ordered]
    tolerance = max(parameters.absolute_tolerance_da, mz_value * parameters.ppm_tolerance * 1e-6)
    index = bisect_left(masses, mz_value - tolerance)
    candidates = []
    while index < len(ordered) and ordered[index].apex_mz <= mz_value + tolerance:
        candidates.append(ordered[index])
        index += 1
    return min(candidates, key=lambda peak: (abs(peak.apex_mz - mz_value), peak.intensity_rank, peak.peak_id)) if candidates else None


def build_discriminating_evidence(
    txt_profile: T1TxtPeakProfile, replicate_result: ReplicateConsistencyAuditResult,
    *, parameters: SourceLinkageParameters | None = None,
) -> tuple[DiscriminatingPeakEvidence, ...]:
    params = parameters or SourceLinkageParameters()
    if len(replicate_result.run_profiles) < 2:
        return ()
    run1, run2 = replicate_result.run_profiles[:2]
    maps = ({peak.peak_id: peak for peak in run1.peaks}, {peak.peak_id: peak for peak in run2.peaks})
    output: list[DiscriminatingPeakEvidence] = []
    for match in replicate_result.matches:
        peak1 = maps[0].get(match.run_1_peak_id or "")
        peak2 = maps[1].get(match.run_2_peak_id or "")
        kind = support = ""
        observed = None
        if peak1 is not None and peak2 is None:
            kind, support, observed = "RUN_1_ONLY", HypothesisType.RUN_1_ONLY_EXPORT.value, peak1.apex_mz
        elif peak2 is not None and peak1 is None:
            kind, support, observed = "RUN_2_ONLY", HypothesisType.RUN_2_ONLY_EXPORT.value, peak2.apex_mz
        elif peak1 is not None and peak2 is not None and match.log2_normalized_intensity_ratio is not None and abs(match.log2_normalized_intensity_ratio) >= params.replicate_parameters.extreme_log2_intensity_ratio:
            if match.log2_normalized_intensity_ratio > 0:
                kind, support = "RUN_2_STRONGER_INTENSITY", HypothesisType.RUN_2_ONLY_EXPORT.value
            else:
                kind, support = "RUN_1_STRONGER_INTENSITY", HypothesisType.RUN_1_ONLY_EXPORT.value
            observed = (peak1.apex_mz + peak2.apex_mz) / 2
        else:
            continue
        present = [peak for peak in (peak1, peak2) if peak is not None]
        if not present:
            continue
        if max(peak.scan_recurrence_fraction for peak in present) < params.replicate_parameters.supportive_recurrence_fraction:
            continue
        if max((peak.relative_prominence or 0) for peak in present) < params.replicate_parameters.minimum_relative_prominence:
            continue
        if all(peak.detection_status in _LOW_QUALITY for peak in present):
            continue
        txt = _nearest_txt_peak(txt_profile.comparison_profile.peaks, float(observed), params.replicate_parameters)
        recurrence = max(peak.scan_recurrence_fraction for peak in present)
        confidence = "HIGH" if recurrence >= 0.10 and txt is not None else "MEDIUM" if txt is not None else "LOW"
        output.append(DiscriminatingPeakEvidence(
            observed_mz=float(observed), discrimination_type=kind,
            run_1_status="DETECTED" if peak1 else "NOT_DETECTED",
            run_2_status="DETECTED" if peak2 else "NOT_DETECTED",
            run_1_normalized_intensity=peak1.normalized_apex_intensity if peak1 else None,
            run_2_normalized_intensity=peak2.normalized_apex_intensity if peak2 else None,
            run_1_rank=peak1.intensity_rank if peak1 else None, run_2_rank=peak2.intensity_rank if peak2 else None,
            run_1_recurrence=peak1.scan_recurrence_fraction if peak1 else None,
            run_2_recurrence=peak2.scan_recurrence_fraction if peak2 else None,
            txt_matched=txt is not None, txt_intensity=txt.raw_apex_intensity if txt else None,
            txt_rank=txt.intensity_rank if txt else None,
            supports_hypothesis=support if txt else "NO_TXT_SUPPORT", discrimination_confidence=confidence,
        ))
    return tuple(sorted(output, key=lambda item: (item.observed_mz, item.discrimination_type)))


def _add_discrimination_counts(
    results: Sequence[TxtMzMLLinkageHypothesisResult],
    evidence: Sequence[DiscriminatingPeakEvidence],
) -> tuple[TxtMzMLLinkageHypothesisResult, ...]:
    run1_support = sum(item.txt_matched and item.supports_hypothesis == HypothesisType.RUN_1_ONLY_EXPORT.value for item in evidence)
    run2_support = sum(item.txt_matched and item.supports_hypothesis == HypothesisType.RUN_2_ONLY_EXPORT.value for item in evidence)
    output = []
    for result in results:
        if result.hypothesis_type is HypothesisType.RUN_1_ONLY_EXPORT:
            support, conflict = run1_support, run2_support
        elif result.hypothesis_type is HypothesisType.RUN_2_ONLY_EXPORT:
            support, conflict = run2_support, run1_support
        elif result.hypothesis_type is HypothesisType.PARTIAL_SCAN_EXPORT_POSSIBLE:
            is_run_2 = any("RUN_2" in label for label in result.reference_run_labels)
            support, conflict = (run2_support, run1_support) if is_run_2 else (run1_support, run2_support)
        else:
            support, conflict = run1_support + run2_support, min(run1_support, run2_support)
        blocks = list(result.block_reasons)
        if support < 2:
            blocks.append("INSUFFICIENT_DISCRIMINATING_PEAKS")
        if conflict:
            blocks.append("CONFLICTING_DISCRIMINATING_PEAKS")
        output.append(replace(result, discriminating_peak_support_count=support,
                              discriminating_peak_conflict_count=conflict,
                              block_reasons=_ordered_blocks(blocks)))
    return tuple(output)


def _replicate_status_by_peak(result: ReplicateConsistencyAuditResult | None) -> dict[str, str]:
    output: dict[str, str] = {}
    if result is None:
        return output
    for match in result.matches:
        for peak_id in (match.run_1_peak_id, match.run_2_peak_id):
            if peak_id:
                output[peak_id] = match.replicate_consistency_status.value
    return output


def build_peak_evidence(
    txt_profile: T1TxtPeakProfile, references: Sequence[SourceLinkageReferenceProfile],
    replicate_result: ReplicateConsistencyAuditResult | None,
    *, parameters: SourceLinkageParameters | None = None,
) -> tuple[TxtPeakLinkageEvidence, ...]:
    params = parameters or SourceLinkageParameters()
    txt_peaks = txt_profile.comparison_profile.peaks
    run_refs = [reference for reference in references if reference.hypothesis_type in {HypothesisType.RUN_1_ONLY_EXPORT, HypothesisType.RUN_2_ONLY_EXPORT}]
    aggregate_refs = [reference for reference in references if reference.hypothesis_type in {HypothesisType.MEAN_OF_RUNS_EXPORT, HypothesisType.SUM_OF_RUNS_EXPORT, HypothesisType.MEDIAN_OF_RUNS_EXPORT}]
    status_map = _replicate_status_by_peak(replicate_result)
    match_maps: list[dict[str, tuple[ReplicateRunPeak, float, float]]] = []
    for reference in run_refs:
        reference_map = {peak.peak_id: peak for peak in reference.run_profile.peaks}
        mapping = {}
        for match in match_replicate_peaks(txt_profile.comparison_profile, reference.run_profile, parameters=params.replicate_parameters):
            txt_id, ref_id, delta, ppm = _match_orientation(match, txt_profile.comparison_profile.run_label)
            if txt_id and ref_id:
                mapping[txt_id] = (reference_map[ref_id], float(delta), float(ppm))
        match_maps.append(mapping)
    aggregate_matches: dict[str, list[str]] = {peak.peak_id: [] for peak in txt_peaks}
    for reference in aggregate_refs:
        for match in match_replicate_peaks(txt_profile.comparison_profile, reference.run_profile, parameters=params.replicate_parameters):
            txt_id, ref_id, _, _ = _match_orientation(match, txt_profile.comparison_profile.run_label)
            if txt_id and ref_id:
                aggregate_matches[txt_id].append(reference.hypothesis_id)
    output = []
    for txt in sorted(txt_peaks, key=lambda peak: (peak.intensity_rank, peak.apex_mz, peak.peak_id)):
        one = match_maps[0].get(txt.peak_id) if len(match_maps) > 0 else None
        two = match_maps[1].get(txt.peak_id) if len(match_maps) > 1 else None
        if one and two:
            linkage, confidence = "MATCHED_BOTH_RUNS", "MEDIUM"
        elif one:
            linkage, confidence = "MATCHED_RUN_1_ONLY", "MEDIUM"
        elif two:
            linkage, confidence = "MATCHED_RUN_2_ONLY", "MEDIUM"
        elif aggregate_matches[txt.peak_id]:
            linkage, confidence = "MATCHED_AGGREGATE_ONLY", "LOW"
        else:
            linkage, confidence = "UNMATCHED", "LOW"
        output.append(TxtPeakLinkageEvidence(
            txt_peak_id=txt.peak_id, txt_apex_mz=txt.apex_mz, txt_intensity=txt.raw_apex_intensity,
            txt_normalized_intensity=txt.normalized_apex_intensity, txt_intensity_rank=txt.intensity_rank,
            run_1_matched=one is not None, run_1_peak_id=one[0].peak_id if one else None,
            run_1_delta_mz=one[1] if one else None, run_1_ppm_error=one[2] if one else None,
            run_1_intensity_rank=one[0].intensity_rank if one else None,
            run_1_scan_recurrence=one[0].scan_recurrence_fraction if one else None,
            run_1_replicate_status=status_map.get(one[0].peak_id, "NOT_MATCHED") if one else "NOT_MATCHED",
            run_2_matched=two is not None, run_2_peak_id=two[0].peak_id if two else None,
            run_2_delta_mz=two[1] if two else None, run_2_ppm_error=two[2] if two else None,
            run_2_intensity_rank=two[0].intensity_rank if two else None,
            run_2_scan_recurrence=two[0].scan_recurrence_fraction if two else None,
            run_2_replicate_status=status_map.get(two[0].peak_id, "NOT_MATCHED") if two else "NOT_MATCHED",
            aggregate_hypothesis_matches=tuple(sorted(aggregate_matches[txt.peak_id])),
            peak_linkage_status=linkage, peak_linkage_confidence=confidence, block_reasons=(),
        ))
    return tuple(output)


def summarize_txt_mzml_source_linkage(
    hypothesis_results: Sequence[TxtMzMLLinkageHypothesisResult], *,
    txt_file: str = "UNKNOWN", source_metadata_records: Sequence[MzMLSourceMetadataRecord] | None = None,
    replicate_audit_result: ReplicateConsistencyAuditResult | None = None,
    discriminating_evidence: Sequence[DiscriminatingPeakEvidence] = (),
    parameters: SourceLinkageParameters | None = None,
) -> TxtMzMLSourceLinkageSummary:
    params = parameters or SourceLinkageParameters()
    ordered = sorted(hypothesis_results, key=lambda result: (-result.composite_linkage_score, result.hypothesis_id))
    records = tuple(source_metadata_records or ())
    if not ordered:
        return TxtMzMLSourceLinkageSummary(
            txt_file=txt_file, rna_identity="UNKNOWN", digest_type="UNKNOWN", context_source="UNKNOWN",
            candidate_run_count=len(records), best_hypothesis="NONE", second_best_hypothesis="NONE",
            best_linkage_status=LinkageStatus.INSUFFICIENT_EVIDENCE, best_linkage_confidence="BLOCKED",
            best_composite_score=0, second_best_composite_score=0, score_margin=0,
            profile_evidence="NOT_AVAILABLE", peak_evidence="NOT_AVAILABLE",
            discriminating_peak_evidence="NOT_AVAILABLE", partial_scan_evidence="NOT_TESTED",
            aggregate_evidence="NOT_AVAILABLE", source_linkage_confirmed=False,
            exact_run_linkage_confirmed=False, source_polarity="UNKNOWN",
            source_polarity_evidence="NOT_AVAILABLE", common_source_polarity_supported=False,
            polarity_propagation_eligible=False, polarity_propagation_applied=False,
            polarity_propagation_block_reasons=("INSUFFICIENT_SOURCE_LINKAGE",),
            overall_block_reasons=("INSUFFICIENT_SOURCE_LINKAGE",),
        )
    best = ordered[0]
    second = ordered[1] if len(ordered) > 1 else None
    second_score = second.composite_linkage_score if second else 0.0
    margin = best.composite_linkage_score - second_score
    similar = second is not None and margin < params.similar_score_margin
    blocks = list(best.block_reasons)
    support = best.discriminating_peak_support_count
    conflict = best.discriminating_peak_conflict_count
    exact = (
        best.hypothesis_type in {HypothesisType.RUN_1_ONLY_EXPORT, HypothesisType.RUN_2_ONLY_EXPORT}
        and best.composite_linkage_score >= params.strong_score and margin >= params.unique_score_margin
        and support >= params.minimum_discriminating_support and conflict == 0
    )
    if similar:
        status, confidence = LinkageStatus.MULTIPLE_HYPOTHESES_SIMILAR, "MEDIUM" if best.composite_linkage_score >= params.supportive_score else "LOW"
        blocks.append("MULTIPLE_HYPOTHESES_SIMILAR")
    elif best.hypothesis_type is HypothesisType.PARTIAL_SCAN_EXPORT_POSSIBLE and best.composite_linkage_score >= params.supportive_score:
        run2 = "RUN_2" in best.reference_run_labels[0]
        status = LinkageStatus.POSSIBLE_PARTIAL_SCAN_EXPORT_RUN_2 if run2 else LinkageStatus.POSSIBLE_PARTIAL_SCAN_EXPORT_RUN_1
        confidence = "MEDIUM"
    elif best.hypothesis_type in {HypothesisType.MEAN_OF_RUNS_EXPORT, HypothesisType.SUM_OF_RUNS_EXPORT, HypothesisType.MEDIAN_OF_RUNS_EXPORT} and best.composite_linkage_score >= params.supportive_score:
        status, confidence = LinkageStatus.BEST_MATCH_AGGREGATE_PROFILE, "MEDIUM"
    elif best.hypothesis_type is HypothesisType.RUN_1_ONLY_EXPORT and best.composite_linkage_score >= params.supportive_score:
        status = LinkageStatus.STRONG_LINK_TO_RUN_1 if exact else LinkageStatus.SUPPORTIVE_LINK_TO_RUN_1
        confidence = "HIGH" if exact else "MEDIUM"
    elif best.hypothesis_type is HypothesisType.RUN_2_ONLY_EXPORT and best.composite_linkage_score >= params.supportive_score:
        status = LinkageStatus.STRONG_LINK_TO_RUN_2 if exact else LinkageStatus.SUPPORTIVE_LINK_TO_RUN_2
        confidence = "HIGH" if exact else "MEDIUM"
    elif best.composite_linkage_score >= params.no_support_score:
        status, confidence = LinkageStatus.DIFFERENT_PROCESSING_EXPORT_POSSIBLE, "LOW"
    else:
        status, confidence = LinkageStatus.NO_SUPPORTED_LINKAGE, "LOW"
        blocks.append("INSUFFICIENT_SOURCE_LINKAGE")
    polarities = {record.polarity_status for record in records}
    all_negative = bool(records) and polarities == {PolarityStatus.NEGATIVE_ONLY}
    manifest_ok = bool(records) and all(record.context_source == "USER_PROVIDED_RUNTIME_MANIFEST" and not record.context_conflict for record in records)
    common_polarity = all_negative and manifest_ok and best.composite_linkage_score >= params.minimum_polarity_support_score
    polarity_blocks: list[str] = []
    if len(polarities) > 1:
        polarity_blocks.append("SOURCE_POLARITY_CONFLICT")
    if any(record.context_conflict for record in records):
        polarity_blocks.append("USER_MANIFEST_METADATA_CONFLICT")
    if not common_polarity:
        polarity_blocks.append("INSUFFICIENT_SOURCE_LINKAGE")
    aggregate = [result for result in ordered if result.hypothesis_type in {HypothesisType.MEAN_OF_RUNS_EXPORT, HypothesisType.SUM_OF_RUNS_EXPORT, HypothesisType.MEDIAN_OF_RUNS_EXPORT}]
    partial = [result for result in ordered if result.hypothesis_type is HypothesisType.PARTIAL_SCAN_EXPORT_POSSIBLE]
    return TxtMzMLSourceLinkageSummary(
        txt_file=txt_file, rna_identity=records[0].rna_identity if records else "UNKNOWN",
        digest_type=records[0].digest_type if records else "UNKNOWN",
        context_source=records[0].context_source if records else "UNKNOWN", candidate_run_count=len(records),
        best_hypothesis=best.hypothesis_id, second_best_hypothesis=second.hypothesis_id if second else "NONE",
        best_linkage_status=status, best_linkage_confidence=confidence,
        best_composite_score=best.composite_linkage_score, second_best_composite_score=second_score,
        score_margin=margin,
        profile_evidence=f"BasePeak_Correlation={best.base_peak_normalized_correlation};Cosine={best.cosine_similarity}",
        peak_evidence=f"Txt_Overlap={best.txt_overlap_fraction};Jaccard={best.peak_jaccard}",
        discriminating_peak_evidence=f"Support={support};Conflict={conflict}",
        partial_scan_evidence=f"Tested={len(partial)};Best={max((item.composite_linkage_score for item in partial), default=0)}",
        aggregate_evidence=f"Tested={len(aggregate)};Best={max((item.composite_linkage_score for item in aggregate), default=0)};NonUnique=True",
        source_linkage_confirmed=exact, exact_run_linkage_confirmed=exact,
        source_polarity="NEGATIVE" if all_negative else "CONFLICT" if len(polarities) > 1 else "UNKNOWN",
        source_polarity_evidence="MZML_INTERNAL_METADATA" if records else "NOT_AVAILABLE",
        common_source_polarity_supported=common_polarity,
        polarity_propagation_eligible=common_polarity, polarity_propagation_applied=False,
        polarity_propagation_block_reasons=_ordered_blocks(polarity_blocks),
        overall_block_reasons=_ordered_blocks(blocks),
    )


def audit_t1_txt_mzml_source_linkage(
    txt_path: Path, replicate_paths: Sequence[Path], *,
    source_metadata_records: Sequence[MzMLSourceMetadataRecord] | None = None,
    replicate_audit_result: ReplicateConsistencyAuditResult | None = None,
    runtime_context: Mapping[str, Any] | None = None,
    parameters: SourceLinkageParameters | None = None,
) -> TxtMzMLSourceLinkageAuditResult:
    params = parameters or SourceLinkageParameters()
    params.validate()
    txt = parse_t1_txt_profile(Path(txt_path))
    txt_peaks = build_txt_profile_peaks(txt)
    metadata = tuple(source_metadata_records or ())
    if replicate_audit_result is None:
        from rna_masshunter.sciex_t1_replicate_consistency_audit import audit_sciex_t1_replicate_consistency
        replicate_audit_result = audit_sciex_t1_replicate_consistency(
            replicate_paths, source_metadata_records=metadata, runtime_context=runtime_context,
            parameters=params.replicate_parameters,
        )
    run_profiles = replicate_audit_result.run_profiles
    references = tuple(_reference_from_run(profile, index) for index, profile in enumerate(run_profiles))
    aggregates = build_replicate_aggregate_profiles(run_profiles)
    references += aggregates
    stage1 = tuple(compare_txt_to_reference_profile(txt_peaks, reference,
        matching_config={"parameters": params}) for reference in references)
    partial_results: tuple[TxtMzMLLinkageHypothesisResult, ...] = ()
    best_stage1 = max((result.composite_linkage_score for result in stage1), default=0)
    partial_refs: list[SourceLinkageReferenceProfile] = []
    if best_stage1 < params.supportive_score:
        metadata_by_name = {record.file_name: record for record in metadata}
        for path in replicate_paths:
            record = metadata_by_name.get(Path(path).name)
            results = evaluate_partial_scan_hypotheses(txt_peaks, Path(path), window_config={
                "parameters": params, "metadata_record": record,
            })
            partial_results += results
    discrimination = build_discriminating_evidence(txt_peaks, replicate_audit_result, parameters=params)
    results = _add_discrimination_counts(stage1 + partial_results, discrimination)
    peak_evidence = build_peak_evidence(txt_peaks, references, replicate_audit_result, parameters=params)
    summary = summarize_txt_mzml_source_linkage(
        results, txt_file=txt.file_name, source_metadata_records=metadata,
        replicate_audit_result=replicate_audit_result, discriminating_evidence=discrimination,
        parameters=params,
    )
    return TxtMzMLSourceLinkageAuditResult(
        params, txt, txt_peaks, references + tuple(partial_refs),
        tuple(sorted(results, key=lambda result: result.hypothesis_id)), peak_evidence,
        discrimination, summary,
    )


def _scalarize(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return ";".join(str(_scalarize(item)) for item in value)
    if isinstance(value, dict):
        return {key: _scalarize(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return "INTERNAL_ARRAY_NOT_EXPORTED"
    return value


def audit_optional_result(result: TxtMzMLSourceLinkageAuditResult) -> dict[str, Any]:
    metadata = asdict(result.txt_profile)
    metadata.pop("coordinates", None)
    metadata.pop("intensities", None)
    return {
        "txt_metadata_records": [_scalarize(metadata)],
        "hypothesis_records": [_scalarize(asdict(item)) for item in result.hypothesis_results],
        "peak_evidence_records": [_scalarize(asdict(item)) for item in result.peak_evidence],
        "discrimination_records": [_scalarize(asdict(item)) for item in result.discriminating_evidence],
        "summary_records": [_scalarize(asdict(result.summary))],
    }
