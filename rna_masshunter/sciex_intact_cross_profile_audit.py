"""Read-only shadow comparison of intact reconstructed peak profiles.

Inputs are established peak-family results. Matching is one-to-one and does not
assign RNA identity, biological origin, or any formal score.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from math import isfinite
from typing import Iterable

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr

from rna_masshunter.sciex_intact_peak_family import PeakFamilyPeak, SciexIntactPeakFamilyResult

ALGORITHM_VERSION = "sciex-intact-cross-profile-audit-v1"


class ComparisonLayer(str, Enum):
    ALL_DETECTED_PEAK_COMPARISON = "ALL_DETECTED_PEAK_COMPARISON"
    SELECTED_MAJOR_PEAK_COMPARISON = "SELECTED_MAJOR_PEAK_COMPARISON"


class CrossProfileMassMatchClass(str, Enum):
    STRICT = "STRICT"
    EXPLORATORY = "EXPLORATORY"


class ShapeSimilarityClass(str, Enum):
    HIGHLY_SIMILAR_PROFILE_PEAK = "HIGHLY_SIMILAR_PROFILE_PEAK"
    MODERATELY_SIMILAR_PROFILE_PEAK = "MODERATELY_SIMILAR_PROFILE_PEAK"
    MASS_MATCH_SHAPE_DIFFERENT = "MASS_MATCH_SHAPE_DIFFERENT"
    MASS_ONLY_MATCH = "MASS_ONLY_MATCH"


class SelectedPeakClassification(str, Enum):
    COMMON_SELECTED_PEAK = "COMMON_SELECTED_PEAK"
    UAA_ONLY_SELECTED_PEAK = "UAA_ONLY_SELECTED_PEAK"
    UAG_ONLY_SELECTED_PEAK = "UAG_ONLY_SELECTED_PEAK"
    AMBIGUOUS_CROSS_PROFILE_MATCH = "AMBIGUOUS_CROSS_PROFILE_MATCH"


@dataclass(frozen=True)
class CrossProfileParameters:
    strict_cross_profile_mass_tolerance_da: float = 0.5
    exploratory_cross_profile_mass_tolerance_da: float = 1.0
    ambiguity_error_margin_da: float = 0.05
    highly_similar_relative_log_ratio_max: float = 0.30
    moderately_similar_relative_log_ratio_max: float = 0.70
    highly_similar_centroid_difference_da: float = 0.5
    moderately_similar_centroid_difference_da: float = 1.0
    minimum_shape_fields: int = 4

    def validate(self) -> None:
        if self.strict_cross_profile_mass_tolerance_da <= 0:
            raise ValueError("strict cross-profile tolerance must be positive")
        if self.exploratory_cross_profile_mass_tolerance_da < self.strict_cross_profile_mass_tolerance_da:
            raise ValueError("exploratory tolerance must be >= strict tolerance")
        if self.ambiguity_error_margin_da < 0:
            raise ValueError("ambiguity margin cannot be negative")
        if self.minimum_shape_fields < 1:
            raise ValueError("minimum_shape_fields must be positive")


@dataclass(frozen=True, kw_only=True)
class CrossProfileSafeguards:
    shadow_analysis_only: bool = True
    mass_evidence_only: bool = True
    rna_identity_confirmed: bool = False
    target_rna_identity_confirmed_by_mass: bool = False
    both_target_trnas_assigned: bool = False
    common_peak_biological_identity_assigned: bool = False
    co_captured_rna_excluded: bool = False
    reconstruction_artifact_excluded: bool = False
    background_component_excluded: bool = False
    sequence_cocapture_ranking_performed: bool = False
    applied_to_formal_score: bool = False
    applied_to_ranking: bool = False
    applied_to_candidate_filtering: bool = False
    applied_to_final_consensus: bool = False


@dataclass(frozen=True, kw_only=True)
class CrossProfileMatch(CrossProfileSafeguards):
    cross_profile_match_id: str
    comparison_layer: ComparisonLayer
    uaa_peak_id: str
    uag_peak_id: str
    uaa_apex_mass: float
    uag_apex_mass: float
    apex_mass_difference_da: float
    uaa_centroid_mass: float | None
    uag_centroid_mass: float | None
    centroid_mass_difference_da: float | None
    uaa_relative_apex_intensity: float
    uag_relative_apex_intensity: float
    relative_apex_intensity_ratio: float | None
    uaa_relative_integrated_intensity: float | None
    uag_relative_integrated_intensity: float | None
    relative_integrated_intensity_ratio: float | None
    uaa_fwhm: float | None
    uag_fwhm: float | None
    fwhm_ratio: float | None
    uaa_peak_width: float | None
    uag_peak_width: float | None
    peak_width_ratio: float | None
    uaa_prominence: float
    uag_prominence: float
    prominence_ratio: float | None
    uaa_quality_class: str
    uag_quality_class: str
    mass_match_class: CrossProfileMassMatchClass
    shape_similarity_class: ShapeSimilarityClass
    ambiguous_assignment: bool
    hypotheses: tuple[str, ...] = (
        "COMMON_PROFILE_COMPONENT_POSSIBLE",
        "COMMON_CO_CAPTURED_RNA_POSSIBLE",
        "COMMON_RECONSTRUCTION_ARTIFACT_POSSIBLE",
        "COMMON_CONTAMINANT_OR_BACKGROUND_POSSIBLE",
        "TARGET_TRNA_SHARED_STATE_UNLIKELY_OR_UNRESOLVED",
    )


@dataclass(frozen=True)
class SelectedPeakStatus:
    profile: str
    peak_id: str
    apex_mass: float
    classification: SelectedPeakClassification
    matched_peak_id: str | None


@dataclass(frozen=True, kw_only=True)
class ProfileCommonPeakSummary(CrossProfileSafeguards):
    profile: str
    selected_peak_count: int
    common_selected_peak_count: int
    ambiguous_selected_peak_count: int
    sample_specific_selected_peak_count: int
    common_selected_apex_intensity_fraction: float
    sample_specific_apex_intensity_fraction: float
    common_selected_integrated_intensity_fraction: float | None
    sample_specific_integrated_intensity_fraction: float | None


@dataclass(frozen=True)
class CrossProfileCorrelations:
    method: str
    apex_mass: float | None
    relative_apex_intensity: float | None
    relative_integrated_intensity: float | None
    prominence: float | None
    fwhm: float | None


@dataclass(frozen=True)
class CrossProfileAuditResult:
    status: str
    parameters: CrossProfileParameters
    all_detected_matches: tuple[CrossProfileMatch, ...]
    selected_major_matches: tuple[CrossProfileMatch, ...]
    selected_peak_statuses: tuple[SelectedPeakStatus, ...]
    uaa_summary: ProfileCommonPeakSummary
    uag_summary: ProfileCommonPeakSummary
    correlations: CrossProfileCorrelations
    algorithm_version: str = ALGORITHM_VERSION


def _ordered(peaks: Iterable[PeakFamilyPeak]) -> tuple[PeakFamilyPeak, ...]:
    return tuple(sorted(peaks, key=lambda peak: (peak.apex_mass, peak.peak_id)))


def _ratio(left: float | None, right: float | None) -> float | None:
    if left is None or right is None or left <= 0 or right <= 0:
        return None
    return right / left


def _log_ratio_distance(ratio: float | None) -> float | None:
    if ratio is None:
        return None
    return abs(float(np.log10(ratio)))


def _shape_class(
    uaa: PeakFamilyPeak, uag: PeakFamilyPeak, centroid_difference: float | None,
    params: CrossProfileParameters,
) -> ShapeSimilarityClass:
    ratios = (
        _ratio(uaa.fwhm_da, uag.fwhm_da),
        _ratio(uaa.peak_width_da, uag.peak_width_da),
        _ratio(uaa.prominence, uag.prominence),
        _ratio(uaa.relative_apex_intensity, uag.relative_apex_intensity),
        _ratio(uaa.relative_integrated_intensity, uag.relative_integrated_intensity),
    )
    distances = [value for value in (_log_ratio_distance(ratio) for ratio in ratios) if value is not None]
    if len(distances) < params.minimum_shape_fields or centroid_difference is None:
        return ShapeSimilarityClass.MASS_ONLY_MATCH
    same_quality = uaa.peak_quality_class is uag.peak_quality_class
    if (same_quality and centroid_difference <= params.highly_similar_centroid_difference_da
            and max(distances) <= params.highly_similar_relative_log_ratio_max):
        return ShapeSimilarityClass.HIGHLY_SIMILAR_PROFILE_PEAK
    if (centroid_difference <= params.moderately_similar_centroid_difference_da
            and sum(value <= params.moderately_similar_relative_log_ratio_max for value in distances)
            >= max(3, len(distances) - 1)):
        return ShapeSimilarityClass.MODERATELY_SIMILAR_PROFILE_PEAK
    return ShapeSimilarityClass.MASS_MATCH_SHAPE_DIFFERENT


def _assign(
    uaa: tuple[PeakFamilyPeak, ...], uag: tuple[PeakFamilyPeak, ...],
    tolerance: float,
) -> tuple[tuple[int, int], ...]:
    """Maximum-cardinality, minimum-total-error deterministic assignment."""
    n, m = len(uaa), len(uag)
    if not n or not m:
        return ()
    size = n + m
    unmatched_cost = tolerance + 1.0
    forbidden_cost = unmatched_cost * (size + 2)
    cost = np.full((size, size), forbidden_cost, dtype=float)
    for i, left in enumerate(uaa):
        for j, right in enumerate(uag):
            error = abs(right.apex_mass - left.apex_mass)
            if error <= tolerance:
                cost[i, j] = error + (i * (m + 1) + j) * 1e-12
        cost[i, m + i] = unmatched_cost
    for j in range(m):
        cost[n + j, j] = unmatched_cost
    cost[n:, m:] = 0.0
    rows, cols = linear_sum_assignment(cost)
    pairs = [(int(i), int(j)) for i, j in zip(rows, cols)
             if i < n and j < m and cost[i, j] < unmatched_cost]
    return tuple(sorted(pairs, key=lambda pair: (uaa[pair[0]].apex_mass, uaa[pair[0]].peak_id,
                                                  uag[pair[1]].peak_id)))


def _ambiguous(
    index: int, assigned_other: int, own: tuple[PeakFamilyPeak, ...],
    other: tuple[PeakFamilyPeak, ...], params: CrossProfileParameters,
) -> bool:
    assigned_error = abs(other[assigned_other].apex_mass - own[index].apex_mass)
    alternatives = sorted(
        abs(candidate.apex_mass - own[index].apex_mass)
        for other_index, candidate in enumerate(other)
        if other_index != assigned_other
        and abs(candidate.apex_mass - own[index].apex_mass)
        <= params.exploratory_cross_profile_mass_tolerance_da
    )
    return bool(alternatives and alternatives[0] <= assigned_error + params.ambiguity_error_margin_da)


def match_cross_profile_peaks(
    uaa_peaks: Iterable[PeakFamilyPeak], uag_peaks: Iterable[PeakFamilyPeak],
    *, comparison_layer: ComparisonLayer,
    parameters: CrossProfileParameters | None = None,
) -> tuple[CrossProfileMatch, ...]:
    params = parameters or CrossProfileParameters(); params.validate()
    uaa, uag = _ordered(uaa_peaks), _ordered(uag_peaks)
    assignments = _assign(uaa, uag, params.exploratory_cross_profile_mass_tolerance_da)
    matches = []
    for uaa_index, uag_index in assignments:
        left, right = uaa[uaa_index], uag[uag_index]
        apex_difference = abs(right.apex_mass - left.apex_mass)
        centroid_difference = (abs(right.centroid_mass - left.centroid_mass)
                               if left.centroid_mass is not None and right.centroid_mass is not None else None)
        ambiguous = (
            _ambiguous(uaa_index, uag_index, uaa, uag, params)
            or _ambiguous(uag_index, uaa_index, uag, uaa, params)
        )
        matches.append(CrossProfileMatch(
            cross_profile_match_id=f"CROSS__{comparison_layer.value}__{left.peak_id}__{right.peak_id}",
            comparison_layer=comparison_layer, uaa_peak_id=left.peak_id, uag_peak_id=right.peak_id,
            uaa_apex_mass=left.apex_mass, uag_apex_mass=right.apex_mass,
            apex_mass_difference_da=apex_difference,
            uaa_centroid_mass=left.centroid_mass, uag_centroid_mass=right.centroid_mass,
            centroid_mass_difference_da=centroid_difference,
            uaa_relative_apex_intensity=left.relative_apex_intensity,
            uag_relative_apex_intensity=right.relative_apex_intensity,
            relative_apex_intensity_ratio=_ratio(left.relative_apex_intensity, right.relative_apex_intensity),
            uaa_relative_integrated_intensity=left.relative_integrated_intensity,
            uag_relative_integrated_intensity=right.relative_integrated_intensity,
            relative_integrated_intensity_ratio=_ratio(left.relative_integrated_intensity, right.relative_integrated_intensity),
            uaa_fwhm=left.fwhm_da, uag_fwhm=right.fwhm_da,
            fwhm_ratio=_ratio(left.fwhm_da, right.fwhm_da),
            uaa_peak_width=left.peak_width_da, uag_peak_width=right.peak_width_da,
            peak_width_ratio=_ratio(left.peak_width_da, right.peak_width_da),
            uaa_prominence=left.prominence, uag_prominence=right.prominence,
            prominence_ratio=_ratio(left.prominence, right.prominence),
            uaa_quality_class=left.peak_quality_class.value,
            uag_quality_class=right.peak_quality_class.value,
            mass_match_class=(CrossProfileMassMatchClass.STRICT
                              if apex_difference <= params.strict_cross_profile_mass_tolerance_da
                              else CrossProfileMassMatchClass.EXPLORATORY),
            shape_similarity_class=_shape_class(left, right, centroid_difference, params),
            ambiguous_assignment=ambiguous,
        ))
    return tuple(matches)


def _statuses(
    uaa: tuple[PeakFamilyPeak, ...], uag: tuple[PeakFamilyPeak, ...],
    matches: tuple[CrossProfileMatch, ...],
) -> tuple[SelectedPeakStatus, ...]:
    uaa_matches = {match.uaa_peak_id: match for match in matches}
    uag_matches = {match.uag_peak_id: match for match in matches}
    rows = []
    for peak in uaa:
        match = uaa_matches.get(peak.peak_id)
        klass = (SelectedPeakClassification.AMBIGUOUS_CROSS_PROFILE_MATCH if match and match.ambiguous_assignment
                 else SelectedPeakClassification.COMMON_SELECTED_PEAK if match
                 else SelectedPeakClassification.UAA_ONLY_SELECTED_PEAK)
        rows.append(SelectedPeakStatus("UAA", peak.peak_id, peak.apex_mass, klass,
                                       match.uag_peak_id if match else None))
    for peak in uag:
        match = uag_matches.get(peak.peak_id)
        klass = (SelectedPeakClassification.AMBIGUOUS_CROSS_PROFILE_MATCH if match and match.ambiguous_assignment
                 else SelectedPeakClassification.COMMON_SELECTED_PEAK if match
                 else SelectedPeakClassification.UAG_ONLY_SELECTED_PEAK)
        rows.append(SelectedPeakStatus("UAG", peak.peak_id, peak.apex_mass, klass,
                                       match.uaa_peak_id if match else None))
    return tuple(rows)


def _fraction(values: list[float | None], selected: list[bool]) -> tuple[float | None, float | None]:
    usable = [(float(value), flag) for value, flag in zip(values, selected) if value is not None]
    total = sum(value for value, _ in usable)
    if not usable or total <= 0:
        return None, None
    common = sum(value for value, flag in usable if flag)
    return common / total, (total - common) / total


def _summary(profile: str, peaks: tuple[PeakFamilyPeak, ...], statuses: tuple[SelectedPeakStatus, ...]) -> ProfileCommonPeakSummary:
    own = [status for status in statuses if status.profile == profile]
    common_flags = [status.classification in {
        SelectedPeakClassification.COMMON_SELECTED_PEAK,
        SelectedPeakClassification.AMBIGUOUS_CROSS_PROFILE_MATCH,
    } for status in own]
    apex_common, apex_specific = _fraction([peak.relative_apex_intensity for peak in peaks], common_flags)
    integrated_common, integrated_specific = _fraction(
        [peak.relative_integrated_intensity for peak in peaks], common_flags,
    )
    ambiguous = sum(status.classification is SelectedPeakClassification.AMBIGUOUS_CROSS_PROFILE_MATCH for status in own)
    common = sum(common_flags)
    return ProfileCommonPeakSummary(
        profile=profile, selected_peak_count=len(peaks), common_selected_peak_count=common,
        ambiguous_selected_peak_count=ambiguous, sample_specific_selected_peak_count=len(peaks)-common,
        common_selected_apex_intensity_fraction=apex_common or 0.0,
        sample_specific_apex_intensity_fraction=apex_specific or 0.0,
        common_selected_integrated_intensity_fraction=integrated_common,
        sample_specific_integrated_intensity_fraction=integrated_specific,
    )


def _spearman(left: list[float | None], right: list[float | None]) -> float | None:
    pairs = [(float(a), float(b)) for a, b in zip(left, right) if a is not None and b is not None]
    if len(pairs) < 3 or len({a for a, _ in pairs}) < 2 or len({b for _, b in pairs}) < 2:
        return None
    value = float(spearmanr([a for a, _ in pairs], [b for _, b in pairs]).statistic)
    return value if isfinite(value) else None


def _correlations(matches: tuple[CrossProfileMatch, ...]) -> CrossProfileCorrelations:
    return CrossProfileCorrelations(
        method="SPEARMAN_RANK_CORRELATION",
        apex_mass=_spearman([m.uaa_apex_mass for m in matches], [m.uag_apex_mass for m in matches]),
        relative_apex_intensity=_spearman([m.uaa_relative_apex_intensity for m in matches], [m.uag_relative_apex_intensity for m in matches]),
        relative_integrated_intensity=_spearman([m.uaa_relative_integrated_intensity for m in matches], [m.uag_relative_integrated_intensity for m in matches]),
        prominence=_spearman([m.uaa_prominence for m in matches], [m.uag_prominence for m in matches]),
        fwhm=_spearman([m.uaa_fwhm for m in matches], [m.uag_fwhm for m in matches]),
    )


def audit_leu_cross_profiles(
    uaa_result: SciexIntactPeakFamilyResult,
    uag_result: SciexIntactPeakFamilyResult,
    *, parameters: CrossProfileParameters | None = None,
) -> CrossProfileAuditResult:
    """Compare all-detected and selected-major layers without redetection."""
    params = parameters or CrossProfileParameters(); params.validate()
    all_matches = match_cross_profile_peaks(
        uaa_result.peaks, uag_result.peaks,
        comparison_layer=ComparisonLayer.ALL_DETECTED_PEAK_COMPARISON, parameters=params,
    )
    selected_matches = match_cross_profile_peaks(
        uaa_result.selected_peaks, uag_result.selected_peaks,
        comparison_layer=ComparisonLayer.SELECTED_MAJOR_PEAK_COMPARISON, parameters=params,
    )
    uaa_selected, uag_selected = _ordered(uaa_result.selected_peaks), _ordered(uag_result.selected_peaks)
    statuses = _statuses(uaa_selected, uag_selected, selected_matches)
    return CrossProfileAuditResult(
        "COMPLETED", params, all_matches, selected_matches, statuses,
        _summary("UAA", uaa_selected, statuses), _summary("UAG", uag_selected, statuses),
        _correlations(selected_matches),
    )
