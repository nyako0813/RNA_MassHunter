"""Shadow-only comparison of two independent T1 txt-to-mzML linkage audits."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from typing import Any

from rna_masshunter.sciex_t1_txt_mzml_source_linkage_audit import TxtMzMLSourceLinkageAuditResult

OPTIONAL_RESULT_KEY = "sciex_t1_source_linkage_comparison_audit"
ALGORITHM_VERSION = "sciex-t1-source-linkage-comparison-audit-v1"

_BLOCK_ORDER = ("UAA_COMPARISON_RESULT_MISSING", "UAA_UAG_NOT_COMPARABLE")


class T1SourceLinkageComparisonStatus(str, Enum):
    UAG_LINKAGE_STRONGER_THAN_UAA = "UAG_LINKAGE_STRONGER_THAN_UAA"
    UAG_LINKAGE_SUPPORTIVE_UAA_UNRESOLVED = "UAG_LINKAGE_SUPPORTIVE_UAA_UNRESOLVED"
    BOTH_LINKAGES_SUPPORTIVE = "BOTH_LINKAGES_SUPPORTIVE"
    BOTH_LINKAGES_WEAK = "BOTH_LINKAGES_WEAK"
    UAA_LINKAGE_STRONGER_THAN_UAG = "UAA_LINKAGE_STRONGER_THAN_UAG"
    NOT_COMPARABLE = "NOT_COMPARABLE"


@dataclass(frozen=True, kw_only=True)
class T1SourceLinkageComparisonRecord:
    sample_label: str
    rna_identity: str
    candidate_run_count: int
    best_hypothesis: str
    best_composite_score: float
    second_best_score: float
    score_margin: float
    profile_correlation: float | None
    spearman_correlation: float | None
    cosine_similarity: float | None
    peak_jaccard: float
    txt_overlap: float
    top_10_match_fraction: float
    top_25_match_fraction: float
    median_absolute_delta_mz: float | None
    source_linkage_status: str
    source_linkage_confidence: str
    source_linkage_confirmed: bool
    exact_run_linkage_confirmed: bool
    common_source_polarity_supported: bool
    polarity_propagation_eligible: bool
    block_reasons: tuple[str, ...]
    shadow_analysis_only: bool = True
    source_linkage_audit_only: bool = True
    formal_propagation: bool = False
    polarity_propagation_applied: bool = False
    chemical_identity_assigned: bool = False
    fragment_identity_assigned: bool = False
    charge_state_confirmed: bool = False
    purity_assigned: bool = False


@dataclass(frozen=True, kw_only=True)
class T1SourceLinkageComparisonSummary:
    left_sample_label: str
    right_sample_label: str
    comparison_status: T1SourceLinkageComparisonStatus
    score_difference_right_minus_left: float
    common_processing_difference_annotations: tuple[str, ...]
    sample_specific_annotations: tuple[str, ...]
    mixture_compatible_interpretation: str
    mixture_confirmed: bool
    purity_assigned: bool
    block_reasons: tuple[str, ...]
    shadow_analysis_only: bool = True
    source_linkage_audit_only: bool = True
    formal_propagation: bool = False
    polarity_propagation_applied: bool = False
    chemical_identity_assigned: bool = False
    fragment_identity_assigned: bool = False
    charge_state_confirmed: bool = False


@dataclass(frozen=True)
class T1SourceLinkageComparisonResult:
    comparison_records: tuple[T1SourceLinkageComparisonRecord, ...]
    summary: T1SourceLinkageComparisonSummary
    algorithm_version: str = ALGORITHM_VERSION
    formal_propagation: bool = False


def _ordered_blocks(values: tuple[str, ...] | list[str]) -> tuple[str, ...]:
    found = set(values)
    return tuple(value for value in _BLOCK_ORDER if value in found) + tuple(sorted(found - set(_BLOCK_ORDER)))


def _best_hypothesis(audit: TxtMzMLSourceLinkageAuditResult) -> Any | None:
    candidates = [item for item in audit.hypothesis_results if item.hypothesis_id == audit.summary.best_hypothesis]
    return sorted(candidates, key=lambda item: item.hypothesis_id)[0] if candidates else None


def _record(audit: TxtMzMLSourceLinkageAuditResult, sample_label: str) -> T1SourceLinkageComparisonRecord:
    best = _best_hypothesis(audit)
    summary = audit.summary
    return T1SourceLinkageComparisonRecord(
        sample_label=sample_label, rna_identity=summary.rna_identity,
        candidate_run_count=summary.candidate_run_count, best_hypothesis=summary.best_hypothesis,
        best_composite_score=summary.best_composite_score,
        second_best_score=summary.second_best_composite_score, score_margin=summary.score_margin,
        profile_correlation=best.base_peak_normalized_correlation if best else None,
        spearman_correlation=best.spearman_rank_correlation if best else None,
        cosine_similarity=best.cosine_similarity if best else None,
        peak_jaccard=best.peak_jaccard if best else 0.0,
        txt_overlap=best.txt_overlap_fraction if best else 0.0,
        top_10_match_fraction=best.top_10_txt_peak_match_fraction if best else 0.0,
        top_25_match_fraction=best.top_25_txt_peak_match_fraction if best else 0.0,
        median_absolute_delta_mz=best.median_absolute_delta_mz if best else None,
        source_linkage_status=summary.best_linkage_status.value,
        source_linkage_confidence=summary.best_linkage_confidence,
        source_linkage_confirmed=summary.source_linkage_confirmed,
        exact_run_linkage_confirmed=summary.exact_run_linkage_confirmed,
        common_source_polarity_supported=summary.common_source_polarity_supported,
        polarity_propagation_eligible=summary.polarity_propagation_eligible,
        block_reasons=summary.overall_block_reasons,
    )


def compare_t1_source_linkage_audits(
    left: TxtMzMLSourceLinkageAuditResult | None,
    right: TxtMzMLSourceLinkageAuditResult | None,
    *,
    left_sample_label: str = "UAA",
    right_sample_label: str = "UAG",
) -> T1SourceLinkageComparisonResult:
    if left is None or right is None:
        blocks = ("UAA_COMPARISON_RESULT_MISSING",) if left is None else ("UAA_UAG_NOT_COMPARABLE",)
        summary = T1SourceLinkageComparisonSummary(
            left_sample_label=left_sample_label, right_sample_label=right_sample_label,
            comparison_status=T1SourceLinkageComparisonStatus.NOT_COMPARABLE,
            score_difference_right_minus_left=0.0, common_processing_difference_annotations=(),
            sample_specific_annotations=(), mixture_compatible_interpretation="POSSIBLE",
            mixture_confirmed=False, purity_assigned=False, block_reasons=_ordered_blocks(blocks),
        )
        return T1SourceLinkageComparisonResult((), summary)
    left_record = _record(left, left_sample_label)
    right_record = _record(right, right_sample_label)
    difference = right_record.best_composite_score - left_record.best_composite_score
    supportive = left.parameters.supportive_score
    unique_margin = left.parameters.unique_score_margin
    left_supportive = left_record.best_composite_score >= supportive
    right_supportive = right_record.best_composite_score >= supportive
    if right_record.source_linkage_confirmed and not left_record.source_linkage_confirmed:
        status = T1SourceLinkageComparisonStatus.UAG_LINKAGE_STRONGER_THAN_UAA
    elif right_supportive and not left_supportive:
        status = (T1SourceLinkageComparisonStatus.UAG_LINKAGE_STRONGER_THAN_UAA
                  if difference >= unique_margin else
                  T1SourceLinkageComparisonStatus.UAG_LINKAGE_SUPPORTIVE_UAA_UNRESOLVED)
    elif left_supportive and right_supportive:
        status = T1SourceLinkageComparisonStatus.BOTH_LINKAGES_SUPPORTIVE
    elif not left_supportive and not right_supportive:
        status = T1SourceLinkageComparisonStatus.BOTH_LINKAGES_WEAK
    elif left_supportive and not right_supportive:
        status = T1SourceLinkageComparisonStatus.UAA_LINKAGE_STRONGER_THAN_UAG
    else:
        status = T1SourceLinkageComparisonStatus.NOT_COMPARABLE
    common_annotations = ()
    specific_annotations = ()
    if status is T1SourceLinkageComparisonStatus.BOTH_LINKAGES_WEAK:
        common_annotations = (
            "COMMON_TXT_MZML_PROCESSING_DIFFERENCE_POSSIBLE",
            "COMMON_EXPORT_PIPELINE_DIFFERENCE_POSSIBLE",
            "MZ_RANGE_OR_SAMPLING_DIFFERENCE_POSSIBLE",
            "PROFILE_AGGREGATION_DIFFERENCE_POSSIBLE",
        )
    elif status in {
        T1SourceLinkageComparisonStatus.UAG_LINKAGE_STRONGER_THAN_UAA,
        T1SourceLinkageComparisonStatus.UAG_LINKAGE_SUPPORTIVE_UAA_UNRESOLVED,
    }:
        specific_annotations = (
            "UAA_SPECIFIC_SOURCE_COMPLEXITY_POSSIBLE",
            "UAA_SPECIFIC_EXPORT_DIFFERENCE_POSSIBLE",
            "UAA_MULTIPLE_RUN_AMBIGUITY_POSSIBLE",
        )
    summary = T1SourceLinkageComparisonSummary(
        left_sample_label=left_sample_label, right_sample_label=right_sample_label,
        comparison_status=status, score_difference_right_minus_left=difference,
        common_processing_difference_annotations=common_annotations,
        sample_specific_annotations=specific_annotations,
        mixture_compatible_interpretation="POSSIBLE", mixture_confirmed=False,
        purity_assigned=False, block_reasons=(),
    )
    records = tuple(sorted((left_record, right_record), key=lambda item: item.sample_label))
    return T1SourceLinkageComparisonResult(records, summary)


def _scalarize(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return ";".join(str(_scalarize(item)) for item in value)
    return value


def audit_optional_result(result: T1SourceLinkageComparisonResult) -> dict[str, Any]:
    return {
        "comparison_records": [_scalarize(asdict(item)) for item in result.comparison_records],
        "comparison_summary_records": [_scalarize(asdict(result.summary))],
    }
