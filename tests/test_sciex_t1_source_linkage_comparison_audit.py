from types import SimpleNamespace

from rna_masshunter.sciex_t1_source_linkage_comparison_audit import (
    T1SourceLinkageComparisonStatus,
    audit_optional_result,
    compare_t1_source_linkage_audits,
)
from rna_masshunter.sciex_t1_txt_mzml_source_linkage_audit import LinkageStatus


def audit(score, *, confirmed=False, status=LinkageStatus.NO_SUPPORTED_LINKAGE, label="BEST"):
    summary = SimpleNamespace(
        rna_identity="tRNA", candidate_run_count=1, best_hypothesis=label,
        best_composite_score=score, second_best_composite_score=0.1,
        score_margin=score - 0.1, best_linkage_status=status,
        best_linkage_confidence="HIGH" if confirmed else "LOW",
        source_linkage_confirmed=confirmed, exact_run_linkage_confirmed=confirmed,
        common_source_polarity_supported=confirmed, polarity_propagation_eligible=confirmed,
        overall_block_reasons=() if confirmed else ("INSUFFICIENT_SOURCE_LINKAGE",),
    )
    best = SimpleNamespace(
        hypothesis_id=label, base_peak_normalized_correlation=score,
        spearman_rank_correlation=score, cosine_similarity=score,
        peak_jaccard=score / 2, txt_overlap_fraction=score / 2,
        top_10_txt_peak_match_fraction=score, top_25_txt_peak_match_fraction=score,
        median_absolute_delta_mz=0.002,
    )
    return SimpleNamespace(
        summary=summary, hypothesis_results=(best,),
        parameters=SimpleNamespace(supportive_score=0.55, unique_score_margin=0.10),
    )


def test_uag_linkage_stronger_than_uaa():
    result = compare_t1_source_linkage_audits(
        audit(0.328), audit(0.90, confirmed=True, status=LinkageStatus.STRONG_LINK_TO_REFERENCE_RUN),
    )
    assert result.summary.comparison_status is T1SourceLinkageComparisonStatus.UAG_LINKAGE_STRONGER_THAN_UAA
    assert result.summary.sample_specific_annotations == (
        "UAA_SPECIFIC_SOURCE_COMPLEXITY_POSSIBLE",
        "UAA_SPECIFIC_EXPORT_DIFFERENCE_POSSIBLE",
        "UAA_MULTIPLE_RUN_AMBIGUITY_POSSIBLE",
    )
    assert not result.summary.mixture_confirmed and not result.summary.purity_assigned


def test_both_linkages_weak_adds_only_possible_common_processing_annotations():
    result = compare_t1_source_linkage_audits(audit(0.30), audit(0.40))
    assert result.summary.comparison_status is T1SourceLinkageComparisonStatus.BOTH_LINKAGES_WEAK
    assert "COMMON_TXT_MZML_PROCESSING_DIFFERENCE_POSSIBLE" in result.summary.common_processing_difference_annotations
    assert all(value.endswith("_POSSIBLE") for value in result.summary.common_processing_difference_annotations)
    assert result.summary.mixture_compatible_interpretation == "POSSIBLE"
    assert not result.summary.mixture_confirmed and not result.summary.purity_assigned


def test_comparison_is_deterministic_and_formally_nonpropagating():
    one = compare_t1_source_linkage_audits(audit(0.30), audit(0.40))
    two = compare_t1_source_linkage_audits(
        audit(0.30), audit(0.40), left_sample_label="UAA", right_sample_label="UAG",
    )
    assert one == two
    payload = audit_optional_result(one)
    assert set(payload) == {"comparison_records", "comparison_summary_records"}
    for records in payload.values():
        for row in records:
            assert row["shadow_analysis_only"] and row["source_linkage_audit_only"]
            for key in (
                "formal_propagation", "polarity_propagation_applied",
                "chemical_identity_assigned", "fragment_identity_assigned",
                "charge_state_confirmed", "purity_assigned",
            ):
                assert not row[key]


def test_missing_uaa_result_is_not_comparable():
    result = compare_t1_source_linkage_audits(None, audit(0.80, confirmed=True))
    assert result.summary.comparison_status is T1SourceLinkageComparisonStatus.NOT_COMPARABLE
    assert result.summary.block_reasons == ("UAA_COMPARISON_RESULT_MISSING",)
