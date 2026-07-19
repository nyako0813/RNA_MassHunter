from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from rna_masshunter.sciex_mzml_source_metadata_audit import (
    MzMLSourceMetadataRecord, PolarityStatus, ReadStatus, RepresentationStatus,
)
from rna_masshunter.sciex_t1_replicate_consistency_audit import (
    ReplicateAuditParameters, ReplicateConsistencyAuditResult,
    match_replicate_peaks, summarize_replicate_consistency,
)
from rna_masshunter.sciex_t1_txt_mzml_source_linkage_audit import (
    HypothesisType, LinkageStatus, SourceLinkageParameters, SourceLinkageReferenceProfile,
    TxtProfileType, _add_discrimination_counts, _derived_profile, _reference_from_run,
    audit_optional_result, build_discriminating_evidence, build_peak_evidence,
    build_replicate_aggregate_profiles, build_txt_profile_peaks,
    compare_txt_to_reference_profile, parse_t1_txt_profile,
    summarize_txt_mzml_source_linkage,
)


def write_profile(tmp_path, name="profile.txt", *, delimiter="\t", header=True,
                  points=None, comments=False):
    points = points if points is not None else [(100 + index * 0.001, 10 * np.exp(-0.5 * ((index - 500) / 8) ** 2)) for index in range(1001)]
    lines = []
    if comments:
        lines += ["# comment", ""]
    if header:
        lines.append(delimiter.join(("Mass/Charge", "Intensity")))
    lines += [delimiter.join((f"{mz:.9e}", f"{intensity:.9e}")) for mz, intensity in points]
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def gaussian_profile(label, centers, *, scale=1.0, low=100.0, high=103.0, step=0.001):
    grid = np.arange(low, high + step / 2, step)
    values = np.zeros_like(grid)
    for center, height in centers:
        values += height * scale * np.exp(-0.5 * ((grid - center) / 0.008) ** 2)
    normalized = values / values.max() if values.max() else values
    return _derived_profile(label, grid, values, normalized)


def with_scan_recurrence(profile, fraction=0.5):
    return replace(
        profile,
        peaks=tuple(
            replace(
                peak,
                supporting_ms1_scan_count=5,
                total_ms1_scan_count=10,
                scan_recurrence_fraction=fraction,
            )
            for peak in profile.peaks
        ),
    )


def txt_peak_profile(tmp_path, centers, *, scale=1.0):
    grid = np.arange(100.0, 103.0001, 0.001)
    values = np.zeros_like(grid)
    for center, height in centers:
        values += height * scale * np.exp(-0.5 * ((grid - center) / 0.008) ** 2)
    path = write_profile(tmp_path, points=list(zip(grid, values, strict=False)))
    return build_txt_profile_peaks(parse_t1_txt_profile(path))


def metadata(label, polarity=PolarityStatus.NEGATIVE_ONLY):
    return MzMLSourceMetadataRecord(
        input_path=f"{label}.mzML", file_name=f"{label}.mzML", read_status=ReadStatus.COMPLETED,
        technical_run_label=label, context_source="USER_PROVIDED_RUNTIME_MANIFEST",
        context_confidence="USER_CONFIRMED", rna_identity="tRNA^Leu-UAA", digest_type="T1_DIGEST",
        polarity_status=polarity, representation_status=RepresentationStatus.PROFILE_ONLY,
    )


@pytest.mark.parametrize("delimiter", [" ", "\t", ","])
@pytest.mark.parametrize("header", [True, False])
def test_txt_parsing_delimiters_headers_scientific_notation(tmp_path, delimiter, header):
    path = write_profile(tmp_path, delimiter=delimiter, header=header, comments=True,
                         points=[(100.0, 1e3), (100.1, 2.5e-4), (100.2, 3.0)])
    result = parse_t1_txt_profile(path)
    assert result.parse_status == "COMPLETED" and result.numeric_row_count == 3
    assert result.mz_column_index == 0 and result.intensity_column_index == 1
    assert result.header_status == ("PRESENT" if header else "ABSENT_TWO_COLUMN_ASSUMED")


def test_invalid_txt_empty_non_numeric_and_unknown_columns(tmp_path):
    empty = tmp_path / "empty.txt"; empty.write_text("", encoding="utf-8")
    bad = tmp_path / "bad.txt"; bad.write_text("foo bar baz\na b c\n", encoding="utf-8")
    short = tmp_path / "short.txt"; short.write_text("Mass/Charge\n100\n", encoding="utf-8")
    assert "TXT_EMPTY" in parse_t1_txt_profile(empty).parse_block_reasons
    assert "TXT_COLUMN_ASSIGNMENT_UNRESOLVED" in parse_t1_txt_profile(bad).parse_block_reasons
    assert "TXT_COLUMN_ASSIGNMENT_UNRESOLVED" in parse_t1_txt_profile(short).parse_block_reasons


def test_dense_profile_classification(tmp_path):
    result = parse_t1_txt_profile(write_profile(tmp_path))
    assert result.profile_or_peaklist_status is TxtProfileType.DENSE_PROFILE
    assert result.duplicate_mz_count == 0 and result.mz_sorted_status == "STRICTLY_INCREASING"


def test_sparse_peak_list_classification(tmp_path):
    points = [(100.0, 10), (101.0, 5), (102.0, 2)]
    result = parse_t1_txt_profile(write_profile(tmp_path, points=points))
    assert result.profile_or_peaklist_status is TxtProfileType.SPARSE_PEAK_LIST
    peaks = build_txt_profile_peaks(result)
    assert peaks.comparison_peak_count == 3


def test_run_1_strong_linkage(tmp_path):
    txt = txt_peak_profile(tmp_path, [(100.5, 10), (101.5, 6)])
    run = gaussian_profile("UAA_T1_RUN_1", [(100.5, 10), (101.5, 6)])
    result = compare_txt_to_reference_profile(txt, _reference_from_run(run, 0))
    assert result.linkage_status is LinkageStatus.STRONG_LINK_TO_RUN_1
    assert result.composite_linkage_score >= 0.75 and result.top_10_txt_peak_match_fraction == 1


def test_run_2_strong_linkage_and_intensity_scaling(tmp_path):
    txt = txt_peak_profile(tmp_path, [(100.5, 10), (101.5, 6)], scale=100)
    run = gaussian_profile("UAA_T1_RUN_2", [(100.5, 10), (101.5, 6)])
    result = compare_txt_to_reference_profile(txt, _reference_from_run(run, 1))
    assert result.linkage_status is LinkageStatus.STRONG_LINK_TO_RUN_2
    assert result.base_peak_normalized_correlation == pytest.approx(1, abs=1e-8)


def test_mz_rounding_within_standard_tolerance(tmp_path):
    txt = txt_peak_profile(tmp_path, [(100.504, 10)])
    run = gaussian_profile("UAA_T1_RUN_1", [(100.500, 10)])
    result = compare_txt_to_reference_profile(txt, _reference_from_run(run, 0))
    assert result.matched_peak_count == 1 and result.median_absolute_delta_mz <= 0.005


def test_similar_runs_are_not_forced_to_one(tmp_path):
    txt = txt_peak_profile(tmp_path, [(100.5, 10), (101.5, 6)])
    run1 = gaussian_profile("UAA_T1_RUN_1", [(100.5, 10), (101.5, 6)])
    run2 = gaussian_profile("UAA_T1_RUN_2", [(100.5, 10), (101.5, 6)])
    results = [compare_txt_to_reference_profile(txt, _reference_from_run(run1, 0)),
               compare_txt_to_reference_profile(txt, _reference_from_run(run2, 1))]
    summary = summarize_txt_mzml_source_linkage(results, txt_file="x", source_metadata_records=[metadata("R1"), metadata("R2")])
    assert summary.best_linkage_status is LinkageStatus.MULTIPLE_HYPOTHESES_SIMILAR
    assert not summary.exact_run_linkage_confirmed


def test_aggregate_profile_and_mean_sum_nonidentifiability(tmp_path):
    txt = txt_peak_profile(tmp_path, [(100.5, 10), (102.5, 10)])
    run1 = gaussian_profile("UAA_T1_RUN_1", [(100.5, 10)])
    run2 = gaussian_profile("UAA_T1_RUN_2", [(102.5, 10)])
    aggregates = build_replicate_aggregate_profiles([run1, run2])
    assert len(aggregates) == 3
    results = [compare_txt_to_reference_profile(txt, reference) for reference in aggregates]
    assert all("AGGREGATE_HYPOTHESIS_NON_UNIQUE" in result.block_reasons for result in results)
    assert all(result.linkage_status is LinkageStatus.BEST_MATCH_AGGREGATE_PROFILE for result in results)
    assert max(result.composite_linkage_score for result in results) >= 0.75
    assert results[0].base_peak_normalized_correlation == pytest.approx(results[1].base_peak_normalized_correlation)


def test_partial_window_result_can_be_best_without_confirming_exact_run(tmp_path):
    txt = txt_peak_profile(tmp_path, [(100.5, 10)])
    whole = gaussian_profile("UAA_T1_RUN_1", [(102.0, 10)])
    window = gaussian_profile("UAA_T1_RUN_1_WINDOW", [(100.5, 10)])
    whole_result = compare_txt_to_reference_profile(txt, _reference_from_run(whole, 0))
    partial_ref = SourceLinkageReferenceProfile("PARTIAL", HypothesisType.PARTIAL_SCAN_EXPORT_POSSIBLE,
        ("UAA_T1_RUN_1",), window, 0, 1, 10)
    partial_result = compare_txt_to_reference_profile(txt, partial_ref)
    summary = summarize_txt_mzml_source_linkage([whole_result, partial_result], txt_file="x",
        source_metadata_records=[metadata("R1")])
    assert summary.best_linkage_status is LinkageStatus.POSSIBLE_PARTIAL_SCAN_EXPORT_RUN_1
    assert not summary.exact_run_linkage_confirmed


def test_discriminating_run_1_peak_support(tmp_path):
    txt = txt_peak_profile(tmp_path, [(100.5, 10)])
    run1 = with_scan_recurrence(gaussian_profile("UAA_T1_RUN_1", [(100.5, 10)]))
    run2 = with_scan_recurrence(gaussian_profile("UAA_T1_RUN_2", [(102.0, 10)]))
    matches = match_replicate_peaks(run1, run2)
    summary = summarize_replicate_consistency(run1, run2, matches)
    replicate = ReplicateConsistencyAuditResult(ReplicateAuditParameters(), (run1, run2), matches, (summary,))
    evidence = build_discriminating_evidence(txt, replicate)
    assert any(item.txt_matched and item.supports_hypothesis == HypothesisType.RUN_1_ONLY_EXPORT.value for item in evidence)


def test_conflicting_discriminating_counts_prevent_unique_support(tmp_path):
    txt = txt_peak_profile(tmp_path, [(100.5, 10), (102.0, 10)])
    run1 = with_scan_recurrence(gaussian_profile("UAA_T1_RUN_1", [(100.5, 10)]))
    run2 = with_scan_recurrence(gaussian_profile("UAA_T1_RUN_2", [(102.0, 10)]))
    matches = match_replicate_peaks(run1, run2)
    replicate = ReplicateConsistencyAuditResult(ReplicateAuditParameters(), (run1, run2), matches,
        (summarize_replicate_consistency(run1, run2, matches),))
    evidence = build_discriminating_evidence(txt, replicate)
    result = compare_txt_to_reference_profile(txt, _reference_from_run(run1, 0))
    updated = _add_discrimination_counts([result], evidence)[0]
    assert updated.discriminating_peak_support_count and updated.discriminating_peak_conflict_count
    assert "CONFLICTING_DISCRIMINATING_PEAKS" in updated.block_reasons


def test_partial_window_discrimination_counts_follow_run_label(tmp_path):
    txt = txt_peak_profile(tmp_path, [(102.0, 10)])
    run1 = with_scan_recurrence(gaussian_profile("UAA_T1_RUN_1", [(100.5, 10)]))
    run2 = with_scan_recurrence(gaussian_profile("UAA_T1_RUN_2", [(102.0, 10)]))
    matches = match_replicate_peaks(run1, run2)
    replicate = ReplicateConsistencyAuditResult(
        ReplicateAuditParameters(), (run1, run2), matches,
        (summarize_replicate_consistency(run1, run2, matches),),
    )
    evidence = build_discriminating_evidence(txt, replicate)
    window = SourceLinkageReferenceProfile(
        "PARTIAL_SCAN__UAA_T1_RUN_1__1", HypothesisType.PARTIAL_SCAN_EXPORT_POSSIBLE,
        ("UAA_T1_RUN_1",), run1, 0, 1, 10,
    )
    updated = _add_discrimination_counts(
        [compare_txt_to_reference_profile(txt, window)], evidence,
    )[0]
    assert updated.discriminating_peak_support_count == 0
    assert updated.discriminating_peak_conflict_count == 1
    assert "CONFLICTING_DISCRIMINATING_PEAKS" in updated.block_reasons


def test_common_negative_polarity_can_be_supported_without_exact_run(tmp_path):
    txt = txt_peak_profile(tmp_path, [(100.5, 10)])
    run1 = gaussian_profile("UAA_T1_RUN_1", [(100.5, 10)])
    run2 = gaussian_profile("UAA_T1_RUN_2", [(100.5, 10)])
    results = [compare_txt_to_reference_profile(txt, _reference_from_run(run1, 0)),
               compare_txt_to_reference_profile(txt, _reference_from_run(run2, 1))]
    summary = summarize_txt_mzml_source_linkage(results, txt_file="x",
        source_metadata_records=[metadata("R1"), metadata("R2")])
    assert not summary.exact_run_linkage_confirmed
    assert summary.common_source_polarity_supported and summary.polarity_propagation_eligible
    assert not summary.polarity_propagation_applied


def test_polarity_conflict_blocks_eligibility(tmp_path):
    txt = txt_peak_profile(tmp_path, [(100.5, 10)])
    run = gaussian_profile("UAA_T1_RUN_1", [(100.5, 10)])
    result = compare_txt_to_reference_profile(txt, _reference_from_run(run, 0))
    summary = summarize_txt_mzml_source_linkage([result], txt_file="x",
        source_metadata_records=[metadata("R1"), metadata("R2", PolarityStatus.POSITIVE_ONLY)])
    assert not summary.common_source_polarity_supported and not summary.polarity_propagation_eligible
    assert "SOURCE_POLARITY_CONFLICT" in summary.polarity_propagation_block_reasons


def test_determinism_formal_nonpropagation_and_optional_payload(tmp_path):
    txt = txt_peak_profile(tmp_path, [(100.5, 10)])
    run1 = gaussian_profile("UAA_T1_RUN_1", [(100.5, 10)])
    run2 = gaussian_profile("UAA_T1_RUN_2", [(101.5, 10)])
    references = [_reference_from_run(run1, 0), _reference_from_run(run2, 1)]
    results = [compare_txt_to_reference_profile(txt, reference) for reference in references]
    one = summarize_txt_mzml_source_linkage(results, txt_file="x", source_metadata_records=[metadata("R1"), metadata("R2")])
    two = summarize_txt_mzml_source_linkage(reversed(results), txt_file="x", source_metadata_records=[metadata("R2"), metadata("R1")])
    assert one == two and not one.formal_propagation and not one.polarity_propagation_applied
    from rna_masshunter.sciex_t1_txt_mzml_source_linkage_audit import TxtMzMLSourceLinkageAuditResult
    audit = TxtMzMLSourceLinkageAuditResult(SourceLinkageParameters(), txt.txt_profile, txt, tuple(references), tuple(results), (), (), one)
    payload = audit_optional_result(audit)
    assert set(payload) == {"txt_metadata_records", "hypothesis_records", "peak_evidence_records", "discrimination_records", "summary_records"}
    safeguards = (
        "shadow_analysis_only", "source_linkage_audit_only", "formal_propagation",
        "polarity_propagation_applied", "chemical_identity_assigned",
        "fragment_identity_assigned", "charge_state_confirmed",
    )
    for records in payload.values():
        for row in records:
            assert all(key in row for key in safeguards)
            assert row["shadow_analysis_only"] and row["source_linkage_audit_only"]
            assert all(not row[key] for key in safeguards[2:])
