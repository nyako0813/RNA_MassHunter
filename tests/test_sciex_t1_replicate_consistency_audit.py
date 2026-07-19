from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from rna_masshunter.sciex_mzml_source_metadata_audit import (
    MzMLSourceMetadataRecord,
    PolarityStatus,
    ReadStatus,
    RepresentationStatus,
)
from rna_masshunter.sciex_t1_replicate_consistency_audit import (
    MatchAmbiguityStatus,
    OverallReplicateStatus,
    ReplicateAuditParameters,
    ReplicateConsistencyStatus,
    ReplicateRunPeak,
    ReplicateRunPeakProfile,
    SystematicMZShiftStatus,
    _validate_pair_context,
    audit_optional_result,
    build_ms1_peak_profile_for_run,
    build_ms1_peak_profile_from_spectra,
    match_replicate_peaks,
    summarize_replicate_consistency,
)


def metadata(label="RUN_1", *, rna="TRNA_LEU_UAA", digest="T1_DIGEST",
             polarity=PolarityStatus.NEGATIVE_ONLY, representation=RepresentationStatus.PROFILE_ONLY):
    return MzMLSourceMetadataRecord(
        input_path=f"{label}.mzML", file_name=f"{label}.mzML", read_status=ReadStatus.COMPLETED,
        technical_run_label=label, context_source="USER_PROVIDED_RUNTIME_MANIFEST",
        context_confidence="USER_CONFIRMED", rna_identity=rna, digest_type=digest,
        polarity_status=polarity, representation_status=representation,
    )


def spectrum(level, peaks=(), *, scale=1.0, rt=1.0, arrays=None):
    if arrays is None:
        x = np.arange(99.8, 101.21, 0.005)
        y = np.zeros_like(x)
        for center, height, width in peaks:
            y += height * scale * np.exp(-0.5 * ((x - center) / width) ** 2)
    else:
        x, y = arrays
    return {
        "ms level": level, "m/z array": x, "intensity array": y,
        "scanList": {"scan": [{"scan start time": rt, "unitName": "minute"}]},
    }


def run_profile(label, peaks):
    return ReplicateRunPeakProfile(
        run_label=label, input_path=f"{label}.mzML", status="COMPLETED",
        aggregation_method="TEST", ms1_spectra_total=10, ms1_spectra_used=10,
        ms1_spectra_excluded=0, ms2_spectra_excluded=0, missing_ms_level_spectra=0,
        mz_grid_method="TEST", intensity_normalization_method="TEST", baseline_method="TEST",
        smoothing_method="TEST", peak_detection_method="TEST", detected_peak_count=len(peaks),
        peaks=tuple(peaks), polarity_status="NEGATIVE_ONLY", representation_status="PROFILE_ONLY",
        block_reasons=(),
    )


def peak(identifier, mz, *, intensity=1.0, rank=1, recurrence=0.1, fwhm=0.02,
         prominence=0.1, quality="MAJOR_SHARP"):
    return ReplicateRunPeak(
        run_label=identifier.split("_")[0], peak_id=identifier, apex_mz=mz, centroid_mz=mz,
        raw_apex_intensity=intensity * 100, normalized_apex_intensity=intensity,
        raw_integrated_intensity=intensity * 10, normalized_integrated_intensity=intensity / 2,
        relative_intensity=intensity, intensity_rank=rank, prominence=prominence,
        relative_prominence=prominence, fwhm=fwhm, left_bound_mz=mz - 0.02,
        right_bound_mz=mz + 0.02, supporting_ms1_scan_count=round(recurrence * 10),
        total_ms1_scan_count=10, scan_recurrence_fraction=recurrence,
        first_supporting_scan_time=1.0, last_supporting_scan_time=2.0,
        detection_status=quality, detection_block_reasons=(),
    )


class TrapArray:
    def decode(self):
        raise AssertionError("MS2 binary array must not be decoded")


def test_ms1_only_extraction_and_ms2_binary_not_decoded():
    spectra = [spectrum(1, [(100.0, 10, 0.01)], rt=1),
               {"ms level": 2, "m/z array": TrapArray(), "intensity array": TrapArray()},
               spectrum(1, [(100.0, 8, 0.01)], rt=2)]
    result = build_ms1_peak_profile_from_spectra(spectra, run_label="RUN_1", metadata_record=metadata())
    assert result.ms1_spectra_total == result.ms1_spectra_used == 2
    assert result.ms2_spectra_excluded == 1
    assert "NON_MS1_SPECTRUM_EXCLUDED" in result.block_reasons


def test_independent_run_processing_and_intensity_normalization():
    one = build_ms1_peak_profile_from_spectra(
        [spectrum(1, [(100.0, 10, 0.01)], scale=1) for _ in range(3)],
        run_label="RUN_1", metadata_record=metadata("RUN_1"))
    two = build_ms1_peak_profile_from_spectra(
        [spectrum(1, [(100.0, 10, 0.01)], scale=10) for _ in range(3)],
        run_label="RUN_2", metadata_record=metadata("RUN_2"))
    assert one.run_label != two.run_label and one.ms1_spectra_used == two.ms1_spectra_used == 3
    assert one.detected_peak_count and two.detected_peak_count
    assert one.peaks[0].normalized_apex_intensity == pytest.approx(two.peaks[0].normalized_apex_intensity)
    assert two.peaks[0].raw_apex_intensity == pytest.approx(one.peaks[0].raw_apex_intensity * 10)


def test_scan_recurrence_distinguishes_multi_and_single_scan_peaks():
    scans = [spectrum(1, [(100.0, 10, 0.01), (101.0, 7 if index == 0 else 0, 0.01)], rt=index) for index in range(5)]
    result = build_ms1_peak_profile_from_spectra(scans, run_label="RUN_1", metadata_record=metadata())
    common = min(result.peaks, key=lambda item: abs(item.apex_mz - 100.0))
    single = min(result.peaks, key=lambda item: abs(item.apex_mz - 101.0))
    assert common.supporting_ms1_scan_count == 5 and common.scan_recurrence_fraction == 1
    assert single.supporting_ms1_scan_count == 1 and single.scan_recurrence_fraction == pytest.approx(0.2)


def test_exact_mz_one_to_one_match():
    matches = match_replicate_peaks(run_profile("A", [peak("A_1", 100)]), run_profile("B", [peak("B_1", 100)]))
    assert len(matches) == 1 and matches[0].absolute_delta_mz == 0
    assert matches[0].match_ambiguity_status is MatchAmbiguityStatus.UNAMBIGUOUS_ONE_TO_ONE


@pytest.mark.parametrize("delta,matched", [(0.01, True), (0.01001, False)])
def test_absolute_tolerance_edge(delta, matched):
    result = match_replicate_peaks(
        run_profile("A", [peak("A_1", 100)]), run_profile("B", [peak("B_1", 100 + delta)]),
        absolute_tolerance_da=0.01, ppm_tolerance=0.0001)
    assert any(item.run_1_peak_id and item.run_2_peak_id for item in result) is matched


def test_ppm_tolerance_at_high_mz():
    result = match_replicate_peaks(
        run_profile("A", [peak("A_1", 2000)]), run_profile("B", [peak("B_1", 2000.015)]),
        absolute_tolerance_da=0.001, ppm_tolerance=10)
    assert sum(bool(item.run_1_peak_id and item.run_2_peak_id) for item in result) == 1


def test_ambiguous_matching_and_alternatives():
    result = match_replicate_peaks(
        run_profile("A", [peak("A_1", 100)]),
        run_profile("B", [peak("B_1", 99.995), peak("B_2", 100.005, rank=2)]),
        absolute_tolerance_da=0.01, ppm_tolerance=0.0001)
    matched = next(item for item in result if item.run_1_peak_id and item.run_2_peak_id)
    assert matched.match_ambiguity_status is MatchAmbiguityStatus.AMBIGUOUS_LEFT_TO_MULTIPLE
    assert matched.alternative_match_peak_ids
    assert matched.replicate_consistency_status is ReplicateConsistencyStatus.REPRODUCED_AMBIGUOUS_MATCH


def test_unmatched_peaks_are_explicit():
    result = match_replicate_peaks(
        run_profile("A", [peak("A_1", 100), peak("A_2", 200, rank=2)]),
        run_profile("B", [peak("B_1", 100), peak("B_2", 300, rank=2)]))
    statuses = {item.replicate_consistency_status for item in result}
    assert ReplicateConsistencyStatus.RUN_1_ONLY in statuses
    assert ReplicateConsistencyStatus.RUN_2_ONLY in statuses


def test_fwhm_missing_is_blocked_without_crash():
    result = match_replicate_peaks(
        run_profile("A", [peak("A_1", 100, fwhm=None)]),
        run_profile("B", [peak("B_1", 100)]))[0]
    assert result.replicate_consistency_status is ReplicateConsistencyStatus.REPRODUCED_WITH_SHAPE_VARIATION
    assert "MISSING_PEAK_SHAPE_METRICS" in result.consistency_block_reasons


def test_prominence_missing_is_blocked_without_crash():
    left_peak = replace(peak("A_1", 100), prominence=None, relative_prominence=None)
    result = match_replicate_peaks(run_profile("A", [left_peak]), run_profile("B", [peak("B_1", 100)]))[0]
    assert result.replicate_consistency_status is ReplicateConsistencyStatus.REPRODUCED_WITH_SHAPE_VARIATION
    assert "MISSING_PEAK_SHAPE_METRICS" in result.consistency_block_reasons


def test_input_file_not_found_is_safe(tmp_path):
    result = build_ms1_peak_profile_for_run(tmp_path / "missing.mzML", metadata_record=metadata())
    assert result.status == "BLOCKED" and result.block_reasons == ("INPUT_FILE_NOT_FOUND",)


def test_audit_does_not_mutate_existing_formal_sentinel():
    sentinel = {"formal_score": 7, "rank": 2, "evidence_tier": "B", "final_consensus": "UNCHANGED"}
    before = dict(sentinel)
    a = run_profile("A", [peak("A_1", 100)])
    b = run_profile("B", [peak("B_1", 100)])
    summarize_replicate_consistency(a, b, match_replicate_peaks(a, b))
    assert sentinel == before


def test_low_quality_match_is_preserved():
    result = match_replicate_peaks(
        run_profile("A", [peak("A_1", 100, quality="LOW_SUPPORT")]),
        run_profile("B", [peak("B_1", 100)]))[0]
    assert result.replicate_consistency_status is ReplicateConsistencyStatus.INSUFFICIENT_PEAK_QUALITY


def test_systematic_mz_shift_summary():
    left = run_profile("A", [peak(f"A_{i}", mz, rank=i) for i, mz in enumerate((100, 200, 300), 1)])
    right = run_profile("B", [peak(f"B_{i}", mz + 0.003, rank=i) for i, mz in enumerate((100, 200, 300), 1)])
    matches = match_replicate_peaks(left, right)
    summary = summarize_replicate_consistency(left, right, matches)
    assert summary.drift.median_delta_mz == pytest.approx(0.003)
    assert summary.drift.mad_delta_mz == pytest.approx(0)
    assert summary.drift.systematic_mz_shift_status is SystematicMZShiftStatus.SMALL_SYSTEMATIC_SHIFT
    assert not summary.drift.drift_adjustment_applied


def test_deterministic_canonical_pair_and_output_order():
    a = run_profile("A", [peak("A_2", 200, rank=2), peak("A_1", 100)])
    b = run_profile("B", [peak("B_2", 200, rank=2), peak("B_1", 100)])
    assert match_replicate_peaks(a, b) == match_replicate_peaks(b, a)


def test_formal_non_propagation_and_optional_rows_are_scalars():
    a = run_profile("A", [peak("A_1", 100)])
    b = run_profile("B", [peak("B_1", 100)])
    matches = match_replicate_peaks(a, b)
    summary = summarize_replicate_consistency(a, b, matches)
    from rna_masshunter.sciex_t1_replicate_consistency_audit import ReplicateConsistencyAuditResult
    result = ReplicateConsistencyAuditResult(ReplicateAuditParameters(), (a, b), matches, (summary,))
    payload = audit_optional_result(result)
    assert set(payload) == {"run_peak_records", "match_records", "summary_records"}
    assert not result.formal_propagation and not summary.formal_propagation
    assert all(not row["formal_propagation"] for rows in payload.values() for row in rows)


def test_metadata_polarity_representation_and_context_mismatch():
    left = metadata("A")
    right = metadata("B", rna="OTHER", polarity=PolarityStatus.POSITIVE_ONLY,
                     representation=RepresentationStatus.CENTROID_ONLY)
    blocks = _validate_pair_context(left, right)
    assert blocks == (
        "REPLICATE_CONTEXT_MISMATCH", "POLARITY_MISMATCH_BETWEEN_REPLICATES",
        "REPRESENTATION_MISMATCH_BETWEEN_REPLICATES",
    )


def test_missing_ms_level_and_invalid_ms1_are_excluded():
    result = build_ms1_peak_profile_from_spectra(
        [{"m/z array": [], "intensity array": []}, spectrum(1, arrays=(np.array([1, 2]), np.array([1, 2])))],
        run_label="RUN_1", metadata_record=metadata())
    assert result.missing_ms_level_spectra == 1 and result.ms1_spectra_excluded == 1
    assert result.status == "BLOCKED" and "MISSING_MS_LEVEL_METADATA" in result.block_reasons


def test_no_peaks_safe_summary():
    empty = run_profile("A", [])
    nonempty = run_profile("B", [peak("B_1", 100)])
    matches = match_replicate_peaks(empty, nonempty)
    summary = summarize_replicate_consistency(empty, nonempty, matches)
    assert summary.matched_peak_pair_count == 0
    assert summary.replicate_consistency_overall_status is OverallReplicateStatus.INSUFFICIENT_COMPARABLE_PEAKS
    assert "NO_DETECTED_PEAKS_RUN_1" in summary.overall_block_reasons
