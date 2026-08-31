from copy import deepcopy
from types import SimpleNamespace

import pandas as pd

from rna_masshunter.models import MS2SpectrumInfo
from rna_masshunter.ms2_annotation import _scan_window
from rna_masshunter.ms2_identity_evidence import build_ms2_modification_identity, physical_observed_peak_key
from rna_masshunter.ms2_unmatched_audit import (
    AUDIT_COLUMNS,
    build_unmatched_ion_audit,
    build_unmatched_ion_summary,
)
from rna_masshunter.review_dashboard import _build_top_candidates


def config(enabled=True):
    return SimpleNamespace(ms2_annotation={"enabled": enabled, "mz_tolerance_ppm": 20})


def spectrum(raw_peaks=None, selected=None, scan_min=50.0, scan_max=150.0, threshold=10.0, threshold_available=True):
    raw = [(100.0, 100.0)] if raw_peaks is None else raw_peaks
    chosen = list(raw if selected is None else selected)
    return MS2SpectrumInfo(
        spectrum_id="S1", scan_index=1, rt=1.25, precursor_mz=500.0,
        precursor_charge=2, precursor_intensity=None, num_peaks=len(chosen),
        base_peak_mz=chosen[0][0] if chosen else None,
        base_peak_intensity=max((x[1] for x in raw), default=None), total_ion_current=sum(x[1] for x in raw),
        peaks=chosen, raw_peaks=raw, scan_mz_min=scan_min, scan_mz_max=scan_max,
        effective_intensity_threshold=threshold, threshold_information_available=threshold_available,
    )


def ion(mz=100.0, ion_id="I1", modification="cnm5U", position=2, charge=1):
    return {
        "Ion_ID": ion_id, "Spectrum_ID": "S1", "Parent_Fragment_ID": "F1",
        "Parent_Start": 35, "Parent_Sequence": "UUU", "Modification_ID": modification,
        "Modification_Name": modification, "Candidate_Modification_Position_In_Parent": position,
        "Candidate_Modification_Base": "U", "Ion_Type": "d", "Ion_Start": 1,
        "Ion_End": 2, "Ion_Length": 2, "Charge": charge, "Theoretical_mz": mz,
        "Ion_Contains_Modification": True,
    }


def match(theoretical=None, modification="cnm5U", position=2, ion_id="I1"):
    theoretical = theoretical or 100.0
    return {
        "Spectrum_ID": "S1", "Scan_Index": 1, "Observed_mz": theoretical,
        "Observed_Intensity": 100.0, "Parent_Fragment_ID": "F1",
        "Modification_ID": modification, "Candidate_Modification_Position_In_Parent": position,
        "Ion_ID": ion_id, "Theoretical_mz": theoretical, "Mass_Error_Da": 0.0,
        "Mass_Error_ppm": 0.0, "Ion_Contains_Modification": True,
    }


def audit_status(spec, current_ion=None, matches=None):
    rows, diagnostics = build_unmatched_ion_audit([spec] if spec is not None else [], [current_ion or ion()], matches or [], config())
    assert list(rows[0]) == AUDIT_COLUMNS
    assert diagnostics[0]["Apply_Unmatched_Audit_To_Final_Score"] is False
    return rows[0]


def test_existing_match_is_matched_without_reevaluation():
    row = audit_status(spectrum(), matches=[match()])
    assert row["Unmatched_Reason_Status"] == "matched"
    assert row["Is_Matched"] is True
    assert row["Existing_Match_IDs"]


def test_theoretical_mz_outside_explicit_scan_range():
    row = audit_status(spectrum(), ion(mz=200.0))
    assert row["Unmatched_Reason_Status"] == "outside_scan_mz_range"


def test_no_raw_observed_peak_in_audit_window():
    row = audit_status(spectrum(raw_peaks=[(101.0, 100.0)]))
    assert row["Unmatched_Reason_Status"] == "no_observed_peak_in_search_window"


def test_nearest_peak_outside_formal_tolerance():
    row = audit_status(spectrum(raw_peaks=[(100.003, 100.0)]))
    assert row["Unmatched_Reason_Status"] == "nearest_peak_outside_tolerance"
    assert row["Nearest_Error_ppm"] > 20


def test_raw_peak_below_effective_intensity_threshold():
    row = audit_status(spectrum(raw_peaks=[(100.001, 5.0)], selected=[], threshold=10.0))
    assert row["Unmatched_Reason_Status"] == "peak_below_intensity_threshold"
    assert row["Below_Threshold"] is True


def test_threshold_unknown_does_not_assert_below_threshold():
    row = audit_status(spectrum(raw_peaks=[(100.001, 5.0)], threshold=None, threshold_available=False))
    assert row["Unmatched_Reason_Status"] == "threshold_information_unavailable"
    assert row["Below_Threshold"] == ""


def test_scan_range_unknown_does_not_assert_outside_range():
    row = audit_status(spectrum(raw_peaks=[(101.0, 20.0)], scan_min=None, scan_max=None))
    assert row["Unmatched_Reason_Status"] == "scan_range_not_available"


def test_spectrum_not_available():
    row = audit_status(None)
    assert row["Unmatched_Reason_Status"] == "spectrum_not_available"


def test_multiple_nearby_raw_peaks_are_ambiguous():
    row = audit_status(spectrum(raw_peaks=[(99.999, 20.0), (100.001, 30.0)]))
    assert row["Unmatched_Reason_Status"] == "ambiguous_multiple_nearby_peaks"
    assert row["Nearby_Peak_Count"] == 2


def test_physical_peak_key_logic_is_reused_for_nearest_peak():
    spec = spectrum(raw_peaks=[(100.003, 25.0)])
    row = audit_status(spec)
    expected = physical_observed_peak_key({"Spectrum_ID": "S1", "Observed_mz": 100.003, "RT": 1.25})
    assert row["Nearest_Physical_Observed_Peak_Key"] == expected


def test_charge_specific_theoretical_mz_is_audited_independently():
    ions = [ion(mz=100.0, ion_id="I1", charge=1), ion(mz=50.0, ion_id="I2", charge=2)]
    rows, _ = build_unmatched_ion_audit([spectrum(raw_peaks=[(100.0, 30.0), (50.003, 30.0)])], ions, [], config())
    assert len(rows) == 2
    assert {row["Ion_Charge"] for row in rows} == {1, 2}


def test_duplicate_semantic_theoretical_ion_is_audited_once():
    duplicate = dict(ion()); duplicate["Ion_ID"] = "DUPLICATE_ID"
    rows, _ = build_unmatched_ion_audit([spectrum()], [ion(), duplicate], [], config())
    assert len(rows) == 1


def test_cross_candidate_assignment_is_not_an_exact_match():
    row = audit_status(spectrum(), matches=[match(modification="other")])
    assert row["Is_Matched"] is False
    assert row["Unmatched_Reason_Status"] != "matched"


def test_zero_theoretical_ions_and_summary():
    rows, diagnostics = build_unmatched_ion_audit([spectrum()], [], [], config())
    assert rows == []
    assert diagnostics[0]["Total_Modified_Theoretical_Ions"] == 0
    summary = build_unmatched_ion_summary([{"Rank": 1, "Modification_ID": "x", "Parent_Fragment_ID": "F", "Candidate_tRNA_Position": 1, "Candidate_Position_In_Parent": 1}], rows)
    assert summary[0]["Total_Modified_Theoretical_Ion_Count"] == 0
    assert summary[0]["Recommended_Followup"] == "insufficient_information"


def test_ms2_disabled_produces_disabled_empty_audit():
    rows, diagnostics = build_unmatched_ion_audit([spectrum()], [ion()], [], config(False), enabled=False)
    assert rows == []
    assert diagnostics[0]["MS2_Unmatched_Ion_Audit_Enabled"] is False


def test_raw_peaks_unavailable_is_information_unavailable():
    spec = spectrum(); spec.raw_peaks = None
    row = audit_status(spec)
    assert row["Unmatched_Reason_Status"] == "insufficient_information"


def test_threshold_passing_peak_removed_by_traced_filter():
    row = audit_status(spectrum(raw_peaks=[(100.001, 20.0)], selected=[], threshold=10.0))
    assert row["Unmatched_Reason_Status"] == "peak_present_but_filtered"


def test_scan_window_comes_from_explicit_mzml_metadata_only():
    raw = {"scanList": {"scan": [{"scanWindowList": {"scanWindow": [{"scan window lower limit": 70.0, "scan window upper limit": 2000.0}]}}]}}
    assert _scan_window(raw) == (70.0, 2000.0)
    assert _scan_window({"m/z array": [100.0, 200.0]}) == (None, None)


def test_cnm5u_canonical_prior_cannot_change_reason_classification():
    spec = spectrum(raw_peaks=[(100.003, 100.0)])
    plain = audit_status(spec)
    cfg = config(); cfg.ms2_annotation["biological_position_prior"] = {"canonical_positions": {"cnm5U": 37}}
    rows, _ = build_unmatched_ion_audit([spec], [ion()], [], cfg)
    assert rows[0]["Unmatched_Reason_Status"] == plain["Unmatched_Reason_Status"]


def test_audit_does_not_mutate_existing_matches_localization_or_identity_inputs():
    matches = [match()]
    localization = [{"Modification_ID": "cnm5U", "Parent_Fragment_ID": "F1", "Candidate_Modification_Position_In_Parent": 2}]
    ranking = [{"Rank": 1, "Final_Score": 7.0, "Final_Confidence": "Medium", "Modification_ID": "cnm5U", "Parent_Fragment_ID": "F1", "Candidate_Position_In_Parent": 2, "Candidate_tRNA_Position": 36}]
    before = deepcopy((matches, localization, ranking))
    build_unmatched_ion_audit([spectrum()], [ion()], matches, config())
    assert (matches, localization, ranking) == before


def test_audit_summary_does_not_change_identity_biological_or_formal_fields():
    ranking = [{"Rank": 1, "Final_Score": 7.0, "Final_Confidence": "Medium", "Modification_ID": "cnm5U", "Modification_Name": "cnm5U", "Parent_Fragment_ID": "F1", "Parent_Sequence": "UUU", "Candidate_Position_In_Parent": 2, "Candidate_tRNA_Position": 36, "Position_Class": "canonical_position", "Biological_Plausibility_Level": "high"}]
    original = deepcopy(ranking)
    audit_rows, _ = build_unmatched_ion_audit([spectrum(raw_peaks=[(101.0, 20.0)])], [ion()], [], config())
    summary = build_unmatched_ion_summary(ranking, audit_rows)
    assert ranking == original
    assert summary[0]["No_Peak_In_Window_Count"] == 1


def test_top_50_order_priority_and_formal_columns_are_unchanged():
    ranking = pd.DataFrame([
        {"Rank": i, "Final_Score": 100 - i, "Final_Confidence": "Medium", "Modification_ID": f"m{i}", "Modification_Name": f"m{i}", "Parent_Fragment_ID": "F", "Parent_Sequence": "UU", "Candidate_Position_In_Parent": i, "Candidate_tRNA_Position": i}
        for i in range(1, 61)
    ])
    plain = _build_top_candidates(ranking, pd.DataFrame(), {"max_top_candidates": 50})
    summaries = pd.DataFrame([{"Modification_ID": f"m{i}", "Parent_Fragment_ID": "F", "Candidate_tRNA_Position": i, "Total_Modified_Theoretical_Ion_Count": 1, "Unmatched_Modified_Theoretical_Ion_Count": 1} for i in range(1, 61)])
    audited = _build_top_candidates(ranking, pd.DataFrame(), {"max_top_candidates": 50}, summaries)
    formal = ["Review_Rank", "Review_Priority", "Modification_ID", "Best_Final_Score", "Best_Final_Confidence"]
    assert len(plain) == len(audited) == 50
    assert plain[formal].equals(audited[formal])


def test_all_new_sheet_names_fit_excel_limit():
    assert all(len(name) <= 31 for name in ["MS2_Unmatched_Ion_Audit", "MS2_Unmatched_Ion_Summary", "MS2_Unmatched_Ion_Diagnostics"])
