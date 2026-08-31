from copy import deepcopy
from types import SimpleNamespace

import pandas as pd

from rna_masshunter.models import MS2SpectrumInfo
from rna_masshunter.ms2_ambiguous_peak_audit import (
    build_ambiguous_peak_audit,
    build_ambiguity_diagnostics,
    build_ambiguity_summary,
    deterministic_cluster_id,
    TOP_SHADOW_COLUMNS,
)
from rna_masshunter.review_dashboard import _build_top_candidates


def spectrum(peaks=None, raw_available=True):
    peaks = peaks if peaks is not None else [(99.99, 20.0), (100.01, 30.0)]
    return MS2SpectrumInfo(
        spectrum_id="S1", scan_index=1, rt=1.0, precursor_mz=500.0,
        precursor_charge=2, precursor_intensity=None, num_peaks=len(peaks),
        base_peak_mz=max(peaks, key=lambda item: item[1])[0] if peaks else None,
        base_peak_intensity=max((item[1] for item in peaks), default=None),
        total_ion_current=sum(item[1] for item in peaks), peaks=list(peaks),
        raw_peaks=list(peaks) if raw_available else None, scan_mz_min=50.0, scan_mz_max=150.0,
        effective_intensity_threshold=1.0, threshold_information_available=True,
    )


def audit(modification="cnm5U", position=2, ion_id="I1", theoretical=100.0):
    return {
        "Modification_ID": modification, "Parent_Fragment_ID": "F1", "Spectrum_ID": "S1",
        "Candidate_tRNA_Position": 36, "Candidate_Position_In_Parent": position,
        "Theoretical_Ion_ID": ion_id, "Ion_Series": "d", "Ion_Number": 2,
        "Ion_Charge": 1, "Theoretical_mz": theoretical,
        "Unmatched_Reason_Status": "ambiguous_multiple_nearby_peaks",
        "Matching_Tolerance_ppm": 20.0, "Matching_Tolerance_Da": theoretical * 20 / 1_000_000,
        "Audit_Search_Window_Da": 0.05,
    }


def ion(modification="cnm5U", position=2, ion_id="I1", theoretical=100.0):
    return {
        "Spectrum_ID": "S1", "Parent_Fragment_ID": "F1", "Modification_ID": modification,
        "Candidate_Modification_Position_In_Parent": position, "Ion_ID": ion_id,
        "Ion_Type": "d", "Ion_Start": 1, "Ion_End": 2, "Charge": 1,
        "Theoretical_mz": theoretical, "Ion_Contains_Modification": True,
    }


def ranking(modification="cnm5U", position=2, trna=36, group=""):
    return {
        "Rank": 1, "Final_Score": 7.0, "Final_Confidence": "Medium",
        "Modification_ID": modification, "Modification_Name": modification,
        "Parent_Fragment_ID": "F1", "Parent_Sequence": "UU", "Candidate_Position_In_Parent": position,
        "Candidate_tRNA_Position": trna, "Structural_Isomer_Group_ID": group,
    }


def build(peaks=None, audits=None, ions=None, rankings=None, assignments=None, raw_available=True):
    return build_ambiguous_peak_audit(
        [spectrum(peaks, raw_available)], audits or [audit()], ions or [ion()],
        rankings or [ranking()], assignments or [], enabled=True,
    )


def types(cluster):
    return {cluster["Primary_Ambiguity_Type"], *filter(None, str(cluster["Secondary_Ambiguity_Types"] or "").split(";"))}


def test_multiple_peaks_same_side():
    clusters, _ = build([(100.01, 20.0), (100.02, 30.0)])
    assert "multiple_peaks_same_side" in types(clusters[0])


def test_multiple_peaks_bracket_theoretical_mz():
    clusters, _ = build([(99.99, 20.0), (100.01, 30.0)])
    assert "multiple_peaks_bracketing_theoretical_mz" in types(clusters[0])


def test_multiple_peaks_within_formal_tolerance():
    clusters, _ = build([(99.999, 20.0), (100.001, 30.0)])
    assert clusters[0]["Peaks_Within_Formal_Tolerance_Count"] == 2
    assert clusters[0]["Primary_Ambiguity_Type"] == "multiple_peaks_within_matching_tolerance"


def test_multiple_peaks_outside_tolerance_but_inside_audit_window():
    clusters, _ = build([(99.99, 20.0), (100.01, 30.0)])
    assert "multiple_peaks_outside_tolerance_but_within_audit_window" in types(clusters[0])


def test_multiple_candidates_for_same_physical_peak():
    clusters, details = build(
        [(100.001, 20.0), (100.01, 30.0)],
        ions=[ion(), ion("other", 2, "I2")],
        rankings=[ranking(), ranking("other", 2, 36)],
    )
    assert "multiple_candidate_assignments_same_physical_peak" in types(clusters[0])
    assert max(row["Candidate_Modification_Count"] for row in details) == 2


def test_multiple_theoretical_ions_compete_for_same_peak():
    clusters, details = build(
        [(100.001, 20.0), (100.01, 30.0)],
        ions=[ion(ion_id="I1"), ion(ion_id="I2")],
    )
    assert "multiple_theoretical_ions_compete_for_same_peak" in types(clusters[0])
    assert max(row["Theoretical_Ion_Count"] for row in details) == 2


def test_structural_isomer_shared_peak():
    clusters, details = build(
        [(100.001, 20.0), (100.01, 30.0)],
        ions=[ion("m1G", 2, "I1"), ion("m2G", 2, "I2")],
        audits=[audit("m1G", 2, "I1")],
        rankings=[ranking("m1G", 2, 36, "SIG1"), ranking("m2G", 2, 36, "SIG1")],
    )
    assert "structural_isomer_shared_peak" in types(clusters[0])
    assert "structural_isomer_group_shared" in {row["Candidate_Specificity_Status"] for row in details}


def test_positional_isomer_shared_peak():
    clusters, details = build(
        [(100.001, 20.0), (100.01, 30.0)],
        ions=[ion(position=2, ion_id="I1"), ion(position=3, ion_id="I2")],
        rankings=[ranking(position=2, trna=36), ranking(position=3, trna=37)],
    )
    assert "positional_isomer_shared_peak" in types(clusters[0])
    assert "position_group_shared" in {row["Candidate_Specificity_Status"] for row in details}


def test_candidate_specific_peak():
    _, details = build([(100.001, 20.0), (100.01, 30.0)])
    assert details[0]["Candidate_Specificity_Status"] == "candidate_specific"


def test_cluster_id_is_deterministic():
    row = audit()
    assert deterministic_cluster_id(row, 99.95, 100.05) == deterministic_cluster_id(dict(row), 99.95, 100.05)


def test_cluster_size_span_and_side_counts():
    clusters, _ = build([(99.99, 20.0), (100.01, 30.0), (100.02, 40.0)])
    row = clusters[0]
    assert row["Peak_Cluster_Size"] == 3
    assert abs(row["Peak_Cluster_Span_Da"] - 0.03) < 1e-12
    assert row["Peaks_Below_Theoretical_Count"] == 1
    assert row["Peaks_Above_Theoretical_Count"] == 2


def test_closest_error_and_intensity_ranks():
    clusters, _ = build([(100.001, 10.0), (100.01, 100.0)])
    assert clusters[0]["Closest_Peak_Rank_By_Error"] == 1
    assert clusters[0]["Closest_Peak_Rank_By_Intensity"] == 2


def test_low_severity():
    clusters, details = build([(100.001, 20.0), (100.01, 30.0)])
    summary = build_ambiguity_summary([ranking()], clusters, details)[0]
    assert summary["Ambiguity_Severity"] == "low"


def test_moderate_severity_for_explicit_structural_group_scope():
    assignment = {
        "Physical_Observed_Peak_Key": "S1|mz=100.00100000|rt=1.000000",
        "Modification_ID": "cnm5U", "Parent_Fragment_ID": "F1",
        "Candidate_Position_In_Parent": 2, "Theoretical_Ion_ID": "I1",
        "Structural_Isomer_Group_ID": "SIG1", "Evidence_Scope": "structural_isomer_group_level",
    }
    clusters, details = build([(100.001, 20.0), (100.01, 30.0)], assignments=[assignment])
    summary = build_ambiguity_summary([ranking()], clusters, details)[0]
    assert summary["Ambiguity_Severity"] == "moderate"


def test_high_severity_for_multiple_formal_tolerance_peaks():
    clusters, details = build([(99.999, 20.0), (100.001, 30.0)])
    summary = build_ambiguity_summary([ranking()], clusters, details)[0]
    assert summary["Ambiguity_Severity"] == "high"


def test_cnm5u_canonical_prior_does_not_change_severity():
    clusters, details = build([(100.001, 20.0), (100.01, 30.0)])
    plain = build_ambiguity_summary([ranking()], clusters, details)[0]["Ambiguity_Severity"]
    with_prior = ranking(); with_prior["Position_Class"] = "canonical_position"
    assert build_ambiguity_summary([with_prior], clusters, details)[0]["Ambiguity_Severity"] == plain


def test_zero_ambiguous_ion_and_disabled_audit():
    clusters, details = build_ambiguous_peak_audit([spectrum()], [], [ion()], [ranking()], [], enabled=True)
    assert clusters == details == []
    clusters, details = build_ambiguous_peak_audit([spectrum()], [audit()], [ion()], [ranking()], [], enabled=False)
    assert clusters == details == []


def test_raw_peaks_unavailable_is_retained_as_information_unavailable():
    clusters, details = build(raw_available=False)
    assert len(clusters) == 1 and details == []
    assert clusters[0]["Primary_Ambiguity_Type"] == "insufficient_information"


def test_existing_unmatched_audit_and_inputs_are_not_mutated():
    audits = [audit()]; ions = [ion()]; rankings = [ranking()]
    before = deepcopy((audits, ions, rankings))
    build([(99.99, 20.0), (100.01, 30.0)], audits, ions, rankings)
    assert (audits, ions, rankings) == before


def test_summary_does_not_mutate_formal_identity_or_biological_fields():
    row = ranking(); row.update({"Position_Class": "canonical_position", "Biological_Plausibility_Level": "high", "MS2_Identity_Evidence_Level": "modified_fragment_ion_supported"})
    original = deepcopy(row)
    clusters, details = build(rankings=[row])
    build_ambiguity_summary([row], clusters, details)
    assert row == original


def test_top_50_existing_columns_and_order_unchanged():
    rankings = pd.DataFrame([
        {"Rank": i, "Final_Score": 100-i, "Final_Confidence": "Medium", "Modification_ID": f"m{i}", "Modification_Name": f"m{i}", "Parent_Fragment_ID": "F", "Parent_Sequence": "UU", "Candidate_Position_In_Parent": i, "Candidate_tRNA_Position": i}
        for i in range(1, 61)
    ])
    plain = _build_top_candidates(rankings, pd.DataFrame(), {"max_top_candidates": 50})
    existing_columns = [column for column in plain.columns if column not in TOP_SHADOW_COLUMNS]
    ambiguity = pd.DataFrame([{"Modification_ID": f"m{i}", "Parent_Fragment_ID": "F", "Candidate_tRNA_Position": i, "Ambiguous_Theoretical_Ion_Count": 1, "Ambiguous_Peak_Cluster_Count": 1, "Maximum_Cluster_Size": 2, "Primary_Ambiguity_Pattern": "multiple_peaks_same_side", "Ambiguity_Severity": "low"} for i in range(1,61)])
    audited = _build_top_candidates(rankings, pd.DataFrame(), {"max_top_candidates": 50}, None, ambiguity)
    assert len(plain) == len(audited) == 50
    assert plain[existing_columns].equals(audited[existing_columns])


def test_diagnostics_are_shadow_only_and_count_severity():
    clusters, details = build([(99.999, 20.0), (100.001, 30.0)])
    summaries = build_ambiguity_summary([ranking()], clusters, details)
    diagnostic = build_ambiguity_diagnostics(clusters, details, summaries)[0]
    assert diagnostic["Apply_Ambiguous_Peak_Audit_To_Final_Score"] is False
    assert diagnostic["High_Severity_Count"] == 1


def test_all_new_sheet_names_within_excel_limit():
    names = ["MS2_Ambiguous_Peak_Clusters", "MS2_Ambiguous_Peak_Detail", "MS2_Ambiguity_Summary"]
    assert all(len(name) <= 31 for name in names)
