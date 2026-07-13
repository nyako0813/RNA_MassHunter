from copy import deepcopy
from math import nan
from types import SimpleNamespace

from rna_masshunter.ms2_identity_evidence import physical_observed_peak_key
from rna_masshunter.ms2_zero_intensity_audit import (
    DETAIL_COLUMNS,
    DIAGNOSTIC_COLUMNS,
    SPECTRA_COLUMNS,
    SUMMARY_COLUMNS,
    TOP_SHADOW_COLUMNS,
    _origin_for_values,
    build_zero_intensity_audit,
    capture_source_spectrum,
    intensity_state,
    record_parsed_spectrum,
)
from rna_masshunter.review_dashboard import TOP_CANDIDATE_COLUMNS, _build_top_candidates


def _key(mz):
    return physical_observed_peak_key({"Spectrum_ID": "S1", "Observed_mz": mz, "RT": 1.0})


def _context(intensities, parsed=None, annotation_indices=None, mz=None):
    mz = mz or [100.0 + index for index in range(len(intensities))]
    source = capture_source_spectrum(
        {"id": "S1", "ms level": 2, "m/z array": mz, "intensity array": intensities, "centroid spectrum": True}, 1,
    )
    if parsed is not None:
        record_parsed_spectrum(source, mz, parsed, range(len(parsed)) if annotation_indices is None else annotation_indices)
    spectrum = SimpleNamespace(
        spectrum_id="S1", scan_index=1, rt=1.0, precursor_mz=500.0,
        precursor_charge=2, peaks=list(zip(mz, parsed or [], strict=False)), raw_peaks=list(zip(mz, parsed or [], strict=False)),
    )
    return {"source_spectra": [source], "spectra": [spectrum]}


def _ranking(mod="M", position=36, parent_position=1):
    return [{
        "Rank": 1, "Modification_ID": mod, "Modification_Name": mod,
        "Parent_Fragment_ID": "P", "Candidate_tRNA_Position": position,
        "Candidate_Position_In_Parent": parent_position, "Final_Score": 5.0,
        "Final_Confidence": "Low",
    }]


def _ambiguity(mz_values, mod="M", position=36, statuses=None):
    cluster = {
        "Peak_Cluster_ID": "CL1", "Spectrum_ID": "S1", "Modification_ID": mod,
        "Parent_Fragment_ID": "P", "Candidate_tRNA_Position": position,
        "Candidate_Position_In_Parent": 1, "Ion_Series": "y", "Ion_Number": 3,
    }
    details = []
    for index, mz in enumerate(mz_values):
        details.append({
            "Peak_Cluster_ID": "CL1", "Physical_Observed_Peak_Key": _key(mz),
            "Error_ppm": index + 1.0, "Within_Formal_Tolerance": index == 0,
            "Candidate_Specificity_Status": (statuses or ["candidate_specific"] * len(mz_values))[index],
        })
    return [cluster], details


def _build(intensities, *, parsed=None, annotation_indices=None, mod="M", ambiguity=True, max_rows=1_048_573,
           modified_matches=None, identity=None, localization=None):
    mz = [100.0 + index for index in range(len(intensities))]
    context = _context(intensities, parsed=intensities if parsed is None else parsed, annotation_indices=annotation_indices, mz=mz)
    clusters, details = _ambiguity(mz, mod=mod) if ambiguity else ([], [])
    return build_zero_intensity_audit(
        context, _ranking(mod=mod), [], modified_matches or [], identity or [], localization or [], clusters, details,
        max_detail_rows=max_rows,
    )


def test_all_positive_peaks():
    spectra, detail, summary, _, _ = _build([1.0, 2.0])
    assert spectra[0]["Positive_Intensity_Count"] == 2
    assert summary[0]["Total_Zero_Intensity_Peaks"] == 0
    assert all(row["Intensity_State"] == "positive" for row in detail)


def test_some_zero_peaks():
    _, _, summary, _, _ = _build([0.0, 2.0])
    assert summary[0]["Total_Zero_Intensity_Peaks"] == 1
    assert summary[0]["Zero_Intensity_Fraction"] == 0.5


def test_all_zero_peaks_and_spectrum():
    spectra, _, summary, _, _ = _build([0.0, 0.0])
    assert spectra[0]["Zero_Intensity_Count"] == 2
    assert summary[0]["All_Zero_Spectra"] == 1


def test_nan_intensity():
    assert intensity_state(nan) == "nan"
    _, detail, summary, _, _ = _build([nan], parsed=[nan])
    assert detail[0]["Intensity_State"] == "nan"
    assert summary[0]["Total_NaN_Peaks"] == 1


def test_none_intensity():
    context = _context([None], parsed=None)
    _, detail, summary, _, _ = build_zero_intensity_audit(context)
    assert detail[0]["Intensity_State"] == "missing"
    assert summary[0]["Total_Missing_Peaks"] == 1


def test_negative_intensity():
    _, detail, summary, _, _ = _build([-1.0])
    assert detail[0]["Intensity_State"] == "negative"
    assert summary[0]["Total_Negative_Peaks"] == 1


def test_non_numeric_intensity():
    context = _context(["bad"], parsed=None)
    _, detail, summary, _, _ = build_zero_intensity_audit(context)
    assert detail[0]["Intensity_State"] == "non_numeric"
    assert summary[0]["Total_Non_Numeric_Peaks"] == 1


def test_raw_and_annotation_input_equal():
    _, detail, _, _, diagnostics = _build([2.0], annotation_indices=[0])
    assert detail[0]["Original_Intensity"] == detail[0]["Annotation_Input_Intensity"] == 2.0
    assert diagnostics[0]["Annotation_Input_Peak_Count"] == 1


def test_parsing_introduced_zero_classification():
    assert _origin_for_values(5.0, 0.0, 0.0) == "introduced_during_parsing"


def test_annotation_preparation_introduced_zero_classification():
    assert _origin_for_values(5.0, 5.0, 0.0) == "introduced_during_annotation_preparation"


def _zero_modified_match():
    return {
        "Spectrum_ID": "S1", "RT": 1.0, "Observed_mz": 100.0, "Observed_Intensity": 0.0,
        "Modification_ID": "M", "Parent_Fragment_ID": "P",
        "Candidate_Modification_Position_In_Parent": 1, "Ion_Type": "y", "Ion_End": 3,
        "Mass_Error_ppm": 1.0, "Discriminates_Position": True,
    }


def test_zero_peak_used_for_theoretical_match():
    _, detail, summary, _, _ = _build([0.0], modified_matches=[_zero_modified_match()])
    assert detail[0]["Used_For_Theoretical_Match"] is True
    assert summary[0]["Matched_Zero_Intensity_Peaks"] == 1


def test_zero_peak_selected_as_best_match():
    _, detail, summary, _, _ = _build([0.0], modified_matches=[_zero_modified_match()])
    assert detail[0]["Selected_As_Best_Match"] is True
    assert summary[0]["Zero_Intensity_Best_Matches"] == 1


def test_zero_peak_used_for_identity():
    identity = [{
        "Physical_Observed_Peak_Key": _key(100.0), "Modification_ID": "M",
        "Parent_Fragment_ID": "P", "Candidate_tRNA_Position": 36,
    }]
    _, detail, summary, _, _ = _build([0.0], identity=identity)
    assert detail[0]["Used_For_Identity"] is True
    assert summary[0]["Zero_Intensity_Identity_Assignments"] == 1


def test_zero_peak_used_for_localization():
    localization = [{
        "Modification_ID": "M", "Parent_Fragment_ID": "P",
        "Candidate_Modification_Position_In_tRNA": 36, "Num_Modified_Ion_Matches": 1,
    }]
    _, detail, summary, _, _ = _build([0.0], modified_matches=[_zero_modified_match()], localization=localization)
    assert detail[0]["Used_For_Localization"] is True
    assert summary[0]["Zero_Intensity_Localization_Uses"] == 1


def test_zero_peak_used_in_ambiguity_cluster():
    _, detail, summary, _, _ = _build([0.0])
    assert detail[0]["Used_In_Ambiguity_Cluster"] is True
    assert summary[0]["Ambiguity_Clusters_With_Zero_Intensity"] == 1


def test_all_zero_cluster():
    _, _, summary, candidates, _ = _build([0.0, 0.0])
    assert summary[0]["All_Zero_Ambiguity_Clusters"] == 1
    assert candidates[0]["All_Zero_Cluster_Count"] == 1


def test_mixed_zero_positive_cluster():
    _, _, summary, candidates, _ = _build([0.0, 2.0])
    assert summary[0]["All_Zero_Ambiguity_Clusters"] == 0
    assert candidates[0]["Zero_Intensity_Audit_Severity"] == "moderate"


def test_candidate_specific_zero_peak():
    _, _, _, candidates, _ = _build([0.0])
    assert candidates[0]["Zero_Intensity_Cluster_Count"] == 1
    assert candidates[0]["Shadow_Nonzero_Candidate_Specific_Count"] == 0


def test_shared_zero_peak():
    context = _context([0.0], parsed=[0.0])
    clusters, details = _ambiguity([100.0], statuses=["cross_candidate_shared"])
    _, _, summary, _, _ = build_zero_intensity_audit(context, _ranking(), ambiguity_clusters=clusters, ambiguity_details=details)
    assert summary[0]["Ambiguity_Clusters_With_Zero_Intensity"] == 1


def test_cnm5u_positions_36_37_38_are_independent_candidates():
    context = _context([0.0], parsed=[0.0])
    ranking = [_ranking("cnm5U", position)[0] for position in (36, 37, 38)]
    clusters = []
    details = []
    for position in (36, 37, 38):
        c, d = _ambiguity([100.0], mod="cnm5U", position=position)
        c[0]["Peak_Cluster_ID"] = f"CL{position}"
        d[0]["Peak_Cluster_ID"] = f"CL{position}"
        clusters += c
        details += d
    _, _, summary, candidates, _ = build_zero_intensity_audit(
        context, ranking, ambiguity_clusters=clusters, ambiguity_details=details,
    )
    assert summary[0]["Affected_cnm5U_Count"] == 3
    assert {row["Candidate_tRNA_Position"] for row in candidates if row["Modification_ID"] == "cnm5U"} == {"36", "37", "38"}


def test_nonzero_shadow_simulation():
    _, _, summary, candidates, _ = _build([0.0, 2.0])
    assert summary[0]["Shadow_Nonzero_Ambiguity_Cluster_Count"] == 0
    assert candidates[0]["Shadow_Nonzero_Conclusion"] == "some_positive_intensity_evidence_remains"


def test_formal_score_confidence_rank_inputs_unchanged():
    ranking = _ranking()
    before = deepcopy(ranking)
    build_zero_intensity_audit(_context([0.0], parsed=[0.0]), ranking)
    assert ranking == before
    assert ranking[0]["Final_Score"] == 5.0
    assert ranking[0]["Final_Confidence"] == "Low"
    assert ranking[0]["Rank"] == 1


def test_top50_existing_columns_unchanged():
    ranking = _ranking()
    base = _build_top_candidates(__import__("pandas").DataFrame(ranking), __import__("pandas").DataFrame(), {}, None, None)
    zero = __import__("pandas").DataFrame(_build([0.0])[3])
    after = _build_top_candidates(__import__("pandas").DataFrame(ranking), __import__("pandas").DataFrame(), {}, None, None, zero)
    old_columns = [column for column in TOP_CANDIDATE_COLUMNS if column not in TOP_SHADOW_COLUMNS]
    assert base[old_columns].equals(after[old_columns])


def test_existing_diagnostics_values_unchanged_when_appending():
    existing = {"Existing": 7}
    before = deepcopy(existing)
    diagnostics = _build([0.0])[4][0]
    existing.update(diagnostics)
    assert {key: existing[key] for key in before} == before
    assert list(diagnostics) == DIAGNOSTIC_COLUMNS


def test_internal_context_name_is_not_an_excel_sheet_name():
    output_names = {"MS2_Zero_Intensity_Spectra", "MS2_Zero_Intensity_Detail", "MS2_Zero_Intensity_Summary"}
    assert "_MS2_Zero_Intensity_Audit_Context" not in output_names
    assert "_MS2_Zero_Intensity_Candidate_Summary" not in output_names


def test_sheet_names_within_excel_limit():
    names = ["MS2_Zero_Intensity_Spectra", "MS2_Zero_Intensity_Detail", "MS2_Zero_Intensity_Summary"]
    assert all(len(name) <= 31 for name in names)


def test_detail_truncation_is_consistent():
    _, detail, summary, _, _ = _build([0.0, 1.0, 2.0], max_rows=2)
    assert len(detail) == 2
    assert summary[0]["Original_Detail_Row_Count"] == 3
    assert summary[0]["Written_Detail_Row_Count"] == 2
    assert summary[0]["Detail_Truncated"] is True


def test_deterministic_row_order():
    _, detail, _, _, _ = _build([1.0, 0.0, 2.0])
    assert [row["Peak_Index"] for row in detail] == [0, 1, 2]


def test_identical_rerun_is_equal():
    first = _build([0.0, 2.0])
    second = _build([0.0, 2.0])
    assert first == second


def test_required_column_contracts():
    assert "Median_Positive_Intensity" in SPECTRA_COLUMNS
    assert "Transformation_History" in DETAIL_COLUMNS
    assert "Applied_To_Final_Score" in SUMMARY_COLUMNS
    assert "Formal_Best_Match_Definition" in SUMMARY_COLUMNS
    assert "Ambiguity_Diagnostic_Best_Peak_Definition" in SUMMARY_COLUMNS
    summary = _build([0.0])[2][0]
    assert "intensity filtering" in summary["Formal_Best_Match_Definition"]
    assert "not a formal match" in summary["Ambiguity_Diagnostic_Best_Peak_Definition"]
    assert TOP_SHADOW_COLUMNS[-1] == "Zero_Intensity_Audit_Applied_To_Final_Score"


def test_exact_zero_origin_is_present_in_decoded_mzml_array():
    spectra, _, summary, _, _ = _build([0.0])
    assert spectra[0]["Origin_Category"] == "present_in_raw_mzml"
    assert summary[0]["Likely_Origin_Category"] == "present_in_raw_mzml"
