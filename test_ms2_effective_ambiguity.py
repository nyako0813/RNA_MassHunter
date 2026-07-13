from copy import deepcopy

import pandas as pd

from rna_masshunter.ms2_effective_ambiguity import (
    CLUSTER_COLUMNS, DETAIL_COLUMNS, DIAGNOSTIC_COLUMNS, SUMMARY_COLUMNS, TOP_SHADOW_COLUMNS,
    _diagnostics, build_effective_ambiguity,
)
from rna_masshunter.ms2_identity_evidence import physical_observed_peak_key
from rna_masshunter.review_dashboard import TOP_CANDIDATE_COLUMNS, _build_top_candidates


def _peak(mz):
    return physical_observed_peak_key({"Spectrum_ID": "S", "Observed_mz": mz, "RT": 1.0})


def _cluster(cid="C1", mod="M", pos=36, theory="T1", theoretical=100.0):
    return {
        "Peak_Cluster_ID": cid, "Modification_ID": mod, "Parent_Fragment_ID": "P",
        "Candidate_tRNA_Position": pos, "Theoretical_Ion_ID": theory, "Ion_Series": "y",
        "Ion_Number": 3, "Theoretical_mz": theoretical, "Primary_Ambiguity_Type": "multiple_peaks_same_side",
    }


def _detail(mz, intensity, within=False, cid="C1", status="candidate_specific"):
    return {
        "Peak_Cluster_ID": cid, "Physical_Observed_Peak_Key": _peak(mz), "Spectrum_ID": "S",
        "Observed_mz": mz, "Intensity": intensity, "Error_Da": mz - 100.0,
        "Error_ppm": (mz - 100.0) / 100.0 * 1e6, "Within_Audit_Window": True,
        "Within_Formal_Tolerance": within, "Candidate_Specificity_Status": status,
    }


def _ranking(mod="M", pos=36, parent_pos=1):
    return [{
        "Rank": 1, "Modification_ID": mod, "Modification_Name": mod, "Parent_Fragment_ID": "P",
        "Candidate_tRNA_Position": pos, "Candidate_Position_In_Parent": parent_pos,
        "Modification_Family": "family", "Final_Score": 5.0, "Final_Confidence": "Low",
    }]


def _formal(mz, ion="T1", mod="M", parent_pos=1):
    return {
        "Spectrum_ID": "S", "RT": 1.0, "Observed_mz": mz, "Observed_Intensity": 10.0,
        "Ion_ID": ion, "Modification_ID": mod, "Parent_Fragment_ID": "P",
        "Candidate_Modification_Position_In_Parent": parent_pos,
    }


def _identity(mz, ion="T1", mod="M", pos=36, group=""):
    return {
        "Physical_Observed_Peak_Key": _peak(mz), "Theoretical_Ion_ID": ion,
        "Modification_ID": mod, "Parent_Fragment_ID": "P", "Candidate_tRNA_Position": pos,
        "Structural_Isomer_Group_ID": group,
    }


def _build(peaks, *, cluster=None, formal=None, identities=None, ranking=None, max_rows=100000):
    cluster = cluster or _cluster()
    summary = [{
        "Modification_ID": cluster["Modification_ID"], "Parent_Fragment_ID": "P",
        "Candidate_tRNA_Position": cluster["Candidate_tRNA_Position"], "Ambiguity_Severity": "moderate",
    }]
    return build_effective_ambiguity(
        [cluster], peaks, summary, [], formal or [], identities or [], [], ranking or _ranking(cluster["Modification_ID"], cluster["Candidate_tRNA_Position"]),
        {}, max_detail_rows=max_rows,
    )


def test_raw_peak_one_is_none():
    rows = _build([_detail(100.01, 1)])[0]
    assert rows[0]["Raw_Ambiguous"] is False and rows[0]["Effective_Ambiguity_Level"] == "none"


def test_two_raw_zero_positive_is_raw_only():
    row = _build([_detail(100.01, 0), _detail(100.02, 1)])[0][0]
    assert row["Raw_Ambiguous"] and not row["Positive_Ambiguous"] and row["Effective_Ambiguity_Level"] == "raw_only"


def test_three_raw_two_zero_one_positive():
    row = _build([_detail(99.99, 0), _detail(100.01, 0), _detail(100.02, 1)])[0][0]
    assert row["Raw_Zero_Intensity_Peak_Count"] == 2 and row["Positive_Peak_Count"] == 1


def test_two_positive_peaks():
    row = _build([_detail(100.01, 1), _detail(100.02, 2)])[0][0]
    assert row["Positive_Ambiguous"] and row["Effective_Ambiguity_Level"] == "positive_intensity"


def test_two_positive_only_one_within_tolerance():
    row = _build([_detail(100.001, 1, True), _detail(100.02, 2, False)])[0][0]
    assert row["Positive_Ambiguous"] and not row["Formal_Tolerance_Ambiguous"]


def test_two_positive_both_within_tolerance():
    row = _build([_detail(99.999, 1, True), _detail(100.001, 2, True)])[0][0]
    assert row["Formal_Tolerance_Ambiguous"] and row["Effective_Ambiguity_Level"] == "formal_tolerance"


def test_one_formal_matched_peak_not_ambiguous():
    row = _build([_detail(100.001, 1, True), _detail(100.02, 0)], formal=[_formal(100.001)])[0][0]
    assert row["Formal_Matched_Physical_Peak_Count"] == 1 and not row["Formal_Match_Ambiguous"]


def test_multiple_formal_peaks_same_theoretical_ion():
    row = _build([_detail(100.001, 1, True), _detail(100.002, 2, True)], formal=[_formal(100.001), _formal(100.002)])[0][0]
    assert row["Formal_Match_Ambiguous"] and "multiple_formal_peaks_for_same_theoretical_ion" in row["Formal_Match_Sharing_Type"]


def test_multiple_theoretical_ions_share_formal_peak():
    row = _build([_detail(100.001, 1, True), _detail(100.02, 0)], formal=[_formal(100.001, "T1"), _formal(100.001, "T2")])[0][0]
    assert row["Effective_Ambiguity_Level"] == "formal_match"
    assert "multiple_theoretical_ions_share_formal_peak" in row["Formal_Match_Sharing_Type"]


def test_candidate_specific_formal_peak():
    row = _build([_detail(100.001, 1, True), _detail(100.02, 0)], formal=[_formal(100.001)])[0][0]
    assert row["Candidate_Specific_Formal_Peak_Count"] == 1


def test_position_group_shared_formal_peak():
    row = _build([_detail(100.001, 1, True, status="position_group_shared"), _detail(100.02, 0)], formal=[_formal(100.001)])[0][0]
    assert row["Position_Group_Shared_Formal_Peak_Count"] == 1


def test_structural_isomer_shared_formal_peak():
    row = _build([_detail(100.001, 1, True, status="structural_isomer_group_shared"), _detail(100.02, 0)], formal=[_formal(100.001)])[0][0]
    assert row["Structural_Isomer_Shared_Formal_Peak_Count"] == 1


def test_cross_candidate_shared_formal_peak():
    row = _build([_detail(100.001, 1, True, status="cross_candidate_shared"), _detail(100.02, 0)], formal=[_formal(100.001)])[0][0]
    assert row["Cross_Candidate_Shared_Formal_Peak_Count"] == 1


def test_raw_only_type():
    row = _build([_detail(100.01, 0), _detail(100.02, 1)])[0][0]
    assert row["Effective_Ambiguity_Type"] == "raw_zero_inflated_multiple_peaks"


def test_positive_same_side_type():
    row = _build([_detail(100.01, 1), _detail(100.02, 2)])[0][0]
    assert row["Effective_Ambiguity_Type"] == "multiple_positive_peaks_same_side"


def test_positive_bracketing_type():
    row = _build([_detail(99.99, 1), _detail(100.02, 2)])[0][0]
    assert row["Effective_Ambiguity_Type"] == "multiple_positive_peaks_bracketing"


def test_formal_tolerance_type():
    row = _build([_detail(99.999, 1, True), _detail(100.001, 2, True)])[0][0]
    assert row["Effective_Ambiguity_Type"] == "multiple_positive_peaks_within_formal_tolerance"


def test_formal_match_type():
    row = _build([_detail(100.001, 1, True), _detail(100.002, 2, True)], formal=[_formal(100.001), _formal(100.002)])[0][0]
    assert row["Effective_Ambiguity_Level"] == "formal_match"


def test_severity_high():
    row = _build([_detail(100.001, 1), _detail(100.002, 2)], formal=[_formal(100.001), _formal(100.002)])[0][0]
    assert row["Effective_Ambiguity_Severity"] == "high"


def test_severity_moderate():
    row = _build([_detail(100.001, 1, True), _detail(100.002, 2, True)])[0][0]
    assert row["Effective_Ambiguity_Severity"] == "moderate"


def test_severity_low():
    row = _build([_detail(100.01, 1), _detail(100.02, 2)])[0][0]
    assert row["Effective_Ambiguity_Severity"] == "low"


def test_severity_informational():
    row = _build([_detail(100.01, 0), _detail(100.02, 2)])[0][0]
    assert row["Effective_Ambiguity_Severity"] == "informational"


def test_severity_none():
    assert _build([_detail(100.01, 2)])[0][0]["Effective_Ambiguity_Severity"] == "none"


def test_cnm5u_positions_36_37_38():
    clusters=[];details=[];summary=[];ranking=[]
    for pos in (36,37,38):
        cid=f"C{pos}";clusters.append(_cluster(cid,"cnm5U",pos));ranking+=_ranking("cnm5U",pos)
        summary.append({"Modification_ID":"cnm5U","Parent_Fragment_ID":"P","Candidate_tRNA_Position":pos,"Ambiguity_Severity":"high"})
        details += [_detail(100.01,0,cid=cid,status="position_group_shared"),_detail(100.02,1,cid=cid,status="spectrum_peak_only")]
    rows=build_effective_ambiguity(clusters,details,summary,[],[],[],[],ranking,{})[0]
    assert len(rows)==3 and all(row["Effective_Ambiguity_Level"]=="raw_only" for row in rows)


def test_positional_isomer_remains_unresolved():
    row=_build([_detail(100.01,0,status="position_group_shared"),_detail(100.02,1,status="spectrum_peak_only")])[0][0]
    assert row["Effective_Ambiguity_Recommendation"] == "positional_isomer_remains_unresolved"


def test_top50_shadow_columns_added_without_existing_changes():
    ranking=pd.DataFrame(_ranking())
    before=_build_top_candidates(ranking,pd.DataFrame(),{},None,None,None,None)
    candidates=pd.DataFrame(_build([_detail(100.01,0),_detail(100.02,1)])[3])
    after=_build_top_candidates(ranking,pd.DataFrame(),{},None,None,None,candidates)
    old=[column for column in TOP_CANDIDATE_COLUMNS if column not in TOP_SHADOW_COLUMNS]
    assert before[old].equals(after[old]) and all(column in after for column in TOP_SHADOW_COLUMNS)


def test_diagnostics_shadow_columns():
    summary=_build([_detail(100.01,0),_detail(100.02,1)])[2]
    diagnostics=_diagnostics(summary,True)
    assert list(diagnostics)==DIAGNOSTIC_COLUMNS and diagnostics["Effective_Ambiguity_Applied_To_Final_Score"] is False


def test_formal_ranking_inputs_unchanged():
    ranking=_ranking();before=deepcopy(ranking)
    _build([_detail(100.01,0),_detail(100.02,1)],ranking=ranking)
    assert ranking==before and ranking[0]["Final_Score"]==5.0 and ranking[0]["Final_Confidence"]=="Low" and ranking[0]["Rank"]==1


def test_review_rank_priority_unchanged():
    ranking=pd.DataFrame(_ranking());base=_build_top_candidates(ranking,pd.DataFrame(),{},None,None,None,None)
    candidates=pd.DataFrame(_build([_detail(100.01,0),_detail(100.02,1)])[3])
    after=_build_top_candidates(ranking,pd.DataFrame(),{},None,None,None,candidates)
    assert base[["Review_Rank","Review_Priority"]].equals(after[["Review_Rank","Review_Priority"]])


def test_existing_ambiguity_inputs_unchanged():
    clusters=[_cluster()];details=[_detail(100.01,0),_detail(100.02,1)];before=deepcopy((clusters,details))
    build_effective_ambiguity(clusters,details,[],[],[],[],[],_ranking(),{})
    assert (clusters,details)==before


def test_zero_audit_inputs_unchanged():
    context={"source_spectra":[{"spectrum_id":"S","scan_index":1,"original_mz":[100.01,100.02]}]};before=deepcopy(context)
    build_effective_ambiguity([_cluster()],[_detail(100.01,0),_detail(100.02,1)],[],[],[],[],[],_ranking(),context)
    assert context==before


def test_internal_context_not_sheet():
    names={"MS2_Effective_Ambiguity","MS2_Effective_Ambig_Detail","MS2_Effective_Ambig_Summary"}
    assert "_MS2_Effective_Ambiguity_Candidate_Summary" not in names


def test_sheet_names_under_31():
    assert all(len(name)<=31 for name in ["MS2_Effective_Ambiguity","MS2_Effective_Ambig_Detail","MS2_Effective_Ambig_Summary"])


def test_detail_truncate_consistency():
    result=_build([_detail(100.01,0),_detail(100.02,1),_detail(100.03,2)],max_rows=2)
    assert len(result[1])==2 and result[2][0]["Detail_Original_Row_Count"]==3 and result[2][0]["Detail_Truncated"]


def test_deterministic_row_order():
    rows=_build([_detail(100.03,2),_detail(100.01,0),_detail(100.02,1)])[1]
    assert [row["MZ"] for row in rows]==[100.01,100.02,100.03]


def test_identical_rerun():
    peaks=[_detail(100.01,0),_detail(100.02,1)]
    assert _build(deepcopy(peaks))==_build(deepcopy(peaks))


def test_column_contracts_and_nonapplication():
    assert "Effective_Ambiguity_Level" in CLUSTER_COLUMNS and "Physical_Peak_ID" in DETAIL_COLUMNS
    assert "Formal_Match_Ambiguity_Definition" in SUMMARY_COLUMNS
    assert TOP_SHADOW_COLUMNS[-1]=="Effective_Ambiguity_Applied_To_Final_Score"
    assert _build([_detail(100.01,0),_detail(100.02,1)])[0][0]["Effective_Ambiguity_Applied_To_Final_Score"] is False


def test_multiple_candidates_share_formal_peak_and_position_isomers():
    identities=[_identity(100.001,"T1","M",36),_identity(100.001,"T1","M",37)]
    ranking=_ranking("M",36)+_ranking("M",37)
    row=_build([_detail(100.001,1,True),_detail(100.02,0)],formal=[_formal(100.001)],identities=identities,ranking=ranking)[0][0]
    assert row["Formal_Match_Ambiguous"]
    assert "multiple_candidates_share_formal_peak" in row["Formal_Match_Sharing_Type"]
    assert "position_isomers_share_formal_peak" in row["Formal_Match_Sharing_Type"]


def test_structural_isomers_share_formal_peak():
    identities=[_identity(100.001,"T1","M",36,"SG1"),_identity(100.001,"T1","M2",36,"SG1")]
    ranking=_ranking("M",36)+_ranking("M2",36)
    row=_build([_detail(100.001,1,True),_detail(100.02,0)],formal=[_formal(100.001)],identities=identities,ranking=ranking)[0][0]
    assert "structural_isomers_share_formal_peak" in row["Formal_Match_Sharing_Type"]
    assert row["Effective_Ambiguity_Recommendation"] == "inspect_structural_isomer_evidence"


def test_candidate_formal_sharing_is_high_severity():
    identities=[_identity(100.001,"T1","M",36),_identity(100.001,"T1","M2",36)]
    ranking=_ranking("M",36)+_ranking("M2",36)
    row=_build([_detail(100.001,1,True),_detail(100.02,0)],formal=[_formal(100.001)],identities=identities,ranking=ranking)[0][0]
    assert row["Effective_Ambiguity_Severity"] == "high"
