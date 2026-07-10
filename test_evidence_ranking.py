from types import SimpleNamespace

from rna_masshunter.evidence_ranking import (
    _confidence_limiting_factor,
    _final_confidence,
    build_ambiguity_groups,
    build_modification_evidence_ranking,
)
from rna_masshunter.models import Fragment, Modification


def _config():
    return SimpleNamespace(
        organism={"group": "archaea", "species": "Test archaeon", "rule_set": "test"},
        sequence={"type": "RNA", "anticodon": "UUC"},
        modification_evidence_ranking={
            "enabled": True, "use_ms1_fragment_evidence": True,
            "use_known_modification_candidates": True, "use_ms2_precursor_evidence": True,
            "use_ms2_modified_ion_evidence": True, "use_localization_evidence": True,
            "use_organism_rules": True, "use_trna_context": True,
            "require_ms2_evidence_for_high_confidence": True,
            "min_final_score_to_report": -10, "max_ranked_candidates": 100,
            "weights": {
                "ms1_fragment_match": 1, "known_modification_candidate": 1,
                "ms2_precursor_rescue": 2, "ms2_modified_ion_match": 2,
                "localization_weak": 1, "localization_moderate": 3, "localization_strong": 5,
                "organism_rule_supported": 1.5, "trna_context_supported": 1,
                "low_information_penalty": -1, "ambiguous_position_penalty": -1,
                "isobaric_precursor_penalty": -2,
            },
        },
    )


def test_integrated_candidate_ranks_first_and_penalties_apply():
    modifications = [
        Modification("mA", "mA", 14.0, "methyl", ["A"], raw={
            "name": "A modification", "source_priority": "user_pdf_for_mass_shift",
            "curation_status": "manually_checked", "candidate_policy": {"include_by_mass_search": True, "include_if_position_rule_exists": True},
            "detectability": {"ms1": True, "ms2": True}, "chemical_group": "methylation",
            "near_isobaric_group": "near_test",
        }),
        Modification("mB", "mB", 28.0, "other", ["C"], raw={"name": "B modification"}),
        Modification("mI", "mI", 0.0, "isobaric", ["U"], raw={"name": "Isobaric"}),
        Modification("mC", "mC", 15.0, "other", ["C"], raw={"name": "C modification"}),
    ]
    fragments = [Fragment("F1", "T", "AAC", 10, 12, None, None, "RNase_T1", 0, "default", 1000.0)]
    ms1 = [{"fragment_id": "F1", "confidence": "High", "intensity": 5000.0}]
    known = [
        {"modification_id": "mA", "source_type": "fragment", "source_id": "F1", "priority_score": 5.0},
        {"modification_id": "mB", "source_type": "fragment", "source_id": "", "priority_score": 1.0},
    ]
    localization = [
        {"Modification_ID": "mA", "Modification_Name": "A modification", "Parent_Fragment_ID": "F1",
         "Candidate_Modification_Position_In_Parent": 1, "Candidate_Modification_Position_In_tRNA": 10,
         "Candidate_Modification_Base": "A", "Parent_Sequence": "AAC", "Parent_Start": 10, "Parent_End": 12,
         "Localization_Level": "Moderate", "Localization_Score": 5.0,
         "Localization_Interpretation": "modification-supported-on-position", "Num_c_Modified_Ions": 1, "Num_y_Modified_Ions": 1,
         "Has_Position_Discriminating_Evidence": True, "Num_Position_Discriminating_Modified_Ions": 1,
         "Num_Informative_Position_Discriminating_Modified_Ions": 1},
        {"Modification_ID": "mC", "Modification_Name": "C modification", "Parent_Fragment_ID": "F1",
         "Candidate_Modification_Position_In_Parent": 3, "Candidate_Modification_Position_In_tRNA": 12,
         "Candidate_Modification_Base": "C", "Parent_Sequence": "AAC", "Parent_Start": 10, "Parent_End": 12,
         "Localization_Level": "Weak", "Localization_Score": 0.25,
         "Localization_Interpretation": "ambiguous-multiple-positions", "Num_c_Modified_Ions": 1, "Num_y_Modified_Ions": 0},
    ]
    precursors = [
        {"Modification_ID": "mA", "Modification_Name": "A modification", "Parent_Fragment_ID": "F1", "Precursor_Error_ppm": 1.0, "Modified_Precursor_Rescue": True},
        {"Modification_ID": "mI", "Modification_Name": "Isobaric", "Parent_Fragment_ID": "F1", "Precursor_Error_ppm": 1.0, "Modified_Precursor_Rescue": False},
        {"Modification_ID": "mC", "Modification_Name": "C modification", "Parent_Fragment_ID": "F1", "Precursor_Error_ppm": 2.0, "Modified_Precursor_Rescue": True},
    ]
    ion_matches = [
        {"Modification_ID": "mA", "Parent_Fragment_ID": "F1", "Candidate_Modification_Position_In_Parent": 1,
         "Ion_Contains_Modification": True, "Informative_Ion": True, "Mass_Error_ppm": 1.0},
        {"Modification_ID": "mA", "Parent_Fragment_ID": "F1", "Candidate_Modification_Position_In_Parent": 1,
         "Ion_Contains_Modification": True, "Informative_Ion": True, "Mass_Error_ppm": 2.0},
        {"Modification_ID": "mC", "Parent_Fragment_ID": "F1", "Candidate_Modification_Position_In_Parent": 3,
         "Ion_Contains_Modification": True, "Informative_Ion": False, "Mass_Error_ppm": 2.0},
    ]
    rows, summary = build_modification_evidence_ranking(
        _config(), modifications, fragments, ms1, known,
        {"MS2_Modification_Localization_Evidence": localization,
         "MS2_Modified_Precursor_Candidates": precursors,
         "MS2_Modified_Ion_Matches": ion_matches},
        rule_set={},
    )
    assert rows[0]["Modification_ID"] == "mA"
    assert rows[0]["Final_Confidence"] in {"Very High", "High"}
    assert rows[0]["Source_Priority"] == "user_pdf_for_mass_shift"
    assert rows[0]["Curation_Status"] == "manually_checked"
    assert rows[0]["Candidate_Policy_By_Mass_Search"] is True
    assert rows[0]["Detectability_MS2"] is True
    assert rows[0]["Chemical_Group"] == "methylation"
    known_only = next(row for row in rows if row["Modification_ID"] == "mB")
    assert known_only["Final_Confidence"] not in {"Very High", "High"}
    ambiguous = next(row for row in rows if row["Modification_ID"] == "mC" and row["Candidate_Position_In_Parent"] == 3)
    assert ambiguous["Ambiguous_Position"] is True
    assert ambiguous["Low_Information_Evidence"] is True
    isobaric = next(row for row in rows if row["Modification_ID"] == "mI")
    assert isobaric["Is_Isobaric"] is True
    assert isobaric["Final_Confidence"] not in {"Very High", "High"}
    assert summary[0]["Total_Ranked_Candidates"] == len(rows)


def test_confidence_calibration_for_weak_and_multi_ion_support():
    config = {"require_ms2_evidence_for_high_confidence": True}
    weak_single = _final_confidence(
        8.0, True, True, "Weak", True, True, False,
        1, 1, 0, 1.0, 1.0, False, 0, "ambiguous", config,
    )
    assert weak_single == "Medium"
    factors = _confidence_limiting_factor(True, True, "Weak", 1, 1, 0, True, True, False, "ambiguous")
    assert "weak-localization" in factors
    assert "single-modified-ion" in factors
    assert "one-sided-ion-series" in factors
    assert "position-ambiguous" in factors
    assert "no-position-discriminating-ion" in factors

    multi_series = _final_confidence(
        7.0, True, True, "Weak", True, True, False,
        2, 1, 1, 2.0, 2.0, True, 2, "resolved", config,
    )
    assert multi_series == "High"

    strong = _final_confidence(
        9.0, True, True, "Strong", True, True, False,
        3, 2, 1, 2.0, 2.0, True, 2, "resolved", config,
    )
    assert strong == "Very High"


def test_ambiguity_groups_distinguish_shared_and_position_specific_evidence():
    localization = []
    for spectrum_id in ("AMB", "RES"):
        for position in (1, 3):
            localization.append({
                "Spectrum_ID": spectrum_id, "Parent_Fragment_ID": "F1", "Modification_ID": "mA",
                "Modification_Name": "A modification", "Parent_Sequence": "AUA", "Parent_Start": 10,
                "Parent_End": 12, "Candidate_Modification_Position_In_Parent": position,
                "Candidate_Modification_Position_In_tRNA": 9 + position, "Candidate_Modification_Base": "A",
                "Localization_Score": 2.0 if spectrum_id == "AMB" else (4.0 if position == 1 else 1.0),
                "RT": 1.0, "Precursor_mz": 500.0,
            })
    common = {"Parent_Fragment_ID": "F1", "Modification_ID": "mA", "Ion_Contains_Modification": True, "Informative_Ion": True}
    matches = [
        {**common, "Spectrum_ID": "AMB", "Candidate_Modification_Position_In_Parent": 1, "Position_Discriminating_Ion": False},
        {**common, "Spectrum_ID": "AMB", "Candidate_Modification_Position_In_Parent": 3, "Position_Discriminating_Ion": False},
        {**common, "Spectrum_ID": "RES", "Candidate_Modification_Position_In_Parent": 1, "Position_Discriminating_Ion": True},
        {**common, "Spectrum_ID": "RES", "Candidate_Modification_Position_In_Parent": 3, "Position_Discriminating_Ion": False},
    ]
    groups = build_ambiguity_groups(localization, matches)
    assert next(row for row in groups if row["Spectrum_ID"] == "AMB")["Position_Ambiguity_Status"] == "ambiguous"
    assert next(row for row in groups if row["Spectrum_ID"] == "RES")["Position_Ambiguity_Status"] == "resolved_by_discriminating_ions"


if __name__ == "__main__":
    test_integrated_candidate_ranks_first_and_penalties_apply()
    test_confidence_calibration_for_weak_and_multi_ion_support()
    test_ambiguity_groups_distinguish_shared_and_position_specific_evidence()
    print("synthetic evidence ranking test: OK")
