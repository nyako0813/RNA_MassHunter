from types import SimpleNamespace

from rna_masshunter.evidence_ranking import build_modification_evidence_ranking
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
        Modification("mA", "mA", 14.0, "methyl", ["A"], raw={"name": "A modification"}),
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
         "Localization_Interpretation": "modification-supported-on-position", "Num_c_Modified_Ions": 1, "Num_y_Modified_Ions": 1},
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
    known_only = next(row for row in rows if row["Modification_ID"] == "mB")
    assert known_only["Final_Confidence"] not in {"Very High", "High"}
    ambiguous = next(row for row in rows if row["Modification_ID"] == "mC" and row["Candidate_Position_In_Parent"] == 3)
    assert ambiguous["Ambiguous_Position"] is True
    assert ambiguous["Low_Information_Evidence"] is True
    isobaric = next(row for row in rows if row["Modification_ID"] == "mI")
    assert isobaric["Is_Isobaric"] is True
    assert isobaric["Final_Confidence"] not in {"Very High", "High"}
    assert summary[0]["Total_Ranked_Candidates"] == len(rows)


if __name__ == "__main__":
    test_integrated_candidate_ranks_first_and_penalties_apply()
    print("synthetic evidence ranking test: OK")
