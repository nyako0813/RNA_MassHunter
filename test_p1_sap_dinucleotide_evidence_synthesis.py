from copy import deepcopy

from rna_masshunter.p1_sap_dinucleotide_evidence_synthesis import (
    EVIDENCE_COLUMNS,
    GROUP_EVIDENCE_COLUMNS,
    SUMMARY_COLUMNS,
    build_p1_sap_dinucleotide_evidence_synthesis,
)


def _group(group_id, composition="C1", linkage="NORMAL_PHOSPHATE", assignments=1):
    return {
        "Dinucleotide_Group_ID": group_id,
        "Final_Elemental_Composition": composition,
        "Linkage_State": linkage,
        "Structural_Assignment_Count": assignments,
        "Group_Interpretation": "DINUCLEOTIDE_GROUP_SUPPORTED",
        "Source_Bond_Resolution_Status": "SOURCE_BOND_RESOLVED",
    }


def _assignment(group_id, number=1, bond="1-2"):
    return {
        "Structural_Assignment_ID": f"A{group_id}_{number}",
        "Dinucleotide_Group_ID": group_id,
        "Possible_Source_Bond": bond,
        "Source_Bond_Resolution_Status": "SOURCE_BOND_RESOLVED",
    }


def _feature(
    group_id, feature_id, physical_id="P1", *, qualified=True,
    candidate_specific=True, structure_specific=True,
):
    return {
        "Dinucleotide_Feature_ID": feature_id,
        "Dinucleotide_Group_ID": group_id,
        "Physical_Feature_ID": physical_id,
        "Feature_Eligible_For_Support": qualified,
        "Feature_Quality_Status": (
            "QUALIFIED_CHROMATOGRAPHIC_FEATURE" if qualified else "PROFILE_ONLY_REJECTED"
        ),
        "Mass_Accuracy_Support_Status": "MASS_COMPATIBLE",
        "Candidate_Specific": candidate_specific,
        "Linkage_Specific": candidate_specific,
        "Composition_Specific": candidate_specific,
        "Structure_Specific": structure_specific,
        "Source_Bond_Resolution_Status": "SOURCE_BOND_RESOLVED",
    }


def _competition(group_id, types="", *, candidate_specific=False):
    return {
        "Physical_Feature_ID": "P1",
        "Dinucleotide_Group_ID": group_id,
        "Competing_Dinucleotide_Group_Count": 0 if candidate_specific else 1,
        "Competition_Types": types,
        "Candidate_Specific": candidate_specific,
        "Linkage_Specific": candidate_specific,
        "Composition_Specific": candidate_specific,
        "Structure_Specific": candidate_specific,
    }


def _ms2(group_id, feature_id, *, applicable=False):
    return {
        "Dinucleotide_Group_ID": group_id,
        "Dinucleotide_Feature_ID": feature_id,
        "Precursor_Compatible_MS2_Spectrum_Count": 1,
        "MS2_Model_Applicable": applicable,
    }


def _run(groups, assignments, features, competition=None, isotopes=None, ms2=None):
    return build_p1_sap_dinucleotide_evidence_synthesis(
        groups, assignments, features, competition or [], isotopes or [], ms2 or [],
    )


def test_different_composition_competition():
    groups = [_group("G1", "C1"), _group("G2", "C2")]
    assignments = [_assignment("G1"), _assignment("G2")]
    features = [
        _feature("G1", "F1", candidate_specific=False, structure_specific=False),
        _feature("G2", "F2", candidate_specific=False, structure_specific=False),
    ]
    competition = [
        _competition("G1", "DIFFERENT_COMPOSITION_WITHIN_TOLERANCE"),
        _competition("G2", "DIFFERENT_COMPOSITION_WITHIN_TOLERANCE"),
    ]

    row = _run(groups, assignments, features, competition).evidence_rows[0]

    assert row["Composition_Resolution_Status"] == "COMPOSITION_UNRESOLVED"
    assert "COMPOSITION_COMPETITION" in row["Unresolved_Issues"]
    assert row["Targeted_MS2_Priority"] == "HIGH"


def test_normal_pt_competition():
    groups = [
        _group("G1", linkage="NORMAL_PHOSPHATE"),
        _group("G2", linkage="PHOSPHOROTHIOATE"),
    ]
    assignments = [_assignment("G1"), _assignment("G2")]
    features = [
        _feature("G1", "F1", candidate_specific=False, structure_specific=False),
        _feature("G2", "F2", candidate_specific=False, structure_specific=False),
    ]
    competition = [
        _competition("G1", "NORMAL_PHOSPHATE_VS_PT_COMPETITION"),
        _competition("G2", "NORMAL_PHOSPHATE_VS_PT_COMPETITION"),
    ]

    row = _run(groups, assignments, features, competition).evidence_rows[0]

    assert row["Linkage_Resolution_Status"] == "LINKAGE_UNRESOLVED"
    assert "LINKAGE_COMPETITION" in row["Unresolved_Issues"]


def test_same_composition_structural_isomerism():
    groups = [_group("G1", assignments=2)]
    assignments = [_assignment("G1", 1), _assignment("G1", 2, "2-3")]
    features = [_feature("G1", "F1", candidate_specific=True, structure_specific=False)]
    competition = [_competition(
        "G1", "SAME_COMPOSITION_STRUCTURAL_ISOMERS", candidate_specific=True,
    )]

    row = _run(groups, assignments, features, competition).evidence_rows[0]

    assert row["Composition_Resolution_Status"] == "COMPOSITION_RESOLVED"
    assert row["Structure_Resolution_Status"] == "STRUCTURE_UNRESOLVED"
    assert "STRUCTURAL_ISOMERISM" in row["Unresolved_Issues"]


def test_confounded_isotope_is_noninformative():
    groups = [_group("G1")]
    assignments = [_assignment("G1")]
    features = [_feature("G1", "F1")]
    isotopes = [{
        "Dinucleotide_Feature_ID": "F1",
        "Envelope_Status": "ENVELOPE_COMPATIBLE",
        "Envelope_Confounded": True,
    }]

    row = _run(groups, assignments, features, isotopes=isotopes).evidence_rows[0]

    assert row["Evidence_Level"] == "QUALIFIED_MS1_FEATURE"
    assert row["Isotope_Evidence_Status"] == "NONINFORMATIVE"
    assert "ISOTOPE_NOT_INFORMATIVE" in row["Unresolved_Issues"]
    assert "ISOTOPE_COMPATIBLE" not in row["Evidence_Basis"]


def test_no_ms2_is_high_priority_for_qualified_feature():
    result = _run([_group("G1")], [_assignment("G1")], [_feature("G1", "F1")])
    row = result.evidence_rows[0]

    assert row["MS2_Provenance_Status"] == "NO_PRECURSOR_COMPATIBLE_MS2"
    assert row["Targeted_MS2_Priority"] == "HIGH"


def test_ms2_presence_does_not_resolve_structure():
    groups = [_group("G1", assignments=2)]
    assignments = [_assignment("G1", 1), _assignment("G1", 2, "2-3")]
    features = [_feature("G1", "F1", structure_specific=False)]

    row = _run(
        groups, assignments, features, ms2=[_ms2("G1", "F1", applicable=False)],
    ).evidence_rows[0]

    assert row["Evidence_Level"] == "QUALIFIED_MS1_WITH_MS2_PROVENANCE"
    assert row["Structure_Resolution_Status"] == "STRUCTURE_UNRESOLVED"
    assert row["Source_Bond_Resolution_Status"] == "SOURCE_BOND_UNRESOLVED"
    assert row["Targeted_MS2_Priority"] == "MEDIUM"
    assert "MS2_FRAGMENT_MODEL_NOT_VALIDATED" in row["Unresolved_Issues"]


def test_no_usable_feature_is_not_applicable():
    group = _group("G1")
    group["Group_Interpretation"] = "NO_MATCH"
    result = _run([group], [_assignment("G1")], [])

    assert result.evidence_rows == []
    group_row = result.group_evidence_rows[0]
    assert group_row["Best_Evidence_Level"] == "NO_USABLE_EVIDENCE"
    assert group_row["Targeted_MS2_Priority"] == "NOT_APPLICABLE"
    assert group_row["Applied_To_Formal_Result"] is False
    assert group_row["Formal_Result_Changed"] is False


def test_input_order_is_deterministic_and_columns_are_stable():
    groups = [_group("G2", "C2"), _group("G1", "C1")]
    assignments = [_assignment("G2"), _assignment("G1")]
    features = [
        _feature("G2", "F2", candidate_specific=False, structure_specific=False),
        _feature("G1", "F1", candidate_specific=False, structure_specific=False),
    ]
    competition = [
        _competition("G2", "DIFFERENT_COMPOSITION_WITHIN_TOLERANCE"),
        _competition("G1", "DIFFERENT_COMPOSITION_WITHIN_TOLERANCE"),
    ]
    isotopes = [
        {"Dinucleotide_Feature_ID": "F2", "Envelope_Status": "NOT_ASSESSED"},
        {"Dinucleotide_Feature_ID": "F1", "Envelope_Status": "NOT_ASSESSED"},
    ]
    ms2 = [_ms2("G2", "F2"), _ms2("G1", "F1")]

    first = _run(groups, assignments, features, competition, isotopes, ms2)
    second = _run(
        list(reversed(deepcopy(groups))),
        list(reversed(deepcopy(assignments))),
        list(reversed(deepcopy(features))),
        list(reversed(deepcopy(competition))),
        list(reversed(deepcopy(isotopes))),
        list(reversed(deepcopy(ms2))),
    )

    assert first.to_jsonable() == second.to_jsonable()
    assert list(first.evidence_rows[0]) == EVIDENCE_COLUMNS
    assert all(list(row) == GROUP_EVIDENCE_COLUMNS for row in first.group_evidence_rows)
    assert list(first.summary_rows[0]) == SUMMARY_COLUMNS
    for rows in (
        first.evidence_rows, first.group_evidence_rows, first.summary_rows,
    ):
        assert all(row["Applied_To_Formal_Result"] is False for row in rows)
        assert all(row["Formal_Change_Ready"] is False for row in rows)
        assert all(row["Formal_Result_Changed"] is False for row in rows)


def test_raw_match_only_group_level():
    group = _group("G1")
    group["Group_Interpretation"] = "RAW_MATCH_ONLY"
    result = _run([group], [_assignment("G1")], [])
    assert result.group_evidence_rows[0]["Best_Evidence_Level"] == "RAW_MATCH_ONLY"


def test_empty_inputs():
    result = _run([], [], [])
    assert result.evidence_rows == []
    assert result.group_evidence_rows == []
    assert result.summary_rows[0]["Physical_Feature_Count"] == 0
    assert result.summary_rows[0]["Applied_To_Formal_Result"] is False
