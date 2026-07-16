from copy import deepcopy

from rna_masshunter.audit_policy import AuditPolicy, included_sheet_names
from rna_masshunter.rnase_ms2_standard_composite_crosswalk import (
    build_rnase_ms2_standard_composite_crosswalk,
)


def standard(mod="M1", parent="F1", position=3, trna=12, mass=10.0):
    return {
        "Modification_ID": mod, "Parent_Fragment_ID": parent,
        "Candidate_Position_In_Parent": position, "Candidate_tRNA_Position": trna,
        "Mass_Shift": mass,
    }


def composite(candidate="C1", structure="S1", position=12, explicit="M1", mass_only="",
              mass=10.0, isomer=False):
    return {
        "Candidate_ID": candidate, "Complete_Structure_ID": structure,
        "Composite_Position": position, "Applied_Transform_IDs": "T1",
        "Explicit_Legacy_Modification_IDs": explicit,
        "Mass_Equivalent_Modification_IDs": mass_only,
        "Canonical_Structure_ID": f"U@{position}|state={structure}",
        "Exact_Mass_Delta": mass, "Is_Isomeric": isomer,
        "Isomer_Group_ID": "ISO1" if isomer else "",
    }


def fragment(parent="F1", start=10, end=20):
    return {"Fragment_ID": parent, "Start_Position": start, "End_Position": end}


def build(standards=(), composites=(), bonds=(), fragments=()):
    return build_rnase_ms2_standard_composite_crosswalk(
        list(standards), list(composites), list(bonds), list(fragments),
    )


def only(result):
    assert len(result.crosswalk_rows) == 1
    return result.crosswalk_rows[0]


def test_exact_match_and_absolute_position_conversion():
    row = only(build([standard()], [composite()], fragments=[fragment()]))
    assert row["Standard_Absolute_Sequence_Position"] == 12
    assert row["Position_Match_Status"] == "MATCH"
    assert row["Modification_Identity_Match_Status"] == "EXPLICIT_MATCH"
    assert row["Crosswalk_Status"] == "EXACT_MATCH"
    assert row["Crosswalk_Cardinality"] == "ONE_TO_ONE"


def test_mass_only_and_isobaric_different_id_never_exact():
    row = only(build([standard(mod="ISOBAR")], [composite(explicit="M1", mass_only="ISOBAR")], fragments=[fragment()]))
    assert row["Mass_Shift_Match_Status"] == "MATCH"
    assert row["Crosswalk_Status"] == "MASS_EQUIVALENT_ONLY"


def test_position_match_identity_unknown():
    row = only(build([standard()], [composite(explicit="", mass_only="")], fragments=[fragment()]))
    assert row["Crosswalk_Status"] == "POSITION_MATCH_IDENTITY_UNRESOLVED"


def test_position_conflict():
    row = only(build([standard(position=2, trna=11)], [composite(position=12)], fragments=[fragment()]))
    assert row["Parent_Fragment_Match_Status"] == "MATCH"
    assert row["Crosswalk_Status"] == "POSITION_CONFLICT"


def test_identity_conflict():
    row = only(build([standard(mod="M2")], [composite(explicit="M1", mass_only="", mass=20.0)], fragments=[fragment()]))
    assert row["Modification_Identity_Match_Status"] == "CONFLICT"
    assert row["Crosswalk_Status"] == "MODIFICATION_IDENTITY_CONFLICT"


def test_parent_fragment_conflict_is_not_ignored():
    row = only(build([standard(position=30, trna=39)], [composite(position=39)], fragments=[fragment()]))
    assert row["Position_Match_Status"] == "MATCH"
    assert row["Parent_Fragment_Match_Status"] == "CONFLICT"
    assert row["Crosswalk_Status"] == "PARENT_FRAGMENT_CONFLICT"


def test_structural_isomer_prevents_exact():
    row = only(build([standard()], [composite(isomer=True)], fragments=[fragment()]))
    assert row["Structural_Isomer_Sharing"] is True
    assert row["Crosswalk_Status"] == "POSITION_MATCH_IDENTITY_UNRESOLVED"


def test_one_to_many():
    result = build([standard()], [composite(), composite(candidate="C2", structure="S2")], fragments=[fragment()])
    assert {row["Crosswalk_Cardinality"] for row in result.crosswalk_rows} == {"ONE_TO_MANY"}


def test_many_to_one():
    standards = [standard(trna=12), standard(trna=112)]
    result = build(standards, [composite()], fragments=[fragment()])
    assert {row["Crosswalk_Cardinality"] for row in result.crosswalk_rows} == {"MANY_TO_ONE"}


def test_many_to_many():
    standards = [standard(mod="M1"), standard(mod="M2")]
    composites = [
        composite(candidate="C1", structure="S1", explicit="M1;M2"),
        composite(candidate="C2", structure="S2", explicit="M1;M2"),
    ]
    result = build(standards, composites, fragments=[fragment()])
    assert {row["Crosswalk_Cardinality"] for row in result.crosswalk_rows} == {"MANY_TO_MANY"}


def test_overlapping_parent_fragments_use_exact_parent_id():
    standards = [standard(parent="F1"), standard(parent="F2")]
    fragments = [fragment("F1", 10, 20), fragment("F2", 10, 20)]
    result = build(standards, [composite()], fragments=fragments)
    assert len(result.crosswalk_rows) == 2
    assert all(row["Parent_Fragment_Match_Status"] == "MATCH" for row in result.crosswalk_rows)


def test_backbone_rows_are_not_crosswalk_candidates():
    bond = {"Candidate_ID": "B1", "Complete_Structure_ID": "BS1", "Bond_ID": "11_12",
            "Applied_Backbone_Transform_IDs": "PT"}
    result = build([standard()], [], bonds=[bond], fragments=[fragment()])
    assert result.crosswalk_rows == []
    assert result.summary_rows[0]["Composite_Position_Component_Count"] == 0


def test_insufficient_provenance_and_empty_input():
    row = only(build([standard(parent="MISSING")], [composite()], fragments=[]))
    assert row["Crosswalk_Status"] == "INSUFFICIENT_PROVENANCE"
    empty = build()
    assert empty.crosswalk_rows == [] and empty.summary_rows[0]["Crosswalk_Row_Count"] == 0


def test_input_order_determinism():
    standards = [standard(mod="M2"), standard(mod="M1")]
    composites = [composite(candidate="C2", structure="S2", explicit="M1;M2"), composite(explicit="M1;M2")]
    first = build(standards, composites, fragments=[fragment()])
    second = build(list(reversed(deepcopy(standards))), list(reversed(deepcopy(composites))), fragments=[fragment()])
    assert first == second


def test_all_formal_flags_false_and_inputs_unchanged():
    standards = [standard()]; composites = [composite()]
    before = deepcopy((standards, composites))
    result = build(standards, composites, fragments=[fragment()])
    for rows in (result.summary_rows, result.crosswalk_rows):
        for row in rows:
            assert row["Applied_To_Formal_Result"] is False
            assert row["Formal_Change_Ready"] is False
            assert row["Formal_Result_Changed"] is False
    assert (standards, composites) == before


def test_sheet_inclusion():
    summary = "RNase_MS2_Standard_Composite_Su"
    detail = "RNase_MS2_Standard_Composite_Cr"
    names = [summary, detail]
    standard_names, _ = included_sheet_names(names, AuditPolicy.from_level("standard"))
    audit_names, _ = included_sheet_names(names, AuditPolicy.from_level("audit"))
    full_names, _ = included_sheet_names(names, AuditPolicy.from_level("full"))
    assert standard_names == []
    assert audit_names == [summary]
    assert full_names == names
