from pathlib import Path
from types import SimpleNamespace
import json

import pytest

from rna_masshunter.audit_policy import AuditPolicy, included_sheet_names
from rna_masshunter.config import load_config
from rna_masshunter.p1_sap_dinucleotide_candidates import (
    CONDENSATION_ADJUSTMENT, DINUCLEOTIDE_MODEL_VERSION, FORMAL_FALSE,
    LINKAGE_COMPOSITIONS, LOCALIZATION_FALSE, NUCLEOSIDE_COMPOSITIONS, PT_DELTA,
    adjacent_bonds, build_position_states, dinucleotide_settings,
    extract_target_candidates, generate_dinucleotide_candidates,
)
from rna_masshunter.p1_sap_dinucleotide_interpretation import (
    build_p1_sap_dinucleotide_audit, build_target_results,
)

ROOT = Path(__file__).resolve().parent


def cfg(*, max_mods=1, charges=(1,), polarity="positive", strict=None, targets=None, search=(100, 1000), organism="generic"):
    return SimpleNamespace(
        p1_sap_dinucleotide={
            "enabled": True,
            "candidate_generation": {
                "max_modifications_per_side": max_mods,
                "max_composite_states_per_position": 64,
                "max_candidate_count": 100000,
                "include_normal_phosphate": True,
                "include_phosphorothioate": True,
                "charges": list(charges), "polarity": polarity,
                "strict_positions": strict or {},
            },
            "search": {"mz_min": search[0], "mz_max": search[1], "tolerance_ppm": 10},
            "mass_accuracy": {"strong_ppm": 2, "moderate_ppm": 5, "search_ppm": 10},
            "feature_quality": {"min_spectrum_count": 2, "min_profile_point_count": 2},
            "isotope": {"enabled": True, "tolerance_ppm": 20, "require_same_scan": True},
            "ms2_provenance": {"enabled": True}, "targets": targets or [],
        },
        sequence={"wobble_position": 34}, organism={"species": organism},
        instrument={}, p1_annotation={},
    )


def assignments(result, left=None, right=None, linkage=None):
    return [row for row in result.assignments
            if (left is None or row["Left_State_ID"] == left)
            and (right is None or row["Right_State_ID"] == right)
            and (linkage is None or row["Linkage_State"] == linkage)]


@pytest.mark.parametrize("sequence,expected", [("", 0), ("A", 0), ("AG", 1), ("AGCU", 3), ("AGCUAGCU", 7)])
def test_arbitrary_sequence_lengths_generate_n_minus_one_bonds(sequence, expected):
    assert len(adjacent_bonds(sequence)) == expected


@pytest.mark.parametrize("sequence,pairs", [("AG", [("A", "G")]), ("CU", [("C", "U")]), ("GCAU", [("G", "C"), ("C", "A"), ("A", "U")])])
def test_all_standard_base_pairs_preserve_direction(sequence, pairs):
    assert [(left, right) for _i, _j, left, right in adjacent_bonds(sequence)] == pairs


def test_no_position_constraint_allows_transform_at_matching_base():
    states, _ = build_position_states("GG", ROOT, config=cfg())
    assert "m2G" in {state.state_id for state in states[1]}
    assert "m2G" in {state.state_id for state in states[2]}


def test_strict_position_constraint_excludes_other_positions():
    states, rejected = build_position_states("GG", ROOT, config=cfg(strict={"m2G": [2]}))
    assert "m2G" not in {state.state_id for state in states[1]}
    assert "m2G" in {state.state_id for state in states[2]}
    assert any(row["Reason"] == "position_disallowed" for row in rejected)


def test_unmatched_tRNA_does_not_inherit_specific_position_hypotheses():
    states, _ = build_position_states("GG", ROOT, config=cfg(organism="unrelated_species"))
    assert "m22G" in {state.state_id for state in states[1]}


def test_organism_context_is_generic_and_not_hardcoded():
    left = generate_dinucleotide_candidates("AG", ROOT, config=cfg(organism="species_a"))
    right = generate_dinucleotide_candidates("AG", ROOT, config=cfg(organism="species_b"))
    assert {(row["Final_Elemental_Composition"], row["Linkage_State"]) for row in left.candidates} == {(row["Final_Elemental_Composition"], row["Linkage_State"]) for row in right.candidates}


def test_left_right_and_composite_modifications_are_generated():
    result = generate_dinucleotide_candidates("GU", ROOT, config=cfg(max_mods=2))
    assert assignments(result, "G", "U", "PHOSPHOROTHIOATE")
    assert assignments(result, "m2G", "U", "PHOSPHOROTHIOATE")
    assert assignments(result, "G", "Um", "PHOSPHOROTHIOATE")
    assert assignments(result, "m2G", "Um", "PHOSPHOROTHIOATE")
    assert any("+" in row["Left_State_ID"] or "+" in row["Right_State_ID"] for row in result.assignments)


def test_normal_phosphate_and_pt_are_separate_groups_with_defined_delta():
    result = generate_dinucleotide_candidates("AG", ROOT, config=cfg(max_mods=0))
    normal = assignments(result, "A", "G", "NORMAL_PHOSPHATE")[0]
    pt = assignments(result, "A", "G", "PHOSPHOROTHIOATE")[0]
    assert normal["Linkage_Composition"] == LINKAGE_COMPOSITIONS["NORMAL_PHOSPHATE"].canonical_string()
    assert pt["Linkage_Composition"] == LINKAGE_COMPOSITIONS["PHOSPHOROTHIOATE"].canonical_string()
    assert pt["Neutral_Mass"] - normal["Neutral_Mass"] == pytest.approx(PT_DELTA.exact_mass)


def test_full_composition_uses_two_nucleosides_linkage_and_condensation():
    result = generate_dinucleotide_candidates("AG", ROOT, config=cfg(max_mods=0))
    row = assignments(result, "A", "G", "NORMAL_PHOSPHATE")[0]
    expected = NUCLEOSIDE_COMPOSITIONS["A"] + NUCLEOSIDE_COMPOSITIONS["G"] + LINKAGE_COMPOSITIONS["NORMAL_PHOSPHATE"] + CONDENSATION_ADJUSTMENT
    assert row["Final_Elemental_Composition"] == expected.canonical_string()
    assert row["Neutral_Mass"] == pytest.approx(expected.exact_mass)


def test_group_has_generic_schema_and_no_confirmed_representative():
    group = generate_dinucleotide_candidates("GUG", ROOT, config=cfg(max_mods=0)).candidates[0]
    required = {"Dinucleotide_Group_ID", "Model_Version", "Final_Elemental_Composition", "Neutral_Mass", "Theoretical_mz", "Charge", "Polarity", "Linkage_State", "Observable", "Search_Enabled", "Structural_Assignment_Count", "Possible_Source_Bonds", "Position_Constraint_Summary", "Chemical_Constraint_Summary"}
    assert required <= group.keys()
    assert group["Representative_Is_Confirmed"] is False
    assert group["Structure_Resolution_Status"] == "STRUCTURE_UNRESOLVED"


def test_assignments_are_separate_rows_with_generic_schema():
    result = generate_dinucleotide_candidates("AG", ROOT, config=cfg(max_mods=0))
    row = result.assignments[0]
    assert row["Structural_Assignment_ID"].startswith(row["Dinucleotide_Group_ID"])
    assert row["Possible_Source_Bond"] == "1-2"
    assert row["Chemical_Constraint_Status"] == "COMPATIBLE"
    assert row["Sequence_Position_Localized"] is False


def test_direction_isomers_and_source_bonds_group_without_loss():
    result = generate_dinucleotide_candidates("GUG", ROOT, config=cfg(max_mods=0))
    group = next(row for row in result.candidates if row["Possible_Source_Bonds"] == "1-2;2-3" and row["Linkage_State"] == "NORMAL_PHOSPHATE")
    assert group["Structural_Assignment_Count"] == 2
    assert "1:G-NORMAL_PHOSPHATE-2:U" in group["Possible_Position_Assignments"]
    assert "2:U-NORMAL_PHOSPHATE-3:G" in group["Possible_Position_Assignments"]


def test_observable_ranges_do_not_use_intact_reconstruction_range():
    config = cfg(max_mods=0, search=(100, 700)); config.reconstruction = {"mz_min": 2000, "mz_max": 3000}
    groups = generate_dinucleotide_candidates("AG", ROOT, config=config).candidates
    assert groups and all(group["Observable_In_Dinucleotide_Search"] for group in groups)
    assert all(group["Observable"] for group in groups)


def test_outside_search_group_is_retained_but_not_searched():
    groups = generate_dinucleotide_candidates("AG", ROOT, config=cfg(max_mods=0, search=(100, 200))).candidates
    assert groups and all(not group["Observable_In_Dinucleotide_Search"] for group in groups)
    assert all(group["Search_Enabled"] is False and "OUTSIDE_DINUCLEOTIDE_SEARCH_RANGE" in group["Not_Observable_Reason"] for group in groups)


def test_z1_z2_and_polarity_are_generic():
    positive = generate_dinucleotide_candidates("AG", ROOT, config=cfg(max_mods=0), charges=(1, 2), polarity="positive")
    negative = generate_dinucleotide_candidates("AG", ROOT, config=cfg(max_mods=0), charges=(1,), polarity="negative")
    assert {row["Charge"] for row in positive.candidates} == {1, 2}
    assert all(row["Polarity"] == "negative" for row in negative.candidates)
    assert positive.candidates[0]["Theoretical_mz"] != negative.candidates[0]["Theoretical_mz"]


def test_candidate_limit_reports_truncation():
    config = cfg(max_mods=2); config.p1_sap_dinucleotide["candidate_generation"]["max_candidate_count"] = 3
    result = generate_dinucleotide_candidates("GUG", ROOT, config=config)
    assert result.truncated and result.summary["Candidate_Generation_Reason"] == "max_candidate_count=3"


def test_target_absence_does_not_prevent_all_group_analysis():
    result = build_p1_sap_dinucleotide_audit(ROOT, "AG", [], cfg(max_mods=0, targets=[]), audit_level="audit")
    assert result["generated"].candidates
    assert result["sheets"]["P1_SAP_Dinuc_Targets"] == []
    assert result["summary"]["Searched_Group_Count"] == len(result["generated"].candidates)


def test_target_list_is_optional_multi_label_post_filter_only():
    base = generate_dinucleotide_candidates("AG", ROOT, config=cfg(max_mods=0))
    values = [base.candidates[0]["Theoretical_mz"], base.candidates[-1]["Theoretical_mz"]]
    targets = [{"label": "anything_a", "theoretical_mz": values[0], "tolerance_ppm": 1}, {"label": "anything_b", "theoretical_mz": values[1], "tolerance_ppm": 1}]
    before = json.dumps(base.candidates, sort_keys=True)
    rows = build_target_results(targets, base.candidates, [])
    assert [row["Target_Label"] for row in rows] == ["anything_a", "anything_b"]
    assert all(row["Matched_Group_Count"] == 1 for row in rows)
    assert json.dumps(base.candidates, sort_keys=True) == before


def test_legacy_flat_config_is_backward_compatible():
    config = SimpleNamespace(p1_sap_dinucleotide={"max_modifications_per_side": 2, "candidate_mass_tolerance_ppm": 15, "polarity": "negative"}, instrument={}, p1_annotation={})
    settings = dinucleotide_settings(config)
    assert settings["max_modifications_per_side"] == 2
    assert settings["search_tolerance_ppm"] == 15
    assert settings["polarity"] == "negative"


def test_fixture_mass_uses_generic_generation_and_preserves_assignments():
    config = load_config(ROOT / "config.yaml")
    result = generate_dinucleotide_candidates(config.sequence["sequence"], ROOT, config=config)
    groups = extract_target_candidates(result.candidates, 634.13269, 634.13272)
    assert len(groups) == 1
    group = groups[0]
    assert group["Structural_Assignment_Count"] == 85
    assert group["Possible_Source_Bonds"] == "7-8;10-11;11-12;12-13;15-16;55-56"
    assert "10:m22G-PHOSPHOROTHIOATE-11:U" in group["Possible_Position_Assignments"]


def test_generic_sheet_policy_and_excel_lengths():
    names = ["P1_SAP_Dinuc_Summary", "P1_SAP_Dinuc_Groups", "P1_SAP_Dinuc_Assignments", "P1_SAP_Dinuc_SpecPeaks", "P1_SAP_Dinuc_Features", "P1_SAP_Dinuc_Isotopes", "P1_SAP_Dinuc_Competition", "P1_SAP_Dinuc_MS2", "P1_SAP_Dinuc_Targets"]
    standard, _ = included_sheet_names(names, AuditPolicy.from_level("standard")); audit, _ = included_sheet_names(names, AuditPolicy.from_level("audit")); full, _ = included_sheet_names(names, AuditPolicy.from_level("full"))
    assert standard == []
    assert set(audit) == {names[0], names[1], names[-1]}
    assert set(full) == set(names)
    assert max(map(len, names)) <= 31


def test_json_generic_object_and_formal_flags():
    result = build_p1_sap_dinucleotide_audit(ROOT, "AG", [], cfg(max_mods=0), audit_level="full")
    payload = result["payload"]; json.dumps(payload)
    assert set(payload) == {"dinucleotide_audit"}
    assert "target_results" in payload["dinucleotide_audit"]
    assert not any("634" in key or "m22" in key for key in payload["dinucleotide_audit"])
    for rows in result["sheets"].values():
        for row in rows:
            assert all(row[key] is False for key in FORMAL_FALSE)
            if "Sequence_Position_Localized" in row:
                assert row["Sequence_Position_Localized"] is LOCALIZATION_FALSE["Sequence_Position_Localized"]


def test_disabled_config_produces_no_shadow_sheets():
    config = cfg(max_mods=0); config.p1_sap_dinucleotide["enabled"] = False
    result = build_p1_sap_dinucleotide_audit(ROOT, "AG", [], config, audit_level="full")
    assert result["sheets"] == {}
    assert result["payload"]["dinucleotide_audit"]["status"] == "DISABLED_BY_CONFIG"


def test_excel_column_schemas_are_unique():
    from rna_masshunter.p1_sap_dinucleotide_candidates import GROUP_COLUMNS, ASSIGNMENT_COLUMNS, SUMMARY_COLUMNS
    from rna_masshunter.p1_sap_dinucleotide_feature_audit import SPECPEAK_COLUMNS, FEATURE_COLUMNS, ISOTOPE_COLUMNS, COMPETITION_COLUMNS, MS2_COLUMNS
    from rna_masshunter.p1_sap_dinucleotide_interpretation import TARGET_COLUMNS
    for columns in (GROUP_COLUMNS, ASSIGNMENT_COLUMNS, SUMMARY_COLUMNS, SPECPEAK_COLUMNS, FEATURE_COLUMNS, ISOTOPE_COLUMNS, COMPETITION_COLUMNS, MS2_COLUMNS, TARGET_COLUMNS):
        assert len(columns) == len(set(columns))
