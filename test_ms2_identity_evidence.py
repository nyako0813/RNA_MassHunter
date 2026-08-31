from copy import deepcopy

import pandas as pd

from rna_masshunter.ms2_identity_evidence import IDENTITY_COLUMNS, IDENTITY_SHADOW_COLUMNS, build_ms2_modification_identity
from rna_masshunter.review_dashboard import _build_top_candidates


def candidate(mod="mA", pos=1, rank=1, score=7.0):
    return {
        "Rank": rank, "Final_Score": score, "Final_Confidence": "High",
        "Modification_ID": mod, "Modification_Name": mod, "Parent_Fragment_ID": "F1",
        "Candidate_tRNA_Position": 10 + pos - 1, "Candidate_Position_In_Parent": pos,
        "Candidate_Base": "A", "Has_MS2_Precursor_Evidence": True,
        "Has_Known_Modification_Candidate": True, "Has_MS1_Fragment_Evidence": True,
        "Has_Localization_Evidence": True, "Position_Discriminating_Evidence": False,
        "Num_Position_Discriminating_Ions": 0, "Ambiguity_Group_ID": "G1",
        "Structural_Isomer_Group_ID": "", "Structure_Ambiguity_Status": "no_structural_alternative_identified",
        "Alternative_Structural_Candidates": "", "Structure_Discriminating_Evidence": False,
        "Biological_Plausibility_Level": "high", "Position_Class": "canonical_position",
        "Biological_Context_Score": 0.0, "Position_Ambiguity_Status": "ambiguous",
    }


def match(mod="mA", pos=1, spectrum="S1", ion="d2", observed=500.0, theoretical=500.001):
    return {
        "Spectrum_ID": spectrum, "Scan_Index": 5, "Parent_Fragment_ID": "F1",
        "Modification_ID": mod, "Candidate_Modification_Position_In_Parent": pos,
        "Ion_Contains_Modification": True, "Observed_mz": observed, "Observed_Intensity": 1000.0,
        "Ion_ID": ion, "Ion_Type": ion[0], "Theoretical_mz": theoretical, "Mass_Error_ppm": -2.0,
    }


def groups(positions="10;11", status="unknown"):
    return [{"Ambiguity_Group_ID": "G1", "Candidate_Positions_In_tRNA": positions, "Position_Ambiguity_Status": status}]


def run(rows=None, matches=None, ambiguity=None, enabled=True, localization=None):
    return build_ms2_modification_identity(
        rows if rows is not None else [candidate()], matches or [], localization or [],
        ambiguity if ambiguity is not None else groups(), enabled=enabled,
    )


def test_modified_fragment_evidence_present_and_absent():
    enriched, _ = run(matches=[match()])
    assert enriched[0]["Has_Modified_Fragment_Ion_Evidence"] is True
    enriched, _ = run(matches=[])
    assert enriched[0]["Has_Modified_Fragment_Ion_Evidence"] is False
    assert enriched[0]["MS2_Identity_Evidence_Level"] == "fragment_mass_shift_supported"


def test_duplicate_match_is_counted_once():
    item = match()
    enriched, _ = run(matches=[item, dict(item)])
    assert enriched[0]["Modified_Fragment_Match_Count"] == 1
    assert enriched[0]["Unique_Modified_Fragment_Ion_Count"] == 1


def test_unique_ion_count_and_series():
    enriched, _ = run(matches=[match(ion="d2"), match(ion="w3", spectrum="S2", observed=600)])
    assert enriched[0]["Unique_Modified_Fragment_Ion_Count"] == 2
    assert enriched[0]["Modified_Fragment_Ion_Series"] == "d;w"


def test_same_observed_peak_multiple_assignments_warns():
    enriched, _ = run(matches=[match(ion="d2"), match(ion="w3", theoretical=500.002)])
    assert "multiple theoretical ion assignments" in enriched[0]["MS2_Identity_Warnings"]


def test_position_evidence_is_separate_from_modified_fragment_identity():
    row = candidate(); row["Position_Discriminating_Evidence"] = True
    row["Num_Position_Discriminating_Ions"] = 2
    enriched, _ = run(rows=[row], matches=[])
    assert enriched[0]["Position_Localization_Status"] == "adjacent_positions_ambiguous"
    assert enriched[0]["Has_Modified_Fragment_Ion_Evidence"] is False
    assert enriched[0]["MS2_Identity_Evidence_Level"] == "fragment_mass_shift_supported"


def test_adjacent_and_multiple_position_ambiguity():
    enriched, _ = run(matches=[match()], ambiguity=groups("10;11"))
    assert enriched[0]["Position_Localization_Status"] == "adjacent_positions_ambiguous"
    enriched, _ = run(matches=[match()], ambiguity=groups("10;13"))
    assert enriched[0]["Position_Localization_Status"] == "multiple_positions_ambiguous"


def test_position_localized_without_isomer_group():
    row = candidate(); row["Position_Discriminating_Evidence"] = True
    localization = [{"Modification_ID": "mA", "Parent_Fragment_ID": "F1", "Candidate_Modification_Position_In_Parent": 1, "Num_Position_Discriminating_Modified_Ions": 1, "Num_Unmodified_Counterpart_Matches": 1}]
    enriched, _ = run(rows=[row], matches=[match()], ambiguity=groups("10;11", "resolved_by_discriminating_ions"), localization=localization)
    assert enriched[0]["MS2_Identity_Evidence_Level"] == "position_localized"
    assert enriched[0]["Shadow_MS2_Identity_Confidence"] == "Moderate"
    assert enriched[0]["Structure_Resolution_Status"] == "no_structural_alternative"


def test_position_localized_structure_unresolved_caps_confidence():
    row = candidate(); row.update({
        "Position_Discriminating_Evidence": True, "Structural_Isomer_Group_ID": "SIG_X",
        "Structure_Ambiguity_Status": "position_resolved_structure_unresolved",
        "Alternative_Structural_Candidates": "mB",
    })
    localization = [{"Modification_ID": "mA", "Parent_Fragment_ID": "F1", "Candidate_Modification_Position_In_Parent": 1, "Num_Position_Discriminating_Modified_Ions": 1, "Num_Unmodified_Counterpart_Matches": 1}]
    enriched, _ = run(rows=[row], matches=[match()], ambiguity=groups("10;11", "resolved_by_discriminating_ions"), localization=localization)
    assert enriched[0]["Position_Localization_Status"] == "position_resolved"
    assert enriched[0]["Structure_Resolution_Status"] == "position_resolved_structure_unresolved"
    assert enriched[0]["MS2_Identity_Evidence_Level"] == "structure_isomer_unresolved"
    assert enriched[0]["Shadow_MS2_Identity_Confidence"] == "Low"


def test_isomer_group_without_structure_evidence_is_unresolved():
    row = candidate(); row.update({
        "Structural_Isomer_Group_ID": "SIG_X", "Structure_Ambiguity_Status": "position_and_structure_unresolved",
        "Alternative_Structural_Candidates": "mB;mC",
    })
    enriched, _ = run(rows=[row], matches=[match()])
    assert enriched[0]["Structure_Resolution_Status"] == "position_and_structure_unresolved"
    assert "do not distinguish structural isomers" in enriched[0]["MS2_Identity_Warnings"]


def test_structure_resolved_requires_structure_discriminating_evidence():
    row = candidate(); row.update({
        "Structural_Isomer_Group_ID": "SIG_X", "Alternative_Structural_Candidates": "mB",
        "Structure_Discriminating_Evidence": True,
    })
    row["Position_Discriminating_Evidence"] = True
    localization = [{"Modification_ID": "mA", "Parent_Fragment_ID": "F1", "Candidate_Modification_Position_In_Parent": 1, "Num_Position_Discriminating_Modified_Ions": 1, "Num_Unmodified_Counterpart_Matches": 1}]
    enriched, _ = run(rows=[row], matches=[match()], ambiguity=groups("10;11", "resolved_by_discriminating_ions"), localization=localization)
    assert enriched[0]["Structure_Resolution_Status"] == "structure_resolved"
    assert enriched[0]["MS2_Identity_Evidence_Level"] == "structure_isomer_resolved"
    assert enriched[0]["Shadow_MS2_Identity_Confidence"] == "Moderate"


def test_canonical_biology_alone_does_not_make_identity_high():
    enriched, _ = run(matches=[])
    assert enriched[0]["Biological_Plausibility_Level"] == "high"
    assert enriched[0]["Shadow_MS2_Identity_Confidence"] == "Low"


def test_zero_candidate_and_ms2_disabled_are_safe():
    assert run(rows=[])[0] == []
    enriched, identity = run(matches=[match()], enabled=False)
    assert len(enriched) == len(identity) == 1
    assert enriched[0]["MS2_Identity_Evidence_Level"] == "unsupported"
    assert "MS2 annotation disabled" in enriched[0]["MS2_Identity_Warnings"]


def test_biological_prior_and_existing_rank_fields_are_unchanged():
    source = [candidate()]
    before = deepcopy(source)
    enriched, _ = run(rows=source, matches=[match()])
    for field in ["Rank", "Final_Score", "Final_Confidence", "Position_Class", "Biological_Plausibility_Level"]:
        assert enriched[0][field] == before[0][field]
    assert source == before


def test_top_candidate_count_order_and_existing_values_are_unchanged():
    source = [candidate("mA", 1, 1, 9), candidate("mB", 2, 2, 8)]
    plain = _build_top_candidates(pd.DataFrame(source), pd.DataFrame(), {"max_top_candidates": 50})
    enriched, _ = run(rows=source, matches=[match("mA", 1), match("mB", 2)])
    shadow = _build_top_candidates(pd.DataFrame(enriched), pd.DataFrame(), {"max_top_candidates": 50})
    old_columns = [column for column in plain.columns if column not in IDENTITY_SHADOW_COLUMNS]
    assert len(plain) == len(shadow) == 2
    assert plain[old_columns].equals(shadow[old_columns])
    assert "MS2_Identity_Evidence_Level" in shadow.columns


def test_localization_rows_are_reused_for_discriminating_count():
    localization = [{
        "Modification_ID": "mA", "Parent_Fragment_ID": "F1",
        "Candidate_Modification_Position_In_Parent": 1,
        "Num_Position_Discriminating_Modified_Ions": 3,
    }]
    enriched, _ = run(matches=[match()], localization=localization)
    assert enriched[0]["Position_Discriminating_Ion_Count"] == 3


def test_sheet_name_and_required_columns_contract():
    assert len("MS2_Modification_Identity") <= 31
    required = {
        "Rank", "Modification_ID", "Has_Modified_Fragment_Ion_Evidence",
        "Position_Localization_Status", "Structure_Resolution_Status",
        "MS2_Identity_Evidence_Level", "Shadow_MS2_Identity_Confidence",
    }
    assert required <= set(IDENTITY_COLUMNS)


def test_physical_key_ignores_theoretical_ion_id():
    from rna_masshunter.ms2_identity_evidence import physical_observed_peak_key
    first, second = match(ion="d2"), match(ion="w9", theoretical=501.0)
    first["RT"] = second["RT"] = 1.234567
    assert physical_observed_peak_key(first) == physical_observed_peak_key(second)


def test_different_observed_peaks_have_different_physical_keys():
    from rna_masshunter.ms2_identity_evidence import physical_observed_peak_key
    first, second = match(observed=500.0), match(observed=500.0001)
    assert physical_observed_peak_key(first) != physical_observed_peak_key(second)


def test_cross_candidate_physical_assignment_is_detected_and_propagated():
    rows = [candidate("mA", 1, 1), candidate("mB", 2, 2)]
    enriched, _, assignments = build_ms2_modification_identity(
        rows, [match("mA", 1), match("mB", 2)], [], groups(), return_assignments=True,
    )
    assert len({item["Physical_Observed_Peak_Key"] for item in assignments}) == 1
    assert all(item["Physical_Peak_Candidate_Count"] == 2 for item in assignments)
    assert all(item["Physical_Peak_Shared_Across_Candidates"] for item in assignments)
    assert all(item["Has_Cross_Candidate_Peak_Sharing"] for item in enriched)


def test_structural_isomer_shared_peak_is_group_level_not_individual():
    rows = [candidate("mA", 1, 1), candidate("mB", 1, 2)]
    for row in rows:
        row["Structural_Isomer_Group_ID"] = "SIG_X"
        row["Alternative_Structural_Candidates"] = "mB" if row["Modification_ID"] == "mA" else "mA"
        row["Structure_Ambiguity_Status"] = "position_and_structure_unresolved"
    enriched, _, assignments = build_ms2_modification_identity(
        rows, [match("mA", 1), match("mB", 1)], [], groups(), return_assignments=True,
    )
    assert all(item["Physical_Peak_Shared_Across_Isomers"] for item in assignments)
    assert all(item["Identity_Evidence_Scope"] == "structural_isomer_group_level" for item in enriched)
    assert all(item["Candidate_Specific_Evidence_Peak_Count"] == 0 for item in enriched)
    assert all(item["Shadow_MS2_Identity_Confidence"] == "Low" for item in enriched)


def test_candidate_specific_peak_counts_for_individual_identity():
    enriched, _, assignments = build_ms2_modification_identity(
        [candidate()], [match()], [], groups(), return_assignments=True,
    )
    assert enriched[0]["Candidate_Specific_Evidence_Peak_Count"] == 1
    assert enriched[0]["Unique_Physical_Peak_Count"] == 1
    assert assignments[0]["Counts_For_Individual_Identity"] is True


def test_shared_peak_does_not_receive_candidate_specific_score_bonus():
    specific, _ = run(matches=[match()])
    rows = [candidate("mA", 1, 1), candidate("mB", 2, 2)]
    shared, _ = build_ms2_modification_identity(rows, [match("mA", 1), match("mB", 2)], [], groups())
    assert specific[0]["Shadow_MS2_Identity_Score"] > shared[0]["Shadow_MS2_Identity_Score"]


def test_unrelated_candidate_shared_peak_is_cross_candidate_ambiguous():
    rows = [candidate("mA", 1, 1), candidate("mB", 2, 2)]
    enriched, _, assignments = build_ms2_modification_identity(
        rows, [match("mA", 1), match("mB", 2)], [], groups(), return_assignments=True,
    )
    assert all(item["Evidence_Scope"] == "cross_candidate_ambiguous" for item in assignments)
    assert all(item["Cross_Candidate_Ambiguous_Peak_Count"] == 1 for item in enriched)
    assert all(item["Shadow_MS2_Identity_Confidence"] == "Low" for item in enriched)


def test_group_unresolved_prevents_candidate_position_resolved():
    row = candidate(); row["Position_Discriminating_Evidence"] = True; row["Num_Position_Discriminating_Ions"] = 1
    enriched, _ = run(rows=[row], matches=[match()], ambiguity=groups("10;12", "ambiguous"))
    assert enriched[0]["Candidate_Position_Resolution_Status"] == "multiple_positions_ambiguous"
    assert enriched[0]["Position_Resolution_Ceiling_Applied"] is True


def test_group_partially_resolved_caps_candidate_status():
    row = candidate(); row["Position_Discriminating_Evidence"] = True; row["Num_Position_Discriminating_Ions"] = 1
    enriched, _ = run(rows=[row], matches=[match()], ambiguity=groups("10;12", "partially_resolved"))
    assert enriched[0]["Candidate_Position_Resolution_Status"] == "partially_resolved"
    assert enriched[0]["Position_Resolution_Ceiling_Applied"] is True
    assert enriched[0]["Shadow_MS2_Identity_Confidence"] == "Low"


def test_single_modified_ion_without_counterpart_has_caveat():
    row = candidate(); row["Position_Discriminating_Evidence"] = True; row["Num_Position_Discriminating_Ions"] = 1
    localization = [{"Modification_ID": "mA", "Parent_Fragment_ID": "F1", "Candidate_Modification_Position_In_Parent": 1, "Num_Position_Discriminating_Modified_Ions": 1, "Num_Unmodified_Counterpart_Matches": 0}]
    enriched, _ = run(rows=[row], matches=[match()], ambiguity=groups("10;12", "resolved_by_discriminating_ions"), localization=localization)
    assert enriched[0]["Candidate_Position_Resolution_Status"] == "position_supported_single_ion"
    assert "single modified ion" in enriched[0]["Position_Resolution_Caveat"]
    assert enriched[0]["Shadow_MS2_Identity_Confidence"] == "Low"


def test_counterpart_evidence_allows_position_resolved_and_moderate():
    row = candidate(); row["Position_Discriminating_Evidence"] = True; row["Num_Position_Discriminating_Ions"] = 1
    localization = [{"Modification_ID": "mA", "Parent_Fragment_ID": "F1", "Candidate_Modification_Position_In_Parent": 1, "Num_Position_Discriminating_Modified_Ions": 1, "Num_Unmodified_Counterpart_Matches": 1}]
    enriched, _ = run(rows=[row], matches=[match()], ambiguity=groups("10;12", "resolved_by_discriminating_ions"), localization=localization)
    assert enriched[0]["Candidate_Position_Resolution_Status"] == "position_resolved"
    assert enriched[0]["Shadow_MS2_Identity_Confidence"] == "Moderate"


def test_identity_shadow_columns_do_not_duplicate_structural_group_column():
    from rna_masshunter.biological_position_prior import SHADOW_RANKING_COLUMNS
    from rna_masshunter.ms2_identity_evidence import IDENTITY_SHADOW_COLUMNS
    combined = SHADOW_RANKING_COLUMNS + IDENTITY_SHADOW_COLUMNS
    assert len(combined) == len(set(combined))
    assert combined.count("Structural_Isomer_Group_ID") == 1


def test_conflicting_canonical_structural_group_value_raises():
    import pytest
    from rna_masshunter.ms2_identity_evidence import _merge_shadow
    row = {"Structural_Isomer_Group_ID": "SIG_A"}
    with pytest.raises(ValueError, match="Conflicting Structural_Isomer_Group_ID"):
        _merge_shadow(row, {"Structural_Isomer_Group_ID": "SIG_B"})


def test_cnm5u_like_no_match_remains_low_and_biology_unchanged():
    row = candidate("cnm5U"); row["Position_Class"] = "canonical_position"; row["Biological_Plausibility_Level"] = "moderate"
    enriched, _ = run(rows=[row], matches=[])
    assert enriched[0]["Shadow_MS2_Identity_Confidence"] == "Low"
    assert enriched[0]["Position_Class"] == "canonical_position"
    assert enriched[0]["Biological_Plausibility_Level"] == "moderate"


def test_new_sheet_names_are_within_excel_limit():
    assert len("MS2_Modification_Identity") <= 31
    assert len("MS2_Identity_Peak_Assignments") <= 31


def test_isomer_sharing_is_retained_when_unrelated_candidate_also_uses_peak():
    rows = [candidate("mA", 1, 1), candidate("mB", 1, 2), candidate("other", 2, 3)]
    for row in rows[:2]:
        row["Structural_Isomer_Group_ID"] = "SIG_X"
        row["Alternative_Structural_Candidates"] = "mB" if row["Modification_ID"] == "mA" else "mA"
    enriched, _, assignments = build_ms2_modification_identity(
        rows, [match("mA", 1), match("mB", 1), match("other", 2)], [], groups(), return_assignments=True,
    )
    isomers = [item for item in enriched if item["Modification_ID"] in {"mA", "mB"}]
    assert all(item["Isomer_Group_Shared_Peak_Count"] == 1 for item in isomers)
    assert all(item["Cross_Candidate_Ambiguous_Peak_Count"] == 1 for item in isomers)
    assert all(item["Identity_Evidence_Scope"] == "cross_candidate_ambiguous" for item in isomers)
    assert all(item["Physical_Peak_Shared_Across_Isomers"] for item in assignments)
