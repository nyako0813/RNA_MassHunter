from copy import deepcopy
from types import SimpleNamespace

from rna_masshunter.biological_position_prior import (
    BIOLOGICAL_PLAUSIBILITY_COLUMNS,
    POSITION_PRIOR_COLUMNS,
    evaluate_biological_position_priors,
)
from rna_masshunter.models import Modification
from rna_masshunter.review_dashboard import _build_top_candidates


def config(landmarks=None, enabled=True):
    return SimpleNamespace(
        sequence={"sequence": "A" * 33 + "U" + "A" * 5 + "A" + "C" * 40, "wobble_position": 34},
        organism={"group": "archaea"},
        ms2_annotation={"biological_position_prior": {
            "enabled": enabled, "apply_to_final_score": False,
            "position_numbering": {"landmarks": landmarks or {}},
        }},
    )


def rules():
    return {"version": "test", "families": [
        {"id": "wobble_U_family", "modification_ids": ["mnm5U", "s2U", "s4U"], "canonical_landmark": "wobble", "chemically_possible_distance": 3},
        {"id": "t6A_family", "modification_ids": ["t6A"], "canonical_landmark": "anticodon_adjacent_3prime", "chemically_possible_distance": 3},
        {"id": "preQ1_family", "modification_ids": ["preQ1"], "canonical_landmark": "queuosine_site", "chemically_possible_distance": 2, "organism_warnings": {"archaea": "review archaeal pathway context"}},
    ]}


def mods():
    return [
        Modification("mnm5U", "mnm5U", 43.0, "x", ["U"]),
        Modification("s2U", "s2U", 15.9772, "x", ["U"]),
        Modification("s4U", "s4U", 15.9772, "x", ["U"]),
        Modification("t6A", "t6A", 145.0, "x", ["A"]),
        Modification("preQ1", "preQ1", 28.0, "x", ["G"]),
    ]


def row(mod="mnm5U", pos=34, base="U", score=7.0, rank=1, parent="F1"):
    return {
        "Rank": rank, "Final_Score": score, "Final_Confidence": "High",
        "Modification_ID": mod, "Modification_Name": mod, "Parent_Fragment_ID": parent,
        "Candidate_tRNA_Position": pos, "Candidate_Base": base, "Candidate_Position_In_Parent": 1,
        "Parent_Start": pos, "Localization_Level": "Strong", "Position_Discriminating_Evidence": True,
        "Biological_Context_Score": 0.0, "Has_MS2_Precursor_Evidence": True,
        "Has_Modified_Ion_Evidence": True, "Position_Ambiguity_Status": "single_candidate_position",
    }


def evaluate(rows, cfg=None):
    return evaluate_biological_position_priors(cfg or config(), rows, mods(), rules())


def test_canonical_and_adjacent_and_inconsistent_classes():
    _, positions, _, _ = evaluate([row(pos=34), row(pos=35, rank=2), row(pos=60, rank=3)])
    assert [x["Position_Class"] for x in positions] == ["canonical_position", "adjacent_to_canonical", "biologically_inconsistent"]


def test_no_landmark_means_unknown_not_fixed_sprinzl():
    cfg = config(); cfg.sequence.pop("wobble_position")
    _, positions, _, _ = evaluate([row()], cfg)
    assert positions[0]["Position_Class"] == "unknown"


def test_explicit_landmark_map_is_input_sequence_based():
    _, positions, _, _ = evaluate([row(mod="t6A", pos=40, base="A")], config({"anticodon_adjacent_3prime": 40}))
    assert positions[0]["Position_Class"] == "canonical_position"
    assert positions[0]["Position_Numbering_System"] == "input_sequence_1_based"


def test_parent_base_compatible_incompatible_ambiguous_unknown():
    _, _, rows, _ = evaluate([row(base="U"), row(base="A", rank=2), row(base="N", rank=3), row(base="", pos="", rank=4)])
    assert [x["Parent_Base_Compatibility"] for x in rows] == ["compatible", "incompatible", "ambiguous", "unknown"]


def test_sequence_supplies_missing_candidate_base():
    _, _, rows, _ = evaluate([row(base="")])
    assert rows[0]["Candidate_Base"] == "U"
    assert rows[0]["Parent_Base_Compatibility"] == "compatible"


def test_ms2_localization_is_separate_from_biology():
    enriched, _, rows, _ = evaluate([row(pos=60)])
    assert rows[0]["MS2_Localization_Evidence"] == "Strong"
    assert rows[0]["Position_Class"] == "biologically_inconsistent"
    assert enriched[0]["Final_Confidence"] == "High"


def test_structural_isobars_reported_at_same_candidate_position():
    _, _, rows, _ = evaluate([row(mod="s2U")])
    assert rows[0]["Structure_Ambiguity_Status"] == "position_resolved_structure_unresolved"
    assert "s4U" in rows[0]["Alternative_Structural_Candidates"]
    assert rows[0]["Structure_Discriminating_Evidence"] is False


def test_shadow_does_not_mutate_final_fields_or_membership_order():
    original = [row(score=3, rank=2), row(mod="t6A", pos=40, base="A", score=9, rank=1)]
    before = deepcopy(original)
    enriched, _, _, diagnostics = evaluate(original, config({"anticodon_adjacent_3prime": 40}))
    assert original == before
    assert [(x["Rank"], x["Final_Score"], x["Final_Confidence"], x["Modification_ID"]) for x in enriched] == [(x["Rank"], x["Final_Score"], x["Final_Confidence"], x["Modification_ID"]) for x in before]
    assert diagnostics[0]["Rows_Changed"] == 0
    assert all(x["Shadow_Only"] for x in enriched)


def test_unknown_family_and_position_are_safe():
    enriched, positions, rows, _ = evaluate([row(mod="unlisted", pos="", base="")])
    assert enriched[0]["Position_Class"] == "unknown"
    assert positions[0]["Modification_Family"] == "unclassified"
    assert rows[0]["Biological_Plausibility_Level"] == "neutral"


def test_empty_candidates_safe_and_columns_contract():
    enriched, positions, rows, diagnostics = evaluate([])
    assert enriched == positions == rows == []
    assert diagnostics[0]["Total_Candidates"] == 0
    assert len("Modification_Position_Priors") <= 31
    assert len("MS2_Biological_Plausibility") <= 31
    assert "Position_Class" in POSITION_PRIOR_COLUMNS
    assert "Shadow_Final_Score" in BIOLOGICAL_PLAUSIBILITY_COLUMNS


def test_disabled_returns_unchanged_copy():
    source = [row()]
    enriched, positions, rows, diagnostics = evaluate(source, config(enabled=False))
    assert enriched == source and enriched is not source
    assert positions == rows == []
    assert diagnostics[0]["Evaluated_Candidates"] == 0


def test_top_candidates_get_shadow_columns_without_priority_or_order_change():
    source = [row(score=9, rank=1), row(mod="t6A", pos=40, base="A", score=8, rank=2, parent="F2")]
    plain = _build_top_candidates(__import__('pandas').DataFrame(source), __import__('pandas').DataFrame(), {"max_top_candidates": 50})
    enriched, _, _, _ = evaluate(source, config({"anticodon_adjacent_3prime": 40}))
    shadow = _build_top_candidates(__import__('pandas').DataFrame(enriched), __import__('pandas').DataFrame(), {"max_top_candidates": 50})
    assert list(plain["Modification_ID"]) == list(shadow["Modification_ID"])
    assert list(plain["Review_Rank"]) == list(shadow["Review_Rank"])
    assert list(plain["Review_Priority"]) == list(shadow["Review_Priority"])
    assert "Shadow_Final_Score" in shadow.columns


def test_preQ1_archaea_warning_can_be_expressed_as_unknown_without_landmark():
    _, positions, _, _ = evaluate([row(mod="preQ1", base="G")])
    assert positions[0]["Position_Class"] == "unknown"
    assert "not configured" in positions[0]["Position_Prior_Reason"]
    assert "Organism warning" in positions[0]["Position_Prior_Reason"]
    assert positions[0]["Position_Prior_Level"] == "warning"


def test_position_falls_back_to_parent_start_plus_local_offset():
    candidate = row(); candidate["Candidate_tRNA_Position"] = ""; candidate["Parent_Start"] = 33; candidate["Candidate_Position_In_Parent"] = 2
    enriched, _, _, _ = evaluate([candidate])
    assert enriched[0]["Candidate_tRNA_Position"] == 34


def test_canonical_position_has_high_prior_score():
    enriched, _, _, _ = evaluate([row()])
    assert enriched[0]["Position_Prior_Level"] == "high"
    assert enriched[0]["Position_Prior_Score"] == 2.0


def test_adjacent_position_has_moderate_prior_score():
    enriched, _, _, _ = evaluate([row(pos=35)])
    assert enriched[0]["Position_Prior_Level"] == "moderate"
    assert enriched[0]["Position_Prior_Score"] == 1.0


def test_incompatible_parent_base_penalizes_shadow_only():
    enriched, _, _, _ = evaluate([row(base="A", pos=60)])
    assert enriched[0]["Parent_Base_Prior_Score"] == -2.0
    assert enriched[0]["Final_Score"] == 7.0
    assert enriched[0]["Shadow_Final_Score"] == 3.0


def test_shadow_confidence_copies_final_confidence():
    enriched, _, _, _ = evaluate([row()])
    assert enriched[0]["Shadow_Final_Confidence"] == enriched[0]["Final_Confidence"]


def test_unclassified_family_has_no_position_reward():
    enriched, _, _, _ = evaluate([row(mod="other", pos=34, base="U")])
    assert enriched[0]["Modification_Family"] == "unclassified"
    assert enriched[0]["Position_Prior_Score"] == 0.0


def test_structural_alternatives_filter_incompatible_parent_bases():
    _, _, rows, _ = evaluate([row(mod="s2U", base="U")])
    assert "t6A" not in str(rows[0]["Alternative_Structural_Candidates"])


def test_nonisobaric_candidate_has_no_structural_alternative():
    _, _, rows, _ = evaluate([row(mod="mnm5U")])
    assert rows[0]["Structure_Ambiguity_Status"] == "no_structural_alternative_identified"


def test_position_discrimination_changes_structure_status_not_alternatives():
    candidate = row(mod="s2U"); candidate["Position_Discriminating_Evidence"] = False
    _, _, rows, _ = evaluate([candidate])
    assert rows[0]["Structure_Ambiguity_Status"] == "position_and_structure_unresolved"
    assert "s4U" in rows[0]["Alternative_Structural_Candidates"]


def test_diagnostics_count_position_classes():
    _, _, _, diagnostics = evaluate([row(), row(pos=35, rank=2), row(pos=60, rank=3)])
    assert diagnostics[0]["Canonical_Position"] == 1
    assert diagnostics[0]["Adjacent_To_Canonical"] == 1
    assert diagnostics[0]["Biologically_Inconsistent"] == 1


def test_diagnostics_count_parent_compatibility():
    _, _, _, diagnostics = evaluate([row(base="U"), row(base="A", rank=2), row(base="N", rank=3)])
    assert diagnostics[0]["Compatible_Parent_Base"] == 1
    assert diagnostics[0]["Incompatible_Parent_Base"] == 1
    assert diagnostics[0]["Ambiguous_Parent_Base"] == 1


def test_plausibility_sheet_rows_follow_ranking_order():
    source = [row(rank=7), row(mod="t6A", pos=40, base="A", rank=2)]
    _, _, rows, _ = evaluate(source, config({"anticodon_adjacent_3prime": 40}))
    assert [item["Rank"] for item in rows] == [7, 2]


def test_top_best_final_score_is_not_shadow_score():
    enriched, _, _, _ = evaluate([row(score=7.0)])
    top = _build_top_candidates(__import__('pandas').DataFrame(enriched), __import__('pandas').DataFrame(), {"max_top_candidates": 50})
    assert top.iloc[0]["Best_Final_Score"] == 7.0
    assert top.iloc[0]["Shadow_Final_Score"] == 10.0


def test_apply_to_final_score_flag_still_reports_zero_mutated_rows():
    cfg = config(); cfg.ms2_annotation["biological_position_prior"]["apply_to_final_score"] = True
    enriched, _, _, diagnostics = evaluate([row()], cfg)
    assert enriched[0]["Final_Score"] == 7.0
    assert diagnostics[0]["Rows_Changed"] == 0
