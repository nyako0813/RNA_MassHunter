from types import SimpleNamespace

from rna_masshunter.biological_context import biological_context_priority_rows, score_biological_context
from rna_masshunter.evidence_ranking import build_modification_evidence_ranking
from rna_masshunter.models import Modification


def _context_config(priority=True):
    return SimpleNamespace(
        organism={"group": "archaea", "species": "Example archaeon"},
        sequence={"type": "RNA"},
        biological_context={
            "enabled": True, "organism_group": "", "organism_species": "",
            "trna_name": "", "trna_type": "", "anticodon": "",
            "focus_positions": [34], "focus_position_window": 2,
            "priority_modifications": ["mA"] if priority else [],
            "priority_keywords": ["methylation"] if priority else [],
            "boost": {"priority_modification": 1.5, "priority_keyword_match": 0.75,
                      "focus_position_match": 1.0, "focus_position_nearby": 0.5,
                      "organism_rule_supported": 1.0, "pathway_supported": 1.0,
                      "trna_context_supported": 1.0},
            "penalties": {"organism_context_conflict": -2.0},
        },
        modification_evidence_ranking={
            "enabled": True, "use_biological_context": True,
            "cap_context_only_confidence": "Medium",
            "require_ms_evidence_for_context_boosted_high": True,
            "use_known_modification_candidates": True,
            "use_ms1_fragment_evidence": True, "use_ms2_precursor_evidence": True,
            "use_ms2_modified_ion_evidence": True, "use_localization_evidence": True,
            "use_organism_rules": True, "use_trna_context": True,
            "require_ms2_evidence_for_high_confidence": True,
            "require_position_discriminating_ions_for_localization_confidence": True,
            "min_final_score_to_report": 0, "max_ranked_candidates": 100,
            "weights": {"known_modification_candidate": 1.0},
        },
    )


def test_empty_and_user_configured_context_scoring():
    modification = Modification("mA", "mA", 14.0, "biological", ["A"], raw={
        "name": "methylated adenosine", "chemical_group": "methylation",
    })
    empty = _context_config(priority=False)
    empty.biological_context["focus_positions"] = []
    result = score_biological_context({"Modification_ID": "mA", "Candidate_tRNA_Position": 34}, modification, empty)
    assert result["Biological_Context_Score"] == 0

    configured = _context_config(priority=True)
    exact = score_biological_context({
        "Modification_ID": "mA", "Modification_Name": "methylated adenosine",
        "Chemical_Group": "methylation", "Candidate_tRNA_Position": 34,
    }, modification, configured)
    assert exact["Biological_Context_Score"] == 3.25
    assert exact["Context_Matched_Priority_Modification"] is True
    assert exact["Context_Matched_Keywords"] == "methylation"
    assert exact["Context_Focus_Position_Match"] == "exact"

    nearby = score_biological_context({
        "Modification_ID": "other", "Modification_Name": "other",
        "Candidate_tRNA_Position": 35,
    }, modification, configured)
    assert nearby["Context_Focus_Position_Match"] == "nearby"
    assert nearby["Context_Focus_Position_Distance"] == 1
    assert biological_context_priority_rows(configured)


def test_context_only_candidate_is_not_high_confidence():
    config = _context_config(priority=True)
    modification = Modification("mA", "mA", 14.0, "biological", ["A"], raw={
        "name": "methylated adenosine", "chemical_group": "methylation",
    })
    rows, summary = build_modification_evidence_ranking(
        config=config, modifications=[modification], theoretical_fragments=[], fragment_ms1_matches=[],
        known_candidates=[{"modification_id": "mA", "source_type": "intact", "source_id": "I1", "priority_score": 1.0}],
        ms2_results={}, rule_set={}, pathways=[],
    )
    assert rows[0]["Biological_Context_Score"] > 0
    assert rows[0]["Final_Confidence"] not in {"High", "Very High"}
    assert summary[0]["Candidates_With_Biological_Context_Support"] == 1


if __name__ == "__main__":
    test_empty_and_user_configured_context_scoring()
    test_context_only_candidate_is_not_high_confidence()
    print("synthetic biological context tests: OK")
