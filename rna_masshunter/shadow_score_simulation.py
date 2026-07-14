"""What-if score/rank simulation that never mutates formal rows."""
from __future__ import annotations
from typing import Any

def simulate_shadow_scores(support_rows: list[dict[str, Any]], comparison_rows: list[dict[str, Any]],
    formal_ranking: list[dict[str, Any]], *, audit_level: str = "full") -> list[dict[str, Any]]:
    output = []
    for support in support_rows:
        cid = support["Candidate_ID"]
        comparison = next((r for r in comparison_rows if r["Candidate_ID"] == cid), {})
        legacy_ids = set(filter(None, str(comparison.get("Legacy_Candidate_IDs") or "").split(";")))
        formal = next((r for r in formal_ranking if str(r.get("Modification_ID") or "") in legacy_ids), None)
        positive = (float(support.get("MS1_Unique_Support_Count") or 0) * 0.5
                    + float(support.get("MS2_Position_Informative_Count") or 0) * 1.0
                    + float(support.get("MS2_Backbone_Informative_Count") or 0) * 1.0
                    + float(support.get("Blocked_Cleavage_Match_Count") or 0) * 1.0)
        penalty = float(support.get("Conflicting_Observation_Count") or 0) * 1.0
        has_formal_baseline = formal is not None
        formal_score = float(formal.get("Final_Score") or 0.0) if formal else ""
        simulated = formal_score + positive - penalty if formal else ""
        formal_rank = formal.get("Rank", "") if formal else ""
        if formal:
            competing = [row for row in formal_ranking if row is not formal]
            simulated_rank = 1 + sum(float(row.get("Final_Score") or 0) >= simulated for row in competing)
        else:
            simulated_rank = ""
        output.append({
            "Candidate_ID": cid, "Baseline_Status": "exact_legacy_equivalent" if formal else "no_exact_legacy_equivalent",
            "Legacy_Formal_Score": formal_score,
            "Composite_Shadow_Support": positive, "Composite_Shadow_Penalty": penalty,
            "Simulated_Shadow_Score": simulated, "Formal_Rank": formal_rank,
            "Simulated_Shadow_Rank": simulated_rank,
            "Score_Delta": (simulated-formal_score) if formal else "",
            "Rank_Delta": (simulated_rank-int(formal_rank)) if formal_rank not in ("", None) else "",
            "Would_Change_Top_Candidate": bool(formal) and simulated_rank == 1 and formal_rank != 1,
            "Would_Change_Confidence": bool(formal) and positive-penalty >= 2,
            "Would_Change_Formal_Result": False, "Audit_Level": audit_level,
            "Applied_To_Formal_Result": False, "Formal_Change_Ready": False,
        })
    return output
