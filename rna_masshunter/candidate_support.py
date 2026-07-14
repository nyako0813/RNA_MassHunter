"""Candidate-level aggregation of composite MS1/MS2 shadow evidence."""
from __future__ import annotations
from typing import Any

def aggregate_candidate_support(structures: list[Any], fragment_rows: list[dict[str, Any]],
    ms1_rows: list[dict[str, Any]], ms2_rows: list[dict[str, Any]],
    blocked_rows: list[dict[str, Any]], *, audit_level: str = "full") -> list[dict[str, Any]]:
    output = []
    for structure in structures:
        cid = structure.candidate_id
        fragments = [r for r in fragment_rows if r.get("Candidate_ID") == cid]
        ms1 = [r for r in ms1_rows if r.get("Candidate_ID") == cid]
        ms2 = [r for r in ms2_rows if r.get("Candidate_ID") == cid]
        blocked = [r for r in blocked_rows if r.get("Candidate_ID") == cid]
        matched_ids = {r["Fragment_ID"] for r in ms1 if r.get("Match_Status") == "matched"}
        observable_ids = {r["Fragment_ID"] for r in ms1 if r.get("Match_Status") != "not_observable"}
        physical_key = lambda r: (str(r.get("Observed_Scan") or ""), r.get("Observed_RT"), r.get("Observed_mz"))
        unique = len({physical_key(r) for r in ms1 if r.get("Support_Class") == "unique_composite_support"})
        shared = len({physical_key(r) for r in ms1 if str(r.get("Support_Class", "")).startswith("shared_")})
        nondiscriminating = len({physical_key(r) for r in ms1 if r.get("Support_Class") == "observation_nondiscriminating"})
        conflicts = sum(r.get("Support_Class") == "mass_match_but_structure_conflict" for r in ms1)
        coverage = len(matched_ids) / len(observable_ids) if observable_ids else 0.0
        positions = sorted(structure.position_states)
        bonds = sorted(k for k, v in structure.bond_states.items() if v.state != "normal_phosphate")
        output.append({
            "Candidate_ID": cid, "Complete_Structure_ID": next((r.get("Complete_Structure_ID") for r in fragments), cid),
            "Parent_Base": ";".join(structure.position_states[p].parent_base for p in positions),
            "Modified_Positions": ";".join(map(str, positions)), "Backbone_Modified_Bonds": ";".join(bonds),
            "Theoretical_Fragment_Count": len(fragments), "Observable_Fragment_Count": len(observable_ids),
            "MS1_Matched_Fragment_Count": len(matched_ids), "MS1_Unique_Support_Count": unique,
            "MS1_Shared_Support_Count": shared, "MS1_Nondiscriminating_Count": nondiscriminating,
            "MS1_Isomeric_Unresolved_Count": sum(r.get("Support_Class") == "isomeric_unresolved" for r in ms1),
            "MS2_Matched_Ion_Count": len(ms2),
            "MS2_Position_Informative_Count": sum(bool(r.get("Position_Informative")) for r in ms2),
            "MS2_Backbone_Informative_Count": sum(bool(r.get("Backbone_Informative")) for r in ms2),
            "Blocked_Cleavage_Match_Count": sum(r.get("Observed_mz") not in ("", None) for r in blocked),
            "Conflicting_Observation_Count": conflicts, "Support_Coverage": coverage,
            "Support_Status": "supported" if unique or ms2 or any(r.get("Observed_mz") not in ("", None) for r in blocked)
                              else "shared_only" if shared else "no_observation",
            "Audit_Level": audit_level, "Applied_To_Formal_Result": False, "Formal_Change_Ready": False,
        })
    return output