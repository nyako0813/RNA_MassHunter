"""Non-propagating comparison of formal legacy and complete-structure candidates."""
from __future__ import annotations
from typing import Any

def classify_legacy_composite(*, mass_equal: bool, structure_equal: bool = False,
    isomeric: bool = False, legacy_parent: bool = False, ms2_discriminated: bool = False,
    legacy_present: bool = True, composite_present: bool = True) -> str:
    if not legacy_present: return "COMPOSITE_ONLY"
    if not composite_present: return "LEGACY_ONLY"
    if ms2_discriminated: return "MS2_DISCRIMINATED"
    if structure_equal: return "BOTH_EQUIVALENT"
    if mass_equal and isomeric: return "BOTH_ISOMERIC"
    if mass_equal: return "MASS_EQUIVALENT_STRUCTURE_EXCLUSIVE"
    if legacy_parent: return "LEGACY_PARENT_OF_COMPOSITE"
    return "COMPOSITE_REFINES_LEGACY"

def _position_state_ids(complete_structure_id: Any) -> set[str]:
    return {part for part in str(complete_structure_id or "").split("|") if "@" in part}

def _split_ids(value: Any) -> set[str]:
    return {item for item in str(value or "").split(";") if item}

def compare_legacy_composite(support_rows: list[dict[str, Any]], phase1_rows: list[dict[str, Any]],
    ranking_rows: list[dict[str, Any]], *, audit_level: str = "full") -> list[dict[str, Any]]:
    output = []
    for support in support_rows:
        cid = support["Candidate_ID"]
        positions = {int(x) for x in str(support.get("Modified_Positions") or "").split(";") if x}
        complete_structure_id = str(support.get("Complete_Structure_ID") or "")
        exact_phase = [r for r in phase1_rows
            if int(r.get("Position") or -1) in positions
            and str(r.get("Complete_Structure_ID") or "")
            and str(r.get("Complete_Structure_ID") or "") in complete_structure_id]
        equivalent_ids = sorted(set().union(*(_split_ids(r.get("Legacy_Equivalent_IDs")) for r in exact_phase))
            if exact_phase else set())
        parent_ids = sorted(set().union(*(_split_ids(r.get("Included_Component_IDs")) for r in exact_phase))
            if exact_phase else set())
        ranks = [r for r in ranking_rows if str(r.get("Modification_ID") or "") in equivalent_ids]
        mass_equal = bool(equivalent_ids)
        legacy_parent = bool(parent_ids) and not mass_equal
        ms2 = bool(support.get("MS2_Position_Informative_Count") or support.get("MS2_Backbone_Informative_Count"))
        classification = classify_legacy_composite(
            mass_equal=mass_equal, isomeric=any(bool(r.get("Is_Isomeric")) for r in exact_phase),
            legacy_parent=legacy_parent, ms2_discriminated=ms2,
            legacy_present=bool(equivalent_ids or parent_ids), composite_present=True,
        )
        output.append({
            "Candidate_ID": cid, "Legacy_Candidate_IDs": ";".join(equivalent_ids),
            "Legacy_Parent_IDs": ";".join(parent_ids),
            "Exact_Phase1_Candidate_IDs": ";".join(str(r.get("Candidate_ID") or "") for r in exact_phase),
            "Comparison_Class": classification, "Neutral_Mass_Equivalent": mass_equal,
            "Position_Compatible": bool(exact_phase), "Chemical_Exclusivity_Checked": False,
            "MS1_Support_Count": support.get("MS1_Unique_Support_Count", 0),
            "MS1_Nondiscriminating_Count": support.get("MS1_Nondiscriminating_Count", 0),
            "MS2_Localization_Support_Count": support.get("MS2_Position_Informative_Count", 0),
            "Legacy_Formal_Ranks": ";".join(str(r.get("Rank")) for r in ranks),
            "Audit_Level": audit_level, "Applied_To_Formal_Result": False, "Formal_Change_Ready": False,
        })
    return output
