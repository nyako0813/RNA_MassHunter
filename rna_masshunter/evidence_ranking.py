"""Integrate modification evidence into review-oriented rankings (MVP-5.4)."""

from dataclasses import asdict, is_dataclass
from typing import Any


RANKING_COLUMNS = [
    "Rank", "Final_Score", "Final_Confidence", "Final_Interpretation",
    "Modification_ID", "Modification_Name", "Modification_Category", "Target_Base", "Mass_Shift", "Is_Isobaric",
    "Candidate_tRNA_Position", "Candidate_Base", "Parent_Fragment_ID", "Parent_Sequence", "Parent_Start", "Parent_End", "Candidate_Position_In_Parent",
    "Has_MS1_Fragment_Evidence", "MS1_Fragment_Best_Confidence", "MS1_Fragment_Total_Intensity",
    "Has_Known_Modification_Candidate", "Known_Modification_Priority_Score",
    "Has_MS2_Precursor_Evidence", "Num_MS2_Precursor_Candidates", "Best_Precursor_Error_ppm", "Modified_Precursor_Rescue",
    "Has_Modified_Ion_Evidence", "Num_Modified_Ion_Matches", "Num_Informative_Modified_Ion_Matches", "Best_Modified_Ion_Error_ppm",
    "Has_Localization_Evidence", "Localization_Level", "Localization_Score", "Localization_Interpretation", "Num_c_Modified_Ions", "Num_y_Modified_Ions",
    "Organism_Group", "Organism_Species", "Rule_Set", "Organism_Rule_Supported", "TRNA_Context_Supported", "Context_Notes",
    "Ambiguous_Position", "Low_Information_Evidence", "Evidence_Warnings", "Notes",
]

SUMMARY_COLUMNS = [
    "Total_Ranked_Candidates", "Very_High", "High", "Medium", "Low", "Very_Low",
    "Candidates_With_MS2_Precursor_Evidence", "Candidates_With_Modified_Ion_Evidence",
    "Candidates_With_Localization_Evidence", "Ambiguous_Candidates", "Top_Modification_IDs", "Notes",
]


def _raw(item: Any) -> dict[str, Any]:
    return asdict(item) if is_dataclass(item) else dict(item)


def _bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _confidence_rank(value: Any) -> int:
    return {"high": 3, "medium": 2, "low": 1}.get(str(value or "").lower(), 0)


def _rule_modification_ids(rule_set: dict[str, Any] | None) -> set[str]:
    ids: set[str] = set()
    def visit(value: Any) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                if str(key).lower() in {"modification_id", "modification", "mod_id"} and isinstance(child, (str, int)):
                    ids.add(str(child))
                visit(child)
        elif isinstance(value, list):
            for child in value:
                visit(child)
    visit(rule_set or {})
    return ids


def build_modification_evidence_ranking(
    config: Any,
    modifications: list[Any],
    theoretical_fragments: list[Any],
    fragment_ms1_matches: list[Any],
    known_candidates: list[Any],
    ms2_results: dict[str, Any],
    rule_set: dict[str, Any] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    ranking = getattr(config, "modification_evidence_ranking", {}) or {}
    if not _bool(ranking.get("enabled"), True):
        return [], []
    weights = ranking.get("weights", {}) or {}
    weight = lambda name, default: _float(weights.get(name), default)
    modification_lookup = {str(getattr(item, "id", "")): item for item in modifications or []}
    fragment_lookup = {str(getattr(item, "fragment_id", "")): item for item in theoretical_fragments or []}
    rule_ids = _rule_modification_ids(rule_set)
    localization = list(ms2_results.get("MS2_Modification_Localization_Evidence", []) or [])
    precursors = list(ms2_results.get("MS2_Modified_Precursor_Candidates", []) or [])
    ion_matches = list(ms2_results.get("MS2_Modified_Ion_Matches", []) or [])

    candidates: dict[tuple[str, str, Any], dict[str, Any]] = {}
    def ensure(mod_id: Any, fragment_id: Any = "", position: Any = "") -> dict[str, Any]:
        key = (str(mod_id or ""), str(fragment_id or ""), position if position not in {None, ""} else "")
        return candidates.setdefault(key, {"mod_id": key[0], "fragment_id": key[1], "position": key[2]})

    for row in localization:
        item = ensure(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), row.get("Candidate_Modification_Position_In_Parent"))
        if not item.get("localization") or _float(row.get("Localization_Score")) > _float(item["localization"].get("Localization_Score")):
            item["localization"] = row
    for row in precursors:
        ensure(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), "").setdefault("precursors", []).append(row)
    for raw_candidate in known_candidates or []:
        row = _raw(raw_candidate)
        fragment_id = row.get("source_id") if str(row.get("source_type", "")).lower() == "fragment" else ""
        ensure(row.get("modification_id"), fragment_id, "").setdefault("known", []).append(row)

    # Attach non-position evidence to every matching localized position; retain non-localized candidates too.
    for key, item in list(candidates.items()):
        if item["position"] == "":
            continue
        base = candidates.get((item["mod_id"], item["fragment_id"], ""), {})
        item.setdefault("precursors", list(base.get("precursors", [])))
        item.setdefault("known", list(base.get("known", [])))

    output = []
    for item in candidates.values():
        mod_id, fragment_id, position = item["mod_id"], item["fragment_id"], item["position"]
        modification = modification_lookup.get(mod_id)
        mod_raw = getattr(modification, "raw", {}) or {}
        loc = item.get("localization", {})
        precursor_rows = item.get("precursors", [])
        known_rows = item.get("known", [])
        fragment = fragment_lookup.get(fragment_id)
        matching_ms1 = [_raw(match) for match in fragment_ms1_matches or [] if str(_raw(match).get("fragment_id") or "") == fragment_id]
        matching_ions = [row for row in ion_matches if str(row.get("Modification_ID") or "") == mod_id and str(row.get("Parent_Fragment_ID") or "") == fragment_id and (position == "" or row.get("Candidate_Modification_Position_In_Parent") == position)]
        modified_ions = [row for row in matching_ions if row.get("Ion_Contains_Modification")]
        informative_ions = [row for row in modified_ions if row.get("Informative_Ion")]
        level = str(loc.get("Localization_Level") or "None")
        ambiguous = loc.get("Localization_Interpretation") == "ambiguous-multiple-positions"
        low_information = bool(modified_ions) and not informative_ions
        shift = _float(getattr(modification, "mass_shift_from_unmodified", None), _float(precursor_rows[0].get("Modification_Mass_Shift") if precursor_rows else 0))
        isobaric = abs(shift) <= 1e-6
        score = 0.0
        if matching_ms1 and _bool(ranking.get("use_ms1_fragment_evidence"), True): score += weight("ms1_fragment_match", 1.0)
        if known_rows and _bool(ranking.get("use_known_modification_candidates"), True): score += weight("known_modification_candidate", 1.0)
        if precursor_rows and _bool(ranking.get("use_ms2_precursor_evidence"), True): score += weight("ms2_precursor_rescue", 2.0)
        if modified_ions and _bool(ranking.get("use_ms2_modified_ion_evidence"), True): score += weight("ms2_modified_ion_match", 2.0)
        if _bool(ranking.get("use_localization_evidence"), True): score += weight(f"localization_{level.lower()}", {"Weak": 1.0, "Moderate": 3.0, "Strong": 5.0}.get(level, 0.0))
        rule_supported = mod_id in rule_ids if rule_ids else False
        if rule_supported and _bool(ranking.get("use_organism_rules"), True): score += weight("organism_rule_supported", 1.5)
        trna_position = loc.get("Candidate_Modification_Position_In_tRNA", "")
        trna_supported = bool(trna_position and (getattr(config, "sequence", {}) or {}).get("type", "").upper() == "RNA")
        if trna_supported and _bool(ranking.get("use_trna_context"), True): score += weight("trna_context_supported", 1.0)
        if low_information: score += weight("low_information_penalty", -1.0)
        if ambiguous: score += weight("ambiguous_position_penalty", -1.0)
        if isobaric and precursor_rows and not modified_ions: score += weight("isobaric_precursor_penalty", -2.0)
        has_ms2 = bool(precursor_rows or modified_ions)
        confidence = _final_confidence(score, bool(precursor_rows), bool(modified_ions), level, bool(known_rows), has_ms2, isobaric, ranking)
        interpretation = _interpretation(bool(precursor_rows), bool(modified_ions), level, bool(matching_ms1), bool(known_rows), ambiguous, low_information)
        best_ms1 = max(matching_ms1, key=lambda row: _confidence_rank(row.get("confidence")), default={})
        best_precursor_error = min((abs(_float(row.get("Precursor_Error_ppm"))) for row in precursor_rows), default="")
        best_ion_error = min((abs(_float(row.get("Mass_Error_ppm"))) for row in modified_ions), default="")
        target_bases = getattr(modification, "target_bases", []) if modification else []
        output.append({
            "Final_Score": score, "Final_Confidence": confidence, "Final_Interpretation": interpretation,
            "Modification_ID": mod_id, "Modification_Name": mod_raw.get("name") or getattr(modification, "symbol", None) or (precursor_rows[0].get("Modification_Name") if precursor_rows else mod_id),
            "Modification_Category": getattr(modification, "category", ""), "Target_Base": ",".join(target_bases), "Mass_Shift": shift, "Is_Isobaric": isobaric,
            "Candidate_tRNA_Position": trna_position, "Candidate_Base": loc.get("Candidate_Modification_Base", ""), "Parent_Fragment_ID": fragment_id,
            "Parent_Sequence": loc.get("Parent_Sequence") or getattr(fragment, "sequence", ""), "Parent_Start": loc.get("Parent_Start", getattr(fragment, "start", "")),
            "Parent_End": loc.get("Parent_End", getattr(fragment, "end", "")), "Candidate_Position_In_Parent": position,
            "Has_MS1_Fragment_Evidence": bool(matching_ms1), "MS1_Fragment_Best_Confidence": best_ms1.get("confidence", ""),
            "MS1_Fragment_Total_Intensity": sum(_float(row.get("intensity")) for row in matching_ms1),
            "Has_Known_Modification_Candidate": bool(known_rows), "Known_Modification_Priority_Score": max((_float(row.get("priority_score")) for row in known_rows), default=""),
            "Has_MS2_Precursor_Evidence": bool(precursor_rows), "Num_MS2_Precursor_Candidates": len(precursor_rows), "Best_Precursor_Error_ppm": best_precursor_error,
            "Modified_Precursor_Rescue": any(_bool(row.get("Modified_Precursor_Rescue")) for row in precursor_rows),
            "Has_Modified_Ion_Evidence": bool(modified_ions), "Num_Modified_Ion_Matches": len(modified_ions),
            "Num_Informative_Modified_Ion_Matches": len(informative_ions), "Best_Modified_Ion_Error_ppm": best_ion_error,
            "Has_Localization_Evidence": bool(loc), "Localization_Level": level, "Localization_Score": _float(loc.get("Localization_Score")),
            "Localization_Interpretation": loc.get("Localization_Interpretation", ""), "Num_c_Modified_Ions": loc.get("Num_c_Modified_Ions", 0), "Num_y_Modified_Ions": loc.get("Num_y_Modified_Ions", 0),
            "Organism_Group": getattr(config, "organism", {}).get("group", ""), "Organism_Species": getattr(config, "organism", {}).get("species", ""),
            "Rule_Set": getattr(config, "organism", {}).get("rule_set", ""), "Organism_Rule_Supported": rule_supported,
            "TRNA_Context_Supported": trna_supported, "Context_Notes": "Rule/context points are applied only when explicit loaded data supports them.",
            "Ambiguous_Position": ambiguous, "Low_Information_Evidence": low_information,
            "Evidence_Warnings": "; ".join(filter(None, ["isobaric precursor evidence is non-specific" if isobaric else "", "ambiguous localization" if ambiguous else "", "1 nt/low-information ion evidence" if low_information else ""])),
            "Notes": "Ranking prioritizes candidates; it does not confirm modification identity or position.",
        })
    minimum = _float(ranking.get("min_final_score_to_report"), 0.0)
    output = [row for row in output if row["Final_Score"] >= minimum]
    output.sort(key=lambda row: (-row["Final_Score"], row["Modification_ID"], row["Parent_Fragment_ID"], str(row["Candidate_Position_In_Parent"])))
    output = output[:int(ranking.get("max_ranked_candidates", 10000) or 10000)]
    for rank, row in enumerate(output, start=1): row["Rank"] = rank
    return output, build_modification_evidence_summary(output)


def _final_confidence(score: float, precursor: bool, ions: bool, level: str, known: bool, has_ms2: bool, isobaric: bool, config: dict[str, Any]) -> str:
    if score >= 8 and precursor and ions and level in {"Moderate", "Strong"}: result = "Very High"
    elif score >= 6 and precursor and ions: result = "High"
    elif score >= 4 and precursor and (level in {"Weak", "Moderate", "Strong"} or known): result = "Medium"
    elif score >= 2: result = "Low"
    else: result = "Very Low"
    if _bool(config.get("require_ms2_evidence_for_high_confidence"), True) and not has_ms2 and result in {"Very High", "High"}: result = "Medium"
    if isobaric and precursor and not ions and result in {"Very High", "High"}: result = "Medium"
    return result


def _interpretation(precursor: bool, ions: bool, level: str, ms1: bool, known: bool, ambiguous: bool, low_information: bool) -> str:
    if ambiguous: return "ambiguous-localization"
    if low_information: return "low-information-ion-only"
    if precursor and ions and level in {"Moderate", "Strong"}: return "strong-modified-ms2-evidence"
    if precursor and ions: return "precursor-and-modified-ion-supported"
    if precursor and level == "Weak": return "precursor-supported-localization-weak"
    if ms1 and not precursor and not ions: return "ms1-only-candidate"
    if known and not precursor and not ions: return "known-modification-only-candidate"
    return "insufficient-evidence"


def build_modification_evidence_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    counts: dict[str, int] = {}
    for row in rows: counts[row["Modification_ID"]] = counts.get(row["Modification_ID"], 0) + 1
    top = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:10]
    return [{
        "Total_Ranked_Candidates": len(rows), "Very_High": sum(row["Final_Confidence"] == "Very High" for row in rows),
        "High": sum(row["Final_Confidence"] == "High" for row in rows), "Medium": sum(row["Final_Confidence"] == "Medium" for row in rows),
        "Low": sum(row["Final_Confidence"] == "Low" for row in rows), "Very_Low": sum(row["Final_Confidence"] == "Very Low" for row in rows),
        "Candidates_With_MS2_Precursor_Evidence": sum(row["Has_MS2_Precursor_Evidence"] for row in rows),
        "Candidates_With_Modified_Ion_Evidence": sum(row["Has_Modified_Ion_Evidence"] for row in rows),
        "Candidates_With_Localization_Evidence": sum(row["Has_Localization_Evidence"] for row in rows),
        "Ambiguous_Candidates": sum(row["Ambiguous_Position"] for row in rows),
        "Top_Modification_IDs": "; ".join(f"{mod_id}/{count}" for mod_id, count in top),
        "Notes": "Evidence ranking is a review priority, not a modification call.",
    }]
