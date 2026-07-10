"""Integrate modification evidence into review-oriented rankings (MVP-5.4)."""

from dataclasses import asdict, is_dataclass
from typing import Any

from rna_masshunter.biological_context import score_biological_context


RANKING_COLUMNS = [
    "Rank", "Final_Score", "Final_Confidence", "Final_Interpretation", "Confidence_Limiting_Factor",
    "Modification_ID", "Modification_Name", "Modification_Category", "Target_Base", "Mass_Shift", "Is_Isobaric",
    "Source_Priority", "Curation_Status", "Candidate_Policy_By_Mass_Search",
    "Candidate_Policy_Position_Rule", "Detectability_MS1", "Detectability_MS2",
    "Chemical_Group", "Near_Isobaric_Group",
    "Biological_Context_Score", "Biological_Context_Level", "Biological_Context_Notes",
    "Context_Matched_Priority_Modification", "Context_Matched_Keywords",
    "Context_Focus_Position_Match", "Context_Focus_Position_Distance",
    "Context_Pathway_Supported", "Context_Organism_Supported", "Context_TRNA_Supported",
    "Context_Conflict",
    "Candidate_tRNA_Position", "Candidate_Base", "Parent_Fragment_ID", "Parent_Sequence", "Parent_Start", "Parent_End", "Candidate_Position_In_Parent",
    "Has_MS1_Fragment_Evidence", "MS1_Fragment_Best_Confidence", "MS1_Fragment_Total_Intensity",
    "Has_Known_Modification_Candidate", "Known_Modification_Priority_Score",
    "Has_MS2_Precursor_Evidence", "Num_MS2_Precursor_Candidates", "Best_Precursor_Error_ppm", "Modified_Precursor_Rescue",
    "Has_Modified_Ion_Evidence", "Num_Modified_Ion_Matches", "Num_Informative_Modified_Ion_Matches", "Best_Modified_Ion_Error_ppm",
    "Has_Localization_Evidence", "Localization_Level", "Localization_Score", "Localization_Interpretation", "Num_c_Modified_Ions", "Num_y_Modified_Ions",
    "Organism_Group", "Organism_Species", "Rule_Set", "Organism_Rule_Supported", "TRNA_Context_Supported", "Context_Notes",
    "Ambiguous_Position", "Low_Information_Evidence", "Evidence_Warnings", "Notes",
    "Ambiguity_Group_ID", "Num_Positions_In_Ambiguity_Group", "Position_Discriminating_Evidence",
    "Num_Position_Discriminating_Ions", "Num_Informative_Position_Discriminating_Ions",
    "Position_Ambiguity_Status", "Group_Best_Position", "Group_Best_Position_Score",
    "Candidate_Is_Group_Best", "Ambiguity_Penalty_Applied",
]

SUMMARY_COLUMNS = [
    "Total_Ranked_Candidates", "Very_High", "High", "Medium", "Low", "Very_Low",
    "Candidates_With_MS2_Precursor_Evidence", "Candidates_With_Modified_Ion_Evidence",
    "Candidates_With_Localization_Evidence", "Ambiguous_Candidates", "Top_Modification_IDs",
    "Total_Ambiguity_Groups", "Resolved_Ambiguity_Groups", "Partially_Resolved_Ambiguity_Groups",
    "Ambiguous_Groups", "Candidates_With_Position_Discriminating_Evidence",
    "Candidates_Without_Position_Discriminating_Evidence", "Notes",
    "Candidates_With_Biological_Context_Support", "Candidates_With_Priority_Modification",
    "Candidates_With_Priority_Keyword", "Candidates_With_Focus_Position_Match",
    "Candidates_With_Context_Conflict", "Top_Context_Supported_Modifications", "Top_Context_Keywords",
]

AMBIGUITY_GROUP_COLUMNS = [
    "Ambiguity_Group_ID", "Spectrum_ID", "RT", "Precursor_mz", "Parent_Fragment_ID",
    "Parent_Sequence", "Parent_Start", "Parent_End", "Modification_ID", "Modification_Name",
    "Candidate_Positions_In_Parent", "Candidate_Positions_In_tRNA", "Candidate_Bases",
    "Num_Candidate_Positions", "Num_Modified_Ion_Matches", "Num_Position_Discriminating_Ions",
    "Num_Informative_Position_Discriminating_Ions", "Num_Non_Discriminating_Ions",
    "Best_Position_By_Score", "Best_Position_Score", "Position_Ambiguity_Status",
    "Group_Interpretation", "Notes",
]


def build_ambiguity_groups(
    localization_rows: list[dict[str, Any]], ion_matches: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for row in localization_rows or []:
        key = (str(row.get("Spectrum_ID") or ""), str(row.get("Parent_Fragment_ID") or ""), str(row.get("Modification_ID") or ""))
        grouped.setdefault(key, []).append(row)
    results = []
    for group_index, (key, positions) in enumerate(sorted(grouped.items()), start=1):
        position_values = sorted({int(row["Candidate_Modification_Position_In_Parent"]) for row in positions})
        group_matches = [row for row in ion_matches or [] if (str(row.get("Spectrum_ID") or ""), str(row.get("Parent_Fragment_ID") or ""), str(row.get("Modification_ID") or "")) == key and row.get("Ion_Contains_Modification")]
        discriminating = [row for row in group_matches if row.get("Position_Discriminating_Ion")]
        informative_discriminating = [row for row in discriminating if row.get("Informative_Ion")]
        scores = {int(row["Candidate_Modification_Position_In_Parent"]): _float(row.get("Localization_Score")) for row in positions}
        best_score = max(scores.values(), default=0.0)
        best_positions = [position for position, score in scores.items() if score == best_score]
        best_position = best_positions[0] if len(best_positions) == 1 else ""
        discriminating_positions = {int(row["Candidate_Modification_Position_In_Parent"]) for row in informative_discriminating}
        if not group_matches:
            status, interpretation = "no_localization_evidence", "no-modified-ion-support"
        elif len(position_values) == 1 or (len(discriminating_positions) == 1 and best_position in discriminating_positions):
            status, interpretation = "resolved", "position-resolved-by-discriminating-ions"
        elif discriminating_positions:
            status, interpretation = "partially_resolved", "position-partially-supported"
        else:
            status, interpretation = "ambiguous", "position-ambiguous-same-evidence"
        if status == "ambiguous" and not discriminating:
            interpretation = "precursor-supported-but-position-unresolved"
        first = positions[0]
        results.append({
            "Ambiguity_Group_ID": f"AMBG_{group_index:06d}", "Spectrum_ID": key[0],
            "RT": first.get("RT", ""), "Precursor_mz": first.get("Precursor_mz", ""),
            "Parent_Fragment_ID": key[1], "Parent_Sequence": first.get("Parent_Sequence", ""),
            "Parent_Start": first.get("Parent_Start", ""), "Parent_End": first.get("Parent_End", ""),
            "Modification_ID": key[2], "Modification_Name": first.get("Modification_Name", ""),
            "Candidate_Positions_In_Parent": ";".join(map(str, position_values)),
            "Candidate_Positions_In_tRNA": ";".join(str(row.get("Candidate_Modification_Position_In_tRNA", "")) for row in positions),
            "Candidate_Bases": ";".join(str(row.get("Candidate_Modification_Base", "")) for row in positions),
            "Num_Candidate_Positions": len(position_values), "Num_Modified_Ion_Matches": len(group_matches),
            "Num_Position_Discriminating_Ions": len(discriminating),
            "Num_Informative_Position_Discriminating_Ions": len(informative_discriminating),
            "Num_Non_Discriminating_Ions": len(group_matches) - len(discriminating),
            "Best_Position_By_Score": best_position, "Best_Position_Score": best_score,
            "Position_Ambiguity_Status": status, "Group_Interpretation": interpretation,
            "Notes": "Candidate positions are retained; the group reports whether MS2 ions separate them.",
        })
    return results


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


def _cap_confidence(value: str, cap: str) -> str:
    order = ["Very Low", "Low", "Medium", "High", "Very High"]
    if value not in order or cap not in order:
        return value
    return order[min(order.index(value), order.index(cap))]


def _empty_context_result() -> dict[str, Any]:
    return {
        "Biological_Context_Score": 0.0, "Biological_Context_Level": "None",
        "Biological_Context_Notes": "", "Context_Matched_Priority_Modification": False,
        "Context_Matched_Keywords": "", "Context_Focus_Position_Match": "",
        "Context_Focus_Position_Distance": "", "Context_Pathway_Supported": False,
        "Context_Organism_Supported": False, "Context_TRNA_Supported": False,
        "Context_Conflict": False,
    }


def _rule_modification_ids(rule_set: dict[str, Any] | None) -> set[str]:
    ids: set[str] = set()
    def visit(value: Any) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                if str(key).lower() in {"modification_id", "modification", "mod_id"} and isinstance(child, (str, int)):
                    ids.add(str(child))
                if str(key).lower() in {"preferred_modifications", "supported_modifications", "modification_ids"} and isinstance(child, list):
                    ids.update(str(item) for item in child if isinstance(item, (str, int)))
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
    pathways: list[dict[str, Any]] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    ranking = getattr(config, "modification_evidence_ranking", {}) or {}
    if not _bool(ranking.get("enabled"), True):
        return [], []
    weights = ranking.get("weights", {}) or {}
    weight = lambda name, default: _float(weights.get(name), default)
    modification_lookup = {str(getattr(item, "id", "")): item for item in modifications or []}
    near_isobaric_counts: dict[str, int] = {}
    for modification_item in modifications or []:
        group = str((getattr(modification_item, "raw", {}) or {}).get("near_isobaric_group") or "")
        if group:
            near_isobaric_counts[group] = near_isobaric_counts.get(group, 0) + 1
    fragment_lookup = {str(getattr(item, "fragment_id", "")): item for item in theoretical_fragments or []}
    rule_ids = _rule_modification_ids(rule_set)
    localization = list(ms2_results.get("MS2_Modification_Localization_Evidence", []) or [])
    precursors = list(ms2_results.get("MS2_Modified_Precursor_Candidates", []) or [])
    ion_matches = list(ms2_results.get("MS2_Modified_Ion_Matches", []) or [])
    ambiguity_groups = list(ms2_results.get("Modification_Ambiguity_Groups", []) or [])

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
        source = mod_raw.get("source", {}) or {}
        source_priority_raw = mod_raw.get("source_priority") or source.get("source_priority") or ""
        source_priority = source_priority_raw.get("mass_shift", "") if isinstance(source_priority_raw, dict) else str(source_priority_raw)
        curation_status = str(mod_raw.get("curation_status") or source.get("curation_status") or (mod_raw.get("curation", {}) or {}).get("status") or "")
        candidate_policy = mod_raw.get("candidate_policy", {}) or {}
        detectability = mod_raw.get("detectability", {}) or {}
        chemical_group = str(mod_raw.get("chemical_group") or "")
        near_isobaric_group = str(mod_raw.get("near_isobaric_group") or "")
        loc = item.get("localization", {})
        precursor_rows = item.get("precursors", [])
        known_rows = item.get("known", [])
        fragment = fragment_lookup.get(fragment_id)
        matching_ms1 = [_raw(match) for match in fragment_ms1_matches or [] if str(_raw(match).get("fragment_id") or "") == fragment_id]
        matching_ions = [row for row in ion_matches if str(row.get("Modification_ID") or "") == mod_id and str(row.get("Parent_Fragment_ID") or "") == fragment_id and (position == "" or row.get("Candidate_Modification_Position_In_Parent") == position)]
        modified_ions = [row for row in matching_ions if row.get("Ion_Contains_Modification")]
        informative_ions = [row for row in modified_ions if row.get("Informative_Ion")]
        level = str(loc.get("Localization_Level") or "None")
        group = next((row for row in ambiguity_groups if str(row.get("Spectrum_ID") or "") == str(loc.get("Spectrum_ID") or "") and str(row.get("Parent_Fragment_ID") or "") == fragment_id and str(row.get("Modification_ID") or "") == mod_id), {})
        ambiguity_status = str(group.get("Position_Ambiguity_Status") or "")
        ambiguous = ambiguity_status == "ambiguous" or loc.get("Localization_Interpretation") in {"ambiguous-multiple-positions", "position-ambiguous-non-discriminating-ions"}
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
        if curation_status == "manually_checked": score += weight("curation_manually_checked", 0.5)
        if source_priority in {"user_pdf_for_mass_shift", "user_pdf"}: score += weight("source_user_pdf", 0.5)
        if modified_ions and _bool(detectability.get("ms2"), False): score += weight("detectability_ms2_supported", 0.5)
        if low_information: score += weight("low_information_penalty", -1.0)
        if ambiguous and not group: score += weight("ambiguous_position_penalty", -1.0)
        ambiguity_penalty = weight("ambiguity_penalty", -1.5) if ambiguous and group else 0.0
        score += ambiguity_penalty
        if isobaric and precursor_rows and not modified_ions: score += weight("isobaric_precursor_penalty", -2.0)
        has_ms2 = bool(precursor_rows or modified_ions)
        best_ms1 = max(matching_ms1, key=lambda row: _confidence_rank(row.get("confidence")), default={})
        best_precursor_error = min((abs(_float(row.get("Precursor_Error_ppm"))) for row in precursor_rows), default="")
        best_ion_error = min((abs(_float(row.get("Mass_Error_ppm"))) for row in modified_ions), default="")
        num_c = int(loc.get("Num_c_Modified_Ions", 0) or 0)
        num_y = int(loc.get("Num_y_Modified_Ions", 0) or 0)
        num_discriminating = int(loc.get("Num_Position_Discriminating_Modified_Ions", 0) or 0)
        num_informative_discriminating = int(loc.get("Num_Informative_Position_Discriminating_Modified_Ions", 0) or 0)
        position_discriminating = bool(loc.get("Has_Position_Discriminating_Evidence"))
        if modified_ions and not position_discriminating:
            score += weight("non_discriminating_ion_penalty", -0.5)
        context_candidate = {
            "Modification_ID": mod_id,
            "Modification_Name": mod_raw.get("name") or getattr(modification, "symbol", None) or mod_id,
            "Chemical_Group": chemical_group, "Near_Isobaric_Group": near_isobaric_group,
            "Source_Priority": source_priority, "Notes": mod_raw.get("notes") or "",
            "Candidate_tRNA_Position": trna_position,
        }
        context_result = score_biological_context(
            context_candidate, modification, config, rule_set=rule_set, pathways=pathways,
            rule_supported=rule_supported,
        ) if _bool(ranking.get("use_biological_context"), True) else _empty_context_result()
        score += _float(context_result.get("Biological_Context_Score"), 0.0)
        confidence = _final_confidence(
            score, bool(precursor_rows), bool(modified_ions), level, bool(known_rows), has_ms2,
            isobaric, len(informative_ions), num_c, num_y, best_precursor_error, best_ion_error,
            position_discriminating, num_informative_discriminating, ambiguity_status, ranking,
        )
        if not _bool(candidate_policy.get("include_by_mass_search"), True) and not rule_supported and not modified_ions and confidence in {"Very High", "High"}:
            confidence = "Medium"
        if context_result.get("Biological_Context_Score", 0) and not has_ms2 and _bool(ranking.get("require_ms_evidence_for_context_boosted_high"), True):
            confidence = _cap_confidence(confidence, str(ranking.get("cap_context_only_confidence") or "Medium"))
        limiting_factor = _confidence_limiting_factor(
            bool(precursor_rows), bool(modified_ions), level, len(informative_ions), num_c, num_y,
            rule_supported, bool(matching_ms1), position_discriminating, ambiguity_status,
        )
        interpretation = _interpretation(
            bool(precursor_rows), bool(modified_ions), level, bool(matching_ms1), bool(known_rows),
            ambiguous, low_information, len(informative_ions), num_c, num_y,
        )
        target_bases = getattr(modification, "target_bases", []) if modification else []
        output.append({
            "Final_Score": score, "Final_Confidence": confidence, "Final_Interpretation": interpretation,
            "Confidence_Limiting_Factor": limiting_factor,
            "Modification_ID": mod_id, "Modification_Name": mod_raw.get("name") or getattr(modification, "symbol", None) or (precursor_rows[0].get("Modification_Name") if precursor_rows else mod_id),
            "Modification_Category": getattr(modification, "category", ""), "Target_Base": ",".join(target_bases), "Mass_Shift": shift, "Is_Isobaric": isobaric,
            "Source_Priority": source_priority, "Curation_Status": curation_status,
            "Candidate_Policy_By_Mass_Search": _bool(candidate_policy.get("include_by_mass_search"), True),
            "Candidate_Policy_Position_Rule": _bool(candidate_policy.get("include_if_position_rule_exists"), False),
            "Detectability_MS1": detectability.get("ms1", ""), "Detectability_MS2": detectability.get("ms2", ""),
            "Chemical_Group": chemical_group, "Near_Isobaric_Group": near_isobaric_group,
            **context_result,
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
            "Localization_Interpretation": loc.get("Localization_Interpretation", ""), "Num_c_Modified_Ions": num_c, "Num_y_Modified_Ions": num_y,
            "Organism_Group": getattr(config, "organism", {}).get("group", ""), "Organism_Species": getattr(config, "organism", {}).get("species", ""),
            "Rule_Set": getattr(config, "organism", {}).get("rule_set", ""), "Organism_Rule_Supported": rule_supported,
            "TRNA_Context_Supported": trna_supported, "Context_Notes": "Rule/context points are applied only when explicit loaded data supports them.",
            "Ambiguous_Position": ambiguous, "Low_Information_Evidence": low_information,
            "Evidence_Warnings": "; ".join(filter(None, ["isobaric precursor evidence is non-specific" if isobaric else "", "near-isobaric alternative" if near_isobaric_group and near_isobaric_counts.get(near_isobaric_group, 0) > 1 else "", "excluded from blind mass search by curated policy" if not _bool(candidate_policy.get("include_by_mass_search"), True) else "", "ambiguous localization" if ambiguous else "", "1 nt/low-information ion evidence" if low_information else ""])),
            "Notes": "Ranking prioritizes candidates; it does not confirm modification identity or position.",
            "Ambiguity_Group_ID": group.get("Ambiguity_Group_ID", ""),
            "Num_Positions_In_Ambiguity_Group": group.get("Num_Candidate_Positions", 1 if position != "" else 0),
            "Position_Discriminating_Evidence": position_discriminating,
            "Num_Position_Discriminating_Ions": num_discriminating,
            "Num_Informative_Position_Discriminating_Ions": num_informative_discriminating,
            "Position_Ambiguity_Status": ambiguity_status,
            "Group_Best_Position": group.get("Best_Position_By_Score", ""),
            "Group_Best_Position_Score": group.get("Best_Position_Score", ""),
            "Candidate_Is_Group_Best": bool(group and group.get("Best_Position_By_Score") == position),
            "Ambiguity_Penalty_Applied": ambiguity_penalty,
        })
    minimum = _float(ranking.get("min_final_score_to_report"), 0.0)
    output = [row for row in output if row["Final_Score"] >= minimum]
    output.sort(key=lambda row: (-row["Final_Score"], row["Modification_ID"], row["Parent_Fragment_ID"], str(row["Candidate_Position_In_Parent"])))
    output = output[:int(ranking.get("max_ranked_candidates", 10000) or 10000)]
    for rank, row in enumerate(output, start=1): row["Rank"] = rank
    return output, build_modification_evidence_summary(output, ambiguity_groups)


def _final_confidence(
    score: float, precursor: bool, ions: bool, level: str, known: bool, has_ms2: bool,
    isobaric: bool, informative_count: int, num_c: int, num_y: int,
    precursor_error: Any, ion_error: Any, position_discriminating: bool,
    informative_discriminating_count: int, ambiguity_status: str, config: dict[str, Any],
) -> str:
    both_series = num_c >= 1 and num_y >= 1
    require_discrimination = _bool(config.get("require_position_discriminating_ions_for_localization_confidence"), True)
    discrimination_gate = position_discriminating or not require_discrimination
    localization_gate = level in {"Moderate", "Strong"} and discrimination_gate
    min_discriminating_for_high = int(config.get("min_informative_discriminating_ions_for_high", 2) or 2)
    multi_ion_gate = informative_count >= 2 and both_series and informative_discriminating_count >= min_discriminating_for_high
    good_errors = precursor_error != "" and ion_error != "" and _float(precursor_error, 999) <= 5.0 and _float(ion_error, 999) <= 5.0
    if score >= 8 and precursor and ions and level == "Strong" and informative_count >= 3 and both_series and good_errors and informative_discriminating_count >= 2 and ambiguity_status in {"", "resolved"}: result = "Very High"
    elif score >= 6 and precursor and ions and (localization_gate or multi_ion_gate): result = "High"
    elif precursor and ions and level == "Weak" and informative_count >= 1: result = "Medium"
    elif score >= 4 and precursor and (level in {"Weak", "Moderate", "Strong"} or known): result = "Medium"
    elif score >= 2: result = "Low"
    else: result = "Very Low"
    if _bool(config.get("require_ms2_evidence_for_high_confidence"), True) and not has_ms2 and result in {"Very High", "High"}: result = "Medium"
    if isobaric and precursor and not ions and result in {"Very High", "High"}: result = "Medium"
    if require_discrimination and not position_discriminating and result in {"Very High", "High"}: result = "Medium"
    if ambiguity_status == "ambiguous" and result in {"Very High", "High"}: result = "Medium"
    return result


def _interpretation(
    precursor: bool, ions: bool, level: str, ms1: bool, known: bool, ambiguous: bool,
    low_information: bool, informative_count: int, num_c: int, num_y: int,
) -> str:
    if ambiguous: return "ambiguous-localization"
    if low_information: return "low-information-ion-only"
    if precursor and ions and level in {"Moderate", "Strong"}: return "strong-modified-ms2-evidence"
    if precursor and ions and informative_count == 1: return "precursor-and-single-modified-ion-supported"
    if precursor and ions and (num_c == 0 or num_y == 0): return "candidate-needs-additional-ms2-support"
    if precursor and ions: return "precursor-and-modified-ion-supported"
    if precursor and level == "Weak": return "precursor-supported-localization-weak"
    if ms1 and not precursor and not ions: return "ms1-only-candidate"
    if known and not precursor and not ions: return "known-modification-only-candidate"
    return "insufficient-evidence"


def _confidence_limiting_factor(
    precursor: bool, ions: bool, level: str, informative_count: int, num_c: int, num_y: int,
    rule_supported: bool, ms1_supported: bool, position_discriminating: bool,
    ambiguity_status: str,
) -> str:
    factors = []
    if precursor and not ions: factors.append("precursor-only")
    if level == "Weak": factors.append("weak-localization")
    if ions and informative_count <= 1: factors.append("single-modified-ion")
    if ions and (num_c == 0 or num_y == 0): factors.append("one-sided-ion-series")
    if not rule_supported: factors.append("no-known-modification-rule")
    if not ms1_supported: factors.append("no-ms1-fragment-evidence")
    if ambiguity_status == "ambiguous": factors.append("position-ambiguous")
    if ions and not position_discriminating: factors.append("no-position-discriminating-ion")
    return "; ".join(factors)


def build_modification_evidence_summary(
    rows: list[dict[str, Any]], ambiguity_groups: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    ambiguity_groups = ambiguity_groups or []
    counts: dict[str, int] = {}
    for row in rows: counts[row["Modification_ID"]] = counts.get(row["Modification_ID"], 0) + 1
    top = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:10]
    context_modifications: dict[str, int] = {}
    context_keywords: dict[str, int] = {}
    for row in rows:
        if _float(row.get("Biological_Context_Score"), 0.0) > 0:
            mod_id = str(row.get("Modification_ID") or "")
            context_modifications[mod_id] = context_modifications.get(mod_id, 0) + 1
        for keyword in str(row.get("Context_Matched_Keywords") or "").split(";"):
            if keyword:
                context_keywords[keyword] = context_keywords.get(keyword, 0) + 1
    return [{
        "Total_Ranked_Candidates": len(rows), "Very_High": sum(row["Final_Confidence"] == "Very High" for row in rows),
        "High": sum(row["Final_Confidence"] == "High" for row in rows), "Medium": sum(row["Final_Confidence"] == "Medium" for row in rows),
        "Low": sum(row["Final_Confidence"] == "Low" for row in rows), "Very_Low": sum(row["Final_Confidence"] == "Very Low" for row in rows),
        "Candidates_With_MS2_Precursor_Evidence": sum(row["Has_MS2_Precursor_Evidence"] for row in rows),
        "Candidates_With_Modified_Ion_Evidence": sum(row["Has_Modified_Ion_Evidence"] for row in rows),
        "Candidates_With_Localization_Evidence": sum(row["Has_Localization_Evidence"] for row in rows),
        "Ambiguous_Candidates": sum(row["Ambiguous_Position"] for row in rows),
        "Top_Modification_IDs": "; ".join(f"{mod_id}/{count}" for mod_id, count in top),
        "Total_Ambiguity_Groups": len(ambiguity_groups),
        "Resolved_Ambiguity_Groups": sum(row.get("Position_Ambiguity_Status") == "resolved" for row in ambiguity_groups),
        "Partially_Resolved_Ambiguity_Groups": sum(row.get("Position_Ambiguity_Status") == "partially_resolved" for row in ambiguity_groups),
        "Ambiguous_Groups": sum(row.get("Position_Ambiguity_Status") == "ambiguous" for row in ambiguity_groups),
        "Candidates_With_Position_Discriminating_Evidence": sum(bool(row.get("Position_Discriminating_Evidence")) for row in rows),
        "Candidates_Without_Position_Discriminating_Evidence": sum(not bool(row.get("Position_Discriminating_Evidence")) for row in rows),
        "Candidates_With_Biological_Context_Support": sum(_float(row.get("Biological_Context_Score"), 0.0) > 0 for row in rows),
        "Candidates_With_Priority_Modification": sum(bool(row.get("Context_Matched_Priority_Modification")) for row in rows),
        "Candidates_With_Priority_Keyword": sum(bool(row.get("Context_Matched_Keywords")) for row in rows),
        "Candidates_With_Focus_Position_Match": sum(bool(row.get("Context_Focus_Position_Match")) for row in rows),
        "Candidates_With_Context_Conflict": sum(bool(row.get("Context_Conflict")) for row in rows),
        "Top_Context_Supported_Modifications": "; ".join(f"{key}/{value}" for key, value in sorted(context_modifications.items(), key=lambda item: (-item[1], item[0]))[:10]),
        "Top_Context_Keywords": "; ".join(f"{key}/{value}" for key, value in sorted(context_keywords.items(), key=lambda item: (-item[1], item[0]))[:10]),
        "Notes": "High confidence requires localization-level support or multi-ion c/y support. Weak localization candidates are treated as Medium review candidates. Candidates sharing the same parent fragment and modification may form ambiguity groups. Position confidence requires position-discriminating ions. Biological context boosts prioritization but does not establish modification identity without mass/MS evidence. Evidence ranking is not a modification call.",
    }]
