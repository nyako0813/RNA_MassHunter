"""Diagnostic-only biological position prior for MS/MS modification candidates.

This module never changes candidate membership, evidence scores, confidence, or rank.
Positions are input-sequence 1-based; no Sprinzl numbering is assumed.
"""
from __future__ import annotations

from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

import yaml

POSITION_PRIOR_COLUMNS = [
    "Modification_ID", "Modification_Family", "Candidate_tRNA_Position", "Candidate_Base",
    "Position_Class", "Position_Prior_Score", "Position_Prior_Level", "Position_Prior_Reason",
    "Landmark_Name", "Landmark_Position", "Position_Numbering_System", "Rule_Source", "Notes",
]
BIOLOGICAL_PLAUSIBILITY_COLUMNS = [
    "Rank", "Modification_ID", "Modification_Name", "Parent_Fragment_ID", "Candidate_tRNA_Position",
    "Candidate_Base", "Modification_Family", "Position_Class", "Position_Prior_Score",
    "Parent_Base_Compatibility", "Parent_Base_Prior_Score", "Parent_Base_Reason",
    "MS2_Localization_Evidence", "Structure_Ambiguity_Status", "Alternative_Structural_Candidates",
    "Structure_Discriminating_Evidence", "Biological_Plausibility_Score", "Biological_Plausibility_Level",
    "Shadow_Final_Score", "Shadow_Final_Confidence", "Shadow_Only", "Notes",
]
SHADOW_RANKING_COLUMNS = [
    "Modification_Family", "Position_Class", "Position_Prior_Score", "Position_Prior_Level",
    "Position_Prior_Reason", "Landmark_Name", "Landmark_Position", "Position_Numbering_System",
    "Parent_Base_Compatibility", "Parent_Base_Prior_Score", "Parent_Base_Reason",
    "MS2_Localization_Evidence", "Structure_Ambiguity_Status", "Alternative_Structural_Candidates",
    "Structure_Discriminating_Evidence", "Biological_Plausibility_Score", "Biological_Plausibility_Level",
    "Shadow_Final_Score", "Shadow_Final_Confidence", "Shadow_Only",
]
DIAGNOSTIC_COLUMNS = [
    "Total_Candidates", "Evaluated_Candidates", "Unknown_Position_Candidates", "Compatible_Parent_Base",
    "Incompatible_Parent_Base", "Ambiguous_Parent_Base", "Unknown_Parent_Base",
    "Canonical_Position", "Adjacent_To_Canonical", "Reported_Noncanonical_Position",
    "Chemically_Possible_But_Unreported", "Biologically_Inconsistent", "Unknown_Position_Class",
    "Structure_Unresolved", "Apply_To_Final_Score", "Rows_Changed", "Notes",
]

DEFAULT_LEVELS = {"high": 2.0, "moderate": 1.0, "neutral": 0.0, "warning": -1.0, "inconsistent": -2.0, "unknown": 0.0}
DEFAULT_PARENT_SCORES = {"compatible": 1.0, "incompatible": -2.0, "ambiguous": 0.0, "unknown": 0.0}


def _raw(value: Any) -> dict[str, Any]:
    if is_dataclass(value):
        return asdict(value)
    return dict(value) if isinstance(value, dict) else {}


def _bool(value: Any, default: bool = False) -> bool:
    if value is None:
        return default
    return value if isinstance(value, bool) else str(value).strip().lower() in {"1", "true", "yes", "on"}


def _float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _position(value: Any) -> int | None:
    try:
        number = int(float(value))
        return number if number > 0 else None
    except (TypeError, ValueError):
        return None


def load_position_prior_rules(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    if not source.exists():
        return {"version": "missing", "families": []}
    with source.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return data if isinstance(data, dict) else {"version": "invalid", "families": []}


def _family(mod_id: str, rules: dict[str, Any]) -> dict[str, Any]:
    key = str(mod_id or "").casefold()
    for rule in rules.get("families", []) or []:
        aliases = [str(item).casefold() for item in rule.get("modification_ids", [])]
        if key in aliases:
            return rule
    return {}


def _landmarks(config: Any) -> dict[str, int]:
    prior = ((getattr(config, "ms2_annotation", {}) or {}).get("biological_position_prior") or {})
    numbering = prior.get("position_numbering") or {}
    raw = numbering.get("landmarks") or prior.get("landmarks") or {}
    result = {str(k): v for k, v in raw.items() if _position(v) is not None}
    # wobble_position is an explicitly configured input-sequence landmark, not a hard-coded convention.
    configured_wobble = (getattr(config, "sequence", {}) or {}).get("wobble_position")
    if "wobble" not in result and _position(configured_wobble) is not None:
        result["wobble"] = configured_wobble
    return {key: int(value) for key, value in result.items()}


def _position_prior(position: int | None, rule: dict[str, Any], landmarks: dict[str, int], scores: dict[str, float]) -> tuple[str, float, str, str, Any]:
    if position is None or not rule:
        return "unknown", scores["unknown"], "unknown", "Position or family rule is unavailable.", "", ""
    landmark_name = str(rule.get("canonical_landmark") or "")
    landmark = landmarks.get(landmark_name)
    reported = {_position(item) for item in rule.get("reported_noncanonical_positions", [])}
    reported.discard(None)
    if landmark is None:
        if position in reported:
            return "reported_noncanonical_position", scores["moderate"], "moderate", "Position is explicitly reported by the configured rule.", landmark_name, ""
        return "unknown", scores["unknown"], "unknown", "Required landmark is not configured; no canonical position was assumed.", landmark_name, ""
    if position == landmark:
        return "canonical_position", scores["high"], "high", "Candidate matches the configured canonical landmark.", landmark_name, landmark
    if abs(position - landmark) == 1:
        return "adjacent_to_canonical", scores["moderate"], "moderate", "Candidate is adjacent to the configured canonical landmark.", landmark_name, landmark
    if position in reported:
        return "reported_noncanonical_position", scores["moderate"], "moderate", "Position is explicitly reported as noncanonical.", landmark_name, landmark
    max_distance = int(rule.get("chemically_possible_distance", 3) or 3)
    if abs(position - landmark) <= max_distance:
        return "chemically_possible_but_unreported", scores["neutral"], "neutral", "Position is near the landmark but is not explicitly reported.", landmark_name, landmark
    return "biologically_inconsistent", scores["inconsistent"], "inconsistent", "Position is distant from the configured family landmark.", landmark_name, landmark


def _parent_base(row: dict[str, Any], modification: dict[str, Any], sequence: str) -> tuple[str, float, str, str]:
    position = _position(row.get("Candidate_tRNA_Position"))
    observed = str(row.get("Candidate_Base") or "").upper()
    if not observed and position and position <= len(sequence):
        observed = sequence[position - 1].upper()
    allowed = [str(item).upper() for item in modification.get("target_bases", []) if item]
    if not observed or not allowed:
        status = "unknown"
        reason = "Candidate base or modification target base is unavailable."
    elif observed in allowed:
        status = "compatible"
        reason = f"Observed parent base {observed} is allowed ({','.join(allowed)})."
    elif observed in {"N", "X", "?"} or len(observed) != 1:
        status = "ambiguous"
        reason = "Observed parent base is ambiguous."
    else:
        status = "incompatible"
        reason = f"Observed parent base {observed} is not allowed ({','.join(allowed)})."
    return status, observed, reason, ";".join(allowed)


def _structure_groups(rows: list[dict[str, Any]], modifications: dict[str, dict[str, Any]], tolerance: float) -> dict[int, tuple[str, str]]:
    result: dict[int, tuple[str, str]] = {}
    for index, row in enumerate(rows):
        alternatives = []
        current_id = str(row.get("Modification_ID") or "")
        current = modifications.get(current_id, {})
        mass = _float(current.get("mass_shift_from_unmodified"), float("nan"))
        group = str(current.get("near_isobaric_group") or "")
        candidate_base = str(row.get("Candidate_Base") or "").upper()
        for other_id, other in modifications.items():
            if other_id == current_id:
                continue
            other_targets = {str(base).upper() for base in other.get("target_bases", []) if base}
            if candidate_base and other_targets and candidate_base not in other_targets:
                continue
            other_mass = _float(other.get("mass_shift_from_unmodified"), float("nan"))
            same_group = bool(group and group == str(other.get("near_isobaric_group") or ""))
            same_mass = mass == mass and other_mass == other_mass and abs(mass - other_mass) <= tolerance
            if same_group or same_mass:
                alternatives.append(other_id)
        alternatives = sorted(set(alternatives))
        if alternatives:
            status = "position_resolved_structure_unresolved" if row.get("Position_Discriminating_Evidence") else "position_and_structure_unresolved"
        else:
            status = "no_structural_alternative_identified"
        result[index] = (status, ";".join(alternatives))
    return result


def evaluate_biological_position_priors(config: Any, ranking_rows: list[dict[str, Any]], modifications: list[Any], rules: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    prior_cfg = ((getattr(config, "ms2_annotation", {}) or {}).get("biological_position_prior") or {})
    enabled = _bool(prior_cfg.get("enabled"), True)
    apply_final = _bool(prior_cfg.get("apply_to_final_score"), False)
    if not enabled:
        return [dict(row) for row in ranking_rows], [], [], [{"Total_Candidates": len(ranking_rows), "Evaluated_Candidates": 0, "Apply_To_Final_Score": apply_final, "Rows_Changed": 0, "Notes": "disabled"}]
    score_cfg = {**DEFAULT_LEVELS, **(prior_cfg.get("position_prior_scores") or {})}
    parent_scores = {**DEFAULT_PARENT_SCORES, **(prior_cfg.get("parent_base_scores") or {})}
    mod_map = {str(_raw(item).get("id") or ""): _raw(item) for item in modifications}
    landmarks = _landmarks(config)
    sequence = str((getattr(config, "sequence", {}) or {}).get("sequence") or "").upper()
    structure = _structure_groups(ranking_rows, mod_map, _float(prior_cfg.get("structural_mass_tolerance_da"), 0.001))
    position_rows, plausibility_rows, enriched = [], [], []
    for index, original in enumerate(ranking_rows):
        row = dict(original)
        mod_id = str(row.get("Modification_ID") or "")
        mod = mod_map.get(mod_id, {})
        family_rule = _family(mod_id, rules)
        family = str(family_rule.get("id") or "unclassified")
        pos = _position(row.get("Candidate_tRNA_Position"))
        if pos is None:
            start, local = _position(row.get("Parent_Start")), _position(row.get("Candidate_Position_In_Parent"))
            pos = start + local - 1 if start and local else None
        row["Candidate_tRNA_Position"] = pos if pos is not None else row.get("Candidate_tRNA_Position", "")
        position_class, position_score, level, reason, landmark_name, landmark_position = _position_prior(pos, family_rule, landmarks, score_cfg)
        organism_group = str((getattr(config, "organism", {}) or {}).get("group") or "").casefold()
        organism_warning = str((family_rule.get("organism_warnings") or {}).get(organism_group) or "")
        if organism_warning:
            position_score += _float(score_cfg.get("warning"), -1.0)
            level = "warning"
            reason = f"{reason} Organism warning: {organism_warning}"
        compatibility, observed_base, base_reason, _ = _parent_base(row, mod, sequence)
        base_score = _float(parent_scores.get(compatibility), 0.0)
        structure_status, alternatives = structure[index]
        ms2_localization = str(row.get("Localization_Level") or "None")
        plausibility_score = position_score + base_score
        plausibility_level = "high" if plausibility_score >= 2 else "moderate" if plausibility_score >= 1 else "warning" if plausibility_score < 0 else "neutral"
        shadow_score = _float(row.get("Final_Score")) + plausibility_score
        shadow_confidence = str(row.get("Final_Confidence") or "")
        shadow = {
            "Modification_Family": family, "Position_Class": position_class,
            "Position_Prior_Score": position_score, "Position_Prior_Level": level,
            "Position_Prior_Reason": reason, "Landmark_Name": landmark_name,
            "Landmark_Position": landmark_position, "Position_Numbering_System": "input_sequence_1_based",
            "Parent_Base_Compatibility": compatibility, "Parent_Base_Prior_Score": base_score,
            "Parent_Base_Reason": base_reason, "MS2_Localization_Evidence": ms2_localization,
            "Structure_Ambiguity_Status": structure_status,
            "Alternative_Structural_Candidates": alternatives,
            "Structure_Discriminating_Evidence": False,
            "Biological_Plausibility_Score": plausibility_score,
            "Biological_Plausibility_Level": plausibility_level,
            "Shadow_Final_Score": shadow_score, "Shadow_Final_Confidence": shadow_confidence,
            "Shadow_Only": not apply_final,
        }
        row.update(shadow)
        enriched.append(row)
        position_rows.append({
            "Modification_ID": mod_id, "Modification_Family": family,
            "Candidate_tRNA_Position": pos if pos is not None else "", "Candidate_Base": observed_base,
            "Position_Class": position_class, "Position_Prior_Score": position_score,
            "Position_Prior_Level": level, "Position_Prior_Reason": reason,
            "Landmark_Name": landmark_name, "Landmark_Position": landmark_position,
            "Position_Numbering_System": "input_sequence_1_based", "Rule_Source": rules.get("version", ""),
            "Notes": "Diagnostic prior only; no Sprinzl numbering is assumed.",
        })
        plausibility_rows.append({
            "Rank": row.get("Rank"), "Modification_ID": mod_id, "Modification_Name": row.get("Modification_Name"),
            "Parent_Fragment_ID": row.get("Parent_Fragment_ID"), "Candidate_tRNA_Position": pos if pos is not None else "",
            "Candidate_Base": observed_base, **{key: shadow[key] for key in BIOLOGICAL_PLAUSIBILITY_COLUMNS if key in shadow},
            "Notes": "Final_Score, Final_Confidence, Rank, and candidate membership are unchanged.",
        })
    counts = {name: 0 for name in ["compatible", "incompatible", "ambiguous", "unknown"]}
    classes = {name: 0 for name in ["canonical_position", "adjacent_to_canonical", "reported_noncanonical_position", "chemically_possible_but_unreported", "biologically_inconsistent", "unknown"]}
    for row in plausibility_rows:
        counts[row["Parent_Base_Compatibility"]] += 1
        classes[row["Position_Class"]] += 1
    diagnostics = [{
        "Total_Candidates": len(ranking_rows), "Evaluated_Candidates": len(plausibility_rows),
        "Unknown_Position_Candidates": classes["unknown"], "Compatible_Parent_Base": counts["compatible"],
        "Incompatible_Parent_Base": counts["incompatible"], "Ambiguous_Parent_Base": counts["ambiguous"],
        "Unknown_Parent_Base": counts["unknown"], "Canonical_Position": classes["canonical_position"],
        "Adjacent_To_Canonical": classes["adjacent_to_canonical"],
        "Reported_Noncanonical_Position": classes["reported_noncanonical_position"],
        "Chemically_Possible_But_Unreported": classes["chemically_possible_but_unreported"],
        "Biologically_Inconsistent": classes["biologically_inconsistent"], "Unknown_Position_Class": classes["unknown"],
        "Structure_Unresolved": sum("unresolved" in row["Structure_Ambiguity_Status"] for row in plausibility_rows),
        "Apply_To_Final_Score": apply_final, "Rows_Changed": 0,
        "Notes": "Shadow evaluation only; existing final score/confidence/rank are not mutated.",
    }]
    return enriched, position_rows, plausibility_rows, diagnostics
