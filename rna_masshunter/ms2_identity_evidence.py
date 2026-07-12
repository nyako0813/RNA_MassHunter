"""Candidate-level shadow summary of existing MS/MS identity evidence.

The module aggregates existing report rows only. It does not rematch spectra or
change ranking, confidence, candidate membership, localization, or biology.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Any

IDENTITY_COLUMNS = [
    "Rank", "Modification_ID", "Modification_Name", "Parent_Fragment_ID",
    "Candidate_tRNA_Position", "Candidate_Base",
    "Has_Modified_Fragment_Ion_Evidence", "Modified_Fragment_Match_Count",
    "Unique_Modified_Fragment_Ion_Count", "Modified_Fragment_Ion_Series",
    "Supporting_Modified_Fragment_Match_IDs", "Best_Modified_Fragment_Error_ppm",
    "Maximum_Modified_Fragment_Intensity", "Position_Localization_Status",
    "Position_Discriminating_Ion_Count", "Structural_Isomer_Group_ID",
    "Structure_Resolution_Status", "Alternative_Modification_IDs",
    "MS2_Identity_Evidence_Level", "Shadow_MS2_Identity_Score",
    "Shadow_MS2_Identity_Confidence", "Shadow_MS2_Identity_Priority",
    "MS2_Identity_Evidence_Reason", "MS2_Identity_Warnings",
]
IDENTITY_SHADOW_COLUMNS = IDENTITY_COLUMNS[6:]


def _bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().lower() in {"1", "true", "yes", "on"}


def _float(value: Any, default: float | None = None) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _position(value: Any) -> int | None:
    number = _float(value)
    return int(number) if number is not None and number > 0 else None


def _key(row: dict[str, Any], position_name: str) -> tuple[str, str, int | None]:
    return (
        str(row.get("Modification_ID") or ""),
        str(row.get("Parent_Fragment_ID") or ""),
        _position(row.get(position_name)),
    )


def _match_id(row: dict[str, Any]) -> str:
    values = (
        row.get("Spectrum_ID"), row.get("Scan_Index"), row.get("Observed_mz"),
        row.get("Ion_ID"), row.get("Theoretical_mz"),
    )
    return ":".join(str(value if value not in (None, "") else "NA") for value in values)


def _deduplicate_matches(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    unique: dict[str, dict[str, Any]] = {}
    for row in rows:
        if not _bool(row.get("Ion_Contains_Modification")):
            continue
        unique.setdefault(_match_id(row), row)
    return list(unique.values())


def _multiple_assignment_warning(rows: list[dict[str, Any]]) -> bool:
    observed: dict[tuple[str, str, str], set[tuple[str, str]]] = defaultdict(set)
    for row in rows:
        peak = (str(row.get("Spectrum_ID") or ""), str(row.get("Scan_Index") or ""), str(row.get("Observed_mz") or ""))
        theoretical = (str(row.get("Ion_ID") or ""), str(row.get("Theoretical_mz") or ""))
        observed[peak].add(theoretical)
    return any(len(assignments) > 1 for assignments in observed.values())


def _ambiguity_positions(group: dict[str, Any] | None) -> list[int]:
    if not group:
        return []
    values = str(group.get("Candidate_Positions_In_tRNA") or "").split(";")
    return sorted({position for value in values if (position := _position(value)) is not None})


def _position_status(row: dict[str, Any], group: dict[str, Any] | None, has_modified: bool) -> str:
    if _bool(row.get("Position_Discriminating_Evidence")):
        return "position_resolved"
    positions = _ambiguity_positions(group)
    if len(positions) > 1:
        if positions[-1] - positions[0] == len(positions) - 1:
            return "adjacent_positions_ambiguous"
        return "multiple_positions_ambiguous"
    if has_modified or _bool(row.get("Has_Localization_Evidence")):
        return "fragment_level_only"
    return "no_localization_evidence"


def _structure_status(row: dict[str, Any]) -> str:
    alternatives = str(row.get("Alternative_Structural_Candidates") or "")
    if _bool(row.get("Structure_Discriminating_Evidence")):
        return "structure_resolved"
    existing = str(row.get("Structure_Ambiguity_Status") or "")
    if alternatives:
        if _bool(row.get("Position_Discriminating_Evidence")):
            return "position_resolved_structure_unresolved"
        return "position_and_structure_unresolved"
    if existing == "no_structural_alternative_identified":
        return "no_structural_alternative"
    return "unknown"


def _level(row: dict[str, Any], has_modified: bool, position_status: str, structure_status: str) -> str:
    if not has_modified:
        if any(_bool(row.get(field)) for field in ("Has_MS2_Precursor_Evidence", "Has_Known_Modification_Candidate", "Has_MS1_Fragment_Evidence")):
            return "fragment_mass_shift_supported"
        return "unsupported"
    if structure_status == "structure_resolved":
        return "structure_isomer_resolved"
    if "structure_unresolved" in structure_status:
        return "structure_isomer_unresolved"
    if position_status == "position_resolved":
        return "position_localized"
    return "modified_fragment_ion_supported"


def _score_confidence(level: str, match_count: int, multiple_assignment: bool) -> tuple[float, str, str]:
    base = {
        "unsupported": 0.0, "fragment_mass_shift_supported": 1.0,
        "modified_fragment_ion_supported": 2.0, "position_localized": 4.0,
        "structure_isomer_unresolved": 4.0, "structure_isomer_resolved": 6.0,
    }[level]
    score = base + min(match_count, 4) * 0.25
    if multiple_assignment:
        score -= 0.5
    if level == "structure_isomer_resolved":
        confidence, priority = "High", "high_identity_review"
    elif level in {"position_localized", "structure_isomer_unresolved"}:
        confidence, priority = "Moderate", "moderate_identity_review"
    elif level == "modified_fragment_ion_supported":
        confidence, priority = "Low", "limited_identity_review"
    else:
        confidence, priority = "Low", "low_identity_information"
    return score, confidence, priority


def build_ms2_modification_identity(
    ranking_rows: list[dict[str, Any]], modified_ion_matches: list[dict[str, Any]] | None,
    localization_rows: list[dict[str, Any]] | None,
    ambiguity_groups: list[dict[str, Any]] | None, enabled: bool = True,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    matches_by_key: dict[tuple[str, str, int | None], list[dict[str, Any]]] = defaultdict(list)
    for match in modified_ion_matches or []:
        matches_by_key[_key(match, "Candidate_Modification_Position_In_Parent")].append(match)
    localization_by_key: dict[tuple[str, str, int | None], list[dict[str, Any]]] = defaultdict(list)
    for item in localization_rows or []:
        localization_by_key[_key(item, "Candidate_Modification_Position_In_Parent")].append(item)
    groups_by_id = {str(item.get("Ambiguity_Group_ID") or ""): item for item in ambiguity_groups or []}

    enriched: list[dict[str, Any]] = []
    identity_rows: list[dict[str, Any]] = []
    for original in ranking_rows or []:
        row = dict(original)
        candidate_key = _key(row, "Candidate_Position_In_Parent")
        matches = _deduplicate_matches(matches_by_key.get(candidate_key, [])) if enabled else []
        localization = localization_by_key.get(candidate_key, [])
        if localization and not _bool(row.get("Has_Localization_Evidence")):
            row["Has_Localization_Evidence"] = True
        has_modified = bool(matches)
        match_ids = sorted(_match_id(item) for item in matches)
        unique_ions = {
            (str(item.get("Spectrum_ID") or ""), str(item.get("Ion_ID") or ""))
            for item in matches
        }
        series = sorted({str(item.get("Ion_Type") or "") for item in matches if item.get("Ion_Type")})
        ppm_values = [abs(value) for item in matches if (value := _float(item.get("Mass_Error_ppm"))) is not None]
        intensities = [value for item in matches if (value := _float(item.get("Observed_Intensity"))) is not None]
        multiple_assignment = _multiple_assignment_warning(matches)
        group = groups_by_id.get(str(row.get("Ambiguity_Group_ID") or ""))
        position_status = _position_status(row, group, has_modified)
        structure_status = _structure_status(row)
        level = _level(row, has_modified, position_status, structure_status) if enabled else "unsupported"
        score, confidence, priority = _score_confidence(level, len(matches), multiple_assignment)
        warnings = []
        if not enabled:
            warnings.append("MS2 annotation disabled")
        if multiple_assignment:
            warnings.append("same observed peak has multiple theoretical ion assignments")
        if "structure_unresolved" in structure_status:
            warnings.append("modified fragment ions do not distinguish structural isomers")
        if position_status in {"adjacent_positions_ambiguous", "multiple_positions_ambiguous"}:
            warnings.append("modified fragment ions do not uniquely localize the candidate position")
        if not has_modified:
            warnings.append("no modified fragment ion match")
        reasons = [
            f"modified fragment matches={len(matches)}",
            f"unique modified ions={len(unique_ions)}",
            f"position status={position_status}",
            f"structure status={structure_status}",
        ]
        shadow = {
            "Has_Modified_Fragment_Ion_Evidence": has_modified,
            "Modified_Fragment_Match_Count": len(matches),
            "Unique_Modified_Fragment_Ion_Count": len(unique_ions),
            "Modified_Fragment_Ion_Series": ";".join(series),
            "Supporting_Modified_Fragment_Match_IDs": ";".join(match_ids),
            "Best_Modified_Fragment_Error_ppm": min(ppm_values) if ppm_values else "",
            "Maximum_Modified_Fragment_Intensity": max(intensities) if intensities else "",
            "Position_Localization_Status": position_status,
            "Position_Discriminating_Ion_Count": max(
                [int(_float(row.get("Num_Position_Discriminating_Ions"), 0) or 0)]
                + [int(_float(item.get("Num_Position_Discriminating_Modified_Ions"), 0) or 0) for item in localization]
            ),
            "Structural_Isomer_Group_ID": row.get("Structural_Isomer_Group_ID", ""),
            "Structure_Resolution_Status": structure_status,
            "Alternative_Modification_IDs": row.get("Alternative_Structural_Candidates", ""),
            "MS2_Identity_Evidence_Level": level,
            "Shadow_MS2_Identity_Score": score,
            "Shadow_MS2_Identity_Confidence": confidence,
            "Shadow_MS2_Identity_Priority": priority,
            "MS2_Identity_Evidence_Reason": "; ".join(reasons),
            "MS2_Identity_Warnings": "; ".join(warnings),
        }
        row.update(shadow)
        enriched.append(row)
        identity_rows.append({
            "Rank": row.get("Rank"), "Modification_ID": row.get("Modification_ID"),
            "Modification_Name": row.get("Modification_Name"),
            "Parent_Fragment_ID": row.get("Parent_Fragment_ID"),
            "Candidate_tRNA_Position": row.get("Candidate_tRNA_Position"),
            "Candidate_Base": row.get("Candidate_Base"), **shadow,
        })
    return enriched, identity_rows
