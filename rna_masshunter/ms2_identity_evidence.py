"""Candidate-level shadow summary of existing MS/MS identity evidence.

Existing match rows are reused without rematching. Candidate assignments are
retained, while scoring and ambiguity use candidate-crossing physical peaks.
Formal ranking, confidence, localization, biology, and inclusion are unchanged.
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
    "Maximum_Modified_Fragment_Intensity", "Physical_Observed_Peak_Keys",
    "Shared_Physical_Peak_Count", "Unique_Physical_Peak_Count",
    "Candidate_Specific_Physical_Peak_Count", "Isomer_Group_Shared_Peak_Count",
    "Has_Cross_Candidate_Peak_Sharing", "Cross_Candidate_Peak_Sharing_Warning",
    "Candidate_Specific_Evidence_Peak_Count", "Group_Shared_Evidence_Peak_Count",
    "Cross_Candidate_Ambiguous_Peak_Count", "Identity_Evidence_Scope",
    "Position_Localization_Status", "Group_Position_Resolution_Status",
    "Candidate_Position_Resolution_Status", "Position_Resolution_Ceiling_Applied",
    "Position_Resolution_Caveat", "Position_Discriminating_Ion_Count",
    "Structural_Isomer_Group_ID", "Structure_Resolution_Status",
    "Alternative_Modification_IDs", "MS2_Identity_Evidence_Level",
    "Shadow_MS2_Identity_Score", "Shadow_MS2_Identity_Confidence",
    "Shadow_MS2_Identity_Priority", "MS2_Identity_Evidence_Reason",
    "MS2_Identity_Warnings",
]
# Structural_Isomer_Group_ID already belongs to the biological shadow columns.
IDENTITY_SHADOW_COLUMNS = [column for column in IDENTITY_COLUMNS[6:] if column != "Structural_Isomer_Group_ID"]

PEAK_ASSIGNMENT_COLUMNS = [
    "Physical_Observed_Peak_Key", "Match_ID", "Spectrum_ID", "RT", "Observed_mz",
    "Observed_Intensity", "Modification_ID", "Parent_Fragment_ID",
    "Candidate_Position_In_Parent", "Candidate_tRNA_Position", "Theoretical_Ion_ID",
    "Theoretical_mz", "Structural_Isomer_Group_ID", "Physical_Peak_Assignment_Count",
    "Physical_Peak_Candidate_Count", "Physical_Peak_Theoretical_Ion_Count",
    "Physical_Peak_Shared_Across_Candidates", "Physical_Peak_Shared_Across_Isomers",
    "Physical_Peak_Assignment_Status", "Evidence_Scope",
    "Counts_For_Individual_Identity", "Assignment_Warning",
]


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
    return (str(row.get("Modification_ID") or ""), str(row.get("Parent_Fragment_ID") or ""), _position(row.get(position_name)))


def _match_id(row: dict[str, Any]) -> str:
    values = (row.get("Spectrum_ID"), row.get("Scan_Index"), row.get("Observed_mz"), row.get("Ion_ID"), row.get("Theoretical_mz"))
    return ":".join(str(value if value not in (None, "") else "NA") for value in values)


def physical_observed_peak_key(row: dict[str, Any], mz_decimals: int = 8, rt_decimals: int = 6) -> str:
    """Return a generic deterministic physical peak key without theoretical-ion metadata."""
    existing = row.get("Observed_Peak_ID") or row.get("Observed_Peak_Index") or row.get("Peak_ID")
    spectrum = str(row.get("Spectrum_ID") or "unknown_spectrum")
    if existing not in (None, ""):
        return f"{spectrum}|peak={existing}"
    mz = _float(row.get("Observed_mz"))
    rt = _float(row.get("RT"))
    mz_text = f"{mz:.{mz_decimals}f}" if mz is not None else "NA"
    rt_text = f"{rt:.{rt_decimals}f}" if rt is not None else "NA"
    return f"{spectrum}|mz={mz_text}|rt={rt_text}"


def _deduplicate_all_matches(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    unique: dict[tuple[tuple[str, str, int | None], str], dict[str, Any]] = {}
    for row in rows:
        if not _bool(row.get("Ion_Contains_Modification")):
            continue
        candidate = _key(row, "Candidate_Modification_Position_In_Parent")
        unique.setdefault((candidate, _match_id(row)), row)
    return list(unique.values())


def _ambiguity_positions(group: dict[str, Any] | None) -> list[int]:
    values = str((group or {}).get("Candidate_Positions_In_tRNA") or "").split(";")
    return sorted({position for value in values if (position := _position(value)) is not None})


def _group_status(group: dict[str, Any] | None) -> str:
    return str((group or {}).get("Position_Ambiguity_Status") or "unknown")


def _position_status(
    row: dict[str, Any], group: dict[str, Any] | None, has_modified: bool,
    discriminating_count: int, counterpart_count: int,
) -> tuple[str, bool, str]:
    group_status = _group_status(group)
    positions = _ambiguity_positions(group)
    candidate_flag = _bool(row.get("Position_Discriminating_Evidence")) and discriminating_count > 0
    if group_status in {"resolved", "resolved_by_discriminating_ions"} and candidate_flag:
        if discriminating_count == 1 and counterpart_count <= 0:
            return "position_supported_single_ion", True, "single modified ion without paired unmodified counterpart"
        return "position_resolved", False, "group and candidate evidence support resolution"
    if group_status == "partially_resolved":
        return "partially_resolved", candidate_flag, "ambiguity group is only partially resolved"
    if len(positions) > 1:
        status = "adjacent_positions_ambiguous" if positions[-1] - positions[0] == len(positions) - 1 else "multiple_positions_ambiguous"
        return status, candidate_flag, f"ambiguity group status={group_status}"
    if candidate_flag:
        return "position_supported_single_ion", True, f"candidate flag exceeds group status={group_status}"
    if has_modified or _bool(row.get("Has_Localization_Evidence")):
        return "fragment_level_only", False, "no group-level position resolution"
    return "no_localization_evidence", False, "no localization evidence"


def _structure_status(row: dict[str, Any], position_status: str) -> str:
    alternatives = str(row.get("Alternative_Structural_Candidates") or "")
    if _bool(row.get("Structure_Discriminating_Evidence")):
        return "structure_resolved"
    if alternatives:
        if position_status in {"position_resolved", "position_supported_single_ion"}:
            return "position_resolved_structure_unresolved"
        return "position_and_structure_unresolved"
    if str(row.get("Structure_Ambiguity_Status") or "") == "no_structural_alternative_identified":
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


def _merge_shadow(row: dict[str, Any], shadow: dict[str, Any]) -> None:
    canonical = "Structural_Isomer_Group_ID"
    old, new = row.get(canonical), shadow.get(canonical)
    if old not in (None, "") and new not in (None, "") and str(old) != str(new):
        raise ValueError(f"Conflicting {canonical}: existing={old!r}, shadow={new!r}")
    row.update(shadow)


def build_ms2_modification_identity(
    ranking_rows: list[dict[str, Any]], modified_ion_matches: list[dict[str, Any]] | None,
    localization_rows: list[dict[str, Any]] | None,
    ambiguity_groups: list[dict[str, Any]] | None, enabled: bool = True,
    return_assignments: bool = False,
):
    ranking_by_key = {_key(row, "Candidate_Position_In_Parent"): row for row in ranking_rows or []}
    all_matches = _deduplicate_all_matches(modified_ion_matches or []) if enabled else []
    matches_by_key: dict[tuple[str, str, int | None], list[dict[str, Any]]] = defaultdict(list)
    physical_groups: dict[str, list[tuple[tuple[str, str, int | None], dict[str, Any]]]] = defaultdict(list)
    for match in all_matches:
        candidate = _key(match, "Candidate_Modification_Position_In_Parent")
        matches_by_key[candidate].append(match)
        physical_groups[physical_observed_peak_key(match)].append((candidate, match))

    peak_meta: dict[str, dict[str, Any]] = {}
    for peak_key, assignments in physical_groups.items():
        candidate_keys = {candidate for candidate, _ in assignments}
        theoretical = {(str(match.get("Ion_ID") or ""), str(match.get("Theoretical_mz") or "")) for _, match in assignments}
        group_ids = {str(ranking_by_key.get(candidate, {}).get("Structural_Isomer_Group_ID") or "") for candidate in candidate_keys}
        group_counts: dict[str, int] = defaultdict(int)
        for candidate in candidate_keys:
            group_id = str(ranking_by_key.get(candidate, {}).get("Structural_Isomer_Group_ID") or "")
            if group_id:
                group_counts[group_id] += 1
        shared_group_ids = {group_id for group_id, count in group_counts.items() if count > 1}
        nonempty_groups = {group for group in group_ids if group}
        same_mod_parent = {(candidate[0], candidate[1]) for candidate in candidate_keys}
        isomer_shared = bool(shared_group_ids)
        exclusively_one_isomer_group = len(candidate_keys) > 1 and len(nonempty_groups) == 1 and len(group_ids) == 1
        if len(candidate_keys) == 1 and len(assignments) == 1:
            status, scope = "unique_candidate_assignment", "candidate_specific"
        elif len(candidate_keys) == 1:
            status, scope = "shared_within_same_candidate", "candidate_specific"
        elif exclusively_one_isomer_group:
            status, scope = "shared_across_structural_isomers", "structural_isomer_group_level"
        elif len(same_mod_parent) == 1:
            status, scope = "shared_across_candidates", "position_group_level"
        elif len(candidate_keys) > 1 and isomer_shared:
            status, scope = "ambiguous_peak_identity", "cross_candidate_ambiguous"
        elif len(candidate_keys) > 1:
            status, scope = "shared_across_candidates", "cross_candidate_ambiguous"
        else:
            status, scope = "ambiguous_peak_identity", "cross_candidate_ambiguous"
        peak_meta[peak_key] = {
            "assignments": len(assignments), "candidates": len(candidate_keys),
            "theoretical": len(theoretical), "candidate_keys": candidate_keys,
            "isomer_shared": isomer_shared, "shared_group_ids": shared_group_ids,
            "status": status, "scope": scope,
        }

    localization_by_key: dict[tuple[str, str, int | None], list[dict[str, Any]]] = defaultdict(list)
    for item in localization_rows or []:
        localization_by_key[_key(item, "Candidate_Modification_Position_In_Parent")].append(item)
    groups_by_id = {str(item.get("Ambiguity_Group_ID") or ""): item for item in ambiguity_groups or []}

    enriched: list[dict[str, Any]] = []
    identity_rows: list[dict[str, Any]] = []
    for original in ranking_rows or []:
        row = dict(original)
        candidate_key = _key(row, "Candidate_Position_In_Parent")
        matches = matches_by_key.get(candidate_key, [])
        localization = localization_by_key.get(candidate_key, [])
        has_modified = bool(matches)
        match_ids = sorted(_match_id(item) for item in matches)
        peak_keys = sorted({physical_observed_peak_key(item) for item in matches})
        candidate_specific = [key for key in peak_keys if peak_meta[key]["scope"] == "candidate_specific"]
        group_shared = [
            key for key in peak_keys
            if peak_meta[key]["scope"] in {"structural_isomer_group_level", "position_group_level"}
            or (str(row.get("Structural_Isomer_Group_ID") or "") in peak_meta[key]["shared_group_ids"] and str(row.get("Structural_Isomer_Group_ID") or ""))
        ]
        cross_ambiguous = [key for key in peak_keys if peak_meta[key]["scope"] == "cross_candidate_ambiguous"]
        shared = [key for key in peak_keys if peak_meta[key]["candidates"] > 1]
        structural_group_for_candidate = str(row.get("Structural_Isomer_Group_ID") or "")
        isomer_shared = [
            key for key in peak_keys
            if structural_group_for_candidate and structural_group_for_candidate in peak_meta[key]["shared_group_ids"]
        ]
        if cross_ambiguous:
            evidence_scope = "cross_candidate_ambiguous"
        elif isomer_shared:
            evidence_scope = "structural_isomer_group_level"
        elif group_shared:
            evidence_scope = "position_group_level"
        elif candidate_specific:
            evidence_scope = "candidate_specific"
        else:
            evidence_scope = "no_modified_peak_evidence"

        unique_ions = {(str(item.get("Spectrum_ID") or ""), str(item.get("Ion_ID") or "")) for item in matches}
        series = sorted({str(item.get("Ion_Type") or "") for item in matches if item.get("Ion_Type")})
        ppm_values = [abs(value) for item in matches if (value := _float(item.get("Mass_Error_ppm"))) is not None]
        intensities = [value for item in matches if (value := _float(item.get("Observed_Intensity"))) is not None]
        group = groups_by_id.get(str(row.get("Ambiguity_Group_ID") or ""))
        discriminating_count = max(
            [int(_float(row.get("Num_Position_Discriminating_Ions"), 0) or 0)]
            + [int(_float(item.get("Num_Position_Discriminating_Modified_Ions"), 0) or 0) for item in localization]
        )
        counterpart_count = max([int(_float(item.get("Num_Unmodified_Counterpart_Matches"), 0) or 0) for item in localization] or [0])
        position_status, ceiling_applied, position_caveat = _position_status(
            row, group, has_modified, discriminating_count, counterpart_count
        )
        structure_status = _structure_status(row, position_status)
        level = _level(row, has_modified, position_status, structure_status) if enabled else "unsupported"

        score = 1.0 if level == "fragment_mass_shift_supported" else 0.0
        if has_modified:
            score = 2.0 + min(len(candidate_specific), 4) * 0.5
            if position_status == "position_resolved":
                score += 1.0
            elif position_status == "position_supported_single_ion":
                score += 0.25
            if structure_status == "structure_resolved":
                score += 2.0
        moderate_gate = (
            bool(candidate_specific) and not shared and not cross_ambiguous
            and position_status == "position_resolved"
            and structure_status in {"no_structural_alternative", "structure_resolved"}
        )
        confidence = "Moderate" if moderate_gate else "Low"
        priority = "moderate_identity_review" if moderate_gate else (
            "limited_identity_review" if has_modified else "low_identity_information"
        )

        warnings = []
        if not enabled:
            warnings.append("MS2 annotation disabled")
        if any(peak_meta[key]["assignments"] > 1 for key in peak_keys):
            warnings.append("same physical observed peak has multiple theoretical ion assignments")
        if shared:
            warnings.append("physical observed peak is shared across candidates")
        if isomer_shared:
            warnings.append("physical observed peak is group-level evidence shared across structural isomers")
        if cross_ambiguous:
            warnings.append("physical observed peak has unrelated candidate assignments")
        if position_caveat:
            warnings.append(position_caveat)
        if "structure_unresolved" in structure_status:
            warnings.append("modified fragment ions do not distinguish structural isomers")
        if position_status in {"adjacent_positions_ambiguous", "multiple_positions_ambiguous", "partially_resolved"}:
            warnings.append("modified fragment ions do not fully resolve the candidate position")
        if not has_modified:
            warnings.append("no modified fragment ion match")

        structural_group = row.get("Structural_Isomer_Group_ID", "")
        shadow = {
            "Has_Modified_Fragment_Ion_Evidence": has_modified,
            "Modified_Fragment_Match_Count": len(matches),
            "Unique_Modified_Fragment_Ion_Count": len(unique_ions),
            "Modified_Fragment_Ion_Series": ";".join(series),
            "Supporting_Modified_Fragment_Match_IDs": ";".join(match_ids),
            "Best_Modified_Fragment_Error_ppm": min(ppm_values) if ppm_values else "",
            "Maximum_Modified_Fragment_Intensity": max(intensities) if intensities else "",
            "Physical_Observed_Peak_Keys": ";".join(peak_keys),
            "Shared_Physical_Peak_Count": len(shared),
            "Unique_Physical_Peak_Count": len(peak_keys),
            "Candidate_Specific_Physical_Peak_Count": len(candidate_specific),
            "Isomer_Group_Shared_Peak_Count": len(isomer_shared),
            "Has_Cross_Candidate_Peak_Sharing": bool(shared),
            "Cross_Candidate_Peak_Sharing_Warning": "physical observed peak shared across candidate assignments" if shared else "",
            "Candidate_Specific_Evidence_Peak_Count": len(candidate_specific),
            "Group_Shared_Evidence_Peak_Count": len(group_shared),
            "Cross_Candidate_Ambiguous_Peak_Count": len(cross_ambiguous),
            "Identity_Evidence_Scope": evidence_scope,
            "Position_Localization_Status": position_status,
            "Group_Position_Resolution_Status": _group_status(group),
            "Candidate_Position_Resolution_Status": position_status,
            "Position_Resolution_Ceiling_Applied": ceiling_applied,
            "Position_Resolution_Caveat": position_caveat,
            "Position_Discriminating_Ion_Count": discriminating_count,
            "Structural_Isomer_Group_ID": structural_group,
            "Structure_Resolution_Status": structure_status,
            "Alternative_Modification_IDs": row.get("Alternative_Structural_Candidates", ""),
            "MS2_Identity_Evidence_Level": level,
            "Shadow_MS2_Identity_Score": score,
            "Shadow_MS2_Identity_Confidence": confidence,
            "Shadow_MS2_Identity_Priority": priority,
            "MS2_Identity_Evidence_Reason": (
                f"candidate matches={len(matches)}; physical peaks={len(peak_keys)}; "
                f"candidate-specific={len(candidate_specific)}; group-shared={len(group_shared)}; "
                f"cross-ambiguous={len(cross_ambiguous)}; position={position_status}; structure={structure_status}"
            ),
            "MS2_Identity_Warnings": "; ".join(warnings),
        }
        _merge_shadow(row, shadow)
        enriched.append(row)
        identity_rows.append({
            "Rank": row.get("Rank"), "Modification_ID": row.get("Modification_ID"),
            "Modification_Name": row.get("Modification_Name"),
            "Parent_Fragment_ID": row.get("Parent_Fragment_ID"),
            "Candidate_tRNA_Position": row.get("Candidate_tRNA_Position"),
            "Candidate_Base": row.get("Candidate_Base"), **shadow,
        })

    assignment_rows: list[dict[str, Any]] = []
    for peak_key, assignments in sorted(physical_groups.items()):
        meta = peak_meta[peak_key]
        warning = "" if meta["status"] == "unique_candidate_assignment" else f"physical peak assignment status={meta['status']}"
        for candidate, match in assignments:
            ranking = ranking_by_key.get(candidate, {})
            assignment_rows.append({
                "Physical_Observed_Peak_Key": peak_key, "Match_ID": _match_id(match),
                "Spectrum_ID": match.get("Spectrum_ID"), "RT": match.get("RT"),
                "Observed_mz": match.get("Observed_mz"), "Observed_Intensity": match.get("Observed_Intensity"),
                "Modification_ID": match.get("Modification_ID"), "Parent_Fragment_ID": match.get("Parent_Fragment_ID"),
                "Candidate_Position_In_Parent": candidate[2], "Candidate_tRNA_Position": ranking.get("Candidate_tRNA_Position", ""),
                "Theoretical_Ion_ID": match.get("Ion_ID"), "Theoretical_mz": match.get("Theoretical_mz"),
                "Structural_Isomer_Group_ID": ranking.get("Structural_Isomer_Group_ID", ""),
                "Physical_Peak_Assignment_Count": meta["assignments"],
                "Physical_Peak_Candidate_Count": meta["candidates"],
                "Physical_Peak_Theoretical_Ion_Count": meta["theoretical"],
                "Physical_Peak_Shared_Across_Candidates": meta["candidates"] > 1,
                "Physical_Peak_Shared_Across_Isomers": meta["isomer_shared"],
                "Physical_Peak_Assignment_Status": meta["status"], "Evidence_Scope": meta["scope"],
                "Counts_For_Individual_Identity": meta["scope"] == "candidate_specific",
                "Assignment_Warning": warning,
            })

    if return_assignments:
        return enriched, identity_rows, assignment_rows
    return enriched, identity_rows
