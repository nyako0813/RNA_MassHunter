"""Shadow-only synthesis of existing RNase MS/MS evidence tables.

The builder is a pure table-integration layer. It does not read mzML, inspect
spectra, consume configuration, rematch ions, or alter formal results.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Any


FORMAL_FALSE = {
    "Applied_To_Formal_Result": False,
    "Formal_Change_Ready": False,
    "Formal_Result_Changed": False,
}

SUMMARY_COLUMNS = [
    "Candidate_Count", "Peak_Evidence_Row_Count", "Ranking_Row_Count",
    "Ambiguity_Group_Count", "Modified_Precursor_Row_Count",
    "Modified_Theoretical_Ion_Row_Count", "Modified_Ion_Match_Row_Count",
    "Localization_Row_Count", "Identity_Row_Count",
    "Identity_Peak_Assignment_Row_Count", "Ambiguous_Cluster_Count",
    "Ambiguous_Peak_Detail_Row_Count", "Effective_Ambiguity_Row_Count",
    "Effective_Ambiguity_Detail_Row_Count",
    "Candidate_Specific_Fragment_Evidence_Count",
    "Precursor_Compatible_Only_Count", "Localization_Review_Required_Count",
    "Existing_Ambiguity_Present_Count", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]

CANDIDATE_EVIDENCE_COLUMNS = [
    "Candidate_Key", "Modification_ID", "Modification_Name",
    "Parent_Fragment_ID", "Candidate_Position_In_Parent",
    "Candidate_tRNA_Position", "Candidate_Base", "Spectrum_IDs",
    "Modification_Identity_Status", "Localization_Status", "Structure_Status",
    "Ambiguity_Status", "Has_Precursor_Compatibility",
    "Modified_Precursor_Count", "Modified_Fragment_Match_Count",
    "Candidate_Specific_Physical_Peak_Count", "Shared_Physical_Peak_Count",
    "Position_Discriminating_Ion_Count", "Competing_Position_Count",
    "Competing_Theoretical_Ion_Count", "Existing_Localization_Status",
    "Existing_Structure_Status", "Existing_Ambiguity_Status",
    "Evidence_Basis", "Limiting_Reasons", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]

PEAK_EVIDENCE_COLUMNS = [
    "Physical_Observed_Peak_Key", "Match_ID", "Spectrum_ID", "RT",
    "Observed_mz", "Observed_Intensity", "Modification_ID",
    "Parent_Fragment_ID", "Candidate_Position_In_Parent",
    "Candidate_tRNA_Position", "Theoretical_Ion_ID", "Theoretical_mz",
    "Structural_Isomer_Group_ID", "Physical_Peak_Assignment_Count",
    "Physical_Peak_Candidate_Count", "Physical_Peak_Theoretical_Ion_Count",
    "Physical_Peak_Shared_Across_Candidates",
    "Physical_Peak_Shared_Across_Isomers", "Physical_Peak_Assignment_Status",
    "Evidence_Scope", "Counts_For_Individual_Identity",
    "Ambiguous_Peak_Cluster_IDs", "Effective_Ambiguity_Level",
    "Formal_Match_Ambiguous", "Used_In_Identity", "Used_In_Localization",
    "Assignment_Warning", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]


@dataclass(frozen=True)
class RNaseMS2EvidenceSynthesisResult:
    summary_rows: list[dict[str, Any]]
    candidate_rows: list[dict[str, Any]]
    peak_rows: list[dict[str, Any]]

    @property
    def sheets(self) -> dict[str, list[dict[str, Any]]]:
        return {
            "RNase_MS2_Evidence_Summary": self.summary_rows,
            "RNase_MS2_Candidate_Evidence": self.candidate_rows,
            "RNase_MS2_Peak_Evidence": self.peak_rows,
        }

    def to_jsonable(self) -> dict[str, Any]:
        return {
            "summary": self.summary_rows[0] if self.summary_rows else {},
            "candidates": self.candidate_rows,
            "peaks": self.peak_rows,
        }


def _position(value: Any) -> str:
    if value in (None, ""):
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return str(int(number)) if number.is_integer() else str(number)


def _candidate_key(modification: Any, parent: Any, position: Any) -> tuple[str, str, str]:
    return str(modification or ""), str(parent or ""), _position(position)


def _candidate_text(key: tuple[str, str, str]) -> str:
    return "|".join(value or "NA" for value in key)


def _split(values: Any) -> set[str]:
    return {value for value in str(values or "").split(";") if value}


def _joined(values: set[str]) -> str:
    return ";".join(sorted(value for value in values if value))


def _bool(value: Any) -> bool:
    return value is True or str(value or "").strip().lower() in {"1", "true", "yes", "on"}


def _float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _match_id(row: dict[str, Any]) -> str:
    values = (
        row.get("Spectrum_ID"), row.get("Scan_Index"), row.get("Observed_mz"),
        row.get("Ion_ID"), row.get("Theoretical_mz"),
    )
    return ":".join(
        str(value if value not in (None, "") else "NA") for value in values
    )


def _positive(row: dict[str, Any]) -> bool:
    value = _float(row.get("Observed_Intensity"))
    return value is None or value > 0


def _ambiguity_status(types: set[str]) -> str:
    return "NONE" if not types else sorted(types)[0] if len(types) == 1 else "MULTIPLE"


def _row_sort_key(row: dict[str, Any]) -> tuple[tuple[str, str], ...]:
    return tuple(sorted((str(key), repr(value)) for key, value in row.items()))


def build_rnase_ms2_evidence_synthesis(
    ranking_rows: list[dict[str, Any]],
    ambiguity_groups: list[dict[str, Any]],
    modified_precursors: list[dict[str, Any]],
    modified_theoretical_ions: list[dict[str, Any]],
    modified_ion_matches: list[dict[str, Any]],
    localization_rows: list[dict[str, Any]],
    identity_rows: list[dict[str, Any]],
    identity_peak_assignments: list[dict[str, Any]],
    ambiguous_clusters: list[dict[str, Any]],
    ambiguous_peak_details: list[dict[str, Any]],
    effective_ambiguity_rows: list[dict[str, Any]],
    effective_ambiguity_details: list[dict[str, Any]],
) -> RNaseMS2EvidenceSynthesisResult:
    """Build conservative skeleton statuses from existing rows only."""
    inputs = [
        ranking_rows, ambiguity_groups, modified_precursors,
        modified_theoretical_ions, modified_ion_matches, localization_rows,
        identity_rows, identity_peak_assignments, ambiguous_clusters,
        ambiguous_peak_details, effective_ambiguity_rows,
        effective_ambiguity_details,
    ]
    (
        ranking_rows, ambiguity_groups, modified_precursors,
        modified_theoretical_ions, modified_ion_matches, localization_rows,
        identity_rows, identity_peak_assignments, ambiguous_clusters,
        ambiguous_peak_details, effective_ambiguity_rows,
        effective_ambiguity_details,
    ) = tuple(list(rows or []) for rows in inputs)

    candidates: dict[tuple[str, str, str], dict[str, Any]] = {}

    def ensure(row: dict[str, Any], position_field: str) -> dict[str, Any]:
        key = _candidate_key(
            row.get("Modification_ID"), row.get("Parent_Fragment_ID"),
            row.get(position_field),
        )
        item = candidates.setdefault(key, {
            "key": key, "ranking": [], "identity": [], "precursors": [],
            "ions": [], "matches": [], "localization": [], "assignments": [],
            "ambiguity_groups": [], "effective": [],
        })
        return item

    ranking_parent_position = {
        (
            str(row.get("Modification_ID") or ""),
            str(row.get("Parent_Fragment_ID") or ""),
            _position(row.get("Candidate_tRNA_Position")),
        ): row.get("Candidate_Position_In_Parent")
        for row in ranking_rows
    }
    for row in ranking_rows:
        ensure(row, "Candidate_Position_In_Parent")["ranking"].append(row)
    for row in identity_rows:
        normalized = dict(row)
        normalized["Candidate_Position_In_Parent"] = ranking_parent_position.get(
            (
                str(row.get("Modification_ID") or ""),
                str(row.get("Parent_Fragment_ID") or ""),
                _position(row.get("Candidate_tRNA_Position")),
            ),
            row.get("Candidate_Position_In_Parent", ""),
        )
        ensure(normalized, "Candidate_Position_In_Parent")["identity"].append(row)
    for row in localization_rows:
        ensure(row, "Candidate_Modification_Position_In_Parent")["localization"].append(row)
    for row in modified_ion_matches:
        ensure(row, "Candidate_Modification_Position_In_Parent")["matches"].append(row)
    for row in modified_theoretical_ions:
        ensure(row, "Candidate_Modification_Position_In_Parent")["ions"].append(row)
    for row in identity_peak_assignments:
        ensure(row, "Candidate_Position_In_Parent")["assignments"].append(row)

    positioned_by_mod_parent: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for item in candidates.values():
        positioned_by_mod_parent[(item["key"][0], item["key"][1])].append(item)
    for row in modified_precursors:
        key = (str(row.get("Modification_ID") or ""), str(row.get("Parent_Fragment_ID") or ""))
        targets = positioned_by_mod_parent.get(key)
        if not targets:
            targets = [ensure(row, "Candidate_Modification_Position_In_Parent")]
        for item in targets:
            item["precursors"].append(row)
    for row in ambiguity_groups:
        key = (str(row.get("Modification_ID") or ""), str(row.get("Parent_Fragment_ID") or ""))
        for item in positioned_by_mod_parent.get(key, []):
            item["ambiguity_groups"].append(row)
    for item in candidates.values():
        for name in (
            "ranking", "identity", "precursors", "ions", "matches",
            "localization", "assignments", "ambiguity_groups", "effective",
        ):
            item[name].sort(key=_row_sort_key)

    cluster_ids_by_peak: dict[str, set[str]] = defaultdict(set)
    ambiguous_detail_by_peak: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in ambiguous_peak_details:
        peak = str(row.get("Physical_Observed_Peak_Key") or "")
        if peak:
            cluster_ids_by_peak[peak].add(str(row.get("Peak_Cluster_ID") or ""))
            ambiguous_detail_by_peak[peak].append(row)
    effective_by_peak: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in effective_ambiguity_details:
        peak = str(row.get("Physical_Peak_ID") or "")
        if peak:
            effective_by_peak[peak].append(row)

    candidate_rows: list[dict[str, Any]] = []
    for key in sorted(candidates):
        item = candidates[key]
        ranking = item["ranking"][0] if item["ranking"] else {}
        identity = item["identity"][0] if item["identity"] else {}
        localization = item["localization"][0] if item["localization"] else {}
        assignments = item["assignments"]
        positive_assignments = [
            row for row in assignments if _positive(row)
        ]
        candidate_specific_assignments = [
            row for row in positive_assignments
            if _bool(row.get("Counts_For_Individual_Identity"))
            and str(row.get("Evidence_Scope") or "") == "candidate_specific"
        ]
        shared_assignments = [
            row for row in positive_assignments
            if row not in candidate_specific_assignments
        ]
        candidate_specific = {
            str(row.get("Physical_Observed_Peak_Key") or "")
            for row in candidate_specific_assignments
        }
        shared = {
            str(row.get("Physical_Observed_Peak_Key") or "")
            for row in shared_assignments
        }
        modified_matches = [
            row for row in item["matches"]
            if _bool(row.get("Ion_Contains_Modification")) and _positive(row)
        ]
        modified_by_id = {_match_id(row): row for row in modified_matches}
        specific_modified_assignments = [
            row for row in candidate_specific_assignments
            if str(row.get("Match_ID") or "") in modified_by_id
            or not row.get("Match_ID")
        ]
        position_support_assignments = []
        partial_position_assignments = []
        for assignment in positive_assignments:
            match = modified_by_id.get(str(assignment.get("Match_ID") or ""))
            if match is None:
                continue
            discriminating = _bool(match.get("Position_Discriminating_Ion"))
            informative = _bool(match.get("Informative_Ion"))
            if discriminating:
                partial_position_assignments.append(assignment)
            if (
                discriminating and informative
                and assignment in candidate_specific_assignments
            ):
                position_support_assignments.append(assignment)

        ambiguity_types: set[str] = set()
        relevant_peaks = {
            str(row.get("Physical_Observed_Peak_Key") or "")
            for row in positive_assignments
        }
        for assignment in positive_assignments:
            scope = str(assignment.get("Evidence_Scope") or "")
            if scope == "position_group_level":
                ambiguity_types.add("POSITION_GROUP")
            elif scope == "structural_isomer_group_level" or _bool(
                assignment.get("Physical_Peak_Shared_Across_Isomers")
            ):
                ambiguity_types.add("STRUCTURAL_ISOMER")
            elif scope == "cross_candidate_ambiguous":
                ambiguity_types.add("CROSS_CANDIDATE")
            if int(assignment.get("Physical_Peak_Theoretical_Ion_Count") or 0) > 1:
                ambiguity_types.add("THEORETICAL_ION_COMPETITION")
        for peak in relevant_peaks:
            for detail in ambiguous_detail_by_peak.get(peak, []):
                specificity = str(detail.get("Candidate_Specificity_Status") or "")
                if specificity == "position_group_shared":
                    ambiguity_types.add("POSITION_GROUP")
                elif specificity == "structural_isomer_group_shared":
                    ambiguity_types.add("STRUCTURAL_ISOMER")
                elif specificity == "cross_candidate_shared":
                    ambiguity_types.add("CROSS_CANDIDATE")
                if int(detail.get("Theoretical_Ion_Count") or 0) > 1:
                    ambiguity_types.add("THEORETICAL_ION_COMPETITION")
            for detail in effective_by_peak.get(peak, []):
                level = str(detail.get("Effective_Ambiguity_Level") or "")
                if level == "raw_only":
                    continue
                if _bool(detail.get("Position_Group_Shared")):
                    ambiguity_types.add("POSITION_GROUP")
                if _bool(detail.get("Structural_Isomer_Shared")):
                    ambiguity_types.add("STRUCTURAL_ISOMER")
                if _bool(detail.get("Cross_Candidate_Shared")):
                    ambiguity_types.add("CROSS_CANDIDATE")
                if _bool(detail.get("Formal_Match_Ambiguous")) and not (
                    _bool(detail.get("Position_Group_Shared"))
                    or _bool(detail.get("Structural_Isomer_Shared"))
                    or _bool(detail.get("Cross_Candidate_Shared"))
                ):
                    ambiguity_types.add("THEORETICAL_ION_COMPETITION")
        group_statuses = {
            str(row.get("Position_Ambiguity_Status") or "")
            for row in item["ambiguity_groups"]
            if row.get("Position_Ambiguity_Status")
        }
        if group_statuses & {"ambiguous", "partially_resolved"}:
            ambiguity_types.add("POSITION_GROUP")
        if _bool(identity.get("Has_Cross_Candidate_Peak_Sharing")):
            ambiguity_types.add("CROSS_CANDIDATE")
        structural_status = str(identity.get("Structure_Resolution_Status") or "")
        if "structure_unresolved" in structural_status or _bool(
            identity.get("Isomer_Group_Shared_Peak_Count")
        ):
            ambiguity_types.add("STRUCTURAL_ISOMER")
        ambiguity_status = _ambiguity_status(ambiguity_types)

        has_precursor = bool(item["precursors"])
        has_shared_modified = bool(modified_matches or _bool(
            identity.get("Has_Modified_Fragment_Ion_Evidence")
        )) and not bool(specific_modified_assignments)
        identity_limited = bool(ambiguity_types & {
            "CROSS_CANDIDATE", "THEORETICAL_ION_COMPETITION",
        })
        if specific_modified_assignments and not identity_limited:
            identity_status = "FRAGMENT_SUPPORTED"
        elif specific_modified_assignments or has_shared_modified:
            identity_status = "AMBIGUOUS"
        elif has_precursor:
            identity_status = "PRECURSOR_COMPATIBLE"
        else:
            identity_status = "UNSUPPORTED"

        localization_limited = bool(ambiguity_types & {
            "POSITION_GROUP", "CROSS_CANDIDATE",
            "THEORETICAL_ION_COMPETITION",
        })
        if position_support_assignments and not localization_limited:
            localization_status = "LOCALIZED"
        elif position_support_assignments or (
            partial_position_assignments and localization_limited
        ):
            localization_status = "AMBIGUOUS"
        elif partial_position_assignments or (
            item["localization"]
            and int(localization.get("Num_Modified_Ion_Matches") or 0) > 0
        ):
            localization_status = "PARTIALLY_LOCALIZED"
        else:
            localization_status = "UNRESOLVED"

        if not specific_modified_assignments:
            structure_status = (
                "AMBIGUOUS" if "STRUCTURAL_ISOMER" in ambiguity_types
                and has_shared_modified else "NOT_EVALUATED"
            )
        elif "STRUCTURAL_ISOMER" in ambiguity_types:
            structure_status = "AMBIGUOUS"
        elif structural_status == "structure_resolved" and ambiguity_status == "NONE":
            structure_status = "SUPPORTED"
        else:
            structure_status = "UNRESOLVED"

        spectra = {
            str(row.get("Spectrum_ID") or "")
            for collection in (
                item["precursors"], item["matches"], item["localization"],
                positive_assignments,
            )
            for row in collection
            if row.get("Spectrum_ID")
        }
        competing_positions = {
            value for row in item["ambiguity_groups"]
            for value in _split(row.get("Candidate_Positions_In_Parent"))
        }
        competing_theories = {
            str(row.get("Theoretical_Ion_ID") or "")
            for row in positive_assignments if row.get("Theoretical_Ion_ID")
        }
        for peak in relevant_peaks:
            for detail in ambiguous_detail_by_peak.get(peak, []):
                competing_theories.update(
                    _split(detail.get("Competing_Theoretical_Ion_IDs"))
                )
        basis = []
        if has_precursor:
            basis.append("PRECURSOR_COMPATIBLE")
        if specific_modified_assignments:
            basis.append("CANDIDATE_SPECIFIC_MODIFIED_FRAGMENT_PEAK")
        if position_support_assignments:
            basis.append("INFORMATIVE_CANDIDATE_SPECIFIC_POSITION_DISCRIMINATING_ION")
        limits = []
        if shared:
            limits.append("SHARED_PHYSICAL_PEAKS_EXCLUDED_FROM_INDIVIDUAL_SUPPORT")
        if has_precursor and not specific_modified_assignments:
            limits.append("PRECURSOR_COMPATIBILITY_DOES_NOT_ESTABLISH_IDENTITY_OR_LOCALIZATION")
        if "single_candidate_position" in group_statuses and not position_support_assignments:
            limits.append("SINGLE_CANDIDATE_POSITION_DOES_NOT_ESTABLISH_LOCALIZATION")
        if "THEORETICAL_ION_COMPETITION" in ambiguity_types:
            limits.append("COMPETING_THEORETICAL_IONS_LIMIT_IDENTITY_AND_LOCALIZATION")
        if not limits:
            limits.append("NONE")

        candidate_rows.append({
            "Candidate_Key": _candidate_text(key),
            "Modification_ID": key[0],
            "Modification_Name": (
                ranking.get("Modification_Name")
                or identity.get("Modification_Name")
                or localization.get("Modification_Name") or ""
            ),
            "Parent_Fragment_ID": key[1],
            "Candidate_Position_In_Parent": key[2],
            "Candidate_tRNA_Position": (
                ranking.get("Candidate_tRNA_Position")
                or identity.get("Candidate_tRNA_Position")
                or localization.get("Candidate_Modification_Position_In_tRNA")
                or ""
            ),
            "Candidate_Base": (
                ranking.get("Candidate_Base")
                or identity.get("Candidate_Base")
                or localization.get("Candidate_Modification_Base") or ""
            ),
            "Spectrum_IDs": _joined(spectra),
            "Modification_Identity_Status": identity_status,
            "Localization_Status": localization_status,
            "Structure_Status": structure_status,
            "Ambiguity_Status": ambiguity_status,
            "Has_Precursor_Compatibility": has_precursor,
            "Modified_Precursor_Count": len(item["precursors"]),
            "Modified_Fragment_Match_Count": len(modified_matches),
            "Candidate_Specific_Physical_Peak_Count": len(candidate_specific),
            "Shared_Physical_Peak_Count": len(shared),
            "Position_Discriminating_Ion_Count": len(position_support_assignments),
            "Competing_Position_Count": len(competing_positions),
            "Competing_Theoretical_Ion_Count": len(competing_theories),
            "Existing_Localization_Status": (
                identity.get("Candidate_Position_Resolution_Status")
                or localization.get("Localization_Interpretation") or ""
            ),
            "Existing_Structure_Status": structural_status,
            "Existing_Ambiguity_Status": _joined(group_statuses),
            "Evidence_Basis": ";".join(basis),
            "Limiting_Reasons": ";".join(limits),
            **FORMAL_FALSE,
        })

    peak_rows: list[dict[str, Any]] = []
    for assignment in sorted(
        identity_peak_assignments,
        key=lambda row: (
            str(row.get("Physical_Observed_Peak_Key") or ""),
            str(row.get("Match_ID") or ""),
        ),
    ):
        peak = str(assignment.get("Physical_Observed_Peak_Key") or "")
        effective = effective_by_peak.get(peak, [])
        peak_rows.append({
            **{
                column: assignment.get(column, "")
                for column in PEAK_EVIDENCE_COLUMNS[:21]
            },
            "Ambiguous_Peak_Cluster_IDs": _joined(cluster_ids_by_peak.get(peak, set())),
            "Effective_Ambiguity_Level": _joined({
                str(row.get("Effective_Ambiguity_Level") or "") for row in effective
            }),
            "Formal_Match_Ambiguous": any(
                _bool(row.get("Formal_Match_Ambiguous")) for row in effective
            ),
            "Used_In_Identity": any(
                _bool(row.get("Used_In_Identity")) for row in effective
            ),
            "Used_In_Localization": any(
                _bool(row.get("Used_In_Localization")) for row in effective
            ),
            "Assignment_Warning": assignment.get("Assignment_Warning", ""),
            **FORMAL_FALSE,
        })

    summary = {
        "Candidate_Count": len(candidate_rows),
        "Peak_Evidence_Row_Count": len(peak_rows),
        "Ranking_Row_Count": len(ranking_rows),
        "Ambiguity_Group_Count": len(ambiguity_groups),
        "Modified_Precursor_Row_Count": len(modified_precursors),
        "Modified_Theoretical_Ion_Row_Count": len(modified_theoretical_ions),
        "Modified_Ion_Match_Row_Count": len(modified_ion_matches),
        "Localization_Row_Count": len(localization_rows),
        "Identity_Row_Count": len(identity_rows),
        "Identity_Peak_Assignment_Row_Count": len(identity_peak_assignments),
        "Ambiguous_Cluster_Count": len(ambiguous_clusters),
        "Ambiguous_Peak_Detail_Row_Count": len(ambiguous_peak_details),
        "Effective_Ambiguity_Row_Count": len(effective_ambiguity_rows),
        "Effective_Ambiguity_Detail_Row_Count": len(effective_ambiguity_details),
        "Candidate_Specific_Fragment_Evidence_Count": sum(
            row["Candidate_Specific_Physical_Peak_Count"] > 0
            for row in candidate_rows
        ),
        "Precursor_Compatible_Only_Count": sum(
            row["Modification_Identity_Status"] == "PRECURSOR_COMPATIBLE"
            for row in candidate_rows
        ),
        "Localization_Review_Required_Count": sum(
            row["Localization_Status"] == "PARTIALLY_LOCALIZED"
            for row in candidate_rows
        ),
        "Existing_Ambiguity_Present_Count": sum(
            row["Ambiguity_Status"] != "NONE"
            for row in candidate_rows
        ),
        **FORMAL_FALSE,
    }
    return RNaseMS2EvidenceSynthesisResult([summary], candidate_rows, peak_rows)
