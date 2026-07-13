from __future__ import annotations

from typing import Any

import pandas as pd

from rna_masshunter.ms2_unmatched_audit import TOP_SHADOW_COLUMNS, primary_unmatched_reason
from rna_masshunter.ms2_ambiguous_peak_audit import TOP_SHADOW_COLUMNS as AMBIGUOUS_PEAK_TOP_COLUMNS


DASHBOARD_NOTE = (
    "Review dashboard prioritizes candidates for manual review; it does not confirm modification identity or position."
)

TOP_CANDIDATE_COLUMNS = [
    "Review_Rank",
    "Review_Priority",
    "Review_Category",
    "Modification_ID",
    "Modification_Name",
    "Parent_Fragment_ID",
    "Parent_Sequence",
    "Candidate_Positions_In_tRNA",
    "Candidate_Positions_In_Parent",
    "Best_Final_Score",
    "Best_Final_Confidence",
    "Best_Biological_Context_Score",
    "Has_MS2_Precursor_Evidence",
    "Has_Modified_Ion_Evidence",
    "Has_Position_Discriminating_Evidence",
    "Position_Ambiguity_Status",
    "Position_Resolution_Basis",
    "Num_Positions_In_Ambiguity_Group",
    "Ambiguity_Group_ID",
    "Evidence_Summary",
    "Key_Warnings",
    "Recommended_Next_Check",
    "Modification_Family",
    "Position_Class",
    "Position_Prior_Score",
    "Parent_Base_Compatibility",
    "MS2_Localization_Evidence",
    "Structural_Isomer_Group_ID",
    "Structure_Ambiguity_Status",
    "Alternative_Structural_Candidates",
    "Biological_Plausibility_Score",
    "Biological_Plausibility_Level",
    "Shadow_Final_Score",
    "Shadow_Final_Confidence",
    "Shadow_Only",
    "Has_Modified_Fragment_Ion_Evidence",
    "Modified_Fragment_Match_Count",
    "Unique_Modified_Fragment_Ion_Count",
    "Modified_Fragment_Ion_Series",
    "Supporting_Modified_Fragment_Match_IDs",
    "Best_Modified_Fragment_Error_ppm",
    "Maximum_Modified_Fragment_Intensity",
    "Physical_Observed_Peak_Keys",
    "Shared_Physical_Peak_Count",
    "Unique_Physical_Peak_Count",
    "Candidate_Specific_Physical_Peak_Count",
    "Isomer_Group_Shared_Peak_Count",
    "Has_Cross_Candidate_Peak_Sharing",
    "Cross_Candidate_Peak_Sharing_Warning",
    "Candidate_Specific_Evidence_Peak_Count",
    "Group_Shared_Evidence_Peak_Count",
    "Cross_Candidate_Ambiguous_Peak_Count",
    "Identity_Evidence_Scope",
    "Position_Localization_Status",
    "Group_Position_Resolution_Status",
    "Candidate_Position_Resolution_Status",
    "Position_Resolution_Ceiling_Applied",
    "Position_Resolution_Caveat",
    "Position_Discriminating_Ion_Count",
    "Structure_Resolution_Status",
    "Alternative_Modification_IDs",
    "MS2_Identity_Evidence_Level",
    "Shadow_MS2_Identity_Score",
    "Shadow_MS2_Identity_Confidence",
    "Shadow_MS2_Identity_Priority",
    "MS2_Identity_Evidence_Reason",
    "MS2_Identity_Warnings",
    *TOP_SHADOW_COLUMNS,
    *AMBIGUOUS_PEAK_TOP_COLUMNS,
    "Notes",
]

DECISION_COLUMNS = [
    "Review_Rank",
    "Modification_ID",
    "Parent_Fragment_ID",
    "Candidate_Positions_In_tRNA",
    "Review_Priority",
    "Decision_Label",
    "Decision_Text",
    "Evidence_For",
    "Evidence_Against",
    "Ambiguity_Text",
    "Context_Text",
    "Next_Action",
]

CHECKLIST_COLUMNS = [
    "Review_Rank",
    "Modification_ID",
    "Parent_Fragment_ID",
    "Candidate_Positions_In_tRNA",
    "MS1_Fragment_Evidence",
    "Known_Modification_Candidate",
    "MS2_Precursor_Evidence",
    "Modified_Ion_Evidence",
    "Informative_Modified_Ion",
    "c_Ion_Support",
    "y_Ion_Support",
    "Both_c_y_Series",
    "Position_Discriminating_Ion",
    "Ambiguity_Group",
    "Ambiguous_Position",
    "Single_Candidate_Position",
    "Curated_Source",
    "Candidate_Policy_Allows_Mass_Search",
    "Biological_Context_Support",
    "Near_Isobaric_Warning",
    "Confidence_Limiting_Factors",
]

DASHBOARD_COLUMNS = [
    "Total_Ranked_Candidates",
    "Very_High",
    "High",
    "Medium",
    "Low",
    "Very_Low",
    "Total_Ambiguity_Groups",
    "Single_Candidate_Position_Groups",
    "Resolved_By_Discriminating_Ion_Groups",
    "Ambiguous_Groups",
    "Candidates_With_Modified_Ion_Evidence",
    "Candidates_With_Position_Discriminating_Evidence",
    "Candidates_With_Biological_Context_Support",
    "Top_Modification_IDs",
    "Top_Review_Priority_Modifications",
    "Key_Warnings",
    "Notes",
]

CONFIDENCE_ORDER = {"VERY_HIGH": 5, "HIGH": 4, "MEDIUM": 3, "LOW": 2, "VERY_LOW": 1}


def _frame(value: Any) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value.copy()
    if isinstance(value, list):
        return pd.DataFrame(value)
    if isinstance(value, dict):
        return pd.DataFrame([value])
    return pd.DataFrame()


def _first_existing(row: pd.Series, names: list[str], default: Any = "") -> Any:
    lowered = {str(key).lower(): key for key in row.index}
    for name in names:
        key = lowered.get(name.lower())
        if key is not None:
            value = row.get(key)
            if pd.notna(value) and value != "":
                return value
    return default


def _joined_existing(row: pd.Series, names: list[str]) -> str:
    lowered = {str(key).lower(): key for key in row.index}
    values = []
    for name in names:
        key = lowered.get(name.lower())
        if key is None:
            continue
        value = row.get(key)
        if pd.isna(value) or value == "":
            continue
        text = str(value)
        if text not in values:
            values.append(text)
    return "; ".join(values)


def _truthy(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value > 0
    text = str(value).strip().lower()
    return text in {"true", "yes", "y", "1", "present", "supported", "high", "medium"} or (
        bool(text) and text not in {"false", "no", "n", "0", "none", "nan", "missing", "absent"}
    )


def _confidence_value(confidence: Any) -> int:
    text = str(confidence or "").strip().replace("-", "_").replace(" ", "_").upper()
    return CONFIDENCE_ORDER.get(text, 0)


def _float(value: Any, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _contains_any(value: Any, needles: list[str]) -> bool:
    text = str(value or "").lower()
    return any(needle.lower() in text for needle in needles)


def _priority(
    confidence: Any,
    score: float,
    has_precursor: bool,
    has_modified_ion: bool,
    has_discriminating: bool,
    ambiguity_status: str,
    thresholds: dict[str, Any],
) -> tuple[str, str]:
    ambiguous = _contains_any(ambiguity_status, ["ambiguous", "unresolved"])
    single = _contains_any(ambiguity_status, ["single_candidate_position", "single candidate"])
    strong_min = _float(thresholds.get("strong_review_min_score"), 8.0)
    medium_min = _float(thresholds.get("medium_review_min_score"), 5.0)
    weak_min = _float(thresholds.get("weak_review_min_score"), 2.0)
    confidence_rank = _confidence_value(confidence)

    if ambiguous and has_precursor and has_modified_ion and not has_discriminating:
        return "C_ambiguous_review", "ambiguous_position_review"
    if confidence_rank >= CONFIDENCE_ORDER["HIGH"] or (
        score >= strong_min and has_precursor and has_modified_ion and has_discriminating
    ):
        return "A_strong_review", "strong_multi_evidence_review"
    if (
        confidence_rank == CONFIDENCE_ORDER["MEDIUM"]
        or (score >= medium_min and has_precursor and has_modified_ion and (has_discriminating or single))
    ):
        return "B_medium_review", "medium_evidence_review"
    if score >= weak_min or has_precursor or has_modified_ion:
        return "D_weak_review", "limited_evidence_review"
    return "E_low_information", "low_information_review"


def _recommended_next_check(priority: str, has_discriminating: bool, ambiguity_status: str, key_warnings: str) -> str:
    if priority == "C_ambiguous_review":
        return "Inspect MS2_Modified_Ion_Matches and Modification_Ambiguity_Groups."
    if not has_discriminating and _contains_any(ambiguity_status, ["single_candidate_position", "ambiguous", "unresolved"]):
        return "Check whether additional MS2 spectra support position-discriminating ions."
    if _contains_any(key_warnings, ["near-isobaric", "isobaric"]):
        return "Compare near-isobaric alternatives."
    return "Review curated source and biological context columns."


def _build_top_candidates(
    ranking: pd.DataFrame, ambiguity: pd.DataFrame, config: dict[str, Any],
    unmatched_audit_summary: pd.DataFrame | None = None,
    ambiguous_peak_summary: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if ranking.empty:
        return pd.DataFrame(columns=TOP_CANDIDATE_COLUMNS)

    max_candidates = int(config.get("max_top_candidates", 50) or 50)
    thresholds = config.get("review_priority_thresholds") or {}
    audit_frame = unmatched_audit_summary if unmatched_audit_summary is not None else pd.DataFrame()
    audit_lookup: dict[tuple[str, str, str], dict[str, Any]] = {}
    for _, audit_row in audit_frame.iterrows():
        audit_key = (
            str(audit_row.get("Modification_ID") or ""),
            str(audit_row.get("Parent_Fragment_ID") or ""),
            str(audit_row.get("Candidate_tRNA_Position") if pd.notna(audit_row.get("Candidate_tRNA_Position")) else ""),
        )
        audit_lookup[audit_key] = audit_row.to_dict()
    ambiguity_frame = ambiguous_peak_summary if ambiguous_peak_summary is not None else pd.DataFrame()
    ambiguity_lookup: dict[tuple[str, str, str], dict[str, Any]] = {}
    for _, ambiguity_row in ambiguity_frame.iterrows():
        ambiguity_key = (
            str(ambiguity_row.get("Modification_ID") or ""),
            str(ambiguity_row.get("Parent_Fragment_ID") or ""),
            str(ambiguity_row.get("Candidate_tRNA_Position") if pd.notna(ambiguity_row.get("Candidate_tRNA_Position")) else ""),
        )
        ambiguity_lookup[ambiguity_key] = ambiguity_row.to_dict()
    rows: list[dict[str, Any]] = []
    for _, row in ranking.iterrows():
        mod_id = _first_existing(row, ["Modification_ID", "Modification", "modification_id"])
        parent_id = _first_existing(row, ["Parent_Fragment_ID", "Fragment_ID", "parent_fragment_id"])
        group_id = _first_existing(row, ["Ambiguity_Group_ID", "Group_ID", "ambiguity_group_id"])
        final_score = _float(_first_existing(row, ["Final_Score", "Score", "Ranking_Score", "Best_Final_Score"], 0))
        final_confidence = _first_existing(row, ["Final_Confidence", "Confidence", "Best_Final_Confidence"])
        candidate_trna_position = _first_existing(row, ["Candidate_Positions_In_tRNA", "Candidate_tRNA_Position", "Positions_In_tRNA", "tRNA_Position"] )
        audit_key = (str(mod_id or ""), str(parent_id or ""), str(candidate_trna_position if candidate_trna_position != "" else ""))
        audit_summary = audit_lookup.get(audit_key, {})
        ambiguity_summary = ambiguity_lookup.get(audit_key, {})
        context_score = _float(_first_existing(row, ["Biological_Context_Score", "Context_Score", "Best_Biological_Context_Score"], 0))
        ambiguity_status = str(_first_existing(row, ["Position_Ambiguity_Status", "Ambiguity_Status", "Position_Status"], ""))
        resolution_basis = str(_first_existing(row, ["Position_Resolution_Basis", "Resolution_Basis"], ""))
        if not ambiguity_status and group_id and not ambiguity.empty:
            matches = ambiguity[ambiguity.astype(str).eq(str(group_id)).any(axis=1)]
            if not matches.empty:
                ambiguity_status = str(_first_existing(matches.iloc[0], ["Position_Ambiguity_Status", "Ambiguity_Status", "Group_Status"], ""))
                resolution_basis = resolution_basis or str(
                    _first_existing(matches.iloc[0], ["Position_Resolution_Basis", "Resolution_Basis"], "")
                )

        has_precursor = _truthy(_first_existing(row, ["Has_MS2_Precursor_Evidence", "MS2_Precursor_Evidence", "Precursor_Evidence"]))
        has_modified_ion = _truthy(_first_existing(row, ["Has_Modified_Ion_Evidence", "Modified_Ion_Evidence", "Modified_Ion_Count"]))
        has_discriminating = _truthy(
            _first_existing(row, ["Has_Position_Discriminating_Evidence", "Position_Discriminating_Evidence", "Position_Discriminating_Ion", "Discriminating_Ion_Count"])
        )
        priority, category = _priority(
            final_confidence,
            final_score,
            has_precursor,
            has_modified_ion,
            has_discriminating,
            ambiguity_status,
            thresholds,
        )
        warnings = _joined_existing(
            row,
            ["Key_Warnings", "Warnings", "Near_Isobaric_Warning", "Confidence_Limiting_Factors", "Limiting_Factors"],
        )
        structure_status = str(row.get("Structure_Ambiguity_Status") or "")
        if "unresolved" in structure_status:
            structure_warnings = [
                "isobaric structural alternatives remain",
                "modified-ion evidence does not distinguish structural isomers",
            ]
            if has_discriminating:
                structure_warnings.insert(0, "position localized but modification structure unresolved")
            warnings = "; ".join(item for item in [warnings, *structure_warnings] if item)
        evidence_summary = str(_first_existing(row, ["Evidence_Summary", "Summary", "Evidence"], ""))
        if not evidence_summary:
            parts = []
            if has_precursor:
                parts.append("precursor evidence")
            if has_modified_ion:
                parts.append("modified ion evidence")
            if has_discriminating:
                parts.append("position-discriminating evidence")
            evidence_summary = "; ".join(parts) or "No strong evidence flags were present in the ranking table."

        rows.append(
            {
                "Review_Priority": priority,
                "Review_Category": category,
                "Modification_ID": mod_id,
                "Modification_Name": _first_existing(row, ["Modification_Name", "Name", "modification_name"]),
                "Parent_Fragment_ID": parent_id,
                "Parent_Sequence": _first_existing(row, ["Parent_Sequence", "Sequence", "parent_sequence"]),
                "Candidate_Positions_In_tRNA": candidate_trna_position,
                "Candidate_Positions_In_Parent": _first_existing(row, ["Candidate_Positions_In_Parent", "Candidate_Position_In_Parent", "Positions_In_Parent", "Parent_Position"]),
                "Best_Final_Score": final_score,
                "Best_Final_Confidence": final_confidence,
                "Best_Biological_Context_Score": context_score,
                "Has_MS2_Precursor_Evidence": has_precursor,
                "Has_Modified_Ion_Evidence": has_modified_ion,
                "Has_Position_Discriminating_Evidence": has_discriminating,
                "Position_Ambiguity_Status": ambiguity_status,
                "Position_Resolution_Basis": resolution_basis,
                "Num_Positions_In_Ambiguity_Group": _first_existing(row, ["Num_Positions_In_Ambiguity_Group", "Group_Size", "Num_Positions"]),
                "Ambiguity_Group_ID": group_id,
                "Evidence_Summary": evidence_summary,
                "Key_Warnings": warnings,
                "Recommended_Next_Check": _recommended_next_check(priority, has_discriminating, ambiguity_status, warnings),
                "Modification_Family": row.get("Modification_Family", ""),
                "Position_Class": row.get("Position_Class", ""),
                "Position_Prior_Score": row.get("Position_Prior_Score", ""),
                "Parent_Base_Compatibility": row.get("Parent_Base_Compatibility", ""),
                "MS2_Localization_Evidence": row.get("MS2_Localization_Evidence", ""),
                "Structural_Isomer_Group_ID": row.get("Structural_Isomer_Group_ID", ""),
                "Structure_Ambiguity_Status": row.get("Structure_Ambiguity_Status", ""),
                "Alternative_Structural_Candidates": row.get("Alternative_Structural_Candidates", ""),
                "Biological_Plausibility_Score": row.get("Biological_Plausibility_Score", ""),
                "Biological_Plausibility_Level": row.get("Biological_Plausibility_Level", ""),
                "Shadow_Final_Score": row.get("Shadow_Final_Score", ""),
                "Shadow_Final_Confidence": row.get("Shadow_Final_Confidence", ""),
                "Shadow_Only": row.get("Shadow_Only", ""),
                "Has_Modified_Fragment_Ion_Evidence": row.get("Has_Modified_Fragment_Ion_Evidence", ""),
                "Modified_Fragment_Match_Count": row.get("Modified_Fragment_Match_Count", ""),
                "Unique_Modified_Fragment_Ion_Count": row.get("Unique_Modified_Fragment_Ion_Count", ""),
                "Modified_Fragment_Ion_Series": row.get("Modified_Fragment_Ion_Series", ""),
                "Supporting_Modified_Fragment_Match_IDs": row.get("Supporting_Modified_Fragment_Match_IDs", ""),
                "Best_Modified_Fragment_Error_ppm": row.get("Best_Modified_Fragment_Error_ppm", ""),
                "Maximum_Modified_Fragment_Intensity": row.get("Maximum_Modified_Fragment_Intensity", ""),
                "Physical_Observed_Peak_Keys": row.get("Physical_Observed_Peak_Keys", ""),
                "Shared_Physical_Peak_Count": row.get("Shared_Physical_Peak_Count", ""),
                "Unique_Physical_Peak_Count": row.get("Unique_Physical_Peak_Count", ""),
                "Candidate_Specific_Physical_Peak_Count": row.get("Candidate_Specific_Physical_Peak_Count", ""),
                "Isomer_Group_Shared_Peak_Count": row.get("Isomer_Group_Shared_Peak_Count", ""),
                "Has_Cross_Candidate_Peak_Sharing": row.get("Has_Cross_Candidate_Peak_Sharing", ""),
                "Cross_Candidate_Peak_Sharing_Warning": row.get("Cross_Candidate_Peak_Sharing_Warning", ""),
                "Candidate_Specific_Evidence_Peak_Count": row.get("Candidate_Specific_Evidence_Peak_Count", ""),
                "Group_Shared_Evidence_Peak_Count": row.get("Group_Shared_Evidence_Peak_Count", ""),
                "Cross_Candidate_Ambiguous_Peak_Count": row.get("Cross_Candidate_Ambiguous_Peak_Count", ""),
                "Identity_Evidence_Scope": row.get("Identity_Evidence_Scope", ""),
                "Position_Localization_Status": row.get("Position_Localization_Status", ""),
                "Group_Position_Resolution_Status": row.get("Group_Position_Resolution_Status", ""),
                "Candidate_Position_Resolution_Status": row.get("Candidate_Position_Resolution_Status", ""),
                "Position_Resolution_Ceiling_Applied": row.get("Position_Resolution_Ceiling_Applied", ""),
                "Position_Resolution_Caveat": row.get("Position_Resolution_Caveat", ""),
                "Position_Discriminating_Ion_Count": row.get("Position_Discriminating_Ion_Count", ""),
                "Structure_Resolution_Status": row.get("Structure_Resolution_Status", ""),
                "Alternative_Modification_IDs": row.get("Alternative_Modification_IDs", ""),
                "MS2_Identity_Evidence_Level": row.get("MS2_Identity_Evidence_Level", ""),
                "Shadow_MS2_Identity_Score": row.get("Shadow_MS2_Identity_Score", ""),
                "Shadow_MS2_Identity_Confidence": row.get("Shadow_MS2_Identity_Confidence", ""),
                "Shadow_MS2_Identity_Priority": row.get("Shadow_MS2_Identity_Priority", ""),
                "MS2_Identity_Evidence_Reason": row.get("MS2_Identity_Evidence_Reason", ""),
                "MS2_Identity_Warnings": row.get("MS2_Identity_Warnings", ""),
                "Total_Modified_Theoretical_Ion_Count": audit_summary.get("Total_Modified_Theoretical_Ion_Count", 0),
                "Matched_Modified_Theoretical_Ion_Count": audit_summary.get("Matched_Modified_Theoretical_Ion_Count", 0),
                "Unmatched_Modified_Theoretical_Ion_Count": audit_summary.get("Unmatched_Modified_Theoretical_Ion_Count", 0),
                "Primary_Unmatched_Reason": primary_unmatched_reason(audit_summary),
                "Outside_Scan_Range_Count": audit_summary.get("Outside_Scan_Range_Count", 0),
                "No_Peak_In_Window_Count": audit_summary.get("No_Peak_In_Window_Count", 0),
                "Nearest_Peak_Outside_Tolerance_Count": audit_summary.get("Nearest_Peak_Outside_Tolerance_Count", 0),
                "Below_Threshold_Count": audit_summary.get("Below_Threshold_Count", 0),
                "Information_Unavailable_Count": audit_summary.get("Information_Unavailable_Count", 0),
                "Best_Unmatched_Error_ppm": audit_summary.get("Best_Unmatched_Error_ppm", ""),
                "Unmatched_Ion_Audit_Warnings": audit_summary.get("Audit_Warnings", ""),
                "Ambiguous_Theoretical_Ion_Count": ambiguity_summary.get("Ambiguous_Theoretical_Ion_Count", 0),
                "Ambiguous_Peak_Cluster_Count": ambiguity_summary.get("Ambiguous_Peak_Cluster_Count", 0),
                "Maximum_Ambiguous_Cluster_Size": ambiguity_summary.get("Maximum_Cluster_Size", 0),
                "Primary_Ambiguity_Pattern": ambiguity_summary.get("Primary_Ambiguity_Pattern", "insufficient_information"),
                "Ambiguity_Severity": ambiguity_summary.get("Ambiguity_Severity", "unknown"),
                "Candidate_Specific_Ambiguous_Peak_Count": ambiguity_summary.get("Candidate_Specific_Peak_Count", 0),
                "Position_Group_Shared_Peak_Count": ambiguity_summary.get("Position_Group_Shared_Peak_Count", 0),
                "Structural_Isomer_Shared_Peak_Count": ambiguity_summary.get("Structural_Isomer_Shared_Peak_Count", 0),
                "Cross_Candidate_Shared_Peak_Count": ambiguity_summary.get("Cross_Candidate_Shared_Peak_Count", 0),
                "Ambiguous_Peak_Recommended_Followup": ambiguity_summary.get("Recommended_Followup", "insufficient_information"),
                "Ambiguous_Peak_Warnings": ambiguity_summary.get("Ambiguity_Warnings", ""),
                "Notes": "Review_Priority is for triage order only; Final_Confidence is unchanged.",
            }
        )

    priority_order = {
        "A_strong_review": 1,
        "B_medium_review": 2,
        "C_ambiguous_review": 3,
        "D_weak_review": 4,
        "E_low_information": 5,
    }
    top = pd.DataFrame(rows)
    top["_Priority_Order"] = top["Review_Priority"].map(priority_order).fillna(99)
    top = top.sort_values(
        by=["_Priority_Order", "Best_Final_Score", "Best_Biological_Context_Score"],
        ascending=[True, False, False],
        kind="mergesort",
    ).head(max_candidates)
    top.insert(0, "Review_Rank", range(1, len(top) + 1))
    return top[TOP_CANDIDATE_COLUMNS]


def _decision_text(row: pd.Series) -> tuple[str, str, str, str]:
    priority = str(row.get("Review_Priority") or "")
    has_precursor = bool(row.get("Has_MS2_Precursor_Evidence"))
    has_modified = bool(row.get("Has_Modified_Ion_Evidence"))
    has_discriminating = bool(row.get("Has_Position_Discriminating_Evidence"))
    ambiguity = str(row.get("Position_Ambiguity_Status") or "")
    context_score = _float(row.get("Best_Biological_Context_Score"))

    if priority == "C_ambiguous_review":
        label = "ambiguous_review_candidate"
        text = "This candidate has precursor and modified ion evidence, but the modification position remains ambiguous."
    elif priority == "E_low_information":
        label = "insufficient_evidence"
        text = "This candidate is supported by low-information evidence and should not be treated as confirmed."
    elif priority == "D_weak_review":
        label = "weak_candidate"
        text = "This candidate is supported mainly by limited evidence and should not be treated as localized."
    elif context_score > 0:
        label = "context_supported_review_candidate"
        text = "This candidate has review evidence with additional biological context support."
    elif _contains_any(ambiguity, ["single_candidate_position"]) and not has_discriminating:
        label = "review_candidate"
        text = "This parent/modification has a single candidate position, but no position-discriminating ion evidence is available."
    else:
        label = "review_candidate"
        text = "This candidate has evidence suitable for manual review, but it is not a confirmed modification call."

    evidence_for = "; ".join(
        item
        for item, present in [
            ("MS2 precursor evidence", has_precursor),
            ("modified ion evidence", has_modified),
            ("position-discriminating ion evidence", has_discriminating),
            ("biological context support", context_score > 0),
        ]
        if present
    )
    evidence_against = str(row.get("Key_Warnings") or "")
    ambiguity_text = str(ambiguity or "No ambiguity status reported.")
    context_text = "Biological context score present." if context_score > 0 else "No biological context support reported."
    return label, text, evidence_for or "No positive evidence flags reported.", evidence_against, ambiguity_text, context_text


def _build_decisions(top: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, row in top.iterrows():
        label, text, evidence_for, evidence_against, ambiguity_text, context_text = _decision_text(row)
        rows.append(
            {
                "Review_Rank": row.get("Review_Rank"),
                "Modification_ID": row.get("Modification_ID"),
                "Parent_Fragment_ID": row.get("Parent_Fragment_ID"),
                "Candidate_Positions_In_tRNA": row.get("Candidate_Positions_In_tRNA"),
                "Review_Priority": row.get("Review_Priority"),
                "Decision_Label": label,
                "Decision_Text": text,
                "Evidence_For": evidence_for,
                "Evidence_Against": evidence_against,
                "Ambiguity_Text": ambiguity_text,
                "Context_Text": context_text,
                "Next_Action": row.get("Recommended_Next_Check"),
            }
        )
    return pd.DataFrame(rows, columns=DECISION_COLUMNS)


def _build_checklist(top: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, row in top.iterrows():
        warnings = str(row.get("Key_Warnings") or "")
        ambiguity_status = str(row.get("Position_Ambiguity_Status") or "")
        has_c = _contains_any(row.get("Evidence_Summary"), ["c-ion", "c ion", "c_ion"])
        has_y = _contains_any(row.get("Evidence_Summary"), ["y-ion", "y ion", "y_ion"])
        rows.append(
            {
                "Review_Rank": row.get("Review_Rank"),
                "Modification_ID": row.get("Modification_ID"),
                "Parent_Fragment_ID": row.get("Parent_Fragment_ID"),
                "Candidate_Positions_In_tRNA": row.get("Candidate_Positions_In_tRNA"),
                "MS1_Fragment_Evidence": _confidence_value(row.get("Best_Final_Confidence")) > 0 or row.get("Best_Final_Score", 0) > 0,
                "Known_Modification_Candidate": bool(row.get("Modification_ID")),
                "MS2_Precursor_Evidence": bool(row.get("Has_MS2_Precursor_Evidence")),
                "Modified_Ion_Evidence": bool(row.get("Has_Modified_Ion_Evidence")),
                "Informative_Modified_Ion": bool(row.get("Has_Modified_Ion_Evidence")) and not _contains_any(warnings, ["low-information"]),
                "c_Ion_Support": has_c,
                "y_Ion_Support": has_y,
                "Both_c_y_Series": has_c and has_y,
                "Position_Discriminating_Ion": bool(row.get("Has_Position_Discriminating_Evidence")),
                "Ambiguity_Group": bool(row.get("Ambiguity_Group_ID")),
                "Ambiguous_Position": _contains_any(ambiguity_status, ["ambiguous", "unresolved"]),
                "Single_Candidate_Position": _contains_any(ambiguity_status, ["single_candidate_position", "single candidate"]),
                "Curated_Source": not _contains_any(warnings, ["uncurated", "no curated"]),
                "Candidate_Policy_Allows_Mass_Search": not _contains_any(warnings, ["policy blocks", "mass search disabled"]),
                "Biological_Context_Support": _float(row.get("Best_Biological_Context_Score")) > 0,
                "Near_Isobaric_Warning": _contains_any(warnings, ["near-isobaric", "isobaric"]),
                "Confidence_Limiting_Factors": warnings,
            }
        )
    return pd.DataFrame(rows, columns=CHECKLIST_COLUMNS)


def _dashboard(top: pd.DataFrame, ambiguity: pd.DataFrame, ranking: pd.DataFrame) -> pd.DataFrame:
    confidence_counts = {name: 0 for name in ["Very_High", "High", "Medium", "Low", "Very_Low"]}
    for value in ranking.get("Final_Confidence", pd.Series(dtype=str)).fillna("").astype(str):
        key = value.strip().replace("-", "_").replace(" ", "_")
        title_key = "_".join(part.capitalize() for part in key.split("_") if part)
        if title_key in confidence_counts:
            confidence_counts[title_key] += 1

    top_mod_ids = "; ".join(str(value) for value in top.get("Modification_ID", pd.Series(dtype=str)).dropna().astype(str).head(10) if value)
    top_priorities = "; ".join(
        f"{row.Modification_ID}:{row.Review_Priority}"
        for row in top[["Modification_ID", "Review_Priority"]].itertuples(index=False)
        if row.Modification_ID
    )
    warnings = "; ".join(str(value) for value in top.get("Key_Warnings", pd.Series(dtype=str)).dropna().astype(str).head(10) if value)
    if not ambiguity.empty and "Position_Ambiguity_Status" in ambiguity.columns:
        ambiguity_statuses = ambiguity["Position_Ambiguity_Status"].fillna("").astype(str)
    else:
        ambiguity_statuses = pd.Series(dtype=str)

    row = {
        "Total_Ranked_Candidates": len(ranking),
        **confidence_counts,
        "Total_Ambiguity_Groups": len(ambiguity),
        "Single_Candidate_Position_Groups": int((ambiguity_statuses == "single_candidate_position").sum()),
        "Resolved_By_Discriminating_Ion_Groups": int(
            ambiguity_statuses.isin(["resolved", "resolved_by_discriminating_ions"]).sum()
        ),
        "Ambiguous_Groups": int((ambiguity_statuses == "ambiguous").sum()),
        "Candidates_With_Modified_Ion_Evidence": int(top.get("Has_Modified_Ion_Evidence", pd.Series(dtype=bool)).fillna(False).sum()),
        "Candidates_With_Position_Discriminating_Evidence": int(
            top.get("Has_Position_Discriminating_Evidence", pd.Series(dtype=bool)).fillna(False).sum()
        ),
        "Candidates_With_Biological_Context_Support": int((top.get("Best_Biological_Context_Score", pd.Series(dtype=float)).fillna(0) > 0).sum()),
        "Top_Modification_IDs": top_mod_ids,
        "Top_Review_Priority_Modifications": top_priorities,
        "Key_Warnings": warnings,
        "Notes": DASHBOARD_NOTE,
    }
    return pd.DataFrame([row], columns=DASHBOARD_COLUMNS)


def build_review_dashboard_results(optional_results: dict[str, Any] | None, config: Any) -> dict[str, pd.DataFrame]:
    review_config = getattr(config, "review_dashboard", {}) or {}
    if not review_config.get("enabled", True):
        return {}

    source = optional_results or {}
    ranking = _frame(source.get("Modification_Evidence_Ranking"))
    ambiguity = _frame(source.get("Modification_Ambiguity_Groups"))
    audit_summary = _frame(source.get("MS2_Unmatched_Ion_Summary"))
    ambiguous_peak_summary = _frame(source.get("MS2_Ambiguity_Summary"))
    top = _build_top_candidates(ranking, ambiguity, review_config, audit_summary, ambiguous_peak_summary)
    return {
        "Review_Dashboard": _dashboard(top, ambiguity, ranking),
        "Top_Modification_Candidates": top,
        "Candidate_Decision_Summary": _build_decisions(top),
        "Evidence_Checklist": _build_checklist(top),
    }
