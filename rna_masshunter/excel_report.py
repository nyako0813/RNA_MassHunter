from dataclasses import asdict, is_dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
from openpyxl.utils import get_column_letter

from rna_masshunter.intact_reconstruction import (
    ASSIGNMENT_DRY_RUN_COLUMNS,
    ASSIGNMENT_DRY_RUN_SUMMARY_COLUMNS,
    ASSIGNMENT_SENSITIVITY_COLUMNS,
    ASSIGNMENT_STABILITY_COLUMNS,
    ASSIGNMENT_CANDIDATE_AUDIT_COLUMNS,
    ASSIGNMENT_AMBIGUOUS_COLUMNS,
    PREASSIGNMENT_COMPARISON_COLUMNS,
    COMPARISON_CANDIDATE_COLUMNS as INTACT_COMPARISON_CANDIDATE_COLUMNS,
    COMPETITION_GROUP_COLUMNS as INTACT_COMPETITION_GROUP_COLUMNS,
    COMPETITION_SCORE_COLUMNS as INTACT_COMPETITION_SCORE_COLUMNS,
    DIAGNOSTIC_COLUMNS as INTACT_DIAGNOSTIC_COLUMNS,
    ENGINE_COMPARISON_COLUMNS,
    GROUP_COLUMNS as INTACT_GROUP_COLUMNS,
    MISSING_CHARGE_DIAGNOSTIC_COLUMNS,
    QC_COLUMNS as INTACT_QC_COLUMNS,
    RECONSTRUCTED_MASS_SPECTRUM_COLUMNS,
    RT_ENGINE_QC_SUMMARY_COLUMNS,
    RT_ENVELOPE_DIAGNOSTIC_COLUMNS,
    TARGET_REVIEW_CANDIDATE_COLUMNS as INTACT_TARGET_REVIEW_CANDIDATE_COLUMNS,
    build_assignment_dry_run_rows,
    build_assignment_dry_run_summary_rows,
    build_assignment_sensitivity_rows,
    build_assignment_stability_rows,
    build_assignment_candidate_audit_rows,
    build_assignment_ambiguous_rows,
    build_preassignment_comparison_rows,
    build_intact_comparison_candidate_rows,
    build_intact_competition_group_rows,
    build_intact_competition_score_rows,
    build_intact_envelope_group_rows,
    build_intact_reconstruction_qc,
    build_reconstructed_mass_spectrum_rows,
    build_rt_engine_qc_summary_rows,
    build_target_review_candidate_rows,
)
from rna_masshunter.ms2_annotation import (
    MS2_FRAGMENT_EVIDENCE_COLUMNS,
    MS2_ION_MATCH_COLUMNS,
    MS2_MODIFIED_PRECURSOR_COLUMNS,
    MS2_MODIFIED_THEORETICAL_ION_COLUMNS,
    MS2_MODIFIED_ION_MATCH_COLUMNS,
    MS2_LOCALIZATION_EVIDENCE_COLUMNS,
    MS2_PARENT_CANDIDATE_COLUMNS,
    MS2_SPECTRA_COLUMNS,
    MS2_SUMMARY_COLUMNS,
    MS2_THEORETICAL_ION_COLUMNS,
    MS2_UNMATCHED_COLUMNS,
)
from rna_masshunter.evidence_ranking import AMBIGUITY_GROUP_COLUMNS, RANKING_COLUMNS, SUMMARY_COLUMNS
from rna_masshunter.biological_context import CONTEXT_PRIORITY_COLUMNS
from rna_masshunter.p1_annotation import (
    P1_ANNOTATION_COLUMNS,
    P1_SUMMARY_COLUMNS,
    P1_THEORETICAL_COLUMNS,
    P1_UNMATCHED_COLUMNS,
)


EXCEL_MAX_ROWS = 1_048_576
DATA_START_ROW = 3
EXCEL_DATA_ROW_LIMIT = EXCEL_MAX_ROWS - DATA_START_ROW


INTACT_COLUMNS = [
    "Cluster_ID",
    "Reconstructed_Mass",
    "Observed_Mass",
    "In_Neutral_Mass_Search_Range",
    "Neutral_Mass_Search_Min_Da",
    "Neutral_Mass_Search_Max_Da",
    "Neutral_Mass_Range_Status",
    "In_Target_Review_Mass_Range",
    "Target_Review_Mass_Range_Status",
    "Target_Review_Priority",
    "Envelope_QC_Eligible",
    "Intact_Review_Eligible",
    "Intact_Strict_Eligible",
    "Intact_Envelope_QC_Score",
    "Intact_Envelope_QC_Rank",
    "Strict_Eligible_Rank",
    "Review_Eligible_Rank",
    "Dominant_Intact_Envelope_Flag",
    "Supporting_Peak_IDs",
    "Supporting_Peak_Count",
    "Supporting_Scan_IDs",
    "Supporting_RT_Values",
    "Supporting_Charge_States",
    "Exact_Peak_Set_Key",
    "Exact_Duplicate_Group_ID",
    "Exact_Duplicate_Count",
    "Is_Exact_Duplicate_Representative",
    "Intact_Envelope_Group_ID",
    "Envelope_Group_Size",
    "Group_Representative",
    "Group_Ambiguity_Status",
    "Comparison_Representative",
    "Comparison_Representative_Reason",
    "Comparison_Representative_Rank",
    "Excluded_From_Comparison_Reason",
    "Target_Review_Group_Representative",
    "Target_Review_Rank",
    "Dominant_Target_Review_Eligible_Flag",
    "Reconstruction_Status",
    "Reconstruction_Confidence",
    "Reconstruction_Engine",
    "RT_Window_ID",
    "RT_Window_Start_Min",
    "RT_Window_End_Min",
    "RT_Window_Center_Min",
    "Num_MS1_Scans_In_Window",
    "Peak_Aggregation_Method",
    "Anchor_MZ",
    "Anchor_Charge",
    "Predicted_Charge_States",
    "Observed_Charge_States",
    "Missing_Charge_States",
    "Missing_Charge_Predicted_MZ",
    "Num_Predicted_Charges",
    "Num_Observed_Charges",
    "Charge_Coverage_Fraction",
    "Consecutive_Charge_Run_Length",
    "Longest_Consecutive_Charge_Run",
    "Charge_Gap_Count",
    "Charge_Continuity_Fraction",
    "Peak_Usage_Count",
    "Shared_Peak_Count",
    "Shared_Peak_Fraction",
    "Local_Window_Max_Intensity",
    "Local_Relative_Peak_Intensity_Percent",
    "Local_Envelope_Relative_Intensity_Percent",
    "Neutral_Mass_Estimator",
    "Neutral_Mass_Unweighted_Mean",
    "Neutral_Mass_Weighted_Mean",
    "Neutral_Mass_Median",
    "Envelope_Internal_Error_Max_ppm",
    "Envelope_Internal_Error_Mean_ppm",
    "Envelope_Internal_Error_Median_ppm",
    "Source_RT_Window_IDs",
    "Num_Source_RT_Windows",
    "Merged_Across_RT_Windows",
    "Extended_Lower_Charges_Evaluated",
    "Extended_Upper_Charges_Evaluated",
    "Extended_Charges_Detected",
    "Extended_Weak_Charges_Detected",
    "Extended_Charges_Not_Detected",
    "Charge_Extension_Improved_Envelope",
    "Original_Charge_States",
    "Final_Charge_States",
    "Split_Envelope_Group_ID",
    "Split_Envelope_Member_Count",
    "Split_Envelope_Merged",
    "Charge_Gaps_Before_Merge",
    "Charge_Gaps_After_Merge",
    "Max_Peak_Usage_Count",
    "Mean_Peak_Usage_Count",
    "Num_Highly_Shared_Peaks",
    "Highly_Shared_Peak_Fraction",
    "Competing_Candidate_Count",
    "Peak_Sharing_Status",
    "Competing_Envelope_Group_ID",
    "Competing_Envelope_Group_Size",
    "Shared_Peak_Competitor_Count",
    "Maximum_Shared_Peak_Fraction",
    "Mean_Shared_Peak_Fraction",
    "Competitor_Cluster_IDs",
    "Is_Noncompeting_Candidate",
    "Envelope_Evidence_Score",
    "Evidence_Score_Rank_In_Competition",
    "Evidence_Score_Components",
    "Evidence_Score_Penalties",
    "Evidence_Score_Config_Version",
    "Direct_Competitor_Count",
    "Direct_Competitor_Cluster_IDs",
    "Direct_Shared_Peak_Count_Max",
    "Direct_Shared_Peak_Fraction_Max",
    "Competition_Component_Size",
    "Dry_Run_Assignment_Status",
    "Dry_Run_Selected",
    "Dry_Run_Selection_Order",
    "Supporting_Peak_Count_Before_Assignment",
    "Independent_Supporting_Peak_Count",
    "Independent_Supporting_Peak_Fraction",
    "Supporting_Charge_Count_Before_Assignment",
    "Independent_Charge_State_Count",
    "Peaks_Already_Assigned_Count",
    "Charges_Already_Assigned_Count",
    "Excluded_By_Cluster_ID",
    "Dry_Run_Exclusion_Reason",
    "Score_Margin_To_Excluding_Candidate",
    "Close_Score_Ambiguity",
    "Assignment_Confidence",
    "Shared_Observed_Peak_Count",
    "Shared_Peak_Charge_Assignment_Count",
    "Independent_Observed_Peak_Count",
    "Pass_Min_Charge_Count",
    "Pass_Min_Consecutive_Charge_Count",
    "Pass_Charge_Continuity",
    "Pass_Internal_Error",
    "Pass_Neutral_Mass_SD",
    "Pass_Neutral_Mass_Range",
    "Pass_RT_Consistency",
    "Pass_Local_Intensity",
    "Pass_Competing_Envelope",
    "Pass_Peak_Sharing",
    "Num_Strict_Criteria_Passed",
    "Num_Review_Criteria_Passed",
    "Strict_Failure_Reasons",
    "Review_Failure_Reasons",
    "Intact_Quality_Tier",
    "Quality_Tier_Reason",
    "Quality_Tier_Rank",
    "Comparison_Ready_Strict",
    "Comparison_Ready_Review",
    "Comparison_Ready",
    "Comparison_Readiness_Reason",
    "Total_Supporting_Intensity",
    "Mean_Supporting_Intensity",
    "Max_Supporting_Intensity",
    "Reconstructed_Envelope_Intensity",
    "Intensity_Method",
    "Relative_Envelope_Intensity_Percent",
    "Relative_Overall_Envelope_Intensity_Percent",
    "Relative_In_Range_Raw_Intensity_Percent",
    "Relative_Intact_Eligible_Intensity_Percent",
    "Supporting_Peak_Classes",
    "Trace_Only_Envelope",
    "Num_Supporting_Charge_States",
    "Charge_State_Count",
    "Charge_States",
    "Charge_State_Range",
    "Charge_State_Continuity",
    "Supporting_Peak_Count",
    "RT_Min",
    "RT_Max",
    "RT_Mean",
    "RT_Range_Min",
    "Max_RT_Difference_Min",
    "RT_Consistency",
    "Neutral_Mass_SD",
    "Neutral_Mass_Range",
    "Envelope_Internal_Error_ppm",
    "Max_Mass_Error_ppm",
    "Theoretical_Mass",
    "Mass_Error_Da",
    "Mass_Error_ppm",
    "Unmodified_Theory_Delta_Da",
    "Unmodified_Theory_Delta_ppm",
    "Best_Reference_Label",
    "Best_Reference_Mass_Da",
    "Reference_Mass_Error_Da",
    "Reference_Mass_Error_ppm",
    "Reference_Mass_Matched",
    "Competing_Envelope_Count",
    "Limiting_Factors",
    "Severe_Limiting_Factors",
    "Num_Limiting_Factors",
    "Primary_Limiting_Factor",
    "Total_Intensity",
    "Assignment",
    "Confidence",
    "Warnings",
]

CHARGE_COLUMNS = ["Cluster_ID", "mz", "Intensity", "RT", "Scan_ID", "Charge", "Neutral_Mass", "Peak_Tier"]

THEORETICAL_FRAGMENT_COLUMNS = [
    "Fragment_ID",
    "Target_ID",
    "Sequence",
    "Length",
    "Start",
    "End",
    "Enzyme",
    "Missed_Cleavages",
    "Terminal_Form",
    "Unmodified_Mass",
    "Warnings",
]

FRAGMENT_MS1_MATCH_COLUMNS = [
    "Match_ID",
    "Fragment_ID",
    "Target_ID",
    "Sequence",
    "Start",
    "End",
    "Enzyme",
    "Missed_Cleavages",
    "Terminal_Form",
    "Fragment_Mass",
    "Charge",
    "Theoretical_mz",
    "Observed_mz",
    "Mass_Error_Da",
    "Mass_Error_ppm",
    "Intensity",
    "RT",
    "Scan_ID",
    "Peak_Tier",
    "Confidence",
    "Warnings",
]

FRAGMENT_MS1_FILTERED_COLUMNS = [
    "Match_ID",
    "Fragment_ID",
    "Target_ID",
    "Sequence",
    "Length",
    "Start",
    "End",
    "Enzyme",
    "Missed_Cleavages",
    "Terminal_Form",
    "Fragment_Mass",
    "Charge",
    "Theoretical_mz",
    "Observed_mz",
    "Mass_Error_Da",
    "Mass_Error_ppm",
    "Intensity",
    "RT",
    "Scan_ID",
    "Peak_Tier",
    "Confidence",
    "Warnings",
]

FRAGMENT_MS1_SUMMARY_COLUMNS = [
    "Fragment_ID",
    "Target_ID",
    "Sequence",
    "Length",
    "Start",
    "End",
    "Enzyme",
    "Missed_Cleavages",
    "Terminal_Form",
    "Best_Charge",
    "Best_Theoretical_mz",
    "Best_Observed_mz",
    "Best_Mass_Error_ppm",
    "Best_Intensity",
    "Best_RT",
    "Best_Peak_Tier",
    "Best_Confidence",
    "Match_Count",
    "Major_Count",
    "Minor_Count",
    "Trace_Count",
    "High_Count",
    "Medium_Count",
    "Low_Count",
]

KNOWN_MODIFICATION_CANDIDATE_COLUMNS = [
    "candidate_id",
    "source_type",
    "source_id",
    "target_id",
    "sequence",
    "start",
    "end",
    "observed_mz",
    "theoretical_mz",
    "observed_mass",
    "unmodified_mass",
    "mass_error_unmodified_da",
    "mass_error_unmodified_ppm",
    "modification_id",
    "modification_symbol",
    "modification_name",
    "target_base",
    "modification_mass_shift",
    "modified_mass",
    "mass_error_modified_da",
    "mass_error_modified_ppm",
    "charge",
    "intensity",
    "rt",
    "peak_tier",
    "confidence",
    "priority_score",
    "notes",
    "warnings",
]

WORKFLOW_SUMMARY_COLUMNS = [
    "Analysis_Mode",
    "Step_Name",
    "Step_Status",
    "Enabled_By_Config",
    "Executed",
    "Skip_Reason",
    "Output_Sheets",
    "Notes",
]


KNOWN_MODIFICATION_SUMMARY_COLUMNS = [
    "Modification_ID",
    "Modification_Name",
    "Symbol",
    "Target_Base",
    "Candidate_Count",
    "Best_Source_ID",
    "Best_Sequence",
    "Best_Mass_Error_Modified_ppm",
    "Best_Intensity",
    "Best_Peak_Tier",
    "Best_Confidence",
    "Best_Priority_Score",
]

SHEET_DESCRIPTIONS = {
    "Run_summary": "Run-level summary for this RNA_MassHunter MVP-3 report.",
    "Workflow_Summary": "Workflow step execution and skip status for the selected analysis mode.",
    "Input_parameters": "Flattened parameters loaded from config.yaml.",
    "mzML_diagnostics": "mzML scan counts, ranges, precursor metadata, and warnings.",
    "Intact_mass_reconstruction": "Reconstructed intact mass clusters, mass errors, and reconstruction QC fields.",
    "Charge_state_peaks": "Peak and charge-state evidence supporting reconstructed masses.",
    "Intact_Reconstruction_QC": "Per-candidate intact mass reconstruction quality diagnostics.",
    "Intact_Reconstruction_Diag": "Run-level intact reconstruction QC settings, status counts, and limiting reasons.",
    "Intact_Envelope_Groups": "Grouped intact envelope candidates and selected group representatives.",
    "Intact_Comparison_Candidates": "Group representative intact candidates suitable for condition comparison.",
    "Target_Review_Candidates": "Optional target review range candidates when configured.",
    "Reconstructed_Mass_Spectrum": "Neutral-mass reconstructed spectrum points with representative envelope intensities.",
    "RT_Envelope_Diagnostics": "RT-localized reconstruction envelope generation diagnostics.",
    "RT_Engine_QC_Summary": "RT-localized engine quality tier, failure reason, missing charge, split envelope, and engine-match summaries.",
    "Missing_Charge_Diagnostics": "Predicted missing charge-state m/z diagnostics for RT-localized envelopes.",
    "Intact_Engine_Comparison": "Optional comparison between legacy_cluster and rt_localized intact engines.",
    "Intact_Competition_Groups": "Diagnostic groups of intact candidates sharing supporting local peaks; no candidate exclusion is applied.",
    "Intact_Competition_Scores": "Envelope-internal evidence scores and rank details within competition groups.",
    "Intact_Assignment_Dry_Run": "Diagnostic-only dry-run peak assignment for competing intact candidates.",
    "Competition_Dry_Run_Summary": "Component-level summary of diagnostic dry-run assignment outcomes.",
    "Assignment_Sensitivity": "Threshold-only competitive assignment sensitivity scenarios.",
    "Assignment_Stability": "Per-candidate selection stability across assignment scenarios.",
    "Assignment_Candidate_Audit": "Optional audit-mass candidate extraction; does not affect assignment or QC.",
    "Assignment_Ambiguous_Candidates": "Ambiguous and threshold-sensitive candidates retained for assignment review.",
    "Preassignment_Comparison": "Comparison representatives before optional assignment eligibility gating.",
    "Theoretical_fragments": "Theoretical RNase digestion fragments and terminal forms.",
    "Fragment_MS1_matches": "MS1 peak matches for unmodified theoretical fragments.",
    "Fragment_MS1_filtered": "Filtered MS1 fragment matches for practical review.",
    "Fragment_MS1_summary": "Best MS1 match per fragment with match counts.",
    "Known_Modification_Candidates": "Known modification candidates explaining fragment or intact mass shifts.",
    "Known_Modification_Summary": "Grouped summary of known modification candidates.",
    "Modification_Evidence_Summary": "Run-level counts for integrated modification evidence ranking.",
    "Modification_Evidence_Ranking": "Integrated evidence scores for prioritizing modification candidates.",
    "Modification_Ambiguity_Groups": "Position ambiguity groups for shared parent-fragment modification candidates.",
    "Biological_Context_Priorities": "Biological context settings used for generic candidate prioritization.",
    "Context_Supported_Candidates": "Ranking candidates receiving a user-configured biological context boost.",
    "P1_Summary": "Summary of P1 observed peak annotation results.",
    "P1_Theoretical_Structures": "P1 monomer and short oligonucleotide theoretical structure candidates.",
    "P1_Peak_Annotations": "Observed P1 peaks matched to theoretical structure candidates, retaining unmatched peaks.",
    "P1_Unmatched_Peaks": "Observed P1 peaks outside tolerance retained for unknown/adduct/phosphate review.",
    "MS2_Summary": "Run-level summary of MS2 c/y ion annotation.",
    "MS2_Spectra": "MS2 spectrum metadata, peak counts, and annotation status.",
    "MS2_Parent_Candidates": "Precursor m/z matches between MS2 spectra and theoretical digestion fragments.",
    "MS2_Theoretical_Ions": "Theoretical c/y RNA fragment ions generated from digestion fragments.",
    "MS2_Ion_Matches": "Matched observed MS2 peaks only; unmatched peaks are reported separately.",
    "MS2_Unmatched_Peaks": "Observed MS2 peaks outside tolerance retained for review.",
    "MS2_Fragment_Evidence": "Spectrum-parent fragment evidence summary from matched MS2 ions.",
    "MS2_Peak_Annotations": "Optional all-peak MS2 annotation sheet, disabled by default.",
    "Warnings": "Warnings and errors recorded during startup, loading, and analysis.",
}


def _flatten_dict(data: dict[str, Any], prefix: str = "") -> list[dict[str, Any]]:
    rows = []
    for key, value in data.items():
        full_key = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            rows.extend(_flatten_dict(value, full_key))
        else:
            rows.append({"Parameter": full_key, "Value": value})
    return rows


def _autosize_and_freeze(writer: pd.ExcelWriter) -> None:
    for worksheet in writer.book.worksheets:
        worksheet.freeze_panes = "A2" if worksheet.title == "Index" else "A4"
        for column_cells in worksheet.columns:
            max_length = 0
            column = get_column_letter(column_cells[0].column)
            for cell in column_cells:
                value = "" if cell.value is None else str(cell.value)
                max_length = max(max_length, min(len(value), 60))
            worksheet.column_dimensions[column].width = max(10, max_length + 2)


def _sheet_link(sheet_name: str, cell: str = "A1") -> str:
    return f"#'{sheet_name}'!{cell}"


def _coerce_to_frame(value: Any) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, list):
        return pd.DataFrame(value)
    if isinstance(value, dict):
        return pd.DataFrame([value])
    return pd.DataFrame([{"Value": value}])


def _analysis_mode(config) -> str:
    workflow_mode = str((getattr(config, "analysis", {}) or {}).get("mode") or "full")
    if workflow_mode == "intact_only":
        return "intact_only"
    reconstruction_enabled = _as_bool((config.reconstruction or {}).get("enabled"), True)
    digestion_enabled = _as_bool((config.digestion or {}).get("enabled"), True)
    if reconstruction_enabled and digestion_enabled:
        return "Intact + digested fragment analysis"
    if reconstruction_enabled:
        return "Intact reconstruction only"
    if digestion_enabled:
        return "Digested fragment MS1 mapping"
    return "No active mass analysis"


def _add_index_and_backlinks(writer: pd.ExcelWriter, sheet_names: list[str]) -> None:
    workbook = writer.book
    index_sheet = workbook["Index"]
    for row_index, sheet_name in enumerate(sheet_names, start=2):
        link_cell = index_sheet.cell(row=row_index, column=1)
        link_cell.value = sheet_name
        link_cell.hyperlink = _sheet_link(sheet_name, "A1")
        link_cell.style = "Hyperlink"

    for sheet_name in sheet_names:
        worksheet = workbook[sheet_name]
        worksheet["A1"] = "← Back to Index"
        worksheet["A1"].hyperlink = _sheet_link("Index", "A1")
        worksheet["A1"].style = "Hyperlink"


def _fragment_rows(theoretical_fragments: list[Any]) -> list[dict[str, Any]]:
    rows = []
    for item in theoretical_fragments:
        raw = asdict(item) if is_dataclass(item) else dict(item)
        fragment_warnings = raw.get("warnings", [])
        if isinstance(fragment_warnings, list):
            fragment_warnings = "; ".join(map(str, fragment_warnings))
        rows.append(
            {
                "Fragment_ID": raw.get("fragment_id"),
                "Target_ID": raw.get("target_id"),
                "Sequence": raw.get("sequence"),
                "Length": len(raw.get("sequence") or ""),
                "Start": raw.get("start"),
                "End": raw.get("end"),
                "Enzyme": raw.get("enzyme"),
                "Missed_Cleavages": raw.get("missed_cleavages"),
                "Terminal_Form": raw.get("terminal_form"),
                "Unmodified_Mass": raw.get("unmodified_mass"),
                "Warnings": fragment_warnings,
            }
        )
    return rows


def _match_raw(item: Any) -> dict[str, Any]:
    return asdict(item) if is_dataclass(item) else dict(item)


def _normalize_filter_values(values: Any) -> set[str]:
    if values is None:
        return set()
    if isinstance(values, str):
        return {values.lower()}
    return {str(value).lower() for value in values}


def _fragment_ms1_match_rows(fragment_ms1_matches: list[Any], include_length: bool = False) -> list[dict[str, Any]]:
    rows = []
    for item in fragment_ms1_matches:
        raw = _match_raw(item)
        match_warnings = raw.get("warnings", [])
        if isinstance(match_warnings, list):
            match_warnings = "; ".join(map(str, match_warnings))
        row = {
            "Match_ID": raw.get("match_id"),
            "Fragment_ID": raw.get("fragment_id"),
            "Target_ID": raw.get("target_id"),
            "Sequence": raw.get("sequence"),
            "Start": raw.get("start"),
            "End": raw.get("end"),
            "Enzyme": raw.get("enzyme"),
            "Missed_Cleavages": raw.get("missed_cleavages"),
            "Terminal_Form": raw.get("terminal_form"),
            "Fragment_Mass": raw.get("fragment_mass"),
            "Charge": raw.get("charge"),
            "Theoretical_mz": raw.get("theoretical_mz"),
            "Observed_mz": raw.get("observed_mz"),
            "Mass_Error_Da": raw.get("mass_error_da"),
            "Mass_Error_ppm": raw.get("mass_error_ppm"),
            "Intensity": raw.get("intensity"),
            "RT": raw.get("rt"),
            "Scan_ID": raw.get("scan_id"),
            "Peak_Tier": raw.get("peak_tier"),
            "Confidence": raw.get("confidence"),
            "Warnings": match_warnings,
        }
        if include_length:
            row["Length"] = len(raw.get("sequence") or "")
        rows.append(row)
    return rows


def _filter_fragment_ms1_matches(fragment_ms1_matches: list[Any], mapping_config: dict[str, Any]) -> list[Any]:
    min_length = _as_positive_int(mapping_config.get("min_fragment_length_for_filtered"), 3)
    allowed_tiers = _normalize_filter_values(mapping_config.get("filtered_peak_tiers", ["Major", "Minor"]))
    allowed_confidence = _normalize_filter_values(mapping_config.get("filtered_confidence", ["High", "Medium"]))
    filtered = []
    for item in fragment_ms1_matches:
        raw = _match_raw(item)
        if len(raw.get("sequence") or "") < min_length:
            continue
        if allowed_tiers and str(raw.get("peak_tier") or "").lower() not in allowed_tiers:
            continue
        if allowed_confidence and str(raw.get("confidence") or "").lower() not in allowed_confidence:
            continue
        filtered.append(item)
    return filtered


def _confidence_rank(value: Any) -> int:
    return {"high": 3, "medium": 2, "low": 1}.get(str(value or "").lower(), 0)


def _peak_tier_rank(value: Any) -> int:
    return {"major": 3, "minor": 2, "trace": 1}.get(str(value or "").lower(), 0)


def _best_match_sort_key(item: Any) -> tuple[int, int, float, float]:
    raw = _match_raw(item)
    return (
        -_confidence_rank(raw.get("confidence")),
        -_peak_tier_rank(raw.get("peak_tier")),
        abs(float(raw.get("mass_error_ppm") or 0.0)),
        -float(raw.get("intensity") or 0.0),
    )


def _fragment_ms1_summary_rows(fragment_ms1_matches: list[Any], mapping_config: dict[str, Any]) -> list[dict[str, Any]]:
    group_key = str(mapping_config.get("summary_best_match_by", "fragment_id") or "fragment_id")
    if group_key != "fragment_id":
        group_key = "fragment_id"

    grouped: dict[str, list[Any]] = {}
    for item in fragment_ms1_matches:
        raw = _match_raw(item)
        key = str(raw.get(group_key) or "")
        if not key:
            continue
        grouped.setdefault(key, []).append(item)

    rows = []
    for fragment_id, matches in grouped.items():
        best = min(matches, key=_best_match_sort_key)
        best_raw = _match_raw(best)
        tier_counts = {"Major": 0, "Minor": 0, "Trace": 0}
        confidence_counts = {"High": 0, "Medium": 0, "Low": 0}
        for item in matches:
            raw = _match_raw(item)
            tier = str(raw.get("peak_tier") or "")
            confidence = str(raw.get("confidence") or "")
            if tier in tier_counts:
                tier_counts[tier] += 1
            if confidence in confidence_counts:
                confidence_counts[confidence] += 1
        rows.append(
            {
                "Fragment_ID": fragment_id,
                "Target_ID": best_raw.get("target_id"),
                "Sequence": best_raw.get("sequence"),
                "Length": len(best_raw.get("sequence") or ""),
                "Start": best_raw.get("start"),
                "End": best_raw.get("end"),
                "Enzyme": best_raw.get("enzyme"),
                "Missed_Cleavages": best_raw.get("missed_cleavages"),
                "Terminal_Form": best_raw.get("terminal_form"),
                "Best_Charge": best_raw.get("charge"),
                "Best_Theoretical_mz": best_raw.get("theoretical_mz"),
                "Best_Observed_mz": best_raw.get("observed_mz"),
                "Best_Mass_Error_ppm": best_raw.get("mass_error_ppm"),
                "Best_Intensity": best_raw.get("intensity"),
                "Best_RT": best_raw.get("rt"),
                "Best_Peak_Tier": best_raw.get("peak_tier"),
                "Best_Confidence": best_raw.get("confidence"),
                "Match_Count": len(matches),
                "Major_Count": tier_counts["Major"],
                "Minor_Count": tier_counts["Minor"],
                "Trace_Count": tier_counts["Trace"],
                "High_Count": confidence_counts["High"],
                "Medium_Count": confidence_counts["Medium"],
                "Low_Count": confidence_counts["Low"],
            }
        )
    return sorted(rows, key=lambda row: (row["Start"] or 0, row["End"] or 0, row["Fragment_ID"]))


def _as_bool(value: Any, default: bool = True) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _as_positive_int(value: Any, default: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed > 0 else default


def _append_excel_warning(
    warnings: list[dict[str, Any]],
    sheet_name: str,
    original_rows: int,
    written_rows: int,
) -> None:
    warnings.append(
        {
            "Timestamp": datetime.now().isoformat(timespec="seconds"),
            "Level": "WARNING",
            "Source": "excel_report",
            "Message": "Excel sheet was truncated because it exceeded max_excel_rows_per_sheet.",
            "Context": {"sheet": sheet_name, "original_rows": original_rows, "written_rows": written_rows},
        }
    )


def _truncate_frame_if_needed(
    sheet_name: str,
    frame: pd.DataFrame,
    max_rows: int,
    truncate_large_sheets: bool,
    warnings: list[dict[str, Any]],
    truncations: list[dict[str, Any]],
) -> pd.DataFrame:
    original_rows = len(frame)
    safe_limit = min(max_rows, EXCEL_DATA_ROW_LIMIT)
    if original_rows <= safe_limit:
        return frame

    if truncate_large_sheets:
        written_rows = safe_limit
    else:
        written_rows = EXCEL_DATA_ROW_LIMIT
    written_rows = min(written_rows, original_rows)
    _append_excel_warning(warnings, sheet_name, original_rows, written_rows)
    truncations.append({"sheet": sheet_name, "original_rows": original_rows, "written_rows": written_rows})
    return frame.head(written_rows).copy()


def _truncation_summary(truncations: list[dict[str, Any]]) -> str:
    if not truncations:
        return "None"
    return "; ".join(
        f"{item['sheet']}: {item['original_rows']} -> {item['written_rows']}" for item in truncations
    )


def write_excel_report(
    output_dir: str | Path,
    config,
    diagnostics: dict[str, Any],
    intact_results: list[Any],
    charge_state_peaks: list[dict[str, Any]],
    warnings: list[dict[str, Any]],
    modifications: list[Any] | None = None,
    rule_set: dict[str, Any] | None = None,
    pathways: list[dict[str, Any]] | None = None,
    theoretical_fragments: list[Any] | None = None,
    fragment_ms1_matches: list[Any] | None = None,
    known_modification_candidates: list[dict[str, Any]] | None = None,
    known_modification_summary: list[dict[str, Any]] | None = None,
    optional_results: dict[str, Any] | None = None,
) -> Path:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / f"RNA_MassHunter_MVP5_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"

    reporting = config.reporting or {}
    max_excel_rows = _as_positive_int(reporting.get("max_excel_rows_per_sheet"), 100000)
    truncate_large_sheets = _as_bool(reporting.get("truncate_large_sheets"), True)
    max_charge_state_peak_rows = _as_positive_int(
        reporting.get("max_charge_state_peak_rows", config.reconstruction.get("max_charge_state_peak_rows")),
        max_excel_rows,
    )
    truncations: list[dict[str, Any]] = []
    reconstruction_enabled = _as_bool(config.reconstruction.get("enabled"), True)
    intact_qc_rows, intact_diagnostic_rows = build_intact_reconstruction_qc(
        intact_results,
        charge_state_peaks,
        config.reconstruction or {},
        reconstruction_enabled=reconstruction_enabled,
    )
    intact_group_rows = build_intact_envelope_group_rows(intact_qc_rows)
    intact_competition_group_rows = build_intact_competition_group_rows(intact_qc_rows)
    intact_competition_score_rows = build_intact_competition_score_rows(intact_qc_rows)
    assignment_dry_run_rows = build_assignment_dry_run_rows(intact_qc_rows)
    assignment_dry_run_summary_rows = build_assignment_dry_run_summary_rows(intact_qc_rows)
    assignment_sensitivity_rows = build_assignment_sensitivity_rows(config.reconstruction or {})
    assignment_stability_rows = build_assignment_stability_rows(intact_qc_rows)
    assignment_candidate_audit_rows = build_assignment_candidate_audit_rows(config.reconstruction or {})
    assignment_ambiguous_rows = build_assignment_ambiguous_rows(intact_qc_rows)
    preassignment_comparison_rows = build_preassignment_comparison_rows(intact_qc_rows)
    intact_comparison_rows = build_intact_comparison_candidate_rows(intact_qc_rows)
    intact_target_review_rows = build_target_review_candidate_rows(intact_qc_rows)
    reconstructed_spectrum_rows = build_reconstructed_mass_spectrum_rows(intact_qc_rows, config.reconstruction or {})
    rt_engine_qc_summary_rows = build_rt_engine_qc_summary_rows(intact_diagnostic_rows)

    intact_rows = []
    for item in intact_results:
        raw = asdict(item) if is_dataclass(item) else dict(item)
        intact_rows.append(
            {
                "Cluster_ID": raw.get("cluster_id"),
                "Reconstructed_Mass": raw.get("observed_mass"),
                "Observed_Mass": raw.get("observed_mass"),
                "In_Neutral_Mass_Search_Range": raw.get("in_neutral_mass_search_range"),
                "Neutral_Mass_Search_Min_Da": raw.get("neutral_mass_search_min_da"),
                "Neutral_Mass_Search_Max_Da": raw.get("neutral_mass_search_max_da"),
                "Neutral_Mass_Range_Status": raw.get("neutral_mass_range_status"),
                "In_Target_Review_Mass_Range": raw.get("in_target_review_mass_range"),
                "Target_Review_Mass_Range_Status": raw.get("target_review_mass_range_status"),
                "Target_Review_Priority": raw.get("target_review_priority"),
                "Envelope_QC_Eligible": raw.get("envelope_qc_eligible"),
                "Intact_Review_Eligible": raw.get("intact_review_eligible"),
                "Intact_Strict_Eligible": raw.get("intact_strict_eligible"),
                "Intact_Envelope_QC_Score": raw.get("intact_envelope_qc_score"),
                "Intact_Envelope_QC_Rank": raw.get("intact_envelope_qc_rank"),
                "Strict_Eligible_Rank": raw.get("strict_eligible_rank"),
                "Review_Eligible_Rank": raw.get("review_eligible_rank"),
                "Dominant_Intact_Envelope_Flag": raw.get("dominant_intact_envelope_flag"),
                "Supporting_Peak_IDs": raw.get("supporting_peak_ids"),
                "Supporting_Peak_Count": raw.get("supporting_peak_count"),
                "Supporting_Scan_IDs": raw.get("supporting_scan_ids"),
                "Supporting_RT_Values": raw.get("supporting_rt_values"),
                "Supporting_Charge_States": raw.get("supporting_charge_states"),
                "Exact_Peak_Set_Key": raw.get("exact_peak_set_key"),
                "Exact_Duplicate_Group_ID": raw.get("exact_duplicate_group_id"),
                "Exact_Duplicate_Count": raw.get("exact_duplicate_count"),
                "Is_Exact_Duplicate_Representative": raw.get("is_exact_duplicate_representative"),
                "Intact_Envelope_Group_ID": raw.get("intact_envelope_group_id"),
                "Envelope_Group_Size": raw.get("envelope_group_size"),
                "Group_Representative": raw.get("group_representative"),
                "Group_Ambiguity_Status": raw.get("group_ambiguity_status"),
                "Comparison_Representative": raw.get("comparison_representative"),
                "Comparison_Representative_Reason": raw.get("comparison_representative_reason"),
                "Comparison_Representative_Rank": raw.get("comparison_representative_rank"),
                "Excluded_From_Comparison_Reason": raw.get("excluded_from_comparison_reason"),
                "Target_Review_Group_Representative": raw.get("target_review_group_representative"),
                "Target_Review_Rank": raw.get("target_review_rank"),
                "Dominant_Target_Review_Eligible_Flag": raw.get("dominant_target_review_eligible_flag"),
                "Reconstruction_Status": raw.get("reconstruction_status"),
                "Reconstruction_Confidence": raw.get("reconstruction_confidence"),
                "Reconstruction_Engine": raw.get("reconstruction_engine"),
                "RT_Window_ID": raw.get("rt_window_id"),
                "RT_Window_Start_Min": raw.get("rt_window_start_min"),
                "RT_Window_End_Min": raw.get("rt_window_end_min"),
                "RT_Window_Center_Min": raw.get("rt_window_center_min"),
                "Num_MS1_Scans_In_Window": raw.get("num_ms1_scans_in_window"),
                "Peak_Aggregation_Method": raw.get("peak_aggregation_method"),
                "Anchor_MZ": raw.get("anchor_mz"),
                "Anchor_Charge": raw.get("anchor_charge"),
                "Predicted_Charge_States": raw.get("predicted_charge_states"),
                "Observed_Charge_States": raw.get("observed_charge_states"),
                "Missing_Charge_States": raw.get("missing_charge_states"),
                "Missing_Charge_Predicted_MZ": raw.get("missing_charge_predicted_mz"),
                "Num_Predicted_Charges": raw.get("num_predicted_charges"),
                "Num_Observed_Charges": raw.get("num_observed_charges"),
                "Charge_Coverage_Fraction": raw.get("charge_coverage_fraction"),
                "Consecutive_Charge_Run_Length": raw.get("consecutive_charge_run_length"),
                "Longest_Consecutive_Charge_Run": raw.get("longest_consecutive_charge_run"),
                "Charge_Gap_Count": raw.get("charge_gap_count"),
                "Charge_Continuity_Fraction": raw.get("charge_continuity_fraction"),
                "Peak_Usage_Count": raw.get("peak_usage_count"),
                "Shared_Peak_Count": raw.get("shared_peak_count"),
                "Shared_Peak_Fraction": raw.get("shared_peak_fraction"),
                "Local_Window_Max_Intensity": raw.get("local_window_max_intensity"),
                "Local_Relative_Peak_Intensity_Percent": raw.get("local_relative_peak_intensity_percent"),
                "Local_Envelope_Relative_Intensity_Percent": raw.get("local_envelope_relative_intensity_percent"),
                "Neutral_Mass_Estimator": raw.get("neutral_mass_estimator"),
                "Neutral_Mass_Unweighted_Mean": raw.get("neutral_mass_unweighted_mean"),
                "Neutral_Mass_Weighted_Mean": raw.get("neutral_mass_weighted_mean"),
                "Neutral_Mass_Median": raw.get("neutral_mass_median"),
                "Envelope_Internal_Error_Max_ppm": raw.get("envelope_internal_error_max_ppm"),
                "Envelope_Internal_Error_Mean_ppm": raw.get("envelope_internal_error_mean_ppm"),
                "Envelope_Internal_Error_Median_ppm": raw.get("envelope_internal_error_median_ppm"),
                "Source_RT_Window_IDs": raw.get("source_rt_window_ids"),
                "Num_Source_RT_Windows": raw.get("num_source_rt_windows"),
                "Merged_Across_RT_Windows": raw.get("merged_across_rt_windows"),
                "Extended_Lower_Charges_Evaluated": raw.get("extended_lower_charges_evaluated"),
                "Extended_Upper_Charges_Evaluated": raw.get("extended_upper_charges_evaluated"),
                "Extended_Charges_Detected": raw.get("extended_charges_detected"),
                "Extended_Weak_Charges_Detected": raw.get("extended_weak_charges_detected"),
                "Extended_Charges_Not_Detected": raw.get("extended_charges_not_detected"),
                "Charge_Extension_Improved_Envelope": raw.get("charge_extension_improved_envelope"),
                "Original_Charge_States": raw.get("original_charge_states"),
                "Final_Charge_States": raw.get("final_charge_states"),
                "Split_Envelope_Group_ID": raw.get("split_envelope_group_id"),
                "Split_Envelope_Member_Count": raw.get("split_envelope_member_count"),
                "Split_Envelope_Merged": raw.get("split_envelope_merged"),
                "Charge_Gaps_Before_Merge": raw.get("charge_gaps_before_merge"),
                "Charge_Gaps_After_Merge": raw.get("charge_gaps_after_merge"),
                "Max_Peak_Usage_Count": raw.get("max_peak_usage_count"),
                "Mean_Peak_Usage_Count": raw.get("mean_peak_usage_count"),
                "Num_Highly_Shared_Peaks": raw.get("num_highly_shared_peaks"),
                "Highly_Shared_Peak_Fraction": raw.get("highly_shared_peak_fraction"),
                "Competing_Candidate_Count": raw.get("competing_candidate_count"),
                "Peak_Sharing_Status": raw.get("peak_sharing_status"),
                "Competing_Envelope_Group_ID": raw.get("competing_envelope_group_id"),
                "Competing_Envelope_Group_Size": raw.get("competing_envelope_group_size"),
                "Shared_Peak_Competitor_Count": raw.get("shared_peak_competitor_count"),
                "Maximum_Shared_Peak_Fraction": raw.get("maximum_shared_peak_fraction"),
                "Mean_Shared_Peak_Fraction": raw.get("mean_shared_peak_fraction"),
                "Competitor_Cluster_IDs": raw.get("competitor_cluster_ids"),
                "Is_Noncompeting_Candidate": raw.get("is_noncompeting_candidate"),
                "Envelope_Evidence_Score": raw.get("envelope_evidence_score"),
                "Evidence_Score_Rank_In_Competition": raw.get("evidence_score_rank_in_competition"),
                "Evidence_Score_Components": raw.get("evidence_score_components"),
                "Evidence_Score_Penalties": raw.get("evidence_score_penalties"),
                "Evidence_Score_Config_Version": raw.get("evidence_score_config_version"),
                "Direct_Competitor_Count": raw.get("direct_competitor_count"),
                "Direct_Competitor_Cluster_IDs": raw.get("direct_competitor_cluster_ids"),
                "Direct_Shared_Peak_Count_Max": raw.get("direct_shared_peak_count_max"),
                "Direct_Shared_Peak_Fraction_Max": raw.get("direct_shared_peak_fraction_max"),
                "Competition_Component_Size": raw.get("competition_component_size"),
                "Dry_Run_Assignment_Status": raw.get("dry_run_assignment_status"),
                "Dry_Run_Selected": raw.get("dry_run_selected"),
                "Dry_Run_Selection_Order": raw.get("dry_run_selection_order"),
                "Supporting_Peak_Count_Before_Assignment": raw.get("supporting_peak_count_before_assignment"),
                "Independent_Supporting_Peak_Count": raw.get("independent_supporting_peak_count"),
                "Independent_Supporting_Peak_Fraction": raw.get("independent_supporting_peak_fraction"),
                "Supporting_Charge_Count_Before_Assignment": raw.get("supporting_charge_count_before_assignment"),
                "Independent_Charge_State_Count": raw.get("independent_charge_state_count"),
                "Peaks_Already_Assigned_Count": raw.get("peaks_already_assigned_count"),
                "Charges_Already_Assigned_Count": raw.get("charges_already_assigned_count"),
                "Excluded_By_Cluster_ID": raw.get("excluded_by_cluster_id"),
                "Dry_Run_Exclusion_Reason": raw.get("dry_run_exclusion_reason"),
                "Score_Margin_To_Excluding_Candidate": raw.get("score_margin_to_excluding_candidate"),
                "Close_Score_Ambiguity": raw.get("close_score_ambiguity"),
                "Assignment_Confidence": raw.get("assignment_confidence"),
                "Shared_Observed_Peak_Count": raw.get("shared_observed_peak_count"),
                "Shared_Peak_Charge_Assignment_Count": raw.get("shared_peak_charge_assignment_count"),
                "Independent_Observed_Peak_Count": raw.get("independent_observed_peak_count"),
                "Pass_Min_Charge_Count": raw.get("pass_min_charge_count"),
                "Pass_Min_Consecutive_Charge_Count": raw.get("pass_min_consecutive_charge_count"),
                "Pass_Charge_Continuity": raw.get("pass_charge_continuity"),
                "Pass_Internal_Error": raw.get("pass_internal_error"),
                "Pass_Neutral_Mass_SD": raw.get("pass_neutral_mass_sd"),
                "Pass_Neutral_Mass_Range": raw.get("pass_neutral_mass_range"),
                "Pass_RT_Consistency": raw.get("pass_rt_consistency"),
                "Pass_Local_Intensity": raw.get("pass_local_intensity"),
                "Pass_Competing_Envelope": raw.get("pass_competing_envelope"),
                "Pass_Peak_Sharing": raw.get("pass_peak_sharing"),
                "Num_Strict_Criteria_Passed": raw.get("num_strict_criteria_passed"),
                "Num_Review_Criteria_Passed": raw.get("num_review_criteria_passed"),
                "Strict_Failure_Reasons": raw.get("strict_failure_reasons"),
                "Review_Failure_Reasons": raw.get("review_failure_reasons"),
                "Intact_Quality_Tier": raw.get("intact_quality_tier"),
                "Quality_Tier_Reason": raw.get("quality_tier_reason"),
                "Quality_Tier_Rank": raw.get("quality_tier_rank"),
                "Comparison_Ready_Strict": raw.get("comparison_ready_strict"),
                "Comparison_Ready_Review": raw.get("comparison_ready_review"),
                "Comparison_Ready": raw.get("comparison_ready"),
                "Comparison_Readiness_Reason": raw.get("comparison_readiness_reason"),
                "Total_Supporting_Intensity": raw.get("total_supporting_intensity"),
                "Mean_Supporting_Intensity": raw.get("mean_supporting_intensity"),
                "Max_Supporting_Intensity": raw.get("max_supporting_intensity"),
                "Reconstructed_Envelope_Intensity": raw.get("reconstructed_envelope_intensity"),
                "Intensity_Method": raw.get("intensity_method"),
                "Relative_Envelope_Intensity_Percent": raw.get("relative_envelope_intensity_percent"),
                "Relative_Overall_Envelope_Intensity_Percent": raw.get("relative_overall_envelope_intensity_percent"),
                "Relative_In_Range_Raw_Intensity_Percent": raw.get("relative_in_range_raw_intensity_percent"),
                "Relative_Intact_Eligible_Intensity_Percent": raw.get("relative_intact_eligible_intensity_percent"),
                "Supporting_Peak_Classes": raw.get("supporting_peak_classes"),
                "Trace_Only_Envelope": raw.get("trace_only_envelope"),
                "Num_Supporting_Charge_States": raw.get("num_supporting_charge_states"),
                "Charge_State_Count": raw.get("charge_state_count"),
                "Charge_States": ",".join(map(str, raw.get("charge_states", []))),
                "Charge_State_Range": raw.get("charge_state_range"),
                "Charge_State_Continuity": raw.get("charge_state_continuity"),
                "Supporting_Peak_Count": raw.get("supporting_peak_count"),
                "RT_Min": raw.get("rt_min"),
                "RT_Max": raw.get("rt_max"),
                "RT_Mean": raw.get("rt_mean"),
                "RT_Range_Min": raw.get("rt_range_min"),
                "Max_RT_Difference_Min": raw.get("max_rt_difference_min"),
                "RT_Consistency": raw.get("rt_consistency"),
                "Neutral_Mass_SD": raw.get("neutral_mass_sd"),
                "Neutral_Mass_Range": raw.get("neutral_mass_range"),
                "Envelope_Internal_Error_ppm": raw.get("envelope_internal_error_ppm"),
                "Max_Mass_Error_ppm": raw.get("max_mass_error_ppm"),
                "Theoretical_Mass": raw.get("theoretical_mass"),
                "Mass_Error_Da": raw.get("mass_error_da"),
                "Mass_Error_ppm": raw.get("mass_error_ppm"),
                "Unmodified_Theory_Delta_Da": raw.get("unmodified_theory_delta_da"),
                "Unmodified_Theory_Delta_ppm": raw.get("unmodified_theory_delta_ppm"),
                "Best_Reference_Label": raw.get("best_reference_label"),
                "Best_Reference_Mass_Da": raw.get("best_reference_mass_da"),
                "Reference_Mass_Error_Da": raw.get("reference_mass_error_da"),
                "Reference_Mass_Error_ppm": raw.get("reference_mass_error_ppm"),
                "Reference_Mass_Matched": raw.get("reference_mass_matched"),
                "Competing_Envelope_Count": raw.get("competing_envelope_count"),
                "Limiting_Factors": raw.get("limiting_factors"),
                "Severe_Limiting_Factors": raw.get("severe_limiting_factors"),
                "Num_Limiting_Factors": raw.get("num_limiting_factors"),
                "Primary_Limiting_Factor": raw.get("primary_limiting_factor"),
                "Total_Intensity": raw.get("total_intensity"),
                "Assignment": raw.get("assignment"),
                "Confidence": raw.get("confidence"),
                "Warnings": raw.get("warnings"),
            }
        )

    charge_state_peak_rows = charge_state_peaks
    if len(charge_state_peaks) > max_charge_state_peak_rows and truncate_large_sheets:
        _append_excel_warning(warnings, "Charge_state_peaks", len(charge_state_peaks), max_charge_state_peak_rows)
        truncations.append(
            {
                "sheet": "Charge_state_peaks",
                "original_rows": len(charge_state_peaks),
                "written_rows": max_charge_state_peak_rows,
            }
        )
        charge_state_peak_rows = charge_state_peaks[:max_charge_state_peak_rows]

    theoretical_fragments = theoretical_fragments or []
    fragment_ms1_matches = fragment_ms1_matches or []
    known_modification_candidates = known_modification_candidates or []
    known_modification_summary = known_modification_summary or []
    fragment_ms1_filtered = _filter_fragment_ms1_matches(fragment_ms1_matches, config.fragment_mapping or {})
    fragment_ms1_summary_rows = _fragment_ms1_summary_rows(fragment_ms1_matches, config.fragment_mapping or {})
    input_parameters = {
        "analysis": getattr(config, "analysis", {}),
        "project": config.project,
        "input": config.input,
        "organism": config.organism,
        "sequence": config.sequence,
        "experiment": config.experiment,
        "instrument": config.instrument,
        "reconstruction": config.reconstruction,
        "digestion": config.digestion,
        "alkaline_phosphatase": config.alkaline_phosphatase,
        "fragment_mapping": config.fragment_mapping,
        "modification_search": config.modification_search,
        "peak_filtering": config.peak_filtering,
        "p1_annotation": config.p1_annotation,
        "ms2_annotation": config.ms2_annotation,
        "modification_evidence_ranking": config.modification_evidence_ranking,
        "biological_context": config.biological_context,
        "performance": config.performance,
        "reporting": config.reporting,
    }
    analysis_mode = str((getattr(config, "analysis", {}) or {}).get("mode") or "full")
    optional_results = optional_results or {}
    workflow_summary_rows = optional_results.get("Workflow_Summary") or [
        {
            "Analysis_Mode": analysis_mode,
            "Step_Name": "workflow_summary",
            "Step_Status": "executed",
            "Enabled_By_Config": True,
            "Executed": True,
            "Skip_Reason": "",
            "Output_Sheets": "Workflow_Summary",
            "Notes": "Default summary row generated by report writer.",
        }
    ]

    data_sheets: dict[str, pd.DataFrame] = {
        "Workflow_Summary": pd.DataFrame(workflow_summary_rows, columns=WORKFLOW_SUMMARY_COLUMNS),
        "Input_parameters": pd.DataFrame(_flatten_dict(input_parameters)),
        "mzML_diagnostics": pd.DataFrame([diagnostics] if diagnostics else [{}]),
        "Theoretical_fragments": pd.DataFrame(_fragment_rows(theoretical_fragments), columns=THEORETICAL_FRAGMENT_COLUMNS),
        "Fragment_MS1_matches": pd.DataFrame(_fragment_ms1_match_rows(fragment_ms1_matches), columns=FRAGMENT_MS1_MATCH_COLUMNS),
        "Fragment_MS1_filtered": pd.DataFrame(_fragment_ms1_match_rows(fragment_ms1_filtered, include_length=True), columns=FRAGMENT_MS1_FILTERED_COLUMNS),
        "Fragment_MS1_summary": pd.DataFrame(fragment_ms1_summary_rows, columns=FRAGMENT_MS1_SUMMARY_COLUMNS),
        "Known_Modification_Candidates": pd.DataFrame(known_modification_candidates, columns=KNOWN_MODIFICATION_CANDIDATE_COLUMNS),
        "Known_Modification_Summary": pd.DataFrame(known_modification_summary, columns=KNOWN_MODIFICATION_SUMMARY_COLUMNS),
    }
    intact_qc_sheets = {
        "Intact_Reconstruction_QC": pd.DataFrame(intact_qc_rows, columns=INTACT_QC_COLUMNS),
        "Intact_Reconstruction_Diag": pd.DataFrame(intact_diagnostic_rows, columns=INTACT_DIAGNOSTIC_COLUMNS),
        "Intact_Envelope_Groups": pd.DataFrame(intact_group_rows, columns=INTACT_GROUP_COLUMNS),
        "Intact_Competition_Groups": pd.DataFrame(intact_competition_group_rows, columns=INTACT_COMPETITION_GROUP_COLUMNS),
        "Intact_Competition_Scores": pd.DataFrame(intact_competition_score_rows, columns=INTACT_COMPETITION_SCORE_COLUMNS),
        "Intact_Assignment_Dry_Run": pd.DataFrame(assignment_dry_run_rows, columns=ASSIGNMENT_DRY_RUN_COLUMNS),
        "Competition_Dry_Run_Summary": pd.DataFrame(assignment_dry_run_summary_rows, columns=ASSIGNMENT_DRY_RUN_SUMMARY_COLUMNS),
        "Assignment_Sensitivity": pd.DataFrame(assignment_sensitivity_rows, columns=ASSIGNMENT_SENSITIVITY_COLUMNS),
        "Assignment_Stability": pd.DataFrame(assignment_stability_rows, columns=ASSIGNMENT_STABILITY_COLUMNS),
        "Assignment_Candidate_Audit": pd.DataFrame(assignment_candidate_audit_rows, columns=ASSIGNMENT_CANDIDATE_AUDIT_COLUMNS),
        "Assignment_Ambiguous_Candidates": pd.DataFrame(assignment_ambiguous_rows, columns=ASSIGNMENT_AMBIGUOUS_COLUMNS),
        "Preassignment_Comparison": pd.DataFrame(preassignment_comparison_rows, columns=PREASSIGNMENT_COMPARISON_COLUMNS),
        "Intact_Comparison_Candidates": pd.DataFrame(intact_comparison_rows, columns=INTACT_COMPARISON_CANDIDATE_COLUMNS),
        "Target_Review_Candidates": pd.DataFrame(intact_target_review_rows, columns=INTACT_TARGET_REVIEW_CANDIDATE_COLUMNS),
        "Reconstructed_Mass_Spectrum": pd.DataFrame(reconstructed_spectrum_rows, columns=RECONSTRUCTED_MASS_SPECTRUM_COLUMNS),
        "RT_Engine_QC_Summary": pd.DataFrame(rt_engine_qc_summary_rows, columns=RT_ENGINE_QC_SUMMARY_COLUMNS),
    }
    if reconstruction_enabled:
        data_sheets = {
            "Workflow_Summary": data_sheets["Workflow_Summary"],
            "Input_parameters": data_sheets["Input_parameters"],
            "mzML_diagnostics": data_sheets["mzML_diagnostics"],
            "Intact_mass_reconstruction": pd.DataFrame(intact_rows, columns=INTACT_COLUMNS),
            "Charge_state_peaks": pd.DataFrame(charge_state_peak_rows, columns=CHARGE_COLUMNS),
            **intact_qc_sheets,
            **{key: value for key, value in data_sheets.items() if key not in {"Input_parameters", "mzML_diagnostics"}},
        }
    else:
        data_sheets = {
            "Workflow_Summary": data_sheets["Workflow_Summary"],
            "Input_parameters": data_sheets["Input_parameters"],
            "mzML_diagnostics": data_sheets["mzML_diagnostics"],
            **intact_qc_sheets,
            **{key: value for key, value in data_sheets.items() if key not in {"Input_parameters", "mzML_diagnostics"}},
        }
    if analysis_mode == "intact_only":
        intact_only_sheet_names = {
            "Workflow_Summary",
            "Input_parameters",
            "mzML_diagnostics",
            "Intact_mass_reconstruction",
            "Charge_state_peaks",
            "Intact_Reconstruction_QC",
            "Intact_Reconstruction_Diag",
            "Intact_Envelope_Groups",
            "Intact_Competition_Groups",
            "Intact_Competition_Scores",
            "Intact_Assignment_Dry_Run",
            "Competition_Dry_Run_Summary",
            "Assignment_Sensitivity",
            "Assignment_Stability",
            "Assignment_Candidate_Audit",
            "Assignment_Ambiguous_Candidates",
            "Preassignment_Comparison",
            "Intact_Comparison_Candidates",
            "Target_Review_Candidates",
            "Reconstructed_Mass_Spectrum",
            "RT_Engine_QC_Summary",
            "RT_Envelope_Diagnostics",
            "Missing_Charge_Diagnostics",
            "Intact_Engine_Comparison",
        }
        data_sheets = {key: value for key, value in data_sheets.items() if key in intact_only_sheet_names}

    optional_columns = {
        "P1_Summary": P1_SUMMARY_COLUMNS,
        "P1_Theoretical_Structures": P1_THEORETICAL_COLUMNS,
        "P1_Peak_Annotations": P1_ANNOTATION_COLUMNS,
        "P1_Unmatched_Peaks": P1_UNMATCHED_COLUMNS,
        "MS2_Summary": MS2_SUMMARY_COLUMNS,
        "MS2_Spectra": MS2_SPECTRA_COLUMNS,
        "MS2_Parent_Candidates": MS2_PARENT_CANDIDATE_COLUMNS,
        "MS2_Modified_Precursor_Candidates": MS2_MODIFIED_PRECURSOR_COLUMNS,
        "MS2_Modified_Theoretical_Ions": MS2_MODIFIED_THEORETICAL_ION_COLUMNS,
        "MS2_Modified_Ion_Matches": MS2_MODIFIED_ION_MATCH_COLUMNS,
        "MS2_Modification_Localization_Evidence": MS2_LOCALIZATION_EVIDENCE_COLUMNS,
        "Modification_Evidence_Summary": SUMMARY_COLUMNS,
        "Modification_Evidence_Ranking": RANKING_COLUMNS,
        "Modification_Ambiguity_Groups": AMBIGUITY_GROUP_COLUMNS,
        "Biological_Context_Priorities": CONTEXT_PRIORITY_COLUMNS,
        "Context_Supported_Candidates": RANKING_COLUMNS,
        "MS2_Theoretical_Ions": MS2_THEORETICAL_ION_COLUMNS,
        "MS2_Ion_Matches": MS2_ION_MATCH_COLUMNS,
        "MS2_Unmatched_Peaks": MS2_UNMATCHED_COLUMNS,
        "MS2_Fragment_Evidence": MS2_FRAGMENT_EVIDENCE_COLUMNS,
        "MS2_Peak_Annotations": MS2_ION_MATCH_COLUMNS,
        "RT_Envelope_Diagnostics": RT_ENVELOPE_DIAGNOSTIC_COLUMNS,
        "RT_Engine_QC_Summary": RT_ENGINE_QC_SUMMARY_COLUMNS,
        "Missing_Charge_Diagnostics": MISSING_CHARGE_DIAGNOSTIC_COLUMNS,
        "Intact_Engine_Comparison": ENGINE_COMPARISON_COLUMNS,
        "Intact_Competition_Groups": INTACT_COMPETITION_GROUP_COLUMNS,
        "Intact_Competition_Scores": INTACT_COMPETITION_SCORE_COLUMNS,
        "Intact_Assignment_Dry_Run": ASSIGNMENT_DRY_RUN_COLUMNS,
        "Competition_Dry_Run_Summary": ASSIGNMENT_DRY_RUN_SUMMARY_COLUMNS,
        "Assignment_Sensitivity": ASSIGNMENT_SENSITIVITY_COLUMNS,
        "Assignment_Stability": ASSIGNMENT_STABILITY_COLUMNS,
        "Assignment_Candidate_Audit": ASSIGNMENT_CANDIDATE_AUDIT_COLUMNS,
        "Assignment_Ambiguous_Candidates": ASSIGNMENT_AMBIGUOUS_COLUMNS,
        "Preassignment_Comparison": PREASSIGNMENT_COMPARISON_COLUMNS,
    }
    for sheet_name, value in optional_results.items():
        if sheet_name in {"Index", "Run_summary", "Warnings", "Workflow_Summary"}:
            continue
        frame = _coerce_to_frame(value)
        columns = optional_columns.get(sheet_name)
        if columns:
            frame = pd.DataFrame(frame, columns=columns)
        data_sheets[sheet_name[:31]] = frame

    truncated_data_sheets = {
        sheet_name: _truncate_frame_if_needed(
            sheet_name,
            frame,
            max_excel_rows,
            truncate_large_sheets,
            warnings,
            truncations,
        )
        for sheet_name, frame in data_sheets.items()
    }

    summary_rows = [
        {"Item": "Project", "Value": config.project.get("name")},
        {"Item": "Generated", "Value": datetime.now().isoformat(timespec="seconds")},
        {"Item": "Analysis mode", "Value": _analysis_mode(config)},
        {"Item": "Modification dictionary entries", "Value": len(modifications or [])},
        {"Item": "Rule set", "Value": config.organism.get("rule_set") or (rule_set or {}).get("id") or (rule_set or {}).get("name")},
        {"Item": "Pathway files", "Value": len(pathways or [])},
        {"Item": "Intact mass candidates", "Value": len(intact_results)},
        {"Item": "Theoretical fragments", "Value": len(theoretical_fragments)},
        {"Item": "Fragment MS1 matches", "Value": len(fragment_ms1_matches)},
        {"Item": "Fragment MS1 filtered", "Value": len(fragment_ms1_filtered)},
        {"Item": "Fragment MS1 summary", "Value": len(fragment_ms1_summary_rows)},
        {"Item": "Known modification candidates", "Value": len(known_modification_candidates)},
        {"Item": "Known modification summary", "Value": len(known_modification_summary)},
        {"Item": "Truncated sheets", "Value": _truncation_summary(truncations)},
        {"Item": "Warnings", "Value": len(warnings)},
    ]

    sheets: dict[str, pd.DataFrame] = {
        "Run_summary": pd.DataFrame(summary_rows),
        **truncated_data_sheets,
        "Warnings": pd.DataFrame(warnings, columns=["Timestamp", "Level", "Source", "Message", "Context"]),
    }
    sheets = {
        sheet_name: _truncate_frame_if_needed(
            sheet_name,
            frame,
            max_excel_rows,
            truncate_large_sheets,
            warnings,
            truncations,
        )
        if sheet_name == "Warnings"
        else frame
        for sheet_name, frame in sheets.items()
    }

    index_rows = [
        {
            "Sheet": sheet_name,
            "Description": SHEET_DESCRIPTIONS.get(sheet_name, "Optional result sheet."),
            "Notes": "Data starts at A3.",
        }
        for sheet_name in sheets
    ]

    with pd.ExcelWriter(report_path, engine="openpyxl") as writer:
        pd.DataFrame(index_rows, columns=["Sheet", "Description", "Notes"]).to_excel(writer, sheet_name="Index", index=False)
        for sheet_name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=sheet_name, index=False, startrow=2)
        _add_index_and_backlinks(writer, list(sheets))
        _autosize_and_freeze(writer)
    return report_path
