from bisect import bisect_left
from statistics import median, pstdev
from time import perf_counter
from typing import Any

from rna_masshunter.masses import mz_from_neutral_mass, neutral_mass_from_mz
from rna_masshunter.models import IntactMassCandidate, PeakTierResult


def _confidence(charge_state_count: int, min_charge_states: int) -> str:
    if charge_state_count >= min_charge_states + 2:
        return "High"
    if charge_state_count >= min_charge_states:
        return "Medium"
    return "Low"


DEFAULT_INTACT_QC_CONFIG = {
    "min_charge_states_for_reliable": 3,
    "min_charge_states_for_review": 2,
    "require_contiguous_charge_states": True,
    "max_neutral_mass_sd_da": 0.5,
    "max_neutral_mass_range_da": 1.5,
    "max_mass_error_ppm": 20,
    "max_envelope_internal_error_ppm": 20,
    "min_relative_intensity_percent": 0.5,
    "min_relative_envelope_intensity_percent_for_reliable": 1.0,
    "min_relative_envelope_intensity_percent_for_review": 0.1,
    "max_competing_envelopes": 3,
    "comparison_ready_statuses": ["Reliable", "Review"],
    "max_rt_range_min_for_reliable": 0.15,
    "max_rt_range_min_for_review": 0.30,
    "allow_trace_only_reliable": False,
    "search_mode": "untargeted",
    "reference_masses": [],
    "reference_mass_tolerance_ppm": 20,
    "neutral_mass_range": {"enabled": True, "min_da": 20000, "max_da": 30000},
    "target_review_mass_range": {"enabled": False, "min_da": None, "max_da": None},
    "envelope_grouping": {
        "enabled": True,
        "mass_tolerance_da": 1.0,
        "rt_tolerance_min": 0.15,
        "min_shared_peak_fraction": 0.5,
        "min_shared_charge_fraction": 0.5,
        "require_peak_overlap": True,
    },
    "engine": "legacy_cluster",
    "compare_with_legacy": False,
    "comparison_ready_tiers": {
        "strict": ["Tier_1_high_quality"],
        "review": ["Tier_1_high_quality", "Tier_2_supported"],
    },
    "quality_tiers": {
        "tier1_min_charge_states": 3,
        "tier1_min_consecutive_charge_states": 3,
        "tier1_min_local_envelope_relative_intensity_percent": 1.0,
        "tier2_min_charge_states": 2,
        "tier2_min_consecutive_charge_states": 2,
        "tier2_min_local_envelope_relative_intensity_percent": 0.1,
    },
    "engine_comparison": {
        "mass_tolerance_ppm": 20,
        "rt_tolerance_min": 0.15,
        "min_shared_charge_fraction": 0.5,
        "require_mass_match": True,
    },
    "competitive_assignment": {
        "enabled": True,
        "rt_tolerance_min": 0.15,
        "close_score_margin": 5.0,
        "dry_run": True,
        "min_independent_peak_fraction": 0.5,
        "min_independent_charge_states": 2,
        "allow_shared_peaks_between_selected": False,
        "minimum_score_margin_for_exclusive_selection": 1.0,
        "evidence_score_config_version": "MVP-5.9.8a-v1",
        "score_weights": {
            "charge_count": 12.0,
            "consecutive_charge_run": 10.0,
            "charge_coverage": 8.0,
            "local_relative_intensity": 5.0,
            "supporting_scan_count": 3.0,
            "rt_consistency": 5.0,
            "extension_support": 2.0,
            "split_merge_support": 2.0,
            "internal_error_penalty": 10.0,
            "neutral_mass_sd_penalty": 5.0,
            "neutral_mass_range_penalty": 5.0,
            "rt_range_penalty": 5.0,
            "peak_sharing_penalty": 12.0,
            "peak_usage_penalty": 5.0,
            "charge_gap_penalty": 8.0,
            "severe_limiting_factor_penalty": 15.0,
        },
    },
    "rt_localized": {
        "enabled": True,
        "rt_window_min": 0.10,
        "rt_step_min": 0.05,
        "min_scans_per_window": 1,
        "peak_aggregation": "max",
        "mz_merge_tolerance_ppm": 10,
        "adjacent_charge_mz_tolerance_ppm": 20,
        "max_charge_gap": 1,
        "min_charge_states": 2,
        "min_consecutive_charge_states": 2,
        "require_consecutive_for_candidate": True,
        "min_local_relative_peak_intensity_percent": 0.1,
        "neutral_mass_estimator": "intensity_weighted_mean",
        "merge_across_windows": {
            "enabled": True,
            "mass_tolerance_ppm": 10,
            "rt_overlap_required": True,
            "min_shared_charge_fraction": 0.5,
        },
        "charge_extension": {
            "enabled": True,
            "max_extension_charges": 2,
            "weak_peak_tolerance_ppm": 30,
            "weak_peak_min_local_relative_percent": 0.01,
            "add_weak_peaks_to_envelope": False,
        },
        "split_envelope_merge": {
            "enabled": True,
            "mass_tolerance_ppm": 20,
            "rt_tolerance_min": 0.10,
            "max_charge_gap": 1,
        },
        "peak_sharing": {
            "high_usage_threshold": 5,
            "max_highly_shared_fraction_for_tier1": 0.25,
            "max_highly_shared_fraction_for_tier2": 0.50,
        },
    },
    "mass_spectrum_output": {
        "enabled": True,
        "representatives_only": True,
        "comparison_ready_only": False,
        "include_qc_ineligible": True,
        "intensity_method": "total_supporting_intensity",
        "normalize_to_percent": True,
        "bin_width_da": None,
        "minimum_quality_tier": "Tier_3_weak",
    },
}

QC_COLUMNS = [
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
    "Shared_Peak_Count_With_Representative",
    "Shared_Peak_Fraction_With_Representative",
    "Shared_Charge_Count_With_Representative",
    "Shared_Charge_Fraction_With_Representative",
    "Mass_Delta_To_Group_Representative_Da",
    "RT_Delta_To_Group_Representative_Min",
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
    "Charge_State_Range",
    "Charge_State_Continuity",
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
]

DIAGNOSTIC_COLUMNS = [
    "Total_Reconstruction_Candidates",
    "Reliable_Count",
    "Review_Count",
    "Insufficient_Count",
    "Failed_Count",
    "Envelope_QC_Eligible_Count",
    "Intact_Strict_Eligible_Count",
    "Intact_Review_Eligible_Count",
    "Comparison_Ready_Strict_Count",
    "Comparison_Ready_Review_Count",
    "Comparison_Ready_Count",
    "Exact_Duplicate_Group_Count",
    "Exact_Duplicate_Candidate_Count",
    "Intact_Envelope_Group_Count",
    "Unique_Envelope_Group_Count",
    "Overlapping_Envelope_Group_Count",
    "Competing_Reconstruction_Group_Count",
    "Comparison_Representative_Count",
    "Target_Review_Representative_Count",
    "Candidates_Removed_As_Exact_Duplicates",
    "Candidates_Removed_As_Group_Nonrepresentatives",
    "Dominant_Comparison_Representative_Mass",
    "Dominant_Comparison_Representative_Intensity",
    "Dominant_Target_Review_Representative_Mass",
    "Dominant_Target_Review_Representative_Intensity",
    "Grouping_Mass_Tolerance_Da",
    "Grouping_RT_Tolerance_Min",
    "Grouping_Min_Shared_Peak_Fraction",
    "Grouping_Min_Shared_Charge_Fraction",
    "Trace_Only_Envelope_Count",
    "Noncontiguous_Envelope_Count",
    "RT_Inconsistent_Count",
    "Internal_Mass_Error_Count",
    "Theory_Near_Match_Count",
    "Reference_Match_Count",
    "Dominant_Envelope_Mass",
    "Dominant_Envelope_Intensity",
    "Dominant_Envelope_Status",
    "Dominant_Envelope_Comparison_Ready",
    "Dominant_Envelope_Overall_Mass",
    "Dominant_Envelope_Overall_Intensity",
    "Dominant_Envelope_In_Mass_Range_Mass",
    "Dominant_Envelope_In_Mass_Range_Intensity",
    "Dominant_Envelope_In_Mass_Range_Status",
    "Dominant_Envelope_In_Mass_Range_Comparison_Ready",
    "Dominant_Envelope_In_Search_Range_Raw_Mass",
    "Dominant_Envelope_In_Search_Range_Raw_Intensity",
    "Dominant_Intact_Strict_Envelope_Mass",
    "Dominant_Intact_Strict_Envelope_Intensity",
    "Dominant_Intact_Strict_QC_Score",
    "Dominant_Intact_Review_Envelope_Mass",
    "Dominant_Intact_Review_Envelope_Intensity",
    "Dominant_Intact_Review_QC_Score",
    "Dominant_Intact_Eligible_Envelope_Mass",
    "Dominant_Intact_Eligible_Envelope_Intensity",
    "Dominant_Intact_Eligible_QC_Score",
    "Dominant_Intact_Eligible_Reference_Label",
    "Failure_Reason_Counts",
    "Reconstruction_Enabled",
    "Reconstruction_Engine",
    "Num_RT_Windows",
    "Num_Local_Peaks",
    "Num_Anchor_Peaks_Evaluated",
    "Num_Raw_Envelope_Candidates",
    "Num_Candidates_After_Charge_Filter",
    "Num_Candidates_After_RT_Window_Merge",
    "Num_Candidates_With_Consecutive_Charges",
    "Num_Candidates_With_Charge_Gaps",
    "Num_Missing_Charges_Evaluated",
    "Num_Missing_Charges_With_Weak_Peaks",
    "Num_Missing_Charges_Not_Detected",
    "Median_RT_Range_Min",
    "Median_Internal_Error_ppm",
    "Median_Charge_Count",
    "Processing_Time_Seconds",
    "Candidate_Count_By_Quality_Tier",
    "Strict_Failure_Reason_Counts",
    "Review_Failure_Reason_Counts",
    "Charge_Count_Distribution",
    "Consecutive_Charge_Distribution",
    "Internal_Error_Distribution",
    "RT_Range_Distribution",
    "Peak_Usage_Distribution",
    "Split_Envelope_Count",
    "Missing_Charge_Status_Count",
    "Engine_Match_Status_Count",
    "Matched_Engine_Candidate_Count",
    "Legacy_Only_Count",
    "RT_Localized_Only_Count",
    "Mass_Mismatch_Count",
    "RT_Mismatch_Count",
    "Charge_Overlap_Failure_Count",
    "Median_Mass_Delta_For_Matched_Only",
    "Competing_Envelope_Group_Count",
    "MultiCandidate_Competition_Group_Count",
    "Noncompeting_Candidate_Count",
    "Largest_Competition_Group_Size",
    "Median_Competition_Group_Size",
    "Mean_Competition_Group_Size",
    "Competition_Group_Size_Distribution",
    "Median_Top_Evidence_Score",
    "Median_Top_vs_Second_Score_Margin",
    "Competition_Groups_With_Close_Scores",
    "High_Peak_Sharing_Candidate_Count",
    "Competitive_Scoring_Time_Seconds",
    "Dry_Run_Selected_Candidate_Count",
    "Dry_Run_Primary_Selected_Count",
    "Dry_Run_Independent_Selected_Count",
    "Dry_Run_Would_Exclude_Count",
    "Dry_Run_Ambiguous_Count",
    "Dry_Run_Noncompeting_Count",
    "Would_Exclude_Peak_Reuse_Count",
    "Would_Exclude_Independent_Peak_Shortage_Count",
    "Would_Exclude_Independent_Charge_Shortage_Count",
    "Median_Independent_Peak_Fraction",
    "Median_Direct_Competitor_Count",
    "Largest_Component_Selected_Count",
    "Largest_Component_Excluded_Count",
    "Largest_Component_Ambiguous_Count",
    "Assignment_Dry_Run_Time_Seconds",
    "Neutral_Mass_Search_Min_Da",
    "Neutral_Mass_Search_Max_Da",
    "Total_Candidates_Before_Mass_Range_Filter",
    "Total_Candidates_In_Mass_Range",
    "Total_Candidates_Outside_Mass_Range",
    "Target_Review_Mass_Range_Settings",
    "Target_Review_Candidate_Count",
    "Search_Mode",
    "Intensity_Normalization_Method",
    "RT_Tolerance_Settings",
    "Reference_Masses_Used",
    "Min_Charge_States_For_Reliable",
    "Min_Charge_States_For_Review",
    "Require_Contiguous_Charge_States",
    "Max_Neutral_Mass_SD_Da",
    "Max_Neutral_Mass_Range_Da",
    "Max_Envelope_Internal_Error_ppm",
    "Min_Relative_Intensity_Percent",
    "Min_Relative_Envelope_Intensity_Percent_For_Reliable",
    "Min_Relative_Envelope_Intensity_Percent_For_Review",
    "Max_Competing_Envelopes",
    "Comparison_Ready_Statuses",
    "Notes",
]

GROUP_COLUMNS = [
    "Intact_Envelope_Group_ID",
    "Group_Size",
    "Exact_Duplicate_Count",
    "Representative_Cluster_ID",
    "Representative_Mass",
    "Representative_Status",
    "Representative_QC_Score",
    "Representative_Comparison_Ready",
    "Representative_Total_Intensity",
    "Representative_Reconstructed_Envelope_Intensity",
    "Intensity_Method",
    "Representative_Charge_States",
    "Representative_RT_Range",
    "Group_Mass_Min",
    "Group_Mass_Max",
    "Group_Mass_Range",
    "Group_RT_Min",
    "Group_RT_Max",
    "Group_Ambiguity_Status",
    "Member_Cluster_IDs",
    "Notes",
]

COMPARISON_CANDIDATE_COLUMNS = [
    "Comparison_Representative_Rank",
    "Cluster_ID",
    "Reconstructed_Mass",
    "Reconstruction_Engine",
    "Intact_Envelope_Group_ID",
    "Reconstruction_Status",
    "Comparison_Ready",
    "Competing_Envelope_Group_ID",
    "Competing_Envelope_Group_Size",
    "Envelope_Evidence_Score",
    "Evidence_Score_Rank_In_Competition",
    "Maximum_Shared_Peak_Fraction",
    "Intact_Envelope_QC_Score",
    "Total_Supporting_Intensity",
    "Reconstructed_Envelope_Intensity",
    "Intensity_Method",
    "Relative_Intact_Eligible_Intensity_Percent",
    "Supporting_Charge_States",
    "Charge_State_Continuity",
    "RT_Range_Min",
    "Envelope_Internal_Error_ppm",
    "Neutral_Mass_SD",
    "Neutral_Mass_Range",
    "Unmodified_Theory_Delta_Da",
    "Best_Reference_Label",
    "Reference_Mass_Error_ppm",
    "Target_Review_Rank",
    "Limiting_Factors",
]

TARGET_REVIEW_CANDIDATE_COLUMNS = [
    "Target_Review_Rank",
    "Cluster_ID",
    "Reconstructed_Mass",
    "Reconstruction_Engine",
    "Intact_Envelope_Group_ID",
    "Comparison_Representative_Rank",
    "Reconstruction_Status",
    "Comparison_Ready",
    "Intact_Envelope_QC_Score",
    "Total_Supporting_Intensity",
    "Reconstructed_Envelope_Intensity",
    "Intensity_Method",
    "Supporting_Charge_States",
    "Charge_State_Continuity",
    "RT_Range_Min",
    "Envelope_Internal_Error_ppm",
    "Neutral_Mass_SD",
    "Neutral_Mass_Range",
    "Best_Reference_Label",
    "Reference_Mass_Error_ppm",
    "Limiting_Factors",
]

RECONSTRUCTED_MASS_SPECTRUM_COLUMNS = [
    "Spectrum_Point_Rank",
    "Reconstructed_Mass_Da",
    "Reconstructed_Envelope_Intensity",
    "Relative_Intensity_Percent",
    "Intensity_Method",
    "Cluster_ID",
    "Reconstruction_Engine",
    "Intact_Envelope_Group_ID",
    "Group_Representative",
    "Comparison_Representative",
    "Reconstruction_Status",
    "Intact_Quality_Tier",
    "Quality_Tier_Rank",
    "Peak_Sharing_Status",
    "Max_Peak_Usage_Count",
    "Competing_Envelope_Group_ID",
    "Competing_Envelope_Group_Size",
    "Envelope_Evidence_Score",
    "Evidence_Score_Rank_In_Competition",
    "Maximum_Shared_Peak_Fraction",
    "Envelope_QC_Eligible",
    "Intact_Strict_Eligible",
    "Intact_Review_Eligible",
    "Comparison_Ready",
    "Num_Supporting_Charge_States",
    "Supporting_Charge_States",
    "RT_Mean",
    "RT_Range_Min",
    "Envelope_Internal_Error_ppm",
    "Neutral_Mass_SD",
    "Neutral_Mass_Range",
    "Best_Reference_Label",
    "Reference_Mass_Error_ppm",
    "Limiting_Factors",
]

RT_ENVELOPE_DIAGNOSTIC_COLUMNS = [
    "Cluster_ID",
    "Reconstruction_Engine",
    "RT_Window_ID",
    "RT_Window_Start_Min",
    "RT_Window_End_Min",
    "RT_Window_Center_Min",
    "Num_MS1_Scans_In_Window",
    "Peak_Aggregation_Method",
    "Anchor_MZ",
    "Anchor_Charge",
    "Reconstructed_Mass",
    "Predicted_Charge_States",
    "Observed_Charge_States",
    "Missing_Charge_States",
    "Num_Predicted_Charges",
    "Num_Observed_Charges",
    "Charge_Coverage_Fraction",
    "Consecutive_Charge_Run_Length",
    "Longest_Consecutive_Charge_Run",
    "Charge_Gap_Count",
    "Charge_Continuity_Fraction",
    "Local_Window_Max_Intensity",
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
    "Split_Envelope_Group_ID",
    "Split_Envelope_Member_Count",
    "Split_Envelope_Merged",
    "Charge_Gaps_Before_Merge",
    "Charge_Gaps_After_Merge",
    "Max_Peak_Usage_Count",
    "Mean_Peak_Usage_Count",
    "Num_Highly_Shared_Peaks",
    "Highly_Shared_Peak_Fraction",
    "Peak_Sharing_Status",
    "Competing_Envelope_Group_ID",
    "Competing_Envelope_Group_Size",
    "Envelope_Evidence_Score",
    "Evidence_Score_Rank_In_Competition",
    "Maximum_Shared_Peak_Fraction",
    "Direct_Competitor_Count",
    "Dry_Run_Assignment_Status",
    "Dry_Run_Selected",
    "Independent_Supporting_Peak_Count",
    "Independent_Supporting_Peak_Fraction",
    "Independent_Charge_State_Count",
    "Assignment_Confidence",
    "Notes",
]

MISSING_CHARGE_DIAGNOSTIC_COLUMNS = [
    "Cluster_ID",
    "RT_Window_ID",
    "Reconstructed_Mass",
    "Missing_Charge",
    "Predicted_MZ",
    "Nearest_Observed_MZ",
    "Error_ppm",
    "Nearest_Intensity",
    "Detection_Status",
    "Notes",
]

ENGINE_COMPARISON_COLUMNS = [
    "Legacy_Cluster_ID",
    "RT_Localized_Cluster_ID",
    "Legacy_Mass",
    "RT_Localized_Mass",
    "Mass_Delta_Da",
    "Legacy_Charge_Count",
    "RT_Localized_Charge_Count",
    "Legacy_RT_Range",
    "RT_Localized_RT_Range",
    "Legacy_Internal_Error_ppm",
    "RT_Localized_Internal_Error_ppm",
    "Peak_Overlap_Fraction",
    "Engine_Match_Status",
    "Match_Failure_Reason",
    "Mass_Matched",
    "RT_Matched",
    "Charge_Matched",
    "Peak_Matched",
    "Notes",
]

RT_ENGINE_QC_SUMMARY_COLUMNS = [
    "Metric",
    "Value",
    "Notes",
]

COMPETITION_GROUP_COLUMNS = [
    "Competing_Envelope_Group_ID",
    "Group_Size",
    "MultiCandidate_Group",
    "Top_Ranked_Cluster_ID",
    "Top_Ranked_Mass",
    "Top_Evidence_Score",
    "Second_Evidence_Score",
    "Score_Margin_Top_vs_Second",
    "Shared_Local_Peak_Count",
    "Maximum_Peak_Usage_Count",
    "Member_Cluster_IDs",
    "Member_Masses",
    "Member_Quality_Tiers",
    "Notes",
]

COMPETITION_SCORE_COLUMNS = [
    "Cluster_ID",
    "Reconstructed_Mass",
    "Competing_Envelope_Group_ID",
    "Competing_Envelope_Group_Size",
    "Envelope_Evidence_Score",
    "Evidence_Score_Rank_In_Competition",
    "Evidence_Score_Components",
    "Evidence_Score_Penalties",
    "Intact_Quality_Tier",
    "Supporting_Charge_Count",
    "Longest_Consecutive_Charge_Run",
    "Charge_Coverage_Fraction",
    "Local_Envelope_Relative_Intensity_Percent",
    "Envelope_Internal_Error_ppm",
    "Neutral_Mass_SD",
    "Neutral_Mass_Range",
    "RT_Range_Min",
    "Maximum_Shared_Peak_Fraction",
    "Max_Peak_Usage_Count",
    "Direct_Competitor_Count",
    "Dry_Run_Assignment_Status",
    "Dry_Run_Selected",
    "Independent_Supporting_Peak_Count",
    "Independent_Supporting_Peak_Fraction",
    "Independent_Charge_State_Count",
    "Assignment_Confidence",
    "Limiting_Factors",
]

ASSIGNMENT_DRY_RUN_COLUMNS = [
    "Cluster_ID",
    "Reconstructed_Mass",
    "Competing_Envelope_Group_ID",
    "Competition_Component_Size",
    "Direct_Competitor_Count",
    "Envelope_Evidence_Score",
    "Evidence_Score_Rank_In_Competition",
    "Dry_Run_Assignment_Status",
    "Dry_Run_Selected",
    "Dry_Run_Selection_Order",
    "Supporting_Peak_Count_Before_Assignment",
    "Independent_Supporting_Peak_Count",
    "Independent_Supporting_Peak_Fraction",
    "Supporting_Charge_Count_Before_Assignment",
    "Independent_Charge_State_Count",
    "Peaks_Already_Assigned_Count",
    "Excluded_By_Cluster_ID",
    "Score_Margin_To_Excluding_Candidate",
    "Close_Score_Ambiguity",
    "Assignment_Confidence",
    "Intact_Quality_Tier",
    "Comparison_Ready",
    "Dry_Run_Exclusion_Reason",
    "Notes",
]

ASSIGNMENT_DRY_RUN_SUMMARY_COLUMNS = [
    "Competing_Envelope_Group_ID",
    "Component_Size",
    "Unique_Local_Peak_Count",
    "Direct_Competition_Edge_Count",
    "Selected_Primary_Count",
    "Selected_Independent_Count",
    "Noncompeting_Count",
    "Would_Exclude_Count",
    "Ambiguous_Count",
    "Candidates_Remaining_After_Dry_Run",
    "Primary_Selected_Cluster_IDs",
    "Main_Exclusion_Reasons",
    "Notes",
]

SEVERE_LIMITING_FACTORS = {
    "reconstruction_disabled",
    "no_charge_state_candidates",
    "insufficient_charge_states",
    "mass_spread_too_large",
    "internal_mass_error_too_large",
    "rt_inconsistent",
    "insufficient_intensity_support",
    "multiple_competing_envelopes",
}


def _as_bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _as_reference_masses(value: Any) -> list[dict[str, Any]]:
    references = []
    if not value:
        return references
    raw_items = value if isinstance(value, list) else [value]
    for index, item in enumerate(raw_items, start=1):
        if isinstance(item, dict):
            mass = item.get("mass_da") or item.get("mass") or item.get("Mass_Da")
            label = item.get("label") or item.get("name") or f"reference_{index}"
        else:
            mass = item
            label = f"reference_{index}"
        try:
            references.append({"label": str(label), "mass_da": float(mass)})
        except (TypeError, ValueError):
            continue
    return references


def _optional_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _qc_config(reconstruction_config: dict[str, Any]) -> dict[str, Any]:
    raw = reconstruction_config.get("intact_reconstruction") or reconstruction_config.get("qc") or {}
    merged = {**DEFAULT_INTACT_QC_CONFIG, **raw}
    merged["min_charge_states_for_reliable"] = int(merged.get("min_charge_states_for_reliable") or 3)
    merged["min_charge_states_for_review"] = int(merged.get("min_charge_states_for_review") or 2)
    merged["require_contiguous_charge_states"] = _as_bool(merged.get("require_contiguous_charge_states"), True)
    merged["max_neutral_mass_sd_da"] = float(merged.get("max_neutral_mass_sd_da") or 0.5)
    merged["max_neutral_mass_range_da"] = float(merged.get("max_neutral_mass_range_da") or 1.5)
    merged["max_mass_error_ppm"] = float(merged.get("max_mass_error_ppm") or 20)
    merged["max_envelope_internal_error_ppm"] = float(
        merged.get("max_envelope_internal_error_ppm") or merged["max_mass_error_ppm"]
    )
    merged["min_relative_intensity_percent"] = float(merged.get("min_relative_intensity_percent") or 0.5)
    merged["min_relative_envelope_intensity_percent_for_reliable"] = float(
        merged.get("min_relative_envelope_intensity_percent_for_reliable") or 1.0
    )
    merged["min_relative_envelope_intensity_percent_for_review"] = float(
        merged.get("min_relative_envelope_intensity_percent_for_review") or 0.1
    )
    merged["max_competing_envelopes"] = int(merged.get("max_competing_envelopes") or 3)
    merged["max_rt_range_min_for_reliable"] = float(merged.get("max_rt_range_min_for_reliable") or 0.15)
    merged["max_rt_range_min_for_review"] = float(merged.get("max_rt_range_min_for_review") or 0.30)
    merged["allow_trace_only_reliable"] = _as_bool(merged.get("allow_trace_only_reliable"), False)
    search_mode = str(merged.get("search_mode") or "untargeted")
    merged["search_mode"] = search_mode if search_mode in {"untargeted", "theoretical_targeted", "reference_targeted"} else "untargeted"
    statuses = merged.get("comparison_ready_statuses") or ["Reliable", "Review"]
    if isinstance(statuses, str):
        statuses = [item.strip() for item in statuses.split(",") if item.strip()]
    merged["comparison_ready_statuses"] = list(statuses)
    merged["reference_masses"] = _as_reference_masses(merged.get("reference_masses"))
    merged["reference_mass_tolerance_ppm"] = float(merged.get("reference_mass_tolerance_ppm") or 20)
    neutral_range = merged.get("neutral_mass_range") or {}
    if not isinstance(neutral_range, dict):
        neutral_range = {}
    merged["neutral_mass_range"] = {
        "enabled": _as_bool(neutral_range.get("enabled"), True),
        "min_da": float(neutral_range.get("min_da", 20000) if neutral_range.get("min_da", None) is not None else 20000),
        "max_da": float(neutral_range.get("max_da", 30000) if neutral_range.get("max_da", None) is not None else 30000),
    }
    if merged["neutral_mass_range"]["min_da"] > merged["neutral_mass_range"]["max_da"]:
        merged["neutral_mass_range"]["min_da"], merged["neutral_mass_range"]["max_da"] = (
            merged["neutral_mass_range"]["max_da"],
            merged["neutral_mass_range"]["min_da"],
        )
    target_range = merged.get("target_review_mass_range") or {}
    if not isinstance(target_range, dict):
        target_range = {}
    merged["target_review_mass_range"] = {
        "enabled": _as_bool(target_range.get("enabled"), False),
        "min_da": _optional_float(target_range.get("min_da")),
        "max_da": _optional_float(target_range.get("max_da")),
    }
    target_min = merged["target_review_mass_range"]["min_da"]
    target_max = merged["target_review_mass_range"]["max_da"]
    if target_min is not None and target_max is not None and target_min > target_max:
        merged["target_review_mass_range"]["min_da"], merged["target_review_mass_range"]["max_da"] = (
            target_max,
            target_min,
        )
    engine = str(merged.get("engine") or "legacy_cluster")
    merged["engine"] = engine if engine in {"legacy_cluster", "rt_localized"} else "legacy_cluster"
    merged["compare_with_legacy"] = _as_bool(merged.get("compare_with_legacy"), False)
    rt_localized = merged.get("rt_localized") or {}
    if not isinstance(rt_localized, dict):
        rt_localized = {}
    peak_aggregation = str(rt_localized.get("peak_aggregation") or "max")
    if peak_aggregation not in {"max", "sum", "mean"}:
        peak_aggregation = "max"
    estimator = str(rt_localized.get("neutral_mass_estimator") or "intensity_weighted_mean")
    if estimator not in {"unweighted_mean", "intensity_weighted_mean", "median"}:
        estimator = "intensity_weighted_mean"
    merge_cfg = rt_localized.get("merge_across_windows") or {}
    if not isinstance(merge_cfg, dict):
        merge_cfg = {}
    extension_cfg = rt_localized.get("charge_extension") or {}
    if not isinstance(extension_cfg, dict):
        extension_cfg = {}
    split_cfg = rt_localized.get("split_envelope_merge") or {}
    if not isinstance(split_cfg, dict):
        split_cfg = {}
    sharing_cfg = rt_localized.get("peak_sharing") or {}
    if not isinstance(sharing_cfg, dict):
        sharing_cfg = {}
    merged["rt_localized"] = {
        "enabled": _as_bool(rt_localized.get("enabled"), True),
        "rt_window_min": float(rt_localized.get("rt_window_min") if rt_localized.get("rt_window_min") is not None else 0.10),
        "rt_step_min": float(rt_localized.get("rt_step_min") if rt_localized.get("rt_step_min") is not None else 0.05),
        "min_scans_per_window": int(rt_localized.get("min_scans_per_window") or 1),
        "peak_aggregation": peak_aggregation,
        "mz_merge_tolerance_ppm": float(rt_localized.get("mz_merge_tolerance_ppm") if rt_localized.get("mz_merge_tolerance_ppm") is not None else 10),
        "adjacent_charge_mz_tolerance_ppm": float(rt_localized.get("adjacent_charge_mz_tolerance_ppm") if rt_localized.get("adjacent_charge_mz_tolerance_ppm") is not None else 20),
        "max_charge_gap": int(rt_localized.get("max_charge_gap") or 1),
        "min_charge_states": int(rt_localized.get("min_charge_states") or 2),
        "min_consecutive_charge_states": int(rt_localized.get("min_consecutive_charge_states") or 2),
        "require_consecutive_for_candidate": _as_bool(rt_localized.get("require_consecutive_for_candidate"), True),
        "min_local_relative_peak_intensity_percent": float(rt_localized.get("min_local_relative_peak_intensity_percent") if rt_localized.get("min_local_relative_peak_intensity_percent") is not None else 0.1),
        "neutral_mass_estimator": estimator,
        "merge_across_windows": {
            "enabled": _as_bool(merge_cfg.get("enabled"), True),
            "mass_tolerance_ppm": float(merge_cfg.get("mass_tolerance_ppm") if merge_cfg.get("mass_tolerance_ppm") is not None else 10),
            "rt_overlap_required": _as_bool(merge_cfg.get("rt_overlap_required"), True),
            "min_shared_charge_fraction": float(merge_cfg.get("min_shared_charge_fraction") if merge_cfg.get("min_shared_charge_fraction") is not None else 0.5),
        },
        "charge_extension": {
            "enabled": _as_bool(extension_cfg.get("enabled"), True),
            "max_extension_charges": int(extension_cfg.get("max_extension_charges") or 2),
            "weak_peak_tolerance_ppm": float(extension_cfg.get("weak_peak_tolerance_ppm") if extension_cfg.get("weak_peak_tolerance_ppm") is not None else 30),
            "weak_peak_min_local_relative_percent": float(extension_cfg.get("weak_peak_min_local_relative_percent") if extension_cfg.get("weak_peak_min_local_relative_percent") is not None else 0.01),
            "add_weak_peaks_to_envelope": _as_bool(extension_cfg.get("add_weak_peaks_to_envelope"), False),
        },
        "split_envelope_merge": {
            "enabled": _as_bool(split_cfg.get("enabled"), True),
            "mass_tolerance_ppm": float(split_cfg.get("mass_tolerance_ppm") if split_cfg.get("mass_tolerance_ppm") is not None else 20),
            "rt_tolerance_min": float(split_cfg.get("rt_tolerance_min") if split_cfg.get("rt_tolerance_min") is not None else 0.10),
            "max_charge_gap": int(split_cfg.get("max_charge_gap") or 1),
        },
        "peak_sharing": {
            "high_usage_threshold": int(sharing_cfg.get("high_usage_threshold") or 5),
            "max_highly_shared_fraction_for_tier1": float(sharing_cfg.get("max_highly_shared_fraction_for_tier1") if sharing_cfg.get("max_highly_shared_fraction_for_tier1") is not None else 0.25),
            "max_highly_shared_fraction_for_tier2": float(sharing_cfg.get("max_highly_shared_fraction_for_tier2") if sharing_cfg.get("max_highly_shared_fraction_for_tier2") is not None else 0.50),
        },
    }
    if merged["rt_localized"]["rt_window_min"] <= 0:
        merged["rt_localized"]["rt_window_min"] = 0.10
    if merged["rt_localized"]["rt_step_min"] <= 0:
        merged["rt_localized"]["rt_step_min"] = merged["rt_localized"]["rt_window_min"]

    grouping = merged.get("envelope_grouping") or {}
    if not isinstance(grouping, dict):
        grouping = {}
    merged["envelope_grouping"] = {
        "enabled": _as_bool(grouping.get("enabled"), True),
        "mass_tolerance_da": float(grouping.get("mass_tolerance_da") if grouping.get("mass_tolerance_da") is not None else 1.0),
        "rt_tolerance_min": float(grouping.get("rt_tolerance_min") if grouping.get("rt_tolerance_min") is not None else 0.15),
        "min_shared_peak_fraction": float(grouping.get("min_shared_peak_fraction") if grouping.get("min_shared_peak_fraction") is not None else 0.5),
        "min_shared_charge_fraction": float(grouping.get("min_shared_charge_fraction") if grouping.get("min_shared_charge_fraction") is not None else 0.5),
        "require_peak_overlap": _as_bool(grouping.get("require_peak_overlap"), True),
    }
    tier_cfg = merged.get("comparison_ready_tiers") or {}
    if not isinstance(tier_cfg, dict):
        tier_cfg = {}
    strict_tiers = tier_cfg.get("strict") or ["Tier_1_high_quality"]
    review_tiers = tier_cfg.get("review") or ["Tier_1_high_quality", "Tier_2_supported"]
    if isinstance(strict_tiers, str):
        strict_tiers = [strict_tiers]
    if isinstance(review_tiers, str):
        review_tiers = [review_tiers]
    merged["comparison_ready_tiers"] = {"strict": list(strict_tiers), "review": list(review_tiers)}
    quality_cfg = merged.get("quality_tiers") or {}
    if not isinstance(quality_cfg, dict):
        quality_cfg = {}
    merged["quality_tiers"] = {
        "tier1_min_charge_states": int(quality_cfg.get("tier1_min_charge_states") or 3),
        "tier1_min_consecutive_charge_states": int(quality_cfg.get("tier1_min_consecutive_charge_states") or 3),
        "tier1_min_local_envelope_relative_intensity_percent": float(quality_cfg.get("tier1_min_local_envelope_relative_intensity_percent") if quality_cfg.get("tier1_min_local_envelope_relative_intensity_percent") is not None else 1.0),
        "tier2_min_charge_states": int(quality_cfg.get("tier2_min_charge_states") or 2),
        "tier2_min_consecutive_charge_states": int(quality_cfg.get("tier2_min_consecutive_charge_states") or 2),
        "tier2_min_local_envelope_relative_intensity_percent": float(quality_cfg.get("tier2_min_local_envelope_relative_intensity_percent") if quality_cfg.get("tier2_min_local_envelope_relative_intensity_percent") is not None else 0.1),
    }
    engine_compare = merged.get("engine_comparison") or {}
    if not isinstance(engine_compare, dict):
        engine_compare = {}
    merged["engine_comparison"] = {
        "mass_tolerance_ppm": float(engine_compare.get("mass_tolerance_ppm") if engine_compare.get("mass_tolerance_ppm") is not None else 20),
        "rt_tolerance_min": float(engine_compare.get("rt_tolerance_min") if engine_compare.get("rt_tolerance_min") is not None else 0.15),
        "min_shared_charge_fraction": float(engine_compare.get("min_shared_charge_fraction") if engine_compare.get("min_shared_charge_fraction") is not None else 0.5),
        "require_mass_match": _as_bool(engine_compare.get("require_mass_match"), True),
    }
    competitive = merged.get("competitive_assignment") or {}
    if not isinstance(competitive, dict):
        competitive = {}
    default_competitive = DEFAULT_INTACT_QC_CONFIG["competitive_assignment"]
    default_weights = default_competitive["score_weights"]
    raw_weights = competitive.get("score_weights") or {}
    if not isinstance(raw_weights, dict):
        raw_weights = {}
    score_weights = {key: float(raw_weights.get(key) if raw_weights.get(key) is not None else value) for key, value in default_weights.items()}
    merged["competitive_assignment"] = {
        "enabled": _as_bool(competitive.get("enabled"), True),
        "rt_tolerance_min": float(competitive.get("rt_tolerance_min") if competitive.get("rt_tolerance_min") is not None else default_competitive["rt_tolerance_min"]),
        "close_score_margin": float(competitive.get("close_score_margin") if competitive.get("close_score_margin") is not None else default_competitive["close_score_margin"]),
        "dry_run": _as_bool(competitive.get("dry_run"), True),
        "min_independent_peak_fraction": float(competitive.get("min_independent_peak_fraction") if competitive.get("min_independent_peak_fraction") is not None else default_competitive.get("min_independent_peak_fraction", 0.5)),
        "min_independent_charge_states": int(competitive.get("min_independent_charge_states") or default_competitive.get("min_independent_charge_states", 2)),
        "allow_shared_peaks_between_selected": _as_bool(competitive.get("allow_shared_peaks_between_selected"), False),
        "minimum_score_margin_for_exclusive_selection": float(competitive.get("minimum_score_margin_for_exclusive_selection") if competitive.get("minimum_score_margin_for_exclusive_selection") is not None else default_competitive.get("minimum_score_margin_for_exclusive_selection", 1.0)),
        "evidence_score_config_version": str(competitive.get("evidence_score_config_version") or default_competitive["evidence_score_config_version"]),
        "score_weights": score_weights,
    }
    spectrum_output = merged.get("mass_spectrum_output") or {}
    if not isinstance(spectrum_output, dict):
        spectrum_output = {}
    intensity_method = str(spectrum_output.get("intensity_method") or "total_supporting_intensity")
    if intensity_method not in {"total_supporting_intensity", "mean_supporting_intensity", "max_supporting_intensity"}:
        raise ValueError("intact_reconstruction.mass_spectrum_output.intensity_method must be one of: total_supporting_intensity, mean_supporting_intensity, max_supporting_intensity")
    merged["mass_spectrum_output"] = {
        "enabled": _as_bool(spectrum_output.get("enabled"), True),
        "representatives_only": _as_bool(spectrum_output.get("representatives_only"), True),
        "comparison_ready_only": _as_bool(spectrum_output.get("comparison_ready_only"), False),
        "include_qc_ineligible": _as_bool(spectrum_output.get("include_qc_ineligible"), True),
        "intensity_method": intensity_method,
        "normalize_to_percent": _as_bool(spectrum_output.get("normalize_to_percent"), True),
        "bin_width_da": _optional_float(spectrum_output.get("bin_width_da")),
        "minimum_quality_tier": str(spectrum_output.get("minimum_quality_tier") or "Tier_3_weak"),
    }
    return merged


def _charge_state_range(charges: list[int]) -> str:
    if not charges:
        return ""
    if len(charges) == 1:
        return str(charges[0])
    return f"{min(charges)}-{max(charges)}"


def _charge_continuity(charges: list[int]) -> str:
    if not charges:
        return "missing"
    expected = set(range(min(charges), max(charges) + 1))
    return "contiguous" if set(charges) == expected else "non_contiguous"


def _primary_factor(factors: list[str]) -> str:
    priority = [
        "reconstruction_disabled",
        "no_charge_state_candidates",
        "insufficient_charge_states",
        "internal_mass_error_too_large",
        "mass_spread_too_large",
        "rt_inconsistent",
        "trace_only_envelope",
        "insufficient_intensity_support",
        "non_contiguous_charge_states",
        "multiple_competing_envelopes",
        "rt_not_available",
    ]
    for item in priority:
        if item in factors:
            return item
    return factors[0] if factors else ""


def _ppm(delta: float | None, reference: float | None) -> float | None:
    if delta is None or not reference:
        return None
    return delta / reference * 1_000_000


def _max_abs_ppm(values: list[float], center: float | None) -> float | None:
    if not values or not center:
        return None
    return max(abs(value - center) / center * 1_000_000 for value in values)


def _rt_metrics(cluster_peaks: list[dict[str, Any]], qc_config: dict[str, Any]) -> tuple[float | None, float | None, float | None, float | None, str]:
    rts = []
    for row in cluster_peaks:
        value = row.get("RT")
        if value is None or value == "":
            continue
        try:
            rts.append(float(value))
        except (TypeError, ValueError):
            continue
    if not rts:
        return None, None, None, None, "not_available"
    rt_min = min(rts)
    rt_max = max(rts)
    rt_mean = sum(rts) / len(rts)
    rt_range = rt_max - rt_min
    if rt_range <= qc_config["max_rt_range_min_for_reliable"]:
        consistency = "consistent"
    elif rt_range <= qc_config["max_rt_range_min_for_review"]:
        consistency = "review"
    else:
        consistency = "inconsistent"
    return rt_min, rt_max, rt_mean, rt_range, consistency


def _neutral_mass_range_status(observed_mass: float, qc_config: dict[str, Any]) -> tuple[bool, float, float, str]:
    neutral_range = qc_config["neutral_mass_range"]
    min_da = float(neutral_range["min_da"])
    max_da = float(neutral_range["max_da"])
    if not neutral_range.get("enabled", True):
        return True, min_da, max_da, "not_applied"
    in_range = min_da <= observed_mass <= max_da
    return in_range, min_da, max_da, "in_range" if in_range else "outside_range"


def _target_review_range_status(observed_mass: float, qc_config: dict[str, Any]) -> tuple[bool, str, str]:
    target_range = qc_config["target_review_mass_range"]
    if not target_range.get("enabled", False):
        return False, "not_configured", "not_configured"
    if target_range.get("min_da") is None or target_range.get("max_da") is None:
        return False, "not_configured", "not_configured"
    min_da = float(target_range["min_da"])
    max_da = float(target_range["max_da"])
    in_range = min_da <= observed_mass <= max_da
    status = "in_range" if in_range else "outside_range"
    priority = "target_review" if in_range else "outside_target_review"
    return in_range, status, priority


def _target_review_settings(qc_config: dict[str, Any]) -> str:
    target_range = qc_config["target_review_mass_range"]
    if not target_range.get("enabled", False):
        return "disabled"
    if target_range.get("min_da") is None or target_range.get("max_da") is None:
        return "not_configured"
    return f"enabled:{target_range['min_da']}-{target_range['max_da']} Da"


def _rt_rank_score(value: str) -> int:
    return {"consistent": 2, "review": 1}.get(str(value or ""), 0)


def _safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _reconstructed_intensity(total_intensity: Any, mean_intensity: Any, max_intensity: Any, method: str) -> float:
    if method == "mean_supporting_intensity":
        return _safe_float(mean_intensity)
    if method == "max_supporting_intensity":
        return _safe_float(max_intensity)
    return _safe_float(total_intensity)


def _small_metric_score(value: Any, limit: float, points: float) -> float:
    if value is None or value == "":
        return points
    numeric = _safe_float(value, limit)
    if limit <= 0:
        return points
    return max(0.0, points * (1.0 - min(numeric / limit, 1.0)))




def _tier_rank(tier: str) -> int:
    return {"Tier_1_high_quality": 1, "Tier_2_supported": 2, "Tier_3_weak": 3, "Tier_4_rejected": 4}.get(str(tier or ""), 4)


def _tier_allowed(row_tier: str, minimum_tier: str) -> bool:
    if str(minimum_tier or "").lower() == "all":
        return True
    return _tier_rank(row_tier) <= _tier_rank(minimum_tier)


def _reason_summary(values: list[str]) -> str:
    counts: dict[str, int] = {}
    for value in values:
        for item in str(value or "").split(";"):
            item = item.strip()
            if item:
                counts[item] = counts.get(item, 0) + 1
    return "; ".join(f"{key}:{counts[key]}" for key in sorted(counts))


def _distribution_summary(values: list[Any]) -> str:
    counts: dict[str, int] = {}
    for value in values:
        key = str(value if value not in {None, ""} else "blank")
        counts[key] = counts.get(key, 0) + 1
    return "; ".join(f"{key}:{counts[key]}" for key in sorted(counts))


def _binned_distribution(values: list[Any], bins: list[tuple[float, str]]) -> str:
    counts = {label: 0 for _, label in bins}
    counts["above"] = 0
    for value in values:
        numeric = _safe_float(value, None)
        if numeric is None:
            continue
        placed = False
        for limit, label in bins:
            if numeric <= limit:
                counts[label] += 1
                placed = True
                break
        if not placed:
            counts["above"] += 1
    return "; ".join(f"{key}:{value}" for key, value in counts.items() if value)


def _quality_matrix_and_tier(
    *,
    charges: list[int],
    continuity: str,
    neutral_sd: float | None,
    neutral_range: float | None,
    envelope_internal_error_ppm: float | None,
    rt_consistency: str,
    local_envelope_relative: float,
    competing: int,
    severe_factors: list[str],
    reconstruction_enabled: bool,
    in_mass_range: bool,
    qc_config: dict[str, Any],
    candidate: IntactMassCandidate,
) -> dict[str, Any]:
    quality_cfg = qc_config["quality_tiers"]
    sharing_cfg = qc_config["rt_localized"]["peak_sharing"]
    charge_count = len(charges)
    longest_run = int(_candidate_extra(candidate, "longest_consecutive_charge_run", 0) or 0)
    if longest_run <= 0:
        longest_run = _longest_consecutive_run(charges)
    highly_shared_fraction = float(_candidate_extra(candidate, "highly_shared_peak_fraction", _candidate_extra(candidate, "shared_peak_fraction", 0.0)) or 0.0)
    pass_min_charge = charge_count >= qc_config["min_charge_states_for_review"]
    pass_min_consecutive = longest_run >= quality_cfg["tier2_min_consecutive_charge_states"]
    pass_continuity = continuity == "contiguous" or not qc_config["require_contiguous_charge_states"]
    pass_internal = envelope_internal_error_ppm is None or envelope_internal_error_ppm <= qc_config["max_envelope_internal_error_ppm"]
    pass_sd = neutral_sd is None or neutral_sd <= qc_config["max_neutral_mass_sd_da"]
    pass_range = neutral_range is None or neutral_range <= qc_config["max_neutral_mass_range_da"]
    pass_rt = rt_consistency in {"consistent", "review"}
    pass_local = local_envelope_relative >= quality_cfg["tier2_min_local_envelope_relative_intensity_percent"]
    pass_competing = competing <= qc_config["max_competing_envelopes"]
    pass_peak_sharing = highly_shared_fraction <= sharing_cfg["max_highly_shared_fraction_for_tier2"]
    strict_checks = {
        "min_charge_count": charge_count >= quality_cfg["tier1_min_charge_states"],
        "min_consecutive_charge_count": longest_run >= quality_cfg["tier1_min_consecutive_charge_states"],
        "charge_continuity": pass_continuity,
        "internal_error": pass_internal,
        "neutral_mass_sd": pass_sd,
        "neutral_mass_range": pass_range,
        "rt_consistency": rt_consistency == "consistent",
        "local_intensity": local_envelope_relative >= quality_cfg["tier1_min_local_envelope_relative_intensity_percent"],
        "competing_envelope": pass_competing,
        "peak_sharing": highly_shared_fraction <= sharing_cfg["max_highly_shared_fraction_for_tier1"],
        "severe_warning": not severe_factors,
        "in_neutral_mass_range": in_mass_range,
        "reconstruction_enabled": reconstruction_enabled,
    }
    review_checks = {
        "min_charge_count": pass_min_charge,
        "min_consecutive_charge_count": pass_min_consecutive,
        "charge_continuity": pass_continuity,
        "internal_error": pass_internal,
        "neutral_mass_sd": pass_sd,
        "neutral_mass_range": pass_range,
        "rt_consistency": pass_rt,
        "local_intensity": pass_local,
        "competing_envelope": pass_competing,
        "peak_sharing": pass_peak_sharing,
        "severe_warning": not severe_factors,
        "in_neutral_mass_range": in_mass_range,
        "reconstruction_enabled": reconstruction_enabled,
    }
    strict_failures = [key for key, ok in strict_checks.items() if not ok]
    review_failures = [key for key, ok in review_checks.items() if not ok]
    tier2_blocking_failures = [
        key
        for key in review_failures
        if key not in {"min_consecutive_charge_count", "charge_continuity"}
    ]
    if not review_checks["reconstruction_enabled"] or not review_checks["in_neutral_mass_range"] or not review_checks["min_charge_count"]:
        tier = "Tier_4_rejected"
        reason = review_failures[0] if review_failures else "not_ready"
    elif not strict_failures:
        tier = "Tier_1_high_quality"
        reason = "strict_criteria_passed"
    elif not tier2_blocking_failures:
        tier = "Tier_2_supported"
        reason = "review_criteria_passed" if not review_failures else "; ".join(review_failures)
    elif review_checks["min_charge_count"] and review_checks["reconstruction_enabled"] and review_checks["in_neutral_mass_range"]:
        tier = "Tier_3_weak"
        reason = "; ".join(review_failures) or "weak_support"
    else:
        tier = "Tier_4_rejected"
        reason = "; ".join(review_failures) or "not_ready"
    return {
        "Pass_Min_Charge_Count": pass_min_charge,
        "Pass_Min_Consecutive_Charge_Count": pass_min_consecutive,
        "Pass_Charge_Continuity": pass_continuity,
        "Pass_Internal_Error": pass_internal,
        "Pass_Neutral_Mass_SD": pass_sd,
        "Pass_Neutral_Mass_Range": pass_range,
        "Pass_RT_Consistency": pass_rt,
        "Pass_Local_Intensity": pass_local,
        "Pass_Competing_Envelope": pass_competing,
        "Pass_Peak_Sharing": pass_peak_sharing,
        "Num_Strict_Criteria_Passed": sum(1 for ok in strict_checks.values() if ok),
        "Num_Review_Criteria_Passed": sum(1 for ok in review_checks.values() if ok),
        "Strict_Failure_Reasons": "; ".join(strict_failures),
        "Review_Failure_Reasons": "; ".join(review_failures),
        "Intact_Quality_Tier": tier,
        "Quality_Tier_Reason": reason,
        "Quality_Tier_Rank": _tier_rank(tier),
    }



def _parse_semicolon_set(value: Any) -> set[str]:
    return {item.strip() for item in str(value or "").split(";") if item.strip()}


def _component_text(items: list[tuple[str, float]]) -> str:
    return "; ".join(f"{name}={value:.6f}" for name, value in items if abs(value) > 1e-12)


def _component_sum(text: str) -> float:
    total = 0.0
    for item in str(text or "").split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        try:
            total += float(item.split("=", 1)[1])
        except ValueError:
            continue
    return total


def _bounded_ratio(value: Any, limit: float, default: float = 0.0) -> float:
    if limit <= 0:
        return default
    return max(0.0, min(_safe_float(value, default) / limit, 1.0))


def _evidence_score(row: dict[str, Any], qc_config: dict[str, Any]) -> tuple[float, str, str]:
    competitive = qc_config["competitive_assignment"]
    weights = competitive["score_weights"]
    charge_count = _safe_float(row.get("Num_Supporting_Charge_States"))
    longest_run = _safe_float(row.get("Longest_Consecutive_Charge_Run"))
    coverage = max(0.0, min(_safe_float(row.get("Charge_Coverage_Fraction")), 1.0))
    local_intensity = max(0.0, min(_safe_float(row.get("Local_Envelope_Relative_Intensity_Percent")) / 100.0, 1.0))
    scan_count = len(_parse_semicolon_set(row.get("Supporting_Scan_IDs"))) or _safe_float(row.get("Supporting_Peak_Count"))
    rt_consistency = {"consistent": 1.0, "review": 0.5}.get(str(row.get("RT_Consistency") or ""), 0.0)
    extension_values = _parse_semicolon_set(row.get("Extended_Charges_Detected")) | _parse_semicolon_set(row.get("Extended_Weak_Charges_Detected"))
    extension_support = min(len(extension_values) / 2.0, 1.0)
    split_support = 1.0 if row.get("Split_Envelope_Merged") is True else 0.0
    components = [
        ("charge_count", weights["charge_count"] * min(charge_count / 5.0, 1.0)),
        ("consecutive_run", weights["consecutive_charge_run"] * min(longest_run / 5.0, 1.0)),
        ("charge_coverage", weights["charge_coverage"] * coverage),
        ("local_intensity", weights["local_relative_intensity"] * local_intensity),
        ("supporting_scan_count", weights["supporting_scan_count"] * min(scan_count / 5.0, 1.0)),
        ("rt_consistency", weights["rt_consistency"] * rt_consistency),
        ("extension_support", weights["extension_support"] * extension_support),
        ("split_merge_support", weights["split_merge_support"] * split_support),
    ]
    severe_count = len(_parse_semicolon_set(row.get("Severe_Limiting_Factors")))
    penalties = [
        ("internal_error", -weights["internal_error_penalty"] * _bounded_ratio(row.get("Envelope_Internal_Error_ppm"), qc_config["max_envelope_internal_error_ppm"])),
        ("neutral_mass_sd", -weights["neutral_mass_sd_penalty"] * _bounded_ratio(row.get("Neutral_Mass_SD"), qc_config["max_neutral_mass_sd_da"])),
        ("neutral_mass_range", -weights["neutral_mass_range_penalty"] * _bounded_ratio(row.get("Neutral_Mass_Range"), qc_config["max_neutral_mass_range_da"])),
        ("rt_range", -weights["rt_range_penalty"] * _bounded_ratio(row.get("RT_Range_Min"), qc_config["max_rt_range_min_for_review"])),
        ("peak_sharing", -weights["peak_sharing_penalty"] * max(0.0, min(_safe_float(row.get("Highly_Shared_Peak_Fraction")), 1.0))),
        ("peak_usage", -weights["peak_usage_penalty"] * min(max(_safe_float(row.get("Max_Peak_Usage_Count")) - 1.0, 0.0) / 10.0, 1.0)),
        ("charge_gap", -weights["charge_gap_penalty"] * min(_safe_float(row.get("Charge_Gap_Count")) / 3.0, 1.0)),
        ("severe_limiting_factor", -weights["severe_limiting_factor_penalty"] * min(severe_count / 3.0, 1.0)),
    ]
    score = round(sum(value for _, value in components) + sum(value for _, value in penalties), 6)
    return score, _component_text(components), _component_text(penalties)


def _shared_peak_fraction(left: set[str], right: set[str]) -> float:
    denominator = min(len(left), len(right))
    if denominator <= 0:
        return 0.0
    return len(left & right) / denominator


def apply_competitive_assignment(qc_rows: list[dict[str, Any]], qc_config: dict[str, Any]) -> dict[str, Any]:
    started = perf_counter()
    competitive = qc_config["competitive_assignment"]
    stats: dict[str, Any] = {"Competitive_Scoring_Time_Seconds": 0.0}
    if not qc_rows:
        stats.update({
            "Competing_Envelope_Group_Count": 0,
            "MultiCandidate_Competition_Group_Count": 0,
            "Noncompeting_Candidate_Count": 0,
            "Largest_Competition_Group_Size": 0,
            "Median_Competition_Group_Size": None,
            "Mean_Competition_Group_Size": None,
            "Competition_Group_Size_Distribution": "",
            "Median_Top_Evidence_Score": None,
            "Median_Top_vs_Second_Score_Margin": None,
            "Competition_Groups_With_Close_Scores": 0,
            "High_Peak_Sharing_Candidate_Count": 0,
        })
        return stats

    for row in qc_rows:
        score, components, penalties = _evidence_score(row, qc_config)
        row["Envelope_Evidence_Score"] = score
        row["Evidence_Score_Components"] = components
        row["Evidence_Score_Penalties"] = penalties
        row["Evidence_Score_Config_Version"] = competitive["evidence_score_config_version"]
        row.setdefault("Competing_Envelope_Group_ID", "")
        row.setdefault("Competing_Envelope_Group_Size", 1)
        row.setdefault("Shared_Peak_Competitor_Count", 0)
        row.setdefault("Maximum_Shared_Peak_Fraction", 0.0)
        row.setdefault("Mean_Shared_Peak_Fraction", 0.0)
        row.setdefault("Competitor_Cluster_IDs", "")
        row.setdefault("Is_Noncompeting_Candidate", True)
        row.setdefault("Evidence_Score_Rank_In_Competition", 1)

    if not competitive.get("enabled", True):
        for index, row in enumerate(sorted(qc_rows, key=lambda item: str(item.get("Cluster_ID") or "")), start=1):
            row["Competing_Envelope_Group_ID"] = f"CG{index:05d}"
            row["Competing_Envelope_Group_Size"] = 1
            row["Is_Noncompeting_Candidate"] = True
            row["Evidence_Score_Rank_In_Competition"] = 1
        stats.update(_competition_stats(qc_rows, competitive, perf_counter() - started))
        return stats

    eligible_rows = [row for row in qc_rows if row.get("In_Neutral_Mass_Search_Range")]
    ids = [str(row.get("Cluster_ID") or index) for index, row in enumerate(eligible_rows)]
    dsu = _DisjointSet(ids)
    by_peak: dict[str, list[dict[str, Any]]] = {}
    for row in eligible_rows:
        for peak_id in row.get("_supporting_local_peak_id_set") or set():
            by_peak.setdefault(str(peak_id), []).append(row)
    rt_tolerance = competitive["rt_tolerance_min"]
    for candidates in by_peak.values():
        ordered = sorted(candidates, key=lambda item: (_safe_float(item.get("RT_Mean")), str(item.get("Cluster_ID") or "")))
        for index, row in enumerate(ordered):
            row_rt = _safe_float(row.get("RT_Mean"), None)
            for other in ordered[index + 1:]:
                other_rt = _safe_float(other.get("RT_Mean"), None)
                if row_rt is not None and other_rt is not None and other_rt - row_rt > rt_tolerance:
                    break
                if row_rt is None or other_rt is None or abs(other_rt - row_rt) <= rt_tolerance:
                    dsu.union(str(row.get("Cluster_ID")), str(other.get("Cluster_ID")))
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in eligible_rows:
        grouped.setdefault(dsu.find(str(row.get("Cluster_ID"))), []).append(row)
    group_members = list(grouped.values()) + [[row] for row in qc_rows if not row.get("In_Neutral_Mass_Search_Range")]
    group_members.sort(key=lambda members: min(str(row.get("Cluster_ID") or "") for row in members))
    for group_index, members in enumerate(group_members, start=1):
        group_id = f"CG{group_index:05d}"
        member_ids = [str(row.get("Cluster_ID") or "") for row in members]
        shared_fractions: dict[str, list[float]] = {member_id: [] for member_id in member_ids}
        competitor_ids: dict[str, set[str]] = {member_id: set() for member_id in member_ids}
        for left_index, left in enumerate(members):
            left_id = str(left.get("Cluster_ID") or "")
            left_peaks = set(left.get("_supporting_local_peak_id_set") or set())
            for right in members[left_index + 1:]:
                right_id = str(right.get("Cluster_ID") or "")
                right_peaks = set(right.get("_supporting_local_peak_id_set") or set())
                fraction = _shared_peak_fraction(left_peaks, right_peaks)
                if fraction > 0.0:
                    shared_fractions[left_id].append(fraction)
                    shared_fractions[right_id].append(fraction)
                    competitor_ids[left_id].add(right_id)
                    competitor_ids[right_id].add(left_id)
        ranked = sorted(members, key=_competition_rank_key)
        for rank, row in enumerate(ranked, start=1):
            row_id = str(row.get("Cluster_ID") or "")
            fractions = shared_fractions.get(row_id, [])
            row["Competing_Envelope_Group_ID"] = group_id
            row["Competing_Envelope_Group_Size"] = len(members)
            row["Shared_Peak_Competitor_Count"] = len(competitor_ids.get(row_id, set()))
            row["Maximum_Shared_Peak_Fraction"] = max(fractions, default=0.0)
            row["Mean_Shared_Peak_Fraction"] = (sum(fractions) / len(fractions)) if fractions else 0.0
            row["Competitor_Cluster_IDs"] = "; ".join(sorted(competitor_ids.get(row_id, set())))
            row["Is_Noncompeting_Candidate"] = len(members) == 1
            row["Evidence_Score_Rank_In_Competition"] = rank
    stats.update(_competition_stats(qc_rows, competitive, perf_counter() - started))
    return stats


def _competition_rank_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -_safe_float(row.get("Envelope_Evidence_Score")),
        _safe_float(row.get("Quality_Tier_Rank"), 4.0),
        -_safe_float(row.get("Num_Supporting_Charge_States")),
        _safe_float(row.get("Envelope_Internal_Error_ppm"), 1_000_000.0),
        _safe_float(row.get("Reconstructed_Mass")),
        str(row.get("Cluster_ID") or ""),
    )


def _competition_stats(qc_rows: list[dict[str, Any]], competitive: dict[str, Any], elapsed: float) -> dict[str, Any]:
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in qc_rows:
        group_id = str(row.get("Competing_Envelope_Group_ID") or "")
        if group_id:
            groups.setdefault(group_id, []).append(row)
    sizes = [len(members) for members in groups.values()]
    top_scores = []
    margins = []
    close = 0
    for members in groups.values():
        ranked = sorted(members, key=_competition_rank_key)
        if ranked:
            top_scores.append(_safe_float(ranked[0].get("Envelope_Evidence_Score")))
        if len(ranked) > 1:
            margin = _safe_float(ranked[0].get("Envelope_Evidence_Score")) - _safe_float(ranked[1].get("Envelope_Evidence_Score"))
            margins.append(margin)
            if margin <= competitive["close_score_margin"]:
                close += 1
    return {
        "Competing_Envelope_Group_Count": len(groups),
        "MultiCandidate_Competition_Group_Count": sum(1 for size in sizes if size > 1),
        "Noncompeting_Candidate_Count": sum(1 for row in qc_rows if row.get("Is_Noncompeting_Candidate")),
        "Largest_Competition_Group_Size": max(sizes, default=0),
        "Median_Competition_Group_Size": median(sizes) if sizes else None,
        "Mean_Competition_Group_Size": (sum(sizes) / len(sizes)) if sizes else None,
        "Competition_Group_Size_Distribution": _distribution_summary(sizes),
        "Median_Top_Evidence_Score": median(top_scores) if top_scores else None,
        "Median_Top_vs_Second_Score_Margin": median(margins) if margins else None,
        "Competition_Groups_With_Close_Scores": close,
        "High_Peak_Sharing_Candidate_Count": sum(1 for row in qc_rows if str(row.get("Peak_Sharing_Status") or "") == "highly_shared"),
        "Competitive_Scoring_Time_Seconds": elapsed,
    }




def _peak_charge_keys(row: dict[str, Any]) -> set[str]:
    """Return observed peak/assumed-charge pairs without inventing combinations."""
    peak_charge_map = row.get("_supporting_peak_charge_map") or {}
    return {
        f"{peak}|z={charge}"
        for peak, charges in peak_charge_map.items()
        for charge in charges
    }


def _charges_for_peaks(row: dict[str, Any], peak_ids: set[str]) -> set[int]:
    peak_charge_map = row.get("_supporting_peak_charge_map") or {}
    if peak_charge_map:
        return {
            int(charge)
            for peak_id in peak_ids
            for charge in peak_charge_map.get(str(peak_id), set())
        }
    # Legacy/synthetic rows may not retain peak-level assumed charges. In that
    # case avoid a fabricated Cartesian product and do not make the charge gate
    # stricter than the observed-peak gate.
    return set(row.get("_supporting_charge_set") or set()) if peak_ids else set()


def _initialize_assignment_fields(row: dict[str, Any]) -> None:
    row.setdefault("Direct_Competitor_Count", 0)
    row.setdefault("Direct_Competitor_Cluster_IDs", "")
    row.setdefault("Direct_Shared_Peak_Count_Max", 0)
    row.setdefault("Direct_Shared_Peak_Fraction_Max", 0.0)
    row.setdefault("Competition_Component_Size", row.get("Competing_Envelope_Group_Size", 1))
    row.setdefault("Dry_Run_Assignment_Status", "not_evaluated")
    row.setdefault("Dry_Run_Selected", False)
    row.setdefault("Dry_Run_Selection_Order", None)
    row.setdefault("Supporting_Peak_Count_Before_Assignment", len(row.get("_supporting_local_peak_id_set") or set()))
    row.setdefault("Independent_Supporting_Peak_Count", 0)
    row.setdefault("Independent_Supporting_Peak_Fraction", 0.0)
    row.setdefault("Supporting_Charge_Count_Before_Assignment", len(row.get("_supporting_charge_set") or set()))
    row.setdefault("Independent_Charge_State_Count", 0)
    row.setdefault("Peaks_Already_Assigned_Count", 0)
    row.setdefault("Charges_Already_Assigned_Count", 0)
    row.setdefault("Excluded_By_Cluster_ID", "")
    row.setdefault("Dry_Run_Exclusion_Reason", "")
    row.setdefault("Score_Margin_To_Excluding_Candidate", None)
    row.setdefault("Close_Score_Ambiguity", False)
    row.setdefault("Assignment_Confidence", "low")
    row.setdefault("Shared_Observed_Peak_Count", 0)
    row.setdefault("Shared_Peak_Charge_Assignment_Count", 0)
    row.setdefault("Independent_Observed_Peak_Count", 0)


def apply_assignment_dry_run(qc_rows: list[dict[str, Any]], qc_config: dict[str, Any]) -> dict[str, Any]:
    started = perf_counter()
    competitive = qc_config["competitive_assignment"]
    for row in qc_rows:
        _initialize_assignment_fields(row)
    if not qc_rows:
        return _assignment_stats(qc_rows, elapsed=perf_counter() - started)
    if not competitive.get("enabled", True) or not competitive.get("dry_run", True):
        for row in qc_rows:
            row["Dry_Run_Assignment_Status"] = "not_evaluated"
            row["Assignment_Confidence"] = "low"
        return _assignment_stats(qc_rows, elapsed=perf_counter() - started)

    groups: dict[str, list[dict[str, Any]]] = {}
    for row in qc_rows:
        groups.setdefault(str(row.get("Competing_Envelope_Group_ID") or row.get("Cluster_ID") or ""), []).append(row)

    for group_id in sorted(groups):
        selection_order = 1
        members = groups[group_id]
        direct: dict[str, set[str]] = {str(row.get("Cluster_ID") or ""): set() for row in members}
        shared_peak_count_max: dict[str, int] = {str(row.get("Cluster_ID") or ""): 0 for row in members}
        shared_fraction_max: dict[str, float] = {str(row.get("Cluster_ID") or ""): 0.0 for row in members}
        shared_assignment_max: dict[str, int] = {str(row.get("Cluster_ID") or ""): 0 for row in members}
        by_peak: dict[str, list[dict[str, Any]]] = {}
        for row in members:
            for peak_id in row.get("_supporting_local_peak_id_set") or set():
                by_peak.setdefault(str(peak_id), []).append(row)
        direct_pairs: set[tuple[str, str]] = set()
        for peak_members in by_peak.values():
            ordered = sorted(peak_members, key=lambda item: str(item.get("Cluster_ID") or ""))
            for left_index, left in enumerate(ordered):
                left_id = str(left.get("Cluster_ID") or "")
                left_peaks = set(left.get("_supporting_local_peak_id_set") or set())
                left_assignments = _peak_charge_keys(left)
                for right in ordered[left_index + 1:]:
                    right_id = str(right.get("Cluster_ID") or "")
                    right_peaks = set(right.get("_supporting_local_peak_id_set") or set())
                    shared_count = len(left_peaks & right_peaks)
                    if shared_count <= 0:
                        continue
                    direct[left_id].add(right_id)
                    direct[right_id].add(left_id)
                    direct_pairs.add(tuple(sorted((left_id, right_id))))
                    fraction = _shared_peak_fraction(left_peaks, right_peaks)
                    shared_peak_count_max[left_id] = max(shared_peak_count_max[left_id], shared_count)
                    shared_peak_count_max[right_id] = max(shared_peak_count_max[right_id], shared_count)
                    shared_fraction_max[left_id] = max(shared_fraction_max[left_id], fraction)
                    shared_fraction_max[right_id] = max(shared_fraction_max[right_id], fraction)
                    shared_assignments = len(left_assignments & _peak_charge_keys(right))
                    shared_assignment_max[left_id] = max(shared_assignment_max[left_id], shared_assignments)
                    shared_assignment_max[right_id] = max(shared_assignment_max[right_id], shared_assignments)
        for row in members:
            row_id = str(row.get("Cluster_ID") or "")
            row["Direct_Competitor_Count"] = len(direct[row_id])
            row["Direct_Competitor_Cluster_IDs"] = "; ".join(sorted(direct[row_id]))
            row["Direct_Shared_Peak_Count_Max"] = shared_peak_count_max[row_id]
            row["Direct_Shared_Peak_Fraction_Max"] = shared_fraction_max[row_id]
            row["Competition_Component_Size"] = len(members)
            row["Shared_Observed_Peak_Count"] = shared_peak_count_max[row_id]
            row["Shared_Peak_Charge_Assignment_Count"] = shared_assignment_max[row_id]
            row["Supporting_Peak_Count_Before_Assignment"] = len(row.get("_supporting_local_peak_id_set") or set())
            row["Supporting_Charge_Count_Before_Assignment"] = len(row.get("_supporting_charge_set") or set())
        if len(members) == 1:
            row = members[0]
            row["Dry_Run_Assignment_Status"] = "noncompeting"
            row["Dry_Run_Selected"] = True
            row["Dry_Run_Selection_Order"] = selection_order
            row["Independent_Supporting_Peak_Count"] = len(row.get("_supporting_local_peak_id_set") or set())
            row["Independent_Observed_Peak_Count"] = row["Independent_Supporting_Peak_Count"]
            row["Independent_Supporting_Peak_Fraction"] = 1.0 if row["Independent_Supporting_Peak_Count"] else 0.0
            row["Independent_Charge_State_Count"] = len(row.get("_supporting_charge_set") or set())
            row["Assignment_Confidence"] = "high"
            selection_order += 1
            continue
        assigned_peaks: set[str] = set()
        assigned_charges: set[int] = set()
        selected_rows: list[dict[str, Any]] = []
        for row in sorted(members, key=_competition_rank_key):
            row_id = str(row.get("Cluster_ID") or "")
            peaks = set(row.get("_supporting_local_peak_id_set") or set())
            charges = set(row.get("_supporting_charge_set") or set())
            direct_selected = [selected for selected in selected_rows if str(selected.get("Cluster_ID") or "") in direct[row_id]]
            assigned_direct_peaks = set().union(*(set(selected.get("_supporting_local_peak_id_set") or set()) for selected in direct_selected)) if direct_selected else set()
            assigned_direct_charges = set().union(*(set(selected.get("_supporting_charge_set") or set()) for selected in direct_selected)) if direct_selected else set()
            allow_shared = competitive.get("allow_shared_peaks_between_selected", False)
            reusable_peaks = assigned_direct_peaks if not allow_shared else set()
            independent_peaks = peaks - reusable_peaks
            independent_charges = _charges_for_peaks(row, independent_peaks)
            if allow_shared:
                independent_charges = charges
            before_count = len(peaks)
            independent_fraction = len(independent_peaks) / before_count if before_count else 0.0
            row["Peaks_Already_Assigned_Count"] = len(peaks & assigned_direct_peaks)
            row["Charges_Already_Assigned_Count"] = len(charges & assigned_direct_charges)
            row["Independent_Supporting_Peak_Count"] = len(independent_peaks)
            row["Independent_Observed_Peak_Count"] = len(independent_peaks)
            row["Independent_Supporting_Peak_Fraction"] = independent_fraction
            row["Independent_Charge_State_Count"] = len(independent_charges)
            excluder = max(direct_selected, key=lambda item: _safe_float(item.get("Envelope_Evidence_Score")), default=None)
            score_margin = (_safe_float(excluder.get("Envelope_Evidence_Score")) - _safe_float(row.get("Envelope_Evidence_Score"))) if excluder else None
            row["Excluded_By_Cluster_ID"] = excluder.get("Cluster_ID") if excluder else ""
            row["Score_Margin_To_Excluding_Candidate"] = score_margin
            close_score = score_margin is not None and score_margin < competitive["minimum_score_margin_for_exclusive_selection"]
            row["Close_Score_Ambiguity"] = close_score
            if not direct_selected:
                row["Dry_Run_Assignment_Status"] = "selected_primary"
                row["Dry_Run_Selected"] = True
                row["Dry_Run_Selection_Order"] = selection_order
                row["Assignment_Confidence"] = "high"
                selected_rows.append(row)
                assigned_peaks.update(peaks)
                assigned_charges.update(charges)
                selection_order += 1
                continue
            peak_ok = independent_fraction >= competitive["min_independent_peak_fraction"]
            charge_ok = len(independent_charges) >= competitive["min_independent_charge_states"]
            if peak_ok and charge_ok:
                row["Dry_Run_Assignment_Status"] = "selected_independent"
                row["Dry_Run_Selected"] = True
                row["Dry_Run_Selection_Order"] = selection_order
                row["Assignment_Confidence"] = "ambiguous" if close_score else "moderate"
                selected_rows.append(row)
                assigned_peaks.update(independent_peaks)
                assigned_charges.update(independent_charges)
                selection_order += 1
            elif close_score:
                row["Dry_Run_Assignment_Status"] = "would_exclude_peak_reuse" if not peak_ok else "would_exclude_independent_charge_shortage"
                row["Dry_Run_Selected"] = False
                row["Dry_Run_Exclusion_Reason"] = "ambiguous_close_score"
                row["Assignment_Confidence"] = "ambiguous"
            elif not peak_ok:
                row["Dry_Run_Assignment_Status"] = "would_exclude_independent_peak_shortage" if len(independent_peaks) else "would_exclude_peak_reuse"
                row["Dry_Run_Selected"] = False
                row["Dry_Run_Exclusion_Reason"] = "insufficient_independent_peak_support"
                row["Assignment_Confidence"] = "high" if score_margin is not None and score_margin >= competitive["minimum_score_margin_for_exclusive_selection"] else "low"
            else:
                row["Dry_Run_Assignment_Status"] = "would_exclude_independent_charge_shortage"
                row["Dry_Run_Selected"] = False
                row["Dry_Run_Exclusion_Reason"] = "insufficient_independent_charge_support"
                row["Assignment_Confidence"] = "high" if score_margin is not None and score_margin >= competitive["minimum_score_margin_for_exclusive_selection"] else "low"
    return _assignment_stats(qc_rows, elapsed=perf_counter() - started)


def _assignment_stats(qc_rows: list[dict[str, Any]], elapsed: float) -> dict[str, Any]:
    selected = [row for row in qc_rows if row.get("Dry_Run_Selected")]
    would_exclude = [row for row in qc_rows if str(row.get("Dry_Run_Assignment_Status") or "").startswith("would_exclude")]
    ambiguous = [row for row in qc_rows if row.get("Close_Score_Ambiguity") or row.get("Assignment_Confidence") == "ambiguous"]
    fractions = [_safe_float(row.get("Independent_Supporting_Peak_Fraction")) for row in qc_rows if row.get("Independent_Supporting_Peak_Fraction") is not None]
    direct_counts = [_safe_float(row.get("Direct_Competitor_Count")) for row in qc_rows if row.get("Direct_Competitor_Count") is not None]
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in qc_rows:
        groups.setdefault(str(row.get("Competing_Envelope_Group_ID") or ""), []).append(row)
    largest = max(groups.values(), key=len, default=[])
    return {
        "Dry_Run_Selected_Candidate_Count": len(selected),
        "Dry_Run_Primary_Selected_Count": sum(1 for row in qc_rows if row.get("Dry_Run_Assignment_Status") == "selected_primary"),
        "Dry_Run_Independent_Selected_Count": sum(1 for row in qc_rows if row.get("Dry_Run_Assignment_Status") == "selected_independent"),
        "Dry_Run_Would_Exclude_Count": len(would_exclude),
        "Dry_Run_Ambiguous_Count": len(ambiguous),
        "Dry_Run_Noncompeting_Count": sum(1 for row in qc_rows if row.get("Dry_Run_Assignment_Status") == "noncompeting"),
        "Would_Exclude_Peak_Reuse_Count": sum(1 for row in qc_rows if row.get("Dry_Run_Assignment_Status") == "would_exclude_peak_reuse"),
        "Would_Exclude_Independent_Peak_Shortage_Count": sum(1 for row in qc_rows if row.get("Dry_Run_Assignment_Status") == "would_exclude_independent_peak_shortage"),
        "Would_Exclude_Independent_Charge_Shortage_Count": sum(1 for row in qc_rows if row.get("Dry_Run_Assignment_Status") == "would_exclude_independent_charge_shortage"),
        "Median_Independent_Peak_Fraction": median(fractions) if fractions else None,
        "Median_Direct_Competitor_Count": median(direct_counts) if direct_counts else None,
        "Largest_Component_Selected_Count": sum(1 for row in largest if row.get("Dry_Run_Selected")),
        "Largest_Component_Excluded_Count": sum(1 for row in largest if str(row.get("Dry_Run_Assignment_Status") or "").startswith("would_exclude")),
        "Largest_Component_Ambiguous_Count": sum(1 for row in largest if row.get("Close_Score_Ambiguity") or row.get("Assignment_Confidence") == "ambiguous"),
        "Assignment_Dry_Run_Time_Seconds": elapsed,
    }


def build_assignment_dry_run_rows(qc_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [_candidate_projection(row, ASSIGNMENT_DRY_RUN_COLUMNS) for row in sorted(qc_rows, key=lambda item: (str(item.get("Competing_Envelope_Group_ID") or ""), _safe_float(item.get("Dry_Run_Selection_Order"), 1_000_000), str(item.get("Cluster_ID") or "")))]


def build_assignment_dry_run_summary_rows(qc_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in qc_rows:
        groups.setdefault(str(row.get("Competing_Envelope_Group_ID") or ""), []).append(row)
    rows = []
    for group_id in sorted(groups):
        members = groups[group_id]
        local_peaks = set().union(*(set(row.get("_supporting_local_peak_id_set") or set()) for row in members)) if members else set()
        reasons = _reason_summary([row.get("Dry_Run_Exclusion_Reason") for row in members])
        rows.append({
            "Competing_Envelope_Group_ID": group_id,
            "Component_Size": len(members),
            "Unique_Local_Peak_Count": len(local_peaks),
            "Direct_Competition_Edge_Count": int(sum(_safe_float(row.get("Direct_Competitor_Count")) for row in members) / 2),
            "Selected_Primary_Count": sum(1 for row in members if row.get("Dry_Run_Assignment_Status") == "selected_primary"),
            "Selected_Independent_Count": sum(1 for row in members if row.get("Dry_Run_Assignment_Status") == "selected_independent"),
            "Noncompeting_Count": sum(1 for row in members if row.get("Dry_Run_Assignment_Status") == "noncompeting"),
            "Would_Exclude_Count": sum(1 for row in members if str(row.get("Dry_Run_Assignment_Status") or "").startswith("would_exclude")),
            "Ambiguous_Count": sum(1 for row in members if row.get("Close_Score_Ambiguity") or row.get("Assignment_Confidence") == "ambiguous"),
            "Candidates_Remaining_After_Dry_Run": sum(1 for row in members if row.get("Dry_Run_Selected")),
            "Primary_Selected_Cluster_IDs": "; ".join(str(row.get("Cluster_ID") or "") for row in members if row.get("Dry_Run_Assignment_Status") == "selected_primary"),
            "Main_Exclusion_Reasons": reasons,
            "Notes": "diagnostic_only_no_candidate_exclusion",
        })
    return rows

def build_intact_competition_group_rows(qc_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in qc_rows:
        group_id = str(row.get("Competing_Envelope_Group_ID") or "")
        if group_id:
            groups.setdefault(group_id, []).append(row)
    rows = []
    for group_id in sorted(groups):
        members = groups[group_id]
        ranked = sorted(members, key=_competition_rank_key)
        top = ranked[0] if ranked else {}
        second = ranked[1] if len(ranked) > 1 else {}
        local_peaks = set()
        for member in members:
            local_peaks.update(member.get("_supporting_local_peak_id_set") or set())
        top_score = top.get("Envelope_Evidence_Score")
        second_score = second.get("Envelope_Evidence_Score")
        rows.append({
            "Competing_Envelope_Group_ID": group_id,
            "Group_Size": len(members),
            "MultiCandidate_Group": len(members) > 1,
            "Top_Ranked_Cluster_ID": top.get("Cluster_ID"),
            "Top_Ranked_Mass": top.get("Reconstructed_Mass"),
            "Top_Evidence_Score": top_score,
            "Second_Evidence_Score": second_score,
            "Score_Margin_Top_vs_Second": (_safe_float(top_score) - _safe_float(second_score)) if second else None,
            "Shared_Local_Peak_Count": len(local_peaks),
            "Maximum_Peak_Usage_Count": max((_safe_float(member.get("Max_Peak_Usage_Count")) for member in members), default=0.0),
            "Member_Cluster_IDs": "; ".join(str(member.get("Cluster_ID") or "") for member in ranked),
            "Member_Masses": "; ".join(_format_float(member.get("Reconstructed_Mass"), 6) for member in ranked),
            "Member_Quality_Tiers": "; ".join(str(member.get("Intact_Quality_Tier") or "") for member in ranked),
            "Notes": "diagnostic_only_no_candidate_exclusion",
        })
    return rows


def build_intact_competition_score_rows(qc_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for row in sorted(qc_rows, key=lambda item: (str(item.get("Competing_Envelope_Group_ID") or ""), _safe_float(item.get("Evidence_Score_Rank_In_Competition")), str(item.get("Cluster_ID") or ""))):
        rows.append({
            "Cluster_ID": row.get("Cluster_ID"),
            "Reconstructed_Mass": row.get("Reconstructed_Mass"),
            "Competing_Envelope_Group_ID": row.get("Competing_Envelope_Group_ID"),
            "Competing_Envelope_Group_Size": row.get("Competing_Envelope_Group_Size"),
            "Envelope_Evidence_Score": row.get("Envelope_Evidence_Score"),
            "Evidence_Score_Rank_In_Competition": row.get("Evidence_Score_Rank_In_Competition"),
            "Evidence_Score_Components": row.get("Evidence_Score_Components"),
            "Evidence_Score_Penalties": row.get("Evidence_Score_Penalties"),
            "Intact_Quality_Tier": row.get("Intact_Quality_Tier"),
            "Supporting_Charge_Count": row.get("Num_Supporting_Charge_States"),
            "Longest_Consecutive_Charge_Run": row.get("Longest_Consecutive_Charge_Run"),
            "Charge_Coverage_Fraction": row.get("Charge_Coverage_Fraction"),
            "Local_Envelope_Relative_Intensity_Percent": row.get("Local_Envelope_Relative_Intensity_Percent"),
            "Envelope_Internal_Error_ppm": row.get("Envelope_Internal_Error_ppm"),
            "Neutral_Mass_SD": row.get("Neutral_Mass_SD"),
            "Neutral_Mass_Range": row.get("Neutral_Mass_Range"),
            "RT_Range_Min": row.get("RT_Range_Min"),
            "Maximum_Shared_Peak_Fraction": row.get("Maximum_Shared_Peak_Fraction"),
            "Max_Peak_Usage_Count": row.get("Max_Peak_Usage_Count"),
            "Direct_Competitor_Count": row.get("Direct_Competitor_Count"),
            "Dry_Run_Assignment_Status": row.get("Dry_Run_Assignment_Status"),
            "Dry_Run_Selected": row.get("Dry_Run_Selected"),
            "Independent_Supporting_Peak_Count": row.get("Independent_Supporting_Peak_Count"),
            "Independent_Supporting_Peak_Fraction": row.get("Independent_Supporting_Peak_Fraction"),
            "Independent_Charge_State_Count": row.get("Independent_Charge_State_Count"),
            "Assignment_Confidence": row.get("Assignment_Confidence"),
            "Limiting_Factors": row.get("Limiting_Factors"),
        })
    return rows

def _intact_qc_score(row: dict[str, Any], qc_config: dict[str, Any]) -> float:
    score = 0.0
    if row.get("Intact_Strict_Eligible"):
        score += 30.0
    elif row.get("Intact_Review_Eligible"):
        score += 22.0
    elif row.get("Envelope_QC_Eligible"):
        score += 12.0
    if row.get("Charge_State_Continuity") == "contiguous":
        score += 12.0
    score += min(_safe_float(row.get("Num_Supporting_Charge_States")), 6.0) * 4.0
    score += {"consistent": 12.0, "review": 6.0}.get(str(row.get("RT_Consistency") or ""), 0.0)
    score += _small_metric_score(row.get("Envelope_Internal_Error_ppm"), qc_config["max_envelope_internal_error_ppm"], 10.0)
    score += _small_metric_score(row.get("Neutral_Mass_SD"), qc_config["max_neutral_mass_sd_da"], 8.0)
    score += _small_metric_score(row.get("Neutral_Mass_Range"), qc_config["max_neutral_mass_range_da"], 8.0)
    score += min(_safe_float(row.get("Relative_Overall_Envelope_Intensity_Percent")), 100.0) / 100.0 * 8.0
    return round(score, 3)


def _dominant_sort_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        bool(row.get("Intact_Strict_Eligible")),
        bool(row.get("Intact_Review_Eligible")),
        row.get("Charge_State_Continuity") == "contiguous",
        _safe_float(row.get("Num_Supporting_Charge_States")),
        _rt_rank_score(str(row.get("RT_Consistency") or "")),
        -_safe_float(row.get("Envelope_Internal_Error_ppm"), 1_000_000.0),
        -_safe_float(row.get("Neutral_Mass_SD"), 1_000_000.0),
        -_safe_float(row.get("Neutral_Mass_Range"), 1_000_000.0),
        _safe_float(row.get("Total_Supporting_Intensity")),
    )


def _dominant_row(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return max(rows, key=_dominant_sort_key, default={})


def _reference_match(observed_mass: float, qc_config: dict[str, Any]) -> tuple[str, float | None, float | None, float | None, bool]:
    references = qc_config.get("reference_masses") or []
    if not references:
        return "not_configured", None, None, None, False
    best = None
    for reference in references:
        mass = float(reference["mass_da"])
        error_da = observed_mass - mass
        error_ppm = _ppm(error_da, mass)
        score = abs(error_ppm) if error_ppm is not None else float("inf")
        if best is None or score < best[0]:
            best = (score, str(reference["label"]), mass, error_da, error_ppm)
    if best is None:
        return "not_configured", None, None, None, False
    matched = best[4] is not None and abs(best[4]) <= qc_config["reference_mass_tolerance_ppm"]
    return best[1], best[2], best[3], best[4], matched


def _class_summary(cluster_peaks: list[dict[str, Any]]) -> tuple[str, bool]:
    classes = []
    for row in cluster_peaks:
        value = str(row.get("Peak_Tier") or "").strip()
        if value and value not in classes:
            classes.append(value)
    trace_only = bool(classes) and all(value.lower() == "trace" for value in classes)
    return "; ".join(classes), trace_only



def _format_float(value: Any, digits: int) -> str:
    try:
        if value is None or value == "":
            return "na"
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "na"


def _supporting_peak_id(row: dict[str, Any]) -> str:
    mz = _format_float(row.get("mz"), 6)
    if row.get("Scan_ID") not in {None, ""}:
        scan = str(row.get("Scan_ID"))
    elif mz != "na":
        scan = "scan_na"
    else:
        scan = str(row.get("Cluster_ID") or "scan_na")
    rt = _format_float(row.get("RT"), 5)
    charge = str(row.get("Charge") or "z_na")
    return f"{scan}|rt={rt}|mz={mz}|z={charge}"


def _shared_counts(row: dict[str, Any], representative: dict[str, Any]) -> tuple[int, float, int, float]:
    peaks = set(row.get("_supporting_peak_id_set") or set())
    rep_peaks = set(representative.get("_supporting_peak_id_set") or set())
    charges = set(row.get("_supporting_charge_set") or set())
    rep_charges = set(representative.get("_supporting_charge_set") or set())
    shared_peak_count = len(peaks & rep_peaks)
    shared_charge_count = len(charges & rep_charges)
    peak_denominator = max(min(len(peaks), len(rep_peaks)), 1)
    charge_denominator = max(min(len(charges), len(rep_charges)), 1)
    return (
        shared_peak_count,
        shared_peak_count / peak_denominator,
        shared_charge_count,
        shared_charge_count / charge_denominator,
    )


def _representative_sort_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        bool(row.get("Intact_Strict_Eligible")),
        bool(row.get("Intact_Review_Eligible")),
        bool(row.get("Comparison_Ready_Strict")),
        bool(row.get("Comparison_Ready_Review")),
        row.get("Charge_State_Continuity") == "contiguous",
        _safe_float(row.get("Num_Supporting_Charge_States")),
        _rt_rank_score(str(row.get("RT_Consistency") or "")),
        -_safe_float(row.get("Envelope_Internal_Error_ppm"), 1_000_000.0),
        -_safe_float(row.get("Neutral_Mass_SD"), 1_000_000.0),
        -_safe_float(row.get("Neutral_Mass_Range"), 1_000_000.0),
        _safe_float(row.get("Intact_Envelope_QC_Score")),
        _safe_float(row.get("Total_Supporting_Intensity")),
        str(row.get("Cluster_ID") or ""),
    )


class _DisjointSet:
    def __init__(self, items: list[str]) -> None:
        self.parent = {item: item for item in items}

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def _rt_delta(left: dict[str, Any], right: dict[str, Any]) -> float:
    return abs(_safe_float(left.get("RT_Mean")) - _safe_float(right.get("RT_Mean")))


def _rows_overlap(left: dict[str, Any], right: dict[str, Any], grouping_config: dict[str, Any]) -> bool:
    mass_delta = abs(_safe_float(left.get("Reconstructed_Mass")) - _safe_float(right.get("Reconstructed_Mass")))
    if mass_delta > grouping_config["mass_tolerance_da"]:
        return False
    if _rt_delta(left, right) > grouping_config["rt_tolerance_min"]:
        return False
    _, peak_fraction, _, charge_fraction = _shared_counts(left, right)
    peak_ok = peak_fraction >= grouping_config["min_shared_peak_fraction"]
    charge_ok = charge_fraction >= grouping_config["min_shared_charge_fraction"]
    if grouping_config.get("require_peak_overlap", True):
        return peak_ok and charge_ok
    return peak_ok or charge_ok


def _comparison_exclusion_reason(row: dict[str, Any]) -> str:
    if not row.get("Group_Representative"):
        return "not_group_representative"
    if not row.get("In_Neutral_Mass_Search_Range"):
        return "outside_neutral_mass_search_range"
    if row.get("Severe_Limiting_Factors"):
        return "severe_limiting_factors"
    if not (row.get("Intact_Strict_Eligible") or row.get("Intact_Review_Eligible")):
        return "not_intact_eligible"
    if not row.get("Comparison_Ready"):
        return "not_comparison_ready"
    return ""


def apply_intact_envelope_grouping(qc_rows: list[dict[str, Any]], qc_config: dict[str, Any]) -> None:
    if not qc_rows:
        return
    exact_groups: dict[str, list[dict[str, Any]]] = {}
    for row in qc_rows:
        key = str(row.get("Exact_Peak_Set_Key") or row.get("Cluster_ID") or "")
        exact_groups.setdefault(key, []).append(row)

    exact_representatives: list[dict[str, Any]] = []
    for index, key in enumerate(sorted(exact_groups), start=1):
        members = exact_groups[key]
        representative = max(members, key=_representative_sort_key)
        group_id = f"ED{index:05d}"
        for member in members:
            member["Exact_Duplicate_Group_ID"] = group_id
            member["Exact_Duplicate_Count"] = len(members)
            member["Is_Exact_Duplicate_Representative"] = member is representative
        exact_representatives.append(representative)

    grouping_config = qc_config["envelope_grouping"]
    if grouping_config.get("enabled", True):
        dsu = _DisjointSet([str(row.get("Cluster_ID")) for row in exact_representatives])
        mass_width = max(grouping_config["mass_tolerance_da"], 0.000001)
        rt_width = max(grouping_config["rt_tolerance_min"], 0.000001)
        buckets: dict[tuple[int, int], list[dict[str, Any]]] = {}
        for row in exact_representatives:
            mass_bucket = int(_safe_float(row.get("Reconstructed_Mass")) // mass_width)
            rt_bucket = int(_safe_float(row.get("RT_Mean")) // rt_width)
            for mass_offset in (-1, 0, 1):
                for rt_offset in (-1, 0, 1):
                    for other in buckets.get((mass_bucket + mass_offset, rt_bucket + rt_offset), []):
                        if _rows_overlap(row, other, grouping_config):
                            dsu.union(str(row.get("Cluster_ID")), str(other.get("Cluster_ID")))
            buckets.setdefault((mass_bucket, rt_bucket), []).append(row)
        grouped_reps: dict[str, list[dict[str, Any]]] = {}
        for row in exact_representatives:
            grouped_reps.setdefault(dsu.find(str(row.get("Cluster_ID"))), []).append(row)
    else:
        grouped_reps = {str(row.get("Cluster_ID")): [row] for row in exact_representatives}

    exact_reps_by_id = {str(row.get("Exact_Duplicate_Group_ID")): row for row in exact_representatives}
    group_index = 1
    for rep_key in sorted(grouped_reps, key=lambda key: min(str(row.get("Cluster_ID")) for row in grouped_reps[key])):
        rep_rows = grouped_reps[rep_key]
        member_rows: list[dict[str, Any]] = []
        for rep in rep_rows:
            exact_id = str(rep.get("Exact_Duplicate_Group_ID"))
            member_rows.extend(exact_groups[str(rep.get("Exact_Peak_Set_Key") or rep.get("Cluster_ID") or "")])
        group_representative = max(member_rows, key=_representative_sort_key)
        group_id = f"IG{group_index:05d}"
        group_index += 1
        exact_count = sum(1 for row in member_rows if row.get("Exact_Duplicate_Count", 1) > 1)
        distinct_exact_groups = {row.get("Exact_Duplicate_Group_ID") for row in member_rows}
        if len(member_rows) == 1:
            ambiguity = "unique"
        elif len(distinct_exact_groups) == 1 and exact_count:
            ambiguity = "exact_duplicates"
        elif any(row.get("Comparison_Ready") for row in member_rows):
            ambiguity = "overlapping_envelopes"
        else:
            ambiguity = "competing_reconstructions"
        for member in member_rows:
            shared_peak_count, shared_peak_fraction, shared_charge_count, shared_charge_fraction = _shared_counts(member, group_representative)
            member["Intact_Envelope_Group_ID"] = group_id
            member["Envelope_Group_Size"] = len(member_rows)
            member["Shared_Peak_Count_With_Representative"] = shared_peak_count
            member["Shared_Peak_Fraction_With_Representative"] = shared_peak_fraction
            member["Shared_Charge_Count_With_Representative"] = shared_charge_count
            member["Shared_Charge_Fraction_With_Representative"] = shared_charge_fraction
            member["Mass_Delta_To_Group_Representative_Da"] = _safe_float(member.get("Reconstructed_Mass")) - _safe_float(group_representative.get("Reconstructed_Mass"))
            member["RT_Delta_To_Group_Representative_Min"] = _safe_float(member.get("RT_Mean")) - _safe_float(group_representative.get("RT_Mean"))
            member["Group_Representative"] = member is group_representative
            member["Group_Ambiguity_Status"] = ambiguity
            reason = _comparison_exclusion_reason(member)
            member["Comparison_Representative"] = reason == ""
            member["Comparison_Representative_Reason"] = "group_representative_intact_ready" if reason == "" else ""
            member["Excluded_From_Comparison_Reason"] = reason
            member["Comparison_Representative_Rank"] = None
            member["Target_Review_Group_Representative"] = False
            member["Target_Review_Rank"] = None
            member["Dominant_Target_Review_Eligible_Flag"] = False

    comparison_rows = sorted([row for row in qc_rows if row.get("Comparison_Representative")], key=_representative_sort_key, reverse=True)
    for rank, row in enumerate(comparison_rows, start=1):
        row["Comparison_Representative_Rank"] = rank
    target_rows = sorted(
        [row for row in comparison_rows if row.get("In_Target_Review_Mass_Range")],
        key=_representative_sort_key,
        reverse=True,
    )
    for rank, row in enumerate(target_rows, start=1):
        row["Target_Review_Group_Representative"] = True
        row["Target_Review_Rank"] = rank
    if target_rows:
        target_rows[0]["Dominant_Target_Review_Eligible_Flag"] = True


def build_intact_envelope_group_rows(qc_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in qc_rows:
        group_id = str(row.get("Intact_Envelope_Group_ID") or "")
        if group_id:
            groups.setdefault(group_id, []).append(row)
    rows = []
    for group_id in sorted(groups):
        members = groups[group_id]
        representative = next((row for row in members if row.get("Group_Representative")), max(members, key=_representative_sort_key))
        masses = [_safe_float(row.get("Reconstructed_Mass")) for row in members]
        rt_values = []
        for row in members:
            if row.get("RT_Min") is not None:
                rt_values.append(_safe_float(row.get("RT_Min")))
            if row.get("RT_Max") is not None:
                rt_values.append(_safe_float(row.get("RT_Max")))
        exact_duplicate_count = sum(1 for row in members if _safe_float(row.get("Exact_Duplicate_Count"), 1) > 1)
        rows.append({
            "Intact_Envelope_Group_ID": group_id,
            "Group_Size": len(members),
            "Exact_Duplicate_Count": exact_duplicate_count,
            "Representative_Cluster_ID": representative.get("Cluster_ID"),
            "Representative_Mass": representative.get("Reconstructed_Mass"),
            "Representative_Status": representative.get("Reconstruction_Status"),
            "Representative_QC_Score": representative.get("Intact_Envelope_QC_Score"),
            "Representative_Comparison_Ready": representative.get("Comparison_Ready"),
            "Representative_Total_Intensity": representative.get("Total_Supporting_Intensity"),
            "Representative_Reconstructed_Envelope_Intensity": representative.get("Reconstructed_Envelope_Intensity"),
            "Intensity_Method": representative.get("Intensity_Method"),
            "Representative_Charge_States": representative.get("Supporting_Charge_States"),
            "Representative_RT_Range": representative.get("RT_Range_Min"),
            "Group_Mass_Min": min(masses) if masses else None,
            "Group_Mass_Max": max(masses) if masses else None,
            "Group_Mass_Range": (max(masses) - min(masses)) if masses else None,
            "Group_RT_Min": min(rt_values) if rt_values else None,
            "Group_RT_Max": max(rt_values) if rt_values else None,
            "Group_Ambiguity_Status": representative.get("Group_Ambiguity_Status"),
            "Member_Cluster_IDs": "; ".join(str(row.get("Cluster_ID")) for row in members),
            "Notes": "",
        })
    return rows


def build_reconstructed_mass_spectrum_rows(
    qc_rows: list[dict[str, Any]],
    reconstruction_config: dict[str, Any],
) -> list[dict[str, Any]]:
    qc_config = _qc_config(reconstruction_config or {})
    output_config = qc_config["mass_spectrum_output"]
    if not output_config.get("enabled", True):
        return []
    rows = []
    for row in qc_rows:
        if not row.get("In_Neutral_Mass_Search_Range"):
            continue
        if output_config.get("representatives_only", True) and not row.get("Group_Representative"):
            continue
        if output_config.get("comparison_ready_only", False) and not row.get("Comparison_Ready"):
            continue
        if not output_config.get("include_qc_ineligible", True) and not (row.get("Intact_Strict_Eligible") or row.get("Intact_Review_Eligible")):
            continue
        if not _tier_allowed(str(row.get("Intact_Quality_Tier") or "Tier_4_rejected"), output_config.get("minimum_quality_tier", "Tier_3_weak")):
            continue
        mass = row.get("Reconstructed_Mass")
        intensity = _safe_float(row.get("Reconstructed_Envelope_Intensity"))
        if mass is None or mass == "" or intensity <= 0:
            continue
        rows.append({
            "Spectrum_Point_Rank": None,
            "Reconstructed_Mass_Da": _safe_float(mass),
            "Reconstructed_Envelope_Intensity": intensity,
            "Relative_Intensity_Percent": None,
            "Intensity_Method": row.get("Intensity_Method") or output_config["intensity_method"],
            "Cluster_ID": row.get("Cluster_ID"),
            "Reconstruction_Engine": row.get("Reconstruction_Engine"),
            "Intact_Envelope_Group_ID": row.get("Intact_Envelope_Group_ID"),
            "Group_Representative": row.get("Group_Representative"),
            "Comparison_Representative": row.get("Comparison_Representative"),
            "Reconstruction_Status": row.get("Reconstruction_Status"),
            "Intact_Quality_Tier": row.get("Intact_Quality_Tier"),
            "Quality_Tier_Rank": row.get("Quality_Tier_Rank"),
            "Peak_Sharing_Status": row.get("Peak_Sharing_Status"),
            "Max_Peak_Usage_Count": row.get("Max_Peak_Usage_Count"),
            "Competing_Envelope_Group_ID": row.get("Competing_Envelope_Group_ID"),
            "Competing_Envelope_Group_Size": row.get("Competing_Envelope_Group_Size"),
            "Envelope_Evidence_Score": row.get("Envelope_Evidence_Score"),
            "Evidence_Score_Rank_In_Competition": row.get("Evidence_Score_Rank_In_Competition"),
            "Maximum_Shared_Peak_Fraction": row.get("Maximum_Shared_Peak_Fraction"),
            "Envelope_QC_Eligible": row.get("Envelope_QC_Eligible"),
            "Intact_Strict_Eligible": row.get("Intact_Strict_Eligible"),
            "Intact_Review_Eligible": row.get("Intact_Review_Eligible"),
            "Comparison_Ready": row.get("Comparison_Ready"),
            "Num_Supporting_Charge_States": row.get("Num_Supporting_Charge_States"),
            "Supporting_Charge_States": row.get("Supporting_Charge_States"),
            "RT_Mean": row.get("RT_Mean"),
            "RT_Range_Min": row.get("RT_Range_Min"),
            "Envelope_Internal_Error_ppm": row.get("Envelope_Internal_Error_ppm"),
            "Neutral_Mass_SD": row.get("Neutral_Mass_SD"),
            "Neutral_Mass_Range": row.get("Neutral_Mass_Range"),
            "Best_Reference_Label": row.get("Best_Reference_Label"),
            "Reference_Mass_Error_ppm": row.get("Reference_Mass_Error_ppm"),
            "Limiting_Factors": row.get("Limiting_Factors"),
        })
    ranked = sorted(rows, key=lambda item: _safe_float(item.get("Reconstructed_Envelope_Intensity")), reverse=True)
    for rank, row in enumerate(ranked, start=1):
        row["Spectrum_Point_Rank"] = rank
    max_intensity = max((_safe_float(row.get("Reconstructed_Envelope_Intensity")) for row in rows), default=0.0)
    for row in rows:
        if output_config.get("normalize_to_percent", True):
            row["Relative_Intensity_Percent"] = (_safe_float(row.get("Reconstructed_Envelope_Intensity")) / max_intensity * 100.0) if max_intensity else 0.0
        else:
            row["Relative_Intensity_Percent"] = None
    return sorted(rows, key=lambda item: _safe_float(item.get("Reconstructed_Mass_Da")))


def _candidate_projection(row: dict[str, Any], columns: list[str]) -> dict[str, Any]:
    return {column: row.get(column) for column in columns}


def build_intact_comparison_candidate_rows(qc_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = sorted([row for row in qc_rows if row.get("Comparison_Representative")], key=lambda row: _safe_float(row.get("Comparison_Representative_Rank")))
    return [_candidate_projection(row, COMPARISON_CANDIDATE_COLUMNS) for row in rows]


def build_target_review_candidate_rows(qc_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = sorted([row for row in qc_rows if row.get("Target_Review_Group_Representative")], key=lambda row: _safe_float(row.get("Target_Review_Rank")))
    return [_candidate_projection(row, TARGET_REVIEW_CANDIDATE_COLUMNS) for row in rows]


def build_intact_reconstruction_qc(
    candidates: list[IntactMassCandidate],
    charge_state_peaks: list[dict[str, Any]],
    reconstruction_config: dict[str, Any],
    reconstruction_enabled: bool = True,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    qc_config = _qc_config(reconstruction_config or {})
    spectrum_output_config = qc_config["mass_spectrum_output"]
    intensity_method = spectrum_output_config["intensity_method"]
    peaks_by_cluster: dict[str, list[dict[str, Any]]] = {}
    for row in charge_state_peaks or []:
        peaks_by_cluster.setdefault(str(row.get("Cluster_ID") or ""), []).append(row)

    max_intensity = max((float(getattr(candidate, "total_intensity", 0.0) or 0.0) for candidate in candidates), default=0.0)
    qc_rows: list[dict[str, Any]] = []
    for candidate in candidates:
        cluster_id = candidate.cluster_id or ""
        cluster_peaks = peaks_by_cluster.get(cluster_id, [])
        supporting_peak_ids = sorted({_supporting_peak_id(row) for row in cluster_peaks})
        supporting_local_peak_ids = sorted({str(row.get("Local_Peak_ID")) for row in cluster_peaks if row.get("Local_Peak_ID") not in {None, ""}})
        if not supporting_local_peak_ids:
            supporting_local_peak_ids = supporting_peak_ids
        supporting_scan_ids = sorted({str(row.get("Scan_ID")) for row in cluster_peaks if row.get("Scan_ID") not in {None, ""}})
        supporting_rt_values = []
        for row in cluster_peaks:
            if row.get("RT") is not None and row.get("RT") != "":
                supporting_rt_values.append(_safe_float(row.get("RT")))
        charges = sorted({int(charge) for charge in candidate.charge_states})
        supporting_charge_states = sorted({int(row.get("Charge")) for row in cluster_peaks if row.get("Charge") is not None} or set(charges))
        supporting_peak_charge_map: dict[str, set[int]] = {}
        use_local_peak_ids = any(row.get("Local_Peak_ID") not in {None, ""} for row in cluster_peaks)
        for peak_row in cluster_peaks:
            charge = peak_row.get("Charge")
            if charge is None:
                continue
            if use_local_peak_ids:
                peak_id = peak_row.get("Local_Peak_ID")
                if peak_id in {None, ""}:
                    continue
            else:
                peak_id = _supporting_peak_id(peak_row)
            supporting_peak_charge_map.setdefault(str(peak_id), set()).add(int(charge))
        exact_peak_set_key = ";".join(supporting_peak_ids) if supporting_peak_ids else f"cluster:{cluster_id}"
        neutral_masses = [float(row.get("Neutral_Mass")) for row in cluster_peaks if row.get("Neutral_Mass") is not None]
        if not neutral_masses and candidate.observed_mass is not None:
            neutral_masses = [float(candidate.observed_mass)]
        reconstructed_mass = float(candidate.observed_mass)
        neutral_sd = pstdev(neutral_masses) if len(neutral_masses) > 1 else 0.0 if neutral_masses else None
        neutral_range = (max(neutral_masses) - min(neutral_masses)) if neutral_masses else None
        envelope_internal_error_ppm = _max_abs_ppm(neutral_masses, reconstructed_mass)
        continuity = _charge_continuity(charges)
        rt_min, rt_max, rt_mean, rt_range, rt_consistency = _rt_metrics(cluster_peaks, qc_config)
        peak_classes, trace_only = _class_summary(cluster_peaks)
        intensities = []
        for row in cluster_peaks:
            try:
                intensities.append(float(row.get("Intensity") or 0.0))
            except (TypeError, ValueError):
                continue
        total_intensity = float(candidate.total_intensity or sum(intensities) or 0.0)
        mean_intensity = sum(intensities) / len(intensities) if intensities else total_intensity / max(len(charges), 1)
        max_supporting_intensity = max(intensities) if intensities else total_intensity
        reconstructed_intensity = _reconstructed_intensity(total_intensity, mean_intensity, max_supporting_intensity, intensity_method)
        relative_intensity = (total_intensity / max_intensity * 100.0) if max_intensity else 0.0
        competing = sum(
            1
            for other in candidates
            if other is not candidate
            and abs(float(other.observed_mass) - reconstructed_mass) <= qc_config["max_neutral_mass_range_da"]
        )
        raw_local_envelope_relative = _candidate_extra(candidate, "local_envelope_relative_intensity_percent", None)
        local_envelope_relative = float(raw_local_envelope_relative or 0.0)
        if local_envelope_relative <= 0.0:
            local_envelope_relative = relative_intensity
        max_peak_usage = int(_candidate_extra(candidate, "max_peak_usage_count", _candidate_extra(candidate, "peak_usage_count", 0)) or 0)
        mean_peak_usage = float(_candidate_extra(candidate, "mean_peak_usage_count", max_peak_usage) or 0.0)
        highly_shared_count = int(_candidate_extra(candidate, "num_highly_shared_peaks", _candidate_extra(candidate, "shared_peak_count", 0)) or 0)
        highly_shared_fraction = float(_candidate_extra(candidate, "highly_shared_peak_fraction", _candidate_extra(candidate, "shared_peak_fraction", 0.0)) or 0.0)
        sharing_cfg = qc_config["rt_localized"]["peak_sharing"]
        if highly_shared_fraction > sharing_cfg["max_highly_shared_fraction_for_tier2"]:
            peak_sharing_status = "highly_shared"
        elif highly_shared_fraction > sharing_cfg["max_highly_shared_fraction_for_tier1"]:
            peak_sharing_status = "moderately_shared"
        else:
            peak_sharing_status = "low_shared"
        unmodified_delta_da = candidate.mass_error_da
        unmodified_delta_ppm = candidate.mass_error_ppm
        reference_label, reference_mass, reference_error_da, reference_error_ppm, reference_matched = _reference_match(reconstructed_mass, qc_config)
        in_mass_range, neutral_search_min, neutral_search_max, neutral_range_status = _neutral_mass_range_status(reconstructed_mass, qc_config)
        in_target_range, target_range_status, target_review_priority = _target_review_range_status(reconstructed_mass, qc_config)

        factors: list[str] = []
        if len(charges) < qc_config["min_charge_states_for_review"]:
            factors.append("insufficient_charge_states")
        if qc_config["require_contiguous_charge_states"] and continuity == "non_contiguous":
            factors.append("non_contiguous_charge_states")
        if neutral_sd is not None and neutral_sd > qc_config["max_neutral_mass_sd_da"]:
            factors.append("mass_spread_too_large")
        if neutral_range is not None and neutral_range > qc_config["max_neutral_mass_range_da"]:
            factors.append("mass_spread_too_large")
        if envelope_internal_error_ppm is not None and envelope_internal_error_ppm > qc_config["max_envelope_internal_error_ppm"]:
            factors.append("internal_mass_error_too_large")
        if rt_consistency == "inconsistent":
            factors.append("rt_inconsistent")
        elif rt_consistency == "not_available":
            factors.append("rt_not_available")
        if relative_intensity < qc_config["min_relative_envelope_intensity_percent_for_review"]:
            factors.append("insufficient_intensity_support")
        if trace_only:
            factors.append("trace_only_envelope")
        if competing > qc_config["max_competing_envelopes"]:
            factors.append("multiple_competing_envelopes")
        if peak_sharing_status == "highly_shared":
            factors.append("high_peak_sharing")
        if not in_mass_range:
            factors.append("outside_neutral_mass_search_range")
        factors = list(dict.fromkeys(factors))
        severe_factors = [factor for factor in factors if factor in SEVERE_LIMITING_FACTORS]

        basic_internal_ok = (
            (neutral_sd is None or neutral_sd <= qc_config["max_neutral_mass_sd_da"])
            and (neutral_range is None or neutral_range <= qc_config["max_neutral_mass_range_da"])
            and (envelope_internal_error_ppm is None or envelope_internal_error_ppm <= qc_config["max_envelope_internal_error_ppm"])
        )
        contiguous_ok = continuity == "contiguous" or not qc_config["require_contiguous_charge_states"]
        reliable_intensity_ok = relative_intensity >= qc_config["min_relative_envelope_intensity_percent_for_reliable"]
        review_intensity_ok = relative_intensity >= qc_config["min_relative_envelope_intensity_percent_for_review"]
        reliable_rt_ok = rt_consistency == "consistent"
        review_rt_ok = rt_consistency in {"consistent", "review"}
        trace_ok_for_reliable = not trace_only or qc_config["allow_trace_only_reliable"]

        if not reconstruction_enabled:
            status = "Failed"
            factors.insert(0, "reconstruction_disabled")
            severe_factors.insert(0, "reconstruction_disabled")
        elif len(charges) < qc_config["min_charge_states_for_review"]:
            status = "Insufficient"
        elif (
            len(charges) >= qc_config["min_charge_states_for_reliable"]
            and contiguous_ok
            and basic_internal_ok
            and reliable_rt_ok
            and reliable_intensity_ok
            and trace_ok_for_reliable
            and competing <= qc_config["max_competing_envelopes"]
        ):
            status = "Reliable"
        else:
            status = "Review"
        factors = list(dict.fromkeys(factors))
        severe_factors = list(dict.fromkeys(severe_factors))
        confidence = {"Reliable": "High", "Review": "Medium", "Insufficient": "Low", "Failed": "None"}.get(status, "Low")
        envelope_qc_eligible = (
            in_mass_range
            and reconstruction_enabled
            and len(charges) >= qc_config["min_charge_states_for_review"]
            and basic_internal_ok
            and review_rt_ok
            and not severe_factors
        )
        intact_strict_eligible = (
            envelope_qc_eligible
            and status == "Reliable"
            and contiguous_ok
            and reliable_rt_ok
            and reliable_intensity_ok
            and trace_ok_for_reliable
        )
        intact_review_eligible = (
            envelope_qc_eligible
            and not intact_strict_eligible
            and "Review" in qc_config["comparison_ready_statuses"]
            and review_intensity_ok
        )
        quality = _quality_matrix_and_tier(
            charges=charges,
            continuity=continuity,
            neutral_sd=neutral_sd,
            neutral_range=neutral_range,
            envelope_internal_error_ppm=envelope_internal_error_ppm,
            rt_consistency=rt_consistency,
            local_envelope_relative=local_envelope_relative,
            competing=competing,
            severe_factors=severe_factors,
            reconstruction_enabled=reconstruction_enabled,
            in_mass_range=in_mass_range,
            qc_config=qc_config,
            candidate=candidate,
        )
        tier = quality["Intact_Quality_Tier"]
        comparison_ready_strict = tier in qc_config["comparison_ready_tiers"]["strict"]
        comparison_ready_review = tier in qc_config["comparison_ready_tiers"]["review"]
        comparison_ready = comparison_ready_strict or comparison_ready_review
        intact_strict_eligible = intact_strict_eligible and comparison_ready_strict
        intact_review_eligible = (intact_review_eligible or comparison_ready_review) and comparison_ready_review
        readiness_reason = "strict" if comparison_ready_strict else "review" if comparison_ready_review else quality["Quality_Tier_Reason"] or _primary_factor(factors) or "not_ready"
        primary_factor = _primary_factor(factors)

        candidate.reconstruction_status = status
        candidate.reconstruction_confidence = confidence
        candidate.num_supporting_charge_states = len(charges)
        candidate.charge_state_range = _charge_state_range(charges)
        candidate.charge_state_continuity = continuity
        candidate.neutral_mass_sd = neutral_sd
        candidate.neutral_mass_range = neutral_range
        candidate.envelope_internal_error_ppm = envelope_internal_error_ppm
        candidate.max_mass_error_ppm = envelope_internal_error_ppm
        candidate.unmodified_theory_delta_da = unmodified_delta_da
        candidate.unmodified_theory_delta_ppm = unmodified_delta_ppm
        candidate.best_reference_label = reference_label
        candidate.best_reference_mass_da = reference_mass
        candidate.reference_mass_error_da = reference_error_da
        candidate.reference_mass_error_ppm = reference_error_ppm
        candidate.reference_mass_matched = reference_matched
        candidate.in_neutral_mass_search_range = in_mass_range
        candidate.neutral_mass_search_min_da = neutral_search_min
        candidate.neutral_mass_search_max_da = neutral_search_max
        candidate.neutral_mass_range_status = neutral_range_status
        candidate.in_target_review_mass_range = in_target_range
        candidate.target_review_mass_range_status = target_range_status
        candidate.target_review_priority = target_review_priority
        candidate.envelope_qc_eligible = envelope_qc_eligible
        candidate.intact_review_eligible = intact_review_eligible
        candidate.intact_strict_eligible = intact_strict_eligible
        candidate.supporting_peak_ids = "; ".join(supporting_peak_ids)
        candidate.supporting_local_peak_ids = "; ".join(supporting_local_peak_ids)
        candidate.supporting_scan_ids = "; ".join(supporting_scan_ids)
        candidate.supporting_rt_values = "; ".join(f"{value:.5f}" for value in supporting_rt_values)
        candidate.supporting_charge_states = "; ".join(map(str, supporting_charge_states))
        candidate.exact_peak_set_key = exact_peak_set_key
        candidate.rt_min = rt_min
        candidate.rt_max = rt_max
        candidate.rt_mean = rt_mean
        candidate.rt_range_min = rt_range
        candidate.max_rt_difference_min = rt_range
        candidate.rt_consistency = rt_consistency
        candidate.total_supporting_intensity = total_intensity
        candidate.mean_supporting_intensity = mean_intensity
        candidate.max_supporting_intensity = max_supporting_intensity
        candidate.reconstructed_envelope_intensity = reconstructed_intensity
        candidate.intensity_method = intensity_method
        candidate.relative_envelope_intensity_percent = relative_intensity
        candidate.supporting_peak_classes = peak_classes
        candidate.trace_only_envelope = trace_only
        candidate.competing_envelope_count = competing
        candidate.limiting_factors = "; ".join(factors)
        candidate.severe_limiting_factors = "; ".join(severe_factors)
        candidate.num_limiting_factors = len(factors)
        candidate.primary_limiting_factor = primary_factor
        candidate.comparison_ready_strict = comparison_ready_strict
        candidate.comparison_ready_review = comparison_ready_review
        candidate.comparison_ready = comparison_ready
        candidate.comparison_readiness_reason = readiness_reason
        candidate.max_peak_usage_count = max_peak_usage
        candidate.mean_peak_usage_count = mean_peak_usage
        candidate.num_highly_shared_peaks = highly_shared_count
        candidate.highly_shared_peak_fraction = highly_shared_fraction
        candidate.competing_candidate_count = competing
        candidate.peak_sharing_status = peak_sharing_status
        candidate.pass_min_charge_count = quality["Pass_Min_Charge_Count"]
        candidate.pass_min_consecutive_charge_count = quality["Pass_Min_Consecutive_Charge_Count"]
        candidate.pass_charge_continuity = quality["Pass_Charge_Continuity"]
        candidate.pass_internal_error = quality["Pass_Internal_Error"]
        candidate.pass_neutral_mass_sd = quality["Pass_Neutral_Mass_SD"]
        candidate.pass_neutral_mass_range = quality["Pass_Neutral_Mass_Range"]
        candidate.pass_rt_consistency = quality["Pass_RT_Consistency"]
        candidate.pass_local_intensity = quality["Pass_Local_Intensity"]
        candidate.pass_competing_envelope = quality["Pass_Competing_Envelope"]
        candidate.pass_peak_sharing = quality["Pass_Peak_Sharing"]
        candidate.num_strict_criteria_passed = quality["Num_Strict_Criteria_Passed"]
        candidate.num_review_criteria_passed = quality["Num_Review_Criteria_Passed"]
        candidate.strict_failure_reasons = quality["Strict_Failure_Reasons"]
        candidate.review_failure_reasons = quality["Review_Failure_Reasons"]
        candidate.intact_quality_tier = quality["Intact_Quality_Tier"]
        candidate.quality_tier_reason = quality["Quality_Tier_Reason"]
        candidate.quality_tier_rank = quality["Quality_Tier_Rank"]

        qc_rows.append({
            "Cluster_ID": cluster_id,
            "Reconstructed_Mass": reconstructed_mass,
            "Observed_Mass": reconstructed_mass,
            "In_Neutral_Mass_Search_Range": in_mass_range,
            "Neutral_Mass_Search_Min_Da": neutral_search_min,
            "Neutral_Mass_Search_Max_Da": neutral_search_max,
            "Neutral_Mass_Range_Status": neutral_range_status,
            "In_Target_Review_Mass_Range": in_target_range,
            "Target_Review_Mass_Range_Status": target_range_status,
            "Target_Review_Priority": target_review_priority,
            "Envelope_QC_Eligible": envelope_qc_eligible,
            "Intact_Review_Eligible": intact_review_eligible,
            "Intact_Strict_Eligible": intact_strict_eligible,
            "Intact_Envelope_QC_Score": None,
            "Intact_Envelope_QC_Rank": None,
            "Strict_Eligible_Rank": None,
            "Review_Eligible_Rank": None,
            "Dominant_Intact_Envelope_Flag": False,
            "Supporting_Peak_IDs": "; ".join(supporting_peak_ids),
            "Supporting_Peak_Count": len(supporting_peak_ids) or len(cluster_peaks),
            "Supporting_Scan_IDs": "; ".join(supporting_scan_ids),
            "Supporting_RT_Values": "; ".join(f"{value:.5f}" for value in supporting_rt_values),
            "Supporting_Charge_States": "; ".join(map(str, supporting_charge_states)),
            "Exact_Peak_Set_Key": exact_peak_set_key,
            "Exact_Duplicate_Group_ID": "",
            "Exact_Duplicate_Count": 1,
            "Is_Exact_Duplicate_Representative": True,
            "Intact_Envelope_Group_ID": "",
            "Envelope_Group_Size": 1,
            "Shared_Peak_Count_With_Representative": 0,
            "Shared_Peak_Fraction_With_Representative": 0.0,
            "Shared_Charge_Count_With_Representative": 0,
            "Shared_Charge_Fraction_With_Representative": 0.0,
            "Mass_Delta_To_Group_Representative_Da": 0.0,
            "RT_Delta_To_Group_Representative_Min": 0.0,
            "Group_Representative": True,
            "Group_Ambiguity_Status": "unique",
            "Comparison_Representative": False,
            "Comparison_Representative_Reason": "",
            "Comparison_Representative_Rank": None,
            "Excluded_From_Comparison_Reason": "not_grouped",
            "Target_Review_Group_Representative": False,
            "Target_Review_Rank": None,
            "Dominant_Target_Review_Eligible_Flag": False,
            "_supporting_peak_id_set": set(supporting_peak_ids),
            "_supporting_local_peak_id_set": set(supporting_local_peak_ids),
            "_supporting_charge_set": set(supporting_charge_states),
            "_supporting_peak_charge_map": supporting_peak_charge_map,
            "_peak_charge_assignment_source": "observed" if supporting_peak_charge_map else "estimated",
            "Reconstruction_Status": status,
            "Reconstruction_Confidence": confidence,
            "Reconstruction_Engine": _candidate_extra(candidate, "reconstruction_engine", "legacy_cluster"),
            "RT_Window_ID": _candidate_extra(candidate, "rt_window_id", ""),
            "RT_Window_Start_Min": _candidate_extra(candidate, "rt_window_start_min"),
            "RT_Window_End_Min": _candidate_extra(candidate, "rt_window_end_min"),
            "RT_Window_Center_Min": _candidate_extra(candidate, "rt_window_center_min"),
            "Num_MS1_Scans_In_Window": _candidate_extra(candidate, "num_ms1_scans_in_window", 0),
            "Peak_Aggregation_Method": _candidate_extra(candidate, "peak_aggregation_method", ""),
            "Anchor_MZ": _candidate_extra(candidate, "anchor_mz"),
            "Anchor_Charge": _candidate_extra(candidate, "anchor_charge"),
            "Predicted_Charge_States": _candidate_extra(candidate, "predicted_charge_states", ""),
            "Observed_Charge_States": _candidate_extra(candidate, "observed_charge_states", ""),
            "Missing_Charge_States": _candidate_extra(candidate, "missing_charge_states", ""),
            "Missing_Charge_Predicted_MZ": _candidate_extra(candidate, "missing_charge_predicted_mz", ""),
            "Num_Predicted_Charges": _candidate_extra(candidate, "num_predicted_charges", 0),
            "Num_Observed_Charges": _candidate_extra(candidate, "num_observed_charges", 0),
            "Charge_Coverage_Fraction": _candidate_extra(candidate, "charge_coverage_fraction", 0.0),
            "Consecutive_Charge_Run_Length": _candidate_extra(candidate, "consecutive_charge_run_length", 0) or _longest_consecutive_run(charges),
            "Longest_Consecutive_Charge_Run": _candidate_extra(candidate, "longest_consecutive_charge_run", 0) or _longest_consecutive_run(charges),
            "Charge_Gap_Count": _candidate_extra(candidate, "charge_gap_count", 0) or _charge_gap_count(charges),
            "Charge_Continuity_Fraction": _candidate_extra(candidate, "charge_continuity_fraction", 0.0) or _charge_continuity_fraction_value(charges),
            "Peak_Usage_Count": _candidate_extra(candidate, "peak_usage_count", 0),
            "Shared_Peak_Count": _candidate_extra(candidate, "shared_peak_count", 0),
            "Shared_Peak_Fraction": _candidate_extra(candidate, "shared_peak_fraction", 0.0),
            "Local_Window_Max_Intensity": _candidate_extra(candidate, "local_window_max_intensity", 0.0),
            "Local_Relative_Peak_Intensity_Percent": _candidate_extra(candidate, "local_relative_peak_intensity_percent", 0.0),
            "Local_Envelope_Relative_Intensity_Percent": local_envelope_relative,
            "Neutral_Mass_Estimator": _candidate_extra(candidate, "neutral_mass_estimator", ""),
            "Neutral_Mass_Unweighted_Mean": _candidate_extra(candidate, "neutral_mass_unweighted_mean"),
            "Neutral_Mass_Weighted_Mean": _candidate_extra(candidate, "neutral_mass_weighted_mean"),
            "Neutral_Mass_Median": _candidate_extra(candidate, "neutral_mass_median"),
            "Envelope_Internal_Error_Max_ppm": _candidate_extra(candidate, "envelope_internal_error_max_ppm"),
            "Envelope_Internal_Error_Mean_ppm": _candidate_extra(candidate, "envelope_internal_error_mean_ppm"),
            "Envelope_Internal_Error_Median_ppm": _candidate_extra(candidate, "envelope_internal_error_median_ppm"),
            "Source_RT_Window_IDs": _candidate_extra(candidate, "source_rt_window_ids", ""),
            "Num_Source_RT_Windows": _candidate_extra(candidate, "num_source_rt_windows", 0),
            "Merged_Across_RT_Windows": _candidate_extra(candidate, "merged_across_rt_windows", False),
            "Extended_Lower_Charges_Evaluated": _candidate_extra(candidate, "extended_lower_charges_evaluated", ""),
            "Extended_Upper_Charges_Evaluated": _candidate_extra(candidate, "extended_upper_charges_evaluated", ""),
            "Extended_Charges_Detected": _candidate_extra(candidate, "extended_charges_detected", ""),
            "Extended_Weak_Charges_Detected": _candidate_extra(candidate, "extended_weak_charges_detected", ""),
            "Extended_Charges_Not_Detected": _candidate_extra(candidate, "extended_charges_not_detected", ""),
            "Charge_Extension_Improved_Envelope": _candidate_extra(candidate, "charge_extension_improved_envelope", False),
            "Original_Charge_States": _candidate_extra(candidate, "original_charge_states", ""),
            "Final_Charge_States": _candidate_extra(candidate, "final_charge_states", ""),
            "Split_Envelope_Group_ID": _candidate_extra(candidate, "split_envelope_group_id", ""),
            "Split_Envelope_Member_Count": _candidate_extra(candidate, "split_envelope_member_count", 1),
            "Split_Envelope_Merged": _candidate_extra(candidate, "split_envelope_merged", False),
            "Charge_Gaps_Before_Merge": _candidate_extra(candidate, "charge_gaps_before_merge", _candidate_extra(candidate, "charge_gap_count", 0)),
            "Charge_Gaps_After_Merge": _candidate_extra(candidate, "charge_gaps_after_merge", _candidate_extra(candidate, "charge_gap_count", 0)),
            "Max_Peak_Usage_Count": max_peak_usage,
            "Mean_Peak_Usage_Count": mean_peak_usage,
            "Num_Highly_Shared_Peaks": highly_shared_count,
            "Highly_Shared_Peak_Fraction": highly_shared_fraction,
            "Competing_Candidate_Count": competing,
            "Peak_Sharing_Status": peak_sharing_status,
            "Competing_Envelope_Group_ID": "",
            "Competing_Envelope_Group_Size": 1,
            "Shared_Peak_Competitor_Count": 0,
            "Maximum_Shared_Peak_Fraction": 0.0,
            "Mean_Shared_Peak_Fraction": 0.0,
            "Competitor_Cluster_IDs": "",
            "Is_Noncompeting_Candidate": True,
            "Envelope_Evidence_Score": None,
            "Evidence_Score_Rank_In_Competition": None,
            "Evidence_Score_Components": "",
            "Evidence_Score_Penalties": "",
            "Evidence_Score_Config_Version": "",
            "Direct_Competitor_Count": 0,
            "Direct_Competitor_Cluster_IDs": "",
            "Direct_Shared_Peak_Count_Max": 0,
            "Direct_Shared_Peak_Fraction_Max": 0.0,
            "Competition_Component_Size": 1,
            "Dry_Run_Assignment_Status": "not_evaluated",
            "Dry_Run_Selected": False,
            "Dry_Run_Selection_Order": None,
            "Supporting_Peak_Count_Before_Assignment": len(supporting_local_peak_ids),
            "Independent_Supporting_Peak_Count": 0,
            "Independent_Supporting_Peak_Fraction": 0.0,
            "Supporting_Charge_Count_Before_Assignment": len(charges),
            "Independent_Charge_State_Count": 0,
            "Peaks_Already_Assigned_Count": 0,
            "Charges_Already_Assigned_Count": 0,
            "Excluded_By_Cluster_ID": "",
            "Dry_Run_Exclusion_Reason": "",
            "Score_Margin_To_Excluding_Candidate": None,
            "Close_Score_Ambiguity": False,
            "Assignment_Confidence": "low",
            "Shared_Observed_Peak_Count": 0,
            "Shared_Peak_Charge_Assignment_Count": 0,
            "Independent_Observed_Peak_Count": 0,
            **quality,
            "Comparison_Ready_Strict": comparison_ready_strict,
            "Comparison_Ready_Review": comparison_ready_review,
            "Comparison_Ready": comparison_ready,
            "Comparison_Readiness_Reason": readiness_reason,
            "Total_Supporting_Intensity": total_intensity,
            "Mean_Supporting_Intensity": mean_intensity,
            "Max_Supporting_Intensity": max_supporting_intensity,
            "Reconstructed_Envelope_Intensity": reconstructed_intensity,
            "Intensity_Method": intensity_method,
            "Relative_Envelope_Intensity_Percent": relative_intensity,
            "Relative_Overall_Envelope_Intensity_Percent": relative_intensity,
            "Relative_In_Range_Raw_Intensity_Percent": None,
            "Relative_Intact_Eligible_Intensity_Percent": None,
            "Supporting_Peak_Classes": peak_classes,
            "Trace_Only_Envelope": trace_only,
            "Num_Supporting_Charge_States": len(charges),
            "Charge_State_Range": candidate.charge_state_range,
            "Charge_State_Continuity": continuity,
            "RT_Min": rt_min,
            "RT_Max": rt_max,
            "RT_Mean": rt_mean,
            "RT_Range_Min": rt_range,
            "Max_RT_Difference_Min": rt_range,
            "RT_Consistency": rt_consistency,
            "Neutral_Mass_SD": neutral_sd,
            "Neutral_Mass_Range": neutral_range,
            "Envelope_Internal_Error_ppm": envelope_internal_error_ppm,
            "Max_Mass_Error_ppm": envelope_internal_error_ppm,
            "Unmodified_Theory_Delta_Da": unmodified_delta_da,
            "Unmodified_Theory_Delta_ppm": unmodified_delta_ppm,
            "Best_Reference_Label": reference_label,
            "Best_Reference_Mass_Da": reference_mass,
            "Reference_Mass_Error_Da": reference_error_da,
            "Reference_Mass_Error_ppm": reference_error_ppm,
            "Reference_Mass_Matched": reference_matched,
            "Competing_Envelope_Count": competing,
            "Limiting_Factors": candidate.limiting_factors,
            "Severe_Limiting_Factors": candidate.severe_limiting_factors,
            "Num_Limiting_Factors": len(factors),
            "Primary_Limiting_Factor": primary_factor,
        })

    max_in_range_intensity = max((_safe_float(row.get("Total_Supporting_Intensity")) for row in qc_rows if row.get("In_Neutral_Mass_Search_Range")), default=0.0)
    eligible_rows = [row for row in qc_rows if row.get("Intact_Strict_Eligible") or row.get("Intact_Review_Eligible")]
    max_eligible_intensity = max((_safe_float(row.get("Total_Supporting_Intensity")) for row in eligible_rows), default=0.0)
    for row in qc_rows:
        total = _safe_float(row.get("Total_Supporting_Intensity"))
        row["Relative_In_Range_Raw_Intensity_Percent"] = (total / max_in_range_intensity * 100.0) if max_in_range_intensity and row.get("In_Neutral_Mass_Search_Range") else 0.0
        row["Relative_Intact_Eligible_Intensity_Percent"] = (total / max_eligible_intensity * 100.0) if max_eligible_intensity and (row.get("Intact_Strict_Eligible") or row.get("Intact_Review_Eligible")) else 0.0
        row["Intact_Envelope_QC_Score"] = _intact_qc_score(row, qc_config)

    ranked_rows = sorted(qc_rows, key=_dominant_sort_key, reverse=True)
    for rank, row in enumerate(ranked_rows, start=1):
        row["Intact_Envelope_QC_Rank"] = rank
    strict_rows = sorted([row for row in qc_rows if row.get("Intact_Strict_Eligible")], key=_dominant_sort_key, reverse=True)
    for rank, row in enumerate(strict_rows, start=1):
        row["Strict_Eligible_Rank"] = rank
    review_rows = sorted([row for row in qc_rows if row.get("Intact_Review_Eligible")], key=_dominant_sort_key, reverse=True)
    for rank, row in enumerate(review_rows, start=1):
        row["Review_Eligible_Rank"] = rank
    dominant_eligible = strict_rows[0] if strict_rows else review_rows[0] if review_rows else None
    if dominant_eligible is not None:
        dominant_eligible["Dominant_Intact_Envelope_Flag"] = True

    apply_intact_envelope_grouping(qc_rows, qc_config)
    competition_stats = apply_competitive_assignment(qc_rows, qc_config)
    assignment_stats = apply_assignment_dry_run(qc_rows, qc_config)
    competition_stats.update(assignment_stats)
    reconstruction_config["_intact_competition_stats"] = competition_stats

    candidates_by_cluster = {candidate.cluster_id or "": candidate for candidate in candidates}
    for row in qc_rows:
        candidate = candidates_by_cluster.get(str(row.get("Cluster_ID") or ""))
        if candidate is None:
            continue
        candidate.relative_overall_envelope_intensity_percent = row["Relative_Overall_Envelope_Intensity_Percent"]
        candidate.relative_in_range_raw_intensity_percent = row["Relative_In_Range_Raw_Intensity_Percent"]
        candidate.relative_intact_eligible_intensity_percent = row["Relative_Intact_Eligible_Intensity_Percent"]
        candidate.reconstructed_envelope_intensity = row["Reconstructed_Envelope_Intensity"]
        candidate.intensity_method = row["Intensity_Method"]
        candidate.intact_envelope_qc_score = row["Intact_Envelope_QC_Score"]
        candidate.intact_envelope_qc_rank = row["Intact_Envelope_QC_Rank"]
        candidate.strict_eligible_rank = row["Strict_Eligible_Rank"]
        candidate.review_eligible_rank = row["Review_Eligible_Rank"]
        candidate.dominant_intact_envelope_flag = row["Dominant_Intact_Envelope_Flag"]
        candidate.exact_duplicate_group_id = row["Exact_Duplicate_Group_ID"]
        candidate.exact_duplicate_count = row["Exact_Duplicate_Count"]
        candidate.is_exact_duplicate_representative = row["Is_Exact_Duplicate_Representative"]
        candidate.intact_envelope_group_id = row["Intact_Envelope_Group_ID"]
        candidate.envelope_group_size = row["Envelope_Group_Size"]
        candidate.group_representative = row["Group_Representative"]
        candidate.group_ambiguity_status = row["Group_Ambiguity_Status"]
        candidate.comparison_representative = row["Comparison_Representative"]
        candidate.comparison_representative_reason = row["Comparison_Representative_Reason"]
        candidate.comparison_representative_rank = row["Comparison_Representative_Rank"]
        candidate.excluded_from_comparison_reason = row["Excluded_From_Comparison_Reason"]
        candidate.target_review_group_representative = row["Target_Review_Group_Representative"]
        candidate.target_review_rank = row["Target_Review_Rank"]
        candidate.dominant_target_review_eligible_flag = row["Dominant_Target_Review_Eligible_Flag"]
        candidate.intact_quality_tier = row.get("Intact_Quality_Tier", candidate.intact_quality_tier)
        candidate.quality_tier_rank = row.get("Quality_Tier_Rank", candidate.quality_tier_rank)
        candidate.peak_sharing_status = row.get("Peak_Sharing_Status", candidate.peak_sharing_status)
        candidate.competing_envelope_group_id = row.get("Competing_Envelope_Group_ID", candidate.competing_envelope_group_id)
        candidate.competing_envelope_group_size = row.get("Competing_Envelope_Group_Size", candidate.competing_envelope_group_size)
        candidate.shared_peak_competitor_count = row.get("Shared_Peak_Competitor_Count", candidate.shared_peak_competitor_count)
        candidate.maximum_shared_peak_fraction = row.get("Maximum_Shared_Peak_Fraction", candidate.maximum_shared_peak_fraction)
        candidate.mean_shared_peak_fraction = row.get("Mean_Shared_Peak_Fraction", candidate.mean_shared_peak_fraction)
        candidate.competitor_cluster_ids = row.get("Competitor_Cluster_IDs", candidate.competitor_cluster_ids)
        candidate.is_noncompeting_candidate = row.get("Is_Noncompeting_Candidate", candidate.is_noncompeting_candidate)
        candidate.envelope_evidence_score = row.get("Envelope_Evidence_Score", candidate.envelope_evidence_score)
        candidate.evidence_score_rank_in_competition = row.get("Evidence_Score_Rank_In_Competition", candidate.evidence_score_rank_in_competition)
        candidate.evidence_score_components = row.get("Evidence_Score_Components", candidate.evidence_score_components)
        candidate.evidence_score_penalties = row.get("Evidence_Score_Penalties", candidate.evidence_score_penalties)
        candidate.evidence_score_config_version = row.get("Evidence_Score_Config_Version", candidate.evidence_score_config_version)
        candidate.direct_competitor_count = row.get("Direct_Competitor_Count", candidate.direct_competitor_count)
        candidate.direct_competitor_cluster_ids = row.get("Direct_Competitor_Cluster_IDs", candidate.direct_competitor_cluster_ids)
        candidate.direct_shared_peak_count_max = row.get("Direct_Shared_Peak_Count_Max", candidate.direct_shared_peak_count_max)
        candidate.direct_shared_peak_fraction_max = row.get("Direct_Shared_Peak_Fraction_Max", candidate.direct_shared_peak_fraction_max)
        candidate.competition_component_size = row.get("Competition_Component_Size", candidate.competition_component_size)
        candidate.dry_run_assignment_status = row.get("Dry_Run_Assignment_Status", candidate.dry_run_assignment_status)
        candidate.dry_run_selected = row.get("Dry_Run_Selected", candidate.dry_run_selected)
        candidate.dry_run_selection_order = row.get("Dry_Run_Selection_Order", candidate.dry_run_selection_order)
        candidate.supporting_peak_count_before_assignment = row.get("Supporting_Peak_Count_Before_Assignment", candidate.supporting_peak_count_before_assignment)
        candidate.independent_supporting_peak_count = row.get("Independent_Supporting_Peak_Count", candidate.independent_supporting_peak_count)
        candidate.independent_supporting_peak_fraction = row.get("Independent_Supporting_Peak_Fraction", candidate.independent_supporting_peak_fraction)
        candidate.supporting_charge_count_before_assignment = row.get("Supporting_Charge_Count_Before_Assignment", candidate.supporting_charge_count_before_assignment)
        candidate.independent_charge_state_count = row.get("Independent_Charge_State_Count", candidate.independent_charge_state_count)
        candidate.peaks_already_assigned_count = row.get("Peaks_Already_Assigned_Count", candidate.peaks_already_assigned_count)
        candidate.charges_already_assigned_count = row.get("Charges_Already_Assigned_Count", candidate.charges_already_assigned_count)
        candidate.excluded_by_cluster_id = row.get("Excluded_By_Cluster_ID", candidate.excluded_by_cluster_id)
        candidate.dry_run_exclusion_reason = row.get("Dry_Run_Exclusion_Reason", candidate.dry_run_exclusion_reason)
        candidate.score_margin_to_excluding_candidate = row.get("Score_Margin_To_Excluding_Candidate", candidate.score_margin_to_excluding_candidate)
        candidate.close_score_ambiguity = row.get("Close_Score_Ambiguity", candidate.close_score_ambiguity)
        candidate.assignment_confidence = row.get("Assignment_Confidence", candidate.assignment_confidence)
        candidate.shared_observed_peak_count = row.get("Shared_Observed_Peak_Count", candidate.shared_observed_peak_count)
        candidate.shared_peak_charge_assignment_count = row.get("Shared_Peak_Charge_Assignment_Count", candidate.shared_peak_charge_assignment_count)
        candidate.independent_observed_peak_count = row.get("Independent_Observed_Peak_Count", candidate.independent_observed_peak_count)


    diagnostic_rows = build_intact_reconstruction_diagnostics(qc_rows, reconstruction_config, reconstruction_enabled)
    return qc_rows, diagnostic_rows


def build_intact_reconstruction_diagnostics(
    qc_rows: list[dict[str, Any]],
    reconstruction_config: dict[str, Any],
    reconstruction_enabled: bool = True,
) -> list[dict[str, Any]]:
    qc_config = _qc_config(reconstruction_config or {})
    engine_stats = (reconstruction_config or {}).get("_intact_engine_stats") or {}
    competition_stats = (reconstruction_config or {}).get("_intact_competition_stats") or {}
    status_counts = {status: 0 for status in ["Reliable", "Review", "Insufficient", "Failed"]}
    reason_counts: dict[str, int] = {}
    for row in qc_rows:
        status = str(row.get("Reconstruction_Status") or "")
        if status in status_counts:
            status_counts[status] += 1
        for reason in str(row.get("Limiting_Factors") or row.get("Primary_Limiting_Factor") or "").split(";"):
            reason = reason.strip()
            if reason:
                reason_counts[reason] = reason_counts.get(reason, 0) + 1
    if not reconstruction_enabled:
        reason_counts["reconstruction_disabled"] = max(1, reason_counts.get("reconstruction_disabled", 0))
    elif not qc_rows:
        reason_counts["no_charge_state_candidates"] = 1
    reason_summary = "; ".join(f"{key}:{value}" for key, value in sorted(reason_counts.items()))
    in_range_rows = [row for row in qc_rows if row.get("In_Neutral_Mass_Search_Range")]
    outside_rows = [row for row in qc_rows if not row.get("In_Neutral_Mass_Search_Range")]
    dominant = max(qc_rows, key=lambda row: _safe_float(row.get("Total_Supporting_Intensity")), default={})
    dominant_in_range = max(in_range_rows, key=lambda row: _safe_float(row.get("Total_Supporting_Intensity")), default={})
    strict_rows = [row for row in qc_rows if row.get("Intact_Strict_Eligible")]
    review_rows = [row for row in qc_rows if row.get("Intact_Review_Eligible")]
    dominant_strict = _dominant_row(strict_rows)
    dominant_review = _dominant_row(review_rows)
    dominant_eligible = dominant_strict or dominant_review
    target_review_rows = [row for row in qc_rows if row.get("In_Target_Review_Mass_Range")]
    exact_duplicate_group_ids = {
        row.get("Exact_Duplicate_Group_ID")
        for row in qc_rows
        if row.get("Exact_Duplicate_Group_ID") and _safe_float(row.get("Exact_Duplicate_Count"), 1) > 1
    }
    duplicate_candidate_count = sum(1 for row in qc_rows if _safe_float(row.get("Exact_Duplicate_Count"), 1) > 1)
    group_ids = {row.get("Intact_Envelope_Group_ID") for row in qc_rows if row.get("Intact_Envelope_Group_ID")}
    ambiguity_counts: dict[str, int] = {}
    for row in qc_rows:
        if row.get("Group_Representative"):
            status = str(row.get("Group_Ambiguity_Status") or "unique")
            ambiguity_counts[status] = ambiguity_counts.get(status, 0) + 1
    comparison_reps = [row for row in qc_rows if row.get("Comparison_Representative")]
    target_review_reps = [row for row in qc_rows if row.get("Target_Review_Group_Representative")]
    dominant_comparison = comparison_reps[0] if comparison_reps else {}
    dominant_target = target_review_reps[0] if target_review_reps else {}
    grouping_config = qc_config["envelope_grouping"]
    references = qc_config.get("reference_masses") or []
    reference_summary = "; ".join(f"{item['label']}={item['mass_da']}" for item in references) or "not_configured"
    rt_settings = (
        f"reliable<={qc_config['max_rt_range_min_for_reliable']} min; "
        f"review<={qc_config['max_rt_range_min_for_review']} min"
    )
    tier_summary = _distribution_summary([row.get("Intact_Quality_Tier") for row in qc_rows])
    strict_failure_summary = _reason_summary([row.get("Strict_Failure_Reasons") for row in qc_rows])
    review_failure_summary = _reason_summary([row.get("Review_Failure_Reasons") for row in qc_rows])
    charge_count_summary = _distribution_summary([row.get("Num_Supporting_Charge_States") for row in qc_rows])
    consecutive_summary = _distribution_summary([row.get("Longest_Consecutive_Charge_Run") for row in qc_rows])
    internal_error_summary = _binned_distribution([row.get("Envelope_Internal_Error_ppm") for row in qc_rows], [(5, "<=5ppm"), (20, "<=20ppm"), (100, "<=100ppm")])
    rt_range_summary = _binned_distribution([row.get("RT_Range_Min") for row in qc_rows], [(0.05, "<=0.05min"), (0.15, "<=0.15min"), (0.30, "<=0.30min")])
    peak_usage_summary = _distribution_summary([row.get("Max_Peak_Usage_Count") for row in qc_rows])
    missing_status_summary = engine_stats.get("Missing_Charge_Status_Count", "")
    engine_match_summary = engine_stats.get("Engine_Match_Status_Count", "")
    return [{
        "Total_Reconstruction_Candidates": len(qc_rows),
        "Reliable_Count": status_counts["Reliable"],
        "Review_Count": status_counts["Review"],
        "Insufficient_Count": status_counts["Insufficient"],
        "Failed_Count": status_counts["Failed"],
        "Envelope_QC_Eligible_Count": sum(1 for row in qc_rows if row.get("Envelope_QC_Eligible")),
        "Intact_Strict_Eligible_Count": len(strict_rows),
        "Intact_Review_Eligible_Count": len(review_rows),
        "Comparison_Ready_Strict_Count": sum(1 for row in qc_rows if row.get("Comparison_Ready_Strict")),
        "Comparison_Ready_Review_Count": sum(1 for row in qc_rows if row.get("Comparison_Ready_Review")),
        "Comparison_Ready_Count": sum(1 for row in qc_rows if row.get("Comparison_Ready")),
        "Exact_Duplicate_Group_Count": len(exact_duplicate_group_ids),
        "Exact_Duplicate_Candidate_Count": duplicate_candidate_count,
        "Intact_Envelope_Group_Count": len(group_ids),
        "Unique_Envelope_Group_Count": ambiguity_counts.get("unique", 0),
        "Overlapping_Envelope_Group_Count": ambiguity_counts.get("overlapping_envelopes", 0) + ambiguity_counts.get("exact_duplicates", 0),
        "Competing_Reconstruction_Group_Count": ambiguity_counts.get("competing_reconstructions", 0) + ambiguity_counts.get("unresolved_group", 0),
        "Comparison_Representative_Count": len(comparison_reps),
        "Target_Review_Representative_Count": len(target_review_reps),
        "Candidates_Removed_As_Exact_Duplicates": sum(1 for row in qc_rows if not row.get("Is_Exact_Duplicate_Representative")),
        "Candidates_Removed_As_Group_Nonrepresentatives": sum(1 for row in qc_rows if not row.get("Group_Representative")),
        "Dominant_Comparison_Representative_Mass": dominant_comparison.get("Reconstructed_Mass"),
        "Dominant_Comparison_Representative_Intensity": dominant_comparison.get("Total_Supporting_Intensity"),
        "Dominant_Target_Review_Representative_Mass": dominant_target.get("Reconstructed_Mass"),
        "Dominant_Target_Review_Representative_Intensity": dominant_target.get("Total_Supporting_Intensity"),
        "Grouping_Mass_Tolerance_Da": grouping_config["mass_tolerance_da"],
        "Grouping_RT_Tolerance_Min": grouping_config["rt_tolerance_min"],
        "Grouping_Min_Shared_Peak_Fraction": grouping_config["min_shared_peak_fraction"],
        "Grouping_Min_Shared_Charge_Fraction": grouping_config["min_shared_charge_fraction"],
        "Trace_Only_Envelope_Count": sum(1 for row in qc_rows if row.get("Trace_Only_Envelope")),
        "Noncontiguous_Envelope_Count": sum(1 for row in qc_rows if row.get("Charge_State_Continuity") == "non_contiguous"),
        "RT_Inconsistent_Count": sum(1 for row in qc_rows if row.get("RT_Consistency") == "inconsistent"),
        "Internal_Mass_Error_Count": sum(1 for row in qc_rows if float(row.get("Envelope_Internal_Error_ppm") or 0.0) > qc_config["max_envelope_internal_error_ppm"]),
        "Theory_Near_Match_Count": sum(1 for row in qc_rows if row.get("Unmodified_Theory_Delta_ppm") is not None and abs(float(row.get("Unmodified_Theory_Delta_ppm") or 0.0)) <= qc_config["max_mass_error_ppm"]),
        "Reference_Match_Count": sum(1 for row in qc_rows if row.get("Reference_Mass_Matched")),
        "Dominant_Envelope_Mass": dominant.get("Reconstructed_Mass"),
        "Dominant_Envelope_Intensity": dominant.get("Total_Supporting_Intensity"),
        "Dominant_Envelope_Status": dominant.get("Reconstruction_Status"),
        "Dominant_Envelope_Comparison_Ready": dominant.get("Comparison_Ready"),
        "Dominant_Envelope_Overall_Mass": dominant.get("Reconstructed_Mass"),
        "Dominant_Envelope_Overall_Intensity": dominant.get("Total_Supporting_Intensity"),
        "Dominant_Envelope_In_Mass_Range_Mass": dominant_in_range.get("Reconstructed_Mass"),
        "Dominant_Envelope_In_Mass_Range_Intensity": dominant_in_range.get("Total_Supporting_Intensity"),
        "Dominant_Envelope_In_Mass_Range_Status": dominant_in_range.get("Reconstruction_Status"),
        "Dominant_Envelope_In_Mass_Range_Comparison_Ready": dominant_in_range.get("Comparison_Ready"),
        "Dominant_Envelope_In_Search_Range_Raw_Mass": dominant_in_range.get("Reconstructed_Mass"),
        "Dominant_Envelope_In_Search_Range_Raw_Intensity": dominant_in_range.get("Total_Supporting_Intensity"),
        "Dominant_Intact_Strict_Envelope_Mass": dominant_strict.get("Reconstructed_Mass"),
        "Dominant_Intact_Strict_Envelope_Intensity": dominant_strict.get("Total_Supporting_Intensity"),
        "Dominant_Intact_Strict_QC_Score": dominant_strict.get("Intact_Envelope_QC_Score"),
        "Dominant_Intact_Review_Envelope_Mass": dominant_review.get("Reconstructed_Mass"),
        "Dominant_Intact_Review_Envelope_Intensity": dominant_review.get("Total_Supporting_Intensity"),
        "Dominant_Intact_Review_QC_Score": dominant_review.get("Intact_Envelope_QC_Score"),
        "Dominant_Intact_Eligible_Envelope_Mass": dominant_eligible.get("Reconstructed_Mass"),
        "Dominant_Intact_Eligible_Envelope_Intensity": dominant_eligible.get("Total_Supporting_Intensity"),
        "Dominant_Intact_Eligible_QC_Score": dominant_eligible.get("Intact_Envelope_QC_Score"),
        "Dominant_Intact_Eligible_Reference_Label": dominant_eligible.get("Best_Reference_Label"),
        "Failure_Reason_Counts": reason_summary,
        "Reconstruction_Enabled": reconstruction_enabled,
        "Reconstruction_Engine": engine_stats.get("Reconstruction_Engine", qc_config.get("engine", "legacy_cluster")),
        "Num_RT_Windows": engine_stats.get("Num_RT_Windows", 0),
        "Num_Local_Peaks": engine_stats.get("Num_Local_Peaks", 0),
        "Num_Anchor_Peaks_Evaluated": engine_stats.get("Num_Anchor_Peaks_Evaluated", 0),
        "Num_Raw_Envelope_Candidates": engine_stats.get("Num_Raw_Envelope_Candidates", len(qc_rows)),
        "Num_Candidates_After_Charge_Filter": engine_stats.get("Num_Candidates_After_Charge_Filter", len(qc_rows)),
        "Num_Candidates_After_RT_Window_Merge": engine_stats.get("Num_Candidates_After_RT_Window_Merge", len(qc_rows)),
        "Num_Candidates_With_Consecutive_Charges": engine_stats.get("Num_Candidates_With_Consecutive_Charges", sum(1 for row in qc_rows if row.get("Longest_Consecutive_Charge_Run", 0) >= 2)),
        "Num_Candidates_With_Charge_Gaps": engine_stats.get("Num_Candidates_With_Charge_Gaps", sum(1 for row in qc_rows if row.get("Charge_Gap_Count", 0))),
        "Num_Missing_Charges_Evaluated": engine_stats.get("Num_Missing_Charges_Evaluated", 0),
        "Num_Missing_Charges_With_Weak_Peaks": engine_stats.get("Num_Missing_Charges_With_Weak_Peaks", 0),
        "Num_Missing_Charges_Not_Detected": engine_stats.get("Num_Missing_Charges_Not_Detected", 0),
        "Median_RT_Range_Min": engine_stats.get("Median_RT_Range_Min"),
        "Median_Internal_Error_ppm": engine_stats.get("Median_Internal_Error_ppm"),
        "Median_Charge_Count": engine_stats.get("Median_Charge_Count"),
        "Processing_Time_Seconds": engine_stats.get("Processing_Time_Seconds"),
        "Candidate_Count_By_Quality_Tier": tier_summary,
        "Strict_Failure_Reason_Counts": strict_failure_summary,
        "Review_Failure_Reason_Counts": review_failure_summary,
        "Charge_Count_Distribution": charge_count_summary,
        "Consecutive_Charge_Distribution": consecutive_summary,
        "Internal_Error_Distribution": internal_error_summary,
        "RT_Range_Distribution": rt_range_summary,
        "Peak_Usage_Distribution": peak_usage_summary,
        "Split_Envelope_Count": engine_stats.get("Split_Envelope_Count", sum(1 for row in qc_rows if row.get("Split_Envelope_Merged"))),
        "Missing_Charge_Status_Count": missing_status_summary,
        "Engine_Match_Status_Count": engine_match_summary,
        "Matched_Engine_Candidate_Count": engine_stats.get("Matched_Engine_Candidate_Count", 0),
        "Legacy_Only_Count": engine_stats.get("Legacy_Only_Count", 0),
        "RT_Localized_Only_Count": engine_stats.get("RT_Localized_Only_Count", 0),
        "Mass_Mismatch_Count": engine_stats.get("Mass_Mismatch_Count", 0),
        "RT_Mismatch_Count": engine_stats.get("RT_Mismatch_Count", 0),
        "Charge_Overlap_Failure_Count": engine_stats.get("Charge_Overlap_Failure_Count", 0),
        "Median_Mass_Delta_For_Matched_Only": engine_stats.get("Median_Mass_Delta_For_Matched_Only"),
        "Competing_Envelope_Group_Count": competition_stats.get("Competing_Envelope_Group_Count", 0),
        "MultiCandidate_Competition_Group_Count": competition_stats.get("MultiCandidate_Competition_Group_Count", 0),
        "Noncompeting_Candidate_Count": competition_stats.get("Noncompeting_Candidate_Count", 0),
        "Largest_Competition_Group_Size": competition_stats.get("Largest_Competition_Group_Size", 0),
        "Median_Competition_Group_Size": competition_stats.get("Median_Competition_Group_Size"),
        "Mean_Competition_Group_Size": competition_stats.get("Mean_Competition_Group_Size"),
        "Competition_Group_Size_Distribution": competition_stats.get("Competition_Group_Size_Distribution", ""),
        "Median_Top_Evidence_Score": competition_stats.get("Median_Top_Evidence_Score"),
        "Median_Top_vs_Second_Score_Margin": competition_stats.get("Median_Top_vs_Second_Score_Margin"),
        "Competition_Groups_With_Close_Scores": competition_stats.get("Competition_Groups_With_Close_Scores", 0),
        "High_Peak_Sharing_Candidate_Count": competition_stats.get("High_Peak_Sharing_Candidate_Count", 0),
        "Competitive_Scoring_Time_Seconds": competition_stats.get("Competitive_Scoring_Time_Seconds", 0.0),
        "Dry_Run_Selected_Candidate_Count": competition_stats.get("Dry_Run_Selected_Candidate_Count", 0),
        "Dry_Run_Primary_Selected_Count": competition_stats.get("Dry_Run_Primary_Selected_Count", 0),
        "Dry_Run_Independent_Selected_Count": competition_stats.get("Dry_Run_Independent_Selected_Count", 0),
        "Dry_Run_Would_Exclude_Count": competition_stats.get("Dry_Run_Would_Exclude_Count", 0),
        "Dry_Run_Ambiguous_Count": competition_stats.get("Dry_Run_Ambiguous_Count", 0),
        "Dry_Run_Noncompeting_Count": competition_stats.get("Dry_Run_Noncompeting_Count", 0),
        "Would_Exclude_Peak_Reuse_Count": competition_stats.get("Would_Exclude_Peak_Reuse_Count", 0),
        "Would_Exclude_Independent_Peak_Shortage_Count": competition_stats.get("Would_Exclude_Independent_Peak_Shortage_Count", 0),
        "Would_Exclude_Independent_Charge_Shortage_Count": competition_stats.get("Would_Exclude_Independent_Charge_Shortage_Count", 0),
        "Median_Independent_Peak_Fraction": competition_stats.get("Median_Independent_Peak_Fraction"),
        "Median_Direct_Competitor_Count": competition_stats.get("Median_Direct_Competitor_Count"),
        "Largest_Component_Selected_Count": competition_stats.get("Largest_Component_Selected_Count", 0),
        "Largest_Component_Excluded_Count": competition_stats.get("Largest_Component_Excluded_Count", 0),
        "Largest_Component_Ambiguous_Count": competition_stats.get("Largest_Component_Ambiguous_Count", 0),
        "Assignment_Dry_Run_Time_Seconds": competition_stats.get("Assignment_Dry_Run_Time_Seconds", 0.0),
        "Neutral_Mass_Search_Min_Da": qc_config["neutral_mass_range"]["min_da"],
        "Neutral_Mass_Search_Max_Da": qc_config["neutral_mass_range"]["max_da"],
        "Total_Candidates_Before_Mass_Range_Filter": len(qc_rows),
        "Total_Candidates_In_Mass_Range": len(in_range_rows),
        "Total_Candidates_Outside_Mass_Range": len(outside_rows),
        "Target_Review_Mass_Range_Settings": _target_review_settings(qc_config),
        "Target_Review_Candidate_Count": len(target_review_rows),
        "Search_Mode": qc_config["search_mode"],
        "Intensity_Normalization_Method": "overall, in-range raw, and intact-eligible relative intensity are reported separately",
        "RT_Tolerance_Settings": rt_settings,
        "Reference_Masses_Used": reference_summary,
        "Min_Charge_States_For_Reliable": qc_config["min_charge_states_for_reliable"],
        "Min_Charge_States_For_Review": qc_config["min_charge_states_for_review"],
        "Require_Contiguous_Charge_States": qc_config["require_contiguous_charge_states"],
        "Max_Neutral_Mass_SD_Da": qc_config["max_neutral_mass_sd_da"],
        "Max_Neutral_Mass_Range_Da": qc_config["max_neutral_mass_range_da"],
        "Max_Envelope_Internal_Error_ppm": qc_config["max_envelope_internal_error_ppm"],
        "Min_Relative_Intensity_Percent": qc_config["min_relative_intensity_percent"],
        "Min_Relative_Envelope_Intensity_Percent_For_Reliable": qc_config["min_relative_envelope_intensity_percent_for_reliable"],
        "Min_Relative_Envelope_Intensity_Percent_For_Review": qc_config["min_relative_envelope_intensity_percent_for_review"],
        "Max_Competing_Envelopes": qc_config["max_competing_envelopes"],
        "Comparison_Ready_Statuses": "; ".join(map(str, qc_config["comparison_ready_statuses"])),
        "Notes": "Reliable emphasizes charge-envelope internal quality, RT consistency, and signal support. Comparison_Ready requires intact eligibility, not only neutral mass range membership. Envelope grouping collapses exact duplicate and overlapping reconstruction candidates before comparison representative export. Reference masses and target review range are annotations only and do not affect global grouping or representative selection.",
    }]





def build_rt_engine_qc_summary_rows(diagnostic_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not diagnostic_rows:
        return []
    diagnostic = diagnostic_rows[0]
    metrics = [
        "Candidate_Count_By_Quality_Tier",
        "Strict_Failure_Reason_Counts",
        "Review_Failure_Reason_Counts",
        "Charge_Count_Distribution",
        "Consecutive_Charge_Distribution",
        "Internal_Error_Distribution",
        "RT_Range_Distribution",
        "Peak_Usage_Distribution",
        "Split_Envelope_Count",
        "Missing_Charge_Status_Count",
        "Engine_Match_Status_Count",
        "Matched_Engine_Candidate_Count",
        "Legacy_Only_Count",
        "RT_Localized_Only_Count",
        "Mass_Mismatch_Count",
        "RT_Mismatch_Count",
        "Charge_Overlap_Failure_Count",
        "Median_Mass_Delta_For_Matched_Only",
        "Competing_Envelope_Group_Count",
        "MultiCandidate_Competition_Group_Count",
        "Noncompeting_Candidate_Count",
        "Largest_Competition_Group_Size",
        "Median_Competition_Group_Size",
        "Mean_Competition_Group_Size",
        "Competition_Group_Size_Distribution",
        "Median_Top_Evidence_Score",
        "Median_Top_vs_Second_Score_Margin",
        "Competition_Groups_With_Close_Scores",
        "High_Peak_Sharing_Candidate_Count",
        "Competitive_Scoring_Time_Seconds",
        "Dry_Run_Selected_Candidate_Count",
        "Dry_Run_Primary_Selected_Count",
        "Dry_Run_Independent_Selected_Count",
        "Dry_Run_Would_Exclude_Count",
        "Dry_Run_Ambiguous_Count",
        "Dry_Run_Noncompeting_Count",
        "Would_Exclude_Peak_Reuse_Count",
        "Would_Exclude_Independent_Peak_Shortage_Count",
        "Would_Exclude_Independent_Charge_Shortage_Count",
        "Median_Independent_Peak_Fraction",
        "Median_Direct_Competitor_Count",
        "Largest_Component_Selected_Count",
        "Largest_Component_Excluded_Count",
        "Largest_Component_Ambiguous_Count",
        "Assignment_Dry_Run_Time_Seconds",
        "Processing_Time_Seconds",
    ]
    return [{"Metric": metric, "Value": diagnostic.get(metric), "Notes": ""} for metric in metrics]


def _ppm_error(observed: float, expected: float) -> float:
    if not expected:
        return 0.0
    return (float(observed) - float(expected)) / float(expected) * 1_000_000


def _within_ppm(observed: float, expected: float, tolerance_ppm: float) -> bool:
    return abs(_ppm_error(observed, expected)) <= tolerance_ppm


def _longest_consecutive_run(charges: list[int]) -> int:
    if not charges:
        return 0
    ordered = sorted(set(charges))
    best = current = 1
    for left, right in zip(ordered, ordered[1:]):
        if right == left + 1:
            current += 1
            best = max(best, current)
        else:
            current = 1
    return best


def _charge_gap_count(charges: list[int]) -> int:
    if len(set(charges)) < 2:
        return 0
    ordered = sorted(set(charges))
    return max(0, (ordered[-1] - ordered[0] + 1) - len(ordered))


def _charge_continuity_fraction_value(charges: list[int]) -> float:
    if not charges:
        return 0.0
    ordered = sorted(set(charges))
    span = ordered[-1] - ordered[0] + 1
    return len(ordered) / span if span else 0.0


def _mean(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def _weighted_mean(values: list[float], weights: list[float]) -> float | None:
    denominator = sum(weights)
    if not values or denominator <= 0:
        return _mean(values)
    return sum(value * weight for value, weight in zip(values, weights)) / denominator


def _internal_error_stats(values: list[float], center: float) -> tuple[float | None, float | None, float | None]:
    if not values or not center:
        return None, None, None
    errors = [abs(value - center) / center * 1_000_000 for value in values]
    return max(errors), _mean(errors), median(errors)


def _candidate_extra(candidate: IntactMassCandidate, field: str, default: Any = None) -> Any:
    return getattr(candidate, field, default)


def _set_candidate_extra(candidate: IntactMassCandidate, **values: Any) -> None:
    for key, value in values.items():
        setattr(candidate, key, value)


def _rt_localized_config(reconstruction_config: dict[str, Any]) -> dict[str, Any]:
    return _qc_config(reconstruction_config or {})["rt_localized"]


def _peak_rows_from_tier_result(tier_result: PeakTierResult, include_below: bool = False) -> list[dict[str, Any]]:
    peaks = list(tier_result.usable_peaks)
    if include_below:
        peaks += list(tier_result.below_threshold)
    rows = []
    for index, peak in enumerate(peaks, start=1):
        rows.append({
            "Peak_ID": f"{peak.scan_id or 'scan_na'}|rt={_format_float(peak.rt, 5)}|mz={_format_float(peak.mz, 6)}|i={index}",
            "mz": float(peak.mz),
            "Intensity": float(peak.intensity),
            "RT": peak.rt,
            "Scan_ID": peak.scan_id,
            "Peak_Tier": peak.tier or "",
        })
    return rows


def _rt_windows(peak_rows: list[dict[str, Any]], rt_config: dict[str, Any]) -> list[dict[str, Any]]:
    rts = sorted({float(row["RT"]) for row in peak_rows if row.get("RT") is not None})
    if not rts:
        return [{"RT_Window_ID": "RT00001", "start": None, "end": None, "center": None, "rows": peak_rows}]
    width = float(rt_config["rt_window_min"])
    step = float(rt_config["rt_step_min"])
    start = min(rts)
    stop = max(rts)
    windows = []
    index = 1
    current = start
    epsilon = 1e-9
    while current <= stop + epsilon:
        end = current + width
        rows = [row for row in peak_rows if row.get("RT") is not None and current - epsilon <= float(row["RT"]) <= end + epsilon]
        scans = {row.get("Scan_ID") or row.get("RT") for row in rows}
        if rows and len(scans) >= int(rt_config["min_scans_per_window"]):
            windows.append({
                "RT_Window_ID": f"RT{index:05d}",
                "start": current,
                "end": end,
                "center": current + width / 2.0,
                "rows": rows,
            })
            index += 1
        current += step
    return windows


def _aggregate_peak_group(group: list[dict[str, Any]], window: dict[str, Any], rt_config: dict[str, Any]) -> dict[str, Any]:
    method = rt_config["peak_aggregation"]
    intensities = [float(row.get("Intensity") or 0.0) for row in group]
    if method == "sum":
        intensity = sum(intensities)
    elif method == "mean":
        intensity = sum(intensities) / len(intensities)
    else:
        intensity = max(intensities)
    mz = _weighted_mean([float(row["mz"]) for row in group], intensities) or float(group[0]["mz"])
    rts = [float(row["RT"]) for row in group if row.get("RT") is not None]
    return {
        "Local_Peak_ID": "",
        "mz": mz,
        "Intensity": intensity,
        "Max_Intensity": max(intensities) if intensities else intensity,
        "RT": sum(rts) / len(rts) if rts else window.get("center"),
        "Scan_ID": ";".join(sorted({str(row.get("Scan_ID")) for row in group if row.get("Scan_ID") not in {None, ""}})),
        "Peak_Tier": "; ".join(sorted({str(row.get("Peak_Tier") or "") for row in group if row.get("Peak_Tier")})),
        "Source_Peak_IDs": sorted({str(row.get("Peak_ID")) for row in group}),
        "Source_RT_Values": sorted({float(row["RT"]) for row in group if row.get("RT") is not None}),
        "Num_Contributing_Scans": len({row.get("Scan_ID") or row.get("RT") for row in group}),
        "RT_Window_ID": window["RT_Window_ID"],
        "RT_Window_Start_Min": window.get("start"),
        "RT_Window_End_Min": window.get("end"),
        "RT_Window_Center_Min": window.get("center"),
    }


def _local_peak_table(window: dict[str, Any], rt_config: dict[str, Any]) -> list[dict[str, Any]]:
    rows = sorted(window.get("rows") or [], key=lambda row: float(row["mz"]))
    if not rows:
        return []
    groups: list[list[dict[str, Any]]] = []
    for row in rows:
        if not groups:
            groups.append([row])
            continue
        current = _weighted_mean([float(item["mz"]) for item in groups[-1]], [float(item.get("Intensity") or 0.0) for item in groups[-1]]) or float(groups[-1][0]["mz"])
        if _within_ppm(float(row["mz"]), current, rt_config["mz_merge_tolerance_ppm"]):
            groups[-1].append(row)
        else:
            groups.append([row])
    local_peaks = []
    for index, group in enumerate(groups, start=1):
        local = _aggregate_peak_group(group, window, rt_config)
        local["Local_Peak_ID"] = f"{window['RT_Window_ID']}_LP{index:05d}"
        local_peaks.append(local)
    return local_peaks


def _nearest_peak(local_peaks: list[dict[str, Any]], mz_values: list[float], target_mz: float, tolerance_ppm: float) -> tuple[dict[str, Any] | None, float | None]:
    if not local_peaks:
        return None, None
    position = bisect_left(mz_values, target_mz)
    best = None
    best_error = None
    for index in (position - 1, position, position + 1):
        if 0 <= index < len(local_peaks):
            candidate = local_peaks[index]
            error = _ppm_error(float(candidate["mz"]), target_mz)
            if best is None or abs(error) < abs(best_error):
                best = candidate
                best_error = error
    if best is not None and abs(best_error) <= tolerance_ppm:
        return best, best_error
    return None, best_error


def _missing_charge_rows(
    cluster_id: str,
    window_id: str,
    reconstructed_mass: float,
    observed_charges: list[int],
    local_peaks_all: list[dict[str, Any]],
    mz_values_all: list[float],
    rt_config: dict[str, Any],
    instrument_config: dict[str, Any],
    usable_peak_ids: set[str],
    usable_source_peak_ids: set[str] | None = None,
) -> list[dict[str, Any]]:
    if len(observed_charges) < 2:
        return []
    rows = []
    polarity = instrument_config.get("polarity", "negative")
    for charge in range(min(observed_charges), max(observed_charges) + 1):
        if charge in observed_charges:
            continue
        predicted_mz = mz_from_neutral_mass(reconstructed_mass, charge, polarity)
        nearest, error = _nearest_peak(local_peaks_all, mz_values_all, predicted_mz, rt_config["adjacent_charge_mz_tolerance_ppm"])
        if nearest is None:
            status = "no_peak_in_tolerance"
            nearest_mz = None
            nearest_intensity = None
        elif usable_source_peak_ids and set(nearest.get("Source_Peak_IDs", [])) & usable_source_peak_ids:
            status = "detected"
            nearest_mz = nearest["mz"]
            nearest_intensity = nearest["Intensity"]
        else:
            status = "below_intensity_threshold"
            nearest_mz = nearest["mz"]
            nearest_intensity = nearest["Intensity"]
        rows.append({
            "Cluster_ID": cluster_id,
            "RT_Window_ID": window_id,
            "Reconstructed_Mass": reconstructed_mass,
            "Missing_Charge": charge,
            "Predicted_MZ": predicted_mz,
            "Nearest_Observed_MZ": nearest_mz,
            "Error_ppm": error if nearest is not None else None,
            "Nearest_Intensity": nearest_intensity,
            "Detection_Status": status,
            "Notes": "",
        })
    return rows



def _extension_charge_rows(
    reconstructed_mass: float,
    observed_charges: list[int],
    local_peaks_all: list[dict[str, Any]],
    mz_values_all: list[float],
    rt_config: dict[str, Any],
    instrument_config: dict[str, Any],
    min_charge: int,
    max_charge: int,
    local_max: float,
    usable_source_peak_ids: set[str],
) -> dict[str, Any]:
    extension_cfg = rt_config["charge_extension"]
    if not extension_cfg.get("enabled", True) or not observed_charges:
        return {
            "lower": [], "upper": [], "detected": [], "weak": [], "not_detected": [], "improved": False, "peaks": {}
        }
    polarity = instrument_config.get("polarity", "negative")
    lower = [charge for charge in range(max(min_charge, min(observed_charges) - extension_cfg["max_extension_charges"]), min(observed_charges))]
    upper = [charge for charge in range(max(observed_charges) + 1, min(max_charge, max(observed_charges) + extension_cfg["max_extension_charges"]) + 1)]
    detected = []
    weak = []
    not_detected = []
    peaks: dict[int, dict[str, Any]] = {}
    for charge in lower + upper:
        predicted_mz = mz_from_neutral_mass(reconstructed_mass, charge, polarity)
        peak, _ = _nearest_peak(local_peaks_all, mz_values_all, predicted_mz, extension_cfg["weak_peak_tolerance_ppm"])
        if peak is None:
            not_detected.append(charge)
            continue
        local_rel = (float(peak.get("Intensity") or 0.0) / local_max * 100.0) if local_max else 0.0
        is_usable = bool(set(peak.get("Source_Peak_IDs", [])) & usable_source_peak_ids)
        if is_usable:
            detected.append(charge)
            peaks[charge] = peak
        elif local_rel >= extension_cfg["weak_peak_min_local_relative_percent"]:
            weak.append(charge)
            peaks[charge] = peak
        else:
            not_detected.append(charge)
    return {
        "lower": lower,
        "upper": upper,
        "detected": detected,
        "weak": weak,
        "not_detected": not_detected,
        "improved": bool(detected or weak),
        "peaks": peaks,
    }


def _recalculate_record(record: dict[str, Any], rt_config: dict[str, Any], polarity: str) -> None:
    charges = sorted(record["observed"])
    record["charges"] = charges
    record["neutral_masses"] = [neutral_mass_from_mz(record["observed"][z]["mz"], z, polarity) for z in charges]
    record["intensities"] = [float(record["observed"][z].get("Intensity") or 0.0) for z in charges]
    record["unweighted"] = _mean(record["neutral_masses"])
    record["weighted"] = _weighted_mean(record["neutral_masses"], record["intensities"])
    record["median"] = median(record["neutral_masses"])
    estimator = rt_config["neutral_mass_estimator"]
    record["mass"] = record["weighted"] if estimator == "intensity_weighted_mean" else record["median"] if estimator == "median" else record["unweighted"]
    record["internal_max"], record["internal_mean"], record["internal_median"] = _internal_error_stats(record["neutral_masses"], record["mass"])
    record["source_peak_ids"] = sorted({peak_id for peak in record["observed"].values() for peak_id in peak.get("Source_Peak_IDs", [])})
    record["local_peak_ids"] = sorted({peak["Local_Peak_ID"] for peak in record["observed"].values()})
    record["total_intensity"] = sum(float(record["observed"][z].get("Intensity") or 0.0) for z in charges)
    record["gap_count"] = _charge_gap_count(charges)
    record["longest_run"] = _longest_consecutive_run(charges)
    record["continuity_fraction"] = _charge_continuity_fraction_value(charges)
    record["local_rel_envelope"] = (record["total_intensity"] / record["local_max"] * 100.0) if record.get("local_max") else 0.0


def _apply_split_envelope_merge(records: list[dict[str, Any]], rt_config: dict[str, Any], polarity: str) -> list[dict[str, Any]]:
    split_cfg = rt_config["split_envelope_merge"]
    if not split_cfg.get("enabled", True) or len(records) < 2:
        return records
    merged: list[dict[str, Any]] = []
    group_index = 1
    for record in sorted(records, key=lambda item: (item["mass"], item["window"].get("center") or 0.0)):
        target = None
        for existing in merged:
            mass_ok = abs(_ppm_error(record["mass"], existing["mass"])) <= split_cfg["mass_tolerance_ppm"]
            rt_delta = abs(float(record["window"].get("center") or 0.0) - float(existing["window"].get("center") or 0.0))
            rt_ok = rt_delta <= split_cfg["rt_tolerance_min"]
            combined = sorted(set(record["charges"]) | set(existing["charges"]))
            gap_before = _charge_gap_count(record["charges"]) + _charge_gap_count(existing["charges"])
            gap_after = _charge_gap_count(combined)
            connectable = gap_after <= split_cfg["max_charge_gap"] and combined != sorted(existing["charges"])
            if mass_ok and rt_ok and connectable:
                target = existing
                target["charge_gaps_before_merge"] = gap_before
                target["charge_gaps_after_merge"] = gap_after
                break
        if target is None:
            record["split_envelope_group_id"] = f"SEG{group_index:05d}"
            record["split_envelope_member_count"] = 1
            record["split_envelope_merged"] = False
            record["charge_gaps_before_merge"] = _charge_gap_count(record["charges"])
            record["charge_gaps_after_merge"] = _charge_gap_count(record["charges"])
            merged.append(record)
            group_index += 1
        else:
            for charge, peak in record["observed"].items():
                if charge not in target["observed"] or float(peak.get("Intensity") or 0.0) > float(target["observed"][charge].get("Intensity") or 0.0):
                    target["observed"][charge] = peak
            target["rt_window_ids"] = sorted(set(target["rt_window_ids"] + record["rt_window_ids"]))
            target["split_envelope_member_count"] = int(target.get("split_envelope_member_count", 1)) + 1
            target["split_envelope_merged"] = True
            _recalculate_record(target, rt_config, polarity)
    return merged

def _build_rt_localized_candidates(
    tier_result: PeakTierResult,
    reconstruction_config: dict[str, Any],
    instrument_config: dict[str, Any],
    theoretical_mass: float | None,
) -> tuple[list[IntactMassCandidate], list[dict[str, Any]], dict[str, Any]]:
    started = perf_counter()
    qc_config = _qc_config(reconstruction_config or {})
    rt_config = qc_config["rt_localized"]
    min_charge = int(reconstruction_config.get("min_charge", 5))
    max_charge = int(reconstruction_config.get("max_charge", 40))
    polarity = instrument_config.get("polarity", "negative")
    usable_rows = _peak_rows_from_tier_result(tier_result, include_below=False)
    all_rows = _peak_rows_from_tier_result(tier_result, include_below=True)
    windows = _rt_windows(usable_rows, rt_config)
    candidates: list[IntactMassCandidate] = []
    charge_state_peaks: list[dict[str, Any]] = []
    rt_diagnostics: list[dict[str, Any]] = []
    missing_diagnostics: list[dict[str, Any]] = []
    raw_candidate_count = 0
    charge_filtered_count = 0
    anchor_count = 0
    local_peak_count = 0
    local_peak_usage: dict[str, int] = {}
    raw_candidate_records: list[dict[str, Any]] = []
    for window in windows:
        local_peaks = _local_peak_table(window, rt_config)
        if window.get("start") is None:
            all_window_rows = all_rows
        else:
            all_window_rows = [
                row for row in all_rows
                if row.get("RT") is not None and float(window["start"]) <= float(row["RT"]) <= float(window["end"])
            ]
        all_window = {**window, "rows": all_window_rows}
        local_peaks_all = _local_peak_table(all_window, rt_config)
        local_peak_count += len(local_peaks)
        local_peaks = sorted(local_peaks, key=lambda row: row["mz"])
        local_peaks_all = sorted(local_peaks_all, key=lambda row: row["mz"])
        mz_values = [float(row["mz"]) for row in local_peaks]
        mz_values_all = [float(row["mz"]) for row in local_peaks_all]
        local_max = max((float(row["Intensity"]) for row in local_peaks), default=0.0)
        usable_peak_ids = {row["Local_Peak_ID"] for row in local_peaks}
        usable_source_peak_ids = {peak_id for row in local_peaks for peak_id in row.get("Source_Peak_IDs", [])}
        for anchor in local_peaks:
            local_rel_peak = (float(anchor["Intensity"]) / local_max * 100.0) if local_max else 0.0
            if local_rel_peak < rt_config["min_local_relative_peak_intensity_percent"]:
                continue
            anchor_count += 1
            for charge in range(min_charge, max_charge + 1):
                anchor_mass = neutral_mass_from_mz(anchor["mz"], charge, polarity)
                in_range, _, _, _ = _neutral_mass_range_status(anchor_mass, qc_config)
                if not in_range:
                    continue
                predicted_charges = sorted(set(range(max(min_charge, charge - rt_config["max_charge_gap"]), min(max_charge, charge + rt_config["max_charge_gap"]) + 1)))
                observed: dict[int, dict[str, Any]] = {charge: anchor}
                for predicted_charge in predicted_charges:
                    if predicted_charge == charge:
                        continue
                    predicted_mz = mz_from_neutral_mass(anchor_mass, predicted_charge, polarity)
                    peak, _ = _nearest_peak(local_peaks, mz_values, predicted_mz, rt_config["adjacent_charge_mz_tolerance_ppm"])
                    if peak is not None:
                        observed[predicted_charge] = peak
                observed_charges = sorted(observed)
                raw_candidate_count += 1
                longest_run = _longest_consecutive_run(observed_charges)
                if len(observed_charges) < rt_config["min_charge_states"]:
                    continue
                if rt_config["require_consecutive_for_candidate"] and longest_run < rt_config["min_consecutive_charge_states"]:
                    continue
                charge_filtered_count += 1
                neutral_masses = [neutral_mass_from_mz(observed[z]["mz"], z, polarity) for z in observed_charges]
                intensities = [float(observed[z]["Intensity"]) for z in observed_charges]
                unweighted = _mean(neutral_masses)
                weighted = _weighted_mean(neutral_masses, intensities)
                med = median(neutral_masses)
                estimator = rt_config["neutral_mass_estimator"]
                reconstructed_mass = weighted if estimator == "intensity_weighted_mean" else med if estimator == "median" else unweighted
                if reconstructed_mass is None:
                    continue
                internal_max, internal_mean, internal_median = _internal_error_stats(neutral_masses, reconstructed_mass)
                source_ids = sorted({peak_id for peak in observed.values() for peak_id in peak.get("Source_Peak_IDs", [])})
                local_ids = sorted({peak["Local_Peak_ID"] for peak in observed.values()})
                total_intensity = sum(float(observed[z]["Intensity"]) for z in observed_charges)
                envelope_local_rel = (total_intensity / local_max * 100.0) if local_max else 0.0
                cluster_id = f"RTL_RAW_{len(raw_candidate_records) + 1:06d}"
                gap_count = _charge_gap_count(observed_charges)
                extension = _extension_charge_rows(reconstructed_mass, observed_charges, local_peaks_all, mz_values_all, rt_config, instrument_config, min_charge, max_charge, local_max, usable_source_peak_ids)
                original_charge_states = sorted(observed_charges)
                if rt_config["charge_extension"].get("add_weak_peaks_to_envelope", False):
                    for ext_charge, ext_peak in extension["peaks"].items():
                        observed.setdefault(ext_charge, ext_peak)
                    observed_charges = sorted(observed)
                    neutral_masses = [neutral_mass_from_mz(observed[z]["mz"], z, polarity) for z in observed_charges]
                    intensities = [float(observed[z]["Intensity"]) for z in observed_charges]
                    unweighted = _mean(neutral_masses)
                    weighted = _weighted_mean(neutral_masses, intensities)
                    med = median(neutral_masses)
                    reconstructed_mass = weighted if estimator == "intensity_weighted_mean" else med if estimator == "median" else unweighted
                    internal_max, internal_mean, internal_median = _internal_error_stats(neutral_masses, reconstructed_mass)
                    source_ids = sorted({peak_id for peak in observed.values() for peak_id in peak.get("Source_Peak_IDs", [])})
                    local_ids = sorted({peak["Local_Peak_ID"] for peak in observed.values()})
                    total_intensity = sum(float(observed[z]["Intensity"]) for z in observed_charges)
                    envelope_local_rel = (total_intensity / local_max * 100.0) if local_max else 0.0
                    gap_count = _charge_gap_count(observed_charges)
                missing_rows = _missing_charge_rows(cluster_id, window["RT_Window_ID"], reconstructed_mass, observed_charges, local_peaks_all, mz_values_all, rt_config, instrument_config, usable_peak_ids, usable_source_peak_ids)
                record = {
                    "cluster_id": cluster_id,
                    "mass": reconstructed_mass,
                    "charges": observed_charges,
                    "observed": observed,
                    "neutral_masses": neutral_masses,
                    "intensities": intensities,
                    "source_peak_ids": source_ids,
                    "local_peak_ids": local_ids,
                    "total_intensity": total_intensity,
                    "rt_window_ids": [window["RT_Window_ID"]],
                    "window": window,
                    "local_max": local_max,
                    "local_rel_peak": local_rel_peak,
                    "local_rel_envelope": envelope_local_rel,
                    "anchor_mz": anchor["mz"],
                    "anchor_charge": charge,
                    "predicted_charges": predicted_charges,
                    "missing_rows": missing_rows,
                    "unweighted": unweighted,
                    "weighted": weighted,
                    "median": med,
                    "internal_max": internal_max,
                    "internal_mean": internal_mean,
                    "internal_median": internal_median,
                    "gap_count": gap_count,
                    "longest_run": longest_run,
                    "continuity_fraction": _charge_continuity_fraction_value(observed_charges),
                    "extension": extension,
                    "original_charge_states": original_charge_states,
                    "final_charge_states": sorted(observed_charges),
                    "charge_extension_improved": extension["improved"],
                }
                raw_candidate_records.append(record)
                for local_id in local_ids:
                    local_peak_usage[local_id] = local_peak_usage.get(local_id, 0) + 1
    merge_cfg = rt_config["merge_across_windows"]
    merged_records: list[dict[str, Any]] = []
    for record in sorted(raw_candidate_records, key=lambda item: (item["mass"], item["window"].get("center") or 0.0)):
        target = None
        if merge_cfg.get("enabled", True):
            for existing in merged_records:
                mass_ppm = abs(_ppm_error(record["mass"], existing["mass"]))
                shared_charges = len(set(record["charges"]) & set(existing["charges"]))
                charge_fraction = shared_charges / max(min(len(record["charges"]), len(existing["charges"])), 1)
                rt_ok = True
                if merge_cfg.get("rt_overlap_required", True):
                    rt_ok = not (record["window"].get("start") is not None and existing["window"].get("end") is not None and (record["window"]["start"] > existing["window"]["end"] or existing["window"]["start"] > record["window"]["end"]))
                if mass_ppm <= merge_cfg["mass_tolerance_ppm"] and charge_fraction >= merge_cfg["min_shared_charge_fraction"] and rt_ok:
                    target = existing
                    break
        if target is None:
            merged_records.append(record)
        else:
            seen = set(target["source_peak_ids"])
            for peak_id in record["source_peak_ids"]:
                if peak_id not in seen:
                    seen.add(peak_id)
            target["source_peak_ids"] = sorted(seen)
            target["rt_window_ids"] = sorted(set(target["rt_window_ids"] + record["rt_window_ids"]))
            if record["total_intensity"] > target["total_intensity"]:
                keep_ids = target["rt_window_ids"]
                record["rt_window_ids"] = keep_ids
                target.update(record)
                target["rt_window_ids"] = keep_ids
            target["merged"] = True
    merged_records = _apply_split_envelope_merge(merged_records, rt_config, polarity)
    for index, record in enumerate(merged_records, start=1):
        cluster_id = f"RTL{index:05d}"
        record["cluster_id"] = cluster_id
        mass_error_da = record["mass"] - theoretical_mass if theoretical_mass is not None else None
        mass_error_ppm = (mass_error_da / theoretical_mass * 1_000_000) if theoretical_mass else None
        candidate = IntactMassCandidate(
            observed_mass=record["mass"],
            charge_state_count=len(record["charges"]),
            charge_states=record["charges"],
            supporting_peak_count=len(record["source_peak_ids"]),
            total_intensity=record["total_intensity"],
            theoretical_mass=theoretical_mass,
            mass_error_da=mass_error_da,
            mass_error_ppm=mass_error_ppm,
            confidence=_confidence(len(record["charges"]), int(reconstruction_config.get("min_charge_states", 3))),
            cluster_id=cluster_id,
        )
        missing_charge_states = [row["Missing_Charge"] for row in record["missing_rows"]]
        peak_usage_values = [local_peak_usage.get(local_id, 1) for local_id in record["local_peak_ids"]]
        max_peak_usage = max(peak_usage_values, default=1)
        mean_peak_usage = _mean(peak_usage_values) or 1.0
        highly_shared_count = sum(1 for value in peak_usage_values if value >= rt_config["peak_sharing"]["high_usage_threshold"])
        highly_shared_fraction = highly_shared_count / max(len(record["local_peak_ids"]), 1)
        if highly_shared_fraction > rt_config["peak_sharing"]["max_highly_shared_fraction_for_tier2"]:
            peak_sharing_status = "highly_shared"
        elif highly_shared_fraction > rt_config["peak_sharing"]["max_highly_shared_fraction_for_tier1"]:
            peak_sharing_status = "moderately_shared"
        else:
            peak_sharing_status = "low_shared"
        _set_candidate_extra(
            candidate,
            reconstruction_engine="rt_localized",
            rt_window_id=record["window"]["RT_Window_ID"],
            rt_window_start_min=record["window"].get("start"),
            rt_window_end_min=record["window"].get("end"),
            rt_window_center_min=record["window"].get("center"),
            num_ms1_scans_in_window=len({row.get("Scan_ID") or row.get("RT") for row in record["window"].get("rows", [])}),
            peak_aggregation_method=rt_config["peak_aggregation"],
            anchor_mz=record["anchor_mz"],
            anchor_charge=record["anchor_charge"],
            predicted_charge_states="; ".join(map(str, record["predicted_charges"])),
            observed_charge_states="; ".join(map(str, record["charges"])),
            missing_charge_states="; ".join(map(str, missing_charge_states)),
            missing_charge_predicted_mz="; ".join(_format_float(row["Predicted_MZ"], 6) for row in record["missing_rows"]),
            num_predicted_charges=len(record["predicted_charges"]),
            num_observed_charges=len(record["charges"]),
            charge_coverage_fraction=len(record["charges"]) / max(len(record["predicted_charges"]), 1),
            consecutive_charge_run_length=record["longest_run"],
            longest_consecutive_charge_run=record["longest_run"],
            charge_gap_count=record["gap_count"],
            charge_continuity_fraction=record["continuity_fraction"],
            peak_usage_count=max_peak_usage,
            shared_peak_count=sum(1 for local_id in record["local_peak_ids"] if local_peak_usage.get(local_id, 1) > 1),
            shared_peak_fraction=sum(1 for local_id in record["local_peak_ids"] if local_peak_usage.get(local_id, 1) > 1) / max(len(record["local_peak_ids"]), 1),
            max_peak_usage_count=max_peak_usage,
            mean_peak_usage_count=mean_peak_usage,
            num_highly_shared_peaks=highly_shared_count,
            highly_shared_peak_fraction=highly_shared_fraction,
            peak_sharing_status=peak_sharing_status,
            local_window_max_intensity=record["local_max"],
            local_relative_peak_intensity_percent=record["local_rel_peak"],
            local_envelope_relative_intensity_percent=record["local_rel_envelope"],
            neutral_mass_estimator=rt_config["neutral_mass_estimator"],
            neutral_mass_unweighted_mean=record["unweighted"],
            neutral_mass_weighted_mean=record["weighted"],
            neutral_mass_median=record["median"],
            envelope_internal_error_max_ppm=record["internal_max"],
            envelope_internal_error_mean_ppm=record["internal_mean"],
            envelope_internal_error_median_ppm=record["internal_median"],
            source_rt_window_ids="; ".join(record["rt_window_ids"]),
            num_source_rt_windows=len(record["rt_window_ids"]),
            merged_across_rt_windows=bool(record.get("merged") or len(record["rt_window_ids"]) > 1),
            extended_lower_charges_evaluated="; ".join(map(str, record.get("extension", {}).get("lower", []))),
            extended_upper_charges_evaluated="; ".join(map(str, record.get("extension", {}).get("upper", []))),
            extended_charges_detected="; ".join(map(str, record.get("extension", {}).get("detected", []))),
            extended_weak_charges_detected="; ".join(map(str, record.get("extension", {}).get("weak", []))),
            extended_charges_not_detected="; ".join(map(str, record.get("extension", {}).get("not_detected", []))),
            charge_extension_improved_envelope=bool(record.get("charge_extension_improved", False)),
            original_charge_states="; ".join(map(str, record.get("original_charge_states", record["charges"]))),
            final_charge_states="; ".join(map(str, record.get("final_charge_states", record["charges"]))),
            split_envelope_group_id=record.get("split_envelope_group_id", ""),
            split_envelope_member_count=record.get("split_envelope_member_count", 1),
            split_envelope_merged=record.get("split_envelope_merged", False),
            charge_gaps_before_merge=record.get("charge_gaps_before_merge", record["gap_count"]),
            charge_gaps_after_merge=record.get("charge_gaps_after_merge", record["gap_count"]),
        )
        candidates.append(candidate)
        for charge in record["charges"]:
            peak = record["observed"][charge]
            charge_state_peaks.append({
                "Cluster_ID": cluster_id,
                "mz": peak["mz"],
                "Intensity": peak["Intensity"],
                "RT": peak.get("RT"),
                "Scan_ID": peak.get("Scan_ID"),
                "Charge": charge,
                "Neutral_Mass": neutral_mass_from_mz(peak["mz"], charge, polarity),
                "Peak_Tier": peak.get("Peak_Tier"),
                "RT_Window_ID": record["window"]["RT_Window_ID"],
                "Local_Peak_ID": peak["Local_Peak_ID"],
            })
        for row in record["missing_rows"]:
            row = dict(row)
            row["Cluster_ID"] = cluster_id
            missing_diagnostics.append(row)
        rt_diagnostics.append({
            "Cluster_ID": cluster_id,
            "Reconstruction_Engine": "rt_localized",
            "RT_Window_ID": record["window"]["RT_Window_ID"],
            "RT_Window_Start_Min": record["window"].get("start"),
            "RT_Window_End_Min": record["window"].get("end"),
            "RT_Window_Center_Min": record["window"].get("center"),
            "Num_MS1_Scans_In_Window": len({row.get("Scan_ID") or row.get("RT") for row in record["window"].get("rows", [])}),
            "Peak_Aggregation_Method": rt_config["peak_aggregation"],
            "Anchor_MZ": record["anchor_mz"],
            "Anchor_Charge": record["anchor_charge"],
            "Reconstructed_Mass": record["mass"],
            "Predicted_Charge_States": "; ".join(map(str, record["predicted_charges"])),
            "Observed_Charge_States": "; ".join(map(str, record["charges"])),
            "Missing_Charge_States": "; ".join(map(str, missing_charge_states)),
            "Num_Predicted_Charges": len(record["predicted_charges"]),
            "Num_Observed_Charges": len(record["charges"]),
            "Charge_Coverage_Fraction": len(record["charges"]) / max(len(record["predicted_charges"]), 1),
            "Consecutive_Charge_Run_Length": record["longest_run"],
            "Longest_Consecutive_Charge_Run": record["longest_run"],
            "Charge_Gap_Count": record["gap_count"],
            "Charge_Continuity_Fraction": record["continuity_fraction"],
            "Local_Window_Max_Intensity": record["local_max"],
            "Local_Envelope_Relative_Intensity_Percent": record["local_rel_envelope"],
            "Neutral_Mass_Estimator": rt_config["neutral_mass_estimator"],
            "Neutral_Mass_Unweighted_Mean": record["unweighted"],
            "Neutral_Mass_Weighted_Mean": record["weighted"],
            "Neutral_Mass_Median": record["median"],
            "Envelope_Internal_Error_Max_ppm": record["internal_max"],
            "Envelope_Internal_Error_Mean_ppm": record["internal_mean"],
            "Envelope_Internal_Error_Median_ppm": record["internal_median"],
            "Source_RT_Window_IDs": "; ".join(record["rt_window_ids"]),
            "Num_Source_RT_Windows": len(record["rt_window_ids"]),
            "Merged_Across_RT_Windows": bool(record.get("merged") or len(record["rt_window_ids"]) > 1),
            "Extended_Lower_Charges_Evaluated": "; ".join(map(str, record.get("extension", {}).get("lower", []))),
            "Extended_Upper_Charges_Evaluated": "; ".join(map(str, record.get("extension", {}).get("upper", []))),
            "Extended_Charges_Detected": "; ".join(map(str, record.get("extension", {}).get("detected", []))),
            "Extended_Weak_Charges_Detected": "; ".join(map(str, record.get("extension", {}).get("weak", []))),
            "Extended_Charges_Not_Detected": "; ".join(map(str, record.get("extension", {}).get("not_detected", []))),
            "Charge_Extension_Improved_Envelope": bool(record.get("charge_extension_improved", False)),
            "Split_Envelope_Group_ID": record.get("split_envelope_group_id", ""),
            "Split_Envelope_Member_Count": record.get("split_envelope_member_count", 1),
            "Split_Envelope_Merged": record.get("split_envelope_merged", False),
            "Charge_Gaps_Before_Merge": record.get("charge_gaps_before_merge", record["gap_count"]),
            "Charge_Gaps_After_Merge": record.get("charge_gaps_after_merge", record["gap_count"]),
            "Max_Peak_Usage_Count": max_peak_usage,
            "Mean_Peak_Usage_Count": mean_peak_usage,
            "Num_Highly_Shared_Peaks": highly_shared_count,
            "Highly_Shared_Peak_Fraction": highly_shared_fraction,
            "Peak_Sharing_Status": peak_sharing_status,
            "Notes": "",
        })
    candidates.sort(key=lambda item: (item.charge_state_count < int(reconstruction_config.get("min_charge_states", 3)), -item.charge_state_count, -item.total_intensity))
    missing_weak = sum(1 for row in missing_diagnostics if row.get("Detection_Status") == "below_intensity_threshold")
    missing_not = sum(1 for row in missing_diagnostics if row.get("Detection_Status") == "no_peak_in_tolerance")
    metadata = {
        "engine": "rt_localized",
        "rt_envelope_diagnostics": rt_diagnostics,
        "missing_charge_diagnostics": missing_diagnostics,
        "engine_comparison": [],
        "stats": {
            "Reconstruction_Engine": "rt_localized",
            "Num_RT_Windows": len(windows),
            "Num_Local_Peaks": local_peak_count,
            "Num_Anchor_Peaks_Evaluated": anchor_count,
            "Num_Raw_Envelope_Candidates": raw_candidate_count,
            "Num_Candidates_After_Charge_Filter": charge_filtered_count,
            "Num_Candidates_After_RT_Window_Merge": len(candidates),
            "Num_Candidates_With_Consecutive_Charges": sum(1 for c in candidates if getattr(c, "longest_consecutive_charge_run", 0) >= rt_config["min_consecutive_charge_states"]),
            "Num_Candidates_With_Charge_Gaps": sum(1 for c in candidates if getattr(c, "charge_gap_count", 0) > 0),
            "Num_Missing_Charges_Evaluated": len(missing_diagnostics),
            "Num_Missing_Charges_With_Weak_Peaks": missing_weak,
            "Num_Missing_Charges_Not_Detected": missing_not,
            "Median_RT_Range_Min": median([c.rt_range_min for c in candidates if c.rt_range_min is not None]) if any(c.rt_range_min is not None for c in candidates) else None,
            "Median_Internal_Error_ppm": median([getattr(c, "envelope_internal_error_max_ppm", 0.0) or 0.0 for c in candidates]) if candidates else None,
            "Median_Charge_Count": median([len(c.charge_states) for c in candidates]) if candidates else None,
            "Processing_Time_Seconds": perf_counter() - started,
            "Split_Envelope_Count": sum(1 for c in candidates if getattr(c, "split_envelope_merged", False)),
            "Missing_Charge_Status_Count": _distribution_summary([row.get("Detection_Status") for row in missing_diagnostics]),
        },
    }
    return candidates, charge_state_peaks, metadata


def _legacy_reconstruct_intact_masses(
    tier_result: PeakTierResult,
    reconstruction_config: dict[str, Any],
    instrument_config: dict[str, Any],
    theoretical_mass: float | None = None,
) -> tuple[list[IntactMassCandidate], list[dict[str, Any]], dict[str, Any]]:
    started = perf_counter()
    min_charge = int(reconstruction_config.get("min_charge", 5))
    max_charge = int(reconstruction_config.get("max_charge", 40))
    min_charge_states = int(reconstruction_config.get("min_charge_states", 3))
    tolerance = float(reconstruction_config.get("mass_cluster_tolerance_da", 1.0))
    polarity = instrument_config.get("polarity", "negative")
    observations = []
    for peak in tier_result.usable_peaks:
        for charge in range(min_charge, max_charge + 1):
            neutral_mass = neutral_mass_from_mz(peak.mz, charge, polarity)
            observations.append({"peak": peak, "charge": charge, "neutral_mass": neutral_mass})
    observations.sort(key=lambda row: row["neutral_mass"])
    clusters: list[list[dict[str, Any]]] = []
    for observation in observations:
        if not clusters:
            clusters.append([observation])
            continue
        current_mean = sum(row["neutral_mass"] for row in clusters[-1]) / len(clusters[-1])
        if abs(observation["neutral_mass"] - current_mean) <= tolerance:
            clusters[-1].append(observation)
        else:
            clusters.append([observation])
    candidates: list[IntactMassCandidate] = []
    charge_state_peaks: list[dict[str, Any]] = []
    for index, cluster in enumerate(clusters, start=1):
        cluster_id = f"C{index:04d}"
        observed_mass = sum(row["neutral_mass"] for row in cluster) / len(cluster)
        charges = sorted({int(row["charge"]) for row in cluster})
        total_intensity = sum(float(row["peak"].intensity) for row in cluster)
        mass_error_da = observed_mass - theoretical_mass if theoretical_mass is not None else None
        mass_error_ppm = (mass_error_da / theoretical_mass * 1_000_000) if theoretical_mass else None
        candidate = IntactMassCandidate(
            observed_mass=observed_mass,
            charge_state_count=len(charges),
            charge_states=charges,
            supporting_peak_count=len(cluster),
            total_intensity=total_intensity,
            theoretical_mass=theoretical_mass,
            mass_error_da=mass_error_da,
            mass_error_ppm=mass_error_ppm,
            confidence=_confidence(len(charges), min_charge_states),
            cluster_id=cluster_id,
        )
        _set_candidate_extra(candidate, reconstruction_engine="legacy_cluster")
        candidates.append(candidate)
        for row in cluster:
            peak = row["peak"]
            charge_state_peaks.append({
                "Cluster_ID": cluster_id,
                "mz": peak.mz,
                "Intensity": peak.intensity,
                "RT": peak.rt,
                "Scan_ID": peak.scan_id,
                "Charge": row["charge"],
                "Neutral_Mass": row["neutral_mass"],
                "Peak_Tier": peak.tier,
            })
    candidates.sort(key=lambda item: (item.charge_state_count < min_charge_states, -item.charge_state_count, -item.total_intensity))
    metadata = {
        "engine": "legacy_cluster",
        "rt_envelope_diagnostics": [],
        "missing_charge_diagnostics": [],
        "engine_comparison": [],
        "stats": {
            "Reconstruction_Engine": "legacy_cluster",
            "Num_RT_Windows": 0,
            "Num_Local_Peaks": 0,
            "Num_Anchor_Peaks_Evaluated": 0,
            "Num_Raw_Envelope_Candidates": len(observations),
            "Num_Candidates_After_Charge_Filter": len(candidates),
            "Num_Candidates_After_RT_Window_Merge": len(candidates),
            "Num_Candidates_With_Consecutive_Charges": sum(1 for c in candidates if _longest_consecutive_run(c.charge_states) >= 2),
            "Num_Candidates_With_Charge_Gaps": sum(1 for c in candidates if _charge_gap_count(c.charge_states) > 0),
            "Num_Missing_Charges_Evaluated": 0,
            "Num_Missing_Charges_With_Weak_Peaks": 0,
            "Num_Missing_Charges_Not_Detected": 0,
            "Median_RT_Range_Min": None,
            "Median_Internal_Error_ppm": None,
            "Median_Charge_Count": median([len(c.charge_states) for c in candidates]) if candidates else None,
            "Processing_Time_Seconds": perf_counter() - started,
        },
    }
    return candidates, charge_state_peaks, metadata



def _finalize_engine_stats(metadata: dict[str, Any], candidates: list[IntactMassCandidate]) -> None:
    stats = metadata.setdefault("stats", {})
    rt_ranges = [candidate.rt_range_min for candidate in candidates if candidate.rt_range_min is not None]
    internal_errors = [candidate.envelope_internal_error_ppm for candidate in candidates if candidate.envelope_internal_error_ppm is not None]
    charge_counts = [len(candidate.charge_states) for candidate in candidates]
    if rt_ranges:
        stats["Median_RT_Range_Min"] = median(rt_ranges)
    if internal_errors:
        stats["Median_Internal_Error_ppm"] = median(internal_errors)
    if charge_counts:
        stats["Median_Charge_Count"] = median(charge_counts)

def build_intact_engine_comparison_rows(
    legacy_candidates: list[IntactMassCandidate],
    legacy_peaks: list[dict[str, Any]],
    rt_candidates: list[IntactMassCandidate],
    rt_peaks: list[dict[str, Any]],
    reconstruction_config: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    qc_config = _qc_config(reconstruction_config or {})
    match_cfg = qc_config["engine_comparison"]
    legacy_peak_sets: dict[str, set[str]] = {}
    rt_peak_sets: dict[str, set[str]] = {}
    for row in legacy_peaks:
        legacy_peak_sets.setdefault(str(row.get("Cluster_ID")), set()).add(_supporting_peak_id(row))
    for row in rt_peaks:
        rt_peak_sets.setdefault(str(row.get("Cluster_ID")), set()).add(_supporting_peak_id(row))

    rt_by_mass = sorted(rt_candidates, key=lambda candidate: float(candidate.observed_mass or 0.0))
    rt_masses = [float(candidate.observed_mass or 0.0) for candidate in rt_by_mass]

    def _rt_candidates_for_legacy(legacy: IntactMassCandidate) -> tuple[list[IntactMassCandidate], bool]:
        if not rt_by_mass:
            return [], False
        if not match_cfg.get("require_mass_match", True) or not legacy.observed_mass:
            return rt_by_mass, True
        tolerance_da = abs(float(legacy.observed_mass)) * match_cfg["mass_tolerance_ppm"] / 1_000_000
        lower = float(legacy.observed_mass) - tolerance_da
        upper = float(legacy.observed_mass) + tolerance_da
        left = bisect_left(rt_masses, lower)
        right = bisect_left(rt_masses, upper)
        while right < len(rt_masses) and rt_masses[right] <= upper:
            right += 1
        return rt_by_mass[left:right], right > left

    rows = []
    used_rt: set[str] = set()
    for legacy in legacy_candidates:
        candidates_to_check, has_mass_window_candidates = _rt_candidates_for_legacy(legacy)
        if not candidates_to_check:
            rows.append({
                "Legacy_Cluster_ID": legacy.cluster_id,
                "Legacy_Mass": legacy.observed_mass,
                "Legacy_Charge_Count": len(legacy.charge_states),
                "Engine_Match_Status": "mass_mismatch" if rt_candidates and not has_mass_window_candidates else "legacy_only",
                "Match_Failure_Reason": "no_rt_candidate_within_mass_tolerance" if rt_candidates and not has_mass_window_candidates else "no_rt_localized_candidate",
                "Mass_Matched": False if rt_candidates else None,
                "RT_Matched": None,
                "Charge_Matched": None,
                "Peak_Matched": None,
                "Notes": "legacy_not_matched" if rt_candidates else "legacy_only",
            })
            continue
        best_row = None
        best_score = None
        for rt in candidates_to_check:
            mass_error_ppm = abs(_ppm_error(rt.observed_mass, legacy.observed_mass)) if legacy.observed_mass else float("inf")
            mass_matched = mass_error_ppm <= match_cfg["mass_tolerance_ppm"]
            legacy_rt = getattr(legacy, "rt_mean", None)
            rt_rt = getattr(rt, "rt_mean", None)
            rt_delta = abs(float(rt_rt) - float(legacy_rt)) if rt_rt is not None and legacy_rt is not None else None
            rt_matched = rt_delta is None or rt_delta <= match_cfg["rt_tolerance_min"]
            shared_charges = len(set(legacy.charge_states) & set(rt.charge_states))
            charge_overlap = shared_charges / max(min(len(legacy.charge_states), len(rt.charge_states)), 1)
            charge_matched = charge_overlap >= match_cfg["min_shared_charge_fraction"]
            left = legacy_peak_sets.get(legacy.cluster_id, set())
            right = rt_peak_sets.get(rt.cluster_id, set())
            peak_overlap = len(left & right) / max(min(len(left), len(right)), 1) if (left or right) else 0.0
            peak_matched = peak_overlap > 0.0
            if not mass_matched and match_cfg.get("require_mass_match", True):
                status = "mass_mismatch"
            elif not rt_matched:
                status = "rt_mismatch"
            elif not charge_matched:
                status = "insufficient_overlap"
            else:
                status = "matched"
            score = (status != "matched", mass_error_ppm, rt_delta if rt_delta is not None else 0.0, -charge_overlap, -peak_overlap)
            row = {
                "Legacy_Cluster_ID": legacy.cluster_id,
                "RT_Localized_Cluster_ID": rt.cluster_id,
                "Legacy_Mass": legacy.observed_mass,
                "RT_Localized_Mass": rt.observed_mass,
                "Mass_Delta_Da": (rt.observed_mass - legacy.observed_mass) if status == "matched" else None,
                "Legacy_Charge_Count": len(legacy.charge_states),
                "RT_Localized_Charge_Count": len(rt.charge_states),
                "Legacy_RT_Range": getattr(legacy, "rt_range_min", None),
                "RT_Localized_RT_Range": getattr(rt, "rt_range_min", None),
                "Legacy_Internal_Error_ppm": getattr(legacy, "envelope_internal_error_ppm", None),
                "RT_Localized_Internal_Error_ppm": getattr(rt, "envelope_internal_error_ppm", None),
                "Peak_Overlap_Fraction": peak_overlap,
                "Engine_Match_Status": status,
                "Match_Failure_Reason": "" if status == "matched" else status,
                "Mass_Matched": mass_matched,
                "RT_Matched": rt_matched,
                "Charge_Matched": charge_matched,
                "Peak_Matched": peak_matched,
                "Notes": "criteria_match" if status == "matched" else "not_matched",
            }
            if best_row is None or score < best_score:
                best_row = row
                best_score = score
                if status == "matched":
                    break
        if best_row is None:
            continue
        if best_row["Engine_Match_Status"] == "matched":
            used_rt.add(str(best_row["RT_Localized_Cluster_ID"]))
            rows.append(best_row)
        else:
            rows.append({
                **best_row,
                "RT_Localized_Cluster_ID": "",
                "RT_Localized_Mass": None,
                "RT_Localized_Charge_Count": None,
                "RT_Localized_RT_Range": None,
                "RT_Localized_Internal_Error_ppm": None,
                "Peak_Overlap_Fraction": None,
                "Notes": "legacy_not_matched",
            })
    for rt in rt_candidates:
        if rt.cluster_id not in used_rt:
            rows.append({
                "Legacy_Cluster_ID": "",
                "RT_Localized_Cluster_ID": rt.cluster_id,
                "RT_Localized_Mass": rt.observed_mass,
                "RT_Localized_Charge_Count": len(rt.charge_states),
                "Engine_Match_Status": "rt_localized_only",
                "Match_Failure_Reason": "no_matched_legacy_candidate",
                "Notes": "rt_localized_only",
            })
    return rows

def _engine_comparison_stats(rows: list[dict[str, Any]]) -> dict[str, Any]:
    status_counts: dict[str, int] = {}
    for row in rows:
        status = str(row.get("Engine_Match_Status") or "unknown")
        status_counts[status] = status_counts.get(status, 0) + 1
    matched_deltas = [abs(float(row["Mass_Delta_Da"])) for row in rows if row.get("Engine_Match_Status") == "matched" and row.get("Mass_Delta_Da") is not None]
    return {
        "Engine_Match_Status_Count": "; ".join(f"{key}:{status_counts[key]}" for key in sorted(status_counts)),
        "Matched_Engine_Candidate_Count": status_counts.get("matched", 0),
        "Legacy_Only_Count": status_counts.get("legacy_only", 0),
        "RT_Localized_Only_Count": status_counts.get("rt_localized_only", 0),
        "Mass_Mismatch_Count": status_counts.get("mass_mismatch", 0),
        "RT_Mismatch_Count": status_counts.get("rt_mismatch", 0),
        "Charge_Overlap_Failure_Count": status_counts.get("insufficient_overlap", 0),
        "Median_Mass_Delta_For_Matched_Only": median(matched_deltas) if matched_deltas else None,
    }




def _sync_competition_fields_to_rt_diagnostics(metadata: dict[str, Any], candidates: list[IntactMassCandidate]) -> None:
    diagnostics = metadata.get("rt_envelope_diagnostics") or []
    by_cluster = {candidate.cluster_id: candidate for candidate in candidates}
    for row in diagnostics:
        candidate = by_cluster.get(row.get("Cluster_ID"))
        if candidate is None:
            continue
        row["Competing_Envelope_Group_ID"] = candidate.competing_envelope_group_id
        row["Competing_Envelope_Group_Size"] = candidate.competing_envelope_group_size
        row["Envelope_Evidence_Score"] = candidate.envelope_evidence_score
        row["Evidence_Score_Rank_In_Competition"] = candidate.evidence_score_rank_in_competition
        row["Maximum_Shared_Peak_Fraction"] = candidate.maximum_shared_peak_fraction
        row["Direct_Competitor_Count"] = candidate.direct_competitor_count
        row["Dry_Run_Assignment_Status"] = candidate.dry_run_assignment_status
        row["Dry_Run_Selected"] = candidate.dry_run_selected
        row["Independent_Supporting_Peak_Count"] = candidate.independent_supporting_peak_count
        row["Independent_Supporting_Peak_Fraction"] = candidate.independent_supporting_peak_fraction
        row["Independent_Charge_State_Count"] = candidate.independent_charge_state_count
        row["Assignment_Confidence"] = candidate.assignment_confidence


def reconstruct_intact_masses(
    tier_result: PeakTierResult,
    reconstruction_config: dict[str, Any],
    instrument_config: dict[str, Any],
    theoretical_mass: float | None = None,
    warnings: list[dict[str, Any]] | None = None,
) -> tuple[list[IntactMassCandidate], list[dict[str, Any]], dict[str, Any]]:
    qc_config = _qc_config(reconstruction_config or {})
    engine = qc_config.get("engine", "legacy_cluster")
    if engine == "rt_localized":
        candidates, charge_state_peaks, metadata = _build_rt_localized_candidates(
            tier_result,
            reconstruction_config,
            instrument_config,
            theoretical_mass,
        )
        reconstruction_config["_intact_engine_stats"] = metadata.get("stats", {})
        build_intact_reconstruction_qc(candidates, charge_state_peaks, reconstruction_config, reconstruction_enabled=True)
        _sync_competition_fields_to_rt_diagnostics(metadata, candidates)
        _finalize_engine_stats(metadata, candidates)
        reconstruction_config["_intact_engine_stats"] = metadata.get("stats", {})
        if qc_config.get("compare_with_legacy", False):
            legacy_candidates, legacy_peaks, legacy_metadata = _legacy_reconstruct_intact_masses(
                tier_result,
                reconstruction_config,
                instrument_config,
                theoretical_mass,
            )
            build_intact_reconstruction_qc(legacy_candidates, legacy_peaks, reconstruction_config, reconstruction_enabled=True)
            metadata["engine_comparison"] = build_intact_engine_comparison_rows(legacy_candidates, legacy_peaks, candidates, charge_state_peaks, reconstruction_config)
            metadata["stats"].update(_engine_comparison_stats(metadata["engine_comparison"]))
            metadata["legacy_candidate_count"] = len(legacy_candidates)
            metadata["legacy_processing_time_seconds"] = legacy_metadata.get("stats", {}).get("Processing_Time_Seconds")
        return candidates, charge_state_peaks, metadata
    candidates, charge_state_peaks, metadata = _legacy_reconstruct_intact_masses(
        tier_result,
        reconstruction_config,
        instrument_config,
        theoretical_mass,
    )
    reconstruction_config["_intact_engine_stats"] = metadata.get("stats", {})
    build_intact_reconstruction_qc(candidates, charge_state_peaks, reconstruction_config, reconstruction_enabled=True)
    _finalize_engine_stats(metadata, candidates)
    reconstruction_config["_intact_engine_stats"] = metadata.get("stats", {})
    return candidates, charge_state_peaks, metadata
