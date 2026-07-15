"""Audit-level policy, sheet registry, and non-propagating status metadata."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

FORMAL_CORE = "FORMAL_CORE"
FORMAL_OPTIONAL = "FORMAL_OPTIONAL"
AUDIT_SUMMARY = "AUDIT_SUMMARY"
AUDIT_GROUP = "AUDIT_GROUP"
AUDIT_DETAIL = "AUDIT_DETAIL"
INTERNAL_ONLY = "INTERNAL_ONLY"
AUDIT_LEVELS = ("standard", "audit", "full")
AUDIT_LEVEL_DEFAULT = "full"


@dataclass(frozen=True)
class AuditPolicy:
    level: str
    run_shadow_audits: bool
    include_summary: bool
    include_group_tables: bool
    include_detail: bool
    include_top_shadow_columns: bool
    include_diagnostics_shadow_columns: bool

    @classmethod
    def from_level(cls, level: str | None) -> "AuditPolicy":
        value = str(level or AUDIT_LEVEL_DEFAULT).lower()
        if value not in AUDIT_LEVELS:
            raise ValueError(f"Unknown audit level: {level}")
        if value == "standard":
            return cls(value, False, False, False, False, False, False)
        if value == "audit":
            return cls(value, True, True, False, False, True, True)
        return cls(value, True, True, True, True, True, True)

    def includes_category(self, category: str) -> bool:
        if category in {FORMAL_CORE, FORMAL_OPTIONAL}:
            return True
        if category == AUDIT_SUMMARY:
            return self.include_summary
        if category == AUDIT_GROUP:
            return self.include_group_tables
        if category == AUDIT_DETAIL:
            return self.include_detail
        return False

    def run_summary_items(self) -> list[dict[str, Any]]:
        return [
            {"Item":"Audit_Level","Value":self.level},
            {"Item":"Shadow_Audits_Run","Value":self.run_shadow_audits},
            {"Item":"Shadow_Summary_Included","Value":self.include_summary},
            {"Item":"Shadow_Detail_Included","Value":self.include_detail},
            {"Item":"Top_Shadow_Columns_Included","Value":self.include_top_shadow_columns},
            {"Item":"Diagnostics_Shadow_Columns_Included","Value":self.include_diagnostics_shadow_columns},
            {"Item":"Formal_Result_Changed_By_Audit_Level","Value":False},
            {"Item":"Audit_Level_Default","Value":AUDIT_LEVEL_DEFAULT},
        ]


FORMAL_CORE_SHEETS = {
    "Run_summary", "Workflow_Summary", "Input_parameters", "mzML_diagnostics",
    "Theoretical_fragments", "Fragment_MS1_matches", "Fragment_MS1_filtered",
    "Fragment_MS1_summary", "Known_Modification_Candidates", "Known_Modification_Summary",
    "Modification_Evidence_Summary", "Modification_Evidence_Ranking",
    "Modification_Ambiguity_Groups", "Modification_Position_Priors",
    "MS2_Biological_Plausibility", "MS2_Modification_Identity",
    "MS2_Identity_Peak_Assignments", "Biological_Prior_Diagnostics",
    "Biological_Context_Priorities", "Context_Supported_Candidates",
    "Review_Dashboard", "Top_Modification_Candidates", "Candidate_Decision_Summary",
    "Evidence_Checklist", "MS2_Unmatched_Ion_Diagnostics", "Warnings",
}
FORMAL_OPTIONAL_SHEETS = {
    "Intact_mass_reconstruction", "Charge_state_peaks", "Intact_Reconstruction_QC",
    "Intact_Reconstruction_Diag", "Intact_Envelope_Groups", "Intact_Comparison_Candidates",
    "Target_Review_Candidates", "Reconstructed_Mass_Spectrum", "RT_Envelope_Diagnostics",
    "RT_Engine_QC_Summary", "Missing_Charge_Diagnostics", "Intact_Engine_Comparison",
    "P1_Summary", "P1_Theoretical_Structures", "P1_Peak_Annotations", "P1_Unmatched_Peaks",
    "MS2_Summary", "MS2_Spectra", "MS2_Parent_Candidates", "MS2_Theoretical_Ions",
    "MS2_Ion_Matches", "MS2_Unmatched_Peaks", "MS2_Fragment_Evidence",
    "MS2_Peak_Annotations", "MS2_Modified_Precursor_Candidates",
    "MS2_Modified_Theoretical_Ions", "MS2_Modified_Ion_Matches",
    "MS2_Modification_Localization_Evidence",
    # Explicit Excel 31-character aliases used by the report writer.
    "MS2_Modified_Precursor_Candidat", "MS2_Modification_Localization_E",
}
AUDIT_SUMMARY_SHEETS = {
    "Audit_Status", "MS2_Unmatched_Ion_Summary", "MS2_Ambiguity_Summary",
    "MS2_Zero_Intensity_Summary", "MS2_Effective_Ambig_Summary",
    "MS1_Truncation_Summary", "MS1_Selection_Summary", "MS1_Top50_Dedup_Summary",
    "MS1_CrossFrag_Summary", "Competition_Dry_Run_Summary",
    "Composite_Mod_Summary", "Cleavage_Block_Audit",
    "Composite_MS1_Summary", "Composite_Support_Summary", "Legacy_Composite_Compare",
    "Composite_Shadow_Score", "Composite_Obs_Summary",
    "PT_Paired_Summary", "PT_Discovery_Candidates",
    "PT_Cross_Run_Runs", "PT_Cross_Run_Summary", "PT_Cross_Run_Neutral",
    "PT_Cross_Run_Pairs", "PT_Cross_Run_Decoy", "PT_Cross_Run_Target",
}
AUDIT_GROUP_SHEETS = {
    "Intact_Competition_Groups", "Intact_Competition_Scores", "Intact_Assignment_Dry_Run",
    "Assignment_Sensitivity", "Assignment_Stability", "Assignment_Candidate_Audit",
    "Assignment_Ambiguous_Candidates", "Preassignment_Comparison",
    "MS2_Ambiguous_Peak_Clusters", "MS2_Zero_Intensity_Spectra",
    "MS2_Effective_Ambiguity", "MS1_Truncation_Audit", "MS1_Selection_Strategy",
    "MS1_Top50_Shadow", "MS1_CrossFrag_Ambiguity",
}
AUDIT_DETAIL_SHEETS = {
    "MS2_Unmatched_Ion_Audit", "MS2_Ambiguous_Peak_Detail", "MS2_Zero_Intensity_Detail",
    "MS2_Effective_Ambig_Detail", "MS1_Truncation_Detail", "MS1_Selection_Detail",
    "MS1_Peak_Dedup_Detail", "MS1_CrossFrag_Detail",
    "Composite_Mod_Candidates", "Composite_Mod_Invalid", "Backbone_Mod_Candidates",
    "Composite_Fragment_Masses", "Composite_MS1_Matches", "Composite_MS2_Ions",
    "Composite_MS2_Matches", "Blocked_Cleavage_Matches", "Composite_Obs_Invalid",
    "PT_Paired_Evidence", "PT_State_Search", "PT_Cross_Run_Detail", "PT_Cross_Run_Decoy_Detail", "PT_Cross_Run_MS2_Detail",
}

SHEET_REGISTRY = {
    **{name:FORMAL_CORE for name in FORMAL_CORE_SHEETS},
    **{name:FORMAL_OPTIONAL for name in FORMAL_OPTIONAL_SHEETS},
    **{name:AUDIT_SUMMARY for name in AUDIT_SUMMARY_SHEETS},
    **{name:AUDIT_GROUP for name in AUDIT_GROUP_SHEETS},
    **{name:AUDIT_DETAIL for name in AUDIT_DETAIL_SHEETS},
}


def sheet_category(name: str) -> str | None:
    if str(name).startswith("_"):
        return INTERNAL_ONLY
    return SHEET_REGISTRY.get(str(name))


def unclassified_sheets(names: Iterable[str]) -> list[str]:
    return sorted({str(name) for name in names if sheet_category(str(name)) is None and str(name) != "Index"})


def included_sheet_names(names: Iterable[str], policy: AuditPolicy) -> tuple[list[str], list[str]]:
    included=[]; unknown=[]
    for name in names:
        category=sheet_category(name)
        if category is None:
            unknown.append(name)
            if policy.level == "full": included.append(name)
        elif policy.includes_category(category):
            included.append(name)
    return included, sorted(set(unknown))


AUDIT_STATUS_COLUMNS = [
    "Audit_Name", "Category", "Enabled_For_Level", "Executed", "Summary_Available",
    "Detail_Available", "Applied_To_Formal_Result", "Runtime_Seconds", "Peak_Memory_MB",
    "Status", "Reason_Not_Run",
]
DIAGNOSTIC_COLUMNS = [
    "Audit_Level", "Shadow_Audits_Run", "Shadow_Audit_Summary_Available",
    "Shadow_Audit_Detail_Available", "Shadow_Audit_Sheet_Count",
    "Shadow_Audit_Runtime_Seconds", "Shadow_Audit_Peak_Memory_MB",
    "Formal_Result_Changed_By_Audit_Level",
]


def append_audit_level_diagnostics(rows: Any, policy: AuditPolicy, status_rows: list[dict[str, Any]], sheet_count: int = 0):
    import pandas as pd
    executed=[row for row in status_rows if row.get("Executed") is True]
    values={
        "Audit_Level":policy.level,
        "Shadow_Audits_Run":policy.run_shadow_audits,
        "Shadow_Audit_Summary_Available":policy.include_summary and bool(executed),
        "Shadow_Audit_Detail_Available":policy.include_detail and bool(executed),
        "Shadow_Audit_Sheet_Count":sheet_count if policy.run_shadow_audits else "not_run",
        "Shadow_Audit_Runtime_Seconds":sum(float(row.get("Runtime_Seconds") or 0) for row in executed) if executed else "not_run",
        "Shadow_Audit_Peak_Memory_MB":max((float(row.get("Peak_Memory_MB") or 0) for row in executed),default=0) if executed else "not_run",
        "Formal_Result_Changed_By_Audit_Level":False,
    }
    is_frame=isinstance(rows,pd.DataFrame)
    source=rows.to_dict("records") if is_frame else list(rows or [{}])
    out=[dict(row,**values) for row in source]
    return pd.DataFrame(out,columns=list(rows.columns)+DIAGNOSTIC_COLUMNS) if is_frame else out


def audit_status_row(name: str, category: str, policy: AuditPolicy, executed: bool, summary: bool, detail: bool, runtime: float = 0.0, peak_mb: float = 0.0, reason: str = "") -> dict[str, Any]:
    return {
        "Audit_Name":name,"Category":category,"Enabled_For_Level":policy.run_shadow_audits,
        "Executed":executed,"Summary_Available":summary and executed,"Detail_Available":detail and executed,
        "Applied_To_Formal_Result":False,"Runtime_Seconds":runtime if executed else "",
        "Peak_Memory_MB":peak_mb if executed else "","Status":"completed" if executed else "not_run",
        "Reason_Not_Run":"" if executed else (reason or f"audit_level={policy.level}"),
    }
