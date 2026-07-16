"""Phase-2 observation connection for explicit complete-structure hypotheses."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from rna_masshunter.backbone_state import load_backbone_transformations
from rna_masshunter.blocked_cleavage_matcher import match_blocked_cleavage_fragments
from rna_masshunter.candidate_support import aggregate_candidate_support
from rna_masshunter.composite_fragment_mass import fragment_mass_row
from rna_masshunter.composite_ms1_matcher import match_composite_fragments_to_peaks
from rna_masshunter.composite_ms2_matcher import match_composite_ms2
from rna_masshunter.composite_ms2_propagation import generate_composite_theoretical_ions
from rna_masshunter.legacy_composite_comparison import compare_legacy_composite
from rna_masshunter.modification_constraints import load_transformations
from rna_masshunter.sample_structure_schema import load_sample_structure_hypotheses
from rna_masshunter.shadow_score_simulation import simulate_shadow_scores
from rna_masshunter.structure_fragment import build_complete_structure_state, extract_fragment_from_structure

COMMON_STATUS_COLUMNS = ["Audit_Level", "Applied_To_Formal_Result", "Formal_Change_Ready"]
FRAGMENT_COLUMNS = ["Candidate_ID","Complete_Structure_ID","Fragment_ID","Fragment_Type","Start_Position","End_Position","Included_Positions","Included_Modified_Positions","Included_Backbone_Bonds","Included_Backbone_Modifications","Terminal_Form","Fragment_Elemental_Composition","Neutral_Exact_Mass"] + COMMON_STATUS_COLUMNS
MS1_COLUMNS = ["Candidate_ID","Complete_Structure_ID","Fragment_ID","Fragment_Type","Start_Position","End_Position","Included_Modified_Positions","Included_Backbone_Bonds","Neutral_Exact_Mass","Charge","Theoretical_mz","Observed_mz","Mass_Error_Da","Mass_Error_ppm","Observed_Intensity","Observed_Scan","Observed_RT","Match_Status","Support_Class","Not_Observable_Reason","Legacy_Competition_Class","Is_Isomeric"] + COMMON_STATUS_COLUMNS
MS1_SUMMARY_COLUMNS = ["Candidate_ID","Matched_Count","No_Observation_Count","Not_Observable_Count","Unique_Support_Count","Shared_Legacy_Count","Shared_Composite_Count","Nondiscriminating_Count","Isomeric_Unresolved_Count"] + COMMON_STATUS_COLUMNS
MS2_ION_COLUMNS = ["Candidate_ID","Complete_Structure_ID","Parent_Fragment_ID","Parent_Sequence","Parent_Start","Parent_End","Parent_Neutral_Mass","Ion_ID","Ion_Series","Ion_Number","Cleavage_Position","Ion_Sequence","Included_Positions","Included_Modified_Positions","Included_Backbone_Bonds","Theoretical_Neutral_Mass","Charge","Theoretical_mz","Position_Informative","Backbone_Informative"] + COMMON_STATUS_COLUMNS
MS2_MATCH_COLUMNS = ["Candidate_ID","Complete_Structure_ID","Spectrum_ID","Precursor_mz","Precursor_Charge","Ion_Series","Ion_Number","Cleavage_Position","Included_Positions","Included_Modified_Positions","Included_Backbone_Bonds","Theoretical_Neutral_Mass","Theoretical_mz","Observed_mz","Mass_Error_Da","Mass_Error_ppm","Observed_Intensity","Position_Informative","Backbone_Informative","Candidate_Discriminating","Isomer_Discriminating","Legacy_Competition_Class"] + COMMON_STATUS_COLUMNS
MS2_ASSIGNMENT_COMPETITION_COLUMNS = [
    "Composite_Match_ID", "Physical_Observed_Peak_Key", "Observed_Peak_Index",
    "Raw_Peak_Index", "Raw_Peak_Index_Missing_Reason", "Scan_Index", "Spectrum_ID",
    "RT", "Observed_mz", "Observed_Intensity", "Observed_Intensity_State", "Ion_ID",
    "Candidate_ID", "Complete_Structure_ID", "Parent_Fragment_ID", "Ion_Series",
    "Ion_Number", "Cleavage_Position", "Charge", "Theoretical_mz", "Mass_Error_Da",
    "Mass_Error_ppm", "Assignment_Rank", "Best_Assignment",
    "Within_Tolerance_Assignment_Count", "Competing_Candidate_Count",
    "Competing_Candidate_IDs", "Competing_Complete_Structure_Count",
    "Competing_Complete_Structure_IDs", "Competing_Theoretical_Ion_Count",
    "Competing_Ion_IDs", "Best_Error_ppm", "Second_Best_Error_ppm",
    "Best_vs_Second_Error_Margin_ppm", "Candidate_Specific",
    "Complete_Structure_Specific", "Theoretical_Ion_Specific", "Position_Specific",
    "Backbone_Bond_Specific", "Included_Positions", "Included_Modified_Positions",
    "Included_Backbone_Bonds", "Position_Informative", "Backbone_Informative",
    "Candidate_Discriminating", "Isomer_Discriminating", "Legacy_Competition_Class",
    "Audit_Level", "Applied_To_Formal_Result", "Formal_Change_Ready",
    "Formal_Result_Changed",
]
SUPPORT_COLUMNS = ["Candidate_ID","Complete_Structure_ID","Parent_Base","Modified_Positions","Backbone_Modified_Bonds","Theoretical_Fragment_Count","Observable_Fragment_Count","MS1_Matched_Fragment_Count","MS1_Unique_Support_Count","MS1_Shared_Support_Count","MS1_Nondiscriminating_Count","MS1_Isomeric_Unresolved_Count","MS2_Matched_Ion_Count","MS2_Position_Informative_Count","MS2_Backbone_Informative_Count","Blocked_Cleavage_Match_Count","Conflicting_Observation_Count","Support_Coverage","Support_Status"] + COMMON_STATUS_COLUMNS
BLOCKED_COLUMNS = ["Candidate_ID","Enzyme","Fragment_ID","Start_Position","End_Position","Blocked_Bond_ID","Blocked_Cleavage_Position","Cleavage_Status","Blocked_Cleavage_Reason","Contains_Phosphorothioate","Backbone_Modification_Count","Backbone_Modification_Positions","Backbone_Mass_Delta","Theoretical_mz","Observed_mz","Mass_Error_ppm","Observed_Intensity","Alternative_Stochastic_Fragment_Exists","Mechanism_Discriminating"] + COMMON_STATUS_COLUMNS
COMPARE_COLUMNS = ["Candidate_ID","Legacy_Candidate_IDs","Legacy_Parent_IDs","Exact_Phase1_Candidate_IDs","Comparison_Class","Neutral_Mass_Equivalent","Position_Compatible","Chemical_Exclusivity_Checked","MS1_Support_Count","MS1_Nondiscriminating_Count","MS2_Localization_Support_Count","Legacy_Formal_Ranks"] + COMMON_STATUS_COLUMNS
SCORE_COLUMNS = ["Candidate_ID","Baseline_Status","Legacy_Formal_Score","Composite_Shadow_Support","Composite_Shadow_Penalty","Simulated_Shadow_Score","Formal_Rank","Simulated_Shadow_Rank","Score_Delta","Rank_Delta","Would_Change_Top_Candidate","Would_Change_Confidence","Would_Change_Formal_Result"] + COMMON_STATUS_COLUMNS
OBS_SUMMARY_COLUMNS = ["Schema_Version","Enabled","Valid_Hypothesis_Count","Invalid_Hypothesis_Count","Fragment_Mass_Count","MS1_Match_Count","MS2_Ion_Count","MS2_Match_Count","Blocked_Cleavage_Match_Count","Formal_Result_Changed","Remaining_Risk"] + COMMON_STATUS_COLUMNS
INVALID_COLUMNS = ["Candidate_ID","Valid","Invalid_Reason","Invalid_Detail","Audit_Level","Applied_To_Formal_Result","Formal_Change_Ready"]

@dataclass(frozen=True)
class CompositeObservationResult:
    sheets: dict[str, list[dict[str, Any]]]
    structures: tuple[Any, ...]
    invalid_rows: tuple[dict[str, Any], ...]

def _ms1_summary(structures, rows, level):
    output = []
    for structure in structures:
        group = [r for r in rows if r["Candidate_ID"] == structure.candidate_id]
        output.append({"Candidate_ID": structure.candidate_id,
            "Matched_Count": sum(r["Match_Status"] == "matched" for r in group),
            "No_Observation_Count": sum(r["Match_Status"] == "no_observation" for r in group),
            "Not_Observable_Count": sum(r["Match_Status"] == "not_observable" for r in group),
            "Unique_Support_Count": sum(r["Support_Class"] == "unique_composite_support" for r in group),
            "Shared_Legacy_Count": sum(r["Support_Class"] == "shared_with_legacy" for r in group),
            "Shared_Composite_Count": sum(r["Support_Class"] == "shared_with_other_composite" for r in group),
            "Nondiscriminating_Count": sum(r["Support_Class"] == "observation_nondiscriminating" for r in group),
            "Isomeric_Unresolved_Count": sum(r["Support_Class"] == "isomeric_unresolved" for r in group),
            "Audit_Level": level, "Applied_To_Formal_Result": False, "Formal_Change_Ready": False})
    return output

def build_composite_observation_audit(project_root: str | Path, sequence: str,
    theoretical_fragments: list[Any], peaks: list[Any], spectra: list[Any], config: Any,
    base_masses: dict, phase1_sheets: dict[str, Any], formal_ms1_matches: list[Any],
    formal_ranking: list[dict[str, Any]], *, audit_level: str = "full",
    fixture_path: str | Path | None = None) -> CompositeObservationResult:
    root = Path(project_root)
    transforms = load_transformations(root / "data/modification_transforms_v2.yaml")
    backbone_transform = load_backbone_transformations(root / "data/backbone_modifications.yaml")[0]
    fixture = Path(fixture_path) if fixture_path else root / "data/sample_structure_hypotheses.yaml"
    loaded = load_sample_structure_hypotheses(fixture, sequence=sequence, transformations=transforms,
        backbone_bond_ids={f"{i}_{i+1}" for i in range(1, len(sequence))}, target_identity={
            "name": (getattr(config, "sequence", {}) or {}).get("name", ""),
            "organism": (getattr(config, "organism", {}) or {}).get("species", ""),
            "rule_set": (getattr(config, "organism", {}) or {}).get("rule_set", ""),
        })
    invalid = [dict(row, Audit_Level=audit_level) for row in loaded.invalid_rows]; structures = []
    for hypothesis in loaded.hypotheses:
        structure, error = build_complete_structure_state(hypothesis, sequence, transforms,
            root / "data/nucleoside_slots.yaml", backbone_transform)
        if error: invalid.append(dict(error, Audit_Level=audit_level))
        elif structure is not None: structures.append(structure)
    fragments = []
    for structure in structures:
        for parent in theoretical_fragments:
            fragments.append(extract_fragment_from_structure(structure, sequence, int(parent.start), int(parent.end),
                fragment_id=f"{structure.candidate_id}|{parent.fragment_id}", fragment_type=str(parent.enzyme),
                terminal_form=str(parent.terminal_form)))
    fragment_rows = [fragment_mass_row(item, audit_level) for item in fragments]
    phase1_candidates = phase1_sheets.get("Composite_Mod_Candidates", [])
    isomer_groups = {}
    for structure in structures:
        canonical_ids = {state.canonical_structure_id for state in structure.position_states.values()}
        matching_isomers = [row for row in phase1_candidates
            if row.get("Complete_Structure_ID") in canonical_ids and row.get("Is_Isomeric")]
        if matching_isomers:
            isomer_groups[structure.candidate_id] = str(matching_isomers[0].get("Isomer_Group_ID") or "isomeric")
    ms1_rows = match_composite_fragments_to_peaks(fragments, peaks, config,
        legacy_matches=formal_ms1_matches, isomer_groups=isomer_groups, audit_level=audit_level)
    ions = generate_composite_theoretical_ions(structures, theoretical_fragments, sequence, config, audit_level=audit_level)
    ms2_rows, ms2_competition = match_composite_ms2(
        ions, spectra, config, audit_level=audit_level, return_competition=True,
    )
    blocked = match_blocked_cleavage_fragments(structures, sequence, peaks, config,
        backbone_transform, base_masses, theoretical_fragments, audit_level=audit_level)
    support = aggregate_candidate_support(structures, fragment_rows, ms1_rows, ms2_rows, blocked, audit_level=audit_level)
    comparison = compare_legacy_composite(support, phase1_sheets.get("Composite_Mod_Candidates", []),
        formal_ranking, audit_level=audit_level)
    scores = simulate_shadow_scores(support, comparison, formal_ranking, audit_level=audit_level)
    matched_blocked = sum(r.get("Observed_mz") not in ("", None) for r in blocked)
    summary = [{
        "Schema_Version": loaded.schema_version, "Enabled": loaded.enabled,
        "Valid_Hypothesis_Count": len(structures), "Invalid_Hypothesis_Count": len(invalid),
        "Fragment_Mass_Count": len(fragment_rows),
        "MS1_Match_Count": sum(r["Match_Status"] == "matched" for r in ms1_rows),
        "MS2_Ion_Count": len(ions), "MS2_Match_Count": len(ms2_rows),
        "Blocked_Cleavage_Match_Count": matched_blocked, "Formal_Result_Changed": False,
        "Remaining_Risk": "Shadow hypotheses require curated structural and orthogonal evidence before formal propagation.",
        "Audit_Level": audit_level, "Applied_To_Formal_Result": False, "Formal_Change_Ready": False,
    }]
    sheets = {
        "Composite_MS1_Summary": _ms1_summary(structures, ms1_rows, audit_level),
        "Composite_Support_Summary": support, "Legacy_Composite_Compare": comparison,
        "Composite_Shadow_Score": scores, "Composite_Obs_Summary": summary,
    }
    if audit_level == "full":
        sheets.update({"Composite_Fragment_Masses": fragment_rows, "Composite_MS1_Matches": ms1_rows,
            "Composite_MS2_Ions": ions, "Composite_MS2_Matches": ms2_rows,
            "Composite_MS2_Assignment_Competition": ms2_competition,
            "Blocked_Cleavage_Matches": blocked})
        if invalid: sheets["Composite_Obs_Invalid"] = invalid
    return CompositeObservationResult(sheets, tuple(structures), tuple(invalid))