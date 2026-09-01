"""Generic evidence interpretation and output assembly for P1+AP dinucleotides."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any
import json
import time
import tracemalloc

from rna_masshunter.p1_sap_dinucleotide_candidates import (
    ASSIGNMENT_COLUMNS, DINUCLEOTIDE_MODEL_VERSION, FORMAL_FALSE, GROUP_COLUMNS,
    LOCALIZATION_FALSE, SUMMARY_COLUMNS, dinucleotide_settings, generate_dinucleotide_candidates,
)
from rna_masshunter.p1_sap_dinucleotide_feature_audit import (
    COMPETITION_COLUMNS, FEATURE_COLUMNS, ISOTOPE_COLUMNS, MS2_COLUMNS, SPECPEAK_COLUMNS,
    DinucleotideFeatureAuditResult, audit_dinucleotide_features,
)
from rna_masshunter.p1_sap_dinucleotide_evidence_synthesis import (
    build_p1_sap_dinucleotide_evidence_synthesis,
)
from rna_masshunter.resource_utils import get_maximum_rss_mib

TARGET_COLUMNS = [
    "Target_Label", "Target_mz", "Tolerance_ppm", "Matched_Group_Count", "Matched_Group_IDs",
    "Matched_Feature_Count", "Matched_Feature_IDs", "Qualified_Feature_Count",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]


def _isotope_by_feature(audit: DinucleotideFeatureAuditResult) -> dict[str, dict[str, Any]]:
    return {str(row["Dinucleotide_Feature_ID"]): row for row in audit.isotopes}


def _competition_by_feature(audit: DinucleotideFeatureAuditResult) -> dict[tuple[str, str], dict[str, Any]]:
    return {(str(row["Physical_Feature_ID"]), str(row["Dinucleotide_Group_ID"])): row for row in audit.competition}


def classify_features(audit: DinucleotideFeatureAuditResult, groups: list[dict[str, Any]], config: Any) -> None:
    settings = dinucleotide_settings(config); group_map = {str(row["Dinucleotide_Group_ID"]): row for row in groups}
    isotopes = _isotope_by_feature(audit); competitions = _competition_by_feature(audit)
    min_spectra = int(settings["feature_quality"].get("min_spectrum_count", 2)); min_points = int(settings["feature_quality"].get("min_profile_point_count", 2))
    for feature in audit.features:
        isotope = isotopes.get(str(feature["Dinucleotide_Feature_ID"]), {})
        competition = competitions.get((str(feature["Physical_Feature_ID"]), str(feature["Dinucleotide_Group_ID"])), {})
        group = group_map[str(feature["Dinucleotide_Group_ID"])]
        spectra = int(feature.get("Unique_Spectrum_Count") or 0); points = int(feature.get("Profile_Point_Count") or 0)
        background = str(feature.get("Background_Status") or "NOT_EVALUABLE")
        envelope = str(isotope.get("Envelope_Status") or "NOT_ASSESSED")
        if not group.get("Search_Enabled"): status = "NOT_EVALUABLE"
        elif points <= 1 or spectra <= 1: status = "SINGLE_POINT_REJECTED"
        elif points < min_points or spectra < min_spectra: status = "PROFILE_ONLY_REJECTED"
        elif background == "PERSISTENT_BACKGROUND_TRACE": status = "BACKGROUND_TRACE_REJECTED"
        elif feature.get("Mass_Accuracy_Class") == "OUTSIDE_SEARCH_TOLERANCE": status = "MASS_INCOMPATIBLE"
        elif envelope == "ENVELOPE_INCOMPATIBLE": status = "ISOTOPE_INCOMPATIBLE"
        elif competition.get("Competing_Dinucleotide_Group_Count", 0):
            if not competition.get("Linkage_Specific", False): status = "QUALIFIED_BUT_LINKAGE_AMBIGUOUS"
            elif not competition.get("Composition_Specific", False): status = "QUALIFIED_BUT_COMPOSITION_AMBIGUOUS"
            else: status = "COMPETITION_UNRESOLVED"
        elif int(group.get("Structural_Assignment_Count") or 0) > 1: status = "QUALIFIED_BUT_STRUCTURALLY_AMBIGUOUS"
        else: status = "QUALIFIED_CHROMATOGRAPHIC_FEATURE"
        qualified = status.startswith("QUALIFIED")
        feature["Feature_Quality_Status"] = status
        feature["Feature_Eligible_For_Support"] = qualified


def interpret_groups(groups: list[dict[str, Any]], audit: DinucleotideFeatureAuditResult) -> None:
    raw_map = {str(row["Dinucleotide_Group_ID"]): row for row in audit.raw_group_rows}
    by_group: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for feature in audit.features: by_group[str(feature["Dinucleotide_Group_ID"])].append(feature)
    competition_by_group: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in audit.competition: competition_by_group[str(row["Dinucleotide_Group_ID"])].append(row)
    for group in groups:
        group_id = str(group["Dinucleotide_Group_ID"]); raw = raw_map.get(group_id, {}); features = by_group.get(group_id, [])
        qualified = [row for row in features if str(row.get("Feature_Quality_Status", "")).startswith("QUALIFIED")]
        group.update({key: raw.get(key, "") for key in ("Raw_Profile_Point_Count", "Unique_MS1_Spectrum_Count", "Observed_Min_mz", "Observed_Max_mz", "Raw_RT_Start", "Raw_RT_End")})
        group["Feature_Count"] = len(features); group["Qualified_Feature_Count"] = len(qualified)
        competitions = competition_by_group.get(group_id, [])
        composition_ambiguous = any(not row.get("Composition_Specific", False) and row.get("Competing_Dinucleotide_Group_Count", 0) for row in competitions)
        linkage_ambiguous = any(not row.get("Linkage_Specific", False) and row.get("Competing_Dinucleotide_Group_Count", 0) for row in competitions)
        if not group.get("Search_Enabled"): interpretation = "NOT_EVALUABLE"
        elif not int(raw.get("Raw_Profile_Point_Count") or 0): interpretation = "NO_MATCH"
        elif not features: interpretation = "RAW_MATCH_ONLY"
        elif not qualified: interpretation = "NO_QUALIFIED_FEATURE"
        elif linkage_ambiguous: interpretation = "QUALIFIED_LINKAGE_AMBIGUOUS_FEATURE"
        elif composition_ambiguous: interpretation = "QUALIFIED_COMPOSITION_COMPATIBLE_FEATURE"
        elif int(group["Structural_Assignment_Count"]) > 1: interpretation = "QUALIFIED_STRUCTURE_UNRESOLVED_FEATURE"
        else: interpretation = "DINUCLEOTIDE_GROUP_SUPPORTED"
        group["Group_Interpretation"] = interpretation
        group["Composition_Resolution_Status"] = "COMPOSITION_AMBIGUOUS" if composition_ambiguous else ("COMPOSITION_COMPATIBLE" if qualified else "COMPOSITION_UNRESOLVED")
        group["Linkage_Resolution_Status"] = "LINKAGE_AMBIGUOUS" if linkage_ambiguous else ("LINKAGE_COMPATIBLE" if qualified else "LINKAGE_UNRESOLVED")
        group["Structure_Resolution_Status"] = "STRUCTURE_UNRESOLVED" if int(group["Structural_Assignment_Count"]) > 1 or not qualified else "STRUCTURE_GROUP_COMPATIBLE"
        group["Source_Bond_Resolution_Status"] = "SOURCE_BOND_UNRESOLVED"
        group.update(LOCALIZATION_FALSE)


def build_target_results(targets: list[dict[str, Any]], groups: list[dict[str, Any]], features: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Optional reporting filter; never mutates groups or features."""
    by_group: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for feature in features: by_group[str(feature["Dinucleotide_Group_ID"])].append(feature)
    results = []
    for number, target in enumerate(targets, 1):
        label = str(target.get("label") or f"target_{number}"); target_mz = float(target["theoretical_mz"]); tolerance = float(target.get("tolerance_ppm", 10.0))
        matched_groups = [row for row in groups if abs(float(row["Theoretical_mz"])-target_mz) <= abs(target_mz)*tolerance/1e6]
        matched_features = [feature for group in matched_groups for feature in by_group.get(str(group["Dinucleotide_Group_ID"]), [])]
        results.append({
            "Target_Label": label, "Target_mz": target_mz, "Tolerance_ppm": tolerance,
            "Matched_Group_Count": len(matched_groups), "Matched_Group_IDs": ";".join(str(row["Dinucleotide_Group_ID"]) for row in matched_groups),
            "Matched_Feature_Count": len(matched_features), "Matched_Feature_IDs": ";".join(str(row["Dinucleotide_Feature_ID"]) for row in matched_features),
            "Qualified_Feature_Count": sum(str(row.get("Feature_Quality_Status", "")).startswith("QUALIFIED") for row in matched_features),
            **FORMAL_FALSE,
        })
    return results


def _summary(candidate: Any, audit: DinucleotideFeatureAuditResult, performance: dict[str, float]) -> dict[str, Any]:
    groups = candidate.candidates; qualified_groups = [row for row in groups if int(row.get("Qualified_Feature_Count") or 0) > 0]
    isotopes = audit.isotopes
    result = dict(candidate.summary)
    result.update({
        "Searched_Group_Count": sum(bool(row.get("Search_Executed")) for row in groups),
        "Raw_Matched_Group_Count": sum(int(row.get("Raw_Profile_Point_Count") or 0) > 0 for row in groups),
        "Qualified_Group_Count": len(qualified_groups),
        "Qualified_Feature_Count": sum(str(row.get("Feature_Quality_Status", "")).startswith("QUALIFIED") for row in audit.features),
        "Qualified_Normal_Phosphate_Group_Count": sum(row["Linkage_State"] == "NORMAL_PHOSPHATE" for row in qualified_groups),
        "Qualified_PT_Group_Count": sum(row["Linkage_State"] == "PHOSPHOROTHIOATE" for row in qualified_groups),
        "Competition_Unresolved_Group_Count": sum("AMBIGUOUS" in str(row.get("Composition_Resolution_Status")) or "AMBIGUOUS" in str(row.get("Linkage_Resolution_Status")) for row in groups),
        "Isotope_Assessed_Feature_Count": sum(bool(row.get("Envelope_Assessed")) for row in isotopes),
        "Isotope_Compatible_Feature_Count": sum(row.get("Envelope_Status") == "ENVELOPE_COMPATIBLE" for row in isotopes),
        "Isotope_Incompatible_Feature_Count": sum(row.get("Envelope_Status") == "ENVELOPE_INCOMPATIBLE" for row in isotopes),
        "Precursor_Compatible_MS2_Count": sum(int(row.get("Precursor_Compatible_MS2_Spectrum_Count") or 0) for row in audit.ms2_provenance),
        **performance, **LOCALIZATION_FALSE, **FORMAL_FALSE,
    })
    return result


def build_p1_sap_dinucleotide_audit(project_root: str | Path, sequence: str, peaks: list[Any], config: Any, *, audit_level: str = "audit", mzml_path: str | Path | None = None) -> dict[str, Any]:
    settings = dinucleotide_settings(config)
    if not settings["enabled"]:
        payload = {"dinucleotide_audit": {"model_version": DINUCLEOTIDE_MODEL_VERSION, "status": "DISABLED_BY_CONFIG", "candidate_generation": {}, "group_summary": {}, "groups": [], "features": [], "isotope_audit": [], "competition": [], "ms2_provenance": [], "target_results": []}}
        return {"generated": None, "audit": None, "sheets": {}, "payload": payload, "features": [], "ms2": [], "summary": {}}
    total_started = time.perf_counter(); owns_tracemalloc = not tracemalloc.is_tracing()
    if owns_tracemalloc: tracemalloc.start()
    started = time.perf_counter(); candidate = generate_dinucleotide_candidates(sequence, project_root, config=config); candidate_runtime = time.perf_counter()-started
    audit = audit_dinucleotide_features(candidate.candidates, peaks, config, mzml_path=mzml_path)
    started = time.perf_counter(); classify_features(audit, candidate.candidates, config); interpret_groups(candidate.candidates, audit); interpretation_runtime = time.perf_counter()-started
    synthesis = build_p1_sap_dinucleotide_evidence_synthesis(
        candidate.candidates, candidate.assignments, audit.features,
        audit.competition, audit.isotopes, audit.ms2_provenance,
    )
    targets = build_target_results(dinucleotide_settings(config)["targets"], candidate.candidates, audit.features)
    _current, peak_bytes = tracemalloc.get_traced_memory()
    if owns_tracemalloc: tracemalloc.stop()
    maximum_rss = get_maximum_rss_mib()
    performance = {
        "Candidate_Generation_Runtime": candidate_runtime, **audit.performance,
        "Interpretation_Runtime": interpretation_runtime,
        "Total_Shadow_Audit_Runtime": time.perf_counter()-total_started,
        "Maximum_RSS_MiB": maximum_rss, "Tracemalloc_Peak_MiB": peak_bytes/(1024*1024),
    }
    summary = _summary(candidate, audit, performance)
    sheets = {
        "P1_SAP_Dinuc_Summary": [summary], "P1_SAP_Dinuc_Groups": candidate.candidates,
        "P1_SAP_Dinuc_Targets": targets, **synthesis.sheets,
    }
    if audit_level == "full":
        sheets.update({
            "P1_SAP_Dinuc_Assignments": candidate.assignments,
            "P1_SAP_Dinuc_SpecPeaks": audit.spectrum_peaks,
            "P1_SAP_Dinuc_Features": audit.features,
            "P1_SAP_Dinuc_Isotopes": audit.isotopes,
            "P1_SAP_Dinuc_Competition": audit.competition,
            "P1_SAP_Dinuc_MS2": audit.ms2_provenance,
        })
    payload = {
        "dinucleotide_audit": {
            "model_version": DINUCLEOTIDE_MODEL_VERSION,
            "candidate_generation": {key: candidate.summary.get(key) for key in candidate.summary},
            "group_summary": summary, "groups": candidate.candidates,
            "features": audit.features, "isotope_audit": audit.isotopes,
            "competition": audit.competition, "ms2_provenance": audit.ms2_provenance,
            "evidence_synthesis": synthesis.to_jsonable(), "target_results": targets,
        }
    }
    json.dumps(payload, ensure_ascii=False, default=str)
    return {
        "generated": candidate, "audit": audit, "synthesis": synthesis, "sheets": sheets,
        "payload": payload, "features": audit.features, "ms2": audit.ms2_provenance,
        "summary": summary,
    }
