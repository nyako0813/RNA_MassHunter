"""Generic measured-feature audit for P1+AP dinucleotide groups."""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any
import math
import time

from rna_masshunter.mzml_diagnostics import _rt_minutes
from rna_masshunter.mzml_reader import iter_spectra
from rna_masshunter.p1_sap_dinucleotide_candidates import (
    FORMAL_FALSE, LOCALIZATION_FALSE, MS2_MODEL_REASON, dinucleotide_settings,
)

SPECPEAK_COLUMNS = [
    "Spectrum_Peak_ID", "Dinucleotide_Group_ID", "Dinucleotide_Feature_ID", "Physical_Feature_ID",
    "Spectrum_ID", "RT", "Charge", "Local_Profile_Point_Count", "Local_Apex_mz",
    "Local_Centroid_mz", "Local_mz_SD", "Local_Apex_Intensity", "Local_Integrated_Intensity",
    "Mass_Error_ppm_Apex", "Mass_Error_ppm_Centroid", "Shared_Physical_Local_Peak",
    "Shared_Dinucleotide_Group_IDs", "Applied_To_Formal_Result", "Formal_Change_Ready",
    "Formal_Result_Changed",
]
FEATURE_COLUMNS = [
    "Dinucleotide_Feature_ID", "Dinucleotide_Group_ID", "Physical_Feature_ID", "Charge",
    "RT_Start", "RT_End", "RT_Apex", "RT_Span", "Unique_Spectrum_Count",
    "Profile_Point_Count", "Observed_Apex_mz", "Observed_Centroid_mz",
    "Mass_Error_ppm_Apex", "Mass_Error_ppm_Centroid", "Mass_Accuracy_Class_Apex",
    "Mass_Accuracy_Class_Centroid", "Mass_Accuracy_Class", "Mass_Accuracy_Support_Status",
    "Mass_Accuracy_Reference", "Apex_Intensity", "Integrated_Intensity", "Peak_Symmetry",
    "Peak_Tailing_Factor", "Scan_Continuity", "Maximum_Scan_Gap", "Local_Maximum_Count",
    "Apex_Local_Contrast", "Baseline_Intensity", "Local_Noise", "Signal_to_Noise",
    "Apex_to_Baseline_Ratio", "Outside_Feature_Intensity_Ratio", "Run_RT_Coverage_Fraction",
    "Matched_Spectrum_Fraction", "Background_Persistence", "Background_Status",
    "Feature_Quality_Status", "Feature_Eligible_For_Support", "Competition_Group_ID",
    "Candidate_Specific", "Linkage_Specific", "Composition_Specific", "Structure_Specific",
    "Sequence_Position_Localized", "Original_Bond_Localized", "Position_Localization_Status",
    "Source_Bond_Resolution_Status", "Applied_To_Formal_Result", "Formal_Change_Ready",
    "Formal_Result_Changed",
]
ISOTOPE_COLUMNS = [
    "Dinucleotide_Group_ID", "Dinucleotide_Feature_ID", "Apex_Spectrum_ID", "Envelope_Assessed",
    "Envelope_Status", "Observed_Mz_M", "Observed_Mz_Mplus1", "Observed_Mz_Mplus2",
    "Observed_Intensity_M", "Observed_Intensity_Mplus1", "Observed_Intensity_Mplus2",
    "Expected_Relative_Mplus1", "Expected_Relative_Mplus2", "Observed_Relative_Mplus1",
    "Observed_Relative_Mplus2", "Envelope_Cosine_Similarity", "Envelope_Point_Count",
    "Sulfur_Count", "Sulfur_Isotope_Contribution_Assessed", "Expected_Mplus2_Contribution",
    "Observed_Mplus2_Contribution", "Sulfur_Envelope_Status", "Isotope_Peak_Shared",
    "Shared_Peak_Group_Count", "Shared_Peak_Group_IDs", "Envelope_Confounded",
    "Isomer_Isotope_Indistinguishable", "Applied_To_Formal_Result", "Formal_Change_Ready",
    "Formal_Result_Changed",
]
COMPETITION_COLUMNS = [
    "Competition_Group_ID", "Physical_Feature_ID", "Dinucleotide_Group_ID",
    "Competing_Dinucleotide_Group_Count", "Competing_Dinucleotide_Group_IDs",
    "Competing_Compositions", "Competing_Linkage_States", "Competing_Structural_Assignment_Count",
    "Competition_Types", "Candidate_Specific", "Linkage_Specific", "Composition_Specific",
    "Structure_Specific", "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]
MS2_COLUMNS = [
    "Dinucleotide_Group_ID", "Dinucleotide_Feature_ID", "Precursor_Compatible_MS2_Spectrum_Count",
    "MS2_Spectrum_IDs", "MS2_RT", "Isolation_Target_mz", "Isolation_Lower_Offset",
    "Isolation_Upper_Offset", "Precursor_Error_ppm", "MS2_Product_Min_mz", "MS2_Product_Max_mz",
    "MS2_Product_Count", "MS2_Product_Count_Below_500", "MS2_Product_Count_At_Or_Above_500",
    "MS2_Model_Applicable", "MS2_Model_Reason", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]
RAW_GROUP_COLUMNS = [
    "Dinucleotide_Group_ID", "Search_Executed", "Raw_Profile_Point_Count",
    "Unique_MS1_Spectrum_Count", "Observed_Min_mz", "Observed_Max_mz", "Raw_RT_Start", "Raw_RT_End",
]


@dataclass
class DinucleotideFeatureAuditResult:
    raw_group_rows: list[dict[str, Any]]
    spectrum_peaks: list[dict[str, Any]]
    features: list[dict[str, Any]]
    isotopes: list[dict[str, Any]]
    competition: list[dict[str, Any]]
    ms2_provenance: list[dict[str, Any]]
    performance: dict[str, float]


def classify_mass_accuracy(error_ppm: Any, strong_ppm: float, moderate_ppm: float, search_ppm: float) -> str:
    try:
        error = abs(float(error_ppm))
    except (TypeError, ValueError):
        return "NOT_EVALUABLE"
    if not math.isfinite(error): return "NOT_EVALUABLE"
    if error <= strong_ppm: return "WITHIN_STRONG_TOLERANCE"
    if error <= moderate_ppm: return "WITHIN_MODERATE_TOLERANCE"
    if error <= search_ppm: return "WITHIN_SEARCH_TOLERANCE"
    return "OUTSIDE_SEARCH_TOLERANCE"


def _raw_peak_rows(peaks: list[Any]) -> list[dict[str, Any]]:
    from rna_masshunter.p1_sap_chemical_state_audit import _raw_peak
    return [_raw_peak(peak, index) for index, peak in enumerate(peaks)]


def match_raw_groups(groups: list[dict[str, Any]], raw_peaks: list[dict[str, Any]], tolerance_ppm: float) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]]]:
    """Binary-search candidate windows; never construct group×point products."""
    ordered = sorted(raw_peaks, key=lambda row: float(row["mz"]))
    mzs = [float(row["mz"]) for row in ordered]
    matched: dict[str, list[dict[str, Any]]] = {}
    rows = []
    for group in groups:
        group_id = str(group["Dinucleotide_Group_ID"])
        executed = bool(group.get("Search_Enabled"))
        points: list[dict[str, Any]] = []
        if executed:
            target = float(group["Theoretical_mz"])
            delta = abs(target) * tolerance_ppm / 1e6
            points = ordered[bisect_left(mzs, target - delta):bisect_right(mzs, target + delta)]
        matched[group_id] = points
        rts = [float(point["rt"]) for point in points if point.get("rt") is not None]
        rows.append({
            "Dinucleotide_Group_ID": group_id, "Search_Executed": executed,
            "Raw_Profile_Point_Count": len(points),
            "Unique_MS1_Spectrum_Count": len({str(point.get("scan_id") or "") for point in points}),
            "Observed_Min_mz": min((float(point["mz"]) for point in points), default=""),
            "Observed_Max_mz": max((float(point["mz"]) for point in points), default=""),
            "Raw_RT_Start": min(rts, default=""), "Raw_RT_End": max(rts, default=""),
        })
        group["Search_Executed"] = executed
    return rows, matched


def _quality_proxy(config: Any, settings: dict[str, Any]) -> Any:
    p1 = dict(getattr(config, "p1_annotation", {}) or {})
    nested = dict(p1.get("feature_quality") or {})
    mapping = {
        "min_spectrum_count": "min_spectrum_count", "min_profile_point_count": "min_profile_point_count",
        "max_rt_gap_min": "max_refined_rt_gap_min", "background_window_rt_min": "background_window_rt_min",
        "background_mz_tolerance_ppm": "background_mz_tolerance_ppm",
    }
    for source, destination in mapping.items():
        if source in settings["feature_quality"]:
            nested[destination] = settings["feature_quality"][source]
    p1["feature_quality"] = nested
    return SimpleNamespace(p1_annotation=p1)


def _mass_status(apex: str, centroid: str) -> tuple[str, str]:
    if "OUTSIDE" in {apex, centroid}: return "OUTSIDE_SEARCH_TOLERANCE", "MASS_INCOMPATIBLE"
    order = ["WITHIN_STRONG_TOLERANCE", "WITHIN_MODERATE_TOLERANCE", "WITHIN_SEARCH_TOLERANCE"]
    for value in reversed(order):
        if value in {apex, centroid}: return value, "MASS_COMPATIBLE"
    return "NOT_EVALUABLE", "NOT_EVALUABLE"


def _background_class(row: dict[str, Any], group_feature_count: int, min_spectra: int, min_points: int) -> str:
    if int(row.get("Profile_Point_Count") or 0) <= 1: return "SINGLE_POINT_EVENT"
    if int(row.get("Profile_Point_Count") or 0) < min_points or int(row.get("Spectrum_Count") or 0) < min_spectra: return "INSUFFICIENT_POINTS"
    if str(row.get("Background_Status")) == "PERSISTENT_BACKGROUND_TRACE": return "PERSISTENT_BACKGROUND_TRACE"
    if group_feature_count > 1: return "MULTIPLE_DISCONNECTED_FEATURES"
    return "LOCALIZED_CHROMATOGRAPHIC_FEATURE"


def _map_spectrum_peaks(rows: list[dict[str, Any]], groups: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    by_physical: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        by_physical[str(row.get("Physical_Feature_ID"))].add(str(row.get("Candidate_ID")))
    output = []
    for row in rows:
        group_id = str(row.get("Candidate_ID")); theoretical = float(groups[group_id]["Theoretical_mz"])
        shared_ids = sorted(by_physical[str(row.get("Physical_Feature_ID"))] - {group_id})
        output.append({
            "Spectrum_Peak_ID": row.get("Spectrum_Peak_ID"), "Dinucleotide_Group_ID": group_id,
            "Dinucleotide_Feature_ID": row.get("Feature_ID"), "Physical_Feature_ID": row.get("Physical_Feature_ID"),
            "Spectrum_ID": row.get("Spectrum_ID"), "RT": row.get("RT"), "Charge": row.get("Charge"),
            "Local_Profile_Point_Count": row.get("Local_Profile_Point_Count"), "Local_Apex_mz": row.get("Local_Apex_mz"),
            "Local_Centroid_mz": row.get("Local_Centroid_mz"), "Local_mz_SD": row.get("Local_mz_SD"),
            "Local_Apex_Intensity": row.get("Local_Apex_Intensity"), "Local_Integrated_Intensity": row.get("Local_Integrated_Intensity"),
            "Mass_Error_ppm_Apex": (float(row["Local_Apex_mz"]) - theoretical) / theoretical * 1e6,
            "Mass_Error_ppm_Centroid": (float(row["Local_Centroid_mz"]) - theoretical) / theoretical * 1e6,
            "Shared_Physical_Local_Peak": bool(shared_ids), "Shared_Dinucleotide_Group_IDs": ";".join(shared_ids),
            **FORMAL_FALSE,
        })
    return output


def _competition_rows(features: list[dict[str, Any]], groups: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    by_physical: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for feature in features: by_physical[str(feature.get("Physical_Feature_ID"))].append(feature)
    output = []
    for serial, (physical_id, members) in enumerate(sorted(by_physical.items()), 1):
        ids = sorted({str(row["Dinucleotide_Group_ID"]) for row in members})
        competition_id = f"P1SAP_DC_{serial:06d}"
        for feature in members:
            group_id = str(feature["Dinucleotide_Group_ID"]); other_ids = [value for value in ids if value != group_id]
            all_groups = [groups[value] for value in ids]
            compositions = sorted({str(row["Final_Elemental_Composition"]) for row in all_groups})
            linkages = sorted({str(row["Linkage_State"]) for row in all_groups})
            types = []
            if len(compositions) == 1 and sum(int(row["Structural_Assignment_Count"]) for row in all_groups) > 1: types.append("SAME_COMPOSITION_STRUCTURAL_ISOMERS")
            if len(compositions) > 1: types.append("DIFFERENT_COMPOSITION_WITHIN_TOLERANCE")
            if len(linkages) > 1: types.append("NORMAL_PHOSPHATE_VS_PT_COMPETITION")
            candidate_specific = not other_ids
            linkage_specific = len(linkages) == 1
            composition_specific = len(compositions) == 1
            structure_specific = candidate_specific and int(groups[group_id]["Structural_Assignment_Count"]) == 1
            row = {
                "Competition_Group_ID": competition_id, "Physical_Feature_ID": physical_id,
                "Dinucleotide_Group_ID": group_id, "Competing_Dinucleotide_Group_Count": len(other_ids),
                "Competing_Dinucleotide_Group_IDs": ";".join(other_ids),
                "Competing_Compositions": ";".join(compositions), "Competing_Linkage_States": ";".join(linkages),
                "Competing_Structural_Assignment_Count": sum(int(groups[value]["Structural_Assignment_Count"]) for value in other_ids),
                "Competition_Types": ";".join(types), "Candidate_Specific": candidate_specific,
                "Linkage_Specific": linkage_specific, "Composition_Specific": composition_specific,
                "Structure_Specific": structure_specific, **FORMAL_FALSE,
            }
            output.append(row)
            feature.update({key: row[key] for key in ("Competition_Group_ID", "Candidate_Specific", "Linkage_Specific", "Composition_Specific", "Structure_Specific")})
    return output


def _map_features(legacy: list[dict[str, Any]], quality_rows: list[dict[str, Any]], groups: dict[str, dict[str, Any]], settings: dict[str, Any]) -> list[dict[str, Any]]:
    quality = {str(row.get("Feature_ID")): row for row in quality_rows}
    counts = Counter(str(row.get("Chemical_State_ID")) for row in legacy)
    output = []
    min_spectra = int(settings["feature_quality"].get("min_spectrum_count", 2))
    min_points = int(settings["feature_quality"].get("min_profile_point_count", 2))
    for row in legacy:
        q = quality.get(str(row.get("Feature_ID")), {}); group_id = str(row["Chemical_State_ID"])
        apex_class = classify_mass_accuracy(row.get("Mass_Error_ppm_at_Apex"), settings["strong_mass_accuracy_ppm"], settings["moderate_mass_accuracy_ppm"], settings["search_mass_accuracy_ppm"])
        centroid_class = classify_mass_accuracy(row.get("Mass_Error_ppm_at_Centroid"), settings["strong_mass_accuracy_ppm"], settings["moderate_mass_accuracy_ppm"], settings["search_mass_accuracy_ppm"])
        mass_class, mass_support = _mass_status(apex_class, centroid_class)
        background = ("PERSISTENT_BACKGROUND_TRACE" if row.get("Feature_Continuity_Status") == "continuous_background_trace" else _background_class(q, counts[group_id], min_spectra, min_points))
        output.append({
            "Dinucleotide_Feature_ID": row.get("Feature_ID"), "Dinucleotide_Group_ID": group_id,
            "Physical_Feature_ID": row.get("Physical_Feature_ID"), "Charge": row.get("Charge"),
            "RT_Start": row.get("RT_Start"), "RT_End": row.get("RT_End"), "RT_Apex": row.get("RT_Apex"), "RT_Span": row.get("RT_Span"),
            "Unique_Spectrum_Count": row.get("Spectrum_Count"), "Profile_Point_Count": row.get("Profile_Point_Count"),
            "Observed_Apex_mz": row.get("Apex_mz"), "Observed_Centroid_mz": row.get("mz_Centroid"),
            "Mass_Error_ppm_Apex": row.get("Mass_Error_ppm_at_Apex"), "Mass_Error_ppm_Centroid": row.get("Mass_Error_ppm_at_Centroid"),
            "Mass_Accuracy_Class_Apex": apex_class, "Mass_Accuracy_Class_Centroid": centroid_class,
            "Mass_Accuracy_Class": mass_class, "Mass_Accuracy_Support_Status": mass_support,
            "Mass_Accuracy_Reference": "THEORETICAL_GROUP_MZ_UNADJUSTED", "Apex_Intensity": row.get("Apex_Intensity"),
            "Integrated_Intensity": row.get("Integrated_Intensity"), "Peak_Symmetry": q.get("Peak_Symmetry", 0.0),
            "Peak_Tailing_Factor": q.get("Peak_Tailing_Factor", 0.0), "Scan_Continuity": q.get("Scan_Continuity", 0.0),
            "Maximum_Scan_Gap": q.get("Maximum_Scan_Gap", 0.0), "Local_Maximum_Count": q.get("Local_Maximum_Count", 0),
            "Apex_Local_Contrast": q.get("Apex_Local_Contrast", 0.0), "Baseline_Intensity": q.get("Baseline_Intensity", 0.0),
            "Local_Noise": q.get("Local_Noise", 0.0), "Signal_to_Noise": q.get("Signal_to_Noise", 0.0),
            "Apex_to_Baseline_Ratio": q.get("Apex_to_Baseline_Ratio", 0.0),
            "Outside_Feature_Intensity_Ratio": q.get("Outside_Feature_Intensity_Ratio", 0.0),
            "Run_RT_Coverage_Fraction": q.get("Run_RT_Coverage_Fraction", 0.0), "Matched_Spectrum_Fraction": q.get("Matched_Spectrum_Fraction", 0.0),
            "Background_Persistence": q.get("Background_Persistence", 0.0), "Background_Status": background,
            "Feature_Quality_Status": q.get("Feature_Quality_Status", "NOT_EVALUABLE"),
            "Feature_Eligible_For_Support": bool(q.get("Feature_Eligible_For_Support", False)),
            "Competition_Group_ID": "", "Candidate_Specific": False, "Linkage_Specific": False,
            "Composition_Specific": False, "Structure_Specific": False, **LOCALIZATION_FALSE, **FORMAL_FALSE,
        })
    return output


def _map_isotopes(rows: list[dict[str, Any]], features: list[dict[str, Any]], groups: dict[str, dict[str, Any]], settings: dict[str, Any]) -> list[dict[str, Any]]:
    feature_map = {str(row["Dinucleotide_Feature_ID"]): row for row in features}
    searchable = sorted((float(row["Theoretical_mz"]), str(row["Dinucleotide_Group_ID"]), int(row["Charge"])) for row in groups.values())
    mono_mzs = [row[0] for row in searchable]
    output = []
    for row in rows:
        feature_id = str(row.get("Feature_ID")); feature = feature_map.get(feature_id, {})
        group_id = str(row.get("Candidate_ID") or feature.get("Dinucleotide_Group_ID")); group = groups.get(group_id, {})
        mono = float(group.get("Theoretical_mz") or 0.0); charge = abs(int(group.get("Charge") or 1)); spacing = 1.00335483507 / charge
        shared_ids: set[str] = set()
        for target in (mono + spacing, mono + 2 * spacing):
            delta = target * settings["isotope_tolerance_ppm"] / 1e6
            for _mz, other_id, _charge in searchable[bisect_left(mono_mzs, target-delta):bisect_right(mono_mzs, target+delta)]:
                if other_id != group_id: shared_ids.add(other_id)
        envelope_status = str(row.get("Envelope_Status") or "NOT_ASSESSED")
        confounded = bool(shared_ids or row.get("Isotope_Peak_Shared_With_Other_Candidate"))
        if confounded and envelope_status not in {"MODEL_NOT_DEFINED", "MODEL_NOT_APPLICABLE", "ENVELOPE_TOO_WEAK"}:
            envelope_status = "ENVELOPE_CONFOUNDED"
        mono_intensity = float(row.get("Observed_Intensity_M") or 0.0)
        observed1 = float(row.get("Observed_Intensity_Mplus1") or 0.0) / mono_intensity if mono_intensity else 0.0
        observed2 = float(row.get("Observed_Intensity_Mplus2") or 0.0) / mono_intensity if mono_intensity else 0.0
        sulfur_count = int(row.get("Sulfur_Count") or 0)
        sulfur_status = str(row.get("Sulfur_Envelope_Compatible") or "NOT_ASSESSED")
        if not sulfur_count: sulfur_status = "NOT_APPLICABLE"
        elif confounded: sulfur_status = "ENVELOPE_CONFOUNDED"
        output.append({
            "Dinucleotide_Group_ID": group_id, "Dinucleotide_Feature_ID": feature_id,
            "Apex_Spectrum_ID": row.get("Envelope_Spectrum_ID"), "Envelope_Assessed": bool(row.get("Envelope_Assessed")),
            "Envelope_Status": envelope_status, "Observed_Mz_M": row.get("Observed_Mz_M"),
            "Observed_Mz_Mplus1": row.get("Observed_Mz_Mplus1"), "Observed_Mz_Mplus2": row.get("Observed_Mz_Mplus2"),
            "Observed_Intensity_M": mono_intensity, "Observed_Intensity_Mplus1": row.get("Observed_Intensity_Mplus1"),
            "Observed_Intensity_Mplus2": row.get("Observed_Intensity_Mplus2"),
            "Expected_Relative_Mplus1": row.get("Expected_Relative_Intensity_Mplus1"),
            "Expected_Relative_Mplus2": row.get("Expected_Relative_Intensity_Mplus2"),
            "Observed_Relative_Mplus1": observed1, "Observed_Relative_Mplus2": observed2,
            "Envelope_Cosine_Similarity": row.get("Envelope_Cosine_Similarity"), "Envelope_Point_Count": row.get("Envelope_Point_Count"),
            "Sulfur_Count": sulfur_count, "Sulfur_Isotope_Contribution_Assessed": bool(row.get("Sulfur_Isotope_Contribution_Assessed")),
            "Expected_Mplus2_Contribution": row.get("Mplus2_Expected_Contribution"), "Observed_Mplus2_Contribution": row.get("Mplus2_Observed_Contribution"),
            "Sulfur_Envelope_Status": sulfur_status, "Isotope_Peak_Shared": confounded,
            "Shared_Peak_Group_Count": len(shared_ids), "Shared_Peak_Group_IDs": ";".join(sorted(shared_ids)),
            "Envelope_Confounded": confounded, "Isomer_Isotope_Indistinguishable": bool(row.get("Isomer_Isotope_Indistinguishable")),
            **FORMAL_FALSE,
        })
    return output


def _precursor_metadata(spectrum: dict[str, Any]) -> tuple[float | None, str, Any, Any]:
    precursors = ((spectrum.get("precursorList") or {}).get("precursor") or [])
    if not precursors: return None, "", "", ""
    precursor = precursors[0]; ions = ((precursor.get("selectedIonList") or {}).get("selectedIon") or [])
    ion = ions[0] if ions else {}; isolation = precursor.get("isolationWindow") or {}
    value = ion.get("selected ion m/z", isolation.get("isolation window target m/z"))
    try: observed = float(value)
    except (TypeError, ValueError): observed = None
    target = isolation.get("isolation window target m/z", value)
    return observed, target, isolation.get("isolation window lower offset", ""), isolation.get("isolation window upper offset", "")


def build_ms2_provenance(mzml_path: str | Path | None, features: list[dict[str, Any]], groups: dict[str, dict[str, Any]], tolerance_ppm: float) -> list[dict[str, Any]]:
    if not mzml_path or not features: return []
    records: dict[str, list[dict[str, Any]]] = defaultdict(list)
    indexed_features = sorted((float(groups[str(feature["Dinucleotide_Group_ID"])]["Theoretical_mz"]), feature) for feature in features)
    feature_mzs = [row[0] for row in indexed_features]
    for spectrum in iter_spectra(mzml_path):
        if int(spectrum.get("ms level", 0) or 0) != 2: continue
        precursor, isolation_target, lower, upper = _precursor_metadata(spectrum)
        if precursor is None: continue
        product_array = spectrum.get("m/z array")
        products = [float(value) for value in product_array] if product_array is not None else []
        rt = _rt_minutes(spectrum); delta = precursor * tolerance_ppm / 1e6
        for theoretical, feature in indexed_features[bisect_left(feature_mzs, precursor-delta):bisect_right(feature_mzs, precursor+delta)]:
            if feature.get("RT_Start") is not None and rt is not None and not (float(feature["RT_Start"]) - 0.08 <= float(rt) <= float(feature["RT_End"]) + 0.08): continue
            records[str(feature["Dinucleotide_Feature_ID"])].append({
                "id": str(spectrum.get("id") or ""), "rt": rt, "target": isolation_target,
                "lower": lower, "upper": upper, "error": (precursor-theoretical)/theoretical*1e6,
                "min": min(products, default=""), "max": max(products, default=""), "count": len(products),
                "below": sum(value < 500 for value in products), "above": sum(value >= 500 for value in products),
            })
    output = []
    for feature in features:
        feature_id = str(feature["Dinucleotide_Feature_ID"]); rows = records.get(feature_id, [])
        output.append({
            "Dinucleotide_Group_ID": feature["Dinucleotide_Group_ID"], "Dinucleotide_Feature_ID": feature_id,
            "Precursor_Compatible_MS2_Spectrum_Count": len(rows), "MS2_Spectrum_IDs": ";".join(row["id"] for row in rows),
            "MS2_RT": ";".join(str(row["rt"]) for row in rows), "Isolation_Target_mz": ";".join(str(row["target"]) for row in rows),
            "Isolation_Lower_Offset": ";".join(str(row["lower"]) for row in rows), "Isolation_Upper_Offset": ";".join(str(row["upper"]) for row in rows),
            "Precursor_Error_ppm": ";".join(str(row["error"]) for row in rows), "MS2_Product_Min_mz": ";".join(str(row["min"]) for row in rows),
            "MS2_Product_Max_mz": ";".join(str(row["max"]) for row in rows), "MS2_Product_Count": sum(row["count"] for row in rows),
            "MS2_Product_Count_Below_500": sum(row["below"] for row in rows), "MS2_Product_Count_At_Or_Above_500": sum(row["above"] for row in rows),
            "MS2_Model_Applicable": False, "MS2_Model_Reason": MS2_MODEL_REASON, **FORMAL_FALSE,
        })
    return output


def audit_dinucleotide_features(groups: list[dict[str, Any]], peaks: list[Any], config: Any, *, mzml_path: str | Path | None = None) -> DinucleotideFeatureAuditResult:
    from rna_masshunter.p1_sap_chemical_state_audit import match_and_group_features
    from rna_masshunter.p1_sap_feature_quality import build_p1_sap_feature_quality
    settings = dinucleotide_settings(config); total_started = time.perf_counter()
    raw_peaks = _raw_peak_rows(peaks); group_map = {str(row["Dinucleotide_Group_ID"]): row for row in groups}
    started = time.perf_counter(); raw_rows, _matches = match_raw_groups(groups, raw_peaks, settings["search_tolerance_ppm"]); raw_runtime = time.perf_counter()-started
    started = time.perf_counter(); legacy_features, _legacy_competition, _raw_counts = match_and_group_features(groups, peaks, tolerance_ppm=settings["search_tolerance_ppm"]); chromatographic_runtime = time.perf_counter()-started
    started = time.perf_counter(); quality = build_p1_sap_feature_quality(groups, legacy_features, raw_peaks, _quality_proxy(config, settings), settings["search_tolerance_ppm"]); spectrum_runtime = time.perf_counter()-started
    spectrum_peaks = _map_spectrum_peaks(quality["spectrum_peaks"], group_map)
    features = _map_features(legacy_features, quality["quality_rows"], group_map, settings)
    started = time.perf_counter(); competition = _competition_rows(features, group_map); competition_runtime = time.perf_counter()-started
    started = time.perf_counter(); isotopes = _map_isotopes(quality["isotope_rows"], features, group_map, settings) if settings["isotope_enabled"] else []; isotope_runtime = time.perf_counter()-started
    started = time.perf_counter(); ms2 = build_ms2_provenance(mzml_path, features, group_map, settings["search_tolerance_ppm"]) if settings["ms2_provenance_enabled"] else []; ms2_runtime = time.perf_counter()-started
    performance = {
        "Raw_Matching_Runtime": raw_runtime, "Spectrum_Level_Grouping_Runtime": spectrum_runtime,
        "Chromatographic_Grouping_Runtime": chromatographic_runtime, "Isotope_Runtime": isotope_runtime,
        "Competition_Runtime": competition_runtime, "MS2_Provenance_Runtime": ms2_runtime,
        "Feature_Audit_Runtime": time.perf_counter()-total_started,
    }
    return DinucleotideFeatureAuditResult(raw_rows, spectrum_peaks, features, isotopes, competition, ms2, performance)
