"""Shadow-only P1+SAP feature-quality and isotope-envelope audit."""
from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from typing import Any

FEATURE_QUALITY_MODEL_VERSION = "1.1"
ISOTOPE_SPACING = 1.00335483507
FORMAL_FALSE = {
    "Applied_To_Formal_Result": False,
    "Formal_Change_Ready": False,
    "Formal_Result_Changed": False,
}

P1_SAP_SPECTRUM_PEAK_COLUMNS = [
    "Spectrum_Peak_ID", "Spectrum_ID", "RT", "Candidate_ID", "Charge",
    "Local_Profile_Point_Count", "Local_Centroid_mz", "Local_Apex_mz",
    "Local_Apex_Intensity", "Local_Integrated_Intensity", "Local_mz_SD",
    "Local_Peak_Boundary_Left", "Local_Peak_Boundary_Right",
    "Chemical_Family", "Physical_Feature_ID", "Feature_ID",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]
P1_SAP_REFINED_FEATURE_COLUMNS = [
    "Chrom_Feature_ID", "Candidate_ID", "Physical_Feature_ID", "Chemical_Family", "Charge",
    "RT_Start", "RT_End", "RT_Apex", "Spectrum_Count", "Integrated_Intensity",
    "Apex_Intensity", "Centroid_mz", "Mass_Error_ppm", "Peak_Width", "Peak_Symmetry",
    "Background_Adjusted_Area", "Legacy_Feature_ID", "Refinement_Status", "Feature_Quality_Status",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]
P1_SAP_FEATURE_QUALITY_COLUMNS = [
    "Feature_ID", "Physical_Feature_ID", "Chemical_State_ID", "Chemical_Family", "Charge",
    "RT_Start", "RT_End", "RT_Apex", "RT_Span", "Spectrum_Count", "Profile_Point_Count",
    "Apex_mz", "Centroid_mz", "Mass_Error_ppm_Apex", "Mass_Error_ppm_Centroid", "mz_SD",
    "Apex_Intensity", "Integrated_Intensity", "Median_Intensity", "Baseline_Intensity",
    "Local_Noise", "Signal_to_Noise", "Apex_to_Baseline_Ratio", "Left_Point_Count",
    "Right_Point_Count", "Left_RT_Span", "Right_RT_Span", "Peak_Symmetry", "Peak_Tailing_Factor",
    "Scan_Continuity", "Maximum_Scan_Gap", "RT_Continuity", "Run_RT_Coverage_Fraction",
    "Matched_Spectrum_Fraction", "Apex_Local_Contrast", "Outside_Feature_Intensity_Ratio",
    "Local_Maximum_Count", "Background_Persistence", "Feature_Width_Status", "Feature_Shape_Status",
    "Background_Status", "Feature_Quality_Status", "Feature_Quality_Score", "Feature_Rejection_Reasons",
    "Raw_Match_Present", "Independent_Feature_Qualified", "Feature_Eligible_For_Support",
    "Competition_Count", "Competing_Candidate_IDs", "Isomer_Isotope_Indistinguishable",
    "Mass_Accuracy_Component", "Spectrum_Continuity_Component", "Chromatographic_Localization_Component",
    "Peak_Shape_Component", "Background_Contrast_Component", "Isotope_Envelope_Component",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]
P1_SAP_ISOTOPE_AUDIT_COLUMNS = [
    "Feature_ID", "Physical_Feature_ID", "Chemical_State_ID", "Candidate_ID", "Charge",
    "Sulfur_Count", "Sulfur_Isotope_Contribution_Assessed", "Mplus2_Expected_Contribution",
    "Mplus2_Observed_Contribution", "Sulfur_Envelope_Compatible", "Envelope_Assessed",
    "Envelope_Status", "Envelope_Spectrum_ID", "Observed_Mz_M", "Observed_Intensity_M",
    "Observed_Mz_Mplus1", "Observed_Intensity_Mplus1", "Observed_Mz_Mplus2",
    "Observed_Intensity_Mplus2", "Expected_Relative_Intensity_Mplus1",
    "Expected_Relative_Intensity_Mplus2", "Relative_Intensity_Error_Mplus1",
    "Relative_Intensity_Error_Mplus2", "Envelope_Correlation", "Envelope_Cosine_Similarity",
    "Isotope_Peak_Shared_With_Other_Candidate", "Isomer_Isotope_Indistinguishable",
    "Envelope_Point_Count", "Applied_To_Formal_Result", "Formal_Change_Ready",
    "Formal_Result_Changed",
]
P1_SAP_QUALITY_SUMMARY_COLUMNS = [
    "Feature_Quality_Model_Version", "Spectrum_Level_Local_Peak_Count", "Refined_Feature_Count",
    "Qualified_Feature_Count", "Rejected_Single_Point_Count", "Rejected_Background_Count",
    "Rejected_Profile_Only_Count", "Isotope_Incompatible_Count", "Competition_Unresolved_Count",
    "Model_Not_Defined_Count", "Raw_Match_Only_Count", "Qualified_PT_Feature_Count",
    "PT_Final_Interpretation", "Feature_Quality_Score_Mean", "Feature_Quality_Score_Median",
    "Feature_Quality_Score_Min", "Feature_Quality_Score_Max",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]

DEFAULT_QUALITY_CONFIG = {
    "min_spectrum_count": 2,
    "min_profile_point_count": 2,
    "background_rt_coverage_fraction": 0.30,
    "background_spectrum_fraction": 0.80,
    "background_contrast_min": 4.0,
    "local_maximum_count_threshold": 2,
    "feature_symmetry_min": 0.4,
    "feature_tailing_max": 5.0,
    "max_refined_rt_gap_min": 0.08,
    "background_window_rt_min": 0.5,
    "background_mz_tolerance_ppm": 20.0,
    "isotope_match_tolerance_ppm": 20.0,
}
NOMINAL_ISOTOPE_PROBABILITIES = {
    "C": (0.9893, 0.0107, 0.0),
    "H": (0.999885, 0.000115, 0.0),
    "N": (0.99636, 0.00364, 0.0),
    "O": (0.99757, 0.00038, 0.00205),
    "S": (0.9499, 0.0075, 0.0425),
    "P": (1.0, 0.0, 0.0),
}


def _quality_config(config: Any) -> dict[str, Any]:
    section = getattr(config, "p1_annotation", {}) or {}
    nested = section.get("feature_quality", {}) or {}
    values: dict[str, Any] = {}
    for key, default in DEFAULT_QUALITY_CONFIG.items():
        value = nested.get(key, section.get(key, default))
        values[key] = type(default)(value)
    return values


def _parse_composition(composition: str | None) -> tuple[dict[str, int] | None, str]:
    text = str(composition or "").strip()
    if not text or text == "MODEL_NOT_DEFINED":
        return None, "NOT_ASSESSED"
    matches = list(re.finditer(r"([A-Z][a-z]?)(-?\d*)", text))
    if not matches or "".join(match.group(0) for match in matches) != text:
        return None, "MODEL_NOT_APPLICABLE"
    counts: dict[str, int] = {}
    for match in matches:
        element, token = match.groups()
        if element not in NOMINAL_ISOTOPE_PROBABILITIES:
            return None, "MODEL_NOT_APPLICABLE"
        count = int(token) if token not in {"", "-"} else 1
        if count < 0:
            return None, "MODEL_NOT_APPLICABLE"
        counts[element] = counts.get(element, 0) + count
    if not counts or all(value == 0 for value in counts.values()):
        return None, "MODEL_NOT_APPLICABLE"
    return counts, "ASSESSED"


def _expected_isotope_abundance(composition: dict[str, int] | None) -> dict[str, float]:
    if composition is None:
        return {"M+1": 0.0, "M+2": 0.0}
    distribution = [1.0, 0.0, 0.0]
    for element, count in composition.items():
        atom = NOMINAL_ISOTOPE_PROBABILITIES[element]
        for _ in range(count):
            distribution = [
                distribution[0] * atom[0],
                distribution[1] * atom[0] + distribution[0] * atom[1],
                distribution[2] * atom[0] + distribution[1] * atom[1] + distribution[0] * atom[2],
            ]
    monoisotopic = distribution[0]
    return {
        "M+1": distribution[1] / monoisotopic if monoisotopic else 0.0,
        "M+2": distribution[2] / monoisotopic if monoisotopic else 0.0,
    }


def _median(values: list[float]) -> float:
    values = sorted(values)
    if not values:
        return 0.0
    middle = len(values) // 2
    return values[middle] if len(values) % 2 else (values[middle - 1] + values[middle]) / 2.0


def _mad(values: list[float], center: float) -> float:
    return _median([abs(value - center) for value in values]) if values else 0.0


def _ppm_tolerance(mz: float, tolerance_ppm: float) -> float:
    return abs(mz) * abs(tolerance_ppm) / 1e6


def _cosine_similarity(left: list[float], right: list[float]) -> float:
    nl = math.sqrt(sum(value * value for value in left))
    nr = math.sqrt(sum(value * value for value in right))
    return sum(a * b for a, b in zip(left, right)) / (nl * nr) if nl and nr else 0.0


def _build_local_peak_groups(points: list[dict[str, Any]], tolerance_ppm: float) -> list[list[dict[str, Any]]]:
    groups: list[list[dict[str, Any]]] = []
    for point in sorted(points, key=lambda row: float(row["mz"])):
        if not groups:
            groups.append([point])
            continue
        centroid = sum(float(row["mz"]) for row in groups[-1]) / len(groups[-1])
        if abs(float(point["mz"]) - centroid) <= _ppm_tolerance(centroid, tolerance_ppm):
            groups[-1].append(point)
        else:
            groups.append([point])
    return groups


def _points_by_feature(raw_peaks: list[dict[str, Any]], features: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for feature in features:
        for index in feature.get("_point_ids", []):
            if 0 <= int(index) < len(raw_peaks):
                result[str(feature["Feature_ID"])].append(raw_peaks[int(index)])
    return result


def build_spectrum_level_peaks(
    features: list[dict[str, Any]],
    raw_peaks: list[dict[str, Any]],
    tolerance_ppm: float,
) -> list[dict[str, Any]]:
    points_by_feature = _points_by_feature(raw_peaks, features)
    rows: list[dict[str, Any]] = []
    for feature in features:
        by_scan: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for point in points_by_feature.get(str(feature["Feature_ID"]), []):
            by_scan[str(point.get("scan_id") or "")].append(point)
        for scan_id, scan_points in by_scan.items():
            for number, group in enumerate(_build_local_peak_groups(scan_points, tolerance_ppm), 1):
                integrated = sum(float(point["intensity"]) for point in group)
                centroid = (
                    sum(float(point["mz"]) * float(point["intensity"]) for point in group) / integrated
                    if integrated else sum(float(point["mz"]) for point in group) / len(group)
                )
                apex = max(group, key=lambda point: float(point["intensity"]))
                variance = (
                    sum(float(point["intensity"]) * (float(point["mz"]) - centroid) ** 2 for point in group) / integrated
                    if integrated else 0.0
                )
                rows.append({
                    "Spectrum_Peak_ID": f"SP_{feature['Feature_ID']}_{scan_id}_{number}",
                    "Spectrum_ID": scan_id,
                    "RT": float(group[0].get("rt") or 0.0),
                    "Candidate_ID": feature["Chemical_State_ID"],
                    "Charge": feature.get("Charge"),
                    "Local_Profile_Point_Count": len(group),
                    "Local_Centroid_mz": centroid,
                    "Local_Apex_mz": float(apex["mz"]),
                    "Local_Apex_Intensity": float(apex["intensity"]),
                    "Local_Integrated_Intensity": integrated,
                    "Local_mz_SD": math.sqrt(variance),
                    "Local_Peak_Boundary_Left": min(float(point["mz"]) for point in group),
                    "Local_Peak_Boundary_Right": max(float(point["mz"]) for point in group),
                    "Chemical_Family": feature.get("Chemical_Family"),
                    "Physical_Feature_ID": feature.get("Physical_Feature_ID"),
                    "Feature_ID": feature["Feature_ID"],
                    **FORMAL_FALSE,
                })
    return rows


def _build_refined_feature(
    group: list[dict[str, Any]],
    candidate: dict[str, Any],
    legacy: list[dict[str, Any]],
    serial: int,
) -> dict[str, Any]:
    group = sorted(group, key=lambda row: (float(row["RT"]), str(row["Spectrum_ID"])))
    rts = [float(row["RT"]) for row in group]
    integrated = sum(float(row["Local_Integrated_Intensity"]) for row in group)
    apex = max(group, key=lambda row: float(row["Local_Apex_Intensity"]))
    centroid = (
        sum(float(row["Local_Centroid_mz"]) * float(row["Local_Integrated_Intensity"]) for row in group) / integrated
        if integrated else float(group[0]["Local_Centroid_mz"])
    )
    legacy_ids = sorted({str(row["Feature_ID"]) for row in legacy})
    if not legacy_ids:
        refinement_status = "REJECTED_AS_PROFILE_ONLY"
        quality_status = "PROFILE_ONLY_REJECTED"
    else:
        refinement_status = "CONFIRMED_BY_REFINED_GROUPING" if len(legacy_ids) == 1 else "MERGED_BY_REFINED_GROUPING"
        quality_status = "PENDING_LEGACY_QUALITY"
    left = float(apex["RT"]) - rts[0]
    right = rts[-1] - float(apex["RT"])
    symmetry = min(left, right) / max(left, right) if left > 0 and right > 0 else 0.0
    theoretical = float(candidate.get("Theoretical_mz") or 0.0)
    return {
        "Chrom_Feature_ID": f"RF_{candidate.get('Chemical_State_ID', 'UNKNOWN')}_{serial}",
        "Candidate_ID": candidate.get("Chemical_State_ID"),
        "Physical_Feature_ID": group[0].get("Physical_Feature_ID"),
        "Chemical_Family": candidate.get("Chemical_Family"),
        "Charge": group[0].get("Charge"),
        "RT_Start": rts[0],
        "RT_End": rts[-1],
        "RT_Apex": float(apex["RT"]),
        "Spectrum_Count": len({str(row["Spectrum_ID"]) for row in group}),
        "Integrated_Intensity": integrated,
        "Apex_Intensity": float(apex["Local_Apex_Intensity"]),
        "Centroid_mz": centroid,
        "Mass_Error_ppm": (centroid - theoretical) / theoretical * 1e6 if theoretical else 0.0,
        "Peak_Width": rts[-1] - rts[0],
        "Peak_Symmetry": symmetry,
        "Background_Adjusted_Area": integrated,
        "Legacy_Feature_ID": ";".join(legacy_ids),
        "Refinement_Status": refinement_status,
        "Feature_Quality_Status": quality_status,
        **FORMAL_FALSE,
    }


def build_refined_chromatographic_features(
    spectrum_peaks: list[dict[str, Any]],
    candidates: dict[str, dict[str, Any]],
    features: list[dict[str, Any]],
    config: Any,
    tolerance_ppm: float,
) -> list[dict[str, Any]]:
    settings = _quality_config(config)
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for row in spectrum_peaks:
        grouped[(str(row["Candidate_ID"]), int(row.get("Charge") or 0))].append(row)
    result: list[dict[str, Any]] = []
    serial = 0
    for (candidate_id, _charge), rows in grouped.items():
        current: list[dict[str, Any]] = []
        for row in sorted(rows, key=lambda item: (float(item["RT"]), float(item["Local_Centroid_mz"]))):
            if current:
                previous = current[-1]
                rt_gap = float(row["RT"]) - float(previous["RT"])
                mz_gap = abs(float(row["Local_Centroid_mz"]) - float(previous["Local_Centroid_mz"]))
                mz_limit = _ppm_tolerance(float(previous["Local_Centroid_mz"]), tolerance_ppm)
                if rt_gap > settings["max_refined_rt_gap_min"] or mz_gap > mz_limit:
                    serial += 1
                    legacy = _matching_legacy_features(current, features, candidate_id)
                    result.append(_build_refined_feature(current, candidates.get(candidate_id, {}), legacy, serial))
                    current = []
            current.append(row)
        if current:
            serial += 1
            legacy = _matching_legacy_features(current, features, candidate_id)
            result.append(_build_refined_feature(current, candidates.get(candidate_id, {}), legacy, serial))
    return result


def _matching_legacy_features(
    group: list[dict[str, Any]],
    features: list[dict[str, Any]],
    candidate_id: str,
) -> list[dict[str, Any]]:
    start = min(float(row["RT"]) for row in group)
    end = max(float(row["RT"]) for row in group)
    physical_ids = {row.get("Physical_Feature_ID") for row in group}
    return [
        feature for feature in features
        if str(feature.get("Chemical_State_ID")) == candidate_id
        and feature.get("Physical_Feature_ID") in physical_ids
        and float(feature.get("RT_End") or 0.0) >= start
        and float(feature.get("RT_Start") or 0.0) <= end
    ]


def _target_mz(feature: dict[str, Any], candidate: dict[str, Any]) -> float:
    return float(
        candidate.get("Theoretical_mz")
        or feature.get("mz_Centroid")
        or feature.get("Centroid_mz")
        or feature.get("Apex_mz")
        or 0.0
    )


def _near_mz(point: dict[str, Any], target: float, ppm: float) -> bool:
    return target > 0 and abs(float(point["mz"]) - target) <= _ppm_tolerance(target, ppm)


def _background_points(
    feature: dict[str, Any],
    candidate: dict[str, Any],
    raw_points: list[dict[str, Any]],
    settings: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    target = _target_mz(feature, candidate)
    start = float(feature.get("RT_Start") or 0.0)
    end = float(feature.get("RT_End") or start)
    window = settings["background_window_rt_min"]
    inside = [
        point for point in raw_points
        if point.get("rt") is not None
        and start <= float(point["rt"]) <= end
        and _near_mz(point, target, settings["background_mz_tolerance_ppm"])
    ]
    outside = [
        point for point in raw_points
        if point.get("rt") is not None
        and ((start - window) <= float(point["rt"]) < start or end < float(point["rt"]) <= (end + window))
        and _near_mz(point, target, settings["background_mz_tolerance_ppm"])
    ]
    return inside, outside


def _detect_local_maxima(rows: list[dict[str, Any]]) -> int:
    if not rows:
        return 0
    values = [float(row["Local_Integrated_Intensity"]) for row in sorted(rows, key=lambda item: float(item["RT"]))]
    if len(values) < 3:
        return 1
    result = int(values[0] > values[1]) + int(values[-1] > values[-2])
    result += sum(values[index] > values[index - 1] and values[index] > values[index + 1] for index in range(1, len(values) - 1))
    return result


def _competition_context(
    features: list[dict[str, Any]],
    spectrum_peaks: list[dict[str, Any]],
    candidates: dict[str, dict[str, Any]],
    tolerance_ppm: float,
) -> dict[str, dict[str, Any]]:
    rows_by_feature: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in spectrum_peaks:
        rows_by_feature[str(row["Feature_ID"])].append(row)
    context: dict[str, dict[str, Any]] = {
        str(feature["Feature_ID"]): {"candidate_ids": set(), "isomer": False}
        for feature in features
    }
    for index, left in enumerate(features):
        for right in features[index + 1:]:
            left_id = str(left["Feature_ID"])
            right_id = str(right["Feature_ID"])
            if left.get("Chemical_State_ID") == right.get("Chemical_State_ID"):
                continue
            same_physical = (
                left.get("Physical_Feature_ID") not in {None, ""}
                and left.get("Physical_Feature_ID") == right.get("Physical_Feature_ID")
            )
            overlaps = False
            for lrow in rows_by_feature[left_id]:
                for rrow in rows_by_feature[right_id]:
                    if str(lrow["Spectrum_ID"]) != str(rrow["Spectrum_ID"]):
                        continue
                    left_mz = float(lrow["Local_Centroid_mz"])
                    right_mz = float(rrow["Local_Centroid_mz"])
                    limit = max(_ppm_tolerance(left_mz, tolerance_ppm), _ppm_tolerance(right_mz, tolerance_ppm))
                    boundary_overlap = (
                        float(lrow["Local_Peak_Boundary_Left"]) <= float(rrow["Local_Peak_Boundary_Right"]) + limit
                        and float(rrow["Local_Peak_Boundary_Left"]) <= float(lrow["Local_Peak_Boundary_Right"]) + limit
                    )
                    if boundary_overlap or abs(left_mz - right_mz) <= limit:
                        overlaps = True
                        break
                if overlaps:
                    break
            if not (same_physical or overlaps):
                continue
            left_candidate = str(left.get("Chemical_State_ID"))
            right_candidate = str(right.get("Chemical_State_ID"))
            context[left_id]["candidate_ids"].add(right_candidate)
            context[right_id]["candidate_ids"].add(left_candidate)
            left_comp, left_state = _parse_composition(candidates.get(left_candidate, {}).get("Elemental_Composition"))
            right_comp, right_state = _parse_composition(candidates.get(right_candidate, {}).get("Elemental_Composition"))
            isomer = left_state == right_state == "ASSESSED" and left_comp == right_comp
            context[left_id]["isomer"] = context[left_id]["isomer"] or isomer
            context[right_id]["isomer"] = context[right_id]["isomer"] or isomer
    return context


def _best_peak(
    scan_points: list[dict[str, Any]],
    target_mz: float,
    tolerance_ppm: float,
) -> dict[str, Any] | None:
    matches = [
        point for point in scan_points
        if abs(float(point["mz"]) - target_mz) <= _ppm_tolerance(target_mz, tolerance_ppm)
    ]
    return max(matches, key=lambda point: float(point["intensity"])) if matches else None


def assess_isotope_envelope(
    feature: dict[str, Any],
    raw_points: list[dict[str, Any]],
    candidate: dict[str, Any],
    tolerance_ppm: float,
) -> dict[str, Any]:
    composition, model_status = _parse_composition(candidate.get("Elemental_Composition"))
    charge = abs(int(candidate.get("Charge") or feature.get("Charge") or 1))
    base = {
        "Feature_ID": feature["Feature_ID"],
        "Physical_Feature_ID": feature.get("Physical_Feature_ID"),
        "Chemical_State_ID": feature.get("Chemical_State_ID"),
        "Candidate_ID": candidate.get("Chemical_State_ID"),
        "Charge": charge,
        "Sulfur_Count": composition.get("S", 0) if composition else 0,
        "Sulfur_Isotope_Contribution_Assessed": bool(composition and composition.get("S", 0)),
        "Mplus2_Expected_Contribution": 0.0,
        "Mplus2_Observed_Contribution": 0.0,
        "Sulfur_Envelope_Compatible": "NOT_ASSESSED",
        "Envelope_Assessed": False,
        "Envelope_Status": model_status,
        "Envelope_Spectrum_ID": "",
        "Observed_Mz_M": "",
        "Observed_Intensity_M": 0.0,
        "Observed_Mz_Mplus1": "",
        "Observed_Intensity_Mplus1": 0.0,
        "Observed_Mz_Mplus2": "",
        "Observed_Intensity_Mplus2": 0.0,
        "Expected_Relative_Intensity_Mplus1": 0.0,
        "Expected_Relative_Intensity_Mplus2": 0.0,
        "Relative_Intensity_Error_Mplus1": 0.0,
        "Relative_Intensity_Error_Mplus2": 0.0,
        "Envelope_Correlation": 0.0,
        "Envelope_Cosine_Similarity": 0.0,
        "Isotope_Peak_Shared_With_Other_Candidate": False,
        "Isomer_Isotope_Indistinguishable": False,
        "Envelope_Point_Count": 0,
        "_isotope_point_indices": set(),
        **FORMAL_FALSE,
    }
    if model_status != "ASSESSED":
        if model_status == "NOT_ASSESSED":
            base["Envelope_Status"] = "NOT_ASSESSED"
        return base

    expected = _expected_isotope_abundance(composition)
    base["Mplus2_Expected_Contribution"] = expected["M+2"]
    base["Expected_Relative_Intensity_Mplus1"] = expected["M+1"]
    base["Expected_Relative_Intensity_Mplus2"] = expected["M+2"]
    theoretical = _target_mz(feature, candidate)
    start = float(feature.get("RT_Start") or 0.0)
    end = float(feature.get("RT_End") or start)
    by_scan: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for point in raw_points:
        if point.get("rt") is not None and start <= float(point["rt"]) <= end:
            by_scan[str(point.get("scan_id") or "")].append(point)

    scan_envelopes: list[tuple[float, str, dict[str, Any] | None, dict[str, Any] | None, dict[str, Any] | None]] = []
    spacing = ISOTOPE_SPACING / charge
    for scan_id, scan_points in by_scan.items():
        mono = _best_peak(scan_points, theoretical, tolerance_ppm)
        plus1 = _best_peak(scan_points, theoretical + spacing, tolerance_ppm)
        plus2 = _best_peak(scan_points, theoretical + 2.0 * spacing, tolerance_ppm)
        scan_envelopes.append((float(mono["intensity"]) if mono else 0.0, scan_id, mono, plus1, plus2))
    if not scan_envelopes:
        base["Envelope_Assessed"] = True
        base["Envelope_Status"] = "ENVELOPE_TOO_WEAK"
        return base
    _intensity, scan_id, mono, plus1, plus2 = max(scan_envelopes, key=lambda item: item[0])
    base["Envelope_Assessed"] = True
    base["Envelope_Spectrum_ID"] = scan_id
    if mono is None or float(mono["intensity"]) <= 0:
        base["Envelope_Status"] = "ENVELOPE_TOO_WEAK"
        return base

    mono_intensity = float(mono["intensity"])
    plus1_intensity = float(plus1["intensity"]) if plus1 else 0.0
    plus2_intensity = float(plus2["intensity"]) if plus2 else 0.0
    observed1 = plus1_intensity / mono_intensity
    observed2 = plus2_intensity / mono_intensity
    similarity = _cosine_similarity([1.0, expected["M+1"], expected["M+2"]], [1.0, observed1, observed2])
    incompatible = (
        similarity < 0.80
        or (expected["M+1"] >= 0.01 and plus1 is None)
        or (expected["M+2"] >= 0.02 and plus2 is None)
    )
    status = "ENVELOPE_INCOMPATIBLE" if incompatible else "ENVELOPE_COMPATIBLE"
    sulfur_count = composition.get("S", 0)
    if sulfur_count:
        sulfur_status = "SULFUR_ENVELOPE_COMPATIBLE" if plus2 is not None and not incompatible else "SULFUR_ENVELOPE_INCOMPATIBLE"
    else:
        sulfur_status = "NOT_ASSESSED"
    selected = [point for point in (mono, plus1, plus2) if point is not None]
    base.update({
        "Mplus2_Observed_Contribution": observed2,
        "Sulfur_Envelope_Compatible": sulfur_status,
        "Envelope_Status": status,
        "Observed_Mz_M": float(mono["mz"]),
        "Observed_Intensity_M": mono_intensity,
        "Observed_Mz_Mplus1": float(plus1["mz"]) if plus1 else "",
        "Observed_Intensity_Mplus1": plus1_intensity,
        "Observed_Mz_Mplus2": float(plus2["mz"]) if plus2 else "",
        "Observed_Intensity_Mplus2": plus2_intensity,
        "Relative_Intensity_Error_Mplus1": abs(observed1 - expected["M+1"]),
        "Relative_Intensity_Error_Mplus2": abs(observed2 - expected["M+2"]),
        "Envelope_Correlation": similarity,
        "Envelope_Cosine_Similarity": similarity,
        "Envelope_Point_Count": len(selected),
        "_isotope_point_indices": {int(point["index"]) for point in selected if "index" in point},
    })
    return base


def _quality_measurements(
    feature: dict[str, Any],
    candidate: dict[str, Any],
    raw_points: list[dict[str, Any]],
    spectrum_peaks: list[dict[str, Any]],
    settings: dict[str, Any],
) -> dict[str, Any]:
    inside, outside = _background_points(feature, candidate, raw_points, settings)
    inside_intensities = [float(point["intensity"]) for point in inside]
    outside_intensities = [float(point["intensity"]) for point in outside]
    median_inside = _median(inside_intensities)
    baseline = _median(outside_intensities) if outside_intensities else 0.0
    noise = _mad(outside_intensities, baseline) if outside_intensities else _mad(inside_intensities, median_inside)
    local_rows = [row for row in spectrum_peaks if row["Feature_ID"] == feature["Feature_ID"]]
    local_rows = sorted(local_rows, key=lambda row: (float(row["RT"]), str(row["Spectrum_ID"])))
    scans = {str(row["Spectrum_ID"]) for row in local_rows}
    profile_count = sum(int(row["Local_Profile_Point_Count"]) for row in local_rows)
    rts = sorted({float(row["RT"]) for row in local_rows})
    maximum_gap = max((right - left for left, right in zip(rts, rts[1:])), default=0.0)
    feature_start = float(feature.get("RT_Start") or (rts[0] if rts else 0.0))
    feature_end = float(feature.get("RT_End") or (rts[-1] if rts else feature_start))
    feature_span = max(0.0, feature_end - feature_start)
    run_rts = [float(point["rt"]) for point in raw_points if point.get("rt") is not None]
    run_span = max(run_rts, default=0.0) - min(run_rts, default=0.0)
    all_scans = {str(point.get("scan_id")) for point in raw_points if point.get("scan_id") not in {None, ""}}
    apex_intensity = float(feature.get("Apex_Intensity") or max(inside_intensities, default=0.0))
    apex_rt = float(feature.get("RT_Apex") or 0.0)
    left_span = max(0.0, apex_rt - feature_start)
    right_span = max(0.0, feature_end - apex_rt)
    symmetry = min(left_span, right_span) / max(left_span, right_span) if left_span and right_span else 0.0
    tailing = max(left_span, right_span) / min(left_span, right_span) if left_span and right_span else (float("inf") if inside else 0.0)
    background_scans = {str(point.get("scan_id")) for point in outside if point.get("scan_id") not in {None, ""}}
    background_fraction = len(background_scans) / max(len(all_scans), 1)
    coverage = feature_span / run_span if run_span > 0 else 0.0
    local_maxima = _detect_local_maxima(local_rows)
    centroid = (
        sum(float(point["mz"]) * float(point["intensity"]) for point in inside) / sum(inside_intensities)
        if sum(inside_intensities) else float(feature.get("mz_Centroid") or feature.get("Centroid_mz") or 0.0)
    )
    return {
        "Feature_ID": feature.get("Feature_ID"),
        "Physical_Feature_ID": feature.get("Physical_Feature_ID"),
        "Chemical_State_ID": feature.get("Chemical_State_ID"),
        "Chemical_Family": feature.get("Chemical_Family"),
        "Charge": feature.get("Charge"),
        "RT_Start": feature_start,
        "RT_End": feature_end,
        "RT_Apex": apex_rt,
        "RT_Span": feature_span,
        "Spectrum_Count": len(scans),
        "Profile_Point_Count": profile_count,
        "Apex_mz": feature.get("Apex_mz"),
        "Centroid_mz": centroid,
        "Mass_Error_ppm_Apex": feature.get("Mass_Error_ppm_at_Apex"),
        "Mass_Error_ppm_Centroid": feature.get("Mass_Error_ppm_at_Centroid"),
        "mz_SD": feature.get("mz_SD"),
        "Apex_Intensity": apex_intensity,
        "Integrated_Intensity": sum(inside_intensities),
        "Median_Intensity": median_inside,
        "Baseline_Intensity": baseline,
        "Local_Noise": noise,
        "Signal_to_Noise": apex_intensity / max(noise, 1.0),
        "Apex_to_Baseline_Ratio": apex_intensity / max(baseline, 1.0),
        "Left_Point_Count": sum(float(point.get("rt") or 0.0) < apex_rt for point in inside),
        "Right_Point_Count": sum(float(point.get("rt") or 0.0) > apex_rt for point in inside),
        "Left_RT_Span": left_span,
        "Right_RT_Span": right_span,
        "Peak_Symmetry": symmetry,
        "Peak_Tailing_Factor": tailing,
        "Scan_Continuity": max(0.0, 1.0 - maximum_gap / feature_span) if feature_span else 0.0,
        "Maximum_Scan_Gap": maximum_gap,
        "RT_Continuity": max(0.0, 1.0 - maximum_gap / feature_span) if feature_span else 0.0,
        "Run_RT_Coverage_Fraction": coverage,
        "Matched_Spectrum_Fraction": len(scans) / max(len(all_scans), 1),
        "Apex_Local_Contrast": apex_intensity / max(baseline, 1.0),
        "Outside_Feature_Intensity_Ratio": baseline / max(median_inside, 1.0),
        "Local_Maximum_Count": local_maxima,
        "Background_Persistence": background_fraction,
        "Feature_Width_Status": "WIDE" if coverage >= settings["background_rt_coverage_fraction"] else "LOCALIZED",
        "Feature_Shape_Status": "GOOD" if symmetry >= settings["feature_symmetry_min"] and tailing <= settings["feature_tailing_max"] else "POOR",
    }


def _classify_quality(
    measurements: dict[str, Any],
    feature: dict[str, Any],
    isotope: dict[str, Any],
    competition: dict[str, Any],
    settings: dict[str, Any],
) -> dict[str, Any]:
    candidate_ids = sorted(competition["candidate_ids"])
    contrast = float(measurements["Apex_Local_Contrast"])
    persistent = (
        float(measurements["Background_Persistence"]) >= settings["background_spectrum_fraction"]
        or float(measurements["Run_RT_Coverage_Fraction"]) >= settings["background_rt_coverage_fraction"]
    )
    background_rejected = (
        persistent and contrast < settings["background_contrast_min"]
    ) or (
        float(measurements["Run_RT_Coverage_Fraction"]) >= settings["background_rt_coverage_fraction"]
        and int(measurements["Local_Maximum_Count"]) > settings["local_maximum_count_threshold"]
    )
    if int(measurements["Spectrum_Count"]) < settings["min_spectrum_count"] or int(measurements["Profile_Point_Count"]) < settings["min_profile_point_count"]:
        status, reason = "SINGLE_POINT_REJECTED", "insufficient_profile_support"
    elif background_rejected:
        status, reason = "BACKGROUND_TRACE_REJECTED", "persistent_background_trace"
    elif feature.get("Chemical_Family") == "MODEL_NOT_DEFINED":
        status, reason = "MODEL_NOT_DEFINED", "model_not_defined"
    elif candidate_ids:
        status, reason = "COMPETITION_UNRESOLVED", "shared_physical_or_local_profile_feature"
    elif isotope["Envelope_Status"] == "ENVELOPE_INCOMPATIBLE":
        status, reason = "ISOTOPE_INCOMPATIBLE", "incompatible_isotope_envelope"
    else:
        status, reason = "QUALIFIED_CHROMATOGRAPHIC_FEATURE", ""
    measurements.update({
        "Background_Status": "PERSISTENT_BACKGROUND_TRACE" if background_rejected else "LOCALIZED_CHROMATOGRAPHIC_FEATURE",
        "Feature_Quality_Status": status,
        "Feature_Rejection_Reasons": reason or "NOT_ASSESSED",
        "Raw_Match_Present": True,
        "Independent_Feature_Qualified": status.startswith("QUALIFIED"),
        "Feature_Eligible_For_Support": status.startswith("QUALIFIED"),
        "Competition_Count": len(candidate_ids),
        "Competing_Candidate_IDs": ";".join(candidate_ids),
        "Isomer_Isotope_Indistinguishable": bool(competition["isomer"]),
        **FORMAL_FALSE,
    })
    return measurements


def _quality_score_components(row: dict[str, Any], isotope: dict[str, Any]) -> dict[str, float]:
    mass = max(0.0, 1.0 - min(abs(float(row.get("Mass_Error_ppm_Centroid") or 0.0)) / 20.0, 1.0))
    isotope_component = 0.0 if isotope["Envelope_Status"] == "ENVELOPE_INCOMPATIBLE" else (
        1.0 if isotope["Envelope_Status"] == "ENVELOPE_COMPATIBLE" else 0.5
    )
    return {
        "Mass_Accuracy_Component": mass,
        "Spectrum_Continuity_Component": max(0.0, min(float(row["Scan_Continuity"]), 1.0)),
        "Chromatographic_Localization_Component": max(0.0, min(float(row["Apex_Local_Contrast"]) / 10.0, 1.0)),
        "Peak_Shape_Component": 1.0 if row["Feature_Shape_Status"] == "GOOD" else 0.5,
        "Background_Contrast_Component": 1.0 if row["Background_Status"] == "LOCALIZED_CHROMATOGRAPHIC_FEATURE" else 0.0,
        "Isotope_Envelope_Component": isotope_component,
    }


def _propagate_refined_quality(
    refined: list[dict[str, Any]],
    quality_rows: list[dict[str, Any]],
) -> None:
    status_by_id = {str(row["Feature_ID"]): str(row["Feature_Quality_Status"]) for row in quality_rows}
    for row in refined:
        legacy_ids = [value for value in str(row.get("Legacy_Feature_ID") or "").split(";") if value]
        statuses = {status_by_id[value] for value in legacy_ids if value in status_by_id}
        if not statuses:
            continue
        row["Feature_Quality_Status"] = next(iter(statuses)) if len(statuses) == 1 else ";".join(sorted(statuses))


def _mark_shared_isotope_peaks(
    isotope_rows: list[dict[str, Any]],
    competition: dict[str, dict[str, Any]],
) -> None:
    owners: dict[int, set[str]] = defaultdict(set)
    for row in isotope_rows:
        for index in row["_isotope_point_indices"]:
            owners[index].add(str(row["Candidate_ID"]))
    for row in isotope_rows:
        shared = any(len(owners[index]) > 1 for index in row["_isotope_point_indices"])
        row["Isotope_Peak_Shared_With_Other_Candidate"] = shared
        row["Isomer_Isotope_Indistinguishable"] = bool(competition[str(row["Feature_ID"])]["isomer"])


def build_p1_sap_feature_quality(
    candidates: list[dict[str, Any]],
    features: list[dict[str, Any]],
    raw_peaks: list[dict[str, Any]],
    config: Any,
    tolerance_ppm: float,
) -> dict[str, Any]:
    candidates_by_id = {str(candidate["Chemical_State_ID"]): candidate for candidate in candidates}
    settings = _quality_config(config)
    spectrum_peaks = build_spectrum_level_peaks(features, raw_peaks, tolerance_ppm)
    refined = build_refined_chromatographic_features(spectrum_peaks, candidates_by_id, features, config, tolerance_ppm)
    competition = _competition_context(features, spectrum_peaks, candidates_by_id, tolerance_ppm)

    isotope_internal: list[dict[str, Any]] = []
    for feature in features:
        candidate = candidates_by_id.get(str(feature["Chemical_State_ID"]), {})
        isotope_internal.append(
            assess_isotope_envelope(feature, raw_peaks, candidate, settings["isotope_match_tolerance_ppm"])
        )
    _mark_shared_isotope_peaks(isotope_internal, competition)
    isotope_by_feature = {str(row["Feature_ID"]): row for row in isotope_internal}

    quality_rows: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    scores: list[float] = []
    for feature in features:
        feature_id = str(feature["Feature_ID"])
        candidate = candidates_by_id.get(str(feature["Chemical_State_ID"]), {})
        isotope = isotope_by_feature[feature_id]
        measurements = _quality_measurements(feature, candidate, raw_peaks, spectrum_peaks, settings)
        quality = _classify_quality(measurements, feature, isotope, competition[feature_id], settings)
        components = _quality_score_components(quality, isotope)
        quality.update(components)
        quality["Feature_Quality_Score"] = sum(components.values()) / len(components)
        scores.append(float(quality["Feature_Quality_Score"]))
        output_row = {column: quality.get(column, "") for column in P1_SAP_FEATURE_QUALITY_COLUMNS}
        quality_rows.append(output_row)

        feature.update(output_row)
        feature["Envelope_Assessed"] = isotope["Envelope_Assessed"]
        feature["Isotope_Status"] = isotope["Envelope_Status"]
        feature["Chemical_State_Supported"] = bool(quality["Feature_Eligible_For_Support"])
        if not quality["Feature_Eligible_For_Support"]:
            feature["Final_Interpretation"] = "NOT_EVALUABLE"
            feature["Feature_Exclusion_Reason"] = quality["Feature_Rejection_Reasons"]

        status = str(quality["Feature_Quality_Status"])
        counts[status] += 1
        if status.startswith("QUALIFIED"):
            counts["Qualified_Feature_Count"] += 1
        if status == "SINGLE_POINT_REJECTED":
            counts["Rejected_Single_Point_Count"] += 1
        elif status == "BACKGROUND_TRACE_REJECTED":
            counts["Rejected_Background_Count"] += 1
        elif status == "PROFILE_ONLY_REJECTED":
            counts["Rejected_Profile_Only_Count"] += 1
        elif status == "ISOTOPE_INCOMPATIBLE":
            counts["Isotope_Incompatible_Count"] += 1
        elif status == "COMPETITION_UNRESOLVED":
            counts["Competition_Unresolved_Count"] += 1
        elif status == "MODEL_NOT_DEFINED":
            counts["Model_Not_Defined_Count"] += 1
        elif status == "RAW_MATCH_ONLY":
            counts["Raw_Match_Only_Count"] += 1
        if feature.get("Chemical_Family") in {"PHOSPHOROTHIOATE", "P1_RESISTANT_PT_OLIGOMER"} and status.startswith("QUALIFIED"):
            counts["Qualified_PT_Feature_Count"] += 1

    _propagate_refined_quality(refined, quality_rows)
    isotope_rows = [
        {column: row.get(column, "") for column in P1_SAP_ISOTOPE_AUDIT_COLUMNS}
        for row in isotope_internal
    ]
    summary = {
        "Feature_Quality_Model_Version": FEATURE_QUALITY_MODEL_VERSION,
        "Spectrum_Level_Local_Peak_Count": len(spectrum_peaks),
        "Refined_Feature_Count": len(refined),
        "Qualified_Feature_Count": counts["Qualified_Feature_Count"],
        "Rejected_Single_Point_Count": counts["Rejected_Single_Point_Count"],
        "Rejected_Background_Count": counts["Rejected_Background_Count"],
        "Rejected_Profile_Only_Count": counts["Rejected_Profile_Only_Count"],
        "Isotope_Incompatible_Count": counts["Isotope_Incompatible_Count"],
        "Competition_Unresolved_Count": counts["Competition_Unresolved_Count"],
        "Model_Not_Defined_Count": counts["Model_Not_Defined_Count"],
        "Raw_Match_Only_Count": counts["Raw_Match_Only_Count"],
        "Qualified_PT_Feature_Count": counts["Qualified_PT_Feature_Count"],
        "PT_Final_Interpretation": "NO_QUALIFIED_PT_FEATURE" if counts["Qualified_PT_Feature_Count"] == 0 else "PT_LIKE_STATE_AMBIGUOUS",
        "Feature_Quality_Score_Mean": sum(scores) / len(scores) if scores else 0.0,
        "Feature_Quality_Score_Median": _median(scores),
        "Feature_Quality_Score_Min": min(scores) if scores else 0.0,
        "Feature_Quality_Score_Max": max(scores) if scores else 0.0,
        **FORMAL_FALSE,
    }
    return {
        "spectrum_peaks": spectrum_peaks,
        "refined_features": refined,
        "quality_rows": quality_rows,
        "isotope_rows": isotope_rows,
        "summary_row": summary,
    }
