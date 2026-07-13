"""Shadow audit for unmatched modified theoretical MS2 ions.

The audit consumes existing matching results and raw spectrum metadata.  It never
changes match acceptance, localization, ranking, confidence, or candidate membership.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any

from rna_masshunter.ms1_mapping import ppm_error
from rna_masshunter.ms2_identity_evidence import physical_observed_peak_key

AUDIT_COLUMNS = [
    "Modification_ID", "Modification_Name", "Parent_Fragment_ID", "Spectrum_ID",
    "Candidate_tRNA_Position", "Candidate_Position_In_Parent", "Candidate_Base",
    "Theoretical_Ion_ID", "Ion_Series", "Ion_Number", "Ion_Charge", "Theoretical_mz",
    "Is_Matched", "Existing_Match_IDs", "Unmatched_Reason_Status", "Scan_mz_Min",
    "Scan_mz_Max", "Matching_Tolerance_ppm", "Matching_Tolerance_Da",
    "Audit_Search_Window_Da", "Nearest_Physical_Observed_Peak_Key", "Nearest_Observed_mz",
    "Nearest_Intensity", "Nearest_Error_Da", "Nearest_Error_ppm", "Nearby_Peak_Count",
    "Nearby_Physical_Peak_Keys", "Intensity_Threshold", "Threshold_Information_Available",
    "Below_Threshold", "Audit_Interpretation", "Audit_Warnings",
]

SUMMARY_COLUMNS = [
    "Rank", "Modification_ID", "Parent_Fragment_ID", "Candidate_tRNA_Position",
    "Total_Modified_Theoretical_Ion_Count", "Matched_Modified_Theoretical_Ion_Count",
    "Unmatched_Modified_Theoretical_Ion_Count", "Outside_Scan_Range_Count",
    "No_Peak_In_Window_Count", "Nearest_Peak_Outside_Tolerance_Count",
    "Below_Threshold_Count", "Filtered_Peak_Count", "Ambiguous_Nearby_Peak_Count",
    "Information_Unavailable_Count", "Best_Unmatched_Error_ppm", "Unmatched_Reason_Summary",
    "Recommended_Followup", "Audit_Warnings",
]

TOP_SHADOW_COLUMNS = [
    "Total_Modified_Theoretical_Ion_Count", "Matched_Modified_Theoretical_Ion_Count",
    "Unmatched_Modified_Theoretical_Ion_Count", "Primary_Unmatched_Reason",
    "Outside_Scan_Range_Count", "No_Peak_In_Window_Count",
    "Nearest_Peak_Outside_Tolerance_Count", "Below_Threshold_Count",
    "Information_Unavailable_Count", "Best_Unmatched_Error_ppm",
    "Unmatched_Ion_Audit_Warnings",
]

DIAGNOSTIC_COLUMNS = [
    "MS2_Unmatched_Ion_Audit_Enabled", "Apply_Unmatched_Audit_To_Final_Score",
    "Total_Modified_Theoretical_Ions", "Matched_Modified_Theoretical_Ions",
    "Unmatched_Modified_Theoretical_Ions", "Outside_Scan_Range_Count",
    "No_Peak_In_Window_Count", "Nearest_Peak_Outside_Tolerance_Count",
    "Below_Threshold_Count", "Filtered_Peak_Count", "Ambiguous_Nearby_Peak_Count",
    "Spectrum_Not_Available_Count", "Scan_Range_Not_Available_Count",
    "Threshold_Information_Unavailable_Count", "Insufficient_Information_Count",
    "Audit_Search_Window_Rule", "Matching_Tolerance_Source", "Intensity_Threshold_Source",
]

INFORMATION_STATUSES = {
    "spectrum_not_available", "scan_range_not_available",
    "threshold_information_unavailable", "insufficient_information",
}

STATUS_TO_COUNT = {
    "outside_scan_mz_range": "Outside_Scan_Range_Count",
    "no_observed_peak_in_search_window": "No_Peak_In_Window_Count",
    "nearest_peak_outside_tolerance": "Nearest_Peak_Outside_Tolerance_Count",
    "peak_below_intensity_threshold": "Below_Threshold_Count",
    "peak_present_but_filtered": "Filtered_Peak_Count",
    "ambiguous_multiple_nearby_peaks": "Ambiguous_Nearby_Peak_Count",
}


def _bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().lower() in {"1", "true", "yes", "on"}


def _float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _position(value: Any) -> int | None:
    number = _float(value)
    return int(number) if number is not None and number > 0 else None


def _candidate_key(row: dict[str, Any], position_field: str) -> tuple[str, str, int | None]:
    return (
        str(row.get("Modification_ID") or ""),
        str(row.get("Parent_Fragment_ID") or ""),
        _position(row.get(position_field)),
    )


def _audit_ion_key(row: dict[str, Any]) -> tuple[Any, ...]:
    """Semantic key prevents duplicate auditing even when duplicate Ion_ID rows exist."""
    return (
        str(row.get("Spectrum_ID") or ""), str(row.get("Parent_Fragment_ID") or ""),
        str(row.get("Modification_ID") or ""),
        _position(row.get("Candidate_Modification_Position_In_Parent")),
        str(row.get("Ion_Type") or ""), _position(row.get("Ion_Start")),
        _position(row.get("Ion_End")), _position(row.get("Charge")),
        _float(row.get("Theoretical_mz")),
    )


def _match_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        str(row.get("Spectrum_ID") or ""), str(row.get("Parent_Fragment_ID") or ""),
        str(row.get("Modification_ID") or ""),
        _position(row.get("Candidate_Modification_Position_In_Parent")),
        str(row.get("Ion_ID") or ""),
    )


def _existing_match_id(row: dict[str, Any]) -> str:
    values = (
        row.get("Spectrum_ID"), row.get("Scan_Index"), row.get("Observed_mz"),
        row.get("Ion_ID"), row.get("Theoretical_mz"),
    )
    return ":".join(str(value if value not in (None, "") else "NA") for value in values)


def _trna_position(ion: dict[str, Any]) -> int | None:
    parent_start = _position(ion.get("Parent_Start"))
    position = _position(ion.get("Candidate_Modification_Position_In_Parent"))
    return parent_start + position - 1 if parent_start is not None and position is not None else None


def _ion_number(ion: dict[str, Any]) -> int | None:
    ion_type = str(ion.get("Ion_Type") or "")
    if ion_type == "c":
        return _position(ion.get("Ion_End"))
    return _position(ion.get("Ion_Length"))


def _peak_row(spectrum: Any, mz: float, intensity: float) -> dict[str, Any]:
    return {
        "Spectrum_ID": getattr(spectrum, "spectrum_id", ""),
        "Observed_mz": mz, "Observed_Intensity": intensity,
        "RT": getattr(spectrum, "rt", None),
    }


def _physical_key(spectrum: Any, mz: float, intensity: float) -> str:
    return physical_observed_peak_key(_peak_row(spectrum, mz, intensity))


def _is_selected_peak(spectrum: Any, peak: tuple[float, float]) -> bool:
    mz, intensity = peak
    return any(float(selected_mz) == mz and float(selected_intensity) == intensity for selected_mz, selected_intensity in (getattr(spectrum, "peaks", None) or []))


def _base_row(ion: dict[str, Any], tolerance_ppm: float) -> dict[str, Any]:
    theoretical_mz = _float(ion.get("Theoretical_mz"))
    tolerance_da = abs(theoretical_mz) * tolerance_ppm / 1_000_000 if theoretical_mz is not None else None
    audit_window = max((tolerance_da or 0.0) * 5.0, 0.05) if theoretical_mz is not None else None
    return {
        "Modification_ID": ion.get("Modification_ID"),
        "Modification_Name": ion.get("Modification_Name"),
        "Parent_Fragment_ID": ion.get("Parent_Fragment_ID"),
        "Spectrum_ID": ion.get("Spectrum_ID"),
        "Candidate_tRNA_Position": _trna_position(ion) or "",
        "Candidate_Position_In_Parent": ion.get("Candidate_Modification_Position_In_Parent"),
        "Candidate_Base": ion.get("Candidate_Modification_Base"),
        "Theoretical_Ion_ID": ion.get("Ion_ID"),
        "Ion_Series": ion.get("Ion_Type"), "Ion_Number": _ion_number(ion) or "",
        "Ion_Charge": ion.get("Charge"), "Theoretical_mz": ion.get("Theoretical_mz"),
        "Is_Matched": False, "Existing_Match_IDs": "", "Unmatched_Reason_Status": "insufficient_information",
        "Scan_mz_Min": "", "Scan_mz_Max": "", "Matching_Tolerance_ppm": tolerance_ppm,
        "Matching_Tolerance_Da": tolerance_da if tolerance_da is not None else "",
        "Audit_Search_Window_Da": audit_window if audit_window is not None else "",
        "Nearest_Physical_Observed_Peak_Key": "", "Nearest_Observed_mz": "",
        "Nearest_Intensity": "", "Nearest_Error_Da": "", "Nearest_Error_ppm": "",
        "Nearby_Peak_Count": 0, "Nearby_Physical_Peak_Keys": "", "Intensity_Threshold": "",
        "Threshold_Information_Available": False, "Below_Threshold": "",
        "Audit_Interpretation": "Required information was unavailable.", "Audit_Warnings": "",
    }


def build_unmatched_ion_audit(
    spectra: list[Any] | None,
    modified_ions: list[dict[str, Any]] | None,
    modified_matches: list[dict[str, Any]] | None,
    config: Any,
    enabled: bool = True,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Return per-ion audit rows and one diagnostics row without rematching."""
    ms2 = getattr(config, "ms2_annotation", {}) or {}
    tolerance_ppm = float(ms2.get("mz_tolerance_ppm", 20) or 20)
    spectra_lookup = {str(spectrum.spectrum_id): spectrum for spectrum in (spectra or [])}
    matches_by_key: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for match in modified_matches or []:
        if _bool(match.get("Ion_Contains_Modification")):
            matches_by_key[_match_key(match)].append(match)

    unique_ions: dict[tuple[Any, ...], dict[str, Any]] = {}
    for ion in modified_ions or []:
        if not _bool(ion.get("Ion_Contains_Modification")):
            continue
        unique_ions.setdefault(_audit_ion_key(ion), ion)

    rows: list[dict[str, Any]] = []
    if enabled:
        for ion in unique_ions.values():
            row = _base_row(ion, tolerance_ppm)
            exact_matches = matches_by_key.get((
                str(ion.get("Spectrum_ID") or ""), str(ion.get("Parent_Fragment_ID") or ""),
                str(ion.get("Modification_ID") or ""),
                _position(ion.get("Candidate_Modification_Position_In_Parent")), str(ion.get("Ion_ID") or ""),
            ), [])
            spectrum = spectra_lookup.get(str(ion.get("Spectrum_ID") or ""))
            theoretical_mz = _float(ion.get("Theoretical_mz"))
            tolerance_da = _float(row["Matching_Tolerance_Da"])
            audit_window = _float(row["Audit_Search_Window_Da"])

            if exact_matches:
                best = min(exact_matches, key=lambda item: abs(_float(item.get("Mass_Error_ppm")) or float("inf")))
                observed_mz = _float(best.get("Observed_mz"))
                intensity = _float(best.get("Observed_Intensity"))
                row.update({
                    "Is_Matched": True, "Existing_Match_IDs": ";".join(sorted(_existing_match_id(item) for item in exact_matches)),
                    "Unmatched_Reason_Status": "matched", "Nearest_Observed_mz": observed_mz if observed_mz is not None else "",
                    "Nearest_Intensity": intensity if intensity is not None else "",
                    "Nearest_Error_Da": best.get("Mass_Error_Da", ""), "Nearest_Error_ppm": best.get("Mass_Error_ppm", ""),
                    "Nearby_Peak_Count": 1, "Audit_Interpretation": "Existing modified-ion match retained without reevaluation.",
                })
                if spectrum is not None and observed_mz is not None and intensity is not None:
                    row["Nearest_Physical_Observed_Peak_Key"] = _physical_key(spectrum, observed_mz, intensity)
                    row["Nearby_Physical_Peak_Keys"] = row["Nearest_Physical_Observed_Peak_Key"]
            elif spectrum is None:
                row.update({"Unmatched_Reason_Status": "spectrum_not_available", "Audit_Interpretation": "Target spectrum was not available for audit."})
            elif theoretical_mz is None or tolerance_da is None or audit_window is None:
                row.update({"Unmatched_Reason_Status": "insufficient_information", "Audit_Interpretation": "Theoretical m/z or tolerance information was unavailable."})
            else:
                scan_min = _float(getattr(spectrum, "scan_mz_min", None))
                scan_max = _float(getattr(spectrum, "scan_mz_max", None))
                threshold = _float(getattr(spectrum, "effective_intensity_threshold", None))
                threshold_available = _bool(getattr(spectrum, "threshold_information_available", False))
                row.update({
                    "Scan_mz_Min": scan_min if scan_min is not None else "", "Scan_mz_Max": scan_max if scan_max is not None else "",
                    "Intensity_Threshold": threshold if threshold is not None else "",
                    "Threshold_Information_Available": threshold_available,
                })
                if scan_min is not None and scan_max is not None and not (scan_min <= theoretical_mz <= scan_max):
                    row.update({"Unmatched_Reason_Status": "outside_scan_mz_range", "Audit_Interpretation": "Theoretical m/z is outside the recorded scan window."})
                else:
                    raw_peaks = getattr(spectrum, "raw_peaks", None)
                    if raw_peaks is None:
                        status = "scan_range_not_available" if scan_min is None or scan_max is None else "insufficient_information"
                        row.update({"Unmatched_Reason_Status": status, "Audit_Interpretation": "Raw spectrum peaks were unavailable for audit."})
                    else:
                        nearby = sorted(
                            [(float(mz), float(intensity)) for mz, intensity in raw_peaks if abs(float(mz) - theoretical_mz) <= audit_window],
                            key=lambda peak: (abs(peak[0] - theoretical_mz), -peak[1]),
                        )
                        row["Nearby_Peak_Count"] = len(nearby)
                        row["Nearby_Physical_Peak_Keys"] = ";".join(_physical_key(spectrum, mz, intensity) for mz, intensity in nearby)
                        if nearby:
                            nearest_mz, nearest_intensity = nearby[0]
                            error_da = nearest_mz - theoretical_mz
                            error_ppm = ppm_error(nearest_mz, theoretical_mz)
                            row.update({
                                "Nearest_Physical_Observed_Peak_Key": _physical_key(spectrum, nearest_mz, nearest_intensity),
                                "Nearest_Observed_mz": nearest_mz, "Nearest_Intensity": nearest_intensity,
                                "Nearest_Error_Da": error_da, "Nearest_Error_ppm": error_ppm,
                                "Below_Threshold": nearest_intensity < threshold if threshold_available and threshold is not None else "",
                            })
                        if not nearby:
                            if scan_min is None or scan_max is None:
                                row.update({"Unmatched_Reason_Status": "scan_range_not_available", "Audit_Interpretation": "No nearby peak was found and scan range was unavailable."})
                            else:
                                row.update({"Unmatched_Reason_Status": "no_observed_peak_in_search_window", "Audit_Interpretation": "No raw observed peak was present in the audit search window."})
                        elif len(nearby) > 1:
                            row.update({"Unmatched_Reason_Status": "ambiguous_multiple_nearby_peaks", "Audit_Interpretation": "Multiple raw peaks in the audit window prevent a unique nearby assignment."})
                        elif abs(nearby[0][0] - theoretical_mz) > tolerance_da:
                            row.update({"Unmatched_Reason_Status": "nearest_peak_outside_tolerance", "Audit_Interpretation": "Nearest raw peak is inside the audit window but outside formal tolerance."})
                        elif not threshold_available or threshold is None:
                            row.update({"Unmatched_Reason_Status": "threshold_information_unavailable", "Audit_Interpretation": "A raw peak is within tolerance but threshold information is unavailable."})
                        elif nearby[0][1] < threshold:
                            row.update({"Unmatched_Reason_Status": "peak_below_intensity_threshold", "Below_Threshold": True, "Audit_Interpretation": "A raw peak is within tolerance but below the effective formal intensity threshold."})
                        elif not _is_selected_peak(spectrum, nearby[0]):
                            row.update({"Unmatched_Reason_Status": "peak_present_but_filtered", "Below_Threshold": False, "Audit_Interpretation": "A threshold-passing raw peak was removed by a traced post-threshold peak filter."})
                        else:
                            row.update({
                                "Unmatched_Reason_Status": "insufficient_information",
                                "Below_Threshold": False,
                                "Audit_Interpretation": "A selected peak is within tolerance but no exact existing candidate-ion match exists; existing matching is not reinterpreted.",
                                "Audit_Warnings": "within-tolerance peak was assigned differently by existing matching",
                            })
            if spectrum is not None:
                if row["Scan_mz_Min"] == "": row["Scan_mz_Min"] = getattr(spectrum, "scan_mz_min", "") or ""
                if row["Scan_mz_Max"] == "": row["Scan_mz_Max"] = getattr(spectrum, "scan_mz_max", "") or ""
                if row["Intensity_Threshold"] == "": row["Intensity_Threshold"] = getattr(spectrum, "effective_intensity_threshold", "")
                row["Threshold_Information_Available"] = _bool(getattr(spectrum, "threshold_information_available", False))
            rows.append(row)

    counts = Counter(str(row.get("Unmatched_Reason_Status") or "insufficient_information") for row in rows)
    diagnostics = [{
        "MS2_Unmatched_Ion_Audit_Enabled": bool(enabled),
        "Apply_Unmatched_Audit_To_Final_Score": False,
        "Total_Modified_Theoretical_Ions": len(rows),
        "Matched_Modified_Theoretical_Ions": counts["matched"],
        "Unmatched_Modified_Theoretical_Ions": len(rows) - counts["matched"],
        "Outside_Scan_Range_Count": counts["outside_scan_mz_range"],
        "No_Peak_In_Window_Count": counts["no_observed_peak_in_search_window"],
        "Nearest_Peak_Outside_Tolerance_Count": counts["nearest_peak_outside_tolerance"],
        "Below_Threshold_Count": counts["peak_below_intensity_threshold"],
        "Filtered_Peak_Count": counts["peak_present_but_filtered"],
        "Ambiguous_Nearby_Peak_Count": counts["ambiguous_multiple_nearby_peaks"],
        "Spectrum_Not_Available_Count": counts["spectrum_not_available"],
        "Scan_Range_Not_Available_Count": counts["scan_range_not_available"],
        "Threshold_Information_Unavailable_Count": counts["threshold_information_unavailable"],
        "Insufficient_Information_Count": counts["insufficient_information"],
        "Audit_Search_Window_Rule": "max(5 * per-ion matching tolerance Da, 0.05 Da)",
        "Matching_Tolerance_Source": "config.ms2_annotation.mz_tolerance_ppm",
        "Intensity_Threshold_Source": "max(config.ms2_annotation.min_peak_intensity, raw base peak * config.ms2_annotation.min_relative_intensity_percent / 100)",
    }]
    return rows, diagnostics


def _followup(total: int, matched: int, counts: Counter[str]) -> str:
    unmatched = total - matched
    if total == 0:
        return "insufficient_information"
    if unmatched == 0 and matched:
        return "no_action_existing_match_present"
    ordered = [
        ("outside_scan_mz_range", "acquire_wider_mz_range"),
        ("peak_below_intensity_threshold", "lower_intensity_threshold_for_audit"),
        ("nearest_peak_outside_tolerance", "relax_tolerance_for_diagnostic_review_only"),
        ("ambiguous_multiple_nearby_peaks", "inspect_raw_spectrum_manually"),
        ("no_observed_peak_in_search_window", "improve_fragmentation_coverage"),
    ]
    status, action = max(ordered, key=lambda item: counts[item[0]])
    if counts[status] > 0:
        return action
    if any(counts[status] for status in INFORMATION_STATUSES):
        return "insufficient_information"
    return "acquire_additional_ms2_spectrum"


def build_unmatched_ion_summary(
    ranking_rows: list[dict[str, Any]] | None,
    audit_rows: list[dict[str, Any]] | None,
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, int | None], list[dict[str, Any]]] = defaultdict(list)
    for row in audit_rows or []:
        grouped[_candidate_key(row, "Candidate_Position_In_Parent")].append(row)
    summaries: list[dict[str, Any]] = []
    for ranking in ranking_rows or []:
        key = _candidate_key(ranking, "Candidate_Position_In_Parent")
        group = grouped.get(key, [])
        counts = Counter(str(row.get("Unmatched_Reason_Status") or "insufficient_information") for row in group)
        total = len(group)
        matched = counts["matched"]
        unmatched_errors = [abs(value) for row in group if row.get("Unmatched_Reason_Status") != "matched" and (value := _float(row.get("Nearest_Error_ppm"))) is not None]
        reason_summary = ";".join(f"{status}={count}" for status, count in sorted(counts.items()) if status != "matched" and count)
        warnings = []
        if total == 0: warnings.append("no modified theoretical ions for candidate")
        if any(counts[status] for status in INFORMATION_STATUSES): warnings.append("some unmatched reasons have unavailable information")
        summary = {
            "Rank": ranking.get("Rank"), "Modification_ID": ranking.get("Modification_ID"),
            "Parent_Fragment_ID": ranking.get("Parent_Fragment_ID"),
            "Candidate_tRNA_Position": ranking.get("Candidate_tRNA_Position"),
            "Total_Modified_Theoretical_Ion_Count": total,
            "Matched_Modified_Theoretical_Ion_Count": matched,
            "Unmatched_Modified_Theoretical_Ion_Count": total - matched,
            "Outside_Scan_Range_Count": counts["outside_scan_mz_range"],
            "No_Peak_In_Window_Count": counts["no_observed_peak_in_search_window"],
            "Nearest_Peak_Outside_Tolerance_Count": counts["nearest_peak_outside_tolerance"],
            "Below_Threshold_Count": counts["peak_below_intensity_threshold"],
            "Filtered_Peak_Count": counts["peak_present_but_filtered"],
            "Ambiguous_Nearby_Peak_Count": counts["ambiguous_multiple_nearby_peaks"],
            "Information_Unavailable_Count": sum(counts[status] for status in INFORMATION_STATUSES),
            "Best_Unmatched_Error_ppm": min(unmatched_errors) if unmatched_errors else "",
            "Unmatched_Reason_Summary": reason_summary,
            "Recommended_Followup": _followup(total, matched, counts),
            "Audit_Warnings": "; ".join(warnings),
        }
        summaries.append(summary)
    return summaries


def primary_unmatched_reason(summary: dict[str, Any]) -> str:
    pairs = [
        ("outside_scan_mz_range", "Outside_Scan_Range_Count"),
        ("no_observed_peak_in_search_window", "No_Peak_In_Window_Count"),
        ("nearest_peak_outside_tolerance", "Nearest_Peak_Outside_Tolerance_Count"),
        ("peak_below_intensity_threshold", "Below_Threshold_Count"),
        ("peak_present_but_filtered", "Filtered_Peak_Count"),
        ("ambiguous_multiple_nearby_peaks", "Ambiguous_Nearby_Peak_Count"),
        ("information_unavailable", "Information_Unavailable_Count"),
    ]
    best = max(pairs, key=lambda item: int(_float(summary.get(item[1])) or 0))
    if int(_float(summary.get(best[1])) or 0) > 0:
        return best[0]
    if int(_float(summary.get("Matched_Modified_Theoretical_Ion_Count")) or 0) > 0:
        return "matched"
    return "insufficient_information"
