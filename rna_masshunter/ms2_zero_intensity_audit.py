"""Shadow audit tracing zero-intensity MS2 peaks without changing analysis.

The audit records the pyteomics-decoded mzML arrays, the values after the
project's float conversion, and the peaks retained as annotation input.  It
then joins existing formal and shadow result rows.  It never filters, rematches,
or mutates an existing result.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from math import isfinite, isnan
from statistics import median
from typing import Any

from rna_masshunter.ms2_identity_evidence import physical_observed_peak_key


SPECTRA_COLUMNS = [
    "Spectrum_ID", "Scan_Index", "Precursor_MZ", "Precursor_Charge",
    "Total_Peak_Count", "Positive_Intensity_Count", "Zero_Intensity_Count",
    "Negative_Intensity_Count", "NaN_Intensity_Count", "Missing_Intensity_Count",
    "Non_Numeric_Intensity_Count", "Finite_Intensity_Count",
    "Positive_Intensity_Fraction", "Zero_Intensity_Fraction", "Min_Intensity",
    "Max_Intensity", "Median_Intensity", "Median_Positive_Intensity",
    "Total_Ion_Intensity", "MZ_Array_Length", "Intensity_Array_Length",
    "Array_Length_Mismatch", "Example_Zero_Intensity_MZ", "Spectrum_Mode",
    "Provenance", "Origin_Category", "Origin_Classification_Basis",
]

DETAIL_COLUMNS = [
    "Spectrum_ID", "Scan_Index", "Peak_Index", "Original_MZ", "Original_Intensity",
    "Parsed_Intensity", "Annotation_Input_Intensity", "Intensity_State",
    "Raw_Spectrum_Derived", "Synthetic_Or_Placeholder", "Provenance",
    "Transformation_History", "Used_For_Theoretical_Match", "Selected_As_Best_Match",
    "Used_For_Identity", "Used_For_Localization", "Used_In_Ambiguity_Cluster",
    "Ambiguity_Cluster_ID", "Candidate_Key", "Modification", "Position",
    "Ion_Series", "Ion_Number", "Error_PPM", "Within_Formal_Tolerance",
]

SUMMARY_COLUMNS = [
    "Total_Spectra", "Total_Raw_Peaks", "Total_Positive_Peaks",
    "Total_Zero_Intensity_Peaks", "Total_Negative_Peaks", "Total_NaN_Peaks",
    "Total_Missing_Peaks", "Total_Non_Numeric_Peaks", "Zero_Intensity_Fraction",
    "Spectra_With_Zero_Intensity", "All_Zero_Spectra", "Total_Matched_Peaks",
    "Matched_Zero_Intensity_Peaks", "Zero_Intensity_Best_Matches",
    "Zero_Intensity_Identity_Assignments", "Zero_Intensity_Localization_Uses",
    "Ambiguity_Clusters_With_Zero_Intensity", "All_Zero_Ambiguity_Clusters",
    "Affected_Candidate_Count", "Affected_Top50_Count", "Affected_cnm5U_Count",
    "Likely_Origin_Category", "Origin_Classification_Confidence",
    "Origin_Classification_Basis", "Recommended_Next_Action", "Audit_Mode",
    "Applied_To_Final_Score", "Shadow_Nonzero_Peak_Match_Count",
    "Shadow_Nonzero_Best_Match_Count", "Shadow_Nonzero_Ambiguity_Cluster_Count",
    "Shadow_Nonzero_Candidate_Specific_Count", "Shadow_Nonzero_Position_Discrimination",
    "Shadow_Nonzero_Conclusion", "Original_Detail_Row_Count",
    "Written_Detail_Row_Count", "Detail_Truncated", "Detail_Truncation_Reason",
    "Formal_Best_Match_Definition", "Ambiguity_Diagnostic_Best_Peak_Definition",
]

DIAGNOSTIC_COLUMNS = [
    "Zero_Intensity_Audit_Available", "Raw_Peak_Count",
    "Raw_Zero_Intensity_Peak_Count", "Raw_Zero_Intensity_Fraction",
    "Raw_Positive_Intensity_Peak_Count", "Raw_Median_Positive_Intensity",
    "Annotation_Input_Peak_Count", "Annotation_Input_Zero_Intensity_Count",
    "Zero_Intensity_Matched_Peak_Count", "Zero_Intensity_Best_Match_Count",
    "Zero_Intensity_Identity_Assignment_Count", "Zero_Intensity_Localization_Count",
    "Zero_Intensity_Ambiguity_Cluster_Count", "All_Zero_Ambiguity_Cluster_Count",
    "Zero_Intensity_Origin_Category", "Zero_Intensity_Audit_Severity",
    "Zero_Intensity_Audit_Recommendation",
    "Zero_Intensity_Audit_Applied_To_Final_Score",
]

TOP_SHADOW_COLUMNS = [
    "Zero_Intensity_Affected", "Zero_Intensity_Cluster_Count",
    "All_Zero_Cluster_Count", "Zero_Intensity_Best_Match_Count",
    "Positive_Intensity_Support_Count", "Zero_Intensity_Dependence",
    "Zero_Intensity_Audit_Severity", "Zero_Intensity_Audit_Recommendation",
    "Zero_Intensity_Audit_Applied_To_Final_Score",
]

CANDIDATE_COLUMNS = [
    "Modification_ID", "Parent_Fragment_ID", "Candidate_tRNA_Position",
    *TOP_SHADOW_COLUMNS, "Shadow_Nonzero_Peak_Match_Count",
    "Shadow_Nonzero_Best_Match_Count", "Shadow_Nonzero_Ambiguity_Cluster_Count",
    "Shadow_Nonzero_Candidate_Specific_Count", "Shadow_Nonzero_Position_Discrimination",
    "Shadow_Nonzero_Conclusion",
]

MAX_EXCEL_DETAIL_ROWS = 1_048_573


def _as_bool(value: Any, default: bool = True) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number


def intensity_state(value: Any, *, missing: bool = False) -> str:
    if missing or value is None:
        return "missing"
    number = _float(value)
    if number is None:
        return "non_numeric"
    if isnan(number):
        return "nan"
    if number == 0.0:
        return "zero"
    if number < 0.0:
        return "negative"
    return "positive"


def _safe_list(value: Any) -> Any:
    if value is None:
        return []
    try:
        len(value)
        return value
    except TypeError:
        return []


def _spectrum_mode(spectrum: dict[str, Any]) -> str:
    if spectrum.get("centroid spectrum") is not None or "centroid spectrum" in spectrum:
        return "centroid"
    if spectrum.get("profile spectrum") is not None or "profile spectrum" in spectrum:
        return "profile"
    return "unknown"


def capture_source_spectrum(spectrum: dict[str, Any], scan_index: int) -> dict[str, Any]:
    """Capture decoded arrays before RNA_MassHunter numeric conversion."""
    return {
        "spectrum_id": str(spectrum.get("id") or f"scan_{scan_index}"),
        "scan_index": scan_index,
        "original_mz": _safe_list(spectrum.get("m/z array", [])),
        "original_intensity": _safe_list(spectrum.get("intensity array", [])),
        "parsed_mz": [],
        "parsed_intensity": [],
        "annotation_indices": [],
        "spectrum_mode": _spectrum_mode(spectrum),
        "parser_status": "captured_pyteomics_decoded_arrays",
        "parser_error": "",
        "precursor_mz": "",
        "precursor_charge": "",
        "provenance": "pyteomics.MzML decoded binary arrays",
    }


def record_parsed_spectrum(
    record: dict[str, Any], mz_values: Any, intensity_values: Any,
    annotation_indices: Any = (), spectrum_info: Any | None = None,
) -> None:
    record["parsed_mz"] = [float(value) for value in mz_values]
    record["parsed_intensity"] = [float(value) for value in intensity_values]
    record["annotation_indices"] = [int(value) for value in annotation_indices]
    record["parser_status"] = "parsed_and_annotation_input_recorded"
    if spectrum_info is not None:
        record["precursor_mz"] = getattr(spectrum_info, "precursor_mz", "")
        record["precursor_charge"] = getattr(spectrum_info, "precursor_charge", "")


def record_parser_error(record: dict[str, Any], error: Any) -> None:
    record["parser_status"] = "numeric_conversion_failed"
    record["parser_error"] = str(error)


def _candidate_position(value: Any) -> str:
    number = _float(value)
    if number is None or not isfinite(number) or number <= 0:
        return ""
    return str(int(number))


def _candidate_key(modification: Any, parent: Any, position: Any) -> tuple[str, str, str]:
    return str(modification or ""), str(parent or ""), _candidate_position(position)


def _candidate_key_text(key: tuple[str, str, str]) -> str:
    return "|".join(value or "NA" for value in key)


def _raw_peak_key(spectrum_id: Any, mz: Any, rt: Any = None) -> str:
    return physical_observed_peak_key({"Spectrum_ID": spectrum_id, "Observed_mz": mz, "RT": rt})


def _joined(values: Any) -> str:
    return ";".join(sorted({str(value) for value in values if value not in (None, "")}))


def _ion_number(row: dict[str, Any]) -> str:
    for name in ("Ion_Number", "Ion_End", "Ion_Start"):
        value = row.get(name)
        if value not in (None, ""):
            return str(value)
    return ""


def _origin_for_values(original: Any, parsed: Any, annotation: Any, *, missing: bool = False) -> str:
    original_state = intensity_state(original, missing=missing)
    parsed_state = intensity_state(parsed, missing=parsed in (None, ""))
    annotation_state = intensity_state(annotation, missing=annotation in (None, ""))
    if original_state == "zero":
        return "present_in_raw_mzml"
    if parsed_state == "zero" and original_state != "zero":
        return "introduced_during_parsing"
    if annotation_state == "zero" and parsed_state != "zero":
        return "introduced_during_annotation_preparation"
    return "unresolved"


def _overall_origin(categories: set[str]) -> str:
    informative = {value for value in categories if value != "unresolved"}
    if len(informative) == 1:
        return next(iter(informative))
    if len(informative) > 1:
        return "mixed_origin"
    return "unresolved"


def _origin_basis(category: str) -> str:
    if category == "present_in_raw_mzml":
        return (
            "Exact zeros are present in pyteomics-decoded mzML intensity arrays before "
            "RNA_MassHunter float conversion or annotation filtering; raw base64 bytes were not independently decoded."
        )
    if category == "introduced_during_parsing":
        return "A nonzero or nonnumeric source value became exactly zero during numeric parsing."
    if category == "introduced_during_annotation_preparation":
        return "A parsed nonzero value became zero in the annotation input list."
    if category == "mixed_origin":
        return "Zero-valued peaks have more than one observed transformation origin."
    if category == "report_display_only":
        return "Zero appears only in report values and not in captured source, parsed, or annotation values."
    return "Captured stages do not provide enough evidence to classify zero origin."


def _usage_maps(
    spectra: list[Any], ranking_rows: list[dict[str, Any]],
    ion_matches: list[dict[str, Any]], modified_matches: list[dict[str, Any]],
    identity_assignments: list[dict[str, Any]], localization_rows: list[dict[str, Any]],
    ambiguity_clusters: list[dict[str, Any]], ambiguity_details: list[dict[str, Any]],
) -> tuple[dict[str, dict[str, Any]], dict[tuple[str, str, str], dict[str, Any]]]:
    usage: dict[str, dict[str, Any]] = defaultdict(lambda: {
        "theoretical": False, "best": False, "identity": False, "localization": False,
        "clusters": set(), "candidates": set(), "mods": set(), "positions": set(),
        "series": set(), "numbers": set(), "errors": [], "within": [],
    })
    candidate_usage: dict[tuple[str, str, str], dict[str, Any]] = defaultdict(lambda: {
        "match_keys": set(), "best_keys": set(), "identity_keys": set(),
        "localization_keys": set(), "cluster_keys": defaultdict(set),
        "candidate_specific_keys": set(), "position_discrimination_keys": set(),
    })
    spectrum_rt = {str(getattr(item, "spectrum_id", "")): getattr(item, "rt", None) for item in spectra}
    parent_position_to_trna: dict[tuple[str, str, str], str] = {}
    for row in ranking_rows:
        mod = str(row.get("Modification_ID") or "")
        parent = str(row.get("Parent_Fragment_ID") or "")
        parent_pos = _candidate_position(row.get("Candidate_Position_In_Parent") or row.get("Candidate_Positions_In_Parent"))
        trna_pos = _candidate_position(row.get("Candidate_tRNA_Position") or row.get("Candidate_Positions_In_tRNA"))
        if mod and parent and parent_pos:
            parent_position_to_trna[(mod, parent, parent_pos)] = trna_pos

    localized_candidates = {
        _candidate_key(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), row.get("Candidate_Modification_Position_In_tRNA"))
        for row in localization_rows if int(_float(row.get("Num_Modified_Ion_Matches")) or 0) > 0
    }

    for row in ion_matches:
        peak_key = _raw_peak_key(row.get("Spectrum_ID"), row.get("Observed_mz"), row.get("RT") or spectrum_rt.get(str(row.get("Spectrum_ID") or "")))
        item = usage[peak_key]
        item["theoretical"] = True
        item["best"] = True
        item["series"].add(row.get("Best_Ion_Type"))
        item["numbers"].add(_ion_number(row))
        item["errors"].append(row.get("Mass_Error_ppm"))
        item["within"].append(True)

    for row in modified_matches:
        spectrum_id = str(row.get("Spectrum_ID") or "")
        peak_key = _raw_peak_key(spectrum_id, row.get("Observed_mz"), row.get("RT") or spectrum_rt.get(spectrum_id))
        parent_pos = _candidate_position(row.get("Candidate_Modification_Position_In_Parent"))
        trna_pos = parent_position_to_trna.get((str(row.get("Modification_ID") or ""), str(row.get("Parent_Fragment_ID") or ""), parent_pos), parent_pos)
        key = _candidate_key(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), trna_pos)
        item = usage[peak_key]
        item["theoretical"] = True
        item["best"] = True
        item["candidates"].add(key)
        item["mods"].add(row.get("Modification_ID"))
        item["positions"].add(trna_pos)
        item["series"].add(row.get("Ion_Type"))
        item["numbers"].add(_ion_number(row))
        item["errors"].append(row.get("Mass_Error_ppm"))
        item["within"].append(True)
        candidate_usage[key]["match_keys"].add(peak_key)
        candidate_usage[key]["best_keys"].add(peak_key)
        if row.get("Discriminates_Position"):
            candidate_usage[key]["position_discrimination_keys"].add(peak_key)
        if key in localized_candidates:
            item["localization"] = True
            candidate_usage[key]["localization_keys"].add(peak_key)

    for row in identity_assignments:
        peak_key = str(row.get("Physical_Observed_Peak_Key") or "")
        key = _candidate_key(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), row.get("Candidate_tRNA_Position"))
        if not peak_key:
            continue
        item = usage[peak_key]
        item["identity"] = True
        item["candidates"].add(key)
        item["mods"].add(row.get("Modification_ID"))
        item["positions"].add(_candidate_position(row.get("Candidate_tRNA_Position")))
        item["errors"].append(row.get("Mass_Error_ppm"))
        candidate_usage[key]["identity_keys"].add(peak_key)

    cluster_lookup: dict[str, dict[str, Any]] = {}
    for row in ambiguity_clusters:
        cluster_id = str(row.get("Peak_Cluster_ID") or "")
        key = _candidate_key(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), row.get("Candidate_tRNA_Position"))
        cluster_lookup[cluster_id] = {"candidate": key, "row": row}
    for row in ambiguity_details:
        peak_key = str(row.get("Physical_Observed_Peak_Key") or "")
        cluster_id = str(row.get("Peak_Cluster_ID") or "")
        cluster = cluster_lookup.get(cluster_id, {})
        key = cluster.get("candidate")
        if not peak_key:
            continue
        item = usage[peak_key]
        item["clusters"].add(cluster_id)
        item["errors"].append(row.get("Error_ppm"))
        item["within"].append(row.get("Within_Formal_Tolerance"))
        if key:
            item["candidates"].add(key)
            item["mods"].add(key[0])
            item["positions"].add(key[2])
            cluster_row = cluster.get("row", {})
            item["series"].add(cluster_row.get("Ion_Series"))
            item["numbers"].add(cluster_row.get("Ion_Number"))
            candidate_usage[key]["cluster_keys"][cluster_id].add(peak_key)
            if row.get("Candidate_Specificity_Status") == "candidate_specific":
                candidate_usage[key]["candidate_specific_keys"].add(peak_key)
    return usage, candidate_usage


def _build_detail_rows(
    source_records: list[dict[str, Any]], spectra: list[Any], usage: dict[str, dict[str, Any]],
    max_rows: int,
) -> tuple[list[dict[str, Any]], dict[str, str], list[dict[str, Any]], dict[str, set[str]], int]:
    spectrum_rt = {str(getattr(item, "spectrum_id", "")): getattr(item, "rt", None) for item in spectra}
    rows: list[dict[str, Any]] = []
    used_rows: list[dict[str, Any]] = []
    state_by_peak: dict[str, str] = {}
    origins_by_spectrum: dict[str, set[str]] = defaultdict(set)
    total_rows = 0
    for record in sorted(source_records, key=lambda item: (int(item.get("scan_index") or 0), str(item.get("spectrum_id") or ""))):
        mz_values = record.get("original_mz", [])
        original_values = record.get("original_intensity", [])
        parsed_values = record.get("parsed_intensity", [])
        annotation_indices = set(record.get("annotation_indices", []))
        total = max(len(mz_values), len(original_values))
        spectrum_id = str(record.get("spectrum_id") or "")
        for index in range(total):
            total_rows += 1
            mz = mz_values[index] if index < len(mz_values) else ""
            missing = index >= len(original_values)
            original = original_values[index] if not missing else ""
            parsed = parsed_values[index] if index < len(parsed_values) else ""
            annotation = parsed if index in annotation_indices and parsed not in (None, "") else ""
            peak_key = _raw_peak_key(spectrum_id, mz, spectrum_rt.get(spectrum_id))
            item = usage.get(peak_key, {})
            state = intensity_state(original, missing=missing)
            origin = _origin_for_values(original, parsed, annotation, missing=missing)
            origins_by_spectrum[spectrum_id].add(origin)
            if item:
                state_by_peak[peak_key] = state
            should_write = len(rows) < max_rows
            if not should_write and not item:
                continue
            errors = [value for value in (_float(value) for value in item.get("errors", [])) if value is not None and isfinite(value)]
            within_values = [value for value in item.get("within", []) if value not in (None, "")]
            transformations = ["pyteomics_decoded_mzml_array"]
            if parsed not in (None, ""):
                transformations.append("numpy_float_conversion")
            transformations.append("retained_for_annotation" if index in annotation_indices else "not_retained_for_annotation")
            detail = {
                "Spectrum_ID": spectrum_id, "Scan_Index": record.get("scan_index"),
                "Peak_Index": index, "_Physical_Observed_Peak_Key": peak_key,
                "Original_MZ": mz, "Original_Intensity": original,
                "Parsed_Intensity": parsed, "Annotation_Input_Intensity": annotation,
                "Intensity_State": state, "Raw_Spectrum_Derived": True,
                "Synthetic_Or_Placeholder": False, "Provenance": record.get("provenance"),
                "Transformation_History": ";".join(transformations),
                "Used_For_Theoretical_Match": bool(item.get("theoretical")),
                "Selected_As_Best_Match": bool(item.get("best")),
                "Used_For_Identity": bool(item.get("identity")),
                "Used_For_Localization": bool(item.get("localization")),
                "Used_In_Ambiguity_Cluster": bool(item.get("clusters")),
                "Ambiguity_Cluster_ID": _joined(item.get("clusters", set())),
                "Candidate_Key": _joined(_candidate_key_text(key) for key in item.get("candidates", set())),
                "Modification": _joined(item.get("mods", set())), "Position": _joined(item.get("positions", set())),
                "Ion_Series": _joined(item.get("series", set())), "Ion_Number": _joined(item.get("numbers", set())),
                "Error_PPM": min(errors, key=abs) if errors else "",
                "Within_Formal_Tolerance": any(_as_bool(value, False) for value in within_values) if within_values else "",
            }
            if should_write:
                rows.append(detail)
            if item:
                used_rows.append(detail)
    return rows, state_by_peak, used_rows, origins_by_spectrum, total_rows


def _spectrum_rows(source_records: list[dict[str, Any]], origins_by_spectrum: dict[str, set[str]], spectra: list[Any]) -> list[dict[str, Any]]:
    spectrum_lookup = {str(getattr(item, "spectrum_id", "")): item for item in spectra}
    rows = []
    for record in sorted(source_records, key=lambda item: (int(item.get("scan_index") or 0), str(item.get("spectrum_id") or ""))):
        spectrum_id = str(record.get("spectrum_id") or "")
        states = Counter()
        numeric: list[float] = []
        positive: list[float] = []
        zero_mz: list[str] = []
        mz_values = record.get("original_mz", [])
        values = record.get("original_intensity", [])
        total = max(len(mz_values), len(values))
        for index in range(total):
            missing = index >= len(values)
            value = values[index] if not missing else None
            state = intensity_state(value, missing=missing)
            states[state] += 1
            number = _float(value)
            if number is not None and isfinite(number):
                numeric.append(number)
                if number > 0:
                    positive.append(number)
            if state == "zero" and index < len(mz_values) and len(zero_mz) < 5:
                zero_mz.append(str(mz_values[index]))
        finite_count = sum(states[name] for name in ("positive", "zero", "negative"))
        spectrum = spectrum_lookup.get(spectrum_id)
        origin = _overall_origin(origins_by_spectrum.get(spectrum_id, {"unresolved"}))
        rows.append({
            "Spectrum_ID": spectrum_id, "Scan_Index": record.get("scan_index"),
            "Precursor_MZ": getattr(spectrum, "precursor_mz", record.get("precursor_mz", "")),
            "Precursor_Charge": getattr(spectrum, "precursor_charge", record.get("precursor_charge", "")),
            "Total_Peak_Count": total, "Positive_Intensity_Count": states["positive"],
            "Zero_Intensity_Count": states["zero"], "Negative_Intensity_Count": states["negative"],
            "NaN_Intensity_Count": states["nan"], "Missing_Intensity_Count": states["missing"],
            "Non_Numeric_Intensity_Count": states["non_numeric"], "Finite_Intensity_Count": finite_count,
            "Positive_Intensity_Fraction": states["positive"] / total if total else 0.0,
            "Zero_Intensity_Fraction": states["zero"] / total if total else 0.0,
            "Min_Intensity": min(numeric) if numeric else "", "Max_Intensity": max(numeric) if numeric else "",
            "Median_Intensity": median(numeric) if numeric else "",
            "Median_Positive_Intensity": median(positive) if positive else "",
            "Total_Ion_Intensity": sum(numeric), "MZ_Array_Length": len(mz_values),
            "Intensity_Array_Length": len(values), "Array_Length_Mismatch": len(mz_values) != len(values),
            "Example_Zero_Intensity_MZ": ";".join(zero_mz), "Spectrum_Mode": record.get("spectrum_mode", "unknown"),
            "Provenance": record.get("provenance"), "Origin_Category": origin,
            "Origin_Classification_Basis": _origin_basis(origin),
        })
    return rows


def _candidate_rows(
    ranking_rows: list[dict[str, Any]], state_by_peak: dict[str, str],
    candidate_usage: dict[tuple[str, str, str], dict[str, Any]],
) -> list[dict[str, Any]]:
    def is_zero(key: str) -> bool:
        return state_by_peak.get(key) == "zero"

    def is_positive(key: str) -> bool:
        return state_by_peak.get(key) == "positive"

    keys = {
        _candidate_key(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), row.get("Candidate_tRNA_Position") or row.get("Candidate_Positions_In_tRNA"))
        for row in ranking_rows
    } | set(candidate_usage)
    rows = []
    for key in sorted(keys):
        item = candidate_usage.get(key, {})
        match_keys = set(item.get("match_keys", set()))
        best_keys = set(item.get("best_keys", set()))
        identity_keys = set(item.get("identity_keys", set()))
        localization_keys = set(item.get("localization_keys", set()))
        cluster_keys = item.get("cluster_keys", {})
        zero_clusters = {cluster for cluster, peaks in cluster_keys.items() if any(is_zero(peak) for peak in peaks)}
        all_zero_clusters = {cluster for cluster, peaks in cluster_keys.items() if peaks and all(is_zero(peak) for peak in peaks)}
        positive_clusters = {cluster for cluster, peaks in cluster_keys.items() if sum(is_positive(peak) for peak in peaks) >= 2}
        all_evidence_keys = match_keys | identity_keys | localization_keys | {peak for peaks in cluster_keys.values() for peak in peaks}
        positive_support = {peak for peak in all_evidence_keys if is_positive(peak)}
        zero_best = {peak for peak in best_keys if is_zero(peak)}
        zero_major = {peak for peak in identity_keys | localization_keys if is_zero(peak)}
        affected = bool(zero_clusters or zero_best or zero_major)
        if not affected:
            dependence, severity, recommendation = "none", "none", "no_action_needed"
        elif (all_zero_clusters and not positive_support) or (zero_major and not positive_support):
            dependence, severity, recommendation = "all_zero_evidence", "high", "evidence_depends_on_zero_intensity_peaks"
        elif positive_support:
            dependence, severity, recommendation = "mixed_zero_and_positive", "moderate", "positive_intensity_evidence_available"
        else:
            dependence, severity, recommendation = "zero_unmatched_or_cluster_only", "low", "retain_for_diagnostics_only"
        if key[0] == "cnm5U" and affected and not positive_support:
            severity = "high"
            recommendation = "positional_isomer_evidence_depends_on_zero_intensity"
        nonzero_matches = sum(is_positive(peak) for peak in match_keys)
        nonzero_best = sum(is_positive(peak) for peak in best_keys)
        nonzero_clusters = len(positive_clusters)
        nonzero_specific = sum(is_positive(peak) for peak in item.get("candidate_specific_keys", set()))
        nonzero_discrimination = sum(is_positive(peak) for peak in item.get("position_discrimination_keys", set()))
        if not affected:
            conclusion = "unchanged_no_zero_intensity_involvement"
        elif nonzero_discrimination:
            conclusion = "positive_intensity_position_discrimination_remains"
        elif key[0] == "cnm5U":
            conclusion = "positional_isomer_remains_unresolved_without_zero_intensity_peaks"
        elif nonzero_matches or positive_support:
            conclusion = "some_positive_intensity_evidence_remains"
        else:
            conclusion = "no_positive_intensity_evidence_remains_in_existing_results"
        rows.append({
            "Modification_ID": key[0], "Parent_Fragment_ID": key[1], "Candidate_tRNA_Position": key[2],
            "Zero_Intensity_Affected": affected, "Zero_Intensity_Cluster_Count": len(zero_clusters),
            "All_Zero_Cluster_Count": len(all_zero_clusters), "Zero_Intensity_Best_Match_Count": len(zero_best),
            "Positive_Intensity_Support_Count": len(positive_support), "Zero_Intensity_Dependence": dependence,
            "Zero_Intensity_Audit_Severity": severity, "Zero_Intensity_Audit_Recommendation": recommendation,
            "Zero_Intensity_Audit_Applied_To_Final_Score": False,
            "Shadow_Nonzero_Peak_Match_Count": nonzero_matches,
            "Shadow_Nonzero_Best_Match_Count": nonzero_best,
            "Shadow_Nonzero_Ambiguity_Cluster_Count": nonzero_clusters,
            "Shadow_Nonzero_Candidate_Specific_Count": nonzero_specific,
            "Shadow_Nonzero_Position_Discrimination": nonzero_discrimination,
            "Shadow_Nonzero_Conclusion": conclusion,
        })
    return rows


def build_zero_intensity_audit(
    context: dict[str, Any] | None,
    ranking_rows: list[dict[str, Any]] | None = None,
    ion_matches: list[dict[str, Any]] | None = None,
    modified_matches: list[dict[str, Any]] | None = None,
    identity_assignments: list[dict[str, Any]] | None = None,
    localization_rows: list[dict[str, Any]] | None = None,
    ambiguity_clusters: list[dict[str, Any]] | None = None,
    ambiguity_details: list[dict[str, Any]] | None = None,
    *, enabled: bool = True, nonzero_simulation: bool = True,
    max_detail_rows: int = MAX_EXCEL_DETAIL_ROWS,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    context = context or {}
    source_records = list(context.get("source_spectra") or [])
    spectra = list(context.get("spectra") or [])
    ranking_rows = list(ranking_rows or [])
    if not enabled:
        summary = _empty_summary("disabled")
        diagnostics = _diagnostics(summary, [], False)
        return [], [], [summary], [], [diagnostics]
    usage, candidate_usage = _usage_maps(
        spectra, ranking_rows, list(ion_matches or []), list(modified_matches or []),
        list(identity_assignments or []), list(localization_rows or []),
        list(ambiguity_clusters or []), list(ambiguity_details or []),
    )
    detail_limit = max(1, int(max_detail_rows or MAX_EXCEL_DETAIL_ROWS))
    detail_rows, state_by_peak, used_rows, origins_by_spectrum, original_count = _build_detail_rows(
        source_records, spectra, usage, detail_limit,
    )
    spectra_rows = _spectrum_rows(source_records, origins_by_spectrum, spectra)
    candidate_rows = _candidate_rows(ranking_rows, state_by_peak, candidate_usage)
    truncated = original_count > len(detail_rows)
    reason = f"deterministic scan-index/peak-index truncation at {detail_limit} configured Excel rows" if truncated else ""
    summary = _summary(
        spectra_rows, used_rows, state_by_peak, usage, candidate_rows, nonzero_simulation,
        original_count, len(detail_rows), truncated, reason,
    )
    diagnostics = _diagnostics(summary, source_records, True, candidate_rows)
    return spectra_rows, detail_rows, [summary], candidate_rows, [diagnostics]


def _empty_summary(mode: str) -> dict[str, Any]:
    row = {column: 0 for column in SUMMARY_COLUMNS}
    row.update({
        "Likely_Origin_Category": "unresolved", "Origin_Classification_Confidence": "low",
        "Origin_Classification_Basis": "Audit was disabled or raw MS2 arrays were unavailable.",
        "Recommended_Next_Action": "insufficient_information", "Audit_Mode": mode,
        "Applied_To_Final_Score": False, "Shadow_Nonzero_Conclusion": "insufficient_information",
        "Detail_Truncated": False, "Detail_Truncation_Reason": "",
        "Formal_Best_Match_Definition": "Closest theoretical-ion assignment for a peak retained in formal annotation input; raw-only peaks are excluded.",
        "Ambiguity_Diagnostic_Best_Peak_Definition": "Closest raw peak to theoretical m/z within the audit window; diagnostic only and not a formal match.",
    })
    return row


def _summary(
    spectra_rows: list[dict[str, Any]], used_rows: list[dict[str, Any]],
    state_by_peak: dict[str, str], usage: dict[str, dict[str, Any]],
    candidate_rows: list[dict[str, Any]], nonzero_simulation: bool,
    original_count: int, written_count: int, truncated: bool, reason: str,
) -> dict[str, Any]:
    state_counts = Counter()
    for row in spectra_rows:
        state_counts["positive"] += int(row.get("Positive_Intensity_Count") or 0)
        state_counts["zero"] += int(row.get("Zero_Intensity_Count") or 0)
        state_counts["negative"] += int(row.get("Negative_Intensity_Count") or 0)
        state_counts["nan"] += int(row.get("NaN_Intensity_Count") or 0)
        state_counts["missing"] += int(row.get("Missing_Intensity_Count") or 0)
        state_counts["non_numeric"] += int(row.get("Non_Numeric_Intensity_Count") or 0)
    matched = [row for row in used_rows if row.get("Used_For_Theoretical_Match")]
    cluster_states: dict[str, list[str]] = defaultdict(list)
    for peak_key, item in usage.items():
        state = state_by_peak.get(peak_key, "unknown")
        for cluster in item.get("clusters", set()):
            cluster_states[cluster].append(state)
    zero_clusters = {cluster for cluster, states in cluster_states.items() if "zero" in states}
    all_zero_clusters = {cluster for cluster, states in cluster_states.items() if states and all(state == "zero" for state in states)}
    origin_categories = {row.get("Origin_Category") for row in spectra_rows if int(row.get("Zero_Intensity_Count") or 0) > 0}
    origin = _overall_origin({str(value) for value in origin_categories if value})
    affected = [row for row in candidate_rows if row.get("Zero_Intensity_Affected")]
    if origin == "present_in_raw_mzml":
        recommendation, confidence = "inspect_mzml_raw_binary_arrays", "moderate"
    elif origin == "introduced_during_parsing":
        recommendation, confidence = "inspect_parser_zero_fill_behavior", "high"
    elif origin == "introduced_during_annotation_preparation":
        recommendation, confidence = "compare_raw_and_annotation_peak_lists", "high"
    else:
        recommendation, confidence = "inspect_missing_intensity_coercion", "low"
    nonzero_match = sum(row.get("Used_For_Theoretical_Match") and row.get("Intensity_State") == "positive" for row in used_rows)
    nonzero_best = sum(row.get("Selected_As_Best_Match") and row.get("Intensity_State") == "positive" for row in used_rows)
    nonzero_clusters = sum(sum(state == "positive" for state in states) >= 2 for states in cluster_states.values())
    nonzero_specific = sum(int(row.get("Shadow_Nonzero_Candidate_Specific_Count") or 0) for row in candidate_rows)
    nonzero_position = sum(int(row.get("Shadow_Nonzero_Position_Discrimination") or 0) for row in candidate_rows)
    conclusion = (
        "diagnostic_nonzero_filter_leaves_positive_evidence" if nonzero_match or any(int(row.get("Positive_Intensity_Support_Count") or 0) > 0 for row in candidate_rows)
        else "diagnostic_nonzero_filter_leaves_no_existing_ms2_evidence"
    ) if nonzero_simulation else "nonzero_shadow_simulation_disabled"
    return {
        "Total_Spectra": len(spectra_rows), "Total_Raw_Peaks": original_count,
        "Total_Positive_Peaks": state_counts["positive"], "Total_Zero_Intensity_Peaks": state_counts["zero"],
        "Total_Negative_Peaks": state_counts["negative"], "Total_NaN_Peaks": state_counts["nan"],
        "Total_Missing_Peaks": state_counts["missing"], "Total_Non_Numeric_Peaks": state_counts["non_numeric"],
        "Zero_Intensity_Fraction": state_counts["zero"] / original_count if original_count else 0.0,
        "Spectra_With_Zero_Intensity": sum(int(row.get("Zero_Intensity_Count") or 0) > 0 for row in spectra_rows),
        "All_Zero_Spectra": sum(int(row.get("Total_Peak_Count") or 0) > 0 and int(row.get("Zero_Intensity_Count") or 0) == int(row.get("Total_Peak_Count") or 0) for row in spectra_rows),
        "Total_Matched_Peaks": len(matched),
        "Matched_Zero_Intensity_Peaks": sum(row.get("Intensity_State") == "zero" for row in matched),
        "Zero_Intensity_Best_Matches": sum(row.get("Selected_As_Best_Match") and row.get("Intensity_State") == "zero" for row in used_rows),
        "Zero_Intensity_Identity_Assignments": sum(row.get("Used_For_Identity") and row.get("Intensity_State") == "zero" for row in used_rows),
        "Zero_Intensity_Localization_Uses": sum(row.get("Used_For_Localization") and row.get("Intensity_State") == "zero" for row in used_rows),
        "Ambiguity_Clusters_With_Zero_Intensity": len(zero_clusters), "All_Zero_Ambiguity_Clusters": len(all_zero_clusters),
        "Affected_Candidate_Count": len(affected), "Affected_Top50_Count": 0,
        "Affected_cnm5U_Count": sum(row.get("Modification_ID") == "cnm5U" for row in affected),
        "Likely_Origin_Category": origin, "Origin_Classification_Confidence": confidence,
        "Origin_Classification_Basis": _origin_basis(origin), "Recommended_Next_Action": recommendation,
        "Audit_Mode": "shadow_diagnostic_exactly_zero", "Applied_To_Final_Score": False,
        "Shadow_Nonzero_Peak_Match_Count": nonzero_match if nonzero_simulation else "",
        "Shadow_Nonzero_Best_Match_Count": nonzero_best if nonzero_simulation else "",
        "Shadow_Nonzero_Ambiguity_Cluster_Count": nonzero_clusters if nonzero_simulation else "",
        "Shadow_Nonzero_Candidate_Specific_Count": nonzero_specific if nonzero_simulation else "",
        "Shadow_Nonzero_Position_Discrimination": nonzero_position if nonzero_simulation else "",
        "Shadow_Nonzero_Conclusion": conclusion, "Original_Detail_Row_Count": original_count,
        "Written_Detail_Row_Count": written_count, "Detail_Truncated": truncated, "Detail_Truncation_Reason": reason,
        "Formal_Best_Match_Definition": "Closest theoretical-ion assignment for a peak retained in spectrum.peaks after configured intensity filtering; raw-only peaks are excluded.",
        "Ambiguity_Diagnostic_Best_Peak_Definition": "Closest raw peak to theoretical m/z within the audit window, ordered by absolute m/z error, then higher intensity, then lower m/z; it may have zero intensity and is not a formal match.",
    }


def _diagnostics(
    summary: dict[str, Any], source_records: list[dict[str, Any]],
    available: bool, candidate_rows: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    positive: list[float] = []
    annotation_count = 0
    annotation_zero = 0
    for record in source_records:
        values = record.get("original_intensity", [])
        for value in values:
            number = _float(value)
            if number is not None and isfinite(number) and number > 0:
                positive.append(number)
        parsed = record.get("parsed_intensity", [])
        for index in record.get("annotation_indices", []):
            if 0 <= int(index) < len(parsed):
                annotation_count += 1
                annotation_zero += intensity_state(parsed[int(index)]) == "zero"
    severity_counts = Counter(row.get("Zero_Intensity_Audit_Severity") for row in (candidate_rows or []))
    severity = "high" if severity_counts["high"] else "moderate" if severity_counts["moderate"] else "low" if summary.get("Total_Zero_Intensity_Peaks") else "none"
    return {
        "Zero_Intensity_Audit_Available": available, "Raw_Peak_Count": summary.get("Total_Raw_Peaks", 0),
        "Raw_Zero_Intensity_Peak_Count": summary.get("Total_Zero_Intensity_Peaks", 0),
        "Raw_Zero_Intensity_Fraction": summary.get("Zero_Intensity_Fraction", 0.0),
        "Raw_Positive_Intensity_Peak_Count": summary.get("Total_Positive_Peaks", 0),
        "Raw_Median_Positive_Intensity": median(positive) if positive else "",
        "Annotation_Input_Peak_Count": annotation_count, "Annotation_Input_Zero_Intensity_Count": annotation_zero,
        "Zero_Intensity_Matched_Peak_Count": summary.get("Matched_Zero_Intensity_Peaks", 0),
        "Zero_Intensity_Best_Match_Count": summary.get("Zero_Intensity_Best_Matches", 0),
        "Zero_Intensity_Identity_Assignment_Count": summary.get("Zero_Intensity_Identity_Assignments", 0),
        "Zero_Intensity_Localization_Count": summary.get("Zero_Intensity_Localization_Uses", 0),
        "Zero_Intensity_Ambiguity_Cluster_Count": summary.get("Ambiguity_Clusters_With_Zero_Intensity", 0),
        "All_Zero_Ambiguity_Cluster_Count": summary.get("All_Zero_Ambiguity_Clusters", 0),
        "Zero_Intensity_Origin_Category": summary.get("Likely_Origin_Category", "unresolved"),
        "Zero_Intensity_Audit_Severity": severity,
        "Zero_Intensity_Audit_Recommendation": summary.get("Recommended_Next_Action", "no_action_needed"),
        "Zero_Intensity_Audit_Applied_To_Final_Score": False,
    }


def update_top50_affected(summary_rows: list[dict[str, Any]], top_rows: Any) -> None:
    if not summary_rows:
        return
    try:
        records = top_rows.to_dict("records")
    except AttributeError:
        records = list(top_rows or [])
    summary_rows[0]["Affected_Top50_Count"] = sum(_as_bool(row.get("Zero_Intensity_Affected"), False) for row in records)

