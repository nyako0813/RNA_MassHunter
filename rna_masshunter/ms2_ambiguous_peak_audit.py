"""Shadow peak-cluster audit for ambiguous nearby MS2 peaks.

This module describes ambiguity already identified by MS2_Unmatched_Ion_Audit.
It does not rematch peaks or mutate any formal or existing shadow result.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha1
from statistics import median
from typing import Any

from rna_masshunter.ms1_mapping import ppm_error
from rna_masshunter.ms2_identity_evidence import physical_observed_peak_key

CLUSTER_COLUMNS = [
    "Modification_ID", "Parent_Fragment_ID", "Candidate_tRNA_Position",
    "Candidate_Position_In_Parent", "Spectrum_ID", "Theoretical_Ion_ID",
    "Ion_Series", "Ion_Number", "Ion_Charge", "Theoretical_mz",
    "Peak_Cluster_ID", "Peak_Cluster_Size", "Peak_Cluster_mz_Min", "Peak_Cluster_mz_Max",
    "Peak_Cluster_Span_Da", "Peak_Cluster_Span_ppm", "Peak_Cluster_Total_Intensity",
    "Peak_Cluster_Max_Intensity", "Peak_Cluster_Base_Peak_Fraction",
    "Peaks_Below_Theoretical_Count", "Peaks_Above_Theoretical_Count",
    "Peaks_Within_Formal_Tolerance_Count", "Closest_Peak_Rank_By_Error",
    "Closest_Peak_Rank_By_Intensity", "Primary_Ambiguity_Type", "Secondary_Ambiguity_Types",
    "Best_Peak_Key", "Best_Peak_mz", "Best_Peak_Intensity", "Best_Error_ppm",
    "Best_Peak_Assignment_Scope", "Ambiguity_Interpretation", "Ambiguity_Warnings",
]

DETAIL_COLUMNS = [
    "Peak_Cluster_ID", "Physical_Observed_Peak_Key", "Spectrum_ID", "Observed_mz",
    "Intensity", "Relative_Intensity", "Error_Da", "Error_ppm",
    "Within_Formal_Tolerance", "Within_Audit_Window", "Candidate_Modification_Count",
    "Candidate_Position_Count", "Theoretical_Ion_Count", "Structural_Isomer_Group_Count",
    "Assignment_Scope", "Candidate_Specificity_Status", "Competing_Modification_IDs",
    "Competing_Positions", "Competing_Theoretical_Ion_IDs", "Peak_Interpretation", "Peak_Warnings",
]

SUMMARY_COLUMNS = [
    "Rank", "Modification_ID", "Parent_Fragment_ID", "Candidate_tRNA_Position",
    "Ambiguous_Theoretical_Ion_Count", "Ambiguous_Peak_Cluster_Count",
    "Total_Ambiguous_Physical_Peak_Count", "Candidate_Specific_Peak_Count",
    "Position_Group_Shared_Peak_Count", "Structural_Isomer_Shared_Peak_Count",
    "Cross_Candidate_Shared_Peak_Count", "Within_Tolerance_Multiple_Peak_Count",
    "Outside_Tolerance_Cluster_Count", "Median_Cluster_Size", "Maximum_Cluster_Size",
    "Best_Error_ppm", "Primary_Ambiguity_Pattern", "Ambiguity_Severity",
    "Recommended_Followup", "Ambiguity_Warnings",
]

DIAGNOSTIC_COLUMNS = [
    "MS2_Ambiguous_Peak_Audit_Enabled", "Apply_Ambiguous_Peak_Audit_To_Final_Score",
    "Ambiguous_Theoretical_Ion_Count", "Ambiguous_Peak_Cluster_Count",
    "Ambiguous_Physical_Peak_Count", "Candidate_Specific_Peak_Count",
    "Position_Group_Shared_Peak_Count", "Structural_Isomer_Shared_Peak_Count",
    "Cross_Candidate_Shared_Peak_Count", "Multiple_Peaks_Within_Tolerance_Count",
    "Outside_Tolerance_Cluster_Count", "Low_Severity_Count", "Moderate_Severity_Count",
    "High_Severity_Count", "Audit_Window_Rule", "Physical_Peak_Key_Rule",
]

TOP_SHADOW_COLUMNS = [
    "Ambiguous_Theoretical_Ion_Count", "Ambiguous_Peak_Cluster_Count",
    "Maximum_Ambiguous_Cluster_Size", "Primary_Ambiguity_Pattern", "Ambiguity_Severity",
    "Candidate_Specific_Ambiguous_Peak_Count", "Position_Group_Shared_Peak_Count",
    "Structural_Isomer_Shared_Peak_Count", "Cross_Candidate_Shared_Peak_Count",
    "Ambiguous_Peak_Recommended_Followup", "Ambiguous_Peak_Warnings",
]

TYPE_PRIORITY = [
    "multiple_peaks_within_matching_tolerance",
    "multiple_candidate_assignments_same_physical_peak",
    "multiple_theoretical_ions_compete_for_same_peak",
    "structural_isomer_shared_peak", "positional_isomer_shared_peak",
    "multiple_peaks_bracketing_theoretical_mz", "multiple_peaks_same_side",
    "multiple_peaks_outside_tolerance_but_within_audit_window",
    "unresolved_peak_cluster", "insufficient_information",
]

SHARED_STATUS_PRIORITY = {
    "unknown": 0, "spectrum_peak_only": 1, "candidate_specific": 2,
    "position_group_shared": 3, "structural_isomer_group_shared": 4,
    "cross_candidate_shared": 5,
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
    value = _float(value)
    return int(value) if value is not None and value > 0 else None


def _candidate_key(row: dict[str, Any], position_field: str) -> tuple[str, str, int | None]:
    return str(row.get("Modification_ID") or ""), str(row.get("Parent_Fragment_ID") or ""), _position(row.get(position_field))


def _physical_key(spectrum: Any, mz: float, intensity: float) -> str:
    return physical_observed_peak_key({
        "Spectrum_ID": getattr(spectrum, "spectrum_id", ""),
        "Observed_mz": mz, "Observed_Intensity": intensity,
        "RT": getattr(spectrum, "rt", None),
    })


def deterministic_cluster_id(audit: dict[str, Any], lower: float, upper: float) -> str:
    key = "|".join(str(value if value not in (None, "") else "NA") for value in (
        audit.get("Spectrum_ID"), audit.get("Modification_ID"), audit.get("Parent_Fragment_ID"),
        audit.get("Candidate_Position_In_Parent"), audit.get("Theoretical_Ion_ID"),
        audit.get("Ion_Series"), audit.get("Ion_Number"), audit.get("Ion_Charge"),
        f"{float(audit.get('Theoretical_mz')):.8f}", f"{lower:.8f}", f"{upper:.8f}",
    ))
    return f"APC_{sha1(key.encode('utf-8')).hexdigest()[:16].upper()}"


def _ranking_maps(ranking_rows: list[dict[str, Any]] | None):
    structural = {}
    trna = {}
    for row in ranking_rows or []:
        key = _candidate_key(row, "Candidate_Position_In_Parent")
        structural[key] = str(row.get("Structural_Isomer_Group_ID") or "")
        trna[key] = _position(row.get("Candidate_tRNA_Position"))
    return structural, trna


def _ion_candidates(
    spectrum_id: str, observed_mz: float, tolerance_ppm: float,
    ions_by_spectrum: dict[str, list[dict[str, Any]]], structural: dict[tuple[str, str, int | None], str],
) -> tuple[set[tuple[str, str, int | None]], set[str], set[int], set[str], set[str]]:
    candidate_keys: set[tuple[str, str, int | None]] = set()
    modifications: set[str] = set()
    positions: set[int] = set()
    ion_ids: set[str] = set()
    groups: set[str] = set()
    for ion in ions_by_spectrum.get(spectrum_id, []):
        theoretical_mz = _float(ion.get("Theoretical_mz"))
        if theoretical_mz is None or abs(ppm_error(observed_mz, theoretical_mz)) > tolerance_ppm:
            continue
        key = _candidate_key(ion, "Candidate_Modification_Position_In_Parent")
        candidate_keys.add(key)
        modifications.add(key[0])
        if key[2] is not None:
            positions.add(key[2])
        ion_ids.add(str(ion.get("Ion_ID") or ""))
        group_id = structural.get(key, "")
        if group_id:
            groups.add(group_id)
    return candidate_keys, modifications, positions, ion_ids, groups


def _specificity(
    candidate_keys: set[tuple[str, str, int | None]], groups: set[str],
    existing_scopes: set[str], structural: dict[tuple[str, str, int | None], str],
) -> tuple[str, str]:
    if existing_scopes:
        if "cross_candidate_ambiguous" in existing_scopes:
            return "cross_candidate_shared", "cross_candidate_ambiguous"
        if "structural_isomer_group_level" in existing_scopes:
            return "structural_isomer_group_shared", "structural_isomer_group_level"
        if "position_group_level" in existing_scopes:
            return "position_group_shared", "position_group_level"
        if "candidate_specific" in existing_scopes:
            return "candidate_specific", "candidate_specific"
    if not candidate_keys:
        return "spectrum_peak_only", "spectrum_peak_only"
    if len(candidate_keys) == 1:
        return "candidate_specific", "candidate_specific"
    nonempty_groups = {structural.get(key, "") for key in candidate_keys if structural.get(key, "")}
    same_group = len(nonempty_groups) == 1 and len(nonempty_groups) > 0 and all(structural.get(key, "") in nonempty_groups for key in candidate_keys)
    same_mod_parent = len({(key[0], key[1]) for key in candidate_keys}) == 1
    if same_group:
        return "structural_isomer_group_shared", "structural_isomer_group_level"
    if same_mod_parent:
        return "position_group_shared", "position_group_level"
    return "cross_candidate_shared", "cross_candidate_ambiguous"


def _ambiguity_types(cluster: dict[str, Any], details: list[dict[str, Any]]) -> tuple[str, str]:
    types: set[str] = set()
    below = int(cluster.get("Peaks_Below_Theoretical_Count") or 0)
    above = int(cluster.get("Peaks_Above_Theoretical_Count") or 0)
    within = int(cluster.get("Peaks_Within_Formal_Tolerance_Count") or 0)
    if below and above:
        types.add("multiple_peaks_bracketing_theoretical_mz")
    elif below + above >= 2:
        types.add("multiple_peaks_same_side")
    if within > 1:
        types.add("multiple_peaks_within_matching_tolerance")
    if details and within == 0:
        types.add("multiple_peaks_outside_tolerance_but_within_audit_window")
    if any(int(row.get("Candidate_Modification_Count") or 0) > 1 or row.get("Candidate_Specificity_Status") == "cross_candidate_shared" for row in details):
        types.add("multiple_candidate_assignments_same_physical_peak")
    if any(int(row.get("Theoretical_Ion_Count") or 0) > 1 for row in details):
        types.add("multiple_theoretical_ions_compete_for_same_peak")
    if any(row.get("Candidate_Specificity_Status") == "structural_isomer_group_shared" for row in details):
        types.add("structural_isomer_shared_peak")
    if any(row.get("Candidate_Specificity_Status") == "position_group_shared" for row in details):
        types.add("positional_isomer_shared_peak")
    if not types:
        types.add("unresolved_peak_cluster" if details else "insufficient_information")
    ordered = [item for item in TYPE_PRIORITY if item in types]
    return ordered[0], ";".join(ordered[1:])


def _severity(clusters: list[dict[str, Any]], details: list[dict[str, Any]]) -> str:
    if not clusters:
        return "unknown"
    max_size = max(int(row.get("Peak_Cluster_Size") or 0) for row in clusters)
    statuses = {str(row.get("Candidate_Specificity_Status") or "unknown") for row in details}
    high = (
        "cross_candidate_shared" in statuses
        or any(int(row.get("Theoretical_Ion_Count") or 0) > 1 for row in details)
        or any(int(row.get("Peaks_Within_Formal_Tolerance_Count") or 0) > 1 for row in clusters)
        or max_size >= 4
    )
    if high:
        return "high"
    if statuses & {"position_group_shared", "structural_isomer_group_shared"} or 2 <= max_size <= 3:
        if max_size == 2 and "candidate_specific" in statuses and not statuses & {"position_group_shared", "structural_isomer_group_shared"}:
            return "low"
        return "moderate"
    if max_size == 2 and "candidate_specific" in statuses:
        return "low"
    return "unknown"


def _followup(severity: str, details: list[dict[str, Any]], clusters: list[dict[str, Any]]) -> str:
    statuses = {str(row.get("Candidate_Specificity_Status") or "") for row in details}
    if severity == "unknown":
        return "insufficient_information"
    if "cross_candidate_shared" in statuses:
        return "candidate_not_distinguishable_with_current_peaks"
    if "structural_isomer_group_shared" in statuses:
        return "structural_isomer_remains_unresolved"
    if "position_group_shared" in statuses:
        return "positional_isomer_remains_unresolved"
    if severity == "high" and any(int(row.get("Peaks_Within_Formal_Tolerance_Count") or 0) > 1 for row in clusters):
        return "acquire_higher_resolution_ms2"
    if severity == "high":
        return "inspect_raw_spectrum_manually"
    if severity == "moderate":
        return "acquire_additional_ms2_spectrum"
    return "no_action_low_severity"


def build_ambiguous_peak_audit(
    spectra: list[Any] | None,
    unmatched_audit_rows: list[dict[str, Any]] | None,
    modified_ions: list[dict[str, Any]] | None,
    ranking_rows: list[dict[str, Any]] | None,
    identity_assignment_rows: list[dict[str, Any]] | None,
    enabled: bool = True,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Build theoretical-ion cluster rows and per-physical-peak detail rows."""
    if not enabled:
        return [], []
    spectra_lookup = {str(item.spectrum_id): item for item in (spectra or [])}
    structural, _ = _ranking_maps(ranking_rows)
    ions_by_spectrum: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for ion in modified_ions or []:
        if _bool(ion.get("Ion_Contains_Modification")):
            ions_by_spectrum[str(ion.get("Spectrum_ID") or "")].append(ion)
    assignments_by_peak: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for assignment in identity_assignment_rows or []:
        assignments_by_peak[str(assignment.get("Physical_Observed_Peak_Key") or "")].append(assignment)

    clusters: list[dict[str, Any]] = []
    all_details: list[dict[str, Any]] = []
    for audit in unmatched_audit_rows or []:
        if audit.get("Unmatched_Reason_Status") != "ambiguous_multiple_nearby_peaks":
            continue
        spectrum = spectra_lookup.get(str(audit.get("Spectrum_ID") or ""))
        theoretical_mz = _float(audit.get("Theoretical_mz"))
        window = _float(audit.get("Audit_Search_Window_Da"))
        tolerance_da = _float(audit.get("Matching_Tolerance_Da"))
        tolerance_ppm = _float(audit.get("Matching_Tolerance_ppm"))
        if spectrum is None or theoretical_mz is None or window is None or tolerance_da is None or tolerance_ppm is None or getattr(spectrum, "raw_peaks", None) is None:
            cluster_id = deterministic_cluster_id(audit, theoretical_mz or 0.0, theoretical_mz or 0.0)
            clusters.append({
                "Modification_ID": audit.get("Modification_ID"), "Parent_Fragment_ID": audit.get("Parent_Fragment_ID"),
                "Candidate_tRNA_Position": audit.get("Candidate_tRNA_Position"),
                "Candidate_Position_In_Parent": audit.get("Candidate_Position_In_Parent"),
                "Spectrum_ID": audit.get("Spectrum_ID"), "Theoretical_Ion_ID": audit.get("Theoretical_Ion_ID"),
                "Ion_Series": audit.get("Ion_Series"), "Ion_Number": audit.get("Ion_Number"),
                "Ion_Charge": audit.get("Ion_Charge"), "Theoretical_mz": audit.get("Theoretical_mz"),
                "Peak_Cluster_ID": cluster_id, "Peak_Cluster_Size": 0, "Peak_Cluster_mz_Min": "", "Peak_Cluster_mz_Max": "",
                "Peak_Cluster_Span_Da": "", "Peak_Cluster_Span_ppm": "", "Peak_Cluster_Total_Intensity": 0,
                "Peak_Cluster_Max_Intensity": 0, "Peak_Cluster_Base_Peak_Fraction": "",
                "Peaks_Below_Theoretical_Count": 0, "Peaks_Above_Theoretical_Count": 0,
                "Peaks_Within_Formal_Tolerance_Count": 0, "Closest_Peak_Rank_By_Error": "",
                "Closest_Peak_Rank_By_Intensity": "", "Primary_Ambiguity_Type": "insufficient_information",
                "Secondary_Ambiguity_Types": "", "Best_Peak_Key": "", "Best_Peak_mz": "",
                "Best_Peak_Intensity": "", "Best_Error_ppm": "", "Best_Peak_Assignment_Scope": "unknown",
                "Ambiguity_Interpretation": "Raw spectrum or tolerance information was unavailable.",
                "Ambiguity_Warnings": "information unavailable for peak cluster reconstruction",
            })
            continue

        peaks = sorted(
            [(float(mz), float(intensity)) for mz, intensity in spectrum.raw_peaks if abs(float(mz) - theoretical_mz) <= window],
            key=lambda item: (abs(item[0] - theoretical_mz), -item[1], item[0]),
        )
        lower = theoretical_mz - window
        upper = theoretical_mz + window
        cluster_id = deterministic_cluster_id(audit, lower, upper)
        details: list[dict[str, Any]] = []
        for mz, intensity in peaks:
            peak_key = _physical_key(spectrum, mz, intensity)
            candidates, modifications, positions, ion_ids, groups = _ion_candidates(
                str(audit.get("Spectrum_ID") or ""), mz, tolerance_ppm, ions_by_spectrum, structural,
            )
            assignments = assignments_by_peak.get(peak_key, [])
            existing_scopes = {str(row.get("Evidence_Scope") or "") for row in assignments if row.get("Evidence_Scope")}
            for assignment in assignments:
                key = _candidate_key(assignment, "Candidate_Position_In_Parent")
                candidates.add(key)
                modifications.add(key[0])
                if key[2] is not None:
                    positions.add(key[2])
                if assignment.get("Theoretical_Ion_ID"):
                    ion_ids.add(str(assignment.get("Theoretical_Ion_ID")))
                group_id = str(assignment.get("Structural_Isomer_Group_ID") or "")
                if group_id:
                    groups.add(group_id)
            specificity, scope = _specificity(candidates, groups, existing_scopes, structural)
            error_da = mz - theoretical_mz
            error_ppm = ppm_error(mz, theoretical_mz)
            relative = intensity / float(getattr(spectrum, "base_peak_intensity", None) or intensity or 1.0)
            details.append({
                "Peak_Cluster_ID": cluster_id, "Physical_Observed_Peak_Key": peak_key,
                "Spectrum_ID": audit.get("Spectrum_ID"), "Observed_mz": mz, "Intensity": intensity,
                "Relative_Intensity": relative, "Error_Da": error_da, "Error_ppm": error_ppm,
                "Within_Formal_Tolerance": abs(error_da) <= tolerance_da, "Within_Audit_Window": True,
                "Candidate_Modification_Count": len(modifications), "Candidate_Position_Count": len(positions),
                "Theoretical_Ion_Count": len(ion_ids), "Structural_Isomer_Group_Count": len(groups),
                "Assignment_Scope": scope, "Candidate_Specificity_Status": specificity,
                "Competing_Modification_IDs": ";".join(sorted(item for item in modifications if item)),
                "Competing_Positions": ";".join(str(item) for item in sorted(positions)),
                "Competing_Theoretical_Ion_IDs": ";".join(sorted(item for item in ion_ids if item)),
                "Peak_Interpretation": "Potential formal-tolerance candidate competition." if candidates else "Raw spectrum peak in audit window without a formal-tolerance theoretical candidate.",
                "Peak_Warnings": "" if candidates else "spectrum peak only; no candidate assignment inferred",
            })
        if not peaks:
            # This should not occur for an existing ambiguous row, but retain a diagnostic row if inputs diverge.
            primary, secondary = "insufficient_information", ""
            best = None
        else:
            provisional = {
                "Peaks_Below_Theoretical_Count": sum(mz < theoretical_mz for mz, _ in peaks),
                "Peaks_Above_Theoretical_Count": sum(mz > theoretical_mz for mz, _ in peaks),
                "Peaks_Within_Formal_Tolerance_Count": sum(abs(mz - theoretical_mz) <= tolerance_da for mz, _ in peaks),
            }
            primary, secondary = _ambiguity_types(provisional, details)
            best = details[0]
        mz_values = [mz for mz, _ in peaks]
        intensities = [intensity for _, intensity in peaks]
        closest_intensity_rank = ""
        if best is not None:
            intensity_order = sorted(details, key=lambda row: (-float(row["Intensity"]), abs(float(row["Error_Da"])), float(row["Observed_mz"])))
            closest_intensity_rank = next(index for index, row in enumerate(intensity_order, start=1) if row["Physical_Observed_Peak_Key"] == best["Physical_Observed_Peak_Key"])
        cluster = {
            "Modification_ID": audit.get("Modification_ID"), "Parent_Fragment_ID": audit.get("Parent_Fragment_ID"),
            "Candidate_tRNA_Position": audit.get("Candidate_tRNA_Position"),
            "Candidate_Position_In_Parent": audit.get("Candidate_Position_In_Parent"),
            "Spectrum_ID": audit.get("Spectrum_ID"), "Theoretical_Ion_ID": audit.get("Theoretical_Ion_ID"),
            "Ion_Series": audit.get("Ion_Series"), "Ion_Number": audit.get("Ion_Number"),
            "Ion_Charge": audit.get("Ion_Charge"), "Theoretical_mz": audit.get("Theoretical_mz"),
            "Peak_Cluster_ID": cluster_id, "Peak_Cluster_Size": len(peaks),
            "Peak_Cluster_mz_Min": min(mz_values) if mz_values else "", "Peak_Cluster_mz_Max": max(mz_values) if mz_values else "",
            "Peak_Cluster_Span_Da": max(mz_values) - min(mz_values) if mz_values else "",
            "Peak_Cluster_Span_ppm": ((max(mz_values) - min(mz_values)) / theoretical_mz * 1_000_000) if mz_values else "",
            "Peak_Cluster_Total_Intensity": sum(intensities), "Peak_Cluster_Max_Intensity": max(intensities) if intensities else 0,
            "Peak_Cluster_Base_Peak_Fraction": (max(intensities) / float(getattr(spectrum, "base_peak_intensity", None) or max(intensities))) if intensities else "",
            "Peaks_Below_Theoretical_Count": sum(mz < theoretical_mz for mz, _ in peaks),
            "Peaks_Above_Theoretical_Count": sum(mz > theoretical_mz for mz, _ in peaks),
            "Peaks_Within_Formal_Tolerance_Count": sum(abs(mz - theoretical_mz) <= tolerance_da for mz, _ in peaks),
            "Closest_Peak_Rank_By_Error": 1 if best is not None else "", "Closest_Peak_Rank_By_Intensity": closest_intensity_rank,
            "Primary_Ambiguity_Type": primary, "Secondary_Ambiguity_Types": secondary,
            "Best_Peak_Key": best.get("Physical_Observed_Peak_Key", "") if best else "",
            "Best_Peak_mz": best.get("Observed_mz", "") if best else "",
            "Best_Peak_Intensity": best.get("Intensity", "") if best else "",
            "Best_Error_ppm": best.get("Error_ppm", "") if best else "",
            "Best_Peak_Assignment_Scope": best.get("Assignment_Scope", "unknown") if best else "unknown",
            "Ambiguity_Interpretation": f"cluster size={len(peaks)}; primary={primary}; formal-tolerance peaks={sum(abs(mz - theoretical_mz) <= tolerance_da for mz, _ in peaks)}",
            "Ambiguity_Warnings": "shadow-only cluster; existing formal match remains unchanged",
        }
        clusters.append(cluster)
        all_details.extend(details)
    return clusters, all_details


def build_ambiguity_summary(
    ranking_rows: list[dict[str, Any]] | None,
    clusters: list[dict[str, Any]] | None,
    details: list[dict[str, Any]] | None,
) -> list[dict[str, Any]]:
    clusters_by_candidate: dict[tuple[str, str, int | None], list[dict[str, Any]]] = defaultdict(list)
    details_by_cluster: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in details or []:
        details_by_cluster[str(row.get("Peak_Cluster_ID") or "")].append(row)
    for row in clusters or []:
        clusters_by_candidate[_candidate_key(row, "Candidate_Position_In_Parent")].append(row)
    result = []
    for ranking in ranking_rows or []:
        key = _candidate_key(ranking, "Candidate_Position_In_Parent")
        group = clusters_by_candidate.get(key, [])
        group_details = [item for cluster in group for item in details_by_cluster.get(str(cluster.get("Peak_Cluster_ID") or ""), [])]
        unique_by_status: dict[str, set[str]] = defaultdict(set)
        for item in group_details:
            unique_by_status[str(item.get("Candidate_Specificity_Status") or "unknown")].add(str(item.get("Physical_Observed_Peak_Key") or ""))
        physical_keys = {str(item.get("Physical_Observed_Peak_Key") or "") for item in group_details if item.get("Physical_Observed_Peak_Key")}
        sizes = [int(item.get("Peak_Cluster_Size") or 0) for item in group]
        errors = [abs(value) for item in group if (value := _float(item.get("Best_Error_ppm"))) is not None]
        patterns = Counter(str(item.get("Primary_Ambiguity_Type") or "insufficient_information") for item in group)
        primary = sorted(patterns, key=lambda item: (-patterns[item], TYPE_PRIORITY.index(item) if item in TYPE_PRIORITY else 999, item))[0] if patterns else "insufficient_information"
        severity = _severity(group, group_details)
        warnings = []
        if not group:
            warnings.append("no ambiguous nearby-peak audit rows for candidate")
        if any(item.get("Primary_Ambiguity_Type") == "insufficient_information" for item in group):
            warnings.append("some peak clusters lack required raw information")
        result.append({
            "Rank": ranking.get("Rank"), "Modification_ID": ranking.get("Modification_ID"),
            "Parent_Fragment_ID": ranking.get("Parent_Fragment_ID"),
            "Candidate_tRNA_Position": ranking.get("Candidate_tRNA_Position"),
            "Ambiguous_Theoretical_Ion_Count": len(group), "Ambiguous_Peak_Cluster_Count": len(group),
            "Total_Ambiguous_Physical_Peak_Count": len(physical_keys),
            "Candidate_Specific_Peak_Count": len(unique_by_status["candidate_specific"]),
            "Position_Group_Shared_Peak_Count": len(unique_by_status["position_group_shared"]),
            "Structural_Isomer_Shared_Peak_Count": len(unique_by_status["structural_isomer_group_shared"]),
            "Cross_Candidate_Shared_Peak_Count": len(unique_by_status["cross_candidate_shared"]),
            "Within_Tolerance_Multiple_Peak_Count": sum(int(item.get("Peaks_Within_Formal_Tolerance_Count") or 0) > 1 for item in group),
            "Outside_Tolerance_Cluster_Count": sum(int(item.get("Peaks_Within_Formal_Tolerance_Count") or 0) == 0 and int(item.get("Peak_Cluster_Size") or 0) > 0 for item in group),
            "Median_Cluster_Size": median(sizes) if sizes else "", "Maximum_Cluster_Size": max(sizes) if sizes else 0,
            "Best_Error_ppm": min(errors) if errors else "", "Primary_Ambiguity_Pattern": primary,
            "Ambiguity_Severity": severity, "Recommended_Followup": _followup(severity, group_details, group),
            "Ambiguity_Warnings": "; ".join(warnings),
        })
    return result


def build_ambiguity_diagnostics(
    clusters: list[dict[str, Any]] | None,
    details: list[dict[str, Any]] | None,
    summaries: list[dict[str, Any]] | None,
    enabled: bool = True,
) -> list[dict[str, Any]]:
    clusters = clusters or []
    details = details or []
    summaries = summaries or []
    physical_by_status: dict[str, set[str]] = defaultdict(set)
    for row in details:
        physical_by_status[str(row.get("Candidate_Specificity_Status") or "unknown")].add(str(row.get("Physical_Observed_Peak_Key") or ""))
    severities = Counter(str(row.get("Ambiguity_Severity") or "unknown") for row in summaries if int(row.get("Ambiguous_Theoretical_Ion_Count") or 0) > 0)
    return [{
        "MS2_Ambiguous_Peak_Audit_Enabled": bool(enabled),
        "Apply_Ambiguous_Peak_Audit_To_Final_Score": False,
        "Ambiguous_Theoretical_Ion_Count": len(clusters), "Ambiguous_Peak_Cluster_Count": len(clusters),
        "Ambiguous_Physical_Peak_Count": len({str(row.get("Physical_Observed_Peak_Key") or "") for row in details if row.get("Physical_Observed_Peak_Key")}),
        "Candidate_Specific_Peak_Count": len(physical_by_status["candidate_specific"]),
        "Position_Group_Shared_Peak_Count": len(physical_by_status["position_group_shared"]),
        "Structural_Isomer_Shared_Peak_Count": len(physical_by_status["structural_isomer_group_shared"]),
        "Cross_Candidate_Shared_Peak_Count": len(physical_by_status["cross_candidate_shared"]),
        "Multiple_Peaks_Within_Tolerance_Count": sum(int(row.get("Peaks_Within_Formal_Tolerance_Count") or 0) > 1 for row in clusters),
        "Outside_Tolerance_Cluster_Count": sum(int(row.get("Peaks_Within_Formal_Tolerance_Count") or 0) == 0 and int(row.get("Peak_Cluster_Size") or 0) > 0 for row in clusters),
        "Low_Severity_Count": severities["low"], "Moderate_Severity_Count": severities["moderate"],
        "High_Severity_Count": severities["high"],
        "Audit_Window_Rule": "reuse MS2_Unmatched_Ion_Audit per-ion audit window",
        "Physical_Peak_Key_Rule": "reuse ms2_identity_evidence.physical_observed_peak_key",
    }]
