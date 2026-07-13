"""Four-stage shadow classification of existing MS2 ambiguity clusters.

This module reuses existing raw-cluster, formal-match, identity, and localization
rows.  It does not rematch peaks or alter any formal or existing shadow result.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from typing import Any

from rna_masshunter.ms2_identity_evidence import physical_observed_peak_key

CLUSTER_COLUMNS = [
    "Ambiguity_Cluster_ID", "Candidate_Key", "Modification", "Modification_Family",
    "Position", "Theoretical_Ion_ID", "Ion_Series", "Ion_Number", "Theoretical_MZ",
    "Existing_Primary_Ambiguity_Type", "Existing_Severity",
    "Raw_Peak_Count", "Raw_Zero_Intensity_Peak_Count", "Raw_Positive_Intensity_Peak_Count", "Raw_Ambiguous",
    "Positive_Peak_Count", "Positive_Ambiguous", "Positive_Peak_MZ_List",
    "Positive_Peak_Intensity_List", "Positive_Peak_Error_PPM_List",
    "Positive_Peak_Within_Formal_Tolerance_Count", "Formal_Tolerance_Ambiguous",
    "Formal_Tolerance_Positive_Peak_MZ_List", "Formal_Tolerance_Positive_Peak_Error_PPM_List",
    "Formal_Tolerance_Positive_Peak_Intensity_List", "Formal_Matched_Physical_Peak_Count",
    "Formal_Matched_Theoretical_Ion_Count", "Formal_Match_Ambiguous",
    "Formal_Matched_Peak_ID_List", "Formal_Matched_Theoretical_Ion_ID_List", "Formal_Match_Sharing_Type",
    "Candidate_Specific_Raw_Peak_Count", "Candidate_Specific_Positive_Peak_Count", "Candidate_Specific_Formal_Peak_Count",
    "Position_Group_Shared_Raw_Peak_Count", "Position_Group_Shared_Positive_Peak_Count", "Position_Group_Shared_Formal_Peak_Count",
    "Structural_Isomer_Shared_Raw_Peak_Count", "Structural_Isomer_Shared_Positive_Peak_Count", "Structural_Isomer_Shared_Formal_Peak_Count",
    "Cross_Candidate_Shared_Raw_Peak_Count", "Cross_Candidate_Shared_Positive_Peak_Count", "Cross_Candidate_Shared_Formal_Peak_Count",
    "Effective_Ambiguity_Level", "Effective_Ambiguity_Type", "Effective_Ambiguity_Severity",
    "Effective_Ambiguity_Reason", "Effective_Ambiguity_Recommendation",
    "Effective_Ambiguity_Applied_To_Final_Score",
]

DETAIL_COLUMNS = [
    "Ambiguity_Cluster_ID", "Candidate_Key", "Theoretical_Ion_ID", "Physical_Peak_ID",
    "Spectrum_ID", "Scan_Index", "Peak_Index", "MZ", "Intensity", "Intensity_State",
    "Error_Da", "Error_PPM", "Within_Audit_Window", "Within_Formal_Tolerance",
    "Is_Raw_Peak", "Is_Positive_Peak", "Is_Formal_Matched_Peak", "Used_In_Identity",
    "Used_In_Localization", "Candidate_Specific", "Position_Group_Shared",
    "Structural_Isomer_Shared", "Cross_Candidate_Shared", "Raw_Ambiguous",
    "Positive_Ambiguous", "Formal_Tolerance_Ambiguous", "Formal_Match_Ambiguous",
    "Effective_Ambiguity_Level",
]

SUMMARY_COLUMNS = [
    "Total_Clusters", "Raw_Ambiguous_Clusters", "Positive_Ambiguous_Clusters",
    "Formal_Tolerance_Ambiguous_Clusters", "Formal_Match_Ambiguous_Clusters",
    "Raw_Only_Clusters", "Effective_None_Clusters", "Affected_Candidate_Count",
    "Raw_Only_Candidate_Count", "Positive_Ambiguity_Candidate_Count",
    "Formal_Tolerance_Ambiguity_Candidate_Count", "Formal_Match_Ambiguity_Candidate_Count",
    "Top50_Affected_Count", "cnm5U_Affected_Count", "High_Count", "Moderate_Count",
    "Low_Count", "Informational_Count", "None_Count", "Candidate_Specific_Formal_Peak_Count",
    "Position_Group_Shared_Formal_Peak_Count", "Structural_Isomer_Shared_Formal_Peak_Count",
    "Cross_Candidate_Shared_Formal_Peak_Count", "Detail_Original_Row_Count",
    "Detail_Written_Row_Count", "Detail_Truncated", "Detail_Truncation_Reason",
    "Audit_Mode", "Apply_To_Final_Score", "Effective_Ambiguity_Definition",
    "Formal_Tolerance_Definition", "Formal_Match_Ambiguity_Definition", "Positive_Intensity_Definition",
]

DIAGNOSTIC_COLUMNS = [
    "Effective_Ambiguity_Audit_Available", "Effective_Ambiguity_Total_Clusters",
    "Raw_Ambiguous_Cluster_Count", "Positive_Ambiguous_Cluster_Count",
    "Formal_Tolerance_Ambiguous_Cluster_Count", "Formal_Match_Ambiguous_Cluster_Count",
    "Raw_Only_Cluster_Count", "Effective_Ambiguity_High_Count",
    "Effective_Ambiguity_Moderate_Count", "Effective_Ambiguity_Low_Count",
    "Effective_Ambiguity_Informational_Count", "Effective_Ambiguity_Recommendation",
    "Effective_Ambiguity_Applied_To_Final_Score",
]

TOP_SHADOW_COLUMNS = [
    "Effective_Ambiguity_Affected", "Raw_Only_Ambiguity_Cluster_Count",
    "Positive_Ambiguity_Cluster_Count", "Formal_Tolerance_Ambiguity_Cluster_Count",
    "Formal_Match_Ambiguity_Cluster_Count", "Effective_Ambiguity_High_Count",
    "Effective_Ambiguity_Moderate_Count", "Effective_Ambiguity_Low_Count",
    "Effective_Ambiguity_Informational_Count", "Effective_Ambiguity_Severity",
    "Effective_Ambiguity_Recommendation", "Effective_Ambiguity_Applied_To_Final_Score",
]

CANDIDATE_COLUMNS = ["Modification_ID", "Parent_Fragment_ID", "Candidate_tRNA_Position", *TOP_SHADOW_COLUMNS]
LEVEL_ORDER = {"none": 0, "raw_only": 1, "positive_intensity": 2, "formal_tolerance": 3, "formal_match": 4}
SEVERITY_ORDER = {"none": 0, "informational": 1, "low": 2, "moderate": 3, "high": 4}


def _float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError, OverflowError):
        return None


def _bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().lower() in {"1", "true", "yes", "on"}


def _position(value: Any) -> str:
    number = _float(value)
    return str(int(number)) if number is not None and number > 0 and number.is_integer() else ""


def _candidate(mod: Any, parent: Any, position: Any) -> tuple[str, str, str]:
    return str(mod or ""), str(parent or ""), _position(position)


def _candidate_text(key: tuple[str, str, str]) -> str:
    return "|".join(value or "NA" for value in key)


def _physical_key(row: dict[str, Any]) -> str:
    return physical_observed_peak_key({
        "Spectrum_ID": row.get("Spectrum_ID"), "Observed_mz": row.get("Observed_mz"),
        "RT": row.get("RT"), "Peak_ID": row.get("Peak_ID"),
    })


def _list(values: list[Any]) -> str:
    return ";".join(str(value) for value in values)


def _family_lookup(ranking_rows: list[dict[str, Any]]) -> tuple[dict[tuple[str, str, str], str], dict[tuple[str, str, str], str]]:
    families, parent_to_trna = {}, {}
    for row in ranking_rows:
        key = _candidate(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), row.get("Candidate_tRNA_Position"))
        families[key] = str(row.get("Modification_Family") or "")
        parent_pos = _position(row.get("Candidate_Position_In_Parent"))
        if parent_pos:
            parent_to_trna[(key[0], key[1], parent_pos)] = key[2]
    return families, parent_to_trna


def _formal_assignments(
    ion_matches: list[dict[str, Any]], modified_matches: list[dict[str, Any]],
    identity_assignments: list[dict[str, Any]], parent_to_trna: dict[tuple[str, str, str], str],
) -> dict[str, dict[str, Any]]:
    assignments: dict[str, dict[str, Any]] = defaultdict(lambda: {
        "theoretical": set(), "candidates": set(), "mods": set(), "positions": set(),
        "structural_groups": set(), "identity": False, "localization_candidates": set(),
    })
    for row in ion_matches:
        key = _physical_key(row)
        assignments[key]["theoretical"].add(str(row.get("Best_Ion_ID") or ""))
    for row in modified_matches:
        peak = _physical_key(row)
        parent_pos = _position(row.get("Candidate_Modification_Position_In_Parent"))
        trna = parent_to_trna.get((str(row.get("Modification_ID") or ""), str(row.get("Parent_Fragment_ID") or ""), parent_pos), parent_pos)
        candidate = _candidate(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), trna)
        item = assignments[peak]
        item["theoretical"].add(str(row.get("Ion_ID") or ""))
        item["candidates"].add(candidate); item["mods"].add(candidate[0]); item["positions"].add(candidate[2])
    for row in identity_assignments:
        peak = str(row.get("Physical_Observed_Peak_Key") or "")
        if not peak:
            continue
        candidate = _candidate(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), row.get("Candidate_tRNA_Position"))
        item = assignments[peak]
        item["identity"] = True; item["candidates"].add(candidate)
        item["mods"].add(candidate[0]); item["positions"].add(candidate[2])
        theory = str(row.get("Theoretical_Ion_ID") or "")
        if theory: item["theoretical"].add(theory)
        group = str(row.get("Structural_Isomer_Group_ID") or "")
        if group: item["structural_groups"].add(group)
    return assignments


def _peak_locations(context: dict[str, Any], wanted: dict[str, set[float]]) -> dict[tuple[str, float], tuple[Any, int]]:
    locations = {}
    for record in context.get("source_spectra", []) or []:
        spectrum_id = str(record.get("spectrum_id") or "")
        wanted_mz = wanted.get(spectrum_id, set())
        if not wanted_mz:
            continue
        mz_values = record.get("original_mz", [])
        if mz_values is None:
            mz_values = []
        for index, mz in enumerate(mz_values):
            rounded = round(float(mz), 8)
            if rounded in wanted_mz:
                locations[(spectrum_id, rounded)] = (record.get("scan_index", ""), index)
    return locations


def _sharing_type(items: list[dict[str, Any]], cluster_theory: str) -> tuple[bool, str, list[str]]:
    types = []
    theory_to_peaks: dict[str, set[str]] = defaultdict(set)
    for item in items:
        peak = item["peak"]
        assignment = item["assignment"]
        for theory in assignment.get("theoretical", set()): theory_to_peaks[theory].add(peak)
        if len(assignment.get("theoretical", set())) > 1: types.append("multiple_theoretical_ions_share_formal_peak")
        candidates = assignment.get("candidates", set())
        if len(candidates) > 1:
            types.append("multiple_candidates_share_formal_peak")
            mods = {key[0] for key in candidates}; positions = {key[2] for key in candidates}
            if len(mods) == 1 and len(positions) > 1: types.append("position_isomers_share_formal_peak")
            if assignment.get("structural_groups"): types.append("structural_isomers_share_formal_peak")
    if len(theory_to_peaks.get(cluster_theory, set())) > 1:
        types.append("multiple_formal_peaks_for_same_theoretical_ion")
    unique = list(dict.fromkeys(types))
    return bool(unique), ";".join(unique) if unique else "none", unique


def _effective_type(level: str, positive: list[dict[str, Any]], formal_types: list[str], theoretical: float) -> str:
    if level == "formal_match":
        return formal_types[0] if len(formal_types) == 1 else "mixed_effective_ambiguity"
    if level == "formal_tolerance": return "multiple_positive_peaks_within_formal_tolerance"
    if level == "positive_intensity":
        below = any(float(row["Observed_mz"]) < theoretical for row in positive)
        above = any(float(row["Observed_mz"]) > theoretical for row in positive)
        return "multiple_positive_peaks_bracketing" if below and above else "multiple_positive_peaks_same_side"
    if level == "raw_only": return "raw_zero_inflated_multiple_peaks"
    return "none"


def _classify(raw: bool, positive: bool, tolerance: bool, formal: bool) -> tuple[str, str, str, str]:
    if formal: return "formal_match", "high", "formal assignment competition exists", "inspect_formal_assignment_competition"
    if tolerance: return "formal_tolerance", "moderate", "multiple positive peaks are within formal tolerance", "inspect_formal_tolerance_competition"
    if positive: return "positive_intensity", "low", "multiple positive raw peaks remain outside formal competition", "inspect_positive_peaks"
    if raw: return "raw_only", "informational", "raw multiplicity is inflated by non-positive peaks", "retain_raw_ambiguity_for_provenance_only"
    return "none", "none", "no effective ambiguity", "no_effective_ambiguity"


def build_effective_ambiguity(
    ambiguity_clusters: list[dict[str, Any]] | None, ambiguity_details: list[dict[str, Any]] | None,
    ambiguity_summary: list[dict[str, Any]] | None, ion_matches: list[dict[str, Any]] | None,
    modified_matches: list[dict[str, Any]] | None, identity_assignments: list[dict[str, Any]] | None,
    localization_rows: list[dict[str, Any]] | None, ranking_rows: list[dict[str, Any]] | None,
    zero_context: dict[str, Any] | None = None, *, enabled: bool = True, max_detail_rows: int = 100000,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    clusters = list(ambiguity_clusters or []); details = list(ambiguity_details or []); ranking = list(ranking_rows or [])
    if not enabled:
        summary = _summary_rows([], [], 0, 0, False, "")[0]
        return [], [], [summary], [], [_diagnostics([summary], False)]
    families, parent_to_trna = _family_lookup(ranking)
    formal = _formal_assignments(list(ion_matches or []), list(modified_matches or []), list(identity_assignments or []), parent_to_trna)
    localized = {_candidate(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), row.get("Candidate_Modification_Position_In_tRNA")) for row in localization_rows or [] if int(_float(row.get("Num_Modified_Ion_Matches")) or 0) > 0}
    for peak, item in formal.items(): item["localization_candidates"] = item["candidates"] & localized
    by_cluster: dict[str, list[dict[str, Any]]] = defaultdict(list)
    wanted: dict[str, set[float]] = defaultdict(set)
    for row in details:
        by_cluster[str(row.get("Peak_Cluster_ID") or "")].append(row)
        mz = _float(row.get("Observed_mz"))
        if mz is not None: wanted[str(row.get("Spectrum_ID") or "")].add(round(mz, 8))
    locations = _peak_locations(zero_context or {}, wanted)
    severity_lookup = {_candidate(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), row.get("Candidate_tRNA_Position")): row.get("Ambiguity_Severity", "unknown") for row in ambiguity_summary or []}
    out_clusters, out_details = [], []
    for cluster in sorted(clusters, key=lambda row: str(row.get("Peak_Cluster_ID") or "")):
        cluster_id = str(cluster.get("Peak_Cluster_ID") or "")
        candidate = _candidate(cluster.get("Modification_ID"), cluster.get("Parent_Fragment_ID"), cluster.get("Candidate_tRNA_Position"))
        rows = sorted(by_cluster.get(cluster_id, []), key=lambda row: (float(row.get("Observed_mz") or 0), str(row.get("Physical_Observed_Peak_Key") or "")))
        positive = [row for row in rows if float(row.get("Intensity") or 0) > 0]
        tolerance_positive = [row for row in positive if _bool(row.get("Within_Formal_Tolerance"))]
        formal_items = []
        for row in rows:
            peak = str(row.get("Physical_Observed_Peak_Key") or "")
            if peak in formal: formal_items.append({"peak": peak, "row": row, "assignment": formal[peak]})
        cluster_theory = str(cluster.get("Theoretical_Ion_ID") or "")
        formal_ambiguous, sharing, formal_types = _sharing_type(formal_items, cluster_theory)
        raw_ambiguous = len(rows) >= 2; positive_ambiguous = len(positive) >= 2; tolerance_ambiguous = len(tolerance_positive) >= 2
        level, severity, reason, recommendation = _classify(raw_ambiguous, positive_ambiguous, tolerance_ambiguous, formal_ambiguous)
        if "position_isomers_share_formal_peak" in formal_types: recommendation = "inspect_position_isomer_evidence"
        elif "structural_isomers_share_formal_peak" in formal_types: recommendation = "inspect_structural_isomer_evidence"
        elif level == "raw_only" and any(str(row.get("Candidate_Specificity_Status") or "") == "position_group_shared" for row in rows):
            recommendation = "positional_isomer_remains_unresolved"
            reason += "; no positive position-discriminating formal evidence remains"
        effective_type = _effective_type(level, positive, formal_types, float(cluster.get("Theoretical_mz") or 0))
        formal_peaks = sorted({item["peak"] for item in formal_items})
        formal_theories = sorted({theory for item in formal_items for theory in item["assignment"].get("theoretical", set()) if theory})
        def count(status: str, subset: list[dict[str, Any]], formal_only: bool = False) -> int:
            return sum(str(row.get("Candidate_Specificity_Status") or "") == status and (not formal_only or str(row.get("Physical_Observed_Peak_Key") or "") in formal) for row in subset)
        record = {
            "Ambiguity_Cluster_ID": cluster_id, "Candidate_Key": _candidate_text(candidate), "Modification": candidate[0],
            "Modification_Family": families.get(candidate, ""), "Position": candidate[2], "Theoretical_Ion_ID": cluster_theory,
            "Ion_Series": cluster.get("Ion_Series"), "Ion_Number": cluster.get("Ion_Number"), "Theoretical_MZ": cluster.get("Theoretical_mz"),
            "Existing_Primary_Ambiguity_Type": cluster.get("Primary_Ambiguity_Type"), "Existing_Severity": severity_lookup.get(candidate, "unknown"),
            "Raw_Peak_Count": len(rows), "Raw_Zero_Intensity_Peak_Count": sum(float(row.get("Intensity") or 0) == 0 for row in rows),
            "Raw_Positive_Intensity_Peak_Count": len(positive), "Raw_Ambiguous": raw_ambiguous,
            "Positive_Peak_Count": len(positive), "Positive_Ambiguous": positive_ambiguous,
            "Positive_Peak_MZ_List": _list([row.get("Observed_mz") for row in positive]),
            "Positive_Peak_Intensity_List": _list([row.get("Intensity") for row in positive]),
            "Positive_Peak_Error_PPM_List": _list([row.get("Error_ppm") for row in positive]),
            "Positive_Peak_Within_Formal_Tolerance_Count": len(tolerance_positive), "Formal_Tolerance_Ambiguous": tolerance_ambiguous,
            "Formal_Tolerance_Positive_Peak_MZ_List": _list([row.get("Observed_mz") for row in tolerance_positive]),
            "Formal_Tolerance_Positive_Peak_Error_PPM_List": _list([row.get("Error_ppm") for row in tolerance_positive]),
            "Formal_Tolerance_Positive_Peak_Intensity_List": _list([row.get("Intensity") for row in tolerance_positive]),
            "Formal_Matched_Physical_Peak_Count": len(formal_peaks), "Formal_Matched_Theoretical_Ion_Count": len(formal_theories),
            "Formal_Match_Ambiguous": formal_ambiguous, "Formal_Matched_Peak_ID_List": ";".join(formal_peaks),
            "Formal_Matched_Theoretical_Ion_ID_List": ";".join(formal_theories), "Formal_Match_Sharing_Type": sharing,
            "Candidate_Specific_Raw_Peak_Count": count("candidate_specific", rows), "Candidate_Specific_Positive_Peak_Count": count("candidate_specific", positive), "Candidate_Specific_Formal_Peak_Count": count("candidate_specific", rows, True),
            "Position_Group_Shared_Raw_Peak_Count": count("position_group_shared", rows), "Position_Group_Shared_Positive_Peak_Count": count("position_group_shared", positive), "Position_Group_Shared_Formal_Peak_Count": count("position_group_shared", rows, True),
            "Structural_Isomer_Shared_Raw_Peak_Count": count("structural_isomer_group_shared", rows), "Structural_Isomer_Shared_Positive_Peak_Count": count("structural_isomer_group_shared", positive), "Structural_Isomer_Shared_Formal_Peak_Count": count("structural_isomer_group_shared", rows, True),
            "Cross_Candidate_Shared_Raw_Peak_Count": count("cross_candidate_shared", rows), "Cross_Candidate_Shared_Positive_Peak_Count": count("cross_candidate_shared", positive), "Cross_Candidate_Shared_Formal_Peak_Count": count("cross_candidate_shared", rows, True),
            "Effective_Ambiguity_Level": level, "Effective_Ambiguity_Type": effective_type,
            "Effective_Ambiguity_Severity": severity, "Effective_Ambiguity_Reason": reason,
            "Effective_Ambiguity_Recommendation": recommendation, "Effective_Ambiguity_Applied_To_Final_Score": False,
        }
        out_clusters.append(record)
        for row in rows:
            peak = str(row.get("Physical_Observed_Peak_Key") or ""); assignment = formal.get(peak, {})
            spectrum_id = str(row.get("Spectrum_ID") or ""); mz = float(row.get("Observed_mz") or 0)
            scan, peak_index = locations.get((spectrum_id, round(mz, 8)), ("", ""))
            status = str(row.get("Candidate_Specificity_Status") or "")
            out_details.append({
                "Ambiguity_Cluster_ID": cluster_id, "Candidate_Key": _candidate_text(candidate), "Theoretical_Ion_ID": cluster_theory,
                "Physical_Peak_ID": peak, "Spectrum_ID": spectrum_id, "Scan_Index": scan, "Peak_Index": peak_index,
                "MZ": row.get("Observed_mz"), "Intensity": row.get("Intensity"), "Intensity_State": "positive" if float(row.get("Intensity") or 0) > 0 else "zero",
                "Error_Da": row.get("Error_Da"), "Error_PPM": row.get("Error_ppm"), "Within_Audit_Window": row.get("Within_Audit_Window", True),
                "Within_Formal_Tolerance": row.get("Within_Formal_Tolerance"), "Is_Raw_Peak": True,
                "Is_Positive_Peak": float(row.get("Intensity") or 0) > 0, "Is_Formal_Matched_Peak": peak in formal,
                "Used_In_Identity": bool(assignment.get("identity")), "Used_In_Localization": candidate in assignment.get("localization_candidates", set()),
                "Candidate_Specific": status == "candidate_specific", "Position_Group_Shared": status == "position_group_shared",
                "Structural_Isomer_Shared": status == "structural_isomer_group_shared", "Cross_Candidate_Shared": status == "cross_candidate_shared",
                "Raw_Ambiguous": raw_ambiguous, "Positive_Ambiguous": positive_ambiguous,
                "Formal_Tolerance_Ambiguous": tolerance_ambiguous, "Formal_Match_Ambiguous": formal_ambiguous,
                "Effective_Ambiguity_Level": level,
            })
    original = len(out_details); limit = max(1, int(max_detail_rows or 100000)); written = out_details[:limit]
    truncated = original > len(written); trunc_reason = f"deterministic cluster/peak ordering truncated at {limit} configured Excel rows" if truncated else ""
    candidate_rows = build_candidate_summary(ranking, out_clusters)
    summary = _summary_rows(out_clusters, candidate_rows, original, len(written), truncated, trunc_reason)
    return out_clusters, written, summary, candidate_rows, [_diagnostics(summary, True)]


def build_candidate_summary(ranking_rows: list[dict[str, Any]], clusters: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in clusters:
        parts = str(row.get("Candidate_Key") or "").split("|")
        if len(parts) == 3: grouped[(parts[0] if parts[0] != "NA" else "", parts[1] if parts[1] != "NA" else "", parts[2] if parts[2] != "NA" else "")].append(row)
    keys = {_candidate(row.get("Modification_ID"), row.get("Parent_Fragment_ID"), row.get("Candidate_tRNA_Position")) for row in ranking_rows} | set(grouped)
    output = []
    for key in sorted(keys):
        rows = grouped.get(key, []); levels = Counter(row.get("Effective_Ambiguity_Level") for row in rows); severities = Counter(row.get("Effective_Ambiguity_Severity") for row in rows)
        severity = max(severities, key=lambda value: SEVERITY_ORDER.get(str(value), 0)) if severities else "none"
        recommendation = next((row.get("Effective_Ambiguity_Recommendation") for row in rows if row.get("Effective_Ambiguity_Severity") == severity), "no_effective_ambiguity")
        output.append({
            "Modification_ID": key[0], "Parent_Fragment_ID": key[1], "Candidate_tRNA_Position": key[2],
            "Effective_Ambiguity_Affected": any(level != "none" for level in levels),
            "Raw_Only_Ambiguity_Cluster_Count": levels["raw_only"], "Positive_Ambiguity_Cluster_Count": levels["positive_intensity"],
            "Formal_Tolerance_Ambiguity_Cluster_Count": levels["formal_tolerance"], "Formal_Match_Ambiguity_Cluster_Count": levels["formal_match"],
            "Effective_Ambiguity_High_Count": severities["high"], "Effective_Ambiguity_Moderate_Count": severities["moderate"],
            "Effective_Ambiguity_Low_Count": severities["low"], "Effective_Ambiguity_Informational_Count": severities["informational"],
            "Effective_Ambiguity_Severity": severity, "Effective_Ambiguity_Recommendation": recommendation,
            "Effective_Ambiguity_Applied_To_Final_Score": False,
        })
    return output


def _summary_rows(clusters: list[dict[str, Any]], candidates: list[dict[str, Any]], original: int, written: int, truncated: bool, reason: str) -> list[dict[str, Any]]:
    levels = Counter(row.get("Effective_Ambiguity_Level") for row in clusters); severity = Counter(row.get("Effective_Ambiguity_Severity") for row in clusters)
    affected = [row for row in candidates if row.get("Effective_Ambiguity_Affected")]
    return [{
        "Total_Clusters": len(clusters), "Raw_Ambiguous_Clusters": sum(_bool(row.get("Raw_Ambiguous")) for row in clusters),
        "Positive_Ambiguous_Clusters": sum(_bool(row.get("Positive_Ambiguous")) for row in clusters),
        "Formal_Tolerance_Ambiguous_Clusters": sum(_bool(row.get("Formal_Tolerance_Ambiguous")) for row in clusters),
        "Formal_Match_Ambiguous_Clusters": sum(_bool(row.get("Formal_Match_Ambiguous")) for row in clusters),
        "Raw_Only_Clusters": levels["raw_only"], "Effective_None_Clusters": levels["none"], "Affected_Candidate_Count": len(affected),
        "Raw_Only_Candidate_Count": sum(int(row.get("Raw_Only_Ambiguity_Cluster_Count") or 0) > 0 for row in candidates),
        "Positive_Ambiguity_Candidate_Count": sum(int(row.get("Positive_Ambiguity_Cluster_Count") or 0) > 0 for row in candidates),
        "Formal_Tolerance_Ambiguity_Candidate_Count": sum(int(row.get("Formal_Tolerance_Ambiguity_Cluster_Count") or 0) > 0 for row in candidates),
        "Formal_Match_Ambiguity_Candidate_Count": sum(int(row.get("Formal_Match_Ambiguity_Cluster_Count") or 0) > 0 for row in candidates),
        "Top50_Affected_Count": 0, "cnm5U_Affected_Count": sum(row.get("Modification_ID") == "cnm5U" for row in affected),
        "High_Count": severity["high"], "Moderate_Count": severity["moderate"], "Low_Count": severity["low"],
        "Informational_Count": severity["informational"], "None_Count": severity["none"],
        "Candidate_Specific_Formal_Peak_Count": sum(int(row.get("Candidate_Specific_Formal_Peak_Count") or 0) for row in clusters),
        "Position_Group_Shared_Formal_Peak_Count": sum(int(row.get("Position_Group_Shared_Formal_Peak_Count") or 0) for row in clusters),
        "Structural_Isomer_Shared_Formal_Peak_Count": sum(int(row.get("Structural_Isomer_Shared_Formal_Peak_Count") or 0) for row in clusters),
        "Cross_Candidate_Shared_Formal_Peak_Count": sum(int(row.get("Cross_Candidate_Shared_Formal_Peak_Count") or 0) for row in clusters),
        "Detail_Original_Row_Count": original, "Detail_Written_Row_Count": written, "Detail_Truncated": truncated, "Detail_Truncation_Reason": reason,
        "Audit_Mode": "shadow_four_stage_effective_ambiguity", "Apply_To_Final_Score": False,
        "Effective_Ambiguity_Definition": "strongest of formal_match > formal_tolerance > positive_intensity > raw_only > none",
        "Formal_Tolerance_Definition": "at least two positive physical peaks within existing per-ion formal tolerance",
        "Formal_Match_Ambiguity_Definition": "competition present in existing formal assignment rows, not raw-window multiplicity alone",
        "Positive_Intensity_Definition": "raw peak intensity > 0 without filtering or rematching",
    }]


def _diagnostics(summary_rows: list[dict[str, Any]], available: bool) -> dict[str, Any]:
    row = summary_rows[0] if summary_rows else {}
    if int(row.get("Formal_Match_Ambiguous_Clusters") or 0): recommendation = "inspect_formal_assignment_competition"
    elif int(row.get("Formal_Tolerance_Ambiguous_Clusters") or 0): recommendation = "inspect_formal_tolerance_competition"
    elif int(row.get("Positive_Ambiguous_Clusters") or 0): recommendation = "inspect_positive_peaks"
    elif int(row.get("Raw_Only_Clusters") or 0): recommendation = "retain_raw_ambiguity_for_provenance_only"
    else: recommendation = "no_effective_ambiguity"
    return {
        "Effective_Ambiguity_Audit_Available": available, "Effective_Ambiguity_Total_Clusters": row.get("Total_Clusters", 0),
        "Raw_Ambiguous_Cluster_Count": row.get("Raw_Ambiguous_Clusters", 0), "Positive_Ambiguous_Cluster_Count": row.get("Positive_Ambiguous_Clusters", 0),
        "Formal_Tolerance_Ambiguous_Cluster_Count": row.get("Formal_Tolerance_Ambiguous_Clusters", 0), "Formal_Match_Ambiguous_Cluster_Count": row.get("Formal_Match_Ambiguous_Clusters", 0),
        "Raw_Only_Cluster_Count": row.get("Raw_Only_Clusters", 0), "Effective_Ambiguity_High_Count": row.get("High_Count", 0),
        "Effective_Ambiguity_Moderate_Count": row.get("Moderate_Count", 0), "Effective_Ambiguity_Low_Count": row.get("Low_Count", 0),
        "Effective_Ambiguity_Informational_Count": row.get("Informational_Count", 0), "Effective_Ambiguity_Recommendation": recommendation,
        "Effective_Ambiguity_Applied_To_Final_Score": False,
    }


def update_top50_affected(summary_rows: list[dict[str, Any]], top_rows: Any) -> None:
    if not summary_rows: return
    try: records = top_rows.to_dict("records")
    except AttributeError: records = list(top_rows or [])
    summary_rows[0]["Top50_Affected_Count"] = sum(_bool(row.get("Effective_Ambiguity_Affected")) for row in records)
