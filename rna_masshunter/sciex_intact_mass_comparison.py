"""Terminal-aware, mass-only shadow comparison for SCIEX intact peaks."""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from types import MappingProxyType
from typing import Any, Iterable, Mapping

from rna_masshunter.intact_rna_mass import (
    FivePrimeState,
    IntactRnaMassParameters,
    RnaTopology,
    ThreePrimeState,
    calculate_intact_rna_mass,
)

ALGORITHM_VERSION = "sciex-intact-mass-comparison-v2"
OBSERVED_MASS_TYPE = "UNKNOWN"
MASS_DEFINITION_COMPATIBILITY = "UNKNOWN"
MASS_TYPE_ASSUMPTION = (
    "SCIEX reconstructed Mass type was not supplied; no monoisotopic, average, "
    "or most-abundant-isotope identity is assumed."
)
FORMAL_FALSE = {
    "Shadow_Only": True,
    "Mass_Match_Only": True,
    "Applied_To_Formal_Score": False,
    "Applied_To_Ranking": False,
    "Applied_To_Candidate_Filtering": False,
    "Applied_To_Final_Consensus": False,
    "SCIEX_Intact_Mass_Matching_Applied_To_Formal_Score": False,
    "SCIEX_Intact_Mass_Matching_Applied_To_Ranking": False,
    "SCIEX_Intact_Mass_Matching_Applied_To_Candidate_Filtering": False,
    "SCIEX_Intact_Mass_Matching_Applied_To_Final_Consensus": False,
    "Structure_Identity_Assigned": False,
    "Position_Assigned": False,
    "Molecular_Identity_Assigned": False,
    "Modification_Assigned": False,
    "Modification_Lookup_Performed": False,
}

# Existing columns remain in their original order. New fields are appended so the
# current Excel serializer and downstream readers retain their established keys.
DETAIL_COLUMNS = [
    "Comparison_ID", "Peak_ID", "Source_File", "Observed_Mass",
    "Theoretical_Unmodified_Mass", "Delta_Mass", "Absolute_Delta_Mass", "Delta_ppm",
    "Nearest_Existing_Intact_Mass", "Nearest_Existing_Intact_Delta",
    "Nearest_Existing_Intact_Absolute_Delta", "Nearest_Existing_Intact_ppm",
    "Existing_Intact_Comparison_Status", "Comparison_Status", "Comparison_Eligible",
    "Strict_Tolerance_Da", "Broad_Tolerance_Da", "Apex_Intensity_Raw", "Detection_Tier",
    "Prominence", "Half_Prominence_Width_Da", "Centroid_Mass", "Centroid_Minus_Apex_Da",
    "Left_Boundary_Mass", "Right_Boundary_Mass", "Peak_Area_Raw",
    "Peak_Area_Baseline_Corrected", "Possible_Shoulder", "Broad_Peak_Flag",
    "Severe_Broad_Peak_Flag", "Edge_Peak_Flag", "Algorithm_Version", "Shadow_Only",
    "Applied_To_Formal_Score", "Applied_To_Ranking", "Applied_To_Candidate_Filtering",
    "Molecular_Identity_Assigned", "Modification_Lookup_Performed", "Notes",
    "Candidate_ID", "Candidate_Category", "Terminal_State_Candidate_ID",
    "Sequence_Source", "Sequence_Length", "Sequence_SHA256", "Ends_With_CCA", "CCA_Policy",
    "Apex_Mass", "Preferred_Mass_Field", "Observed_Mass_Raw", "Observed_Mass_Calibrated",
    "Theoretical_Mass", "Theoretical_Formula", "Five_Prime_State", "Three_Prime_State",
    "Topology", "Terminal_State_Confirmed", "Apex_Delta_Da", "Apex_Delta_ppm",
    "Centroid_Delta_Da", "Centroid_Delta_ppm", "Match_Tolerance_Class", "Match_Status",
    "Strict_Tolerance_Passed", "Exploratory_Tolerance_Passed", "Exploratory_Tolerance_Da",
    "Peak_Eligibility", "Peak_Eligibility_Reason", "Strict_Threshold_Passed",
    "Peak_Area_Complete", "Centroid_Complete", "Mass_Equivalent_Candidate_Group_ID",
    "Candidate_Ambiguity_Count", "Terminal_State_Ambiguous",
    "Mass_Equivalent_Candidate_IDs", "Identity_Gate_Status", "Identity_Gate_Passed",
    "Observed_Mass_Type", "Theoretical_Mass_Type", "Mass_Definition_Compatibility",
    "Mass_Type_Assumption", "SCIEX_Reconstruction_Settings_Available",
    "Calibration_Applied", "Calibration_Method", "Calibration_Reference",
    "Calibration_Offset_Da", "Calibration_Reference_Externally_Assigned",
    "Interpretation_Warnings", "Mass_Match_Only", "Structure_Identity_Assigned",
    "Position_Assigned", "Modification_Assigned", "Applied_To_Final_Consensus",
    "SCIEX_Intact_Mass_Matching_Applied_To_Formal_Score",
    "SCIEX_Intact_Mass_Matching_Applied_To_Ranking",
    "SCIEX_Intact_Mass_Matching_Applied_To_Candidate_Filtering",
    "SCIEX_Intact_Mass_Matching_Applied_To_Final_Consensus",
]
SUMMARY_COLUMNS = [
    "Source_File", "Parser_Status", "Detector_Status", "Comparison_Status",
    "Comparison_Eligible", "Detected_Peak_Count", "Compared_Peak_Count",
    "Strict_Match_Count", "Broad_Match_Count", "No_Match_Count",
    "Theoretical_Unmodified_Mass", "Strict_Tolerance_Da", "Broad_Tolerance_Da",
    "Closest_Peak_ID", "Closest_Observed_Mass", "Closest_Delta_Mass",
    "Closest_Absolute_Delta_Mass", "Closest_Delta_ppm", "Strongest_Peak_ID",
    "Strongest_Observed_Mass", "Strongest_Apex_Intensity_Raw",
    "Existing_Intact_Result_Available", "Existing_Intact_Mass_Count",
    "Input_Identity_Audit_Status", "Input_Identity_Conflict",
    "Input_Identity_Warning_Code", "Biological_Interpretation_Eligible",
    "Algorithm_Version", "Shadow_Only", "Applied_To_Formal_Score", "Applied_To_Ranking",
    "Applied_To_Candidate_Filtering", "Molecular_Identity_Assigned",
    "Modification_Lookup_Performed", "Notes",
    "Sequence_Source", "Sequence_Length", "Sequence_SHA256", "Ends_With_CCA", "CCA_Policy",
    "Identity_Gate_Status", "Identity_Gate_Passed", "Identity_Warnings",
    "Observed_Mass_Type", "Theoretical_Mass_Type", "Mass_Definition_Compatibility",
    "Mass_Type_Assumption", "SCIEX_Reconstruction_Settings_Available",
    "Calibration_Applied", "Calibration_Method", "Calibration_Reference",
    "Calibration_Offset_Da", "Calibration_Reference_Externally_Assigned",
    "Terminal_Candidate_Count", "Mass_Equivalent_Group_Count", "Eligible_Peak_Count",
    "Primary_Peak_Count", "Exploratory_Peak_Count", "Shape_Warning_Peak_Count",
    "Nonprimary_Peak_Count", "Exploratory_Match_Count", "Ambiguous_Match_Count",
    "No_Match_Peak_Count", "Exploratory_Tolerance_Da", "Mass_Match_Only",
    "Structure_Identity_Assigned", "Position_Assigned", "Modification_Assigned",
    "Applied_To_Final_Consensus", "SCIEX_Intact_Mass_Matching_Applied_To_Formal_Score",
    "SCIEX_Intact_Mass_Matching_Applied_To_Ranking",
    "SCIEX_Intact_Mass_Matching_Applied_To_Candidate_Filtering",
    "SCIEX_Intact_Mass_Matching_Applied_To_Final_Consensus",
]
PEAK_METADATA_COLUMNS = (
    "Apex_Intensity_Raw", "Detection_Tier", "Prominence", "Half_Prominence_Width_Da",
    "Centroid_Mass", "Centroid_Minus_Apex_Da", "Left_Boundary_Mass",
    "Right_Boundary_Mass", "Peak_Area_Raw", "Peak_Area_Baseline_Corrected",
    "Possible_Shoulder", "Broad_Peak_Flag", "Severe_Broad_Peak_Flag", "Edge_Peak_Flag",
)


@dataclass(frozen=True)
class SciexIntactMassComparisonResult:
    detail_rows: tuple[Mapping[str, Any], ...]
    summary_rows: tuple[Mapping[str, Any], ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "detail_rows", tuple(MappingProxyType(dict(row)) for row in self.detail_rows),
        )
        object.__setattr__(
            self, "summary_rows", tuple(MappingProxyType(dict(row)) for row in self.summary_rows),
        )

    def details(self) -> list[dict[str, Any]]:
        return [dict(row) for row in self.detail_rows]

    def summaries(self) -> list[dict[str, Any]]:
        return [dict(row) for row in self.summary_rows]


def _finite_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if isfinite(result) else None


def _existing_mass(value: Any) -> float | None:
    if isinstance(value, Mapping):
        for key in ("observed_mass", "Observed_Mass", "Reconstructed_Mass"):
            if key in value:
                return _finite_float(value.get(key))
        return None
    observed = getattr(value, "observed_mass", None)
    return _finite_float(observed if observed is not None else value)


def _peak_rows(detection_result: Any) -> list[dict[str, Any]]:
    values = (
        detection_result.peak_rows()
        if hasattr(detection_result, "peak_rows")
        else getattr(detection_result, "peaks", ())
    )
    rows = []
    for value in values or ():
        if isinstance(value, Mapping):
            rows.append(dict(value))
        elif hasattr(value, "to_dict"):
            rows.append(dict(value.to_dict()))
    return rows


def _diagnostics(detection_result: Any) -> dict[str, Any]:
    value = (
        detection_result.diagnostics_row()
        if hasattr(detection_result, "diagnostics_row")
        else getattr(detection_result, "diagnostics", {})
    )
    if isinstance(value, Mapping):
        return dict(value)
    if hasattr(value, "to_dict"):
        return dict(value.to_dict())
    return {}


def _intensity(row: Mapping[str, Any]) -> float:
    return _finite_float(row.get("Apex_Intensity_Raw")) or 0.0


def _identity_mapping(value: Any) -> dict[str, Any]:
    if hasattr(value, "row"):
        value = value.row()
    return dict(value) if isinstance(value, Mapping) else {}


def _identity_gate(
    identity: Mapping[str, Any],
    explicit_status: str | None,
    sequence_source: str,
) -> tuple[str, bool, list[str]]:
    warnings: list[str] = []
    if explicit_status is not None:
        status = str(explicit_status).strip().upper()
        if status not in {"CONFIRMED", "AMBIGUOUS", "CONFLICT", "UNKNOWN"}:
            raise ValueError("identity_status must be CONFIRMED, AMBIGUOUS, CONFLICT, or UNKNOWN")
    else:
        audit_status = str(identity.get("Audit_Status") or "").upper()
        if bool(identity.get("Identity_Conflict")) or audit_status == "CONFLICT":
            status = "CONFLICT"
        elif audit_status in {"MATCH", "PARTIAL_MATCH"}:
            status = "AMBIGUOUS"
            warnings.append("FILENAME_IDENTITY_DOES_NOT_CONFIRM_SEQUENCE_IDENTITY")
        else:
            status = "UNKNOWN"
    if status == "CONFIRMED" and not str(sequence_source or "").strip():
        status = "AMBIGUOUS"
        warnings.append("SEQUENCE_SOURCE_NOT_PROVIDED")
    if status == "AMBIGUOUS":
        warnings.append("IDENTITY_ASSUMPTION_DEPENDENT")
    elif status == "UNKNOWN":
        warnings.append("IDENTITY_UNKNOWN")
    elif status == "CONFLICT":
        warnings.append("IDENTITY_CONFLICT_CANDIDATE_GENERATION_SKIPPED")
    return status, status == "CONFIRMED", list(dict.fromkeys(warnings))


def _terminal_candidates(
    sequence: str,
    sequence_source: str,
    terminal_state_confirmed: bool,
) -> list[dict[str, Any]]:
    specs = (
        ("TERM_LINEAR_5P_3OH", FivePrimeState.MONOPHOSPHATE, ThreePrimeState.OH),
        ("TERM_LINEAR_5OH_3OH", FivePrimeState.OH, ThreePrimeState.OH),
        ("TERM_LINEAR_5OH_3P", FivePrimeState.OH, ThreePrimeState.MONOPHOSPHATE),
        ("TERM_LINEAR_5OH_3CYCLIC", FivePrimeState.OH, ThreePrimeState.CYCLIC_PHOSPHATE),
    )
    candidates: list[dict[str, Any]] = []
    for order, (candidate_id, five_prime, three_prime) in enumerate(specs):
        result = calculate_intact_rna_mass(
            sequence,
            parameters=IntactRnaMassParameters(
                five_prime_state=five_prime,
                three_prime_state=three_prime,
                topology=RnaTopology.LINEAR,
                terminal_state_confirmed=terminal_state_confirmed,
            ),
        )
        candidates.append({
            "Candidate_ID": candidate_id,
            "Candidate_Category": "UNMODIFIED_TERMINAL_STATE",
            "Terminal_State_Candidate_ID": candidate_id,
            "Sequence_Source": str(sequence_source),
            "Sequence_Length": result.sequence_length,
            "Sequence_SHA256": result.sequence_sha256,
            "Ends_With_CCA": result.ends_with_cca,
            "CCA_Policy": result.cca_policy,
            "Theoretical_Mass": result.monoisotopic_neutral_mass,
            "Theoretical_Formula": result.formula,
            "Five_Prime_State": result.five_prime_state,
            "Three_Prime_State": result.three_prime_state,
            "Topology": result.topology,
            "Terminal_State_Confirmed": result.terminal_state_confirmed,
            "Theoretical_Mass_Type": result.theoretical_mass_type,
            "_Candidate_Order": order,
            "_Warnings": result.warnings,
        })
    grouped: dict[str, list[dict[str, Any]]] = {}
    for candidate in candidates:
        grouped.setdefault(str(candidate["Theoretical_Formula"]), []).append(candidate)
    group_ids = {
        formula: f"TERM_MASS_EQ_{index:03d}"
        for index, formula in enumerate(
            sorted(grouped, key=lambda key: (grouped[key][0]["Theoretical_Mass"], key)), 1,
        )
    }
    for formula, members in grouped.items():
        member_ids = ";".join(sorted(str(member["Candidate_ID"]) for member in members))
        for member in members:
            member["Mass_Equivalent_Candidate_Group_ID"] = group_ids[formula]
            member["Candidate_Ambiguity_Count"] = len(members)
            member["Terminal_State_Ambiguous"] = len(members) > 1
            member["Mass_Equivalent_Candidate_IDs"] = member_ids
    return candidates


def _legacy_candidate(theory: float | None) -> dict[str, Any]:
    return {
        "Candidate_ID": "LEGACY_THEORETICAL_MASS",
        "Candidate_Category": "LEGACY_THEORETICAL_MASS_COMPATIBILITY",
        "Terminal_State_Candidate_ID": "",
        "Sequence_Source": "",
        "Sequence_Length": None,
        "Sequence_SHA256": "",
        "Ends_With_CCA": None,
        "CCA_Policy": "",
        "Theoretical_Mass": theory,
        "Theoretical_Formula": "",
        "Five_Prime_State": "UNSPECIFIED",
        "Three_Prime_State": "UNSPECIFIED",
        "Topology": "UNSPECIFIED",
        "Terminal_State_Confirmed": False,
        "Theoretical_Mass_Type": "UNKNOWN_LEGACY",
        "Mass_Equivalent_Candidate_Group_ID": "",
        "Candidate_Ambiguity_Count": 1,
        "Terminal_State_Ambiguous": False,
        "Mass_Equivalent_Candidate_IDs": "LEGACY_THEORETICAL_MASS",
        "_Candidate_Order": 0,
        "_Warnings": ("LEGACY_THEORETICAL_MASS_COMPATIBILITY_MODE",),
    }


def _peak_eligibility(peak: Mapping[str, Any], centroid: float | None) -> tuple[str, str]:
    strict = bool(peak.get("Strict_Threshold_Passed", str(peak.get("Detection_Tier") or "").upper() == "STRICT"))
    edge = bool(peak.get("Edge_Peak_Flag", False))
    area_complete = bool(peak.get("Peak_Area_Complete", True))
    centroid_complete = bool(peak.get("Centroid_Complete", centroid is not None))
    broad = bool(peak.get("Broad_Peak_Flag", False))
    shoulder = bool(peak.get("Possible_Shoulder", False))
    incomplete_reasons = []
    if edge:
        incomplete_reasons.append("edge_peak")
    if not area_complete:
        incomplete_reasons.append("peak_area_incomplete")
    if not centroid_complete:
        incomplete_reasons.append("centroid_incomplete")
    if incomplete_reasons:
        return "NOT_PRIMARY_ELIGIBLE", ";".join(incomplete_reasons)
    shape_reasons = []
    if broad:
        shape_reasons.append("broad_peak")
    if shoulder:
        shape_reasons.append("possible_shoulder")
    if shape_reasons:
        return "MATCHABLE_WITH_SHAPE_WARNING", ";".join(shape_reasons)
    if strict:
        return "PRIMARY", "strict_nonedge_complete"
    return "EXPLORATORY", "sensitive_only"


def _joined_warnings(values: Iterable[str]) -> str:
    return ";".join(dict.fromkeys(str(value) for value in values if str(value)))


def compare_sciex_intact_masses(
    detection_result: Any,
    theoretical_unmodified_mass: float | None,
    existing_intact_results: Iterable[Any] | None = None,
    *,
    source_file: str = "",
    strict_tolerance_da: float = 1.0,
    broad_tolerance_da: float = 5.0,
    exploratory_tolerance_da: float | None = None,
    input_identity_audit: Any = None,
    sequence: str | None = None,
    sequence_source: str = "",
    identity_status: str | None = None,
    terminal_state_confirmed: bool = False,
) -> SciexIntactMassComparisonResult:
    """Compare apex masses to explicit terminal candidates without chemical assignment."""
    strict = float(strict_tolerance_da)
    exploratory = float(
        broad_tolerance_da if exploratory_tolerance_da is None else exploratory_tolerance_da
    )
    if not isfinite(strict) or strict <= 0:
        raise ValueError("strict_tolerance_da must be finite and positive")
    if not isfinite(exploratory) or exploratory < strict:
        raise ValueError("broad/exploratory tolerance must be finite and >= strict_tolerance_da")
    if not isinstance(terminal_state_confirmed, bool):
        raise ValueError("terminal_state_confirmed must be boolean")

    legacy_theory = _finite_float(theoretical_unmodified_mass)
    if legacy_theory is not None and legacy_theory <= 0:
        legacy_theory = None
    peaks = _peak_rows(detection_result)
    diagnostics = _diagnostics(detection_result)
    detector_completed = diagnostics.get("Detection_Status") == "DETECTION_COMPLETED"
    mz_profile = diagnostics.get("Profile_Type") == "MZ_PROFILE"
    identity = _identity_mapping(input_identity_audit)

    explicit_sequence = sequence is not None
    resolved_sequence = sequence
    resolved_source = str(sequence_source or "")
    if explicit_sequence and not resolved_source:
        resolved_source = "CALLER_PROVIDED_SEQUENCE"
    if resolved_sequence is None and str(identity.get("Configured_Sequence") or "").strip():
        resolved_sequence = str(identity["Configured_Sequence"])
        if not resolved_source:
            resolved_source = "IDENTITY_AUDIT_CONFIGURED_SEQUENCE"

    gate_status, gate_passed, identity_warnings = _identity_gate(
        identity, identity_status, resolved_source,
    )
    candidate_mode = "TERMINAL_AWARE" if resolved_sequence is not None else "LEGACY_COMPATIBILITY"
    if gate_status == "CONFLICT" or mz_profile:
        candidates: list[dict[str, Any]] = []
    elif resolved_sequence is not None:
        candidates = _terminal_candidates(
            resolved_sequence, resolved_source, terminal_state_confirmed,
        )
    else:
        candidates = [_legacy_candidate(legacy_theory)]

    existing_masses = sorted({
        mass for item in (existing_intact_results or ())
        if (mass := _existing_mass(item)) is not None
    })
    indexed_peaks = sorted(
        (
            (index, peak, mass)
            for index, peak in enumerate(peaks)
            if (mass := _finite_float(peak.get("Apex_Mass"))) is not None
        ),
        key=lambda item: (item[2], item[0]),
    )
    peak_classes = {
        original_index: _peak_eligibility(peak, _finite_float(peak.get("Centroid_Mass")))
        for original_index, peak, _mass in indexed_peaks
    }

    details: list[dict[str, Any]] = []
    for original_index, peak, observed_mass in indexed_peaks:
        centroid = _finite_float(peak.get("Centroid_Mass"))
        peak_eligibility, eligibility_reason = peak_classes[original_index]
        strict_threshold_passed = bool(
            peak.get(
                "Strict_Threshold_Passed",
                str(peak.get("Detection_Tier") or "").upper() == "STRICT",
            )
        )
        area_complete = bool(peak.get("Peak_Area_Complete", True))
        centroid_complete = bool(peak.get("Centroid_Complete", centroid is not None))
        nearest = (
            min(existing_masses, key=lambda mass: (abs(observed_mass - mass), mass))
            if existing_masses else None
        )
        existing_delta = observed_mass - nearest if nearest is not None else None
        for candidate in candidates:
            theory = _finite_float(candidate.get("Theoretical_Mass"))
            eligible = detector_completed and not mz_profile and theory is not None
            apex_delta = observed_mass - theory if eligible else None
            absolute_delta = abs(apex_delta) if apex_delta is not None else None
            centroid_delta = centroid - theory if eligible and centroid is not None else None
            if not detector_completed or mz_profile:
                comparison_status = "NOT_ELIGIBLE"
                tolerance_class = "NOT_ELIGIBLE"
                match_status = "NOT_ELIGIBLE"
            elif theory is None:
                comparison_status = "NO_THEORETICAL_MASS"
                tolerance_class = "NO_THEORETICAL_MASS"
                match_status = "NO_THEORETICAL_MASS"
            elif absolute_delta <= strict:
                comparison_status = "STRICT_MATCH"
                tolerance_class = "STRICT"
                match_status = "MATCH_STRICT_MASS_ONLY"
            elif absolute_delta <= exploratory:
                comparison_status = "BROAD_MATCH"
                tolerance_class = "EXPLORATORY"
                match_status = "MATCH_EXPLORATORY_MASS_ONLY"
            else:
                comparison_status = "NO_MATCH"
                tolerance_class = "NO_MATCH"
                match_status = "NO_MATCH"

            warnings = [*identity_warnings, *candidate.get("_Warnings", ())]
            warnings.append("OBSERVED_MASS_TYPE_UNKNOWN")
            if candidate.get("Terminal_State_Ambiguous"):
                warnings.append("TERMINAL_STATE_AMBIGUOUS")
            if peak_eligibility == "MATCHABLE_WITH_SHAPE_WARNING":
                warnings.append("PEAK_SHAPE_WARNING")
            elif peak_eligibility == "NOT_PRIMARY_ELIGIBLE":
                warnings.append("PEAK_NOT_PRIMARY_ELIGIBLE")
            if centroid is None or not centroid_complete:
                warnings.append("CENTROID_UNAVAILABLE")

            row = {
                "Comparison_ID": "",
                "Peak_ID": peak.get("Peak_ID", ""),
                "Source_File": str(source_file),
                "Observed_Mass": observed_mass,
                "Theoretical_Unmodified_Mass": theory,
                "Delta_Mass": apex_delta,
                "Absolute_Delta_Mass": absolute_delta,
                "Delta_ppm": apex_delta / theory * 1_000_000 if apex_delta is not None else None,
                "Nearest_Existing_Intact_Mass": nearest,
                "Nearest_Existing_Intact_Delta": existing_delta,
                "Nearest_Existing_Intact_Absolute_Delta": (
                    abs(existing_delta) if existing_delta is not None else None
                ),
                "Nearest_Existing_Intact_ppm": (
                    existing_delta / nearest * 1_000_000
                    if existing_delta is not None and nearest else None
                ),
                "Existing_Intact_Comparison_Status": (
                    "AVAILABLE" if nearest is not None else "NO_EXISTING_INTACT_RESULT"
                ),
                "Comparison_Status": comparison_status,
                "Comparison_Eligible": eligible,
                "Strict_Tolerance_Da": strict,
                "Broad_Tolerance_Da": exploratory,
                **{column: peak.get(column) for column in PEAK_METADATA_COLUMNS},
                "Centroid_Mass": centroid,
                "Algorithm_Version": ALGORITHM_VERSION,
                **FORMAL_FALSE,
                "Notes": "Mass-only terminal-state comparison; no chemical identity assignment.",
                **{key: value for key, value in candidate.items() if not key.startswith("_")},
                "Apex_Mass": observed_mass,
                "Preferred_Mass_Field": "Apex_Mass",
                "Observed_Mass_Raw": observed_mass,
                "Observed_Mass_Calibrated": observed_mass,
                "Apex_Delta_Da": apex_delta,
                "Apex_Delta_ppm": apex_delta / theory * 1_000_000 if apex_delta is not None else None,
                "Centroid_Delta_Da": centroid_delta,
                "Centroid_Delta_ppm": (
                    centroid_delta / theory * 1_000_000 if centroid_delta is not None else None
                ),
                "Match_Tolerance_Class": tolerance_class,
                "Match_Status": match_status,
                "Strict_Tolerance_Passed": tolerance_class == "STRICT",
                "Exploratory_Tolerance_Passed": tolerance_class in {"STRICT", "EXPLORATORY"},
                "Exploratory_Tolerance_Da": exploratory,
                "Peak_Eligibility": peak_eligibility,
                "Peak_Eligibility_Reason": eligibility_reason,
                "Strict_Threshold_Passed": strict_threshold_passed,
                "Peak_Area_Complete": area_complete,
                "Centroid_Complete": centroid_complete,
                "Identity_Gate_Status": gate_status,
                "Identity_Gate_Passed": gate_passed,
                "Observed_Mass_Type": OBSERVED_MASS_TYPE,
                "Mass_Definition_Compatibility": MASS_DEFINITION_COMPATIBILITY,
                "Mass_Type_Assumption": MASS_TYPE_ASSUMPTION,
                "SCIEX_Reconstruction_Settings_Available": False,
                "Calibration_Applied": False,
                "Calibration_Method": "NONE",
                "Calibration_Reference": "",
                "Calibration_Offset_Da": 0.0,
                "Calibration_Reference_Externally_Assigned": False,
                "Interpretation_Warnings": _joined_warnings(warnings),
                "_Original_Index": original_index,
                "_Candidate_Order": candidate.get("_Candidate_Order", 0),
            }
            details.append(row)

    details.sort(key=lambda row: (row["Observed_Mass"], row["_Original_Index"], row["_Candidate_Order"]))
    for position, row in enumerate(details, 1):
        row["Comparison_ID"] = f"SCIEX_CMP_{position:05d}"

    eligible_details = [row for row in details if row["Absolute_Delta_Mass"] is not None]
    closest = min(
        eligible_details,
        key=lambda row: (
            row["Absolute_Delta_Mass"], -_intensity(row), row["Observed_Mass"],
            row["_Original_Index"], row["_Candidate_Order"],
        ),
    ) if eligible_details else None
    strongest = min(
        details,
        key=lambda row: (
            -_intensity(row), row["Observed_Mass"], row["_Original_Index"], row["_Candidate_Order"],
        ),
    ) if details else None

    if gate_status == "CONFLICT":
        summary_status = "IDENTITY_CONFLICT"
    elif mz_profile:
        summary_status = "NOT_ELIGIBLE"
    elif not detector_completed or not indexed_peaks:
        summary_status = "NOT_ELIGIBLE"
    elif not any(_finite_float(candidate.get("Theoretical_Mass")) for candidate in candidates):
        summary_status = "NO_THEORETICAL_MASS"
    else:
        summary_status = "COMPARISON_COMPLETED"

    matched_by_peak: dict[int, list[dict[str, Any]]] = {}
    for row in details:
        if row["Match_Tolerance_Class"] in {"STRICT", "EXPLORATORY"}:
            matched_by_peak.setdefault(row["_Original_Index"], []).append(row)
    no_match_peaks = sum(index not in matched_by_peak for index, _peak, _mass in indexed_peaks)
    ambiguous_match_peaks = sum(
        any(row["Terminal_State_Ambiguous"] for row in rows)
        for rows in matched_by_peak.values()
    )
    terminal_candidates = [
        candidate for candidate in candidates
        if candidate.get("Candidate_Category") == "UNMODIFIED_TERMINAL_STATE"
    ]
    group_ids = {
        candidate["Mass_Equivalent_Candidate_Group_ID"]
        for candidate in terminal_candidates
    }
    reference = next(
        (candidate for candidate in terminal_candidates if candidate["Candidate_ID"] == "TERM_LINEAR_5P_3OH"),
        candidates[0] if candidates else {},
    )
    raw_identity_status = str(identity.get("Audit_Status") or "NOT_RUN")
    biological_eligible = bool(identity.get("Biological_Interpretation_Eligible", False))
    summary = {
        "Source_File": str(source_file),
        "Parser_Status": str(diagnostics.get("Input_Status") or ""),
        "Detector_Status": str(diagnostics.get("Detection_Status") or ""),
        "Comparison_Status": summary_status,
        "Comparison_Eligible": summary_status == "COMPARISON_COMPLETED",
        "Detected_Peak_Count": len(peaks),
        "Compared_Peak_Count": len(indexed_peaks) if details else 0,
        "Strict_Match_Count": sum(row["Comparison_Status"] == "STRICT_MATCH" for row in details),
        "Broad_Match_Count": sum(row["Comparison_Status"] == "BROAD_MATCH" for row in details),
        "No_Match_Count": sum(row["Comparison_Status"] == "NO_MATCH" for row in details),
        "Theoretical_Unmodified_Mass": _finite_float(reference.get("Theoretical_Mass")),
        "Strict_Tolerance_Da": strict,
        "Broad_Tolerance_Da": exploratory,
        "Closest_Peak_ID": closest.get("Peak_ID") if closest else "",
        "Closest_Observed_Mass": closest.get("Observed_Mass") if closest else None,
        "Closest_Delta_Mass": closest.get("Delta_Mass") if closest else None,
        "Closest_Absolute_Delta_Mass": closest.get("Absolute_Delta_Mass") if closest else None,
        "Closest_Delta_ppm": closest.get("Delta_ppm") if closest else None,
        "Strongest_Peak_ID": strongest.get("Peak_ID") if strongest else "",
        "Strongest_Observed_Mass": strongest.get("Observed_Mass") if strongest else None,
        "Strongest_Apex_Intensity_Raw": strongest.get("Apex_Intensity_Raw") if strongest else None,
        "Existing_Intact_Result_Available": bool(existing_masses),
        "Existing_Intact_Mass_Count": len(existing_masses),
        "Input_Identity_Audit_Status": raw_identity_status,
        "Input_Identity_Conflict": gate_status == "CONFLICT",
        "Input_Identity_Warning_Code": str(identity.get("Warning_Code") or ""),
        "Biological_Interpretation_Eligible": biological_eligible and gate_status != "CONFLICT",
        "Algorithm_Version": ALGORITHM_VERSION,
        **FORMAL_FALSE,
        "Notes": (
            f"{candidate_mode}; mass-only terminal-state comparison; no modification, "
            "structure, position, calibration, or formal-result assignment."
        ),
        "Sequence_Source": str(reference.get("Sequence_Source") or resolved_source),
        "Sequence_Length": reference.get("Sequence_Length"),
        "Sequence_SHA256": str(reference.get("Sequence_SHA256") or ""),
        "Ends_With_CCA": reference.get("Ends_With_CCA"),
        "CCA_Policy": str(reference.get("CCA_Policy") or ""),
        "Identity_Gate_Status": gate_status,
        "Identity_Gate_Passed": gate_passed,
        "Identity_Warnings": _joined_warnings(identity_warnings),
        "Observed_Mass_Type": OBSERVED_MASS_TYPE,
        "Theoretical_Mass_Type": str(reference.get("Theoretical_Mass_Type") or "UNKNOWN"),
        "Mass_Definition_Compatibility": MASS_DEFINITION_COMPATIBILITY,
        "Mass_Type_Assumption": MASS_TYPE_ASSUMPTION,
        "SCIEX_Reconstruction_Settings_Available": False,
        "Calibration_Applied": False,
        "Calibration_Method": "NONE",
        "Calibration_Reference": "",
        "Calibration_Offset_Da": 0.0,
        "Calibration_Reference_Externally_Assigned": False,
        "Terminal_Candidate_Count": len(terminal_candidates),
        "Mass_Equivalent_Group_Count": len(group_ids),
        "Eligible_Peak_Count": sum(
            value[0] != "NOT_PRIMARY_ELIGIBLE" for value in peak_classes.values()
        ),
        "Primary_Peak_Count": sum(value[0] == "PRIMARY" for value in peak_classes.values()),
        "Exploratory_Peak_Count": sum(value[0] == "EXPLORATORY" for value in peak_classes.values()),
        "Shape_Warning_Peak_Count": sum(
            value[0] == "MATCHABLE_WITH_SHAPE_WARNING" for value in peak_classes.values()
        ),
        "Nonprimary_Peak_Count": sum(
            value[0] == "NOT_PRIMARY_ELIGIBLE" for value in peak_classes.values()
        ),
        "Exploratory_Match_Count": sum(
            row["Match_Tolerance_Class"] == "EXPLORATORY" for row in details
        ),
        "Ambiguous_Match_Count": ambiguous_match_peaks,
        "No_Match_Peak_Count": no_match_peaks,
        "Exploratory_Tolerance_Da": exploratory,
    }
    for row in details:
        row.pop("_Original_Index", None)
        row.pop("_Candidate_Order", None)
    return SciexIntactMassComparisonResult(tuple(details), (summary,))
