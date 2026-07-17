"""Shadow-only mass-proximity comparison for detected SCIEX intact peaks."""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from types import MappingProxyType
from typing import Any, Iterable, Mapping

ALGORITHM_VERSION = "sciex-intact-mass-comparison-v1"
FORMAL_FALSE = {
    "Shadow_Only": True,
    "Applied_To_Formal_Score": False,
    "Applied_To_Ranking": False,
    "Applied_To_Candidate_Filtering": False,
    "Molecular_Identity_Assigned": False,
    "Modification_Lookup_Performed": False,
}
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


def compare_sciex_intact_masses(
    detection_result: Any,
    theoretical_unmodified_mass: float | None,
    existing_intact_results: Iterable[Any] | None = None,
    *,
    source_file: str = "",
    strict_tolerance_da: float = 1.0,
    broad_tolerance_da: float = 5.0,
    input_identity_audit: Any = None,
) -> SciexIntactMassComparisonResult:
    """Compare apex masses without modification lookup or identity assignment."""
    strict = float(strict_tolerance_da)
    broad = float(broad_tolerance_da)
    if not isfinite(strict) or strict <= 0:
        raise ValueError("strict_tolerance_da must be finite and positive")
    if not isfinite(broad) or broad < strict:
        raise ValueError("broad_tolerance_da must be finite and >= strict_tolerance_da")

    theory = _finite_float(theoretical_unmodified_mass)
    if theory is not None and theory <= 0:
        theory = None
    peaks = _peak_rows(detection_result)
    diagnostics = _diagnostics(detection_result)
    detector_completed = diagnostics.get("Detection_Status") == "DETECTION_COMPLETED"
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

    details: list[dict[str, Any]] = []
    for position, (original_index, peak, observed_mass) in enumerate(indexed_peaks, 1):
        eligible = detector_completed and theory is not None
        delta = observed_mass - theory if eligible else None
        absolute_delta = abs(delta) if delta is not None else None
        if not detector_completed:
            status = "NOT_ELIGIBLE"
        elif theory is None:
            status = "NO_THEORETICAL_MASS"
        elif absolute_delta <= strict:
            status = "STRICT_MATCH"
        elif absolute_delta <= broad:
            status = "BROAD_MATCH"
        else:
            status = "NO_MATCH"

        nearest = (
            min(existing_masses, key=lambda mass: (abs(observed_mass - mass), mass))
            if existing_masses else None
        )
        existing_delta = observed_mass - nearest if nearest is not None else None
        row = {
            "Comparison_ID": f"SCIEX_CMP_{position:05d}",
            "Peak_ID": peak.get("Peak_ID", ""),
            "Source_File": str(source_file),
            "Observed_Mass": observed_mass,
            "Theoretical_Unmodified_Mass": theory,
            "Delta_Mass": delta,
            "Absolute_Delta_Mass": absolute_delta,
            "Delta_ppm": delta / theory * 1_000_000 if delta is not None else None,
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
            "Comparison_Status": status,
            "Comparison_Eligible": eligible,
            "Strict_Tolerance_Da": strict,
            "Broad_Tolerance_Da": broad,
            **{column: peak.get(column) for column in PEAK_METADATA_COLUMNS},
            "Algorithm_Version": ALGORITHM_VERSION,
            **FORMAL_FALSE,
            "Notes": "Mass proximity only; no modification or molecular identity assignment.",
            "_Original_Index": original_index,
        }
        details.append(row)

    closest = None
    if detector_completed and theory is not None and details:
        closest = min(
            details,
            key=lambda row: (
                row["Absolute_Delta_Mass"], -_intensity(row),
                row["Observed_Mass"], row["_Original_Index"],
            ),
        )
    strongest = min(
        details,
        key=lambda row: (-_intensity(row), row["Observed_Mass"], row["_Original_Index"]),
    ) if details else None
    for row in details:
        row.pop("_Original_Index", None)

    if not detector_completed or not details:
        summary_status = "NOT_ELIGIBLE"
    elif theory is None:
        summary_status = "NO_THEORETICAL_MASS"
    else:
        summary_status = "COMPARISON_COMPLETED"
    identity_value = input_identity_audit
    if hasattr(identity_value, "row"):
        identity_value = identity_value.row()
    identity = dict(identity_value) if isinstance(identity_value, Mapping) else {}
    identity_status = str(identity.get("Audit_Status") or "NOT_RUN")
    biological_eligible = bool(identity.get("Biological_Interpretation_Eligible", False))
    summary = {
        "Source_File": str(source_file),
        "Parser_Status": str(diagnostics.get("Input_Status") or ""),
        "Detector_Status": str(diagnostics.get("Detection_Status") or ""),
        "Comparison_Status": summary_status,
        "Comparison_Eligible": detector_completed and theory is not None and bool(details),
        "Detected_Peak_Count": len(peaks),
        "Compared_Peak_Count": len(details),
        "Strict_Match_Count": sum(row["Comparison_Status"] == "STRICT_MATCH" for row in details),
        "Broad_Match_Count": sum(row["Comparison_Status"] == "BROAD_MATCH" for row in details),
        "No_Match_Count": sum(row["Comparison_Status"] == "NO_MATCH" for row in details),
        "Theoretical_Unmodified_Mass": theory,
        "Strict_Tolerance_Da": strict,
        "Broad_Tolerance_Da": broad,
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
        "Input_Identity_Audit_Status": identity_status,
        "Input_Identity_Conflict": bool(identity.get("Identity_Conflict", False)),
        "Input_Identity_Warning_Code": str(identity.get("Warning_Code") or ""),
        "Biological_Interpretation_Eligible": biological_eligible,
        "Algorithm_Version": ALGORITHM_VERSION,
        **FORMAL_FALSE,
        "Notes": (
            "Mass proximity only; no modification lookup, chemical assignment, "
            "or formal-result propagation."
        ),
    }
    return SciexIntactMassComparisonResult(tuple(details), (summary,))
