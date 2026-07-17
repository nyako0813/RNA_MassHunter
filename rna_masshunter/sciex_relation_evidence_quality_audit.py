"""Shadow evidence-quality audit for exported SCIEX delta-mass relations."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import isfinite
from types import MappingProxyType
from typing import Any, Mapping, Sequence

AUDIT_RESULT_KEY = "sciex_relation_evidence_quality_audit"
DETAIL_SHEET = "SCIEX_Relation_Evidence"
SUMMARY_SHEET = "SCIEX_Relation_Evidence_Summary"
ERROR_CODE = "SCIEX_RELATION_EVIDENCE_QUALITY_AUDIT_ERROR"
ALGORITHM_VERSION = "sciex-relation-evidence-quality-v1"

FORMAL_FALSE = {
    "Shadow_Only": True,
    "Applied_To_Formal_Score": False,
    "Applied_To_Ranking": False,
    "Applied_To_Candidate_Filtering": False,
    "Molecular_Identity_Assigned": False,
}
DEFAULT_PARAMETERS = {
    "enabled": True,
    "high_error_fraction_threshold": 0.25,
    "low_error_fraction_threshold": 0.75,
    "minimum_recurrent_support_pairs": 2,
    "minimum_interpretable_resolution_margin": 2.0,
}

DETAIL_COLUMNS = [
    "Relation_Evidence_ID", "Source_Relation_ID", "Peak_Row_A", "Peak_Row_B",
    "Observed_Mass_A", "Observed_Mass_B", "Pair_Mass_Spacing_Da",
    "Pair_Delta_Spacing_Da", "Integer_Spacing_Candidate",
    "Nearest_Integer_Spacing", "Integer_Spacing_Error_Da",
    "Integer_Error_Fraction", "Isotope_Spacing_Candidate",
    "Nearest_Isotope_Multiple", "Isotope_Spacing_Error_Da",
    "Isotope_Error_Fraction", "Best_Numerical_Relation", "Best_Error_Fraction",
    "Numerical_Fit_Quality", "Resolution_Status", "Resolution_Ambiguous",
    "Resolution_Supported", "Effective_Resolution_Margin",
    "Recurrent_Support", "Recurrent_Numerical_Support", "Recurrent_Group_ID",
    "Recurrent_Group_Pair_Count", "Recurrent_Support_Level",
    "Recurrent_Chemical_Interpretation_Eligible",
    "Input_Identity_Audit_Status", "Input_Identity_Conflict",
    "Biological_Interpretation_Eligible", "Numerical_Evidence_Available",
    "Integer_Numerical_Evidence", "Integer_Numerical_Interpretation_Eligible",
    "Integer_Chemical_Interpretation_Eligible", "Isotope_Numerical_Proximity",
    "Isotope_Assignment_Eligible", "Isotope_Evidence_Quality",
    "Identity_Blocked", "Biological_Interpretation_Blocked",
    "Chemical_Interpretation_Blocked", "Chemical_Interpretation_Eligible",
    "Evidence_Tier", "Evidence_Quality_Label", "Interpretation_Block_Reasons",
    "High_Error_Fraction_Threshold", "Low_Error_Fraction_Threshold",
    "Minimum_Recurrent_Support_Pairs", "Minimum_Interpretable_Resolution_Margin",
    "Algorithm_Version", "Shadow_Only", "Applied_To_Formal_Score",
    "Applied_To_Ranking", "Applied_To_Candidate_Filtering",
    "Molecular_Identity_Assigned", "Notes",
]

SUMMARY_COLUMNS = [
    "Audit_Status", "Audit_Eligible", "SCIEX_Source_File", "Total_Relation_Count",
    "Tier_0_Count", "Tier_1_Count", "Tier_2_Count", "Tier_3_Count",
    "Tier_4_Count", "Numerical_Evidence_Count", "Resolution_Ambiguous_Count",
    "Resolution_Supported_Count", "Identity_Blocked_Count",
    "Biological_Interpretation_Blocked_Count",
    "Chemical_Interpretation_Blocked_Count", "Integer_Candidate_Count",
    "Isotope_Candidate_Count", "Dual_Candidate_Count",
    "Isotope_Assignment_Eligible_Count", "Recurrent_Support_Count",
    "Excellent_Fit_Count", "Good_Fit_Count", "Weak_Fit_Count",
    "Outside_Tolerance_Count", "Not_Applicable_Fit_Count",
    "Estimated_Effective_Grid_Da", "Grid_Confidence",
    "Spacing_Resolution_Status", "Spacing_Theoretically_Distinguishable",
    "Input_Identity_Audit_Status", "Input_Identity_Conflict",
    "Biological_Interpretation_Eligible", "Chemical_Interpretation_Eligible",
    "Highest_Available_Evidence_Tier", "Formal_Propagation_Enabled",
    "High_Error_Fraction_Threshold", "Low_Error_Fraction_Threshold",
    "Minimum_Recurrent_Support_Pairs", "Minimum_Interpretable_Resolution_Margin",
    "Algorithm_Version", "Shadow_Only", "Applied_To_Formal_Score",
    "Applied_To_Ranking", "Applied_To_Candidate_Filtering",
    "Molecular_Identity_Assigned", "Notes",
]

TIER_ORDER = {
    "TIER_0_UNUSABLE": 0,
    "TIER_1_NUMERICAL_ONLY": 1,
    "TIER_2_RESOLUTION_LIMITED": 2,
    "TIER_3_RESOLUTION_SUPPORTED": 3,
    "TIER_4_INTERPRETATION_ELIGIBLE": 4,
}


@dataclass(frozen=True)
class RelationEvidenceQualityParameters:
    high_error_fraction_threshold: float = 0.25
    low_error_fraction_threshold: float = 0.75
    minimum_recurrent_support_pairs: int = 2
    minimum_interpretable_resolution_margin: float = 2.0

    @classmethod
    def from_mapping(
        cls, value: Mapping[str, Any] | None,
    ) -> "RelationEvidenceQualityParameters":
        source = dict(value or {})
        source.pop("enabled", None)
        allowed = set(cls.__dataclass_fields__)
        result = cls(**{key: item for key, item in source.items() if key in allowed})
        result.validate()
        return result

    def validate(self) -> None:
        for name in (
            "high_error_fraction_threshold", "low_error_fraction_threshold",
            "minimum_interpretable_resolution_margin",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"{name} must be numeric")
            if not isfinite(float(value)) or float(value) <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if self.high_error_fraction_threshold >= self.low_error_fraction_threshold:
            raise ValueError(
                "high_error_fraction_threshold must be less than low_error_fraction_threshold"
            )
        if self.low_error_fraction_threshold > 1:
            raise ValueError("low_error_fraction_threshold must be <= 1")
        if (
            isinstance(self.minimum_recurrent_support_pairs, bool)
            or not isinstance(self.minimum_recurrent_support_pairs, int)
            or self.minimum_recurrent_support_pairs < 2
        ):
            raise ValueError("minimum_recurrent_support_pairs must be an integer >= 2")


@dataclass(frozen=True)
class SciexRelationEvidenceQualityAuditResult:
    detail_rows: tuple[Mapping[str, Any], ...]
    summary_rows: tuple[Mapping[str, Any], ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "detail_rows",
            tuple(MappingProxyType(dict(row)) for row in self.detail_rows),
        )
        object.__setattr__(
            self, "summary_rows",
            tuple(MappingProxyType(dict(row)) for row in self.summary_rows),
        )

    def details(self) -> list[dict[str, Any]]:
        return [dict(row) for row in self.detail_rows]

    def summaries(self) -> list[dict[str, Any]]:
        return [dict(row) for row in self.summary_rows]


def _records(value: Any, method: str, attribute: str) -> list[dict[str, Any]]:
    if value is None:
        return []
    source = getattr(value, method)() if hasattr(value, method) else getattr(value, attribute, value)
    if isinstance(source, Mapping):
        source = [source]
    return [dict(row) for row in source or () if isinstance(row, Mapping)]


def _finite(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if isfinite(number) else None


def _error_fraction(candidate: bool, error: Any, tolerance: Any) -> float | None:
    if not candidate:
        return None
    error_value = _finite(error)
    tolerance_value = _finite(tolerance)
    if error_value is None or tolerance_value is None or tolerance_value <= 0:
        return None
    return abs(error_value) / tolerance_value


def _fit_quality(
    fraction: float | None, params: RelationEvidenceQualityParameters,
) -> str:
    if fraction is None:
        return "NOT_APPLICABLE"
    if fraction <= params.high_error_fraction_threshold:
        return "EXCELLENT"
    if fraction <= params.low_error_fraction_threshold:
        return "GOOD"
    if fraction <= 1.0:
        return "WEAK"
    return "OUTSIDE_TOLERANCE"


def _best_relation(
    integer_fraction: float | None, isotope_fraction: float | None,
) -> tuple[str, float | None]:
    candidates = []
    if integer_fraction is not None:
        candidates.append((integer_fraction, 0, "INTEGER"))
    if isotope_fraction is not None:
        candidates.append((isotope_fraction, 1, "ISOTOPE_LIKE"))
    if not candidates:
        return "NONE", None
    fraction, _, label = min(candidates)
    return label, fraction


def _resolution_metadata(
    resolution_result: Any,
) -> tuple[dict[str, Any], dict[int, dict[str, Any]]]:
    summaries = _records(resolution_result, "summaries", "summary_rows")
    details = _records(resolution_result, "details", "detail_rows")
    if not summaries:
        raise ValueError("spacing resolution summary is missing")
    summary = summaries[0]
    required = (
        "Estimated_Effective_Grid_Da", "Grid_Confidence", "Resolution_Status",
        "Theoretically_Distinguishable", "Isotope_Interpretation_Eligible",
        "Numerical_Spacing_Interpretation_Eligible",
        "Chemical_Interpretation_Eligible", "Input_Identity_Audit_Status",
        "Input_Identity_Conflict", "Biological_Interpretation_Eligible",
    )
    missing = [key for key in required if key not in summary]
    if missing:
        raise ValueError(f"spacing resolution schema missing fields: {', '.join(missing)}")
    by_multiple = {
        int(row["Spacing_Multiple"]): row
        for row in details
        if _finite(row.get("Spacing_Multiple")) is not None
    }
    return summary, by_multiple


def _relation_resolution(
    relation: Mapping[str, Any],
    resolution_summary: Mapping[str, Any],
    resolution_details: Mapping[int, Mapping[str, Any]],
    params: RelationEvidenceQualityParameters,
) -> tuple[str, bool, bool, float | None]:
    multiples = []
    if relation.get("Integer_Spacing_Candidate"):
        value = _finite(relation.get("Nearest_Integer_Spacing"))
        if value is not None:
            multiples.append(int(value))
    if relation.get("Isotope_Spacing_Candidate"):
        value = _finite(relation.get("Nearest_Isotope_Multiple"))
        if value is not None:
            multiples.append(int(value))
    detail = next(
        (resolution_details[multiple] for multiple in multiples if multiple in resolution_details),
        None,
    )
    status = str(
        (detail or resolution_summary).get(
            "Resolution_Status",
            resolution_summary.get("Resolution_Status") or "INSUFFICIENT_INFORMATION",
        )
    )
    distinguishable = bool(
        (detail or resolution_summary).get(
            "Theoretically_Distinguishable",
            resolution_summary.get("Theoretically_Distinguishable", False),
        )
    )
    margin = _finite(
        (detail or resolution_summary).get(
            "Grid_Steps_Per_Target_Separation",
            resolution_summary.get("Single_Step_Grid_Steps_Per_Separation"),
        )
    )
    ambiguous = status.startswith("NOT_DISTINGUISHABLE")
    supported = bool(
        distinguishable
        and status in {"DISTINGUISHABLE", "MARGINALLY_DISTINGUISHABLE"}
        and margin is not None
        and margin >= params.minimum_interpretable_resolution_margin
    )
    return status, ambiguous, supported, margin


def _recurrent_level(group_id: str, count: int, minimum: int) -> tuple[bool, str]:
    if not group_id:
        return False, "NONE"
    if count < minimum:
        return False, "LOW"
    if count < minimum * 2:
        return True, "MODERATE"
    return True, "HIGH"


def _quality_label(
    unusable: bool, chemical_blocked: bool, biological_blocked: bool,
    identity_blocked: bool, ambiguous: bool, supported: bool,
) -> str:
    if unusable:
        return "INSUFFICIENT_INFORMATION"
    if chemical_blocked or biological_blocked:
        return "INTERPRETATION_BLOCKED"
    if identity_blocked:
        return "IDENTITY_BLOCKED"
    if ambiguous:
        return "RESOLUTION_AMBIGUOUS"
    if supported:
        return "RESOLUTION_SUPPORTED"
    return "NUMERICAL_ONLY"


def audit_sciex_relation_evidence_quality(
    cluster_result: Any,
    resolution_result: Any,
    cluster_parameters: Mapping[str, Any],
    parameters: RelationEvidenceQualityParameters | Mapping[str, Any] | None = None,
) -> SciexRelationEvidenceQualityAuditResult:
    params = (
        parameters if isinstance(parameters, RelationEvidenceQualityParameters)
        else RelationEvidenceQualityParameters.from_mapping(parameters)
    )
    params.validate()
    relations = _records(cluster_result, "relations", "relation_rows")
    if not relations:
        raise ValueError("exported SCIEX relation rows are missing")
    resolution, resolution_details = _resolution_metadata(resolution_result)
    integer_tolerance = cluster_parameters.get("integer_spacing_tolerance_da")
    isotope_tolerance = cluster_parameters.get("isotope_spacing_tolerance_da")
    identity_status = str(resolution.get("Input_Identity_Audit_Status") or "NOT_RUN")
    identity_conflict = bool(resolution.get("Input_Identity_Conflict", False))
    biological_eligible = bool(resolution.get("Biological_Interpretation_Eligible", False))
    group_counts = Counter(
        str(row.get("Recurrent_Spacing_Group_ID") or "")
        for row in relations if row.get("Recurrent_Spacing_Group_ID")
    )

    detail_rows = []
    for index, source in enumerate(relations, 1):
        relation_id = str(source.get("Relation_ID") or "")
        mass_a = _finite(source.get("Observed_Mass_A"))
        mass_b = _finite(source.get("Observed_Mass_B"))
        pair_spacing = _finite(source.get("Pair_Mass_Spacing_Da"))
        delta_spacing = _finite(source.get("Pair_Delta_Spacing_Da"))
        unusable = any(value is None for value in (mass_a, mass_b, pair_spacing, delta_spacing))
        integer_candidate = bool(source.get("Integer_Spacing_Candidate", False))
        isotope_candidate = bool(source.get("Isotope_Spacing_Candidate", False))
        integer_fraction = _error_fraction(
            integer_candidate, source.get("Integer_Spacing_Error_Da"), integer_tolerance,
        )
        isotope_fraction = _error_fraction(
            isotope_candidate, source.get("Isotope_Spacing_Error_Da"), isotope_tolerance,
        )
        best_relation, best_fraction = _best_relation(integer_fraction, isotope_fraction)
        fit_quality = _fit_quality(best_fraction, params)
        status, resolution_ambiguous, resolution_supported, margin = _relation_resolution(
            source, resolution, resolution_details, params,
        )
        numerical_available = bool(not unusable and pair_spacing is not None)
        candidate = integer_candidate or isotope_candidate
        ambiguous = bool(candidate and resolution_ambiguous)
        supported = bool(candidate and resolution_supported)
        identity_blocked = identity_conflict or identity_status != "MATCH"
        biological_blocked = not biological_eligible
        chemical_blocked = True
        group_id = str(source.get("Recurrent_Spacing_Group_ID") or "")
        group_count = group_counts.get(group_id, 0) if group_id else 0
        recurrent_support, recurrent_level = _recurrent_level(
            group_id, group_count, params.minimum_recurrent_support_pairs,
        )

        if unusable:
            tier = "TIER_0_UNUSABLE"
        elif candidate and ambiguous:
            tier = "TIER_2_RESOLUTION_LIMITED"
        elif candidate and supported:
            tier = "TIER_3_RESOLUTION_SUPPORTED"
        else:
            tier = "TIER_1_NUMERICAL_ONLY"

        reasons = []
        if unusable:
            reasons.append("MISSING_OR_NONFINITE_RELATION_VALUES")
        if status == "INSUFFICIENT_INFORMATION":
            reasons.append("RESOLUTION_METADATA_INSUFFICIENT")
        if ambiguous:
            reasons.append("RESOLUTION_NOT_DISTINGUISHABLE")
        if identity_status == "NOT_RUN":
            reasons.append("INPUT_IDENTITY_UNAVAILABLE")
        elif identity_conflict:
            reasons.append("INPUT_IDENTITY_CONFLICT")
        if biological_blocked:
            reasons.append("BIOLOGICAL_INTERPRETATION_BLOCKED")
        reasons.append("CHEMICAL_ASSIGNMENT_DISABLED")

        detail_rows.append({
            "Relation_Evidence_ID": f"SCIEX_REL_EVID_{index:05d}",
            "Source_Relation_ID": relation_id,
            "Peak_Row_A": source.get("Peak_Row_A"),
            "Peak_Row_B": source.get("Peak_Row_B"),
            "Observed_Mass_A": source.get("Observed_Mass_A"),
            "Observed_Mass_B": source.get("Observed_Mass_B"),
            "Pair_Mass_Spacing_Da": source.get("Pair_Mass_Spacing_Da"),
            "Pair_Delta_Spacing_Da": source.get("Pair_Delta_Spacing_Da"),
            "Integer_Spacing_Candidate": integer_candidate,
            "Nearest_Integer_Spacing": source.get("Nearest_Integer_Spacing"),
            "Integer_Spacing_Error_Da": source.get("Integer_Spacing_Error_Da"),
            "Integer_Error_Fraction": integer_fraction,
            "Isotope_Spacing_Candidate": isotope_candidate,
            "Nearest_Isotope_Multiple": source.get("Nearest_Isotope_Multiple"),
            "Isotope_Spacing_Error_Da": source.get("Isotope_Spacing_Error_Da"),
            "Isotope_Error_Fraction": isotope_fraction,
            "Best_Numerical_Relation": best_relation,
            "Best_Error_Fraction": best_fraction,
            "Numerical_Fit_Quality": fit_quality,
            "Resolution_Status": status,
            "Resolution_Ambiguous": ambiguous,
            "Resolution_Supported": supported,
            "Effective_Resolution_Margin": margin,
            "Recurrent_Support": recurrent_support,
            "Recurrent_Numerical_Support": recurrent_support,
            "Recurrent_Group_ID": group_id,
            "Recurrent_Group_Pair_Count": group_count,
            "Recurrent_Support_Level": recurrent_level,
            "Recurrent_Chemical_Interpretation_Eligible": False,
            "Input_Identity_Audit_Status": identity_status,
            "Input_Identity_Conflict": identity_conflict,
            "Biological_Interpretation_Eligible": biological_eligible,
            "Numerical_Evidence_Available": numerical_available,
            "Integer_Numerical_Evidence": integer_candidate,
            "Integer_Numerical_Interpretation_Eligible": bool(
                numerical_available and integer_candidate
            ),
            "Integer_Chemical_Interpretation_Eligible": False,
            "Isotope_Numerical_Proximity": isotope_candidate,
            "Isotope_Assignment_Eligible": False,
            "Isotope_Evidence_Quality": (
                "RESOLUTION_AMBIGUOUS" if isotope_candidate and ambiguous
                else "RESOLUTION_SUPPORTED" if isotope_candidate and supported
                else "NOT_APPLICABLE"
            ),
            "Identity_Blocked": identity_blocked,
            "Biological_Interpretation_Blocked": biological_blocked,
            "Chemical_Interpretation_Blocked": chemical_blocked,
            "Chemical_Interpretation_Eligible": False,
            "Evidence_Tier": tier,
            "Evidence_Quality_Label": _quality_label(
                unusable, chemical_blocked, biological_blocked,
                identity_blocked, ambiguous, supported,
            ),
            "Interpretation_Block_Reasons": "; ".join(reasons),
            "High_Error_Fraction_Threshold": params.high_error_fraction_threshold,
            "Low_Error_Fraction_Threshold": params.low_error_fraction_threshold,
            "Minimum_Recurrent_Support_Pairs": params.minimum_recurrent_support_pairs,
            "Minimum_Interpretable_Resolution_Margin": (
                params.minimum_interpretable_resolution_margin
            ),
            "Algorithm_Version": ALGORITHM_VERSION,
            **FORMAL_FALSE,
            "Notes": (
                "Derived from the existing exported relation row; source relation flags, "
                "values, order, and truncation are unchanged."
            ),
        })

    tier_counts = Counter(row["Evidence_Tier"] for row in detail_rows)
    fit_counts = Counter(row["Numerical_Fit_Quality"] for row in detail_rows)
    highest = max(
        tier_counts, key=lambda tier: TIER_ORDER[tier],
    ) if tier_counts else "TIER_0_UNUSABLE"
    total = len(detail_rows)
    if sum(tier_counts.values()) != total:
        raise RuntimeError("evidence tier aggregation count mismatch")
    cluster_summaries = _records(cluster_result, "summaries", "summary_rows")
    cluster_summary = cluster_summaries[0] if cluster_summaries else {}
    summary = {
        "Audit_Status": "AUDIT_COMPLETED",
        "Audit_Eligible": True,
        "SCIEX_Source_File": str(
            resolution.get("SCIEX_Source_File")
            or cluster_summary.get("SCIEX_Source_File") or ""
        ),
        "Total_Relation_Count": total,
        "Tier_0_Count": tier_counts["TIER_0_UNUSABLE"],
        "Tier_1_Count": tier_counts["TIER_1_NUMERICAL_ONLY"],
        "Tier_2_Count": tier_counts["TIER_2_RESOLUTION_LIMITED"],
        "Tier_3_Count": tier_counts["TIER_3_RESOLUTION_SUPPORTED"],
        "Tier_4_Count": tier_counts["TIER_4_INTERPRETATION_ELIGIBLE"],
        "Numerical_Evidence_Count": sum(
            row["Numerical_Evidence_Available"] for row in detail_rows
        ),
        "Resolution_Ambiguous_Count": sum(
            row["Resolution_Ambiguous"] for row in detail_rows
        ),
        "Resolution_Supported_Count": sum(
            row["Resolution_Supported"] for row in detail_rows
        ),
        "Identity_Blocked_Count": sum(row["Identity_Blocked"] for row in detail_rows),
        "Biological_Interpretation_Blocked_Count": sum(
            row["Biological_Interpretation_Blocked"] for row in detail_rows
        ),
        "Chemical_Interpretation_Blocked_Count": sum(
            row["Chemical_Interpretation_Blocked"] for row in detail_rows
        ),
        "Integer_Candidate_Count": sum(
            row["Integer_Spacing_Candidate"] for row in detail_rows
        ),
        "Isotope_Candidate_Count": sum(
            row["Isotope_Spacing_Candidate"] for row in detail_rows
        ),
        "Dual_Candidate_Count": sum(
            row["Integer_Spacing_Candidate"] and row["Isotope_Spacing_Candidate"]
            for row in detail_rows
        ),
        "Isotope_Assignment_Eligible_Count": 0,
        "Recurrent_Support_Count": sum(row["Recurrent_Support"] for row in detail_rows),
        "Excellent_Fit_Count": fit_counts["EXCELLENT"],
        "Good_Fit_Count": fit_counts["GOOD"],
        "Weak_Fit_Count": fit_counts["WEAK"],
        "Outside_Tolerance_Count": fit_counts["OUTSIDE_TOLERANCE"],
        "Not_Applicable_Fit_Count": fit_counts["NOT_APPLICABLE"],
        "Estimated_Effective_Grid_Da": resolution.get("Estimated_Effective_Grid_Da"),
        "Grid_Confidence": resolution.get("Grid_Confidence"),
        "Spacing_Resolution_Status": resolution.get("Resolution_Status"),
        "Spacing_Theoretically_Distinguishable": bool(
            resolution.get("Theoretically_Distinguishable", False)
        ),
        "Input_Identity_Audit_Status": identity_status,
        "Input_Identity_Conflict": identity_conflict,
        "Biological_Interpretation_Eligible": biological_eligible,
        "Chemical_Interpretation_Eligible": False,
        "Highest_Available_Evidence_Tier": highest,
        "Formal_Propagation_Enabled": False,
        "High_Error_Fraction_Threshold": params.high_error_fraction_threshold,
        "Low_Error_Fraction_Threshold": params.low_error_fraction_threshold,
        "Minimum_Recurrent_Support_Pairs": params.minimum_recurrent_support_pairs,
        "Minimum_Interpretable_Resolution_Margin": (
            params.minimum_interpretable_resolution_margin
        ),
        "Algorithm_Version": ALGORITHM_VERSION,
        **FORMAL_FALSE,
        "Notes": (
            "Relation-level shadow evidence quality only; isotope, chemical, biological, "
            "and molecular identity assignment remain disabled."
        ),
    }
    return SciexRelationEvidenceQualityAuditResult(tuple(detail_rows), (summary,))


def annotate_cluster_summary(cluster_result: Any, evidence_result: Any) -> Any:
    from rna_masshunter.sciex_delta_mass_cluster_audit import (
        SciexDeltaMassClusterAuditResult,
    )
    summaries = _records(cluster_result, "summaries", "summary_rows")
    evidence = _records(evidence_result, "summaries", "summary_rows")
    metadata = evidence[0] if evidence else {}
    updated = []
    for row in summaries:
        copy = dict(row)
        copy.update({
            "Relation_Evidence_Audit_Status": metadata.get("Audit_Status", "NOT_RUN"),
            "Highest_Available_Evidence_Tier": metadata.get(
                "Highest_Available_Evidence_Tier", ""
            ),
            "Resolution_Ambiguous_Relation_Count": metadata.get(
                "Resolution_Ambiguous_Count", 0
            ),
            "Interpretation_Eligible_Relation_Count": metadata.get("Tier_4_Count", 0),
        })
        updated.append(copy)
    return SciexDeltaMassClusterAuditResult(
        tuple(_records(cluster_result, "clusters", "cluster_rows")),
        tuple(updated),
        tuple(_records(cluster_result, "relations", "relation_rows")),
    )


def annotate_resolution_summary(resolution_result: Any, evidence_result: Any) -> Any:
    from rna_masshunter.sciex_spacing_resolution_audit import (
        SciexSpacingResolutionAuditResult,
    )
    summaries = _records(resolution_result, "summaries", "summary_rows")
    details = _records(resolution_result, "details", "detail_rows")
    evidence = _records(evidence_result, "summaries", "summary_rows")
    metadata = evidence[0] if evidence else {}
    updated = []
    for row in summaries:
        copy = dict(row)
        copy.update({
            "Relation_Evidence_Audit_Status": metadata.get("Audit_Status", "NOT_RUN"),
            "Tier_2_Resolution_Limited_Count": metadata.get("Tier_2_Count", 0),
            "Tier_3_Resolution_Supported_Count": metadata.get("Tier_3_Count", 0),
            "Tier_4_Interpretation_Eligible_Count": metadata.get("Tier_4_Count", 0),
        })
        updated.append(copy)
    return SciexSpacingResolutionAuditResult(tuple(updated), tuple(details))
