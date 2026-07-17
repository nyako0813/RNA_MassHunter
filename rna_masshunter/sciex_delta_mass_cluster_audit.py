"""Numerical shadow clustering and spacing audit for SCIEX intact delta masses."""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from statistics import mean, median
from types import MappingProxyType
from typing import Any, Mapping, Sequence
import re

AUDIT_RESULT_KEY = "sciex_delta_mass_cluster_audit"
CLUSTER_SHEET = "SCIEX_Delta_Mass_Clusters"
SUMMARY_SHEET = "SCIEX_Delta_Mass_Clust_Summary"
RELATION_SHEET = "SCIEX_Delta_Mass_Relations"
ERROR_CODE = "SCIEX_DELTA_MASS_CLUSTER_AUDIT_ERROR"
ALGORITHM_VERSION = "sciex-delta-mass-cluster-v1"
FORMAL_FALSE = {
    "Shadow_Only": True,
    "Applied_To_Formal_Score": False,
    "Applied_To_Ranking": False,
    "Applied_To_Candidate_Filtering": False,
    "Molecular_Identity_Assigned": False,
}

DEFAULT_PARAMETERS = {
    "enabled": True,
    "cluster_tolerance_da": 0.5,
    "duplicate_apex_tolerance_da": 0.25,
    "isotope_spacing_da": 1.003355,
    "isotope_spacing_tolerance_da": 0.15,
    "integer_spacing_tolerance_da": 0.15,
    "minimum_cluster_size": 2,
    "max_pair_spacing_da": 200.0,
    "max_pair_rows": 20000,
}

CLUSTER_COLUMNS = [
    "Cluster_ID", "Cluster_Label", "Cluster_Size", "Cluster_Min_Delta_Da",
    "Cluster_Max_Delta_Da", "Cluster_Span_Da", "Cluster_Mean_Delta_Da",
    "Cluster_Median_Delta_Da", "Cluster_Weighted_Mean_Delta_Da",
    "Cluster_Weighted_Mean_Fallback_Used", "Cluster_Total_Intensity",
    "Cluster_Max_Intensity", "Cluster_Strongest_Observed_Mass",
    "Cluster_Closest_To_Theoretical_Mass", "Cluster_Closest_Absolute_Delta_Da",
    "Cluster_Detection_Modes", "Cluster_Comparison_Statuses", "Cluster_Is_Singleton",
    "Cluster_Is_Multi_Peak", "Duplicate_Group_ID", "Duplicate_Group_Size",
    "Duplicate_Like", "Duplicate_Mass_Span_Da", "Duplicate_Strongest_Row",
    "Duplicate_Representative_Row", "Integer_Spacing_Series",
    "Isotope_Spacing_Candidate", "Recurrent_Spacing_Candidate",
    "Integer_Spacing_Relation_Count", "Isotope_Spacing_Relation_Count",
    "Recurrent_Spacing_Relation_Count", "Member_Row_Indices",
    "Member_Observed_Masses", "Member_Delta_Masses", "Algorithm_Version",
    "Shadow_Only", "Applied_To_Formal_Score", "Applied_To_Ranking",
    "Applied_To_Candidate_Filtering", "Molecular_Identity_Assigned", "Notes",
]

SUMMARY_COLUMNS = [
    "Audit_Status", "Audit_Eligible", "SCIEX_Source_File",
    "Theoretical_Unmodified_Mass", "Detected_Peak_Count", "Cluster_Count",
    "Singleton_Cluster_Count", "Multi_Peak_Cluster_Count",
    "Duplicate_Like_Cluster_Count", "Integer_Spacing_Cluster_Count",
    "Isotope_Spacing_Candidate_Count", "Recurrent_Spacing_Group_Count",
    "Recurrent_Spacing_Cluster_Count", "Largest_Cluster_Size", "Largest_Cluster_ID",
    "Strongest_Cluster_ID", "Closest_Cluster_ID", "Closest_Cluster_Median_Delta_Da",
    "Closest_Peak_Row", "Closest_Peak_Cluster_ID", "Strongest_Peak_Row",
    "Strongest_Peak_Cluster_ID", "Pair_Analysis_Performed",
    "Total_Eligible_Pair_Count", "Total_Relevant_Pair_Count", "Exported_Pair_Count",
    "Pair_Rows_Truncated", "Input_Identity_Audit_Status", "Input_Identity_Conflict",
    "Biological_Interpretation_Eligible", "Cluster_Tolerance_Da",
    "Duplicate_Apex_Tolerance_Da", "Isotope_Spacing_Da",
    "Isotope_Spacing_Tolerance_Da", "Integer_Spacing_Tolerance_Da",
    "Minimum_Cluster_Size", "Max_Pair_Spacing_Da", "Max_Pair_Rows",
    "Algorithm_Version", "Shadow_Only", "Applied_To_Formal_Score",
    "Applied_To_Ranking", "Applied_To_Candidate_Filtering",
    "Molecular_Identity_Assigned", "Notes",
]

RELATION_COLUMNS = [
    "Relation_ID", "Relation_Types", "Peak_Row_A", "Peak_Row_B", "Observed_Mass_A",
    "Observed_Mass_B", "Delta_Mass_A", "Delta_Mass_B", "Pair_Mass_Spacing_Da",
    "Pair_Delta_Spacing_Da", "Combined_Apex_Intensity", "Duplicate_Like",
    "Integer_Spacing_Candidate", "Nearest_Integer_Spacing",
    "Integer_Spacing_Error_Da", "Isotope_Spacing_Candidate",
    "Nearest_Isotope_Multiple", "Isotope_Spacing_Error_Da",
    "Recurrent_Spacing_Group_ID", "Shadow_Only", "Applied_To_Formal_Score",
    "Applied_To_Ranking", "Applied_To_Candidate_Filtering",
    "Molecular_Identity_Assigned", "Notes",
]


@dataclass(frozen=True)
class DeltaMassClusterParameters:
    cluster_tolerance_da: float = 0.5
    duplicate_apex_tolerance_da: float = 0.25
    isotope_spacing_da: float = 1.003355
    isotope_spacing_tolerance_da: float = 0.15
    integer_spacing_tolerance_da: float = 0.15
    minimum_cluster_size: int = 2
    max_pair_spacing_da: float = 200.0
    max_pair_rows: int = 20000

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any] | None) -> "DeltaMassClusterParameters":
        source = dict(value or {})
        source.pop("enabled", None)
        allowed = {name for name in cls.__dataclass_fields__}
        result = cls(**{key: item for key, item in source.items() if key in allowed})
        result.validate()
        return result

    def validate(self) -> None:
        positive = (
            "cluster_tolerance_da", "isotope_spacing_da",
            "isotope_spacing_tolerance_da", "integer_spacing_tolerance_da",
            "max_pair_spacing_da",
        )
        for name in positive:
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"{name} must be numeric")
            if not isfinite(float(value)) or float(value) <= 0:
                raise ValueError(f"{name} must be finite and positive")
        duplicate = self.duplicate_apex_tolerance_da
        if isinstance(duplicate, bool) or not isinstance(duplicate, (int, float)):
            raise ValueError("duplicate_apex_tolerance_da must be numeric")
        if not isfinite(float(duplicate)) or float(duplicate) < 0:
            raise ValueError("duplicate_apex_tolerance_da must be finite and nonnegative")
        if isinstance(self.minimum_cluster_size, bool) or not isinstance(self.minimum_cluster_size, int) or self.minimum_cluster_size < 2:
            raise ValueError("minimum_cluster_size must be an integer >= 2")
        if isinstance(self.max_pair_rows, bool) or not isinstance(self.max_pair_rows, int) or self.max_pair_rows < 1:
            raise ValueError("max_pair_rows must be an integer >= 1")


@dataclass(frozen=True)
class SciexDeltaMassClusterAuditResult:
    cluster_rows: tuple[Mapping[str, Any], ...]
    summary_rows: tuple[Mapping[str, Any], ...]
    relation_rows: tuple[Mapping[str, Any], ...]

    def __post_init__(self) -> None:
        for name in ("cluster_rows", "summary_rows", "relation_rows"):
            object.__setattr__(
                self, name, tuple(MappingProxyType(dict(row)) for row in getattr(self, name)),
            )

    def clusters(self) -> list[dict[str, Any]]:
        return [dict(row) for row in self.cluster_rows]

    def summaries(self) -> list[dict[str, Any]]:
        return [dict(row) for row in self.summary_rows]

    def relations(self) -> list[dict[str, Any]]:
        return [dict(row) for row in self.relation_rows]


def _finite(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if isfinite(number) else None


def _records(value: Any, method: str, attribute: str) -> list[dict[str, Any]]:
    if value is None:
        return []
    source = getattr(value, method)() if hasattr(value, method) else getattr(value, attribute, value)
    if isinstance(source, Mapping):
        source = [source]
    return [dict(row) for row in source or () if isinstance(row, Mapping)]


def _source_index(row: Mapping[str, Any], fallback: int) -> int:
    explicit = _finite(row.get("Source_Row_Index"))
    if explicit is not None and explicit >= 0:
        return int(explicit)
    identifier = str(row.get("Comparison_ID") or row.get("Peak_ID") or "")
    match = re.search(r"(\d+)$", identifier)
    return int(match.group(1)) if match else fallback


def _normalize_rows(detail_rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    temporary = []
    for position, source in enumerate(detail_rows, 1):
        observed = _finite(source.get("Observed_Mass"))
        delta = _finite(source.get("Delta_Mass", source.get("Observed_Delta_Mass")))
        if observed is None or delta is None:
            continue
        intensity_raw = _finite(source.get("Apex_Intensity_Raw"))
        temporary.append({
            "Observed_Mass": observed,
            "Delta_Mass": delta,
            "Absolute_Delta_Mass": abs(delta),
            "Intensity": max(0.0, intensity_raw) if intensity_raw is not None else 0.0,
            "Intensity_Available": intensity_raw is not None and intensity_raw >= 0,
            "Detection_Mode": str(source.get("Detection_Tier") or source.get("Detection_Mode") or "UNKNOWN"),
            "Peak_Width": _finite(source.get("Half_Prominence_Width_Da", source.get("Peak_Width"))),
            "Comparison_Status": str(source.get("Comparison_Status") or ""),
            "Source_Row_Index": _source_index(source, position),
            "Comparison_ID": str(source.get("Comparison_ID") or ""),
        })
    temporary.sort(key=lambda row: (
        row["Delta_Mass"], row["Observed_Mass"], row["Source_Row_Index"], row["Comparison_ID"],
    ))
    for position, row in enumerate(temporary, 1):
        row["Audit_Row_ID"] = f"SCIEX_DELTA_ROW_{position:05d}"
    return temporary


def _complete_link_groups(rows: Sequence[Any], value_key, tolerance: float) -> list[list[Any]]:
    groups: list[list[Any]] = []
    current: list[Any] = []
    minimum = 0.0
    for row in rows:
        value = float(value_key(row))
        if not current or value - minimum <= tolerance:
            if not current:
                minimum = value
            current.append(row)
        else:
            groups.append(current)
            current = [row]
            minimum = value
    if current:
        groups.append(current)
    return groups


def _detection_rank(value: str) -> int:
    return {"STRICT": 0, "SENSITIVE": 1}.get(str(value).upper(), 2)


def _representative_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    width = row.get("Peak_Width")
    return (
        -float(row["Intensity"]), _detection_rank(str(row["Detection_Mode"])),
        float(width) if width is not None else float("inf"),
        float(row["Observed_Mass"]), int(row["Source_Row_Index"]),
    )


def _strongest_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    return (-float(row["Intensity"]), float(row["Observed_Mass"]), int(row["Source_Row_Index"]))


def _format_numbers(values: Sequence[float]) -> str:
    return ";".join(f"{value:.9f}" for value in values)


def _relation_priority(row: Mapping[str, Any]) -> tuple[int, float]:
    if row["Duplicate_Like"]:
        return 0, float(row["Pair_Mass_Spacing_Da"])
    if row["Isotope_Spacing_Candidate"]:
        return 1, float(row["Isotope_Spacing_Error_Da"])
    if row["Integer_Spacing_Candidate"]:
        return 2, float(row["Integer_Spacing_Error_Da"])
    return 3, 0.0


def audit_sciex_delta_mass_clusters(
    comparison_result: Any,
    parameters: DeltaMassClusterParameters | Mapping[str, Any] | None = None,
) -> SciexDeltaMassClusterAuditResult:
    """Cluster delta masses and diagnose spacing without chemical assignment."""
    params = parameters if isinstance(parameters, DeltaMassClusterParameters) else DeltaMassClusterParameters.from_mapping(parameters)
    params.validate()
    detail_rows = _records(comparison_result, "details", "detail_rows")
    comparison_summaries = _records(comparison_result, "summaries", "summary_rows")
    comparison_summary = comparison_summaries[0] if comparison_summaries else {}
    theoretical_mass = _finite(comparison_summary.get("Theoretical_Unmodified_Mass"))
    rows = _normalize_rows(detail_rows)
    if theoretical_mass is None or theoretical_mass <= 0 or not rows:
        raise ValueError("comparison must contain theoretical mass and at least one valid detail row")

    duplicate_groups: list[list[dict[str, Any]]] = []
    mass_sorted = sorted(rows, key=lambda row: (
        row["Observed_Mass"], row["Source_Row_Index"], row["Audit_Row_ID"],
    ))
    for group in _complete_link_groups(
        mass_sorted, lambda row: row["Observed_Mass"], params.duplicate_apex_tolerance_da,
    ):
        if len(group) >= 2:
            duplicate_groups.append(group)
    duplicate_by_row: dict[str, dict[str, Any]] = {}
    for position, group in enumerate(duplicate_groups, 1):
        group_id = f"SCIEX_DUP_{position:05d}"
        strongest = min(group, key=_strongest_key)
        representative = min(group, key=_representative_key)
        info = {
            "Duplicate_Group_ID": group_id,
            "Duplicate_Group_Size": len(group),
            "Duplicate_Mass_Span_Da": max(row["Observed_Mass"] for row in group) - min(row["Observed_Mass"] for row in group),
            "Duplicate_Strongest_Row": strongest["Source_Row_Index"],
            "Duplicate_Representative_Row": representative["Source_Row_Index"],
        }
        for row in group:
            duplicate_by_row[row["Audit_Row_ID"]] = info

    eligible_pairs: list[dict[str, Any]] = []
    for left_index, left in enumerate(mass_sorted):
        for right in mass_sorted[left_index + 1:]:
            spacing = right["Observed_Mass"] - left["Observed_Mass"]
            if spacing > params.max_pair_spacing_da:
                break
            nearest_integer = max(1, int(round(spacing)))
            integer_error = abs(spacing - nearest_integer)
            integer_candidate = integer_error <= params.integer_spacing_tolerance_da
            nearest_isotope = max(1, int(round(spacing / params.isotope_spacing_da)))
            isotope_error = abs(spacing - nearest_isotope * params.isotope_spacing_da)
            isotope_candidate = isotope_error <= params.isotope_spacing_tolerance_da
            duplicate_like = spacing <= params.duplicate_apex_tolerance_da
            eligible_pairs.append({
                "Peak_Row_A": left["Source_Row_Index"],
                "Peak_Row_B": right["Source_Row_Index"],
                "Audit_Row_ID_A": left["Audit_Row_ID"],
                "Audit_Row_ID_B": right["Audit_Row_ID"],
                "Observed_Mass_A": left["Observed_Mass"],
                "Observed_Mass_B": right["Observed_Mass"],
                "Delta_Mass_A": left["Delta_Mass"],
                "Delta_Mass_B": right["Delta_Mass"],
                "Pair_Mass_Spacing_Da": spacing,
                "Pair_Delta_Spacing_Da": abs(right["Delta_Mass"] - left["Delta_Mass"]),
                "Combined_Apex_Intensity": left["Intensity"] + right["Intensity"],
                "Duplicate_Like": duplicate_like,
                "Integer_Spacing_Candidate": integer_candidate,
                "Nearest_Integer_Spacing": nearest_integer if integer_candidate else "",
                "Integer_Spacing_Error_Da": integer_error if integer_candidate else "",
                "Isotope_Spacing_Candidate": isotope_candidate,
                "Nearest_Isotope_Multiple": nearest_isotope if isotope_candidate else "",
                "Isotope_Spacing_Error_Da": isotope_error if isotope_candidate else "",
                "Recurrent_Spacing_Group_ID": "",
            })

    spacing_sorted = sorted(eligible_pairs, key=lambda row: (
        row["Pair_Mass_Spacing_Da"], row["Peak_Row_A"], row["Peak_Row_B"],
    ))
    recurrent_tolerance = min(
        params.integer_spacing_tolerance_da, params.isotope_spacing_tolerance_da,
    )
    recurrent_groups = []
    for group in _complete_link_groups(
        spacing_sorted, lambda row: row["Pair_Mass_Spacing_Da"], recurrent_tolerance,
    ):
        if len(group) >= params.minimum_cluster_size:
            recurrent_groups.append(group)
    for position, group in enumerate(recurrent_groups, 1):
        group_id = f"SCIEX_RECUR_{position:05d}"
        for relation in group:
            relation["Recurrent_Spacing_Group_ID"] = group_id

    relevant_relations = [
        row for row in eligible_pairs
        if row["Duplicate_Like"] or row["Integer_Spacing_Candidate"]
        or row["Isotope_Spacing_Candidate"] or row["Recurrent_Spacing_Group_ID"]
    ]
    for relation in relevant_relations:
        types = []
        if relation["Duplicate_Like"]:
            types.append("DUPLICATE_LIKE")
        if relation["Isotope_Spacing_Candidate"]:
            types.append("ISOTOPE_SPACING_CANDIDATE")
        if relation["Integer_Spacing_Candidate"]:
            types.append("INTEGER_SPACING_CANDIDATE")
        if relation["Recurrent_Spacing_Group_ID"]:
            types.append("RECURRENT_SPACING")
        relation["Relation_Types"] = ";".join(types)
    relevant_relations.sort(key=lambda row: (
        *_relation_priority(row), -row["Combined_Apex_Intensity"],
        row["Observed_Mass_A"], row["Observed_Mass_B"],
        row["Peak_Row_A"], row["Peak_Row_B"],
    ))
    pair_rows_truncated = len(relevant_relations) > params.max_pair_rows
    exported_relations = relevant_relations[:params.max_pair_rows]
    relation_rows = []
    for position, relation in enumerate(exported_relations, 1):
        relation_rows.append({
            "Relation_ID": f"SCIEX_REL_{position:05d}",
            **{column: relation.get(column, "") for column in RELATION_COLUMNS if column not in {
                "Relation_ID", "Shadow_Only", "Applied_To_Formal_Score", "Applied_To_Ranking",
                "Applied_To_Candidate_Filtering", "Molecular_Identity_Assigned", "Notes",
            }},
            **FORMAL_FALSE,
            "Notes": "Numerical spacing candidate only; no isotope, adduct, or chemical identity assigned.",
        })

    incident: dict[str, list[dict[str, Any]]] = {row["Audit_Row_ID"]: [] for row in rows}
    for relation in relevant_relations:
        incident[relation["Audit_Row_ID_A"]].append(relation)
        incident[relation["Audit_Row_ID_B"]].append(relation)

    delta_groups = _complete_link_groups(
        rows, lambda row: row["Delta_Mass"], params.cluster_tolerance_da,
    )
    cluster_rows = []
    cluster_by_audit_row: dict[str, str] = {}
    for position, group in enumerate(delta_groups, 1):
        cluster_id = f"SCIEX_DELTA_C{position:05d}"
        for member in group:
            cluster_by_audit_row[member["Audit_Row_ID"]] = cluster_id
        deltas = [row["Delta_Mass"] for row in group]
        intensities = [row["Intensity"] for row in group]
        total_intensity = sum(intensities)
        weighted_fallback = total_intensity <= 0
        weighted_mean = (
            sum(delta * weight for delta, weight in zip(deltas, intensities, strict=False)) / total_intensity
            if total_intensity > 0 else mean(deltas)
        )
        strongest = min(group, key=_strongest_key)
        closest = min(group, key=lambda row: (
            row["Absolute_Delta_Mass"], -row["Intensity"],
            row["Observed_Mass"], row["Source_Row_Index"],
        ))
        relations = [relation for member in group for relation in incident[member["Audit_Row_ID"]]]
        unique_relations = {(relation["Peak_Row_A"], relation["Peak_Row_B"]): relation for relation in relations}.values()
        integer_count = sum(bool(row["Integer_Spacing_Candidate"]) for row in unique_relations)
        isotope_count = sum(bool(row["Isotope_Spacing_Candidate"]) for row in unique_relations)
        recurrent_count = sum(bool(row["Recurrent_Spacing_Group_ID"]) for row in unique_relations)
        duplicate_infos = {
            info["Duplicate_Group_ID"]: info
            for member in group
            if (info := duplicate_by_row.get(member["Audit_Row_ID"])) is not None
        }
        duplicate_like = bool(duplicate_infos)
        integer_series = integer_count > 0
        isotope_candidate = isotope_count > 0
        recurrent_candidate = recurrent_count > 0
        relation_flags = sum((duplicate_like, integer_series, isotope_candidate, recurrent_candidate))
        span = max(deltas) - min(deltas)
        if relation_flags > 1:
            label = "MIXED_RELATION"
        elif duplicate_like:
            label = "DUPLICATE_LIKE"
        elif isotope_candidate:
            label = "ISOTOPE_SPACING_CANDIDATE"
        elif integer_series:
            label = "INTEGER_SPACING_SERIES"
        elif recurrent_candidate:
            label = "RECURRENT_SPACING_CANDIDATE"
        elif len(group) == 1:
            label = "SINGLETON"
        elif span <= params.duplicate_apex_tolerance_da:
            label = "TIGHT_MULTI_PEAK"
        else:
            label = "BROAD_MULTI_PEAK"
        duplicate_ids = sorted(duplicate_infos)
        duplicate_info_values = [duplicate_infos[key] for key in duplicate_ids]
        cluster_rows.append({
            "Cluster_ID": cluster_id,
            "Cluster_Label": label,
            "Cluster_Size": len(group),
            "Cluster_Min_Delta_Da": min(deltas),
            "Cluster_Max_Delta_Da": max(deltas),
            "Cluster_Span_Da": span,
            "Cluster_Mean_Delta_Da": mean(deltas),
            "Cluster_Median_Delta_Da": median(deltas),
            "Cluster_Weighted_Mean_Delta_Da": weighted_mean,
            "Cluster_Weighted_Mean_Fallback_Used": weighted_fallback,
            "Cluster_Total_Intensity": total_intensity,
            "Cluster_Max_Intensity": max(intensities),
            "Cluster_Strongest_Observed_Mass": strongest["Observed_Mass"],
            "Cluster_Closest_To_Theoretical_Mass": closest["Observed_Mass"],
            "Cluster_Closest_Absolute_Delta_Da": closest["Absolute_Delta_Mass"],
            "Cluster_Detection_Modes": ";".join(sorted({row["Detection_Mode"] for row in group})),
            "Cluster_Comparison_Statuses": ";".join(sorted({row["Comparison_Status"] for row in group})),
            "Cluster_Is_Singleton": len(group) == 1,
            "Cluster_Is_Multi_Peak": len(group) >= params.minimum_cluster_size,
            "Duplicate_Group_ID": ";".join(duplicate_ids),
            "Duplicate_Group_Size": max((info["Duplicate_Group_Size"] for info in duplicate_info_values), default=0),
            "Duplicate_Like": duplicate_like,
            "Duplicate_Mass_Span_Da": max((info["Duplicate_Mass_Span_Da"] for info in duplicate_info_values), default=0.0),
            "Duplicate_Strongest_Row": ";".join(str(info["Duplicate_Strongest_Row"]) for info in duplicate_info_values),
            "Duplicate_Representative_Row": ";".join(str(info["Duplicate_Representative_Row"]) for info in duplicate_info_values),
            "Integer_Spacing_Series": integer_series,
            "Isotope_Spacing_Candidate": isotope_candidate,
            "Recurrent_Spacing_Candidate": recurrent_candidate,
            "Integer_Spacing_Relation_Count": integer_count,
            "Isotope_Spacing_Relation_Count": isotope_count,
            "Recurrent_Spacing_Relation_Count": recurrent_count,
            "Member_Row_Indices": ";".join(str(row["Source_Row_Index"]) for row in group),
            "Member_Observed_Masses": _format_numbers([row["Observed_Mass"] for row in group]),
            "Member_Delta_Masses": _format_numbers(deltas),
            "Algorithm_Version": ALGORITHM_VERSION,
            **FORMAL_FALSE,
            "Notes": "Complete-link span-bounded numerical delta-mass cluster; no chemical identity assigned.",
        })

    largest = min(cluster_rows, key=lambda row: (-row["Cluster_Size"], row["Cluster_Span_Da"], row["Cluster_ID"]))
    strongest_cluster = min(cluster_rows, key=lambda row: (-row["Cluster_Total_Intensity"], -row["Cluster_Max_Intensity"], row["Cluster_ID"]))
    closest_cluster = min(cluster_rows, key=lambda row: (row["Cluster_Closest_Absolute_Delta_Da"], -row["Cluster_Max_Intensity"], row["Cluster_ID"]))
    closest_peak = min(rows, key=lambda row: (row["Absolute_Delta_Mass"], -row["Intensity"], row["Observed_Mass"], row["Source_Row_Index"]))
    strongest_peak = min(rows, key=_strongest_key)
    identity_status = str(comparison_summary.get("Input_Identity_Audit_Status") or "NOT_RUN")
    identity_conflict = bool(comparison_summary.get("Input_Identity_Conflict", False))
    biological_eligible = bool(comparison_summary.get("Biological_Interpretation_Eligible", False))
    notes = "Delta-mass clusters are numerical diagnostics only; no modification, isotope, adduct, charge, or molecular identity is assigned."
    if identity_conflict or identity_status in {"INSUFFICIENT_INFORMATION", "NOT_ELIGIBLE", "NOT_RUN"}:
        notes += " SCIEX input identity conflicts with or cannot be confirmed against the configured target; results must not be interpreted biologically."
    summary = {
        "Audit_Status": "AUDIT_COMPLETED",
        "Audit_Eligible": True,
        "SCIEX_Source_File": str(comparison_summary.get("Source_File") or detail_rows[0].get("Source_File") or ""),
        "Theoretical_Unmodified_Mass": theoretical_mass,
        "Detected_Peak_Count": len(rows),
        "Cluster_Count": len(cluster_rows),
        "Singleton_Cluster_Count": sum(row["Cluster_Is_Singleton"] for row in cluster_rows),
        "Multi_Peak_Cluster_Count": sum(row["Cluster_Is_Multi_Peak"] for row in cluster_rows),
        "Duplicate_Like_Cluster_Count": sum(row["Duplicate_Like"] for row in cluster_rows),
        "Integer_Spacing_Cluster_Count": sum(row["Integer_Spacing_Series"] for row in cluster_rows),
        "Isotope_Spacing_Candidate_Count": sum(row["Isotope_Spacing_Candidate"] for row in cluster_rows),
        "Recurrent_Spacing_Group_Count": len(recurrent_groups),
        "Recurrent_Spacing_Cluster_Count": sum(row["Recurrent_Spacing_Candidate"] for row in cluster_rows),
        "Largest_Cluster_Size": largest["Cluster_Size"],
        "Largest_Cluster_ID": largest["Cluster_ID"],
        "Strongest_Cluster_ID": strongest_cluster["Cluster_ID"],
        "Closest_Cluster_ID": closest_cluster["Cluster_ID"],
        "Closest_Cluster_Median_Delta_Da": closest_cluster["Cluster_Median_Delta_Da"],
        "Closest_Peak_Row": closest_peak["Source_Row_Index"],
        "Closest_Peak_Cluster_ID": cluster_by_audit_row[closest_peak["Audit_Row_ID"]],
        "Strongest_Peak_Row": strongest_peak["Source_Row_Index"],
        "Strongest_Peak_Cluster_ID": cluster_by_audit_row[strongest_peak["Audit_Row_ID"]],
        "Pair_Analysis_Performed": True,
        "Total_Eligible_Pair_Count": len(eligible_pairs),
        "Total_Relevant_Pair_Count": len(relevant_relations),
        "Exported_Pair_Count": len(relation_rows),
        "Pair_Rows_Truncated": pair_rows_truncated,
        "Input_Identity_Audit_Status": identity_status,
        "Input_Identity_Conflict": identity_conflict,
        "Biological_Interpretation_Eligible": biological_eligible,
        "Cluster_Tolerance_Da": params.cluster_tolerance_da,
        "Duplicate_Apex_Tolerance_Da": params.duplicate_apex_tolerance_da,
        "Isotope_Spacing_Da": params.isotope_spacing_da,
        "Isotope_Spacing_Tolerance_Da": params.isotope_spacing_tolerance_da,
        "Integer_Spacing_Tolerance_Da": params.integer_spacing_tolerance_da,
        "Minimum_Cluster_Size": params.minimum_cluster_size,
        "Max_Pair_Spacing_Da": params.max_pair_spacing_da,
        "Max_Pair_Rows": params.max_pair_rows,
        "Algorithm_Version": ALGORITHM_VERSION,
        **FORMAL_FALSE,
        "Notes": notes,
    }
    return SciexDeltaMassClusterAuditResult(
        tuple(cluster_rows), (summary,), tuple(relation_rows),
    )
