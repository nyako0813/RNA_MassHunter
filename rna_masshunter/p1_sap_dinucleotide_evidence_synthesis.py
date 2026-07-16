"""Shadow-only evidence synthesis for existing P1+SAP dinucleotide audit rows.

This module is deliberately a pure table integration layer.  It does not read
mzML, extract peaks, consume configuration, or mutate formal results.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Any

from rna_masshunter.p1_sap_dinucleotide_candidates import FORMAL_FALSE


EVIDENCE_COLUMNS = [
    "Physical_Feature_ID", "Dinucleotide_Feature_Count", "Dinucleotide_Feature_IDs",
    "Candidate_Group_Count", "Candidate_Group_IDs", "Qualified_Candidate_Group_Count",
    "Evidence_Level", "Evidence_Basis", "Mass_Accuracy_Status",
    "Chromatographic_Quality_Status", "Isotope_Evidence_Status", "MS2_Provenance_Status",
    "Composition_Resolution_Status", "Linkage_Resolution_Status",
    "Structure_Resolution_Status", "Source_Bond_Resolution_Status",
    "Unresolved_Issue_Count", "Unresolved_Issues",
    "Targeted_MS2_Priority", "Targeted_MS2_Rationale", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]

GROUP_EVIDENCE_COLUMNS = [
    "Dinucleotide_Group_ID", "Physical_Feature_Count", "Qualified_Physical_Feature_Count",
    "Best_Evidence_Level", "Group_Evidence_Status", "Supporting_Physical_Feature_IDs",
    "Structural_Assignment_Count", "Candidate_Specific_Feature_Count",
    "Linkage_Specific_Feature_Count", "Composition_Specific_Feature_Count",
    "Isotope_Compatible_Feature_Count", "Isotope_Confounded_Feature_Count",
    "Precursor_Compatible_MS2_Spectrum_Count", "Unresolved_Issue_Count",
    "Unresolved_Issues", "Targeted_MS2_Priority", "Targeted_MS2_Rationale",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]

SUMMARY_COLUMNS = [
    "Physical_Feature_Count", "Group_Count", "Group_With_Feature_Count",
    "Qualified_Physical_Feature_Count", "Candidate_Specific_Physical_Feature_Count",
    "Competition_Unresolved_Physical_Feature_Count", "Isotope_Compatible_Physical_Feature_Count",
    "Precursor_Compatible_MS2_Physical_Feature_Count", "High_Targeted_MS2_Priority_Count",
    "Medium_Targeted_MS2_Priority_Count", "Low_Targeted_MS2_Priority_Count",
    "Not_Applicable_Targeted_MS2_Priority_Count", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]

_EVIDENCE_RANK = {
    "NO_USABLE_EVIDENCE": 0,
    "RAW_MATCH_ONLY": 1,
    "CHROMATOGRAPHIC_FEATURE": 2,
    "QUALIFIED_MS1_FEATURE": 3,
    "QUALIFIED_MS1_WITH_ISOTOPE_SUPPORT": 4,
    "QUALIFIED_MS1_WITH_MS2_PROVENANCE": 5,
    "QUALIFIED_MULTI_AXIS_EVIDENCE": 6,
}
_PRIORITY_RANK = {"NOT_APPLICABLE": 0, "LOW": 1, "MEDIUM": 2, "HIGH": 3}


@dataclass(frozen=True)
class DinucleotideEvidenceSynthesisResult:
    evidence_rows: list[dict[str, Any]]
    group_evidence_rows: list[dict[str, Any]]
    summary_rows: list[dict[str, Any]]

    @property
    def sheets(self) -> dict[str, list[dict[str, Any]]]:
        return {
            "P1_SAP_Dinuc_Evidence": self.evidence_rows,
            "P1_SAP_Dinuc_Group_Evidence": self.group_evidence_rows,
            "P1_SAP_Dinuc_Evidence_Summary": self.summary_rows,
        }

    def to_jsonable(self) -> dict[str, Any]:
        return {
            "physical_features": self.evidence_rows,
            "groups": self.group_evidence_rows,
            "summary": self.summary_rows[0] if self.summary_rows else {},
        }


def _unique(values: list[Any]) -> list[str]:
    return sorted({str(value) for value in values if value not in (None, "")})


def _joined(values: list[Any]) -> str:
    return ";".join(_unique(values))


def _best(values: list[str], ranking: dict[str, int]) -> str:
    return max(values, key=lambda value: ranking.get(value, -1)) if values else ""


def _is_qualified(row: dict[str, Any]) -> bool:
    return bool(row.get("Feature_Eligible_For_Support"))


def _has_compatible_ms2(rows: list[dict[str, Any]]) -> bool:
    return any(int(row.get("Precursor_Compatible_MS2_Spectrum_Count") or 0) > 0 for row in rows)


def _isotope_informative(row: dict[str, Any]) -> bool:
    noninformative = {
        "ENVELOPE_CONFOUNDED", "ENVELOPE_TOO_WEAK", "NOT_ASSESSED", "INSUFFICIENT_DATA",
    }
    return (
        str(row.get("Envelope_Status") or "NOT_ASSESSED") not in noninformative
        and not bool(row.get("Envelope_Confounded"))
        and not bool(row.get("Isotope_Peak_Shared"))
        and not bool(row.get("Isomer_Isotope_Indistinguishable"))
    )


def _has_informative_isotope_support(rows: list[dict[str, Any]]) -> bool:
    return any(
        _isotope_informative(row) and row.get("Envelope_Status") == "ENVELOPE_COMPATIBLE"
        for row in rows
    )


def _physical_evidence_level(
    feature_rows: list[dict[str, Any]], isotope_rows: list[dict[str, Any]],
    ms2_rows: list[dict[str, Any]],
) -> str:
    if not feature_rows:
        return "NO_USABLE_EVIDENCE"
    if not any(_is_qualified(row) for row in feature_rows):
        return "CHROMATOGRAPHIC_FEATURE"
    isotope = _has_informative_isotope_support(isotope_rows)
    ms2 = _has_compatible_ms2(ms2_rows)
    if isotope and ms2:
        return "QUALIFIED_MULTI_AXIS_EVIDENCE"
    if ms2:
        return "QUALIFIED_MS1_WITH_MS2_PROVENANCE"
    if isotope:
        return "QUALIFIED_MS1_WITH_ISOTOPE_SUPPORT"
    return "QUALIFIED_MS1_FEATURE"


def _competition_types(rows: list[dict[str, Any]]) -> set[str]:
    return {
        value
        for row in rows
        for value in str(row.get("Competition_Types") or "").split(";")
        if value
    }


def _status_is_resolved(value: Any) -> bool:
    normalized = str(value or "").strip().upper()
    return bool(normalized) and "RESOLVED" in normalized and "UNRESOLVED" not in normalized


def _source_bond_candidates(
    group_ids: list[str], assignments_by_group: dict[str, list[dict[str, Any]]],
) -> list[str]:
    return _unique([
        row.get("Possible_Source_Bond")
        for group_id in group_ids
        for row in assignments_by_group.get(group_id, [])
    ])


def _source_bond_resolved(
    group_rows: list[dict[str, Any]], feature_rows: list[dict[str, Any]],
    assignment_rows: list[dict[str, Any]], source_bonds: list[str],
) -> bool:
    existing_rows = [*group_rows, *feature_rows, *assignment_rows]
    existing_resolved = any(
        bool(row.get("Original_Bond_Localized"))
        or _status_is_resolved(row.get("Source_Bond_Resolution_Status"))
        or _status_is_resolved(row.get("Localization_Status"))
        for row in existing_rows
    )
    return len(source_bonds) == 1 and existing_resolved


def _physical_unresolved(
    feature_rows: list[dict[str, Any]], competition_rows: list[dict[str, Any]],
    isotope_rows: list[dict[str, Any]], ms2_rows: list[dict[str, Any]],
    composition_resolved: bool, linkage_resolved: bool, structure_resolved: bool,
    source_bond_resolved: bool,
) -> list[str]:
    issues: set[str] = set()
    statuses = {str(row.get("Feature_Quality_Status") or "") for row in feature_rows}
    if any(value in statuses for value in {"SINGLE_POINT_REJECTED", "PROFILE_ONLY_REJECTED"}):
        issues.add("INSUFFICIENT_PROFILE_SUPPORT")
    if "BACKGROUND_TRACE_REJECTED" in statuses:
        issues.add("PERSISTENT_BACKGROUND")
    if "MASS_INCOMPATIBLE" in statuses:
        issues.add("MASS_ACCURACY_INCOMPATIBLE")

    informative_isotopes = [row for row in isotope_rows if _isotope_informative(row)]
    isotope_statuses = {str(row.get("Envelope_Status") or "NOT_ASSESSED") for row in isotope_rows}
    if not informative_isotopes:
        issues.add("ISOTOPE_NOT_INFORMATIVE")
    if "ENVELOPE_INCOMPATIBLE" in {str(row.get("Envelope_Status")) for row in informative_isotopes}:
        issues.add("ISOTOPE_INCOMPATIBLE")
    if isotope_rows and not informative_isotopes and (
        "ENVELOPE_CONFOUNDED" in isotope_statuses
        or any(
            row.get("Envelope_Confounded")
            or row.get("Isotope_Peak_Shared")
            or row.get("Isomer_Isotope_Indistinguishable")
            for row in isotope_rows
        )
    ):
        issues.add("ISOTOPE_CONFOUNDED")

    types = _competition_types(competition_rows)
    if "DIFFERENT_COMPOSITION_WITHIN_TOLERANCE" in types or not composition_resolved:
        issues.add("COMPOSITION_COMPETITION")
    if "NORMAL_PHOSPHATE_VS_PT_COMPETITION" in types or not linkage_resolved:
        issues.add("LINKAGE_COMPETITION")
    if "SAME_COMPOSITION_STRUCTURAL_ISOMERS" in types or not structure_resolved:
        issues.add("STRUCTURAL_ISOMERISM")
    if not source_bond_resolved:
        issues.add("SOURCE_BOND_UNRESOLVED")

    if not _has_compatible_ms2(ms2_rows):
        issues.add("NO_PRECURSOR_COMPATIBLE_MS2")
    elif any(not bool(row.get("MS2_Model_Applicable")) for row in ms2_rows):
        issues.add("MS2_FRAGMENT_MODEL_NOT_VALIDATED")
    return sorted(issues)


def _targeted_priority(
    qualified: bool, candidate_specific: bool, ambiguity: bool,
    compatible_ms2: bool, fragment_model_unvalidated: bool,
) -> tuple[str, str]:
    if not qualified:
        return "NOT_APPLICABLE", "No qualified usable feature is available."
    if ambiguity and compatible_ms2 and fragment_model_unvalidated:
        return "MEDIUM", "Compatible MS2 exists, but the fragment model is unvalidated and ambiguity remains."
    if ambiguity or not compatible_ms2:
        return "HIGH", "Qualified evidence retains ambiguity or lacks compatible MS2."
    if candidate_specific:
        return "LOW", "Candidate-specific evidence has no major unresolved ambiguity."
    return "MEDIUM", "Qualified evidence is not candidate-specific despite no classified major ambiguity."

def build_p1_sap_dinucleotide_evidence_synthesis(
    groups: list[dict[str, Any]], assignments: list[dict[str, Any]],
    features: list[dict[str, Any]], competition: list[dict[str, Any]],
    isotopes: list[dict[str, Any]], ms2_provenance: list[dict[str, Any]],
) -> DinucleotideEvidenceSynthesisResult:
    """Integrate existing dinucleotide tables without reading source data."""
    group_map = {str(row.get("Dinucleotide_Group_ID")): row for row in groups}
    assignments_by_group: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in assignments:
        assignments_by_group[str(row.get("Dinucleotide_Group_ID"))].append(row)
    features_by_physical: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in features:
        features_by_physical[str(row.get("Physical_Feature_ID"))].append(row)
    competition_by_key = {
        (str(row.get("Physical_Feature_ID")), str(row.get("Dinucleotide_Group_ID"))): row
        for row in competition
    }
    isotopes_by_feature: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in isotopes:
        isotopes_by_feature[str(row.get("Dinucleotide_Feature_ID"))].append(row)
    ms2_by_feature: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in ms2_provenance:
        ms2_by_feature[str(row.get("Dinucleotide_Feature_ID"))].append(row)

    evidence_rows: list[dict[str, Any]] = []
    for physical_id, physical_features in sorted(features_by_physical.items()):
        feature_ids = _unique([row.get("Dinucleotide_Feature_ID") for row in physical_features])
        group_ids = _unique([row.get("Dinucleotide_Group_ID") for row in physical_features])
        physical_isotopes = [row for feature_id in feature_ids for row in isotopes_by_feature.get(feature_id, [])]
        physical_ms2 = [row for feature_id in feature_ids for row in ms2_by_feature.get(feature_id, [])]
        physical_competition = [competition_by_key[key] for key in ((physical_id, group_id) for group_id in group_ids) if key in competition_by_key]
        physical_groups = [group_map[group_id] for group_id in group_ids if group_id in group_map]
        physical_assignments = [
            row for group_id in group_ids for row in assignments_by_group.get(group_id, [])
        ]
        assignment_count = len(physical_assignments)
        qualified = [row for row in physical_features if _is_qualified(row)]

        compositions = _unique([
            row.get("Final_Elemental_Composition") or row.get("Elemental_Composition")
            for row in physical_groups
        ])
        linkages = _unique([row.get("Linkage_State") for row in physical_groups])
        composition_resolved = len(compositions) == 1
        linkage_resolved = len(linkages) == 1

        specificity_rows = physical_competition or physical_features
        candidate_specific = (
            len(group_ids) == 1
            and bool(specificity_rows)
            and all(bool(row.get("Candidate_Specific")) for row in specificity_rows)
        )
        existing_structure_specific = (
            bool(physical_features)
            and all(bool(row.get("Structure_Specific")) for row in physical_features)
            and (
                not physical_competition
                or all(bool(row.get("Structure_Specific")) for row in physical_competition)
            )
        )
        structure_resolved = candidate_specific and assignment_count == 1 and existing_structure_specific
        source_bonds = _source_bond_candidates(group_ids, assignments_by_group)
        source_bond_resolved = _source_bond_resolved(
            physical_groups, physical_features, physical_assignments, source_bonds,
        )

        evidence_level = _physical_evidence_level(
            physical_features, physical_isotopes, physical_ms2,
        )
        issues = _physical_unresolved(
            physical_features, physical_competition, physical_isotopes, physical_ms2,
            composition_resolved, linkage_resolved, structure_resolved,
            source_bond_resolved,
        )
        compatible_ms2 = _has_compatible_ms2(physical_ms2)
        fragment_model_unvalidated = compatible_ms2 and any(
            not bool(row.get("MS2_Model_Applicable")) for row in physical_ms2
        )
        ambiguity = bool(set(issues) & {
            "COMPOSITION_COMPETITION", "LINKAGE_COMPETITION",
            "STRUCTURAL_ISOMERISM", "SOURCE_BOND_UNRESOLVED",
        })
        priority, rationale = _targeted_priority(
            bool(qualified), candidate_specific, ambiguity,
            compatible_ms2, fragment_model_unvalidated,
        )
        isotope_statuses = [
            str(row.get("Envelope_Status") or "NOT_ASSESSED")
            for row in physical_isotopes
        ]
        informative_isotope = _has_informative_isotope_support(physical_isotopes)
        ms2_count = sum(
            int(row.get("Precursor_Compatible_MS2_Spectrum_Count") or 0)
            for row in physical_ms2
        )
        evidence_rows.append({
            "Physical_Feature_ID": physical_id,
            "Dinucleotide_Feature_Count": len(feature_ids),
            "Dinucleotide_Feature_IDs": ";".join(feature_ids),
            "Candidate_Group_Count": len(group_ids),
            "Candidate_Group_IDs": ";".join(group_ids),
            "Qualified_Candidate_Group_Count": len({
                str(row.get("Dinucleotide_Group_ID")) for row in qualified
            }),
            "Evidence_Level": evidence_level,
            "Evidence_Basis": ";".join(filter(None, [
                "QUALIFIED_MS1" if qualified else "MS1_FEATURE_ONLY",
                "ISOTOPE_COMPATIBLE" if informative_isotope else "",
                "PRECURSOR_COMPATIBLE_MS2" if ms2_count else "",
            ])),
            "Mass_Accuracy_Status": _joined([
                row.get("Mass_Accuracy_Support_Status") for row in physical_features
            ]),
            "Chromatographic_Quality_Status": _joined([
                row.get("Feature_Quality_Status") for row in physical_features
            ]),
            "Isotope_Evidence_Status": (
                _joined(isotope_statuses) if any(
                    _isotope_informative(row) for row in physical_isotopes
                ) else "NONINFORMATIVE"
            ),
            "MS2_Provenance_Status": (
                "PRECURSOR_COMPATIBLE_MS2_PRESENT"
                if ms2_count else "NO_PRECURSOR_COMPATIBLE_MS2"
            ),
            "Composition_Resolution_Status": (
                "COMPOSITION_RESOLVED"
                if composition_resolved else "COMPOSITION_UNRESOLVED"
            ),
            "Linkage_Resolution_Status": (
                "LINKAGE_RESOLVED" if linkage_resolved else "LINKAGE_UNRESOLVED"
            ),
            "Structure_Resolution_Status": (
                "STRUCTURE_RESOLVED" if structure_resolved else "STRUCTURE_UNRESOLVED"
            ),
            "Source_Bond_Resolution_Status": (
                "SOURCE_BOND_RESOLVED"
                if source_bond_resolved else "SOURCE_BOND_UNRESOLVED"
            ),
            "Unresolved_Issue_Count": len(issues),
            "Unresolved_Issues": ";".join(issues),
            "Targeted_MS2_Priority": priority,
            "Targeted_MS2_Rationale": rationale,
            **FORMAL_FALSE,
        })

    evidence_by_group: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in evidence_rows:
        for group_id in str(row.get("Candidate_Group_IDs") or "").split(";"):
            if group_id:
                evidence_by_group[group_id].append(row)
    feature_by_group: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in features:
        feature_by_group[str(row.get("Dinucleotide_Group_ID"))].append(row)
    group_evidence_rows: list[dict[str, Any]] = []
    for group_id in sorted(group_map):
        group = group_map[group_id]
        group_features = feature_by_group.get(group_id, [])
        physical_rows = evidence_by_group.get(group_id, [])
        qualified_physical_ids = {str(row.get("Physical_Feature_ID")) for row in group_features if row.get("Feature_Eligible_For_Support")}
        feature_ids = _unique([row.get("Dinucleotide_Feature_ID") for row in group_features])
        group_isotopes = [row for feature_id in feature_ids for row in isotopes_by_feature.get(feature_id, [])]
        group_ms2 = [row for feature_id in feature_ids for row in ms2_by_feature.get(feature_id, [])]
        issues = _unique([issue for row in physical_rows for issue in str(row.get("Unresolved_Issues") or "").split(";")])
        if not physical_rows:
            issues = ["NO_QUALIFIED_PHYSICAL_FEATURE"]
        priority = _best([str(row.get("Targeted_MS2_Priority")) for row in physical_rows], _PRIORITY_RANK) or "NOT_APPLICABLE"
        priority_rows = [row for row in physical_rows if row.get("Targeted_MS2_Priority") == priority]
        best_level = _best(
            [str(row.get("Evidence_Level")) for row in physical_rows], _EVIDENCE_RANK,
        )
        if not best_level:
            best_level = (
                "RAW_MATCH_ONLY"
                if str(group.get("Group_Interpretation")) == "RAW_MATCH_ONLY"
                else "NO_USABLE_EVIDENCE"
            )
        group_evidence_rows.append({
            "Dinucleotide_Group_ID": group_id, "Physical_Feature_Count": len(physical_rows),
            "Qualified_Physical_Feature_Count": len(qualified_physical_ids), "Best_Evidence_Level": best_level,
            "Group_Evidence_Status": str(group.get("Group_Interpretation") or "NOT_EVALUABLE"),
            "Supporting_Physical_Feature_IDs": _joined([row.get("Physical_Feature_ID") for row in physical_rows]),
            "Structural_Assignment_Count": len(assignments_by_group.get(group_id, [])) or int(group.get("Structural_Assignment_Count") or 0),
            "Candidate_Specific_Feature_Count": sum(bool(row.get("Candidate_Specific")) for row in group_features),
            "Linkage_Specific_Feature_Count": sum(bool(row.get("Linkage_Specific")) for row in group_features),
            "Composition_Specific_Feature_Count": sum(bool(row.get("Composition_Specific")) for row in group_features),
            "Isotope_Compatible_Feature_Count": sum(
                row.get("Envelope_Status") == "ENVELOPE_COMPATIBLE"
                and _isotope_informative(row)
                for row in group_isotopes
            ),
            "Isotope_Confounded_Feature_Count": sum(
                not _isotope_informative(row) for row in group_isotopes
            ),
            "Precursor_Compatible_MS2_Spectrum_Count": sum(int(row.get("Precursor_Compatible_MS2_Spectrum_Count") or 0) for row in group_ms2),
            "Unresolved_Issue_Count": len(issues), "Unresolved_Issues": ";".join(issues),
            "Targeted_MS2_Priority": priority,
            "Targeted_MS2_Rationale": _joined([row.get("Targeted_MS2_Rationale") for row in priority_rows]) or "No qualified physical feature available for targeted MS2 prioritization.",
            **FORMAL_FALSE,
        })

    summary = {
        "Physical_Feature_Count": len(evidence_rows), "Group_Count": len(group_evidence_rows),
        "Group_With_Feature_Count": sum(int(row.get("Physical_Feature_Count") or 0) > 0 for row in group_evidence_rows),
        "Qualified_Physical_Feature_Count": sum(int(row.get("Qualified_Candidate_Group_Count") or 0) > 0 for row in evidence_rows),
        "Candidate_Specific_Physical_Feature_Count": sum(
            int(row.get("Candidate_Group_Count") or 0) == 1
            and row.get("Structure_Resolution_Status") == "STRUCTURE_RESOLVED"
            for row in evidence_rows
        ),
        "Competition_Unresolved_Physical_Feature_Count": sum(
            any(
                row.get(column) == unresolved
                for column, unresolved in (
                    ("Composition_Resolution_Status", "COMPOSITION_UNRESOLVED"),
                    ("Linkage_Resolution_Status", "LINKAGE_UNRESOLVED"),
                    ("Structure_Resolution_Status", "STRUCTURE_UNRESOLVED"),
                )
            )
            for row in evidence_rows
        ),
        "Isotope_Compatible_Physical_Feature_Count": sum(
            "ISOTOPE_COMPATIBLE" in str(row.get("Evidence_Basis"))
            for row in evidence_rows
        ),
        "Precursor_Compatible_MS2_Physical_Feature_Count": sum(row.get("MS2_Provenance_Status") == "PRECURSOR_COMPATIBLE_MS2_PRESENT" for row in evidence_rows),
        "High_Targeted_MS2_Priority_Count": sum(row.get("Targeted_MS2_Priority") == "HIGH" for row in evidence_rows),
        "Medium_Targeted_MS2_Priority_Count": sum(row.get("Targeted_MS2_Priority") == "MEDIUM" for row in evidence_rows),
        "Low_Targeted_MS2_Priority_Count": sum(row.get("Targeted_MS2_Priority") == "LOW" for row in evidence_rows),
        "Not_Applicable_Targeted_MS2_Priority_Count": sum(row.get("Targeted_MS2_Priority") == "NOT_APPLICABLE" for row in evidence_rows),
        **FORMAL_FALSE,
    }
    return DinucleotideEvidenceSynthesisResult(evidence_rows, group_evidence_rows, [summary])
