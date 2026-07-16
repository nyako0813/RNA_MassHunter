"""Final shadow-only consensus across standard and complete-structure RNase MS/MS evidence."""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Any

FORMAL_FALSE = {
    "Applied_To_Formal_Result": False,
    "Formal_Change_Ready": False,
    "Formal_Result_Changed": False,
}

EVIDENCE_COLUMNS = [
    "Modification_ID", "Parent_Fragment_ID", "Candidate_Position_In_Parent",
    "Candidate_tRNA_Position", "Candidate_IDs", "Complete_Structure_IDs",
    "Exact_Crosswalk_Count", "Crosswalk_Cardinality", "Standard_Identity_Status",
    "Standard_Localization_Status", "Standard_Ambiguity_Status",
    "Composite_Identity_Status", "Composite_Localization_Status",
    "Composite_Backbone_Status", "Composite_Structure_Status",
    "Composite_Ambiguity_Status", "Standard_Composite_Consistency_Status",
    "Modification_Identity_Consensus", "Localization_Consensus",
    "Backbone_Consensus", "Structure_Consensus", "Consensus_Ambiguity_Status",
    "Review_Priority", "Consensus_Basis", "Limiting_Reasons",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]

SUMMARY_COLUMNS = [
    "Consensus_Row_Count", "Standard_Linked_Row_Count", "Standard_Only_Row_Count",
    "Composite_Only_Row_Count", "Exact_Crosswalk_Row_Count", "CONSISTENT_Count",
    "COMPOSITE_ADDS_CONTEXT_Count", "POSITION_CONFLICT_Count",
    "IDENTITY_CONFLICT_Count", "PARENT_FRAGMENT_CONFLICT_Count",
    "INSUFFICIENT_PROVENANCE_Count", "NOT_LINKED_Count",
    "Identity_SUPPORTED_Count", "Identity_PROVISIONAL_Count",
    "Identity_AMBIGUOUS_Count", "Identity_UNSUPPORTED_Count",
    "Localization_LOCALIZED_Count", "Localization_PARTIALLY_LOCALIZED_Count",
    "Localization_POSITION_COMPATIBLE_Count", "Localization_AMBIGUOUS_Count",
    "Localization_UNRESOLVED_Count", "Review_HIGH_Count", "Review_MEDIUM_Count",
    "Review_LOW_Count", "Review_NONE_Count", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]


@dataclass(frozen=True)
class RNaseMS2ConsensusSynthesisResult:
    summary_rows: list[dict[str, Any]]
    evidence_rows: list[dict[str, Any]]

    @property
    def sheets(self) -> dict[str, list[dict[str, Any]]]:
        return {
            "RNase_MS2_Consensus_Summary": self.summary_rows,
            "RNase_MS2_Consensus_Evidence": self.evidence_rows,
        }


def _text(value: Any) -> str:
    return "" if value is None else str(value)


def _standard_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    return tuple(_text(row.get(name)) for name in (
        "Modification_ID", "Parent_Fragment_ID", "Candidate_Position_In_Parent",
        "Candidate_tRNA_Position",
    ))


def _composite_key(row: dict[str, Any]) -> tuple[str, str]:
    return _text(row.get("Candidate_ID")), _text(row.get("Complete_Structure_ID"))


def _crosswalk_composite_key(row: dict[str, Any]) -> tuple[str, str]:
    return _text(row.get("Candidate_ID")), _text(row.get("Complete_Structure_ID"))


def _deduplicate(rows: list[dict[str, Any]], key_function) -> list[dict[str, Any]]:
    unique: dict[tuple[str, ...], dict[str, Any]] = {}
    for row in rows or []:
        unique.setdefault(key_function(row), dict(row))
    return [unique[key] for key in sorted(unique)]


def _joined(values) -> str:
    return ";".join(sorted({_text(value) for value in values if _text(value)}))


def _aggregate(rows: list[dict[str, Any]], field: str, empty: str = "") -> str:
    values = {_text(row.get(field)) for row in rows if _text(row.get(field))}
    if not values:
        return empty
    if len(values) == 1:
        return next(iter(values))
    return "AMBIGUOUS"


def _ambiguity(types: set[str]) -> str:
    return "NONE_OBSERVED" if not types else next(iter(types)) if len(types) == 1 else "MULTIPLE"


def _identity_consensus(standard: str, composites: list[dict[str, Any]], exact: bool, multiplicity: bool,
                        conflict: bool) -> str:
    if conflict or multiplicity:
        return "AMBIGUOUS"
    composite_values = {_text(row.get("Composite_Identity_Status")) for row in composites}
    if exact and standard == "FRAGMENT_SUPPORTED" and composite_values == {"PROVISIONAL_CANDIDATE_SUPPORT"}:
        return "SUPPORTED"
    if standard == "AMBIGUOUS" or "AMBIGUOUS" in composite_values:
        return "AMBIGUOUS"
    if standard == "FRAGMENT_SUPPORTED":
        return "PROVISIONAL" if not exact else "SUPPORTED"
    if standard == "PRECURSOR_COMPATIBLE" or "PROVISIONAL_CANDIDATE_SUPPORT" in composite_values:
        return "PROVISIONAL"
    return "UNSUPPORTED"


def _localization_consensus(standard: str, composites: list[dict[str, Any]], conflict: bool,
                            multiplicity: bool) -> str:
    values = {_text(row.get("Composite_Localization_Status")) for row in composites}
    if conflict or multiplicity or standard == "AMBIGUOUS" or "AMBIGUOUS" in values:
        return "AMBIGUOUS"
    if standard == "LOCALIZED" and values & {"POSITION_COMPATIBLE", "PARTIALLY_SUPPORTED"}:
        return "LOCALIZED"
    if standard == "LOCALIZED":
        return "LOCALIZED"
    if standard == "PARTIALLY_LOCALIZED" or "PARTIALLY_SUPPORTED" in values:
        return "PARTIALLY_LOCALIZED"
    if "POSITION_COMPATIBLE" in values:
        return "POSITION_COMPATIBLE"
    return "UNRESOLVED"


def _backbone_consensus(composites: list[dict[str, Any]]) -> str:
    values = {_text(row.get("Composite_Backbone_Status")) for row in composites}
    if "AMBIGUOUS" in values or len(values - {"", "NOT_EVALUATED"}) > 1:
        return "AMBIGUOUS"
    if "SUPPORTED" in values:
        return "SUPPORTED"
    if "PARTIALLY_SUPPORTED" in values:
        return "PARTIALLY_SUPPORTED"
    if "BACKBONE_COMPATIBLE" in values:
        return "BOND_COMPATIBLE"
    return "NOT_EVALUATED"


def _structure_consensus(standard_structure: str, composites: list[dict[str, Any]]) -> str:
    values = {_text(row.get("Composite_Structure_Status")) for row in composites}
    if "AMBIGUOUS" in values or standard_structure == "AMBIGUOUS":
        return "AMBIGUOUS"
    if values == {"SUPPORTED"} and standard_structure == "SUPPORTED":
        return "SUPPORTED"
    if values & {"UNRESOLVED", "SUPPORTED"}:
        return "UNRESOLVED"
    return "NOT_EVALUATED"


def _conflict_edge(rows: list[dict[str, Any]]) -> str:
    if any(row.get("Crosswalk_Status") == "MODIFICATION_IDENTITY_CONFLICT" and
           row.get("Position_Match_Status") == "MATCH" and
           row.get("Parent_Fragment_Match_Status") == "MATCH" for row in rows):
        return "IDENTITY_CONFLICT"
    if any(row.get("Crosswalk_Status") == "PARENT_FRAGMENT_CONFLICT" and
           row.get("Position_Match_Status") == "MATCH" for row in rows):
        return "PARENT_FRAGMENT_CONFLICT"
    if any(row.get("Crosswalk_Status") == "POSITION_CONFLICT" and
           row.get("Parent_Fragment_Match_Status") == "MATCH" for row in rows):
        return "POSITION_CONFLICT"
    if any(row.get("Crosswalk_Status") in {"INSUFFICIENT_PROVENANCE", "MASS_EQUIVALENT_ONLY"} for row in rows):
        return "INSUFFICIENT_PROVENANCE"
    return "NOT_LINKED"


def build_rnase_ms2_consensus_synthesis(
    standard_candidate_rows: list[dict[str, Any]],
    composite_candidate_rows: list[dict[str, Any]],
    crosswalk_rows: list[dict[str, Any]],
) -> RNaseMS2ConsensusSynthesisResult:
    """Build final conservative consensus without mutating any source table."""
    standards = _deduplicate(list(standard_candidate_rows or []), _standard_key)
    composites = _deduplicate(list(composite_candidate_rows or []), _composite_key)
    crosswalk = [dict(row) for row in crosswalk_rows or []]
    composite_by_key = {_composite_key(row): row for row in composites}
    crosswalk_by_standard: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in crosswalk:
        crosswalk_by_standard[_standard_key(row)].append(row)

    evidence: list[dict[str, Any]] = []
    exactly_linked_composites: set[tuple[str, str]] = set()
    for standard in standards:
        key = _standard_key(standard)
        edges = crosswalk_by_standard.get(key, [])
        exact_edges = [row for row in edges if row.get("Crosswalk_Status") == "EXACT_MATCH"]
        exact_keys = sorted({_crosswalk_composite_key(row) for row in exact_edges})
        linked_composites = [composite_by_key[item] for item in exact_keys if item in composite_by_key]
        exactly_linked_composites.update(exact_keys)
        cardinalities = {_text(row.get("Crosswalk_Cardinality")) for row in exact_edges if row.get("Crosswalk_Cardinality")}
        multiplicity = len(exact_keys) > 1 or bool(cardinalities & {"ONE_TO_MANY", "MANY_TO_ONE", "MANY_TO_MANY"})
        standard_identity = _text(standard.get("Modification_Identity_Status"))
        standard_localization = _text(standard.get("Localization_Status"))
        standard_structure = _text(standard.get("Structure_Status"))
        standard_ambiguity = _text(standard.get("Ambiguity_Status"))
        composite_identity = _aggregate(linked_composites, "Composite_Identity_Status")
        composite_localization = _aggregate(linked_composites, "Composite_Localization_Status")
        composite_backbone = _aggregate(linked_composites, "Composite_Backbone_Status")
        composite_structure = _aggregate(linked_composites, "Composite_Structure_Status")
        composite_ambiguity = _aggregate(linked_composites, "Composite_Ambiguity_Status")

        evidence_conflict = bool(exact_edges) and (
            (standard_identity == "FRAGMENT_SUPPORTED" and composite_identity in {"UNSUPPORTED", "AMBIGUOUS"})
            or (standard_localization == "LOCALIZED" and composite_localization == "AMBIGUOUS")
        )
        if exact_edges:
            if evidence_conflict and standard_identity == "FRAGMENT_SUPPORTED" and composite_identity in {"UNSUPPORTED", "AMBIGUOUS"}:
                consistency = "IDENTITY_CONFLICT"
            elif evidence_conflict:
                consistency = "POSITION_CONFLICT"
            elif composite_backbone not in {"", "NOT_EVALUATED", "UNRESOLVED"}:
                consistency = "COMPOSITE_ADDS_CONTEXT"
            else:
                consistency = "CONSISTENT"
        elif edges:
            consistency = _conflict_edge(edges)
        else:
            consistency = "NOT_LINKED"

        identity = _identity_consensus(standard_identity, linked_composites, bool(exact_edges), multiplicity, evidence_conflict or consistency == "IDENTITY_CONFLICT")
        localization = _localization_consensus(standard_localization, linked_composites, consistency == "POSITION_CONFLICT", multiplicity)
        backbone = _backbone_consensus(linked_composites)
        structure = _structure_consensus(standard_structure, linked_composites)
        ambiguity_types: set[str] = set()
        if not exact_edges: ambiguity_types.add("STANDARD_ONLY")
        if multiplicity: ambiguity_types.add("CROSSWALK_MULTIPLICITY")
        if evidence_conflict or consistency in {"POSITION_CONFLICT", "IDENTITY_CONFLICT", "PARENT_FRAGMENT_CONFLICT"}:
            ambiguity_types.add("EVIDENCE_CONFLICT")
        if standard_ambiguity not in {"", "NONE"} or composite_ambiguity not in {"", "NONE"}:
            ambiguity_types.add("EVIDENCE_CONFLICT")
        ambiguity = _ambiguity(ambiguity_types)
        if "EVIDENCE_CONFLICT" in ambiguity_types or (multiplicity and len(linked_composites) > 1) or structure == "AMBIGUOUS":
            priority = "HIGH"
        elif identity == "PROVISIONAL" or localization == "POSITION_COMPATIBLE" or backbone == "BOND_COMPATIBLE" or consistency == "INSUFFICIENT_PROVENANCE":
            priority = "MEDIUM"
        elif standard_identity in {"", "UNSUPPORTED"} and not linked_composites:
            priority = "NONE"
        else:
            priority = "LOW"
        reasons = []
        if not exact_edges: reasons.append("no_exact_crosswalk")
        if multiplicity: reasons.append("crosswalk_multiplicity_limits_consensus")
        if evidence_conflict: reasons.append("standard_composite_evidence_conflict")
        if composite_structure != "SUPPORTED": reasons.append("composite_structure_not_supported")
        if consistency == "INSUFFICIENT_PROVENANCE": reasons.append("crosswalk_provenance_insufficient")
        evidence.append({
            "Modification_ID": key[0], "Parent_Fragment_ID": key[1],
            "Candidate_Position_In_Parent": standard.get("Candidate_Position_In_Parent", ""),
            "Candidate_tRNA_Position": standard.get("Candidate_tRNA_Position", ""),
            "Candidate_IDs": _joined(row.get("Candidate_ID") for row in exact_edges),
            "Complete_Structure_IDs": _joined(row.get("Complete_Structure_ID") for row in exact_edges),
            "Exact_Crosswalk_Count": len(exact_edges),
            "Crosswalk_Cardinality": _joined(cardinalities),
            "Standard_Identity_Status": standard_identity,
            "Standard_Localization_Status": standard_localization,
            "Standard_Ambiguity_Status": standard_ambiguity,
            "Composite_Identity_Status": composite_identity,
            "Composite_Localization_Status": composite_localization,
            "Composite_Backbone_Status": composite_backbone,
            "Composite_Structure_Status": composite_structure,
            "Composite_Ambiguity_Status": composite_ambiguity,
            "Standard_Composite_Consistency_Status": consistency,
            "Modification_Identity_Consensus": identity,
            "Localization_Consensus": localization, "Backbone_Consensus": backbone,
            "Structure_Consensus": structure, "Consensus_Ambiguity_Status": ambiguity,
            "Review_Priority": priority,
            "Consensus_Basis": ";".join(filter(None, [
                f"standard_identity={standard_identity or 'missing'}",
                f"exact_crosswalks={len(exact_edges)}",
                f"composite_identity={composite_identity or 'none'}",
            ])),
            "Limiting_Reasons": ";".join(reasons), **FORMAL_FALSE,
        })

    for key in sorted(set(composite_by_key) - exactly_linked_composites):
        composite = composite_by_key[key]
        identity_value = _text(composite.get("Composite_Identity_Status"))
        identity = "PROVISIONAL" if identity_value == "PROVISIONAL_CANDIDATE_SUPPORT" else "AMBIGUOUS" if identity_value == "AMBIGUOUS" else "UNSUPPORTED"
        localization_value = _text(composite.get("Composite_Localization_Status"))
        localization = {"PARTIALLY_SUPPORTED": "PARTIALLY_LOCALIZED", "POSITION_COMPATIBLE": "POSITION_COMPATIBLE", "AMBIGUOUS": "AMBIGUOUS"}.get(localization_value, "UNRESOLVED")
        backbone = _backbone_consensus([composite])
        structure = _structure_consensus("", [composite])
        evidence.append({
            "Modification_ID": "", "Parent_Fragment_ID": "",
            "Candidate_Position_In_Parent": "", "Candidate_tRNA_Position": "",
            "Candidate_IDs": key[0], "Complete_Structure_IDs": key[1],
            "Exact_Crosswalk_Count": 0, "Crosswalk_Cardinality": "",
            "Standard_Identity_Status": "", "Standard_Localization_Status": "",
            "Standard_Ambiguity_Status": "", "Composite_Identity_Status": identity_value,
            "Composite_Localization_Status": localization_value,
            "Composite_Backbone_Status": composite.get("Composite_Backbone_Status", ""),
            "Composite_Structure_Status": composite.get("Composite_Structure_Status", ""),
            "Composite_Ambiguity_Status": composite.get("Composite_Ambiguity_Status", ""),
            "Standard_Composite_Consistency_Status": "NOT_LINKED",
            "Modification_Identity_Consensus": identity,
            "Localization_Consensus": localization, "Backbone_Consensus": backbone,
            "Structure_Consensus": structure, "Consensus_Ambiguity_Status": "COMPOSITE_ONLY",
            "Review_Priority": "MEDIUM" if identity == "PROVISIONAL" or backbone == "BOND_COMPATIBLE" else "NONE",
            "Consensus_Basis": "composite_candidate_without_exact_standard_crosswalk",
            "Limiting_Reasons": "no_exact_crosswalk;standard_candidate_unavailable", **FORMAL_FALSE,
        })

    evidence.sort(key=lambda row: (
        row["Modification_ID"], row["Parent_Fragment_ID"],
        _text(row["Candidate_Position_In_Parent"]), row["Candidate_IDs"],
    ))
    consistency_counts = Counter(row["Standard_Composite_Consistency_Status"] for row in evidence)
    identity_counts = Counter(row["Modification_Identity_Consensus"] for row in evidence)
    localization_counts = Counter(row["Localization_Consensus"] for row in evidence)
    priority_counts = Counter(row["Review_Priority"] for row in evidence)
    summary = [{
        "Consensus_Row_Count": len(evidence),
        "Standard_Linked_Row_Count": sum(bool(row["Exact_Crosswalk_Count"]) for row in evidence),
        "Standard_Only_Row_Count": sum(row["Consensus_Ambiguity_Status"] == "STANDARD_ONLY" for row in evidence),
        "Composite_Only_Row_Count": sum(row["Consensus_Ambiguity_Status"] == "COMPOSITE_ONLY" for row in evidence),
        "Exact_Crosswalk_Row_Count": sum(int(row["Exact_Crosswalk_Count"]) for row in evidence),
        **{f"{status}_Count": consistency_counts[status] for status in (
            "CONSISTENT", "COMPOSITE_ADDS_CONTEXT", "POSITION_CONFLICT", "IDENTITY_CONFLICT",
            "PARENT_FRAGMENT_CONFLICT", "INSUFFICIENT_PROVENANCE", "NOT_LINKED",
        )},
        **{f"Identity_{status}_Count": identity_counts[status] for status in ("SUPPORTED", "PROVISIONAL", "AMBIGUOUS", "UNSUPPORTED")},
        **{f"Localization_{status}_Count": localization_counts[status] for status in ("LOCALIZED", "PARTIALLY_LOCALIZED", "POSITION_COMPATIBLE", "AMBIGUOUS", "UNRESOLVED")},
        **{f"Review_{status}_Count": priority_counts[status] for status in ("HIGH", "MEDIUM", "LOW", "NONE")},
        **FORMAL_FALSE,
    }]
    return RNaseMS2ConsensusSynthesisResult(summary, evidence)
