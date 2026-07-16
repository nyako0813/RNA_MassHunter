"""Shadow-only crosswalk between standard RNase and complete-structure candidates."""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, is_dataclass
from typing import Any

FORMAL_FALSE = {
    "Applied_To_Formal_Result": False,
    "Formal_Change_Ready": False,
    "Formal_Result_Changed": False,
}

CROSSWALK_COLUMNS = [
    "Modification_ID", "Parent_Fragment_ID", "Candidate_Position_In_Parent",
    "Candidate_tRNA_Position", "Standard_Absolute_Sequence_Position",
    "Candidate_ID", "Complete_Structure_ID", "Composite_Position",
    "Applied_Transform_IDs", "Explicit_Legacy_Modification_IDs",
    "Mass_Equivalent_Modification_IDs", "Canonical_Structure_ID",
    "Position_Match_Status", "Modification_Identity_Match_Status",
    "Mass_Shift_Match_Status", "Parent_Fragment_Match_Status",
    "Structural_Isomer_Sharing", "Crosswalk_Status", "Crosswalk_Cardinality",
    "Crosswalk_Basis", "Limiting_Reasons", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]

SUMMARY_COLUMNS = [
    "Standard_Candidate_Count", "Composite_Position_Component_Count",
    "Crosswalk_Row_Count", "Related_Edge_Count", "EXACT_MATCH_Count",
    "POSITION_MATCH_IDENTITY_UNRESOLVED_Count", "MASS_EQUIVALENT_ONLY_Count",
    "POSITION_CONFLICT_Count", "MODIFICATION_IDENTITY_CONFLICT_Count",
    "PARENT_FRAGMENT_CONFLICT_Count", "NOT_MAPPABLE_Count",
    "INSUFFICIENT_PROVENANCE_Count", "ONE_TO_ONE_Count", "ONE_TO_MANY_Count",
    "MANY_TO_ONE_Count", "MANY_TO_MANY_Count", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]

RELATED_STATUSES = {
    "EXACT_MATCH", "POSITION_MATCH_IDENTITY_UNRESOLVED", "MASS_EQUIVALENT_ONLY",
}


@dataclass(frozen=True)
class RNaseMS2StandardCompositeCrosswalkResult:
    summary_rows: list[dict[str, Any]]
    crosswalk_rows: list[dict[str, Any]]

    @property
    def sheets(self) -> dict[str, list[dict[str, Any]]]:
        return {
            "RNase_MS2_Standard_Composite_Summary": self.summary_rows,
            "RNase_MS2_Standard_Composite_Crosswalk": self.crosswalk_rows,
        }


def _row(value: Any) -> dict[str, Any]:
    if isinstance(value, dict):
        return dict(value)
    if is_dataclass(value):
        return asdict(value)
    return dict(vars(value))


def _text(value: Any) -> str:
    return "" if value is None else str(value)


def _number(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _integer(value: Any) -> int | None:
    number = _number(value)
    return int(number) if number is not None and number.is_integer() else None


def _ids(value: Any) -> set[str]:
    return {item for item in _text(value).split(";") if item}


def _bool(value: Any) -> bool:
    return value is True or _text(value).strip().lower() in {"1", "true", "yes", "on"}


def _standard_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    return (
        _text(row.get("Modification_ID")), _text(row.get("Parent_Fragment_ID")),
        _text(row.get("Candidate_Position_In_Parent")),
        _text(row.get("Candidate_tRNA_Position")),
    )


def _composite_key(row: dict[str, Any]) -> tuple[str, str, str]:
    return (
        _text(row.get("Candidate_ID")), _text(row.get("Complete_Structure_ID")),
        _text(row.get("Composite_Position")),
    )


def _fragment_row(value: Any) -> dict[str, Any]:
    row = _row(value)
    return {
        "Parent_Fragment_ID": row.get("Parent_Fragment_ID", row.get("fragment_id", row.get("Fragment_ID", ""))),
        "Parent_Start": row.get("Parent_Start", row.get("start", row.get("Start_Position", ""))),
        "Parent_End": row.get("Parent_End", row.get("end", row.get("End_Position", ""))),
    }


def _deduplicate(rows: list[Any], key_function) -> list[dict[str, Any]]:
    unique: dict[tuple[str, ...], dict[str, Any]] = {}
    for value in rows or []:
        row = _row(value)
        key = key_function(row)
        unique.setdefault(key, row)
    return [unique[key] for key in sorted(unique)]


def _mass_status(standard: dict[str, Any], composite: dict[str, Any], modification_id: str) -> str:
    standard_mass = _number(standard.get("Mass_Shift"))
    composite_mass = _number(composite.get("Exact_Mass_Delta"))
    if modification_id in _ids(composite.get("Mass_Equivalent_Modification_IDs")):
        return "MATCH"
    if standard_mass is None or composite_mass is None:
        return "NOT_ASSESSED"
    return "MATCH" if abs(standard_mass - composite_mass) <= 1e-4 else "CONFLICT"


def build_rnase_ms2_standard_composite_crosswalk(
    standard_candidate_rows: list[dict[str, Any]],
    composite_position_rows: list[dict[str, Any]],
    composite_bond_rows: list[dict[str, Any]],
    theoretical_fragment_rows: list[Any],
) -> RNaseMS2StandardCompositeCrosswalkResult:
    """Build a non-propagating candidate crosswalk from existing tables only."""
    standards = _deduplicate(list(standard_candidate_rows or []), _standard_key)
    composites = _deduplicate(list(composite_position_rows or []), _composite_key)
    # Bond rows are intentionally accepted but never used as nucleoside components.
    _ = list(composite_bond_rows or [])
    fragments = [_fragment_row(value) for value in theoretical_fragment_rows or []]
    fragment_by_id: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for fragment in fragments:
        fragment_by_id[_text(fragment["Parent_Fragment_ID"])].append(fragment)

    provisional: list[dict[str, Any]] = []
    for standard in standards:
        modification_id, parent_id, parent_position_text, trna_position = _standard_key(standard)
        parent_position = _integer(parent_position_text)
        parent_fragments = fragment_by_id.get(parent_id, [])
        parent_start_values = {_integer(row.get("Parent_Start")) for row in parent_fragments}
        parent_start_values.discard(None)
        absolute = (
            next(iter(parent_start_values)) + parent_position - 1
            if len(parent_start_values) == 1 and parent_position is not None else None
        )
        for composite in composites:
            composite_position = _integer(composite.get("Composite_Position"))
            explicit = _ids(composite.get("Explicit_Legacy_Modification_IDs"))
            mass_only = _ids(composite.get("Mass_Equivalent_Modification_IDs"))
            position_status = (
                "MATCH" if absolute is not None and composite_position == absolute else
                "CONFLICT" if absolute is not None and composite_position is not None else "INSUFFICIENT_PROVENANCE"
            )
            if not parent_fragments:
                parent_status = "INSUFFICIENT_PROVENANCE"
            else:
                contains = [
                    row for row in parent_fragments
                    if _integer(row.get("Parent_Start")) is not None
                    and _integer(row.get("Parent_End")) is not None
                    and composite_position is not None
                    and _integer(row["Parent_Start"]) <= composite_position <= _integer(row["Parent_End"])
                ]
                parent_status = "MATCH" if contains else "CONFLICT"
            identity_status = (
                "EXPLICIT_MATCH" if modification_id and modification_id in explicit else
                "MASS_EQUIVALENT_ONLY" if modification_id and modification_id in mass_only else
                "CONFLICT" if modification_id and explicit else "UNRESOLVED"
            )
            mass_status = _mass_status(standard, composite, modification_id)
            isomer = _bool(composite.get("Is_Isomeric")) or bool(_text(composite.get("Isomer_Group_ID")))
            missing = []
            if not modification_id: missing.append("missing_standard_modification_id")
            if not parent_id: missing.append("missing_parent_fragment_id")
            if parent_position is None: missing.append("missing_candidate_position_in_parent")
            if absolute is None: missing.append("absolute_position_not_reconstructable")
            if composite_position is None: missing.append("missing_composite_position")
            if not _text(composite.get("Canonical_Structure_ID")): missing.append("missing_canonical_structure_id")
            if parent_status == "INSUFFICIENT_PROVENANCE": missing.append("parent_fragment_provenance_unavailable")

            if missing:
                status = "INSUFFICIENT_PROVENANCE"
            elif parent_status == "CONFLICT":
                status = "PARENT_FRAGMENT_CONFLICT"
            elif position_status == "CONFLICT":
                status = "POSITION_CONFLICT"
            elif identity_status == "EXPLICIT_MATCH":
                if mass_status == "CONFLICT":
                    status = "MODIFICATION_IDENTITY_CONFLICT"
                elif isomer:
                    status = "POSITION_MATCH_IDENTITY_UNRESOLVED"
                else:
                    status = "EXACT_MATCH"
            elif identity_status == "MASS_EQUIVALENT_ONLY":
                status = "MASS_EQUIVALENT_ONLY"
            elif position_status == "MATCH" and identity_status == "UNRESOLVED":
                status = "POSITION_MATCH_IDENTITY_UNRESOLVED"
            elif position_status == "MATCH" and identity_status == "CONFLICT":
                status = "MODIFICATION_IDENTITY_CONFLICT"
            else:
                status = "NOT_MAPPABLE"

            basis = [
                f"parent_fragment={parent_status.lower()}",
                f"absolute_position={position_status.lower()}",
                f"identity={identity_status.lower()}",
                f"mass_shift={mass_status.lower()}",
            ]
            if isomer: basis.append("structural_isomer_sharing")
            limiting = list(missing)
            if parent_status == "CONFLICT": limiting.append("parent_fragment_conflict")
            if position_status == "CONFLICT": limiting.append("absolute_position_conflict")
            if identity_status == "CONFLICT": limiting.append("explicit_identity_conflict")
            if identity_status == "MASS_EQUIVALENT_ONLY": limiting.append("mass_equivalence_is_not_identity")
            if mass_status == "CONFLICT": limiting.append("mass_shift_conflict")
            if isomer: limiting.append("structural_isomer_prevents_exact_match")
            provisional.append({
                "Modification_ID": modification_id, "Parent_Fragment_ID": parent_id,
                "Candidate_Position_In_Parent": standard.get("Candidate_Position_In_Parent", ""),
                "Candidate_tRNA_Position": trna_position,
                "Standard_Absolute_Sequence_Position": absolute if absolute is not None else "",
                "Candidate_ID": composite.get("Candidate_ID", ""),
                "Complete_Structure_ID": composite.get("Complete_Structure_ID", ""),
                "Composite_Position": composite.get("Composite_Position", ""),
                "Applied_Transform_IDs": composite.get("Applied_Transform_IDs", ""),
                "Explicit_Legacy_Modification_IDs": composite.get("Explicit_Legacy_Modification_IDs", ""),
                "Mass_Equivalent_Modification_IDs": composite.get("Mass_Equivalent_Modification_IDs", ""),
                "Canonical_Structure_ID": composite.get("Canonical_Structure_ID", ""),
                "Position_Match_Status": position_status,
                "Modification_Identity_Match_Status": identity_status,
                "Mass_Shift_Match_Status": mass_status,
                "Parent_Fragment_Match_Status": parent_status,
                "Structural_Isomer_Sharing": isomer, "Crosswalk_Status": status,
                "Crosswalk_Cardinality": "", "Crosswalk_Basis": ";".join(basis),
                "Limiting_Reasons": ";".join(dict.fromkeys(limiting)), **FORMAL_FALSE,
                "_standard_key": _standard_key(standard), "_composite_key": _composite_key(composite),
            })

    related = [row for row in provisional if row["Crosswalk_Status"] in RELATED_STATUSES]
    composites_by_standard: dict[tuple[str, ...], set[tuple[str, ...]]] = defaultdict(set)
    standards_by_composite: dict[tuple[str, ...], set[tuple[str, ...]]] = defaultdict(set)
    for row in related:
        composites_by_standard[row["_standard_key"]].add(row["_composite_key"])
        standards_by_composite[row["_composite_key"]].add(row["_standard_key"])
    for row in provisional:
        if row["Crosswalk_Status"] not in RELATED_STATUSES:
            continue
        standard_count = len(composites_by_standard[row["_standard_key"]])
        composite_count = len(standards_by_composite[row["_composite_key"]])
        row["Crosswalk_Cardinality"] = (
            "MANY_TO_MANY" if standard_count > 1 and composite_count > 1 else
            "ONE_TO_MANY" if standard_count > 1 else
            "MANY_TO_ONE" if composite_count > 1 else "ONE_TO_ONE"
        )
    output = []
    for row in provisional:
        output.append({column: row.get(column, "") for column in CROSSWALK_COLUMNS})
    output.sort(key=lambda row: (
        row["Modification_ID"], row["Parent_Fragment_ID"],
        _text(row["Candidate_Position_In_Parent"]), row["Candidate_ID"],
        row["Complete_Structure_ID"], _text(row["Composite_Position"]),
    ))
    statuses = Counter(row["Crosswalk_Status"] for row in output)
    cardinalities = Counter(row["Crosswalk_Cardinality"] for row in output)
    summary = [{
        "Standard_Candidate_Count": len(standards),
        "Composite_Position_Component_Count": len(composites),
        "Crosswalk_Row_Count": len(output), "Related_Edge_Count": len(related),
        **{f"{status}_Count": statuses[status] for status in (
            "EXACT_MATCH", "POSITION_MATCH_IDENTITY_UNRESOLVED", "MASS_EQUIVALENT_ONLY",
            "POSITION_CONFLICT", "MODIFICATION_IDENTITY_CONFLICT", "PARENT_FRAGMENT_CONFLICT",
            "NOT_MAPPABLE", "INSUFFICIENT_PROVENANCE",
        )},
        **{f"{value}_Count": cardinalities[value] for value in (
            "ONE_TO_ONE", "ONE_TO_MANY", "MANY_TO_ONE", "MANY_TO_MANY",
        )},
        **FORMAL_FALSE,
    }]
    return RNaseMS2StandardCompositeCrosswalkResult(summary, output)
