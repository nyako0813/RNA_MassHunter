"""Normalized component provenance derived directly from complete structure states."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Any

from rna_masshunter.structure_fragment import complete_structure_id

FORMAL_FALSE = {
    "Applied_To_Formal_Result": False,
    "Formal_Change_Ready": False,
    "Formal_Result_Changed": False,
}

POSITION_MAP_COLUMNS = [
    "Candidate_ID", "Complete_Structure_ID", "Composite_Position", "Parent_Base",
    "Applied_Transform_IDs", "Explicit_Legacy_Modification_IDs",
    "Mass_Equivalent_Modification_IDs", "Canonical_Structure_ID",
    "Elemental_Composition_Delta", "Exact_Mass_Delta", "Chemical_Status",
    "Is_Isomeric", "Isomer_Group_ID", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]

BOND_MAP_COLUMNS = [
    "Candidate_ID", "Complete_Structure_ID", "Bond_ID", "Left_Position",
    "Right_Position", "Backbone_State", "Applied_Backbone_Transform_IDs",
    "Elemental_Composition_Delta", "Exact_Mass_Delta", "Stereochemistry",
    "Evidence_Status", "Applied_To_Formal_Result", "Formal_Change_Ready",
    "Formal_Result_Changed",
]


@dataclass(frozen=True)
class CompositeStructureProvenanceResult:
    position_rows: list[dict[str, Any]]
    bond_rows: list[dict[str, Any]]

    @property
    def sheets(self) -> dict[str, list[dict[str, Any]]]:
        return {
            "Composite_Structure_Position_Map": self.position_rows,
            "Composite_Structure_Bond_Map": self.bond_rows,
        }


def _joined(values: set[str] | tuple[str, ...]) -> str:
    return ";".join(sorted(str(value) for value in values if str(value)))


def _modification_value(item: Any, name: str, default: Any = None) -> Any:
    if isinstance(item, dict):
        return item.get(name, default)
    return getattr(item, name, default)


def _mass_equivalent_ids(state: Any, legacy_modifications: list[Any], explicit: set[str], tolerance: float) -> set[str]:
    output: set[str] = set()
    for modification in legacy_modifications:
        modification_id = str(
            _modification_value(modification, "id", "")
            or _modification_value(modification, "Modification_ID", "")
        )
        targets = (
            _modification_value(modification, "target_bases", None)
            or _modification_value(modification, "Target_Bases", None)
            or [_modification_value(modification, "Target_Base", "")]
        )
        if isinstance(targets, str):
            targets = [value for value in targets.replace(",", ";").split(";") if value]
        raw_mass = _modification_value(modification, "mass_shift_from_unmodified", None)
        if raw_mass is None:
            raw_mass = _modification_value(modification, "Mass_Shift", None)
        try:
            mass = float(raw_mass)
        except (TypeError, ValueError):
            continue
        normalized_targets = {str(value).upper().replace("T", "U") for value in targets or ()}
        if modification_id and modification_id not in explicit and state.parent_base in normalized_targets:
            if abs(mass - float(state.exact_mass_delta)) <= tolerance:
                output.add(modification_id)
    return output


def build_composite_structure_provenance(
    structures: list[Any] | tuple[Any, ...],
    legacy_modifications: list[Any] | tuple[Any, ...] = (),
    *,
    mass_tolerance_da: float = 1e-4,
) -> CompositeStructureProvenanceResult:
    """Build position and bond maps without parsing Complete_Structure_ID strings."""
    structures = sorted(
        list(structures or ()),
        key=lambda item: (str(item.candidate_id), complete_structure_id(item)),
    )
    legacy_modifications = list(legacy_modifications or ())

    state_groups: dict[tuple[int, str, str], set[str]] = defaultdict(set)
    for structure in structures:
        for position, state in structure.position_states.items():
            composition = state.elemental_composition_delta.canonical_string()
            state_groups[(int(position), str(state.parent_base), composition)].add(
                str(state.canonical_structure_id)
            )
    isomer_groups: dict[tuple[int, str, str], str] = {}
    isomer_index = 0
    for key in sorted(state_groups):
        if len(state_groups[key]) > 1:
            isomer_index += 1
            isomer_groups[key] = f"CMP_STRUCT_ISO_{isomer_index:06d}"

    position_rows: list[dict[str, Any]] = []
    bond_rows: list[dict[str, Any]] = []
    for structure in structures:
        structure_id = complete_structure_id(structure)
        for position, state in sorted(structure.position_states.items()):
            explicit: set[str] = set()
            transform_map = state.transform_map or {}
            for transform_id in state.applied_transform_ids:
                transform = transform_map.get(transform_id)
                if transform is not None:
                    explicit.update(str(value) for value in transform.legacy_ids if str(value))
            composition = state.elemental_composition_delta.canonical_string()
            group_key = (int(position), str(state.parent_base), composition)
            position_rows.append({
                "Candidate_ID": structure.candidate_id,
                "Complete_Structure_ID": structure_id,
                "Composite_Position": int(position),
                "Parent_Base": state.parent_base,
                "Applied_Transform_IDs": _joined(state.applied_transform_ids),
                "Explicit_Legacy_Modification_IDs": _joined(explicit),
                "Mass_Equivalent_Modification_IDs": _joined(_mass_equivalent_ids(
                    state, legacy_modifications, explicit, mass_tolerance_da,
                )),
                "Canonical_Structure_ID": state.canonical_structure_id,
                "Elemental_Composition_Delta": composition,
                "Exact_Mass_Delta": state.exact_mass_delta,
                "Chemical_Status": state.chemical_status,
                "Is_Isomeric": group_key in isomer_groups,
                "Isomer_Group_ID": isomer_groups.get(group_key, ""),
                **FORMAL_FALSE,
            })
        for bond_id, state in sorted(structure.bond_states.items()):
            bond_rows.append({
                "Candidate_ID": structure.candidate_id,
                "Complete_Structure_ID": structure_id,
                "Bond_ID": bond_id,
                "Left_Position": state.left_position,
                "Right_Position": state.right_position,
                "Backbone_State": state.state,
                "Applied_Backbone_Transform_IDs": _joined(state.applied_transform_ids),
                "Elemental_Composition_Delta": state.composition_delta.canonical_string(),
                "Exact_Mass_Delta": state.exact_mass_delta,
                "Stereochemistry": state.stereochemistry,
                "Evidence_Status": state.evidence_status,
                **FORMAL_FALSE,
            })
    return CompositeStructureProvenanceResult(position_rows, bond_rows)
