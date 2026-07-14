"""Complete structure states and exact fragment composition extraction."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping
from rna_masshunter.backbone_state import BondState, BackboneTransformation, normal_bond
from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.modification_composer import apply_transform_ids
from rna_masshunter.sample_structure_schema import SampleStructureHypothesis

RNA_RESIDUE_COMPOSITIONS = {
    "A": {"C": 10, "H": 12, "N": 5, "O": 6, "P": 1},
    "G": {"C": 10, "H": 12, "N": 5, "O": 7, "P": 1},
    "C": {"C": 9, "H": 12, "N": 3, "O": 7, "P": 1},
    "U": {"C": 9, "H": 11, "N": 2, "O": 8, "P": 1},
}
TERMINAL_ADJUSTMENTS = {
    "default": {"H": 2, "O": 1}, "inherited": {"H": 2, "O": 1},
    "dephosphorylated": {"H": 2, "O": 1},
    "residual_phosphate": {"H": 3, "O": 4, "P": 1},
    "cyclic_phosphate": {"H": 1, "O": 3, "P": 1},
}

@dataclass(frozen=True)
class CompleteStructureState:
    candidate_id: str
    position_states: Mapping[int, Any]
    bond_states: Mapping[str, BondState]
    terminal_five_prime: str
    terminal_three_prime: str

@dataclass(frozen=True)
class StructureFragment:
    candidate_id: str
    complete_structure_id: str
    fragment_id: str
    fragment_type: str
    start: int
    end: int
    sequence: str
    included_positions: tuple[int, ...]
    included_modified_positions: tuple[int, ...]
    included_backbone_bonds: tuple[str, ...]
    included_backbone_modifications: tuple[str, ...]
    terminal_form: str
    elemental_composition: ElementalComposition
    neutral_exact_mass: float
    observable: bool = True
    not_observable_reason: str = ""
    @property
    def elemental_composition_canonical(self) -> str:
        return self.elemental_composition.canonical_string()

def build_complete_structure_state(hypothesis: SampleStructureHypothesis, sequence: str,
    transformations: list[Any], slot_schema_path: str | Path,
    backbone_transform: BackboneTransformation) -> tuple[CompleteStructureState | None, dict[str, Any] | None]:
    states: dict[int, Any] = {}
    for position in hypothesis.positions:
        try:
            state, result, _ = apply_transform_ids(position.parent_base, position.position,
                position.transformation_ids, transformations, slot_schema_path)
        except (KeyError, ValueError) as exc:
            return None, {"Candidate_ID": hypothesis.hypothesis_id, "Valid": False,
                "Invalid_Reason": "state_construction_error", "Invalid_Detail": str(exc),
                "Applied_To_Formal_Result": False, "Formal_Change_Ready": False}
        if not result.valid:
            return None, {"Candidate_ID": hypothesis.hypothesis_id, "Valid": False,
                "Invalid_Reason": result.reason_code, "Invalid_Detail": result.reason,
                "Applied_To_Formal_Result": False, "Formal_Change_Ready": False}
        states[position.position] = state
    bonds: dict[str, BondState] = {}
    for item in hypothesis.backbone:
        left, right = (int(value) for value in item.bond_id.split("_", 1))
        bond = normal_bond(left, right)
        if item.state == "phosphorothioate":
            bond = bond.apply(backbone_transform)
        bonds[item.bond_id] = bond
    return CompleteStructureState(hypothesis.hypothesis_id, states, bonds,
        hypothesis.five_prime, hypothesis.three_prime), None

def complete_structure_id(structure: CompleteStructureState) -> str:
    parts = [structure.position_states[p].canonical_structure_id for p in sorted(structure.position_states)]
    parts += [f"{b}={structure.bond_states[b].state}" for b in sorted(structure.bond_states)]
    return structure.candidate_id + "|" + "|".join(parts)

def unmodified_fragment_composition(sequence: str, terminal_form: str = "default") -> ElementalComposition:
    counts: dict[str, int] = {}
    for base in sequence.upper().replace("T", "U"):
        if base not in RNA_RESIDUE_COMPOSITIONS:
            raise ValueError(f"unsupported_base:{base}")
        for element, count in RNA_RESIDUE_COMPOSITIONS[base].items():
            counts[element] = counts.get(element, 0) + count
    adjustment = TERMINAL_ADJUSTMENTS.get(terminal_form)
    if adjustment is None:
        raise ValueError(f"unsupported_terminal_form:{terminal_form}")
    for element, count in adjustment.items():
        counts[element] = counts.get(element, 0) + count
    return ElementalComposition(counts)

def extract_fragment_from_structure(structure: CompleteStructureState, sequence: str,
    start: int, end: int, *, fragment_id: str = "", fragment_type: str = "internal",
    terminal_form: str = "default") -> StructureFragment:
    if start < 1 or end < start or end > len(sequence):
        raise ValueError("fragment_position_out_of_range")
    subsequence = sequence[start - 1:end]
    composition = unmodified_fragment_composition(subsequence, terminal_form)
    positions = tuple(range(start, end + 1))
    modified_positions = tuple(p for p in positions if p in structure.position_states)
    internal_bonds = tuple(f"{p}_{p + 1}" for p in range(start, end))
    modified_bonds = tuple(b for b in internal_bonds if b in structure.bond_states
                           and structure.bond_states[b].state != "normal_phosphate")
    for position in modified_positions:
        composition = composition + structure.position_states[position].elemental_composition_delta
    for bond_id in modified_bonds:
        composition = composition + structure.bond_states[bond_id].composition_delta
    if any(value < 0 for value in composition.to_dict().values()):
        raise ValueError("negative_final_element_count")
    return StructureFragment(structure.candidate_id, complete_structure_id(structure),
        fragment_id or f"{structure.candidate_id}_{start}_{end}", fragment_type, start, end,
        subsequence, positions, modified_positions, internal_bonds, modified_bonds,
        terminal_form, composition, composition.exact_mass)