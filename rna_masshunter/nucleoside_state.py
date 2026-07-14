"""Immutable slot/state representation of an RNA nucleoside."""
from __future__ import annotations
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Mapping
import yaml
from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.modification_constraints import ConstraintResult, Transformation, validate_transformation

@dataclass(frozen=True)
class NucleosideState:
    position: int
    parent_base: str
    slot_states: tuple[tuple[str,str], ...]
    applied_transform_ids: tuple[str,...]=()
    occupied_slots: tuple[str,...]=()
    elemental_composition_delta: ElementalComposition=ElementalComposition.delta()
    chemical_status: str="curated"
    pathway_tags: tuple[str,...]=()
    warnings: tuple[str,...]=()
    transform_map: Mapping[str,Transformation] | None=None
    @property
    def current_base_identity(self): return self.slot_state_dict["base_identity"]
    @property
    def slot_state_dict(self): return dict(self.slot_states)
    @property
    def exact_mass_delta(self): return self.elemental_composition_delta.exact_mass
    @property
    def canonical_structure_id(self):
        changed=";".join(f"{k}={v}" for k,v in self.slot_states)
        return f"{self.parent_base}@{self.position}|{changed}"
    def apply(self, transform: Transformation, *, pathway_context: str|None=None, organism_context: str|None=None):
        result=validate_transformation(self,transform,position=self.position,pathway_context=pathway_context,organism_context=organism_context)
        if not result.valid: return self,result
        slots=self.slot_state_dict; slots.update(transform.state_updates)
        status="hypothetical" if "hypothetical" in {self.chemical_status,transform.evidence_status} else ("inferred_valid" if transform.evidence_status=="inferred_valid" else self.chemical_status)
        state=replace(self,slot_states=tuple(sorted(slots.items())),applied_transform_ids=self.applied_transform_ids+(transform.id,),
            occupied_slots=tuple(sorted(set(self.occupied_slots)|set(transform.occupies_slots))),
            elemental_composition_delta=self.elemental_composition_delta+transform.composition_delta,chemical_status=status)
        return state,result


def load_initial_states(path: str|Path) -> dict[str,dict[str,str]]:
    with Path(path).open(encoding="utf-8") as handle: raw=yaml.safe_load(handle) or {}
    states=raw.get("initial_states") or {}
    for base in "AGCU":
        if base not in states or states[base].get("base_identity")!=base: raise ValueError(f"unsupported_schema:missing {base} initial state")
    return states


def initial_nucleoside_state(parent_base: str,position: int,schema_path: str|Path,transforms: list[Transformation]|None=None):
    base=parent_base.upper().replace("T","U"); schemas=load_initial_states(schema_path)
    if base not in schemas: raise ValueError(f"Unsupported parent base: {parent_base}")
    mapping={x.id:x for x in (transforms or [])}
    return NucleosideState(position,base,tuple(sorted((str(k),str(v)) for k,v in schemas[base].items())),transform_map=mapping)
