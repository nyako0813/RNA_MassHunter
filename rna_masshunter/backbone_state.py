"""Immutable inter-nucleotide bond states for shadow-only modification audits."""
from __future__ import annotations
from dataclasses import dataclass, replace
from pathlib import Path
import yaml
from rna_masshunter.elemental_composition import ElementalComposition

@dataclass(frozen=True)
class BackboneTransformation:
    id: str
    name: str
    from_state: str
    to_state: str
    composition_delta: ElementalComposition
    stereochemistry: str
    evidence_status: str
    blocking_rules: tuple[tuple[str,str,str],...]
    notes: str=""
    @property
    def exact_mass_delta(self): return self.composition_delta.exact_mass
    def rule_for(self,enzyme):
        for name,status,evidence in self.blocking_rules:
            if name==enzyme: return status,evidence
        return "unknown","unknown"

@dataclass(frozen=True)
class BondState:
    left_position: int
    right_position: int
    state: str="normal_phosphate"
    applied_transform_ids: tuple[str,...]=()
    composition_delta: ElementalComposition=ElementalComposition.delta()
    stereochemistry: str="not_applicable"
    confidence: str="curated"
    evidence_status: str="curated"
    @property
    def bond_id(self): return f"{self.left_position}_{self.right_position}"
    @property
    def exact_mass_delta(self): return self.composition_delta.exact_mass
    def apply(self,transform: BackboneTransformation):
        if transform.id in self.applied_transform_ids: raise ValueError("duplicate_component")
        if self.state!=transform.from_state: raise ValueError("from_state_mismatch")
        return replace(self,state=transform.to_state,applied_transform_ids=self.applied_transform_ids+(transform.id,),
            composition_delta=self.composition_delta+transform.composition_delta,stereochemistry=transform.stereochemistry,
            confidence="hypothesis",evidence_status=transform.evidence_status)


def normal_bond(left_position:int,right_position:int):
    if right_position!=left_position+1: raise ValueError("Bond positions must be adjacent")
    return BondState(left_position,right_position)


def load_backbone_transformations(path: str|Path):
    with Path(path).open(encoding="utf-8") as handle: raw=yaml.safe_load(handle) or {}
    out=[]
    for row in raw.get("transforms",[]):
        rules=tuple((name,str(rule.get("status","unknown")),str(rule.get("evidence_status","unknown"))) for name,rule in sorted((row.get("blocking_rules") or {}).items()))
        out.append(BackboneTransformation(str(row["id"]),str(row["name"]),str(row["from_state"]),str(row["to_state"]),
            ElementalComposition.delta(row.get("composition_delta")),str(row.get("stereochemistry","unknown")),str(row.get("evidence_status","unknown")),rules,str(row.get("notes",""))))
    return out
