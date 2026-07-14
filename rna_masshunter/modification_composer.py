"""Deterministic bounded search over schema-defined nucleoside transformations."""
from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations, permutations
from pathlib import Path
from typing import Any
from rna_masshunter.modification_constraints import ConstraintResult, Transformation
from rna_masshunter.nucleoside_state import NucleosideState, initial_nucleoside_state

@dataclass(frozen=True)
class CompositeCandidate:
    candidate_id: str
    state: NucleosideState
    transforms: tuple[Transformation,...]
    pathway_compatible: bool
    position_compatible: bool
    isomer_group_id: str=""
    is_isomeric: bool=False
    @property
    def component_count(self): return len(self.transforms)

@dataclass(frozen=True)
class InvalidAttempt:
    attempt_id: str
    position: int
    parent_base: str
    transform_ids: tuple[str,...]
    result: ConstraintResult
    attempted_slots: tuple[str,...]

@dataclass(frozen=True)
class CompositionResult:
    valid_candidates: tuple[CompositeCandidate,...]
    invalid_attempts: tuple[InvalidAttempt,...]
    isomer_groups: tuple[tuple[str,tuple[str,...]],...]


def _try_combo(initial: NucleosideState, combo: tuple[Transformation,...], pathway_context: str|None, organism_context: str|None=None):
    best=None
    for order in permutations(combo):
        state=initial; position_ok=True; pathway_ok=True; failure=None
        for transform in order:
            state2,result=state.apply(transform,pathway_context=pathway_context,organism_context=organism_context)
            position_ok &= result.position_compatible; pathway_ok &= result.pathway_compatible
            if not result.valid: failure=result; break
            state=state2
        if failure is None: return state,ConstraintResult(True,pathway_compatible=pathway_ok,position_compatible=position_ok),order
        if best is None: best=(failure,order)
    return initial,best[0],best[1]


def compose_modifications(parent_base: str,position: int,transformations: list[Transformation],schema_path: str|Path,*,max_components: int=3,pathway_context: str|None=None,organism_context: str|None=None) -> CompositionResult:
    if max_components<1: return CompositionResult((),(),())
    eligible=sorted((t for t in transformations if t.target_scope=="nucleoside" and parent_base.upper().replace("T","U") in t.parent_bases),key=lambda t:t.id)
    initial=initial_nucleoside_state(parent_base,position,schema_path,transformations)
    valid=[]; invalid=[]; seen={}
    attempt_index=0
    for size in range(1,min(max_components,len(eligible))+1):
        for combo in combinations(eligible,size):
            attempt_index+=1
            state,result,order=_try_combo(initial,combo,pathway_context,organism_context)
            ids=tuple(t.id for t in combo); slots=tuple(sorted({s for t in combo for s in t.occupies_slots}))
            if not result.valid:
                invalid.append(InvalidAttempt(f"INV-{parent_base}-{position}-{attempt_index:05d}",position,parent_base,ids,result,slots)); continue
            canonical=(state.slot_states,state.elemental_composition_delta)
            if canonical in seen: continue
            ordered=tuple(sorted(combo,key=lambda t:t.id))
            candidate=CompositeCandidate(f"CMP-{parent_base}-{position}-{len(valid)+1:05d}",state,ordered,result.pathway_compatible,result.position_compatible)
            seen[canonical]=candidate; valid.append(candidate)
    mass_groups={}
    for c in valid: mass_groups.setdefault(c.state.elemental_composition_delta.canonical_string(),[]).append(c)
    isomer_groups=[]; replacements={}
    group_i=0
    for composition,members in sorted(mass_groups.items()):
        structures={m.state.canonical_structure_id for m in members}
        if len(structures)<2: continue
        group_i+=1; gid=f"ISO-{parent_base}-{position}-{group_i:04d}"
        isomer_groups.append((gid,tuple(m.candidate_id for m in members)))
        for m in members: replacements[m.candidate_id]=(gid,True)
    valid=[CompositeCandidate(c.candidate_id,c.state,c.transforms,c.pathway_compatible,c.position_compatible,*replacements.get(c.candidate_id,("",False))) for c in valid]
    return CompositionResult(tuple(valid),tuple(invalid),tuple(isomer_groups))

def apply_transform_ids(parent_base: str,position: int,transform_ids:list[str]|tuple[str,...],transformations:list[Transformation],schema_path:str|Path,*,pathway_context:str|None=None):
    """Apply an explicitly ordered hypothesis fixture without candidate deduplication."""
    mapping={t.id:t for t in transformations}; initial=initial_nucleoside_state(parent_base,position,schema_path,transformations)
    selected=tuple(mapping[x] for x in transform_ids)
    state,result,order=_try_combo(initial,selected,pathway_context)
    return state,result,tuple(t.id for t in order)
