"""Schema-driven validation for non-propagating composite modification states."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping
import yaml
from rna_masshunter.elemental_composition import ElementalComposition

@dataclass(frozen=True)
class Transformation:
    id: str
    name: str
    parent_bases: tuple[str, ...]
    target_scope: str
    target_slot: str
    from_state: str
    to_state: str
    composition_delta: ElementalComposition
    occupies_slots: tuple[str, ...]
    sets_states: tuple[tuple[str, str], ...] = ()
    requires: tuple[tuple[str, tuple[str, ...]], ...] = ()
    forbids: tuple[tuple[str, tuple[str, ...]], ...] = ()
    supersedes: tuple[str, ...] = ()
    included_components: tuple[str, ...] = ()
    pathway_tags: tuple[str, ...] = ()
    allowed_positions: tuple[int, ...] = ()
    allowed_organisms: tuple[str, ...] = ()
    allowed_enzymes: tuple[str, ...] = ()
    evidence_status: str = "unknown"
    legacy_ids: tuple[str, ...] = ()
    notes: str = ""

    @property
    def exact_mass_delta(self) -> float: return self.composition_delta.exact_mass
    @property
    def state_updates(self) -> dict[str, str]:
        updates=dict(self.sets_states); updates[self.target_slot]=self.to_state; return updates

@dataclass(frozen=True)
class ConstraintResult:
    valid: bool
    reason_code: str = "valid"
    reason: str = ""
    conflicting_slots: tuple[str, ...] = ()
    missing_requirements: tuple[str, ...] = ()
    forbidden_states: tuple[str, ...] = ()
    duplicate_components: tuple[str, ...] = ()
    superseded_components: tuple[str, ...] = ()
    invalid_composition: bool = False
    chemically_valid: bool = True
    pathway_compatible: bool = True
    position_compatible: bool = True
    pathway_warning: str = ""
    position_warning: str = ""


def _pairs(value: Mapping[str, Any] | None):
    return tuple((str(k),tuple(str(x) for x in (v if isinstance(v,list) else [v]))) for k,v in (value or {}).items())


def load_transformations(path: str | Path) -> list[Transformation]:
    with Path(path).open(encoding="utf-8") as handle: raw=yaml.safe_load(handle) or {}
    out=[]
    for row in raw.get("transforms",[]):
        required=("id","name","parent_bases","target_scope","target_slot","from_state","to_state","composition_delta")
        missing=[key for key in required if key not in row]
        if missing: raise ValueError(f"unsupported_schema:{row.get('id','?')}:{','.join(missing)}")
        delta=ElementalComposition.delta(row.get("composition_delta"))
        supplied=row.get("mass_delta")
        if supplied is not None and abs(float(supplied)-delta.exact_mass)>1e-6:
            raise ValueError(f"mass_delta mismatch for {row['id']}")
        out.append(Transformation(
            id=str(row["id"]),name=str(row["name"]),parent_bases=tuple(row["parent_bases"]),
            target_scope=str(row["target_scope"]),target_slot=str(row["target_slot"]),
            from_state=str(row["from_state"]),to_state=str(row["to_state"]),composition_delta=delta,
            occupies_slots=tuple(row.get("occupies_slots") or [row["target_slot"]]),
            sets_states=tuple((str(k),str(v)) for k,v in (row.get("sets_states") or {}).items()),
            requires=_pairs(row.get("requires")),forbids=_pairs(row.get("forbids")),
            supersedes=tuple(row.get("supersedes") or ()),included_components=tuple(row.get("included_components") or ()),
            pathway_tags=tuple(row.get("pathway_tags") or ()),allowed_positions=tuple(int(x) for x in (row.get("allowed_positions") or ())),
            allowed_organisms=tuple(row.get("allowed_organisms") or ()),allowed_enzymes=tuple(row.get("allowed_enzymes") or ()),
            evidence_status=str(row.get("evidence_status","unknown")),legacy_ids=tuple(row.get("legacy_ids") or ()),notes=str(row.get("notes", "")),
        ))
    ids=[x.id for x in out]
    if len(ids)!=len(set(ids)): raise ValueError("unsupported_schema:duplicate transformation id")
    return out


def validate_transformation(state: Any, transform: Transformation, *, position: int | None = None, pathway_context: str | None = None, organism_context: str | None = None) -> ConstraintResult:
    applied=set(state.applied_transform_ids)
    if transform.id in applied:
        return ConstraintResult(False,"duplicate_component",f"{transform.id} already applied",duplicate_components=(transform.id,),chemically_valid=False)
    if state.parent_base not in transform.parent_bases:
        return ConstraintResult(False,"wrong_parent_base",f"{transform.id} does not apply to {state.parent_base}",chemically_valid=False)
    overlap=sorted(set(transform.included_components)&applied)
    if overlap:
        return ConstraintResult(False,"parent_child_double_count","completed transformation includes an applied parent",duplicate_components=tuple(overlap),chemically_valid=False)
    superseded=sorted((set(transform.supersedes)&applied) | {x for x in applied if transform.id in state.transform_map[x].supersedes}) if getattr(state,"transform_map",None) else sorted(set(transform.supersedes)&applied)
    if superseded:
        return ConstraintResult(False,"superseded_component","parent/derived transformations cannot be added together",superseded_components=tuple(superseded),chemically_valid=False)
    slots=state.slot_state_dict
    if transform.target_slot not in slots:
        return ConstraintResult(False,"unsupported_schema",f"slot {transform.target_slot} is unavailable",chemically_valid=False)
    if slots.get(transform.target_slot)!=transform.from_state:
        code="impossible_oxidation_state" if "oxo" in transform.id else "from_state_mismatch"
        return ConstraintResult(False,code,f"{transform.target_slot} is {slots.get(transform.target_slot)}, expected {transform.from_state}",missing_requirements=(f"{transform.target_slot}={transform.from_state}",),chemically_valid=False)
    conflicts=sorted((set(transform.occupies_slots)-{transform.target_slot})&set(state.occupied_slots))
    if conflicts:
        return ConstraintResult(False,"slot_conflict","occupied slot conflict",conflicting_slots=tuple(conflicts),chemically_valid=False)
    missing=[]
    for slot,allowed in transform.requires:
        if slots.get(slot) not in allowed: missing.append(f"{slot} in {list(allowed)}")
    if missing:
        return ConstraintResult(False,"missing_requirement","required state is absent",missing_requirements=tuple(missing),chemically_valid=False)
    forbidden=[]
    for slot,values in transform.forbids:
        if slots.get(slot) in values: forbidden.append(f"{slot}={slots.get(slot)}")
    if forbidden:
        return ConstraintResult(False,"forbidden_state","forbidden state present",forbidden_states=tuple(forbidden),chemically_valid=False)
    position_ok=not transform.allowed_positions or (position if position is not None else state.position) in transform.allowed_positions
    pathway_ok=(not transform.pathway_tags or pathway_context is None or pathway_context in transform.pathway_tags) and (not transform.allowed_organisms or organism_context is None or organism_context in transform.allowed_organisms)
    compatibility_code="position_disallowed" if not position_ok else ("pathway_disallowed" if not pathway_ok else "valid")
    return ConstraintResult(True,compatibility_code,"chemically valid",pathway_compatible=pathway_ok,position_compatible=position_ok,
        pathway_warning="" if pathway_ok else "pathway_disallowed",position_warning="" if position_ok else "position_disallowed")
