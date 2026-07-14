"""Non-propagating enzyme cleavage evaluation against hypothetical bond states."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Mapping
from rna_masshunter.backbone_state import BackboneTransformation, BondState
from rna_masshunter.enzymes import find_cleavage_sites, normalize_enzyme_name
from rna_masshunter.masses import calculate_unmodified_rna_mass

@dataclass(frozen=True)
class CleavageSiteEvaluation:
    enzyme: str
    position: int
    bond_id: str
    sequence_derived: bool
    status: str
    blocking_rule: str
    blocking_reason: str
    evidence_status: str
    formal_cleavage_allowed: bool
    shadow_cleavage_allowed: bool

@dataclass(frozen=True)
class ShadowFragment:
    fragment_id: str
    start: int
    end: int
    sequence: str
    terminal_form: str
    missed_cleavage_count: int
    stochastic_missed_cleavage_count: int
    blocked_cleavage_count: int
    blocked_cleavage_positions: tuple[int,...]
    blocked_cleavage_bond_ids: tuple[str,...]
    blocked_cleavage_reasons: tuple[str,...]
    contains_phosphorothioate: bool
    backbone_modification_count: int
    backbone_mass_delta: float
    mass: float
    cleavage_origin: str

@dataclass(frozen=True)
class CleavageShadowResult:
    candidate_sites: tuple[int,...]
    allowed_sites: tuple[int,...]
    blocked_sites: tuple[int,...]
    evaluations: tuple[CleavageSiteEvaluation,...]
    fragments: tuple[ShadowFragment,...]


def evaluate_cleavage(sequence:str,enzyme:str,bonds:Mapping[str,BondState],transform:BackboneTransformation,base_masses:dict|None=None,*,stochastic_missed_sites=(),terminal_form="default") -> CleavageShadowResult:
    sequence=sequence.upper().replace("T","U"); enzyme=normalize_enzyme_name(enzyme)
    candidates=tuple(sorted(x for x in set(find_cleavage_sites(sequence,enzyme)) if x<len(sequence)))
    candidate_set=set(candidates); allowed=set(candidates); blocked=set(); evaluations=[]
    for pos in range(1,len(sequence)):
        bond_id=f"{pos}_{pos+1}"; bond=bonds.get(bond_id)
        derived=pos in candidate_set
        if not bond or bond.state!="phosphorothioate":
            status="allowed" if derived else "not_applicable"; rule="normal_phosphate"
            reason="sequence-derived site is unaffected" if derived else "not a sequence-derived cleavage site"; evidence="curated"
        else:
            configured_status,evidence=transform.rule_for(enzyme); rule=f"phosphorothioate:{configured_status}"
            status=configured_status if derived else "not_applicable"
            reason=(f"phosphorothioate bond {bond_id} is {configured_status} for {enzyme} in the Phase-1 rule schema" if derived else "bond is not a sequence-derived cleavage site")
            if derived and configured_status=="blocked": blocked.add(pos); allowed.discard(pos)
        evaluations.append(CleavageSiteEvaluation(enzyme,pos,bond_id,derived,status,rule,reason,evidence,derived,derived and pos in allowed))
    stochastic=set(int(x) for x in stochastic_missed_sites)&allowed
    cut_sites=sorted(allowed-stochastic)
    boundaries=[0]+cut_sites+[len(sequence)]; fragments=[]
    for i,(left,right) in enumerate(zip(boundaries,boundaries[1:]),1):
        internal_blocked=tuple(sorted(x for x in blocked if left<x<right))
        internal_stochastic=tuple(sorted(x for x in stochastic if left<x<right))
        fragment_bonds=[b for b in bonds.values() if left < b.left_position and b.right_position <= right and b.state=="phosphorothioate"]
        delta=sum(b.exact_mass_delta for b in fragment_bonds)
        if internal_blocked and internal_stochastic: origin="mixed_stochastic_and_blocked"
        elif internal_blocked: origin="phosphorothioate_blocked"
        elif internal_stochastic: origin="stochastic_missed"
        else: origin="normal"
        seq=sequence[left:right]
        base_mass=calculate_unmodified_rna_mass(seq,base_masses,terminal_form=terminal_form) if base_masses else 0.0
        reasons=tuple(f"phosphorothioate blocked {enzyme} cleavage after {x}" for x in internal_blocked)
        fragments.append(ShadowFragment(f"SH-{enzyme}-{left+1}-{right}-{i}",left+1,right,seq,terminal_form,
            len(internal_blocked)+len(internal_stochastic),len(internal_stochastic),len(internal_blocked),internal_blocked,
            tuple(f"{x}_{x+1}" for x in internal_blocked),reasons,bool(fragment_bonds),len(fragment_bonds),delta,float(base_mass or 0)+delta,origin))
    return CleavageShadowResult(candidates,tuple(sorted(allowed)),tuple(sorted(blocked)),tuple(evaluations),tuple(fragments))
