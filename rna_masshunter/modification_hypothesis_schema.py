"""Schema and identity/chemical validation for shadow-only position hypotheses."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import hashlib
import yaml
from rna_masshunter.backbone_state import load_backbone_transformations
from rna_masshunter.cleavage_site_discovery import discover_candidate_cleavage_bonds
from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.enzymes import get_enzyme_rule, normalize_enzyme_name
from rna_masshunter.modification_composer import apply_transform_ids
from rna_masshunter.modification_constraints import load_transformations

FALSE_FLAGS={"Applied_To_Formal_Result":False,"Formal_Change_Ready":False,"Formal_Result_Changed":False}
HYPOTHESIS_TYPES={"nucleoside_modification","backbone_modification","terminal_state","composite_structure","cleavage_behavior","absence_hypothesis"}
PRIOR_STATUSES={"known","expected","suspected","possible","exploratory","negative_control"}
PRIOR_STRENGTHS={"strong","moderate","weak","unknown"}
PRIORITIES={"critical","high","normal","low"}
EVIDENCE_CODES={"known_database_annotation","literature_supported","organism_specific_expectation","tRNA_position_conservation","pathway_supported","enzyme_supported","prior_MS1_observation","prior_MS2_observation","cross_run_recurrence","manual_expert_hypothesis","negative_control","unknown"}

@dataclass(frozen=True)
class ModificationPositionHypothesis:
    target_id:str; hypothesis_id:str; hypothesis_type:str; priority:str; prior_status:str; prior_strength:str
    evidence_basis:tuple[str,...]; positions:tuple[int,...]; exact_position:int|None; position_range:tuple[int,int]|None
    bonds:tuple[str,...]; exact_bond:str; parent_base:str; modification_ids:tuple[str,...]; enzyme_context:tuple[str,...]
    terminal_state:tuple[tuple[str,str],...]; required_observation:tuple[str,...]; preferred_observation:tuple[str,...]
    source_notes:tuple[str,...]; chemical_validity:str="valid"; position_compatibility:str="compatible"; pathway_compatibility:str="compatible"
    component_domains:tuple[str,...]=(); canonical_nucleoside_states:tuple[tuple[int,str],...]=(); canonical_backbone_states:tuple[tuple[str,str],...]=()
    elemental_composition_delta:str="0"; exact_mass_delta:float=0.0; allow_additional_backbone_states:bool=False; alias_of:str=""
    modification_family:str=""; oxidation_state:str="not_applicable"; chemical_model_status:str="standard"; chemical_model_note:str=""

@dataclass(frozen=True)
class ModificationHypothesisLoadResult:
    schema_version:int; hypotheses:tuple[ModificationPositionHypothesis,...]; invalid_rows:tuple[dict[str,Any],...]; target_rows:tuple[dict[str,Any],...]

def _invalid(target_id,hid,reason,detail=""):
    return {"Target_ID":target_id,"Hypothesis_ID":hid,"Valid":False,"Invalid_Reason":reason,"Invalid_Detail":str(detail),**FALSE_FLAGS}

def _items(value):
    if value in (None,""):return ()
    return tuple(str(x) for x in (value if isinstance(value,list) else [value]))

def _positions(payload,length):
    exact=payload.get("position");position_range=payload.get("position_range") or {};candidates=payload.get("candidate_positions") or []
    values=[];reasons=[];rng=None
    if exact not in (None,""):
        try:exact=int(exact);values.append(exact)
        except (TypeError,ValueError):reasons.append("invalid_position");exact=None
    else:exact=None
    if position_range:
        try:
            start=int(position_range["start"]);end=int(position_range["end"]);rng=(start,end)
            if start<1 or end<start:reasons.append("invalid_position_range")
            else:values.extend(range(start,end+1))
        except (KeyError,TypeError,ValueError):reasons.append("invalid_position_range")
    for raw in candidates:
        try:values.append(int(raw))
        except (TypeError,ValueError):reasons.append("invalid_candidate_position")
    if any(x<1 or x>length for x in values):reasons.append("position_out_of_range")
    return tuple(sorted(set(values))),exact,rng,reasons

def _bonds(payload,length):
    exact=str(payload.get("bond_id") or "");raw=([exact] if exact else [])+list(payload.get("candidate_bonds") or []);valid=[];reasons=[]
    for value in raw:
        try:left,right=(int(x) for x in str(value).split("_",1))
        except (TypeError,ValueError):reasons.append("invalid_bond");continue
        if left<1 or right!=left+1 or right>length:reasons.append("invalid_bond")
        else:valid.append(f"{left}_{right}")
    return tuple(dict.fromkeys(valid)),exact,reasons

def load_modification_position_hypotheses(path:str|Path,*,project_root:str|Path,sequence:str,sequence_id:str,sequence_name:str,organism:str,rule_set:str):
    source=Path(path)
    if not source.is_file():return ModificationHypothesisLoadResult(1,(),(_invalid("","","hypothesis_file_not_found",source),),())
    raw=yaml.safe_load(source.read_text(encoding="utf-8")) or {};version=int(raw.get("schema_version",0) or 0)
    if version!=1:return ModificationHypothesisLoadResult(version,(),(_invalid("","","unsupported_schema_version",version),),())
    root=Path(project_root);transforms=load_transformations(root/"data/modification_transforms_v2.yaml");transform_map={x.id:x for x in transforms}
    backbone={x.id:x for x in load_backbone_transformations(root/"data/backbone_modifications.yaml")};seq=str(sequence).upper().replace("T","U");sha=hashlib.sha256(seq.encode()).hexdigest()
    valid=[];invalid=[];targets=[];seen=set()
    for target_index,target in enumerate(raw.get("targets") or (),1):
        tid=str(target.get("target_id") or f"TARGET_{target_index}");aliases=set(_items(target.get("sequence_name_aliases")))|set(_items(target.get("sequence_id_aliases")))
        checks=[("target_sequence_id_mismatch",sequence_id in {str(target.get("sequence_id") or ""),*aliases}),
            ("target_sequence_name_mismatch",sequence_name in {str(target.get("sequence_name") or ""),*aliases}),
            ("target_sequence_length_mismatch",int(target.get("sequence_length") or -1)==len(seq)),
            ("target_sequence_sha256_mismatch",str(target.get("sequence_sha256") or "")==sha),
            ("target_organism_mismatch",str(target.get("organism") or "").replace("_"," ").casefold()==str(organism).replace("_"," ").casefold()),
            ("target_rule_set_mismatch",str(target.get("rule_set") or "")==str(rule_set))]
        failures=[name for name,ok in checks if not ok]
        targets.append({"Target_ID":tid,"Target_Valid":not failures,"Invalid_Reason":";".join(failures),**FALSE_FLAGS})
        if failures:
            invalid.extend(_invalid(tid,str((p or {}).get("hypothesis_id") or ""),";".join(failures),"target_identity_mismatch") for p in target.get("hypotheses") or [{}]);continue
        for index,payload in enumerate(target.get("hypotheses") or (),1):
            payload=payload or {};hid=str(payload.get("hypothesis_id") or f"HYP_{index:04d}");reasons=[];details=[]
            if hid in seen:reasons.append("duplicate_hypothesis_id")
            seen.add(hid);kind=str(payload.get("hypothesis_type") or "")
            if kind not in HYPOTHESIS_TYPES:reasons.append("invalid_hypothesis_type")
            status=str(payload.get("prior_status") or "possible");strength=str(payload.get("prior_strength") or "unknown");priority=str(payload.get("priority") or "normal")
            if status not in PRIOR_STATUSES:reasons.append("invalid_prior_status")
            if strength not in PRIOR_STRENGTHS:reasons.append("invalid_prior_strength")
            if priority not in PRIORITIES:reasons.append("invalid_priority")
            basis=_items(payload.get("evidence_basis")) or ("unknown",)
            unknown=sorted(set(basis)-EVIDENCE_CODES)
            if unknown:reasons.append("invalid_evidence_code");details.extend(unknown)
            positions,exact,rng,pos_reasons=_positions(payload,len(seq));reasons.extend(pos_reasons)
            bonds,exact_bond,bond_reasons=_bonds(payload,len(seq));reasons.extend(bond_reasons)
            parent=str(payload.get("parent_base") or "").upper().replace("T","U")
            mods=_items(payload.get("modification_ids") or payload.get("modification_id"))
            nucleoside_ids=tuple(x for x in mods if x in transform_map);backbone_ids=tuple(x for x in mods if x in backbone)
            terminal_payload=payload.get("terminal_state") or {}
            component_domains=tuple(dict.fromkeys((["nucleoside"] if nucleoside_ids else [])+(["backbone"] if backbone_ids else [])+(["terminal"] if terminal_payload else [])+(["cleavage_behavior"] if kind=="cleavage_behavior" else [])))
            if nucleoside_ids and not positions:reasons.append("missing_position")
            if backbone_ids and not bonds:reasons.append("missing_bond")
            if kind in {"nucleoside_modification","absence_hypothesis"} and not positions:reasons.append("missing_position")
            if kind in {"backbone_modification","cleavage_behavior"} and not bonds:reasons.append("missing_bond")
            component_count=len(nucleoside_ids)+(len(backbone_ids)*max(1,len(bonds)))+len(terminal_payload)
            if kind=="composite_structure" and component_count<2:reasons.append("composite_requires_multiple_components")
            if positions and parent and ((rng and not any(seq[x-1]==parent for x in positions)) or (not rng and any(seq[x-1]!=parent for x in positions))):reasons.append("parent_base_mismatch")
            missing=[x for x in mods if x not in transform_map and x not in backbone]
            if missing:reasons.append("transform_not_found");details.extend(missing)
            chemical="valid";position_ok="compatible";pathway="compatible";canonical_states=[];composition=ElementalComposition.delta()
            if positions and nucleoside_ids and not missing:
                for pos in positions:
                    if parent and seq[pos-1]!=parent and rng:continue
                    try:state,result,_=apply_transform_ids(parent or seq[pos-1],pos,nucleoside_ids,transforms,root/"data/nucleoside_slots.yaml",pathway_context=rule_set)
                    except (KeyError,ValueError) as exc:reasons.append("transform_constraint_violation");details.append(str(exc));continue
                    if not result.valid:reasons.append(result.reason_code);details.append(result.reason);chemical="invalid"
                    else:canonical_states.append((pos,state.canonical_structure_id));composition=composition+state.elemental_composition_delta
                    if not result.position_compatible:position_ok="incompatible"
                    if not result.pathway_compatible:pathway="incompatible"
            enzymes=tuple(normalize_enzyme_name(x) for x in _items(payload.get("enzyme_context")))
            for enzyme in enzymes:
                try:get_enzyme_rule(enzyme)
                except ValueError:reasons.append("unknown_enzyme");continue
                if kind=="cleavage_behavior" and bonds and bonds[0] not in {x.bond_id for x in discover_candidate_cleavage_bonds(seq,sequence_id,enzyme)}:reasons.append("bond_not_normal_cleavage_site")
                for mod in mods:
                    if mod in backbone and backbone[mod].rule_for(enzyme)[0]=="unknown":reasons.append("backbone_enzyme_rule_not_found")
            if bonds:
                left,right=(int(x) for x in bonds[0].split("_"));expected_left=str(payload.get("left_base") or "").upper();expected_right=str(payload.get("right_base") or "").upper()
                if expected_left and seq[left-1]!=expected_left:reasons.append("left_base_mismatch")
                if expected_right and seq[right-1]!=expected_right:reasons.append("right_base_mismatch")
            canonical_backbone=[]
            for bond_id in bonds:
                for mod in backbone_ids:
                    canonical_backbone.append((bond_id,backbone[mod].to_state));composition=composition+backbone[mod].composition_delta
            terminal=tuple(sorted((str(k),str(v)) for k,v in terminal_payload.items()))
            if any(k not in {"five_prime","three_prime"} or v not in {"default","inherited","dephosphorylated","residual_phosphate","cyclic_phosphate"} for k,v in terminal):reasons.append("invalid_terminal_state")
            if reasons:
                invalid.append(_invalid(tid,hid,";".join(dict.fromkeys(reasons)),";".join(details)));continue
            selected=[transform_map[x] for x in nucleoside_ids]
            oxidation_transform=next((t for t in selected if "side_chain_carbonyl" in t.target_slot),None)
            family=oxidation_transform.target_slot if oxidation_transform else ""
            state_name=oxidation_transform.to_state if oxidation_transform else ""
            oxidation={"thioamide_sulfur":"unoxidized","oxidized_sulfur_1":"monooxide","oxidized_sulfur_2":"dioxide","carbonyl_oxygen":"precursor"}.get(state_name,"not_applicable")
            shadow_transforms=[t.id for t in selected if any("side_chain_carbonyl" in slot and value=="carbonyl_oxygen" for slot,value in t.sets_states)]
            shadow_model=bool(shadow_transforms);model_status="hypothesis_shadow_model" if shadow_model else "standard"
            model_note=((", ".join(shadow_transforms)+" transform uses a carbonyl-bearing shadow state to permit a side-chain thioamide transformation; this is not by itself a confirmed structural assignment.") if shadow_model else "")
            valid.append(ModificationPositionHypothesis(tid,hid,kind,priority,status,strength,basis,positions,exact,rng,bonds,exact_bond,parent,mods,enzymes,terminal,
                _items(payload.get("required_observation")),_items(payload.get("preferred_observation")),_items(payload.get("source_notes")),chemical,position_ok,pathway,
                component_domains,tuple(canonical_states),tuple(sorted(canonical_backbone)),composition.canonical_string(),composition.exact_mass,
                bool(payload.get("allow_additional_backbone_states",False)),str(payload.get("alias_of") or ""),family,oxidation,model_status,model_note))
    order={"critical":0,"high":1,"normal":2,"low":3}
    valid.sort(key=lambda x:(order[x.priority],x.hypothesis_id))
    return ModificationHypothesisLoadResult(version,tuple(valid),tuple(invalid),tuple(targets))
