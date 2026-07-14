"""Phase-1 composite/backbone/cleavage shadow audit; never mutates formal results."""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import json
import time
import tracemalloc
from rna_masshunter.backbone_state import load_backbone_transformations, normal_bond
from rna_masshunter.cleavage_constraints import evaluate_cleavage
from rna_masshunter.modification_composer import compose_modifications
from rna_masshunter.modification_constraints import load_transformations

COMPOSITE_CANDIDATE_COLUMNS=["Candidate_ID","Position","Parent_Base","Complete_Structure_ID","Chemical_Status","Applied_Transform_IDs","Applied_Transform_Names","Component_Count","Slot_State_Summary","Occupied_Slots","Elemental_Composition_Delta","Exact_Mass_Delta","Legacy_Equivalent_IDs","Included_Component_IDs","Is_Composite","Is_Isomeric","Isomer_Group_ID","Chemically_Valid","Pathway_Compatible","Position_Compatible","Evidence_Status","Notes","Applied_To_Formal_Result"]
COMPOSITE_INVALID_COLUMNS=["Attempt_ID","Position","Parent_Base","Attempted_Transform_IDs","Attempted_Slots","Valid","Reason_Code","Invalid_Reason","Conflicting_Slots","Missing_Requirements","Forbidden_States","Duplicate_Components","Superseded_Components","Invalid_Composition","Pathway_Compatible","Position_Compatible","Applied_To_Formal_Result"]
COMPOSITE_SUMMARY_COLUMNS=["Total_Transform_Definitions","Legacy_Mapped_Transform_Count","Legacy_Only_Modification_Count","Valid_Single_Candidate_Count","Valid_Composite_Candidate_Count","Invalid_Combination_Count","Slot_Conflict_Count","Missing_Requirement_Count","Parent_Child_Double_Count_Count","Duplicate_Component_Count","Impossible_Oxidation_Count","Isomer_Group_Count","A_Candidate_Count","G_Candidate_Count","C_Candidate_Count","U_Candidate_Count","Phosphorothioate_Candidate_Count","Formal_Result_Changed","Audit_Mode","Applied_To_Formal_Result","Formal_Change_Ready","Remaining_Risk"]
BACKBONE_COLUMNS=["Backbone_Candidate_ID","Bond_ID","Left_Position","Right_Position","Backbone_State","Applied_Transform_IDs","Elemental_Composition_Delta","Exact_Mass_Delta","Stereochemistry","Evidence_Status","Blocks_RNase_T1","Blocks_RNase_A","Blocks_Nuclease_P1","Blocks_RNase_T2","Chemical_Status","Applied_To_Formal_Result"]
CLEAVAGE_COLUMNS=["Enzyme","Bond_ID","Left_Position","Right_Position","Sequence_Derived_Cleavage_Site","Backbone_State","Cleavage_Status","Blocking_Rule","Blocking_Reason","Evidence_Status","Formal_Cleavage_Allowed","Shadow_Cleavage_Allowed","Resulting_Fragment_ID","Fragment_Start","Fragment_End","Fragment_Sequence","Cleavage_Origin","Stochastic_Missed_Cleavage_Count","Blocked_Cleavage_Count","Blocked_Cleavage_Bond_IDs","Contains_Phosphorothioate","Backbone_Mass_Delta","Shadow_Fragment_Mass","Applied_To_Formal_Result"]
DIAGNOSTIC_COLUMNS=["Composite_Mod_Audit_Available","Composite_Mod_Valid_Count","Composite_Mod_Invalid_Count","Composite_Mod_Isomer_Group_Count","Composite_Mod_Legacy_Overlap_Count","Backbone_Mod_Audit_Available","Phosphorothioate_Candidate_Count","Cleavage_Block_Audit_Available","RNase_T1_Blocked_Cleavage_Count","Nuclease_P1_Blocked_Cleavage_Count","Composite_Mod_Formal_Change_Ready","Composite_Mod_Applied_To_Formal_Result","Backbone_Mod_Applied_To_Formal_Result"]

@dataclass(frozen=True)
class CompositeAuditResult:
    sheets: dict
    diagnostics: dict
    runtime_seconds: float
    peak_memory_mb: float


def _joined(values): return ";".join(str(x) for x in values)
def _json(value): return json.dumps(value,sort_keys=True,separators=(",",":"))


def build_composite_modification_audit(project_root: str|Path,sequence: str,legacy_modifications:list,base_masses:dict,audit_mode="full",max_components=3) -> CompositeAuditResult:
    root=Path(project_root); sequence=sequence.upper().replace("T","U")
    tracemalloc.start(); started=time.perf_counter()
    transforms=load_transformations(root/"data/modification_transforms_v2.yaml")
    legacy_ids={str(getattr(x,"id","") or getattr(x,"symbol","")) for x in legacy_modifications}
    valid_rows=[]; invalid_rows=[]; isomer_ids=set(); reason_counts=Counter(); base_counts=Counter(); overlap_ids=set()
    for position,base in enumerate(sequence,1):
        result=compose_modifications(base,position,transforms,root/"data/nucleoside_slots.yaml",max_components=max_components)
        isomer_ids.update(gid for gid,_ in result.isomer_groups)
        for candidate in result.valid_candidates:
            t_ids=[t.id for t in candidate.transforms]
            directly_mapped={x for t in candidate.transforms for x in t.legacy_ids if x in legacy_ids}
            mass_equivalent={str(getattr(item,"id","") or getattr(item,"symbol","")) for item in legacy_modifications if base in (getattr(item,"target_bases",[]) or []) and getattr(item,"mass_shift_from_unmodified",float("nan"))==getattr(item,"mass_shift_from_unmodified",float("nan")) and abs(float(item.mass_shift_from_unmodified)-candidate.state.exact_mass_delta)<=1e-4}
            legacy=sorted(directly_mapped|mass_equivalent); overlap_ids.update(legacy); base_counts[base]+=1
            included=sorted({x for t in candidate.transforms for x in t.included_components})
            evidence="hypothetical" if any(t.evidence_status=="hypothetical" for t in candidate.transforms) else ("inferred_valid" if any(t.evidence_status=="inferred_valid" for t in candidate.transforms) else "curated")
            valid_rows.append({"Candidate_ID":candidate.candidate_id,"Position":position,"Parent_Base":base,"Complete_Structure_ID":candidate.state.canonical_structure_id,"Chemical_Status":candidate.state.chemical_status,
                "Applied_Transform_IDs":_joined(t_ids),"Applied_Transform_Names":_joined(t.name for t in candidate.transforms),"Component_Count":candidate.component_count,
                "Slot_State_Summary":_joined(f"{k}={v}" for k,v in candidate.state.slot_states),"Occupied_Slots":_joined(candidate.state.occupied_slots),
                "Elemental_Composition_Delta":candidate.state.elemental_composition_delta.canonical_string(),"Exact_Mass_Delta":candidate.state.exact_mass_delta,
                "Legacy_Equivalent_IDs":_joined(legacy),"Included_Component_IDs":_joined(included),"Is_Composite":candidate.component_count>1,"Is_Isomeric":candidate.is_isomeric,
                "Isomer_Group_ID":candidate.isomer_group_id,"Chemically_Valid":True,"Pathway_Compatible":candidate.pathway_compatible,"Position_Compatible":candidate.position_compatible,
                "Evidence_Status":evidence,"Notes":_joined(t.notes for t in candidate.transforms if t.notes),"Applied_To_Formal_Result":False})
        for invalid in result.invalid_attempts:
            r=invalid.result; reason_counts[r.reason_code]+=1
            invalid_rows.append({"Attempt_ID":invalid.attempt_id,"Position":position,"Parent_Base":base,"Attempted_Transform_IDs":_joined(invalid.transform_ids),"Attempted_Slots":_joined(invalid.attempted_slots),"Valid":False,
                "Reason_Code":r.reason_code,"Invalid_Reason":r.reason,"Conflicting_Slots":_joined(r.conflicting_slots),"Missing_Requirements":_joined(r.missing_requirements),"Forbidden_States":_joined(r.forbidden_states),
                "Duplicate_Components":_joined(r.duplicate_components),"Superseded_Components":_joined(r.superseded_components),"Invalid_Composition":r.invalid_composition,
                "Pathway_Compatible":r.pathway_compatible,"Position_Compatible":r.position_compatible,"Applied_To_Formal_Result":False})
    backbone_transform=load_backbone_transformations(root/"data/backbone_modifications.yaml")[0]
    backbone_rows=[]; cleavage_rows=[]; blocked_counts=Counter()
    for left in range(1,len(sequence)):
        bond=normal_bond(left,left+1).apply(backbone_transform); status={enzyme:backbone_transform.rule_for(enzyme)[0] for enzyme in ("RNase_T1","RNase_A","Nuclease_P1","RNase_T2")}
        backbone_rows.append({"Backbone_Candidate_ID":f"BB-PS-{bond.bond_id}","Bond_ID":bond.bond_id,"Left_Position":left,"Right_Position":left+1,"Backbone_State":bond.state,"Applied_Transform_IDs":_joined(bond.applied_transform_ids),
            "Elemental_Composition_Delta":bond.composition_delta.canonical_string(),"Exact_Mass_Delta":bond.exact_mass_delta,"Stereochemistry":bond.stereochemistry,"Evidence_Status":bond.evidence_status,
            "Blocks_RNase_T1":status["RNase_T1"],"Blocks_RNase_A":status["RNase_A"],"Blocks_Nuclease_P1":status["Nuclease_P1"],"Blocks_RNase_T2":status["RNase_T2"],"Chemical_Status":"hypothetical","Applied_To_Formal_Result":False})
        for enzyme in ("RNase_T1","RNase_A","Nuclease_P1","RNase_T2"):
            shadow=evaluate_cleavage(sequence,enzyme,{bond.bond_id:bond},backbone_transform,base_masses)
            ev=next(x for x in shadow.evaluations if x.bond_id==bond.bond_id)
            fragment=next((x for x in shadow.fragments if x.start<=left and x.end>=left+1),None)
            if ev.status=="blocked" and ev.sequence_derived: blocked_counts[enzyme]+=1
            cleavage_rows.append({"Enzyme":enzyme,"Bond_ID":bond.bond_id,"Left_Position":left,"Right_Position":left+1,"Sequence_Derived_Cleavage_Site":ev.sequence_derived,
                "Backbone_State":bond.state,"Cleavage_Status":ev.status,"Blocking_Rule":ev.blocking_rule,"Blocking_Reason":ev.blocking_reason,"Evidence_Status":ev.evidence_status,
                "Formal_Cleavage_Allowed":ev.formal_cleavage_allowed,"Shadow_Cleavage_Allowed":ev.shadow_cleavage_allowed,"Resulting_Fragment_ID":fragment.fragment_id if fragment else "",
                "Fragment_Start":fragment.start if fragment else "","Fragment_End":fragment.end if fragment else "","Fragment_Sequence":fragment.sequence if fragment else "","Cleavage_Origin":fragment.cleavage_origin if fragment else "",
                "Stochastic_Missed_Cleavage_Count":fragment.stochastic_missed_cleavage_count if fragment else 0,"Blocked_Cleavage_Count":fragment.blocked_cleavage_count if fragment else 0,
                "Blocked_Cleavage_Bond_IDs":_joined(fragment.blocked_cleavage_bond_ids) if fragment else "","Contains_Phosphorothioate":fragment.contains_phosphorothioate if fragment else False,
                "Backbone_Mass_Delta":fragment.backbone_mass_delta if fragment else 0.0,"Shadow_Fragment_Mass":fragment.mass if fragment else "","Applied_To_Formal_Result":False})
    mapped={x for t in transforms for x in t.legacy_ids if x in legacy_ids}
    summary={"Total_Transform_Definitions":len(transforms),"Legacy_Mapped_Transform_Count":len(mapped),"Legacy_Only_Modification_Count":len(legacy_ids-mapped),
        "Valid_Single_Candidate_Count":sum(x["Component_Count"]==1 for x in valid_rows),"Valid_Composite_Candidate_Count":sum(x["Component_Count"]>1 for x in valid_rows),"Invalid_Combination_Count":len(invalid_rows),
        "Slot_Conflict_Count":reason_counts["slot_conflict"],"Missing_Requirement_Count":reason_counts["missing_requirement"],"Parent_Child_Double_Count_Count":reason_counts["parent_child_double_count"]+reason_counts["superseded_component"],
        "Duplicate_Component_Count":reason_counts["duplicate_component"],"Impossible_Oxidation_Count":reason_counts["impossible_oxidation_state"],"Isomer_Group_Count":len(isomer_ids),
        "A_Candidate_Count":base_counts["A"],"G_Candidate_Count":base_counts["G"],"C_Candidate_Count":base_counts["C"],"U_Candidate_Count":base_counts["U"],
        "Phosphorothioate_Candidate_Count":len(backbone_rows),"Formal_Result_Changed":False,"Audit_Mode":audit_mode,"Applied_To_Formal_Result":False,"Formal_Change_Ready":False,
        "Remaining_Risk":"Phase-1 transformations and cleavage rules require structural/evidence curation before formal use."}
    sheets={"Composite_Mod_Candidates":valid_rows,"Composite_Mod_Invalid":invalid_rows,"Composite_Mod_Summary":[summary],"Backbone_Mod_Candidates":backbone_rows,"Cleavage_Block_Audit":cleavage_rows}
    elapsed=time.perf_counter()-started; _,peak=tracemalloc.get_traced_memory(); tracemalloc.stop()
    diagnostics={"Composite_Mod_Audit_Available":True,"Composite_Mod_Valid_Count":len(valid_rows),"Composite_Mod_Invalid_Count":len(invalid_rows),"Composite_Mod_Isomer_Group_Count":len(isomer_ids),
        "Composite_Mod_Legacy_Overlap_Count":len(overlap_ids),"Backbone_Mod_Audit_Available":True,"Phosphorothioate_Candidate_Count":len(backbone_rows),"Cleavage_Block_Audit_Available":True,
        "RNase_T1_Blocked_Cleavage_Count":blocked_counts["RNase_T1"],"Nuclease_P1_Blocked_Cleavage_Count":blocked_counts["Nuclease_P1"],"Composite_Mod_Formal_Change_Ready":False,
        "Composite_Mod_Applied_To_Formal_Result":False,"Backbone_Mod_Applied_To_Formal_Result":False}
    return CompositeAuditResult(sheets,diagnostics,elapsed,peak/(1024*1024))


def append_composite_diagnostics(rows,audit_result:CompositeAuditResult|None):
    import pandas as pd
    is_frame=isinstance(rows,pd.DataFrame); source=rows.to_dict("records") if is_frame else list(rows or [{}])
    values=audit_result.diagnostics if audit_result else {column:(False if "Applied_To_Formal" in column or "Ready" in column else "not_run") for column in DIAGNOSTIC_COLUMNS}
    out=[dict(row,**values) for row in source]
    if is_frame: return pd.DataFrame(out,columns=list(rows.columns)+DIAGNOSTIC_COLUMNS)
    return out
