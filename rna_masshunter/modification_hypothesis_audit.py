"""Data-driven, shadow-only audit of explicit modification position hypotheses."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
import time
import tracemalloc
from pathlib import Path
from typing import Any
from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.enzymes import normalize_enzyme_name
from rna_masshunter.modification_composer import apply_transform_ids
from rna_masshunter.modification_constraints import load_transformations
from rna_masshunter.modification_hypothesis_schema import FALSE_FLAGS,ModificationHypothesisLoadResult,ModificationPositionHypothesis
from rna_masshunter.phosphorothioate_evidence import build_pt_evidence

@dataclass(frozen=True)
class ModificationHypothesisAuditResult:
    sheets:dict[str,list[dict[str,Any]]];metrics:dict[str,Any]

def _pair_positions(pair):return {x.position for x in pair.spec.position_states}
def _pair_mods(pair):return {t for x in pair.spec.position_states for t in x.transform_ids}
def _matches_pair(h,pair):
    position=not h.positions or bool(set(h.positions)&_pair_positions(pair));bond=not h.bonds or pair.spec.bond_id in h.bonds;mods=not h.modification_ids or bool(set(h.modification_ids)&(_pair_mods(pair)|{"phosphorothioate"}))
    if h.hypothesis_type=="nucleoside_modification":return position and bool(set(h.modification_ids)&_pair_mods(pair))
    if h.hypothesis_type in {"backbone_modification","cleavage_behavior"}:return bond
    if h.hypothesis_type=="composite_structure":return bool(h.canonical_backbone_states) and position and bond and set(h.modification_ids)-{"phosphorothioate"}<=_pair_mods(pair)
    if h.hypothesis_type=="absence_hypothesis":return position and bool(set(h.modification_ids)&_pair_mods(pair))
    return False

def _terminal_value(structure,end):
    return getattr(structure,f"terminal_{end}",getattr(structure,end,"inherited"))

def _structure_composition(structure):
    value=ElementalComposition.delta()
    for state in structure.position_states.values():value=value+getattr(state,"elemental_composition_delta",ElementalComposition.delta())
    for bond in structure.bond_states.values():
        if getattr(bond,"state","normal_phosphate")!="normal_phosphate":value=value+getattr(bond,"composition_delta",ElementalComposition.delta())
    return value

def structure_mapping_row(h,structure):
    expected_nuc=dict(h.canonical_nucleoside_states)
    actual_nuc={int(pos):getattr(state,"canonical_structure_id","") for pos,state in structure.position_states.items()}
    expected_backbone=dict(h.canonical_backbone_states)
    actual_backbone={str(bond_id):getattr(state,"state","") for bond_id,state in structure.bond_states.items() if getattr(state,"state","normal_phosphate")!="normal_phosphate"}
    nucleoside_match=expected_nuc==actual_nuc
    backbone_match=(all(actual_backbone.get(k)==v for k,v in expected_backbone.items()) if h.allow_additional_backbone_states else expected_backbone==actual_backbone)
    actual_comp=_structure_composition(structure).canonical_string();composition_match=h.elemental_composition_delta==actual_comp
    expected_terminal={"five_prime":"inherited","three_prime":"inherited"};expected_terminal.update(dict(h.terminal_state))
    actual_terminal={"five_prime":_terminal_value(structure,"five_prime"),"three_prime":_terminal_value(structure,"three_prime")}
    terminal_match=expected_terminal==actual_terminal
    additional_nuc=sum(k not in expected_nuc or expected_nuc.get(k)!=v for k,v in actual_nuc.items())
    additional_backbone=sum(k not in expected_backbone or expected_backbone.get(k)!=v for k,v in actual_backbone.items())
    exact=nucleoside_match and backbone_match and composition_match and terminal_match
    reasons=[]
    if not nucleoside_match:reasons.append("canonical_nucleoside_state_mismatch")
    if not backbone_match:reasons.append("unexpected_backbone_state" if additional_backbone else "canonical_backbone_state_mismatch")
    if not composition_match:reasons.append("elemental_composition_mismatch")
    if not terminal_match:reasons.append("terminal_state_mismatch")
    if exact:status="EXACT_MATCH"
    elif nucleoside_match and not backbone_match:status="NUCLEOSIDE_MATCH_BACKBONE_MISMATCH"
    elif composition_match:status="COMPOSITION_MATCH_STRUCTURE_MISMATCH"
    elif nucleoside_match or backbone_match:status="PARTIAL_MATCH"
    else:status="NO_MATCH"
    return {"Position_Hypothesis_ID":h.hypothesis_id,"Sample_Structure_ID":structure.candidate_id,
        "Canonical_Nucleoside_State_Match":nucleoside_match,"Canonical_Backbone_State_Match":backbone_match,
        "Composition_Match":composition_match,"Terminal_State_Match":terminal_match,"Exact_Structure_Match":exact,
        "Additional_Nucleoside_State_Count":additional_nuc,"Additional_Backbone_State_Count":additional_backbone,
        "Unexpected_Backbone_State":";".join(f"{k}={v}" for k,v in actual_backbone.items() if expected_backbone.get(k)!=v),
        "Unexpected_State_Reason":";".join(reasons),"Mapping_Status":status,"Mapping_Reason":";".join(reasons) or "canonical states, composition, and terminal state match",**FALSE_FLAGS}

def _matches_structure(h,structure):return bool(structure_mapping_row(h,structure)["Exact_Structure_Match"])

def _interpret_oxidation_family(unoxidized_matched,monooxide_matched,unoxidized_specific,monooxide_specific,same_peak,ambiguous,evaluable):
    if not evaluable:return "NOT_EVALUABLE"
    if same_peak or ambiguous:return "STRUCTURE_AMBIGUOUS"
    if unoxidized_specific and monooxide_specific:return "BOTH_SUPPORTED_AS_MIXTURE"
    if unoxidized_specific and not monooxide_matched:return "UNOXIDIZED_SUPPORTED"
    if monooxide_specific and not unoxidized_matched:return "MONOOXIDE_SUPPORTED"
    if unoxidized_matched or monooxide_matched:return "STRUCTURE_AMBIGUOUS"
    return "NEITHER_SUPPORTED"

def _cross_rows(h,pair_ids,cross_sheets):
    rows=(cross_sheets or {}).get("PT_Cross_Run_Summary",[]);out=[]
    for row in rows:
        if row.get("Candidate_ID") not in pair_ids:continue
        if h.bonds and row.get("Bond_ID") not in h.bonds:continue
        if (h.hypothesis_type in {"backbone_modification","cleavage_behavior"} or h.canonical_backbone_states) and row.get("Backbone_State")!="phosphorothioate":continue
        if h.modification_ids and set(h.modification_ids)-{"phosphorothioate"} and not any(x in str(row.get("Nucleoside_State") or "") for x in set(h.modification_ids)-{"phosphorothioate"}):continue
        out.append(row)
    return out

def _interpret(observable,matched,specific,ambiguous,level,alternative_supported,absence=False):
    if not observable:return "NOT_EVALUABLE_WITH_CURRENT_DATA","NOT_EVALUABLE","HYPOTHESIS_NOT_EVALUABLE"
    if absence:
        if specific:return "CONTRADICTED_BY_CURRENT_DATA","CONTRADICTED","HYPOTHESIS_CONTRADICTED"
        return "SUPPORTED_BY_CURRENT_DATA" if not matched else "SUPPORTED_BUT_STRUCTURE_AMBIGUOUS","LOW_DATA_SUPPORT" if not matched else "AMBIGUOUS","HYPOTHESIS_SUPPORTED" if not matched else "HYPOTHESIS_PARTIALLY_SUPPORTED"
    if not matched:
        if alternative_supported:return "DISCOVERY_SUPPORTS_ALTERNATIVE","AMBIGUOUS","DISCOVERY_SUPPORTS_ALTERNATIVE"
        return "NOT_SUPPORTED_IN_CURRENT_DATA","HYPOTHESIS_ONLY","HYPOTHESIS_NOT_SUPPORTED"
    if ambiguous or not specific:return "SUPPORTED_BUT_STRUCTURE_AMBIGUOUS","AMBIGUOUS","HYPOTHESIS_PARTIALLY_SUPPORTED"
    if level>=5:return "SUPPORTED_BY_CURRENT_DATA","STRONG_DATA_SUPPORT","HYPOTHESIS_SUPPORTED"
    if level>=3:return "SUPPORTED_BY_CURRENT_DATA","MODERATE_DATA_SUPPORT","HYPOTHESIS_SUPPORTED"
    return "PARTIALLY_SUPPORTED","LOW_DATA_SUPPORT","HYPOTHESIS_PARTIALLY_SUPPORTED"

def _recommend(h,ms1,ms2,cross,ambiguous):
    choices=[]
    if not ms1:choices.append("higher_MS1_accumulation")
    if not ms2:choices.append("targeted_MS2")
    if cross=="NOT_EVALUABLE":choices.append("independent_T1_digestion")
    if ambiguous:choices.append("isotope_envelope_confirmation")
    if "phosphorothioate" in h.modification_ids:choices.append("normal_PT_control")
    choices.append("synthetic_standard")
    return ";".join(dict.fromkeys(choices))

def _build_oxidation_family_rows(hypotheses,structures,composite_sheets,project_root,config):
    groups=defaultdict(list)
    for h in hypotheses:
        if h.modification_family:groups[(h.target_id,h.positions,h.parent_base,h.modification_family)].append(h)
    if not groups:return []
    support={row.get("Candidate_ID"):row for row in composite_sheets.get("Composite_Support_Summary",[])}
    ms1=defaultdict(list)
    for row in composite_sheets.get("Composite_MS1_Matches",[]):ms1[row.get("Candidate_ID")].append(row)
    transforms=[];slot_schema=None
    if project_root:
        root=Path(project_root);transforms=load_transformations(root/"data/modification_transforms_v2.yaml");slot_schema=root/"data/nucleoside_slots.yaml"
    rows=[]
    for (target_id,positions,parent,slot),members in groups.items():
        family_id=f"{target_id}|{','.join(map(str,positions))}|{slot}"
        oxidation_ids={t.id for t in transforms if t.target_slot==slot}
        common=tuple(x for x in members[0].modification_ids if x not in oxidation_ids)
        variants=[]
        if transforms:
            variants=[(None,common)]+[(t,common+(t.id,)) for t in transforms if t.target_slot==slot]
        else:
            variants=[(None,h.modification_ids) for h in members]
        candidates=[];seen=set()
        for transform,ids in variants:
            if ids in seen:continue
            seen.add(ids);explicit=next((h for h in members if h.modification_ids==ids),None)
            if explicit:
                state_id=dict(explicit.canonical_nucleoside_states).get(positions[0],"");comp=explicit.elemental_composition_delta;mass=explicit.exact_mass_delta;oxidation=explicit.oxidation_state
            elif transforms:
                try:state,result,_=apply_transform_ids(parent,positions[0],ids,transforms,slot_schema)
                except (KeyError,ValueError):continue
                if not result.valid:continue
                state_id=state.canonical_structure_id;comp=state.elemental_composition_delta.canonical_string();mass=state.exact_mass_delta
                oxidation={"thioamide_sulfur":"unoxidized","oxidized_sulfur_1":"monooxide","oxidized_sulfur_2":"dioxide"}.get(getattr(transform,"to_state",None),"precursor")
            else:continue
            matching=[]
            for structure in structures:
                actual={int(pos):getattr(st,"canonical_structure_id","") for pos,st in structure.position_states.items()}
                modified_bonds={k for k,v in structure.bond_states.items() if getattr(v,"state","normal_phosphate")!="normal_phosphate"}
                if actual=={positions[0]:state_id} and not modified_bonds and _structure_composition(structure).canonical_string()==comp:matching.append(structure)
            structure_ids=[x.candidate_id for x in matching];support_rows=[support[x] for x in structure_ids if x in support]
            detail=[r for x in structure_ids for r in ms1[x]]
            position_tokens={str(position) for position in positions}
            informative_detail=[r for r in detail if position_tokens & set(str(r.get("Included_Modified_Positions") or "").split(";"))]
            matched_detail=[r for r in informative_detail if r.get("Match_Status")=="matched"]
            matched=bool(matched_detail) if detail else any(int(r.get("MS1_Matched_Fragment_Count") or 0)>0 for r in support_rows)
            observable=any(r.get("Match_Status")!="not_observable" for r in informative_detail) if detail else any(int(r.get("Observable_Fragment_Count") or 0)>0 for r in support_rows)
            specific=any(r.get("Support_Class")=="unique_composite_support" for r in matched_detail) if detail else any(int(r.get("MS1_Unique_Support_Count") or 0)>0 for r in support_rows)
            physical={f"{r.get('Observed_Scan','')}|{r.get('Observed_RT','')}|{r.get('Observed_mz','')}" for r in matched_detail}
            first=matched_detail[0] if matched_detail else (informative_detail[0] if informative_detail else {})
            candidates.append({"Modification_Family_ID":family_id,"Hypothesis_ID":explicit.hypothesis_id if explicit else "",
                "Candidate_ID":explicit.hypothesis_id if explicit else f"{family_id}|{oxidation}","Oxidation_State":oxidation,
                "Modification_IDs":";".join(ids),"Canonical_Nucleoside_State":state_id,"Composition_Delta":comp,"Exact_Mass_Delta":mass,
                "Theoretical_mz":first.get("Theoretical_mz",""),"Observed_mz":first.get("Observed_mz",""),"Mass_Error_ppm":first.get("Mass_Error_ppm",""),
                "Candidate_Specific":specific,"Observable":observable,"Matched":matched,"Physical_Peak_IDs":";".join(sorted(physical)),
                "Sample_Structure_IDs":";".join(structure_ids),"Exact_Structure_Match":bool(structure_ids),"Unexpected_Backbone_State":"",
                "Chemical_Model_Status":explicit.chemical_model_status if explicit else "hypothesis_shadow_model",
                "Chemical_Model_Note":explicit.chemical_model_note if explicit else "Precursor and oxidation variants are schema-derived shadow states, not confirmed structural assignments.",**FALSE_FLAGS})
        by_state={x["Oxidation_State"]:x for x in candidates};un=by_state.get("unoxidized");mono=by_state.get("monooxide")
        un_matched=bool(un and un["Matched"]);mono_matched=bool(mono and mono["Matched"]);un_specific=bool(un and un["Candidate_Specific"]);mono_specific=bool(mono and mono["Candidate_Specific"])
        same_peak=bool(un and mono and set(filter(None,un["Physical_Peak_IDs"].split(";")))&set(filter(None,mono["Physical_Peak_IDs"].split(";"))))
        ambiguous=bool((un_matched and not un_specific) or (mono_matched and not mono_specific));evaluable=bool(un and mono and un["Observable"] and mono["Observable"])
        interpretation=_interpret_oxidation_family(un_matched,mono_matched,un_specific,mono_specific,same_peak,ambiguous,evaluable)
        delta=(mono["Exact_Mass_Delta"]-un["Exact_Mass_Delta"]) if un and mono else ""
        un_mz=un.get("Theoretical_mz","") if un else "";mono_mz=mono.get("Theoretical_mz","") if mono else ""
        mz_delta=(float(mono_mz)-float(un_mz)) if un_mz not in ("",None) and mono_mz not in ("",None) else ""
        fragment_mapping=getattr(config,"fragment_mapping",{}) or {};min_charge=int(fragment_mapping.get("min_charge",1));max_charge=int(fragment_mapping.get("max_charge",min_charge))
        charge_deltas=";".join(f"z={charge}:{float(delta)/charge:.12f}" for charge in range(min_charge,max_charge+1)) if delta not in ("",None) else ""
        for row in candidates:
            row.update({"Unoxidized_Candidate_ID":un["Candidate_ID"] if un else "","Monooxide_Candidate_ID":mono["Candidate_ID"] if mono else "",
                "Unoxidized_Theoretical_Mass":un["Exact_Mass_Delta"] if un else "","Monooxide_Theoretical_Mass":mono["Exact_Mass_Delta"] if mono else "",
                "Delta_Oxidation_Da":delta,"Oxidation_Delta_Da":delta,"Theoretical_mz_Difference":mz_delta,"Theoretical_mz_Difference_By_Charge":charge_deltas,
                "Unoxidized_Matched":un_matched,"Monooxide_Matched":mono_matched,"Both_Observable":evaluable,
                "Observation_Discriminating":evaluable and not same_peak and not ambiguous,"Observed_State_Comparison":interpretation,
                "Oxidation_Origin_Assessable":False,"Possible_In_Vivo_Oxidation":"unknown","Possible_Ex_Vivo_Oxidation":"unknown",
                "Control_Required":"reducing_or_oxygen_control;fresh_vs_stored_sample;independent_sample_preparation;time_course_oxidation_check;targeted_MS2;synthetic_or_reference_standard",
                "Final_Family_Interpretation":interpretation})
        rows.extend(candidates)
    order={"precursor":0,"unoxidized":1,"monooxide":2,"dioxide":3}
    return sorted(rows,key=lambda x:(x["Modification_Family_ID"],order.get(x["Oxidation_State"],99),x["Candidate_ID"]))

def build_modification_hypothesis_audit(loaded:ModificationHypothesisLoadResult,*,pt_pairs=(),peaks=(),config=None,composite_observation=None,cross_run_sheets=None,audit_level="audit",hypothesis_mode="both",project_root=None):
    started=time.perf_counter();tracemalloc.start();hypotheses=list(loaded.hypotheses)
    configured_enzyme=normalize_enzyme_name(((getattr(config,"digestion",{}) or {}).get("enzyme", "")))
    enzyme_applicable={h.hypothesis_id:not h.enzyme_context or not configured_enzyme or configured_enzyme in h.enzyme_context for h in hypotheses}
    selected_by_h={h.hypothesis_id:[p for p in pt_pairs if _matches_pair(h,p)] if enzyme_applicable[h.hypothesis_id] else [] for h in hypotheses}
    relevant=[]
    for h in hypotheses:
        if not enzyme_applicable[h.hypothesis_id]:continue
        context=selected_by_h[h.hypothesis_id]
        if hypothesis_mode in {"both","discovery"}:
            context=context+[p for p in pt_pairs if (h.bonds and p.spec.bond_id in h.bonds) or (h.positions and set(h.positions)&_pair_positions(p))]
        for pair in context:
            if pair not in relevant:relevant.append(pair)
    evidence,states=build_pt_evidence(relevant,list(peaks or ()),config,audit_level=audit_level,
        include_detail=audit_level=="full",include_compact_states=audit_level!="full") if relevant else ([],[])
    evidence_by_id=defaultdict(list);states_by_id=defaultdict(list)
    for row in evidence:evidence_by_id[row.get("Candidate_ID")].append(row)
    for row in states:states_by_id[row.get("Candidate_ID")].append(row)
    composite_sheets=getattr(composite_observation,"sheets",{}) if composite_observation else {};structures=getattr(composite_observation,"structures",()) if composite_observation else ()
    mapping_rows=[structure_mapping_row(h,structure) for h in hypotheses for structure in structures]
    mapping_by_h=defaultdict(list)
    for row in mapping_rows:mapping_by_h[row["Position_Hypothesis_ID"]].append(row)
    id_rows=[];collision_rows=[];pair_groups=defaultdict(list)
    for pair in pt_pairs:pair_groups[pair.spec.candidate_id].append(pair)
    cross_ids={str(row.get("Candidate_ID") or "") for row in (cross_run_sheets or {}).get("PT_Cross_Run_Summary",[])}
    for h in hypotheses:
        same=[row for row in mapping_by_h[h.hypothesis_id] if row["Sample_Structure_ID"]==h.hypothesis_id]
        same_pairs=pair_groups.get(h.hypothesis_id,[]);pair_consistent=True
        if same_pairs:
            expected_mods=set(h.modification_ids)-{"phosphorothioate"};expected_backbone=dict(h.canonical_backbone_states)
            pair_consistent=all(set(h.positions)==_pair_positions(pair) and expected_mods==_pair_mods(pair) and expected_backbone=={pair.spec.bond_id:"phosphorothioate"} for pair in same_pairs)
        collision=bool(((same and not same[0]["Exact_Structure_Match"]) or not pair_consistent) and not h.alias_of)
        id_rows.append({"Hypothesis_ID":h.hypothesis_id,"Position_Source_Present":True,"Sample_Structure_Source_Present":bool(same),"PT_Pair_Source_Present":bool(same_pairs),
            "Cross_Run_Source_Present":h.hypothesis_id in cross_ids,"Alias_Of":h.alias_of,"Canonical_Structure_Consistent":not collision,
            "ID_Audit_Status":"COLLISION" if collision else "UNIQUE_OR_CONSISTENT",**FALSE_FLAGS})
        if collision:collision_rows.append({"Target_ID":h.target_id,"Hypothesis_ID":h.hypothesis_id,"Valid":False,"Invalid_Reason":"hypothesis_id_structure_collision",
            "Invalid_Detail":same[0]["Mapping_Reason"] if same and not same[0]["Exact_Structure_Match"] else "PT pair canonical state differs",**FALSE_FLAGS})
    for candidate_id,pairs_for_id in pair_groups.items():
        signatures={(tuple(sorted(_pair_positions(pair))),tuple(sorted(_pair_mods(pair))),pair.spec.bond_id) for pair in pairs_for_id}
        if len(signatures)>1:
            collision_rows.append({"Target_ID":"","Hypothesis_ID":candidate_id,"Valid":False,"Invalid_Reason":"hypothesis_id_structure_collision","Invalid_Detail":"PT pair ID maps to multiple structures",**FALSE_FLAGS})
            id_rows.append({"Hypothesis_ID":candidate_id,"Position_Source_Present":False,"Sample_Structure_Source_Present":False,"PT_Pair_Source_Present":True,
                "Cross_Run_Source_Present":candidate_id in cross_ids,"Alias_Of":"","Canonical_Structure_Consistent":False,"ID_Audit_Status":"COLLISION",**FALSE_FLAGS})
    support_by_id={r.get("Candidate_ID"):r for r in composite_sheets.get("Composite_Support_Summary",[])}
    composite_ms1=defaultdict(list);composite_ms2=defaultdict(list)
    for r in composite_sheets.get("Composite_MS1_Matches",[]):composite_ms1[r.get("Candidate_ID")].append(r)
    for r in composite_sheets.get("Composite_MS2_Matches",[]):composite_ms2[r.get("Candidate_ID")].append(r)
    summary=[];details=[];alternatives=[];cross_rows=[]
    for h in hypotheses:
        applicable=enzyme_applicable[h.hypothesis_id]
        pairs=selected_by_h[h.hypothesis_id];pair_ids={p.spec.candidate_id for p in pairs};matching_structures=[s for s in structures if _matches_structure(h,s)] if applicable else [];structure_ids={s.candidate_id for s in matching_structures}
        selected_evidence=[r for cid in pair_ids for r in evidence_by_id[cid]];selected_states=[r for cid in pair_ids for r in states_by_id[cid]];selected_support=[support_by_id[x] for x in structure_ids if x in support_by_id]
        target_backbone=h.hypothesis_type in {"backbone_modification","cleavage_behavior"} or bool(h.canonical_backbone_states)
        target_states=[r for r in selected_states if (not target_backbone or r.get("Backbone_State")=="phosphorothioate") and (not h.modification_ids or not set(h.modification_ids)-{"phosphorothioate"} or any(x in str(r.get("Shared_Nucleoside_States") or "") for x in set(h.modification_ids)-{"phosphorothioate"}))]
        observable=any(bool(r.get("Observable")) for r in target_states) or any(int(r.get("Observable_Fragment_Count") or 0)>0 for r in selected_support)
        matched_states=[r for r in target_states if r.get("Matched")];matched=bool(matched_states) or any(int(r.get("MS1_Matched_Fragment_Count") or 0)>0 for r in selected_support)
        specific=any(r.get("Matched") and int(r.get("Competition_Count") or 0)==0 for r in target_states) or any(int(r.get("MS1_Unique_Support_Count") or 0)>0 for r in selected_support)
        ambiguous=any(r.get("Matched") and int(r.get("Competition_Count") or 0)>0 for r in target_states) or any(int(r.get("MS1_Shared_Support_Count") or 0)+int(r.get("MS1_Isomeric_Unresolved_Count") or 0)>0 for r in selected_support)
        ms2_position=sum(int(r.get("MS2_Position_Informative_Count") or 0) for r in selected_support);ms2_backbone=sum(int(r.get("MS2_Backbone_Informative_Count") or 0) for r in selected_support)
        blocked=sum(int(r.get("Blocked_Cleavage_Match_Count") or 0) for r in selected_support)
        cr=_cross_rows(h,pair_ids,cross_run_sheets);evaluable_runs=max([int(r.get("Evaluable_Run_Count") or 0) for r in cr] or [0]);detected_runs=max([int(r.get("Detected_Run_Count") or 0) for r in cr] or [0]);specific_runs=max([int(r.get("Candidate_Specific_Run_Count") or 0) for r in cr] or [0])
        independent=max([int(r.get("Independent_Digestion_Detected_Count") or 0)+int(r.get("Independent_Preparation_Detected_Count") or 0)+int(r.get("Biological_Replicate_Detected_Count") or 0) for r in cr] or [0]);recurrence=next((r.get("Recurrence_Evidence_Class") for r in cr if str(r.get("Recurrence_Evidence_Class")).startswith("RECURRENT")),"NOT_EVALUABLE" if evaluable_runs<2 else "NO_RECURRENT_SUPPORT")
        level=6 if ms2_backbone else 5 if ms2_position else 4 if independent else 3 if str(recurrence).startswith("RECURRENT") else 2 if specific else 1 if matched else 0
        level_labels={0:"PRIOR_ONLY",1:"SINGLE_RUN_MS1",2:"CANDIDATE_SPECIFIC_MS1",3:"CROSS_RUN_RECURRENT_MS1",4:"INDEPENDENT_RECURRENCE",5:"POSITION_INFORMATIVE_MS2",6:"BOND_LOCALIZING_MS2",7:"CONTROL_OR_PAIRED_EVIDENCE"}
        if any(r.get("Evidence_Class")=="BOTH_PRESENT" for r in selected_evidence):level=7
        all_same=[p for p in pt_pairs if (h.bonds and p.spec.bond_id in h.bonds) or (h.positions and set(h.positions)&_pair_positions(p))]
        alt_pairs=[] if hypothesis_mode=="targeted" else [p for p in all_same if p.spec.candidate_id not in pair_ids];alt_supported=False
        primary_fragments=[p.modified_fragment if target_backbone else p.normal_fragment for p in pairs]
        for p in alt_pairs:
            rows=evidence_by_id.get(p.spec.candidate_id,[]);support=any(r.get("Matched") or "SUPPORT" in str(r.get("Evidence_Class")) for r in rows);alt_supported|=support
            alternative_fragment=p.modified_fragment if target_backbone else p.normal_fragment
            mass_equivalent=any(abs(alternative_fragment.neutral_exact_mass-x.neutral_exact_mass)<=0.001 for x in primary_fragments)
            same_composition=any(alternative_fragment.elemental_composition==x.elemental_composition for x in primary_fragments)
            alternatives.append({"Hypothesis_ID":h.hypothesis_id,"Primary_Candidate_ID":";".join(sorted(pair_ids|structure_ids)),"Alternative_Candidate_ID":p.spec.candidate_id,
                "Alternative_Type":"discovery_candidate","Same_Position":bool(set(h.positions)&_pair_positions(p)),"Same_Bond":p.spec.bond_id in h.bonds,
                "Same_Composition":same_composition,"Mass_Equivalent":mass_equivalent,"Isomeric":same_composition and p.spec.candidate_id not in pair_ids,"Chemically_Exclusive":True,
                "Legacy_or_Composite":"composite","MS1_Discriminating":specific,"MS2_Discriminating":bool(ms2_position or ms2_backbone),"Observed_Support_Class":"supported" if support else "not_supported",**FALSE_FLAGS})
        normal_alternative_count=0
        if target_backbone:
            seen_normal=set()
            for pair in pairs:
                state=";".join(f"{x.position}:{'+'.join(x.transform_ids)}" for x in pair.spec.position_states) or "unmodified"
                key=(pair.spec.bond_id,state)
                if key in seen_normal:continue
                seen_normal.add(key);normal_alternative_count+=1
                erows=evidence_by_id.get(pair.spec.candidate_id,[]);support=any(r.get("Normal_Observed_mz") not in (None,"") for r in erows);alt_supported|=support and recurrence=="NORMAL_DOMINANT"
                alternatives.append({"Hypothesis_ID":h.hypothesis_id,"Primary_Candidate_ID":";".join(sorted(pair_ids|structure_ids)),"Alternative_Candidate_ID":f"{pair.spec.candidate_id}|normal_phosphate",
                    "Alternative_Type":"normal_phosphate_or_stochastic_missed_cleavage","Same_Position":True,"Same_Bond":True,"Same_Composition":False,
                    "Mass_Equivalent":abs(pair.modified_fragment.neutral_exact_mass-pair.normal_fragment.neutral_exact_mass)<=0.001,"Isomeric":False,"Chemically_Exclusive":True,
                    "Legacy_or_Composite":"composite","MS1_Discriminating":specific,"MS2_Discriminating":bool(ms2_backbone),"Observed_Support_Class":"supported" if support else "not_supported",**FALSE_FLAGS})
        interpretation,confidence,hypothesis_result=_interpret(observable,matched,specific,ambiguous,level,alt_supported,h.hypothesis_type=="absence_hypothesis")
        localization="BOND_LOCALIZED" if ms2_backbone else "EXACT_POSITION_LOCALIZED" if ms2_position and h.exact_position else "POSITION_RANGE_SUPPORTED" if matched and h.position_range else "BOND_UNRESOLVED" if h.bonds else "POSITION_UNRESOLVED"
        observed="supporting_observation" if matched else "no_observation" if observable else "not_observable"
        reason=f"observable={observable}; matched={matched}; candidate_specific={specific}; ambiguous={ambiguous}; MS2_position={ms2_position}; MS2_backbone={ms2_backbone}; recurrence={recurrence}"
        not_applicable_reason=""
        if not applicable:
            observable=matched=specific=ambiguous=False;level=0
            interpretation=confidence=hypothesis_result=localization="NOT_APPLICABLE_TO_ENZYME";observed="not_applicable"
            not_applicable_reason="RNase_T1_hypothesis_on_Nuclease_P1_run" if configured_enzyme=="Nuclease_P1" else "hypothesis_enzyme_context_mismatch"
            reason=f"configured_enzyme={configured_enzyme}; required_enzyme={';'.join(h.enzyme_context)}; {not_applicable_reason}"
        source="hypothesis_and_discovery" if any(p.spec.search_mode=="discovery" for p in pairs) else "hypothesis"
        hypothesis_alternatives=[x for x in alternatives if x["Hypothesis_ID"]==h.hypothesis_id]
        supported_alternatives=[x for x in hypothesis_alternatives if x["Observed_Support_Class"]=="supported"]
        summary.append({"Audit_Level":audit_level,"Target_ID":h.target_id,"Hypothesis_ID":h.hypothesis_id,"Hypothesis_Type":h.hypothesis_type,"Priority":h.priority,
            "Prior_Status":h.prior_status,"Prior_Strength":h.prior_strength,"Prior_Hypothesis_Status":h.prior_status,"Prior_Hypothesis_Strength":h.prior_strength,
            "Evidence_Basis":";".join(h.evidence_basis),"Position_or_Bond":";".join(map(str,h.positions)) or ";".join(h.bonds),
            "Parent_Base":h.parent_base,"Modification_ID":";".join(h.modification_ids),"Component_Domains":";".join(h.component_domains),"Composition_Delta":h.elemental_composition_delta,"Exact_Mass_Delta":h.exact_mass_delta,
            "Chemical_Validity":h.chemical_validity,"Position_Compatibility":h.position_compatibility,"Pathway_Compatibility":h.pathway_compatibility,
            "Modification_Family":h.modification_family,"Oxidation_State":h.oxidation_state,"Chemical_Model_Status":h.chemical_model_status,"Chemical_Model_Note":h.chemical_model_note,
            "Exact_Structure_Match":any(x["Exact_Structure_Match"] for x in mapping_by_h[h.hypothesis_id]),
            "Unexpected_Backbone_State":";".join(x["Unexpected_Backbone_State"] for x in mapping_by_h[h.hypothesis_id] if x["Unexpected_Backbone_State"]),
            "Candidate_Source":source,"Hypothesis_Mode":hypothesis_mode,"Configured_Enzyme":configured_enzyme,"Enzyme_Context_Applicable":applicable,"Evidence_Not_Applicable_Reason":not_applicable_reason,"Observable":observable,"Candidate_Observable":observable,"Expected_mz_In_Range":applicable and any(r.get("Theoretical_mz") not in (None,"") and (getattr(config,"reconstruction",{}) or {}).get("mz_min",0)<=float(r.get("Theoretical_mz"))<=(getattr(config,"reconstruction",{}) or {}).get("mz_max",float("inf")) for r in target_states),"Expected_Charge_In_Range":applicable and bool(target_states),"Eligible_MS1_Spectra_Present":bool(peaks),
            "Precursor_Compatible_MS2_Present":bool(ms2_position or ms2_backbone),"Expected_Fragment_Generated":bool(pairs or matching_structures),"Detection_Sensitivity_Adequate":bool(peaks),"Cross_Run_Evaluable":evaluable_runs>=2,
            "Observed_Evidence_Status":observed,"Observed_Evidence_Level":level,"Highest_Evidence_Level":level,"Evidence_Level_Label":level_labels[level],"Evidence_Level_Reason":reason,"Hypothesis_Result_Class":hypothesis_result,
            "Hypothesis_Confidence":confidence,"Cross_Run_Status":recurrence,"Position_Localization_Status":localization,
            "Alternative_Hypothesis_IDs":";".join(x["Alternative_Candidate_ID"] for x in hypothesis_alternatives),"Alternative_Candidate_Count":len(hypothesis_alternatives),
            "Mass_Equivalent_Alternative_Count":sum(bool(x["Mass_Equivalent"]) for x in hypothesis_alternatives),"Isomeric_Alternative_Count":sum(bool(x["Isomeric"]) for x in hypothesis_alternatives),
            "Chemically_Exclusive_Alternative_Count":sum(bool(x["Chemically_Exclusive"]) for x in hypothesis_alternatives),"Observation_Discriminating":specific and not ambiguous,
            "MS2_Discriminating":bool(ms2_position or ms2_backbone),"Best_Alternative_Explanation":supported_alternatives[0]["Alternative_Candidate_ID"] if supported_alternatives else "",
            "MS1_Data_Sufficient":bool(peaks) and observable,"MS2_Data_Sufficient":bool(ms2_position or ms2_backbone),"Cross_Run_Data_Sufficient":evaluable_runs>=2,"Control_Data_Sufficient":level>=7,
            "Overall_Data_Sufficiency":"not_applicable" if not applicable else "sufficient" if level>=5 or independent else "limited","Final_Shadow_Interpretation":interpretation,"Evidence_Reason":reason,
            "Recommended_Next_Evidence":"evaluate_with_RNase_T1_run" if not applicable else _recommend(h,matched,bool(ms2_position or ms2_backbone),"EVALUABLE" if evaluable_runs>=2 else "NOT_EVALUABLE",ambiguous),**FALSE_FLAGS})
        cross_rows.append({"Hypothesis_ID":h.hypothesis_id,"Evaluable_Run_Count":evaluable_runs,"Detected_Run_Count":detected_runs,"Candidate_Specific_Run_Count":specific_runs,"Independent_Run_Count":independent,
            "Detection_Rate":detected_runs/evaluable_runs if evaluable_runs else 0.0,"Independent_Detection_Rate":independent/evaluable_runs if evaluable_runs else 0.0,"Recurrence_Evidence_Class":recurrence,
            "Mass_Error_Consistency":next((r.get("Error_Sign_Consistency") for r in cr),"NOT_EVALUABLE"),"Charge_Consistency":next((r.get("Charge_Consistency_Status") for r in cr),"NOT_EVALUABLE"),
            "RT_Consistency":next((r.get("RT_Consistency_Status") for r in cr),"NOT_EVALUABLE"),"Cross_Run_Status":"EVALUABLE" if evaluable_runs>=2 else "NOT_EVALUABLE",**FALSE_FLAGS})
        if audit_level=="full":
            for r in target_states:
                details.append({"Run_ID":"current_run","Hypothesis_ID":h.hypothesis_id,"Candidate_ID":r.get("Candidate_ID"),"Fragment_ID":f"{r.get('Candidate_ID')}|{r.get('Fragment_Start')}_{r.get('Fragment_End')}",
                    "Fragment_Start":r.get("Fragment_Start"),"Fragment_End":r.get("Fragment_End"),"Included_Positions":r.get("Shared_Modified_Positions"),"Included_Bonds":r.get("Bond_ID"),"Charge":r.get("Charge"),
                    "Elemental_Composition":r.get("Elemental_Composition"),"Neutral_Mass":r.get("Neutral_Mass"),"Theoretical_mz":r.get("Theoretical_mz"),"Observed_mz":r.get("Nearest_Observed_mz") if r.get("Matched") else "",
                    "Mass_Error_ppm":r.get("Mass_Error_ppm") if r.get("Matched") else "","Intensity":r.get("Intensity") if r.get("Matched") else "","Peak_Tier":"","Scan":r.get("Scan") if r.get("Matched") else "","RT":r.get("RT") if r.get("Matched") else "",
                    "Physical_Peak_ID":r.get("Physical_Peak_ID") if r.get("Matched") else "","Candidate_Specific":r.get("Matched") and int(r.get("Competition_Count") or 0)==0,"Competition_Count":r.get("Competition_Count"),
                    "Position_Informative":False,"Backbone_Informative":False,"Evidence_Type":"MS1","Evidence_Status":"supporting_observation" if r.get("Matched") else "no_observation" if r.get("Observable") else "not_observable","Evidence_Level":2 if r.get("Matched") and int(r.get("Competition_Count") or 0)==0 else 1 if r.get("Matched") else 0,**FALSE_FLAGS})
            for cid in structure_ids:
                for r in composite_ms1.get(cid,[]):
                    matched_row=r.get("Match_Status")=="matched";candidate_specific=r.get("Support_Class")=="unique_composite_support"
                    details.append({"Run_ID":"current_run","Hypothesis_ID":h.hypothesis_id,"Candidate_ID":cid,"Fragment_ID":r.get("Fragment_ID"),"Fragment_Start":r.get("Start_Position"),"Fragment_End":r.get("End_Position"),
                        "Included_Positions":r.get("Included_Modified_Positions"),"Included_Bonds":r.get("Included_Backbone_Bonds"),"Charge":r.get("Charge"),"Elemental_Composition":"","Neutral_Mass":r.get("Neutral_Exact_Mass"),
                        "Theoretical_mz":r.get("Theoretical_mz"),"Observed_mz":r.get("Observed_mz"),"Mass_Error_ppm":r.get("Mass_Error_ppm"),"Intensity":r.get("Observed_Intensity"),"Peak_Tier":"","Scan":r.get("Observed_Scan"),"RT":r.get("Observed_RT"),
                        "Physical_Peak_ID":f"{r.get('Observed_Scan','')}|{r.get('Observed_RT','')}|{r.get('Observed_mz','')}" if matched_row else "","Candidate_Specific":candidate_specific,"Competition_Count":0 if candidate_specific else 1 if matched_row else 0,
                        "Position_Informative":False,"Backbone_Informative":False,"Evidence_Type":"MS1","Evidence_Status":"supporting_observation" if matched_row else r.get("Match_Status"),"Evidence_Level":2 if candidate_specific else 1 if matched_row else 0,**FALSE_FLAGS})
                for r in composite_ms2.get(cid,[]):
                    details.append({"Run_ID":"current_run","Hypothesis_ID":h.hypothesis_id,"Candidate_ID":cid,"Fragment_ID":r.get("Spectrum_ID"),"Fragment_Start":"","Fragment_End":"","Included_Positions":r.get("Included_Modified_Positions"),"Included_Bonds":r.get("Included_Backbone_Bonds"),
                        "Charge":r.get("Precursor_Charge"),"Elemental_Composition":"","Neutral_Mass":r.get("Theoretical_Neutral_Mass"),"Theoretical_mz":r.get("Theoretical_mz"),"Observed_mz":r.get("Observed_mz"),"Mass_Error_ppm":r.get("Mass_Error_ppm"),
                        "Intensity":r.get("Observed_Intensity"),"Peak_Tier":"","Scan":r.get("Spectrum_ID"),"RT":"","Physical_Peak_ID":f"{r.get('Spectrum_ID')}|{r.get('Observed_mz')}","Candidate_Specific":r.get("Candidate_Discriminating"),"Competition_Count":0 if r.get("Candidate_Discriminating") else 1,
                        "Position_Informative":r.get("Position_Informative"),"Backbone_Informative":r.get("Backbone_Informative"),"Evidence_Type":"MS2","Evidence_Status":"supporting_observation","Evidence_Level":6 if r.get("Backbone_Informative") else 5 if r.get("Position_Informative") else 2,**FALSE_FLAGS})
    current,peak=tracemalloc.get_traced_memory();tracemalloc.stop();metrics={"Hypothesis_Audit_Runtime":time.perf_counter()-started,"Hypothesis_Audit_Tracemalloc_Peak_MiB":peak/1024/1024,
        "Summary_Row_Count":len(summary),"Detail_Row_Count":len(details) if audit_level=="full" else 0,"Alternative_Row_Count":len(alternatives) if audit_level=="full" else 0,"Cross_Run_Row_Count":len(cross_rows)}
    family_rows=_build_oxidation_family_rows(hypotheses,structures,composite_sheets,project_root,config)
    family_by_h={row.get("Hypothesis_ID"):row for row in family_rows if row.get("Hypothesis_ID")}
    for row in summary:
        family=family_by_h.get(row["Hypothesis_ID"],{})
        row.update({k:family.get(k,"") for k in ("Unoxidized_Candidate_ID","Monooxide_Candidate_ID","Oxidation_Delta_Da","Observed_State_Comparison","Oxidation_Origin_Assessable","Final_Family_Interpretation")})
    sheets={"Mod_Hypothesis_Summary":summary,"Mod_Hypothesis_Cross_Run":cross_rows,
        "Mod_Hypothesis_Invalid":list(loaded.invalid_rows)+collision_rows,"Mod_Hypothesis_Structure_Map":mapping_rows,
        "Mod_Hypothesis_ID_Audit":id_rows,"Mod_Oxidation_Family":family_rows}
    if audit_level=="full":sheets.update({"Mod_Hypothesis_Detail":details,"Mod_Hypothesis_Alternatives":alternatives})
    return ModificationHypothesisAuditResult(sheets,metrics)
