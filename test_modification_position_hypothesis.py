from __future__ import annotations
from copy import deepcopy
from pathlib import Path
import hashlib
import yaml
import pytest
from rna_masshunter.config import load_config
from rna_masshunter.modification_hypothesis_schema import load_modification_position_hypotheses
from rna_masshunter.modification_hypothesis_audit import _interpret,_interpret_oxidation_family,_matches_pair,build_modification_hypothesis_audit,structure_mapping_row
from rna_masshunter.models import Peak
from rna_masshunter.masses import mz_from_neutral_mass

ROOT=Path(__file__).resolve().parent;CFG=load_config(ROOT/"config.yaml");SEQ=CFG.sequence["sequence"]
def target(hypotheses):return {"target_id":"T","sequence_id":CFG.sequence["name"],"sequence_name":CFG.sequence["name"],"sequence_length":len(SEQ),"sequence_sha256":hashlib.sha256(SEQ.encode()).hexdigest(),"organism":CFG.organism["species"],"rule_set":CFG.organism["rule_set"],"hypotheses":hypotheses}
def hyp(hid="H",kind="nucleoside_modification",**extra):
 d={"hypothesis_id":hid,"hypothesis_type":kind,"priority":"normal","prior_status":"suspected","prior_strength":"moderate","evidence_basis":["manual_expert_hypothesis"]};d.update(extra);return d
def load(tmp_path,hypotheses,mutate=None):
 payload={"schema_version":1,"targets":[target(hypotheses)]}
 if mutate:mutate(payload)
 path=tmp_path/"h.yaml";path.write_text(yaml.safe_dump(payload));return load_modification_position_hypotheses(path,project_root=ROOT,sequence=SEQ,sequence_id=CFG.sequence["name"],sequence_name=CFG.sequence["name"],organism=CFG.organism["species"],rule_set=CFG.organism["rule_set"])

def test_example_is_generic_and_valid():
 text=(ROOT/"data/modification_position_hypotheses.example.yaml").read_text(encoding="utf-8")
 assert "Mac_" not in text and "Methanosarcina" not in text
 r=load_modification_position_hypotheses(ROOT/"data/modification_position_hypotheses.example.yaml",project_root=ROOT,sequence="ACGU",sequence_id="example_rna",sequence_name="example_rna",organism="example_organism",rule_set="example_rule_set")
 assert len(r.hypotheses)==1 and not r.invalid_rows
@pytest.mark.parametrize("kind,extra",[
 ("nucleoside_modification",{"position":10,"parent_base":"G","modification_id":"m22G"}),
 ("backbone_modification",{"bond_id":"10_11","modification_id":"phosphorothioate","enzyme_context":["RNase_T1"]}),
 ("composite_structure",{"position":10,"parent_base":"G","modification_ids":["m22G","phosphorothioate"],"bond_id":"10_11"}),
 ("cleavage_behavior",{"bond_id":"10_11","enzyme_context":["RNase_T1"]}),
 ("absence_hypothesis",{"position":10,"parent_base":"G","modification_id":"m22G"}),
 ("terminal_state",{"terminal_state":{"five_prime":"inherited"}}),
])
def test_hypothesis_types(tmp_path,kind,extra):assert len(load(tmp_path,[hyp(kind=kind,**extra)]).hypotheses)==1
def test_duplicate_id(tmp_path):assert "duplicate_hypothesis_id" in load(tmp_path,[hyp(position=10,parent_base="G",modification_id="m22G"),hyp(position=10,parent_base="G",modification_id="m22G")]).invalid_rows[0]["Invalid_Reason"]
@pytest.mark.parametrize("field,value,reason",[("organism","wrong","target_organism_mismatch"),("sequence_sha256","bad","target_sequence_sha256_mismatch"),("sequence_length",77,"target_sequence_length_mismatch"),("sequence_id","wrong","target_sequence_id_mismatch")])
def test_target_guard(tmp_path,field,value,reason):
 def mutate(p):p["targets"][0][field]=value
 r=load(tmp_path,[hyp(position=10,parent_base="G",modification_id="m22G")],mutate);assert reason in r.invalid_rows[0]["Invalid_Reason"]
def test_missing_target_identity(tmp_path):
 def mutate(p):p["targets"][0].pop("sequence_sha256")
 assert "target_sequence_sha256_mismatch" in load(tmp_path,[hyp(position=10,parent_base="G",modification_id="m22G")],mutate).invalid_rows[0]["Invalid_Reason"]
@pytest.mark.parametrize("payload,reason",[
 ({"position":0,"parent_base":"G","modification_id":"m22G"},"position_out_of_range"),
 ({"position":10,"parent_base":"U","modification_id":"m22G"},"parent_base_mismatch"),
 ({"position":10,"parent_base":"G","modification_id":"missing"},"transform_not_found"),
 ({"bond_id":"78_79","modification_id":"phosphorothioate"},"invalid_bond"),
])
def test_invalid_coordinates_and_transform(tmp_path,payload,reason):assert reason in load(tmp_path,[hyp(kind="backbone_modification" if "bond_id" in payload else "nucleoside_modification",**payload)]).invalid_rows[0]["Invalid_Reason"]
def test_invalid_prior_and_evidence(tmp_path):
 r=load(tmp_path,[hyp(position=10,parent_base="G",modification_id="m22G",prior_status="certain",evidence_basis=["invented"])]);assert "invalid_prior_status" in r.invalid_rows[0]["Invalid_Reason"] and "invalid_evidence_code" in r.invalid_rows[0]["Invalid_Reason"]
def test_exact_range_candidate_positions_and_bonds(tmp_path):
 hs=[hyp("exact",position=10,parent_base="G",modification_id="m22G"),hyp("range",position_range={"start":8,"end":12},parent_base="G",modification_id="m22G"),hyp("candidates",candidate_positions=[10,12],parent_base="G",modification_id="m22G"),hyp("bonds",kind="backbone_modification",candidate_bonds=["10_11","11_12"],modification_id="phosphorothioate")]
 r=load(tmp_path,hs);assert {x.hypothesis_id for x in r.hypotheses}=={"exact","range","candidates","bonds"}

def test_prior_only_never_supports():assert _interpret(True,False,False,False,0,False)[0]=="NOT_SUPPORTED_IN_CURRENT_DATA"
def test_not_observable_is_not_contradiction():assert _interpret(False,False,False,False,0,True)[0]=="NOT_EVALUABLE_WITH_CURRENT_DATA"
def test_candidate_specific_and_ambiguous_decisions():
 assert _interpret(True,True,True,False,2,False)[0]=="PARTIALLY_SUPPORTED"
 assert _interpret(True,True,False,True,1,False)[0]=="SUPPORTED_BUT_STRUCTURE_AMBIGUOUS"
def test_independent_or_ms2_support():assert _interpret(True,True,True,False,4,False)[0]=="SUPPORTED_BY_CURRENT_DATA"
def test_discovery_alternative():assert _interpret(True,False,False,False,0,True)[0]=="DISCOVERY_SUPPORTS_ALTERNATIVE"
def test_absence_contradicted_only_by_specific_observation():
 assert _interpret(True,True,True,False,2,False,True)[0]=="CONTRADICTED_BY_CURRENT_DATA"
 assert _interpret(True,False,False,False,0,False,True)[0]=="SUPPORTED_BY_CURRENT_DATA"

def test_data_driven_audit_and_prior_separation(tmp_path):
 from dataclasses import replace
 from test_pt_paired_evidence import pair,cfg
 p=pair(("m22G",));p=replace(p,spec=replace(p.spec,candidate_id="C",search_mode="hypothesis_driven"))
 r=load(tmp_path,[hyp("weak",position=10,parent_base="G",modification_id="m22G",prior_strength="weak"),hyp("strong",position=10,parent_base="G",modification_id="m22G",prior_strength="strong")])
 mz=mz_from_neutral_mass(p.modified_fragment.neutral_exact_mass,1,"negative");config=cfg(max_charge=1)
 a=build_modification_hypothesis_audit(r,pt_pairs=[p],peaks=[Peak(mz,2000,1.0,"scan")],config=config,audit_level="full")
 rows=a.sheets["Mod_Hypothesis_Summary"];assert len(rows)==2 and rows[0]["Observed_Evidence_Status"]==rows[1]["Observed_Evidence_Status"]
 assert all(x["Applied_To_Formal_Result"] is False and x["Formal_Result_Changed"] is False for values in a.sheets.values() for x in values)
 assert a.sheets["Mod_Hypothesis_Detail"] and a.metrics["Detail_Row_Count"]>0

def test_cross_run_integration(tmp_path):
 from test_pt_paired_evidence import pair,cfg
 p=pair(("m22G",));r=load(tmp_path,[hyp(position=10,parent_base="G",modification_id="m22G")])
 cross={"PT_Cross_Run_Summary":[{"Candidate_ID":"C","Nucleoside_State":"10:m22G","Evaluable_Run_Count":3,"Detected_Run_Count":2,"Candidate_Specific_Run_Count":2,"Independent_Digestion_Detected_Count":2,"Independent_Preparation_Detected_Count":0,"Biological_Replicate_Detected_Count":0,"Recurrence_Evidence_Class":"RECURRENT_INDEPENDENT_PT_MS1_SUPPORT","Error_Sign_Consistency":"CONSISTENT","Charge_Consistency_Status":"CONSISTENT","RT_Consistency_Status":"RAW_CONSISTENT"}]}
 a=build_modification_hypothesis_audit(r,pt_pairs=[p],peaks=[],config=cfg(max_charge=1),cross_run_sheets=cross,audit_level="audit")
 row=a.sheets["Mod_Hypothesis_Cross_Run"][0];assert row["Independent_Run_Count"]==2 and row["Cross_Run_Status"]=="EVALUABLE"

def test_standard_cli_does_not_require_file():
 from main import parse_args
 a=parse_args(["--audit-level","standard","--position-hypotheses","/missing"]);assert a.position_hypotheses=="/missing"

def test_transform_requires_and_parent_derived_constraints(tmp_path):
    r=load(tmp_path,[hyp(position=37,parent_base="U",modification_id="side_chain_thioamide_oxo1")])
    assert any(code in r.invalid_rows[0]["Invalid_Reason"] for code in ("missing_requirement","from_state_mismatch","impossible_oxidation_state"))
def test_unknown_enzyme_and_bond_base_guard(tmp_path):
    r=load(tmp_path,[hyp("e",kind="backbone_modification",bond_id="10_11",left_base="U",modification_id="phosphorothioate",enzyme_context=["missing"]),hyp("b",kind="backbone_modification",bond_id="10_11",left_base="U",modification_id="phosphorothioate")])
    reasons=";".join(x["Invalid_Reason"] for x in r.invalid_rows);assert "unknown_enzyme" in reasons and "left_base_mismatch" in reasons
def test_invalid_terminal_state(tmp_path):
    r=load(tmp_path,[hyp(kind="terminal_state",terminal_state={"five_prime":"invented"})]);assert "invalid_terminal_state" in r.invalid_rows[0]["Invalid_Reason"]

def test_position_and_bond_localizing_ms2_levels(tmp_path):
    from types import SimpleNamespace
    from test_pt_paired_evidence import cfg
    expected=load(tmp_path,[hyp(position=10,parent_base="G",modification_id="m22G")]).hypotheses[0]
    from rna_masshunter.elemental_composition import ElementalComposition
    state=SimpleNamespace(applied_transform_ids=("m22G",),parent_base="G",canonical_structure_id=dict(expected.canonical_nucleoside_states)[10],elemental_composition_delta=ElementalComposition.delta({"C":2,"H":4}))
    normal_bond=SimpleNamespace(state="normal_phosphate",composition_delta=ElementalComposition.delta())
    structure=SimpleNamespace(candidate_id="S",position_states={10:state},bond_states={"10_11":normal_bond},five_prime="inherited",three_prime="inherited")
    sheets={"Composite_Support_Summary":[{"Candidate_ID":"S","Observable_Fragment_Count":1,"MS1_Matched_Fragment_Count":1,"MS1_Unique_Support_Count":1,"MS1_Shared_Support_Count":0,"MS1_Isomeric_Unresolved_Count":0,"MS2_Position_Informative_Count":1,"MS2_Backbone_Informative_Count":0,"Blocked_Cleavage_Match_Count":0}]}
    obs=SimpleNamespace(structures=(structure,),sheets=sheets)
    r=load(tmp_path,[hyp(position=10,parent_base="G",modification_id="m22G")]);a=build_modification_hypothesis_audit(r,config=cfg(),composite_observation=obs,audit_level="audit")
    row=a.sheets["Mod_Hypothesis_Summary"][0];assert row["Highest_Evidence_Level"]==5 and row["Position_Localization_Status"]=="EXACT_POSITION_LOCALIZED"
    ptbond=SimpleNamespace(state="phosphorothioate",composition_delta=ElementalComposition.delta({"O":-1,"S":1}));structure=SimpleNamespace(candidate_id="B",position_states={},bond_states={"10_11":ptbond},five_prime="inherited",three_prime="inherited")
    sheets={"Composite_Support_Summary":[{"Candidate_ID":"B","Observable_Fragment_Count":1,"MS1_Matched_Fragment_Count":1,"MS1_Unique_Support_Count":1,"MS1_Shared_Support_Count":0,"MS1_Isomeric_Unresolved_Count":0,"MS2_Position_Informative_Count":0,"MS2_Backbone_Informative_Count":1,"Blocked_Cleavage_Match_Count":1}]}
    obs=SimpleNamespace(structures=(structure,),sheets=sheets);r=load(tmp_path,[hyp(kind="backbone_modification",bond_id="10_11",modification_id="phosphorothioate")]);a=build_modification_hypothesis_audit(r,config=cfg(),composite_observation=obs,audit_level="audit")
    row=a.sheets["Mod_Hypothesis_Summary"][0];assert row["Highest_Evidence_Level"]==6 and row["Position_Localization_Status"]=="BOND_LOCALIZED"

def test_priority_and_audit_do_not_mutate_tolerance_candidates_or_formal_rows(tmp_path):
    from test_pt_paired_evidence import pair,cfg
    pairs=[pair(("m22G",)),pair(())];before_pairs=deepcopy(pairs);formal=[{"Final_Score":7,"Final_Confidence":"High","Rank":1}];before=deepcopy(formal);config=cfg();tol=config.fragment_mapping["mz_tolerance_ppm"]
    r=load(tmp_path,[hyp("critical",position=10,parent_base="G",modification_id="m22G",priority="critical"),hyp("low",position=10,parent_base="G",modification_id="m22G",priority="low")])
    a=build_modification_hypothesis_audit(r,pt_pairs=pairs,peaks=[],config=config,audit_level="audit")
    assert config.fragment_mapping["mz_tolerance_ppm"]==tol and formal==before and len(pairs)==len(before_pairs)
    assert [x["Hypothesis_ID"] for x in a.sheets["Mod_Hypothesis_Summary"]]==["critical","low"]
def test_alternative_rows_have_mass_and_isomer_audit(tmp_path):
    from test_pt_paired_evidence import pair,cfg
    from dataclasses import replace
    pairs=[]
    for cid,item in (("m22",pair(("m22G",))),("m1",pair(("m1G",))),("unmod",pair(()))):pairs.append(replace(item,spec=replace(item.spec,candidate_id=cid)))
    r=load(tmp_path,[hyp(position=10,parent_base="G",modification_id="m22G")])
    a=build_modification_hypothesis_audit(r,pt_pairs=pairs,peaks=[],config=cfg(),audit_level="full")
    assert a.sheets["Mod_Hypothesis_Alternatives"] and all("Mass_Equivalent" in x and "Isomeric" in x for x in a.sheets["Mod_Hypothesis_Alternatives"])
def test_audit_mode_does_not_emit_peak_or_alternative_detail(tmp_path):
    from test_pt_paired_evidence import pair,cfg
    r=load(tmp_path,[hyp(position=10,parent_base="G",modification_id="m22G")]);a=build_modification_hypothesis_audit(r,pt_pairs=[pair(("m22G",))],peaks=[],config=cfg(),audit_level="audit")
    assert "Mod_Hypothesis_Detail" not in a.sheets and "Mod_Hypothesis_Alternatives" not in a.sheets and a.metrics["Detail_Row_Count"]==0


def test_composite_component_domain_validation(tmp_path):
    nuc=load(tmp_path,[hyp(kind="composite_structure",position=37,parent_base="U",modification_ids=["s2U","cnm5U","side_chain_thioamide"])])
    assert len(nuc.hypotheses)==1 and nuc.hypotheses[0].component_domains==("nucleoside",) and not nuc.hypotheses[0].bonds
    mixed=load(tmp_path,[hyp(kind="composite_structure",position=10,parent_base="G",modification_ids=["m22G","phosphorothioate"],bond_id="10_11")])
    assert len(mixed.hypotheses)==1 and mixed.hypotheses[0].component_domains==("nucleoside","backbone")
    single=load(tmp_path,[hyp(kind="composite_structure",position=10,parent_base="G",modification_id="m22G")])
    assert "composite_requires_multiple_components" in single.invalid_rows[0]["Invalid_Reason"]
    no_position=load(tmp_path,[hyp(kind="composite_structure",modification_ids=["m22G","m1G"])])
    assert "missing_position" in no_position.invalid_rows[0]["Invalid_Reason"]
    no_bond=load(tmp_path,[hyp(kind="composite_structure",position=10,parent_base="G",modification_ids=["m22G","phosphorothioate"])])
    assert "missing_bond" in no_bond.invalid_rows[0]["Invalid_Reason"]
    backbone_only=load(tmp_path,[hyp(kind="composite_structure",candidate_bonds=["10_11","11_12"],modification_id="phosphorothioate")])
    assert len(backbone_only.hypotheses)==1 and backbone_only.hypotheses[0].component_domains==("backbone",)

def test_u37_chemistry_and_no_backbone_pt(tmp_path):
    r=load(tmp_path,[
        hyp("U",kind="composite_structure",position=37,parent_base="U",modification_ids=["s2U","cnm5U","side_chain_thioamide"]),
        hyp("O",kind="composite_structure",position=37,parent_base="U",modification_ids=["s2U","cnm5U","side_chain_thioamide_oxo1"])])
    by={x.hypothesis_id:x for x in r.hypotheses}
    assert by["U"].elemental_composition_delta=="C2H2N2O-2S2"
    assert by["O"].elemental_composition_delta=="C2H2N2O-1S2"
    from rna_masshunter.elemental_composition import ElementalComposition
    assert by["O"].exact_mass_delta-by["U"].exact_mass_delta==pytest.approx(ElementalComposition.delta({"O":1}).exact_mass)
    assert by["O"].canonical_backbone_states==() and "phosphorothioate" not in by["O"].modification_ids
    double=load(tmp_path,[hyp(kind="composite_structure",position=37,parent_base="U",modification_ids=["s2U","cnm5U","side_chain_thioamide","side_chain_thioamide_oxo1"])])
    assert any(x in double.invalid_rows[0]["Invalid_Reason"] for x in ("parent_child_double_count","superseded_component"))

def test_strict_mapping_rejects_extra_backbone_and_same_composition_structure(tmp_path):
    from types import SimpleNamespace
    from rna_masshunter.elemental_composition import ElementalComposition
    h=load(tmp_path,[hyp(kind="composite_structure",position=37,parent_base="U",modification_ids=["s2U","cnm5U","side_chain_thioamide_oxo1"])]).hypotheses[0]
    state=SimpleNamespace(canonical_structure_id=dict(h.canonical_nucleoside_states)[37],elemental_composition_delta=ElementalComposition.delta({"C":2,"H":2,"N":2,"O":-1,"S":2}))
    normal=SimpleNamespace(state="normal_phosphate",composition_delta=ElementalComposition.delta())
    exact=SimpleNamespace(candidate_id="exact",position_states={37:state},bond_states={},five_prime="inherited",three_prime="inherited")
    assert structure_mapping_row(h,exact)["Exact_Structure_Match"]
    pt=SimpleNamespace(state="phosphorothioate",composition_delta=ElementalComposition.delta({"O":-1,"S":1}))
    extra=SimpleNamespace(candidate_id="extra",position_states={37:state},bond_states={"37_38":pt},five_prime="inherited",three_prime="inherited")
    row=structure_mapping_row(h,extra);assert not row["Exact_Structure_Match"] and row["Additional_Backbone_State_Count"]==1
    other_state=SimpleNamespace(canonical_structure_id="different-slot-state",elemental_composition_delta=state.elemental_composition_delta)
    isomer=SimpleNamespace(candidate_id="isomer",position_states={37:other_state},bond_states={},five_prime="inherited",three_prime="inherited")
    row=structure_mapping_row(h,isomer);assert row["Composition_Match"] and not row["Exact_Structure_Match"] and row["Mapping_Status"]=="COMPOSITION_MATCH_STRUCTURE_MISMATCH"

def test_same_id_different_structure_collision(tmp_path):
    from types import SimpleNamespace
    from rna_masshunter.elemental_composition import ElementalComposition
    from test_pt_paired_evidence import cfg
    r=load(tmp_path,[hyp("same",position=10,parent_base="G",modification_id="m22G")]);h=r.hypotheses[0]
    wrong=SimpleNamespace(canonical_structure_id="wrong",elemental_composition_delta=ElementalComposition.delta({"C":2,"H":4}))
    structure=SimpleNamespace(candidate_id="same",position_states={10:wrong},bond_states={},five_prime="inherited",three_prime="inherited")
    obs=SimpleNamespace(structures=(structure,),sheets={})
    a=build_modification_hypothesis_audit(r,config=cfg(),composite_observation=obs)
    assert any(x["Invalid_Reason"]=="hypothesis_id_structure_collision" for x in a.sheets["Mod_Hypothesis_Invalid"])

def test_oxidation_family_interpretations():
    f=_interpret_oxidation_family
    assert f(True,False,True,False,False,False,True)=="UNOXIDIZED_SUPPORTED"
    assert f(False,True,False,True,False,False,True)=="MONOOXIDE_SUPPORTED"
    assert f(True,True,True,True,False,False,True)=="BOTH_SUPPORTED_AS_MIXTURE"
    assert f(True,True,True,True,True,False,True)=="STRUCTURE_AMBIGUOUS"
    assert f(False,False,False,False,False,False,True)=="NEITHER_SUPPORTED"
    assert f(False,False,False,False,False,False,False)=="NOT_EVALUABLE"


def test_oxidation_family_rows_keep_candidates_separate(tmp_path):
    from types import SimpleNamespace
    from rna_masshunter.backbone_state import load_backbone_transformations
    from rna_masshunter.modification_constraints import load_transformations
    from rna_masshunter.sample_structure_schema import SampleStructureHypothesis,PositionHypothesis
    from rna_masshunter.structure_fragment import build_complete_structure_state
    from test_pt_paired_evidence import cfg
    r=load(tmp_path,[
        hyp("U",kind="composite_structure",position=37,parent_base="U",modification_ids=["s2U","cnm5U","side_chain_thioamide"]),
        hyp("O",kind="composite_structure",position=37,parent_base="U",modification_ids=["s2U","cnm5U","side_chain_thioamide_oxo1"])])
    transforms=load_transformations(ROOT/"data/modification_transforms_v2.yaml");bt=load_backbone_transformations(ROOT/"data/backbone_modifications.yaml")[0]
    structures=[]
    for candidate,ids in (("U",("s2U","cnm5U","side_chain_thioamide")),("O",("s2U","cnm5U","side_chain_thioamide_oxo1"))):
        sh=SampleStructureHypothesis(candidate,(PositionHypothesis(37,"U",ids),),())
        state,error=build_complete_structure_state(sh,SEQ,transforms,ROOT/"data/nucleoside_slots.yaml",bt);assert error is None;structures.append(state)
    support=[{"Candidate_ID":"U","Observable_Fragment_Count":1,"MS1_Matched_Fragment_Count":1,"MS1_Unique_Support_Count":1},
             {"Candidate_ID":"O","Observable_Fragment_Count":1,"MS1_Matched_Fragment_Count":1,"MS1_Unique_Support_Count":1}]
    obs=SimpleNamespace(structures=tuple(structures),sheets={"Composite_Support_Summary":support})
    a=build_modification_hypothesis_audit(r,config=cfg(),composite_observation=obs,project_root=ROOT)
    rows=a.sheets["Mod_Oxidation_Family"];states={x["Oxidation_State"] for x in rows}
    assert {"precursor","unoxidized","monooxide","dioxide"}<=states
    u=next(x for x in rows if x["Oxidation_State"]=="unoxidized");o=next(x for x in rows if x["Oxidation_State"]=="monooxide")
    assert u["Candidate_ID"]!=o["Candidate_ID"] and u["Composition_Delta"]!=o["Composition_Delta"]
    assert u["Final_Family_Interpretation"]=="BOTH_SUPPORTED_AS_MIXTURE"
    assert u["Delta_Oxidation_Da"]==pytest.approx(15.99491461957)
    assert u["Oxidation_Origin_Assessable"] is False and u["Possible_Ex_Vivo_Oxidation"]=="unknown"

def test_nucleoside_only_composite_does_not_select_pt_pair(tmp_path):
    from test_pt_paired_evidence import pair
    h=load(tmp_path,[hyp(kind="composite_structure",position=10,parent_base="G",modification_ids=["m22G","Gm"])]).hypotheses[0]
    assert not _matches_pair(h,pair(("m22G","Gm")))
