from copy import deepcopy
from pathlib import Path
import pytest
from rna_masshunter.backbone_state import load_backbone_transformations
from rna_masshunter.blocked_cleavage_matcher import match_blocked_cleavage_fragments
from rna_masshunter.composite_ms1_matcher import match_composite_fragments_to_peaks
from rna_masshunter.composite_ms2_matcher import match_composite_ms2
from rna_masshunter.composite_ms2_propagation import generate_composite_theoretical_ions
from rna_masshunter.models import Fragment, Peak, RunConfig, MS2SpectrumInfo
from rna_masshunter.modification_constraints import load_transformations
from rna_masshunter.sample_structure_schema import load_sample_structure_hypotheses
from rna_masshunter.structure_fragment import build_complete_structure_state, extract_fragment_from_structure
from rna_masshunter.masses import load_base_masses, mz_from_neutral_mass

ROOT = Path(__file__).resolve().parent
SEQ = "GCUCCGGUAGUGUAGUCCGGCCAAUCAUUCCGGCCUUUCGAGCCGAAGACUCGGGUUCGAAUCCCGGCCGGAGCACCA"
TARGET_IDENTITY = {"name": "Mac_tRNA-Glu-UUC", "organism": "Methanosarcina acetivorans", "rule_set": "methanosarcina_acetivorans"}

def state_and_fragment():
    transforms = load_transformations(ROOT/"data/modification_transforms_v2.yaml")
    loaded = load_sample_structure_hypotheses(ROOT/"data/sample_structure_hypotheses.yaml", sequence=SEQ,
        transformations=transforms, backbone_bond_ids={f"{i}_{i+1}" for i in range(1,len(SEQ))},
        target_identity=TARGET_IDENTITY)
    bt = load_backbone_transformations(ROOT/"data/backbone_modifications.yaml")[0]
    state, error = build_complete_structure_state(loaded.hypotheses[0], SEQ, transforms,
        ROOT/"data/nucleoside_slots.yaml", bt)
    assert error is None
    parent = Fragment("F35_40","T",SEQ[34:40],35,40,35,40,"Nuclease_P1",0,"default",0.0)
    fragment = extract_fragment_from_structure(state,SEQ,35,40,fragment_id=parent.fragment_id,
        fragment_type=parent.enzyme)
    return state,parent,fragment,bt

def cfg(mz_min=0,mz_max=10000):
    return RunConfig(sequence={"name":"Mac_tRNA-Glu-UUC"},
        organism={"species":"Methanosarcina acetivorans","rule_set":"methanosarcina_acetivorans"},
        instrument={"polarity":"negative"},reconstruction={"mz_min":mz_min,"mz_max":mz_max},
        fragment_mapping={"min_charge":1,"max_charge":1,"mz_tolerance_ppm":10,"polarity":"negative"},
        ms2_annotation={"mz_tolerance_ppm":20,"default_charge":1,"modified_fragment_max_rows":10000})

def test_ms1_tolerance_match_and_no_observation():
    _,_,fragment,_=state_and_fragment(); theoretical=mz_from_neutral_mass(fragment.neutral_exact_mass,1,"negative")
    matched=match_composite_fragments_to_peaks([fragment],[Peak(theoretical,1000)],cfg())
    assert matched[0]["Match_Status"]=="matched"
    assert matched[0]["Support_Class"]=="unique_composite_support"
    isomeric=match_composite_fragments_to_peaks([fragment],[Peak(theoretical,1000)],cfg(),isomer_groups={fragment.candidate_id:"ISO"})
    assert isomeric[0]["Support_Class"]=="isomeric_unresolved"
    missed=match_composite_fragments_to_peaks([fragment],[Peak(theoretical+1,1000)],cfg())
    assert missed[0]["Match_Status"]=="no_observation"

def test_ms1_not_observable_is_distinct():
    _,_,fragment,_=state_and_fragment()
    rows=match_composite_fragments_to_peaks([fragment],[],cfg(mz_max=100))
    assert rows[0]["Match_Status"]=="not_observable"
    assert rows[0]["Not_Observable_Reason"]=="theoretical_mz_above_acquisition_range"

def test_ms2_state_propagates_only_to_containing_ions():
    state,parent,_,_=state_and_fragment()
    ions=generate_composite_theoretical_ions([state],[parent],SEQ,cfg())
    assert any(i["Position_Informative"] for i in ions)
    assert any(i["Backbone_Informative"] for i in ions)
    assert any(not i["Position_Informative"] for i in ions)

def test_ms2_match_is_shadow_only():
    state,parent,fragment,_=state_and_fragment(); ions=generate_composite_theoretical_ions([state],[parent],SEQ,cfg())
    ion=next(i for i in ions if i["Position_Informative"])
    spectrum=MS2SpectrumInfo("S1",1,1.0,mz_from_neutral_mass(fragment.neutral_exact_mass,1,"negative"),1,1000.0,1,ion["Theoretical_mz"],100.0,100.0,
        peaks=[(ion["Theoretical_mz"],100.0)])
    rows=match_composite_ms2(ions,[spectrum],cfg())
    assert rows and all(r["Applied_To_Formal_Result"] is False for r in rows)
    assert all(r["Formal_Change_Ready"] is False for r in rows)
    spectrum.precursor_mz += 1.0
    assert match_composite_ms2(ions,[spectrum],cfg()) == []

def test_phosphorothioate_blocked_is_not_stochastic_missed():
    state,parent,_,bt=state_and_fragment()
    rows=match_blocked_cleavage_fragments([state],SEQ,[],cfg(),bt,
        load_base_masses(ROOT/"data/base_masses.yaml"),[parent])
    assert rows
    assert {r["Cleavage_Status"] for r in rows}=={"phosphorothioate_blocked"}
    assert all(r["Blocked_Bond_ID"]=="37_38" for r in rows)

def test_invalid_hypothesis_never_reaches_matching(tmp_path):
    fixture=tmp_path/"bad.yaml"
    fixture.write_text("""schema_version: 1
enabled: true
hypotheses:
  - hypothesis_id: bad
    positions:
      37: {parent_base: U, transformations: [m5U, cm5U]}
""",encoding="utf-8")
    transforms=load_transformations(ROOT/"data/modification_transforms_v2.yaml")
    loaded=load_sample_structure_hypotheses(fixture,sequence=SEQ,transformations=transforms,
        backbone_bond_ids={f"{i}_{i+1}" for i in range(1,len(SEQ))})
    state,error=build_complete_structure_state(loaded.hypotheses[0],SEQ,transforms,
        ROOT/"data/nucleoside_slots.yaml",load_backbone_transformations(ROOT/"data/backbone_modifications.yaml")[0])
    assert state is None and error["Valid"] is False

def test_formal_inputs_are_not_mutated():
    state,parent,fragment,_=state_and_fragment()
    formal=[{"Modification_ID":"s2U","Final_Score":4.0,"Rank":1}]
    before=deepcopy(formal)
    match_composite_fragments_to_peaks([fragment],[],cfg(),legacy_matches=[])
    generate_composite_theoretical_ions([state],[parent],SEQ,cfg())
    assert formal==before
from rna_masshunter.audit_policy import AuditPolicy, included_sheet_names
from rna_masshunter.composite_observation_audit import build_composite_observation_audit

def test_phase2_orchestrator_audit_and_full_are_formal_immutable():
    state,parent,fragment,_=state_and_fragment()
    peak=Peak(mz_from_neutral_mass(fragment.neutral_exact_mass,1,"negative"),1000.0,rt=1.0,scan_id="MS1")
    formal_ms1=[{"sentinel":"formal-ms1"}]
    formal_ranking=[{"Modification_ID":"s2U","Final_Score":4.0,"Rank":1}]
    before=(deepcopy(formal_ms1),deepcopy(formal_ranking))
    phase1={"Composite_Mod_Candidates":[{"Position":37,"Legacy_Equivalent_IDs":"s2U",
        "Is_Isomeric":False}]}
    audit=build_composite_observation_audit(ROOT,SEQ,[parent],[peak],[],cfg(),
        load_base_masses(ROOT/"data/base_masses.yaml"),phase1,formal_ms1,formal_ranking,audit_level="audit")
    full=build_composite_observation_audit(ROOT,SEQ,[parent],[peak],[],cfg(),
        load_base_masses(ROOT/"data/base_masses.yaml"),phase1,formal_ms1,formal_ranking,audit_level="full")
    assert "Composite_Obs_Summary" in audit.sheets
    assert "Composite_MS1_Matches" not in audit.sheets
    assert "Composite_MS1_Matches" in full.sheets
    assert full.sheets["Composite_MS1_Matches"][0]["Match_Status"]=="matched"
    assert before==(formal_ms1,formal_ranking)
    assert all(row["Applied_To_Formal_Result"] is False for rows in full.sheets.values() for row in rows)
    assert all(row["Formal_Change_Ready"] is False for rows in full.sheets.values() for row in rows)

def test_phase2_sheet_policy_standard_audit_full():
    names=["Run_summary","Composite_Obs_Summary","Composite_MS1_Matches"]
    standard,_=included_sheet_names(names,AuditPolicy.from_level("standard"))
    audit,_=included_sheet_names(names,AuditPolicy.from_level("audit"))
    full,_=included_sheet_names(names,AuditPolicy.from_level("full"))
    assert standard==["Run_summary"]
    assert audit==["Run_summary","Composite_Obs_Summary"]
    assert full==names


def test_unmodified_fragment_observation_is_nondiscriminating():
    state, _, _, _ = state_and_fragment()
    outside = extract_fragment_from_structure(state, SEQ, 1, 6, fragment_id="outside", fragment_type="RNase_T1")
    theoretical = mz_from_neutral_mass(outside.neutral_exact_mass, 1, "negative")
    rows = match_composite_fragments_to_peaks([outside], [Peak(theoretical, 1000)], cfg())
    assert rows[0]["Support_Class"] == "observation_nondiscriminating"
    assert rows[0]["Legacy_Competition_Class"] == "OBSERVATION_NONDISCRIMINATING"


def test_legacy_comparison_uses_exact_complete_state_only():
    from rna_masshunter.legacy_composite_comparison import compare_legacy_composite
    canonical = "U@37|slot=a|slot2=b"
    support = [{"Candidate_ID":"C","Complete_Structure_ID":f"C|{canonical}|37_38=phosphorothioate",
        "Modified_Positions":"37","MS1_Unique_Support_Count":1,"MS1_Nondiscriminating_Count":2}]
    phase = [
        {"Candidate_ID":"wrong","Position":37,"Complete_Structure_ID":"U@37|slot=x","Legacy_Equivalent_IDs":"cnm5U"},
        {"Candidate_ID":"exact","Position":37,"Complete_Structure_ID":canonical,"Included_Component_IDs":"cnm5U;s2U"},
    ]
    rows = compare_legacy_composite(support, phase, [], audit_level="full")
    assert rows[0]["Exact_Phase1_Candidate_IDs"] == "exact"
    assert rows[0]["Legacy_Candidate_IDs"] == ""
    assert rows[0]["Comparison_Class"] == "LEGACY_PARENT_OF_COMPOSITE"
