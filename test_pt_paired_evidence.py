from dataclasses import replace
from pathlib import Path
import pytest
from rna_masshunter.audit_policy import AuditPolicy, included_sheet_names
from rna_masshunter.backbone_state import load_backbone_transformations
from rna_masshunter.cleavage_site_discovery import discover_candidate_cleavage_bonds
from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.models import Peak, RunConfig
from rna_masshunter.modification_composer import apply_transform_ids
from rna_masshunter.modification_constraints import load_transformations
from rna_masshunter.phosphorothioate_evidence import build_pt_evidence
from rna_masshunter.phosphorothioate_pairing import (
    PTPairSpec, PositionStateSpec, build_pt_pair, composition_difference,
    load_pt_pair_hypotheses,
)
from rna_masshunter.structure_fragment import extract_fragment_from_structure

ROOT = Path(__file__).resolve().parent
SEQ = "GCUCCGGUAGUGUAGUCCGGCCAAUCAUUCCGGCCUUUCGAGCCGAAGACUCGGGUUCGAAUCCCGGCCGGAGCACCA"
IDENTITY = {"sequence_id":"Mac_tRNA-Glu-UUC","organism":"Methanosarcina acetivorans","rule_set":"methanosarcina_acetivorans"}

def transforms(): return load_transformations(ROOT / "data/modification_transforms_v2.yaml")
def backbone(): return load_backbone_transformations(ROOT / "data/backbone_modifications.yaml")[0]
def cfg(mz_min=0, mz_max=5000, max_charge=3):
    return RunConfig(instrument={"polarity":"negative"}, reconstruction={"mz_min":mz_min,"mz_max":mz_max},
        fragment_mapping={"min_charge":1,"max_charge":max_charge,"mz_tolerance_ppm":10,
            "polarity":"negative","use_peak_tiers":False}, ms2_annotation={"default_charge":1})

def spec(transform_ids=(), bond_id="10_11", sequence=SEQ, sequence_id="Mac_tRNA-Glu-UUC", enzyme="RNase_T1"):
    candidate = next(x for x in discover_candidate_cleavage_bonds(sequence, sequence_id, enzyme) if x.bond_id == bond_id)
    positions = () if not transform_ids else (PositionStateSpec(candidate.left_position,candidate.left_base,tuple(transform_ids)),)
    return PTPairSpec("C", "H", "hypothesis_driven", sequence_id, sequence, enzyme, bond_id,
        positions, candidate.fragment_start, candidate.fragment_end, "default")

def pair(transform_ids=(), **kwargs):
    return build_pt_pair(spec(transform_ids, **kwargs), transforms(), ROOT/"data/nucleoside_slots.yaml", backbone())

def evidence_for(pair_obj, normal=False, modified=False, custom_cfg=None, legacy=()):
    peaks=[]
    if not normal and not modified:
        peaks.append(Peak(100.0,1000,rt=0.5,scan_id="X"))
    if normal:
        peaks.append(Peak(mz_from_neutral_mass(pair_obj.normal_fragment.neutral_exact_mass,1,"negative"),2000,rt=1.0,scan_id="N"))
    if modified:
        peaks.append(Peak(mz_from_neutral_mass(pair_obj.modified_fragment.neutral_exact_mass,1,"negative"),3000,rt=2.0,scan_id="P"))
    rows,states=build_pt_evidence([pair_obj],peaks,custom_cfg or cfg(max_charge=1),legacy_matches=legacy)
    return rows[0],states

def test_fixture_target_and_scientific_coordinates():
    loaded=load_pt_pair_hypotheses(ROOT/"data/sample_pt_pair_hypotheses.yaml",sequence=SEQ,
        transformations=transforms(),**IDENTITY)
    assert loaded.enabled and len(loaded.specs)==1 and not loaded.invalid_rows
    item=loaded.specs[0]
    assert SEQ[9:11]=="GU" and item.bond_id=="10_11" and item.position_states[0].transform_ids==("m22G",)
    assert "10_11" in {x.bond_id for x in discover_candidate_cleavage_bonds(SEQ,item.sequence_id,item.enzyme)}

def test_fixture_accepts_explicit_sequence_name_alias():
    loaded=load_pt_pair_hypotheses(ROOT/"data/sample_pt_pair_hypotheses.yaml",sequence=SEQ,
        sequence_id="MA_tRNA^Glu-UUC",organism=IDENTITY["organism"],rule_set=IDENTITY["rule_set"],transformations=transforms())
    assert len(loaded.specs)==1 and not loaded.invalid_rows

def test_fixture_target_guard_rejects_mismatch():
    loaded=load_pt_pair_hypotheses(ROOT/"data/sample_pt_pair_hypotheses.yaml",sequence=SEQ,
        sequence_id="wrong",organism=IDENTITY["organism"],rule_set=IDENTITY["rule_set"],transformations=transforms())
    assert not loaded.specs and loaded.invalid_rows[0]["Invalid_Reason"]=="target_identity_mismatch"

def test_m22g_transform_is_generic_exact_and_rejects_invalid_application():
    ts=transforms(); m22=next(x for x in ts if x.id=="m22G")
    assert m22.parent_bases==("G",) and m22.target_slot=="G_N2"
    assert m22.composition_delta==ElementalComposition.delta({"C":2,"H":4})
    state,result,_=apply_transform_ids("G",10,("m22G",),ts,ROOT/"data/nucleoside_slots.yaml")
    assert result.valid and state.elemental_composition_delta.canonical_string()=="C2H4"
    _,wrong_base,_=apply_transform_ids("C",10,("m22G",),ts,ROOT/"data/nucleoside_slots.yaml")
    assert not wrong_base.valid and wrong_base.reason_code=="wrong_parent_base"
    _,conflict,_=apply_transform_ids("G",10,("m2G","m22G"),ts,ROOT/"data/nucleoside_slots.yaml")
    assert not conflict.valid and conflict.reason_code in {"parent_child_double_count","superseded_component","from_state_mismatch"}

def test_normal_pt_pair_cancels_shared_m22g_and_exact_delta():
    p=pair(("m22G",)); expected=ElementalComposition.delta({"O":-1,"S":1})
    assert p.shared_modification_composition==ElementalComposition.delta({"C":2,"H":4})
    assert p.composition_delta==expected and p.expected_backbone_delta==expected
    assert abs((p.modified_fragment.neutral_exact_mass-p.normal_fragment.neutral_exact_mass)-expected.exact_mass)<1e-9
    combined=composition_difference(p.modified_fragment.elemental_composition,
        pair(()).normal_fragment.elemental_composition)
    assert combined==ElementalComposition.delta({"C":2,"H":4,"O":-1,"S":1})
    for z in (1,2,3):
        delta=mz_from_neutral_mass(p.modified_fragment.neutral_exact_mass,z,"negative")-mz_from_neutral_mass(p.normal_fragment.neutral_exact_mass,z,"negative")
        assert abs(delta-expected.exact_mass/z)<1e-10

def test_state_propagation_only_when_fragment_contains_position_or_bond():
    p=pair(("m22G",))
    outside=extract_fragment_from_structure(p.modified_structure,SEQ,11,15)
    assert not outside.included_modified_positions and not outside.included_backbone_modifications
    assert outside.elemental_composition==extract_fragment_from_structure(p.normal_structure,SEQ,11,15).elemental_composition

def test_discovery_is_position_independent_and_omits_terminal_bond():
    seq="AGUGG"; bonds=discover_candidate_cleavage_bonds(seq,"short","RNase_T1")
    assert [x.bond_id for x in bonds]==["2_3","4_5"]
    assert all(x.right_position<=len(seq) for x in bonds)
    other=pair(("m22G",),bond_id="6_7")
    assert other.spec.bond_id=="6_7" and other.composition_delta.canonical_string()=="O-1S1"

def test_enzyme_rule_discovery_and_non_site_rejection():
    assert discover_candidate_cleavage_bonds("ACUA","x","RNase_A")
    assert len(discover_candidate_cleavage_bonds("ACUA","x","Nuclease_P1"))==3
    with pytest.raises(ValueError,match="not_normal_cleavage_site"):
        build_pt_pair(PTPairSpec("x","","discovery","x",SEQ,"RNase_T1","9_10",(),1,12),
            transforms(),ROOT/"data/nucleoside_slots.yaml",backbone())

def test_pt_strong_normal_only_both_and_neither_classes():
    p=pair(("m22G",))
    strong=evidence_for(p,modified=True)[0]
    assert strong["Evidence_Class"]=="PT_CANDIDATE_SPECIFIC_MS1_SUPPORT"
    assert not strong["Mechanism_Discriminating"]
    assert "no observed normal/PT peak pair" in strong["Evidence_Reason"]
    assert strong["Normal_Cleavage_Mechanism"]=="stochastic_missed_cleavage"
    assert strong["Modified_Cleavage_Mechanism"]=="phosphorothioate_blocked"
    assert strong["Nucleoside_Blocking_Status"]=="unknown"
    assert evidence_for(p,normal=True)[0]["Evidence_Class"]=="NORMAL_ONLY_SUPPORT"
    assert evidence_for(p,normal=True,modified=True)[0]["Evidence_Class"]=="BOTH_PRESENT"
    neither,states=evidence_for(p)
    assert neither["Evidence_Class"]=="NEITHER_PRESENT"
    assert all(row["Nearest_Observed_mz"]==100.0 and row["Intensity"]==1000 for row in states)

def test_not_observable_and_mass_shift_inconsistent():
    p=pair(("m22G",))
    assert evidence_for(p,custom_cfg=cfg(mz_min=1,mz_max=2,max_charge=1))[0]["Evidence_Class"]=="NOT_OBSERVABLE"
    bad=replace(p,delta_consistency_error=0.1)
    assert evidence_for(bad,modified=True)[0]["Evidence_Class"]=="MASS_SHIFT_INCONSISTENT"

def test_legacy_competition_is_ambiguous_and_not_candidate_specific():
    p=pair(("m22G",)); mz=mz_from_neutral_mass(p.modified_fragment.neutral_exact_mass,1,"negative")
    legacy=[{"Observed_mz":mz,"Observed_Scan":"P","Observed_RT":2.0}]
    row,_=evidence_for(p,modified=True,legacy=legacy)
    assert row["Evidence_Class"]=="AMBIGUOUS_PEAK_ASSIGNMENT" and not row["Candidate_Specific"]

def test_isotope_competition_is_detected():
    p=pair(("m22G",)); shifted=replace(p.modified_fragment,
        neutral_exact_mass=p.modified_fragment.neutral_exact_mass-1.00335483507)
    competitor=replace(p,spec=replace(p.spec,candidate_id="isotope_competitor"),modified_fragment=shifted)
    mz=mz_from_neutral_mass(p.modified_fragment.neutral_exact_mass,1,"negative")
    rows,_=build_pt_evidence([p,competitor],[Peak(mz,3000,rt=2.0,scan_id="P")],cfg(max_charge=1))
    target=next(x for x in rows if x["Candidate_ID"]=="C")
    assert target["Evidence_Class"]=="AMBIGUOUS_PEAK_ASSIGNMENT" and target["Modified_Competition_Count"]>0

def test_physical_peak_is_not_duplicated_across_identical_search_modes():
    p=pair(("m22G",)); duplicate=replace(p,spec=replace(p.spec,candidate_id="same_chemistry",search_mode="discovery"))
    mz=mz_from_neutral_mass(p.modified_fragment.neutral_exact_mass,1,"negative")
    rows,_=build_pt_evidence([p,duplicate],[Peak(mz,3000,rt=2.0,scan_id="P")],cfg(max_charge=1))
    assert all(r["Modified_Competition_Count"]==0 for r in rows)

def test_audit_mode_does_not_retain_detail_rows():
    p=pair(("m22G",))
    mz=mz_from_neutral_mass(p.modified_fragment.neutral_exact_mass,1,"negative")
    rows,states=build_pt_evidence([p],[Peak(mz,3000,rt=2.0,scan_id="P")],
        cfg(max_charge=1),audit_level="audit",include_detail=False)
    assert states==[]
    assert set(rows[0])=={"Search_Mode","Candidate_ID","Evidence_Class","Observable",
        "Candidate_Specific","Modified_Physical_Peak_ID","Normal_Competition_Count",
        "Modified_Competition_Count"}


def test_pt_sheet_policy_standard_audit_full():
    names=["Run_summary","PT_Paired_Summary","PT_Discovery_Candidates","PT_Paired_Evidence","PT_State_Search"]
    standard,_=included_sheet_names(names,AuditPolicy.from_level("standard"))
    audit,_=included_sheet_names(names,AuditPolicy.from_level("audit"))
    full,_=included_sheet_names(names,AuditPolicy.from_level("full"))
    assert standard==["Run_summary"]
    assert audit==["Run_summary","PT_Paired_Summary","PT_Discovery_Candidates"]
    assert full==names


def test_pt_orchestrator_modes_and_formal_isolation():
    from copy import deepcopy
    from rna_masshunter.pt_paired_audit import build_pt_paired_audit
    formal=[{"observed_mz":999.0,"scan_id":"formal","rt":1.0}]; before=deepcopy(formal)
    run=RunConfig(sequence={"name":IDENTITY["sequence_id"]},
        organism={"species":IDENTITY["organism"],"rule_set":IDENTITY["rule_set"]},
        instrument={"polarity":"negative"},reconstruction={"mz_min":500,"mz_max":3000},
        fragment_mapping={"min_charge":1,"max_charge":2,"mz_tolerance_ppm":10,"polarity":"negative","use_peak_tiers":False},
        ms2_annotation={"default_charge":1,"precursor_match_tolerance_ppm":20,"constrain_by_precursor":True})
    audit=build_pt_paired_audit(ROOT,SEQ,IDENTITY["sequence_id"],[],[],run,legacy_matches=formal,audit_level="audit")
    full=build_pt_paired_audit(ROOT,SEQ,IDENTITY["sequence_id"],[],[],run,legacy_matches=formal,audit_level="full")
    assert "PT_Paired_Summary" in audit.sheets and "PT_Paired_Evidence" not in audit.sheets
    assert "PT_Paired_Evidence" in full.sheets and "PT_State_Search" in full.sheets
    assert {r["Search_Mode"] for r in full.sheets["PT_Paired_Summary"]}=={"hypothesis_driven","discovery"}
    assert all("PT_Candidate_Specific_MS1_Support_Count" in r for r in full.sheets["PT_Paired_Summary"])
    assert any(r["Bond_ID"]=="10_11" for r in full.sheets["PT_Discovery_Candidates"])
    assert formal==before
    assert all(row.get("Applied_To_Formal_Result") is False and row.get("Formal_Change_Ready") is False
        and row.get("Formal_Result_Changed") is False
        for rows in full.sheets.values() for row in rows)


def test_same_composition_structural_isomers_compete():
    a=pair(("Gm",)); b=pair(("m1G",))
    assert a.modified_fragment.elemental_composition==b.modified_fragment.elemental_composition
    mz=mz_from_neutral_mass(a.modified_fragment.neutral_exact_mass,1,"negative")
    rows,_=build_pt_evidence([a,b],[Peak(mz,3000,rt=2.0,scan_id="P")],cfg(max_charge=1))
    assert all(r["Evidence_Class"]=="AMBIGUOUS_PEAK_ASSIGNMENT" for r in rows)
    assert all(r["Modified_Competition_Count"]>0 for r in rows)
