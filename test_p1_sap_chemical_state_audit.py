from pathlib import Path
from types import SimpleNamespace
import pytest

from rna_masshunter.audit_policy import AuditPolicy, included_sheet_names
from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.models import Modification, Peak, RunConfig
from rna_masshunter.p1_sap_chemical_state_audit import (
    FALSE_FLAGS, LOCALIZATION, O_TO_S, _family_rows, _terminal_rows,
    build_p1_sap_chemical_state_audit, generate_chemical_state_candidates,
    match_and_group_features, compatible_ms2_provenance,
)

ROOT=Path(__file__).resolve().parent

def mod(mid="m22G",base="G",shift=28.03130012892):
    return Modification(mid,mid,shift,"biological",[base])

def cfg():
    return RunConfig(
        digestion={"enzyme":"Nuclease_P1"},
        alkaline_phosphatase={"enabled":True,"assume_complete":False},
        instrument={"polarity":"positive","ms1_tolerance_ppm":10},
        reconstruction={"mz_min":50,"mz_max":1000},
        p1_annotation={"mz_tolerance_ppm":10,"charge_states":[1]},
    )

def candidate(cid,family,mz=100.0,charge=1,composition="C1"):
    return {"Chemical_State_ID":cid,"Chemical_Family":family,"Product_Type":"monomer",
        "Elemental_Composition":composition,"Charge":charge,"Theoretical_mz":mz,
        "Search_Enabled":True,"Observable":True,**LOCALIZATION,**FALSE_FLAGS}

def feature(fid,cid,family,pid="PF_1",specific=True):
    return {"Feature_ID":fid,"Physical_Feature_ID":pid,"Chemical_State_ID":cid,
        "Chemical_Family":family,"Integrated_Intensity":1000,"Candidate_Specific":specific,
        "Apex_mz":100,"Mass_Error_ppm_at_Apex":0,"RT_Apex":1.0}

def test_candidate_dedup_ignores_source_positions_and_bonds():
    rows,_=generate_chemical_state_candidates("GGAG",[mod()],ROOT)
    m22=[r for r in rows if r["Nucleoside_Modification_State"]=="m22G" and r["Chemical_Family"]=="PHOSPHOROTHIOATE"]
    assert len(m22)==1
    assert m22[0]["Possible_Source_Position_Count"]==3
    assert m22[0]["Sequence_Position_Localized"] is False
    assert m22[0]["Original_Bond_Localized"] is False
    assert m22[0]["Localization_Status"]=="POSITION_NOT_RETAINED_BY_PREPARATION"

def test_charge_is_candidate_specific_and_composition_is_distinct():
    rows,_=generate_chemical_state_candidates("AG",[],ROOT,charges=(1,2))
    a=[r for r in rows if r["Base_or_Oligomer_Composition"]=="A" and r["Chemical_Family"]=="PHOSPHOROTHIOATE"]
    assert {r["Charge"] for r in a}=={1,2}
    g=next(r for r in rows if r["Base_or_Oligomer_Composition"]=="G" and r["Chemical_Family"]=="PHOSPHOROTHIOATE" and r["Charge"]==1)
    assert a[0]["Elemental_Composition"]!=g["Elemental_Composition"]

def test_normal_pt_pair_uses_elemental_mass_delta():
    rows,families=generate_chemical_state_candidates("A",[],ROOT)
    f=next(x for x in families if x["Chemical_Base_State"]=="A")
    assert f["Expected_Delta_Da"]==pytest.approx(O_TO_S.exact_mass)
    assert f["PT_Theoretical_mz"]-f["Normal_Theoretical_mz"]==pytest.approx(O_TO_S.exact_mass)

def test_undefined_thiophosphate_and_oxidized_pt_do_not_invent_mass():
    rows,_=generate_chemical_state_candidates("A",[],ROOT)
    undefined=[r for r in rows if r["Chemical_Family"] in {"THIOPHOSPHATE_LIKE","OXIDIZED_PT_DERIVATIVE"}]
    assert undefined and all(r["Model_Status"]=="MODEL_NOT_DEFINED" for r in undefined)
    assert all(r["Neutral_Mass"] is None and r["Theoretical_mz"] is None for r in undefined)

def test_28_profile_points_group_to_one_feature():
    c=candidate("PT","PHOSPHOROTHIOATE")
    peaks=[Peak(100+(i-14)*1e-6,100+i,1.0+i*0.002,f"scan{i}") for i in range(28)]
    features,_,counts=match_and_group_features([c],peaks,tolerance_ppm=20,max_rt_gap=.01)
    assert counts["PT"]==28 and len(features)==1
    row=features[0]
    assert row["Profile_Point_Count"]==28 and row["Spectrum_Count"]==28
    assert row["Apex_Intensity"]==127 and row["Integrated_Intensity"]==sum(range(100,128))
    assert row["Feature_Continuity_Status"]=="continuous_profile_feature"
    assert row["Feature_Eligible_For_Support"] is True

def test_disconnected_rt_and_charge_are_separate_features():
    candidates=[candidate("z1","PHOSPHOROTHIOATE",charge=1),candidate("z2","PHOSPHOROTHIOATE",charge=2)]
    peaks=[Peak(100,100,1.0,"s1"),Peak(100,200,1.02,"s2"),Peak(100,300,2.0,"s3")]
    features,_,_=match_and_group_features(candidates,peaks,max_rt_gap=.08)
    assert len([x for x in features if x["Chemical_State_ID"]=="z1"])==2
    assert {x["Charge"] for x in features}=={1,2}
    assert len({x["Physical_Feature_ID"] for x in features})==4

def test_outside_tolerance_is_not_a_match():
    features,_,counts=match_and_group_features([candidate("x","PHOSPHOROTHIOATE")],[Peak(100.1,100,1,"s")],tolerance_ppm=10)
    assert not features and counts["x"]==0

def test_competing_pt_and_non_pt_sulfur_is_ambiguous():
    candidates=[candidate("pt","PHOSPHOROTHIOATE"),candidate("s2u","SULFUR_CONTAINING_NON_PT_ALTERNATIVE")]
    features,competition,_=match_and_group_features(candidates,[Peak(100,500,1,"s1"),Peak(100,600,1.02,"s2")])
    pt=next(x for x in features if x["Chemical_State_ID"]=="pt")
    assert pt["Sulfur_Non_PT_Competition_Count"]==1
    assert pt["Candidate_Specific"] is False
    assert pt["Final_Interpretation"]=="PT_LIKE_STATE_AMBIGUOUS"
    assert competition

def test_isotope_status_is_provisional_without_envelope_fit():
    features,_,_=match_and_group_features([candidate("x","PHOSPHOROTHIOATE")],[Peak(100,500,1,"s")])
    assert features[0]["Envelope_Assessed"] is False
    assert features[0]["Isotope_Status"]=="provisional"

def test_family_interpretations_cover_terminal_states():
    spec={"Family_ID":"F","Chemical_Base_State":"A","Normal_State_ID":"normal_ref",
        "Dephosphorylated_State_ID":"d","Residual_State_ID":"n","PT_State_ID":"p",
        "Thiophosphate_State_ID":"t","Oxidized_PT_State_ID":"o","Normal_Composition":"x",
        "PT_Composition":"y","Normal_Theoretical_mz":1,"PT_Theoretical_mz":2,
        "Expected_Delta_Da":1,"Expected_Delta_mz":1}
    candidates=[candidate("d","DEPHOSPHORYLATED"),candidate("n","RESIDUAL_NORMAL_PHOSPHATE"),candidate("p","PHOSPHOROTHIOATE"),candidate("t","THIOPHOSPHATE_LIKE"),candidate("o","OXIDIZED_PT_DERIVATIVE")]
    assert _family_rows([spec],candidates,[feature("fd","d","DEPHOSPHORYLATED")])[0]["Family_Interpretation"]=="DEPHOSPHORYLATED_DOMINANT"
    mixed=_family_rows([spec],candidates,[feature("fd","d","DEPHOSPHORYLATED"),feature("fn","n","RESIDUAL_NORMAL_PHOSPHATE",pid="PF_2")])[0]
    assert mixed["Family_Interpretation"]=="MIXED_TERMINAL_STATES"
    pt=_family_rows([spec],candidates,[feature("fp","p","PHOSPHOROTHIOATE",specific=False)])[0]
    assert pt["Family_Interpretation"]=="PT_LIKE_STATE_AMBIGUOUS"

def test_sap_model_does_not_assert_complete_or_pt_removal():
    spec={"Family_ID":"F","Chemical_Base_State":"A","Normal_State_ID":"normal_ref",
        "Dephosphorylated_State_ID":"d","Residual_State_ID":"n","PT_State_ID":"p",
        "Thiophosphate_State_ID":"t","Oxidized_PT_State_ID":"o"}
    candidates=[candidate("d","DEPHOSPHORYLATED"),candidate("n","RESIDUAL_NORMAL_PHOSPHATE"),candidate("p","PHOSPHOROTHIOATE")]
    rows=_terminal_rows([spec],candidates,[])
    pt=next(x for x in rows if x["SAP_Substrate_State"]=="phosphorothioate monoester")
    assert pt["SAP_Removal_Unknown"] is True and pt["SAP_Removal_Confirmed"] is False
    normal=next(x for x in rows if x["SAP_Substrate_State"]=="normal phosphate monoester")
    assert normal["SAP_Removal_Expected"] is True and normal["SAP_Removal_Confirmed"] is False

def test_p1_audit_runs_pt_chemistry_and_never_localizes():
    result=build_p1_sap_chemical_state_audit(ROOT,"AG",[],cfg(),[mod()],audit_level="audit")
    assert result.metrics["P1_PT_Cleavage_Behavior"]=="unknown"
    assert any(x["Chemical_Family"]=="PHOSPHOROTHIOATE" for x in result.sheets["P1_SAP_Chemical_State"])
    assert all(x["Sequence_Position_Localized"] is False and x["Original_Bond_Localized"] is False for rows in result.sheets.values() for x in rows if "Sequence_Position_Localized" in x)
    assert result.summary_payload["sap_reaction_model"]["PT_removal_unknown"] is True
    assert result.summary_payload["sap_reaction_model"]["complete_reaction_asserted"] is False

def test_all_rows_keep_formal_false_flags():
    result=build_p1_sap_chemical_state_audit(ROOT,"AG",[Peak(100,100,1,"s")],cfg(),[],audit_level="full")
    for rows in result.sheets.values():
        for row in rows:
            assert row["Applied_To_Formal_Result"] is False
            assert row["Formal_Change_Ready"] is False
            assert row["Formal_Result_Changed"] is False

def test_audit_level_sheet_policy():
    names=["P1_SAP_Chemical_State","P1_SAP_PT_Family","P1_SAP_Terminal_Audit","P1_SAP_Features","P1_SAP_Competition","Cross_Enzyme_Chemistry","P1_SAP_MS2_Provenance"]
    standard,_=included_sheet_names(names,AuditPolicy.from_level("standard"))
    audit,_=included_sheet_names(names,AuditPolicy.from_level("audit"))
    full,_=included_sheet_names(names,AuditPolicy.from_level("full"))
    assert standard==[]
    assert set(audit)==set(names[:3])
    assert set(full)==set(names)


def test_multiple_profile_points_in_same_scan_and_rt_are_grouped():
    c=candidate("x","PHOSPHOROTHIOATE")
    features,_,_=match_and_group_features([c],[Peak(99.9999,100,1.0,"same"),Peak(100.0001,200,1.0,"same")],tolerance_ppm=2)
    assert len(features)==1
    assert features[0]["Profile_Point_Count"]==2
    assert features[0]["Spectrum_Count"]==1


def test_ms2_provenance_keeps_same_mz_competitors(monkeypatch):
    import rna_masshunter.p1_sap_chemical_state_audit as module
    spectrum={"ms level":2,"id":"ms2","scanList":{"scan":[{"scan start time":1.0}]},
        "precursorList":{"precursor":[{"selectedIonList":{"selectedIon":[{"selected ion m/z":100.0,"charge state":1}]},
        "isolationWindow":{"isolation window target m/z":100.0},"activation":{"collision energy":20}}]}}
    monkeypatch.setattr(module,"iter_spectra",lambda _:iter([spectrum]))
    rows=compatible_ms2_provenance("unused",[candidate("a","PHOSPHOROTHIOATE"),candidate("b","PHOSPHOROTHIOATE")])
    assert {row["Chemical_State_ID"] for row in rows}=={"a","b"}
    assert all(row["MS2_Model_Applicable"] is False for row in rows)


def test_single_point_and_long_background_trace_are_not_support():
    c=candidate("pt","PHOSPHOROTHIOATE")
    single,_,_=match_and_group_features([c],[Peak(100,100,1.0,"s")])
    assert single[0]["Feature_Eligible_For_Support"] is False
    assert single[0]["Chemical_State_Supported"] is False
    assert single[0]["Final_Interpretation"]=="NOT_EVALUABLE"
    trace=[Peak(100,100+i,1.0+i*.05,f"s{i}") for i in range(25)]
    long,_,_=match_and_group_features([c],trace,max_rt_gap=.08)
    assert len(long)==1 and long[0]["RT_Span"]>1.0
    assert long[0]["Feature_Continuity_Status"]=="continuous_background_trace"
    assert long[0]["Feature_Eligible_For_Support"] is False


def test_canonical_m22g_is_generated_without_dictionary_id():
    rows,_=generate_chemical_state_candidates("GG",[],ROOT)
    m22=[r for r in rows if r["Nucleoside_Modification_State"]=="m22G"]
    assert m22
    assert all(r["Possible_Source_Position_Count"]==2 for r in m22)
    assert all(r["Sequence_Position_Localized"] is False for r in m22)
