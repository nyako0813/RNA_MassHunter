from __future__ import annotations
from pathlib import Path
import yaml
import pytest
from rna_masshunter.cross_run_manifest import (
    ManifestValidationError, classify_run_independence, load_cross_run_manifest, strongest_independence,
)
from rna_masshunter.models import Peak, RunConfig
from rna_masshunter.pt_cross_run_recurrence import (
    aggregate_candidates, aggregate_pairs, candidate_key, neutral_candidate_key, recurrence_class,
)

BASE={"sample_id":"sample","biological_replicate_id":"bio1","sample_preparation_id":"prep1",
"digestion_id":"dig1","technical_replicate_id":"tech1","acquisition_batch_id":"batch1",
"instrument_method_id":"method1","condition":"WT","enzyme":"RNase_T1","sequence_id":"seq","organism":"org","notes":""}
def run(run_id="r1",**changes):
    value=dict(BASE,run_id=run_id,mzml_path=f"/{run_id}.mzML");value.update(changes);return value

def write_manifest(tmp_path,runs):
    for r in runs: (tmp_path/Path(r["mzml_path"]).name).write_text("x")
    local=[dict(r,mzml_path=Path(r["mzml_path"]).name) for r in runs]
    path=tmp_path/"runs.yaml";path.write_text(yaml.safe_dump({"schema_version":1,"runs":local}));return path

def detail(rid="r1",matched=True,ppm=1.0,specific=True,backbone="phosphorothioate",**extra):
    row={"Run_ID":rid,"Sequence_ID":"seq","Enzyme":"RNase_T1","Fragment_Start":8,"Fragment_End":12,
    "Terminal_Form":"default","Charge":2,"Nucleoside_State":"unmodified","Backbone_State":backbone,
    "Bond_ID":"10_11","Elemental_Composition":"C1H1","Theoretical_mz":1000.0,"Candidate_ID":"C",
    "Search_Mode":"discovery","Observable":True,"Matched":matched,"Mass_Error_ppm":ppm if matched else "",
    "Candidate_Specific":specific,"Isotope_Ambiguity":False,"Charge_Ambiguity":False,
    "Legacy_Competition_Count":0,"Composite_Competition_Count":0,"RT":1.0,"Intensity_Percentile":.2,
    "Continuity_Status":"Single_Scan"}
    row.update(extra);row["Cross_Run_Candidate_Key"]=candidate_key(row);row["Neutral_Candidate_Key"]=neutral_candidate_key(row);return row

def test_valid_manifest_and_relative_paths(tmp_path):
    loaded=load_cross_run_manifest(write_manifest(tmp_path,[run()]))
    assert loaded.schema_version==1 and Path(loaded.runs[0]["mzml_path"]).is_absolute()
@pytest.mark.parametrize("mutation,match",[
    (lambda rs:rs+[dict(rs[0],mzml_path="/other.mzML")],"duplicate_run_id"),
    (lambda rs:rs+[dict(rs[0],run_id="r2")],"duplicate_mzml_path"),
])
def test_manifest_duplicates(tmp_path,mutation,match):
    rows=mutation([run()]);path=tmp_path/"x.yaml";path.write_text(yaml.safe_dump({"schema_version":1,"runs":rows}))
    with pytest.raises(ManifestValidationError,match=match):load_cross_run_manifest(path,require_files=False)
def test_manifest_missing_file(tmp_path):
    path=tmp_path/"x.yaml";path.write_text(yaml.safe_dump({"schema_version":1,"runs":[run()]}))
    with pytest.raises(ManifestValidationError,match="missing_mzml_file"):load_cross_run_manifest(path)
def test_manifest_missing_metadata(tmp_path):
    bad=run();bad.pop("digestion_id");path=tmp_path/"x.yaml";path.write_text(yaml.safe_dump({"schema_version":1,"runs":[bad]}))
    with pytest.raises(ManifestValidationError,match="missing_metadata:digestion_id"):load_cross_run_manifest(path,require_files=False)
@pytest.mark.parametrize("changes,expected",[
    ({"technical_replicate_id":"tech2"},"TECHNICAL_REPLICATE"),
    ({"digestion_id":"dig2"},"INDEPENDENT_DIGESTION"),
    ({"sample_preparation_id":"prep2"},"INDEPENDENT_SAMPLE_PREPARATION"),
    ({"biological_replicate_id":"bio2"},"BIOLOGICAL_REPLICATE"),
    ({},"INDEPENDENT_INJECTION"),
])
def test_independence_levels(changes,expected):assert classify_run_independence(run(),run("r2",**changes))==expected
def test_unknown_independence():assert classify_run_independence(run(),run("r2",digestion_id=""))=="UNKNOWN_INDEPENDENCE"
def test_strongest_mixed_independence():assert strongest_independence([run(),run("r2",technical_replicate_id="t2"),run("r3",digestion_id="d2")])=="INDEPENDENT_DIGESTION"
def test_candidate_key_stable_across_run_fields():
    a=detail();b=dict(a,Run_ID="r9",Observed_mz=1000.1,Intensity=2)
    assert candidate_key(a)==candidate_key(b)
def test_charge_specific_and_neutral_key():
    a=detail();b=dict(a,Charge=3,Theoretical_mz=666.0)
    assert candidate_key(a)!=candidate_key(b) and neutral_candidate_key(a)==neutral_candidate_key(b)
def test_different_fragment_or_state_has_different_key():
    a=detail();assert candidate_key(a)!=candidate_key(dict(a,Fragment_End=13));assert candidate_key(a)!=candidate_key(dict(a,Nucleoside_State="m22G"))
@pytest.mark.parametrize("rows,normal,expected",[
    ([detail("r1"),detail("r2")],0,"RECURRENT_PT_MS1_SUPPORT"),
    ([detail("r1"),detail("r2",matched=False)],0,"SINGLE_RUN_PT_MS1_SUPPORT"),
    ([detail("r1",matched=False),detail("r2",matched=False)],0,"NO_RECURRENT_SUPPORT"),
    ([detail("r1"),detail("r2",ppm=20)],0,"RECURRENT_MASS_INCONSISTENT"),
    ([detail("r1"),detail("r2",specific=False)],0,"RECURRENT_BUT_AMBIGUOUS"),
    ([detail("r1",matched=False),detail("r2",matched=False)],2,"NORMAL_DOMINANT"),
])
def test_recurrence_classes(rows,normal,expected):
    metas={"r1":run(),"r2":run("r2",technical_replicate_id="t2")}
    for row in rows:row["_run_metadata"]=metas[row["Run_ID"]]
    assert recurrence_class(rows,normal_detected_runs=normal)[0]==expected

def test_independent_recurrence_class():
    rows=[detail("r1"),detail("r2")];metas={"r1":run(),"r2":run("r2",digestion_id="d2")}
    for row in rows:row["_run_metadata"]=metas[row["Run_ID"]]
    assert recurrence_class(rows)[0]=="RECURRENT_INDEPENDENT_PT_MS1_SUPPORT"
def test_not_evaluable_class():
    row=detail();row["_run_metadata"]=run();assert recurrence_class([row])[0]=="NOT_EVALUABLE"
def test_mixed_normal_and_pt_aggregation():
    rows=[detail("r1"),detail("r2"),detail("r1",backbone="normal_phosphate"),detail("r2",backbone="normal_phosphate")]
    summaries=aggregate_candidates(rows,[run(),run("r2",technical_replicate_id="t2")])
    pt=next(x for x in summaries if x["Backbone_State"]=="phosphorothioate")
    assert pt["Recurrence_Evidence_Class"]=="MIXED_NORMAL_AND_PT"
def test_pair_states_summary():
    rows=[]
    for rid,state in [("r1","PT_ONLY"),("r2","BOTH_PRESENT")]:rows.append({"Pair_Key":"p","Run_ID":rid,"Pair_State":state,"Normal_Candidate_Key":"n","PT_Candidate_Key":"p"})
    out=aggregate_pairs(rows,[run(),run("r2",digestion_id="d2")])[0]
    assert out["PT_Only_Run_Count"]==1 and out["Both_Present_Run_Count"]==1 and out["Pair_Recurrence_Class"]=="MIXED_NORMAL_AND_PT"
def test_formal_flags_are_false():
    out=aggregate_candidates([detail("r1"),detail("r2")],[run(),run("r2",technical_replicate_id="t2")])[0]
    assert out["Applied_To_Formal_Result"] is False and out["Formal_Change_Ready"] is False and out["Formal_Result_Changed"] is False

def test_standard_mode_does_not_read_manifest():
    from rna_masshunter.pt_cross_run_audit import build_pt_cross_run_audit
    result=build_pt_cross_run_audit("/definitely/missing.yaml",[],RunConfig(),audit_level="standard")
    assert result.sheets=={} and result.metrics["Manifest_Loaded"] is False

def test_cross_run_orchestrator_audit_and_full(tmp_path):
    from dataclasses import replace
    from rna_masshunter.masses import mz_from_neutral_mass
    from rna_masshunter.pt_cross_run_audit import build_pt_cross_run_audit
    from test_pt_paired_evidence import pair, cfg
    unmodified=pair(());unmodified=replace(unmodified,spec=replace(unmodified.spec,candidate_id="target|unmodified",search_mode="hypothesis_driven"))
    modified=pair(("m22G",));modified=replace(modified,spec=replace(modified.spec,candidate_id="target|modified",search_mode="hypothesis_driven"))
    runs=[run(),run("r2",digestion_id="dig2")];manifest=write_manifest(tmp_path,runs)
    mz=mz_from_neutral_mass(modified.modified_fragment.neutral_exact_mass,1,"negative")
    def loader(meta):
        offset=1e-6 if meta["run_id"]=="r1" else 2e-6
        return [Peak(mz*(1+offset),2000,rt=1.0,scan_id="scan=1"),Peak(mz*(1+offset),1500,rt=1.1,scan_id="scan=2")],{"MS1_Spectrum_Count":2,"MS2_Spectrum_Count":0}
    config=cfg(max_charge=1);config.peak_filtering={"trace_intensity_threshold":1000,"minor_intensity_threshold":5000,"major_intensity_threshold":25000}
    audit=build_pt_cross_run_audit(manifest,[unmodified,modified],config,audit_level="audit",peak_loader=loader)
    full=build_pt_cross_run_audit(manifest,[unmodified,modified],config,audit_level="full",peak_loader=loader)
    assert "PT_Cross_Run_Detail" not in audit.sheets and "PT_Cross_Run_Detail" in full.sheets
    assert len(full.sheets["PT_Cross_Run_Runs"])==2
    h4=next(x for x in full.sheets["PT_Cross_Run_Summary"] if x["Hypothesis_ID"]=="H4")
    assert h4["Detected_Run_Count"]==2 and h4["Recurrence_Evidence_Class"]=="RECURRENT_INDEPENDENT_PT_MS1_SUPPORT"
    assert any(r["Continuity_Status"]=="Multi_Scan_Continuous" for r in full.sheets["PT_Cross_Run_Detail"] if r["Hypothesis_ID"]=="H4")
    assert all(not r["Formal_Result_Changed"] for rows in full.sheets.values() for r in rows)

def test_cross_run_sheet_policy():
    from rna_masshunter.audit_policy import AuditPolicy,included_sheet_names
    names=["Run_summary","PT_Cross_Run_Runs","PT_Cross_Run_Summary","PT_Cross_Run_Pairs","PT_Cross_Run_Detail"]
    assert included_sheet_names(names,AuditPolicy.from_level("standard"))[0]==["Run_summary"]
    assert included_sheet_names(names,AuditPolicy.from_level("audit"))[0]==names[:-1]
    assert included_sheet_names(names,AuditPolicy.from_level("full"))[0]==names

def test_orchestrator_reports_missing_mzml_as_invalid_run(tmp_path):
    from rna_masshunter.pt_cross_run_audit import build_pt_cross_run_audit
    missing=run();missing["mzml_path"]=str(tmp_path/"missing.mzML")
    manifest=tmp_path/"manifest.yaml";manifest.write_text(yaml.safe_dump({"schema_version":1,"runs":[missing]}))
    result=build_pt_cross_run_audit(manifest,[],RunConfig(fragment_mapping={},reconstruction={},peak_filtering={}),audit_level="audit")
    row=result.sheets["PT_Cross_Run_Runs"][0]
    assert not row["Run_Valid"] and row["Invalid_Reason"].startswith("missing_mzml_file:")

def test_cross_run_ms2_reuses_precursor_and_composite_matcher(tmp_path):
    from dataclasses import replace
    from rna_masshunter.masses import mz_from_neutral_mass
    from rna_masshunter.models import MS2SpectrumInfo
    from rna_masshunter.pt_cross_run_audit import _build_ms2_ions,build_pt_cross_run_audit
    from test_pt_paired_evidence import pair,cfg
    target=pair(("m22G",));target=replace(target,spec=replace(target.spec,search_mode="hypothesis_driven"))
    config=cfg(max_charge=1);config.peak_filtering={};config.ms2_annotation.update({"mz_tolerance_ppm":20,"modified_fragment_min_ion_length":1})
    ion=next(x for x in _build_ms2_ions([target],config,"full") if x["Backbone_Informative"])
    precursor=mz_from_neutral_mass(ion["Parent_Neutral_Mass"],1,"negative")
    spectrum=MS2SpectrumInfo("ms2",1,1.0,precursor,1,None,1,ion["Theoretical_mz"],1000,1000,[(ion["Theoretical_mz"],1000)])
    manifest=write_manifest(tmp_path,[run()])
    def loader(meta):return [],{"MS1_Spectrum_Count":0,"MS2_Spectrum_Count":1,"_MS2_Spectra":[spectrum]}
    result=build_pt_cross_run_audit(manifest,[target],config,audit_level="full",peak_loader=loader)
    summary=next(x for x in result.sheets["PT_Cross_Run_Summary"] if x["Candidate_ID"]=="C")
    assert summary["MS2_Precursor_Compatible_Run_Count"]==1 and summary["MS2_Matched_Run_Count"]==1
    assert summary["Backbone_Localizing_Run_Count"]==1 and result.sheets["PT_Cross_Run_MS2_Detail"]
