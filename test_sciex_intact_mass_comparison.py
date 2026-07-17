from dataclasses import dataclass
from types import SimpleNamespace

import pandas as pd
import pytest

import main as main_module
from rna_masshunter.audit_policy import AUDIT_DETAIL, AUDIT_SUMMARY, AuditPolicy, included_sheet_names, sheet_category
from rna_masshunter.config import validate_config
from rna_masshunter.excel_report import (
    SCIEX_INTACT_OPTIONAL_RESULT_KEY, SCIEX_MASS_COMPARISON_DETAIL_SHEET,
    SCIEX_MASS_COMPARISON_OPTIONAL_RESULT_KEY, SCIEX_MASS_COMPARISON_SUMMARY_SHEET,
    _sciex_mass_comparison_excel_sheets, write_excel_report,
)
from rna_masshunter.models import IntactMassCandidate, RunConfig
from rna_masshunter.sciex_intact_mass_comparison import (
    DETAIL_COLUMNS, SUMMARY_COLUMNS, compare_sciex_intact_masses,
)

class Detection:
    def __init__(self, peaks, status="DETECTION_COMPLETED", input_status="SUPPORTED_INPUT"):
        self._peaks=[dict(x) for x in peaks]; self._diagnostics={"Detection_Status":status,"Input_Status":input_status}
    def peak_rows(self): return [dict(x) for x in self._peaks]
    def diagnostics_row(self): return dict(self._diagnostics)

def peak(mass, intensity=10, peak_id="P", **extra):
    return {"Peak_ID":peak_id,"Apex_Mass":mass,"Apex_Intensity_Raw":intensity,"Detection_Tier":"STRICT","Prominence":5.0,"Half_Prominence_Width_Da":2.0,**extra}

def compare(masses, theory=100.0, existing=(), **kwargs):
    return compare_sciex_intact_masses(Detection([peak(m, peak_id=f"P{i}") for i,m in enumerate(masses)]),theory,existing,**kwargs)

@pytest.mark.parametrize("mass,status",[(99.0,"STRICT_MATCH"),(101.0,"STRICT_MATCH"),(95.0,"BROAD_MATCH"),(105.0,"BROAD_MATCH"),(105.0001,"NO_MATCH")])
def test_status_boundaries_are_inclusive(mass,status):
    assert compare([mass]).details()[0]["Comparison_Status"]==status

@pytest.mark.parametrize("mass,sign",[(99.0,-1),(100.0,0),(101.0,1)])
def test_delta_and_ppm_sign_follow_observed_minus_theoretical(mass,sign):
    row=compare([mass]).details()[0]
    assert row["Delta_Mass"]==mass-100
    assert row["Absolute_Delta_Mass"]==abs(mass-100)
    assert (row["Delta_ppm"]>0)-(row["Delta_ppm"]<0)==sign

def test_nearest_existing_mass_and_signed_deltas_are_independent_of_theory_status():
    row=compare([104.0],existing=[102.0,110.0]).details()[0]
    assert row["Comparison_Status"]=="BROAD_MATCH"
    assert row["Nearest_Existing_Intact_Mass"]==102
    assert row["Nearest_Existing_Intact_Delta"]==2

def test_nearest_existing_tie_breaks_to_lower_mass():
    assert compare([105.0],existing=[110.0,100.0]).details()[0]["Nearest_Existing_Intact_Mass"]==100

def test_no_existing_result_retains_theory_comparison():
    row=compare([100.2]).details()[0]
    assert row["Comparison_Status"]=="STRICT_MATCH"
    assert row["Existing_Intact_Comparison_Status"]=="NO_EXISTING_INTACT_RESULT"
    assert row["Nearest_Existing_Intact_Mass"] is None

def test_dataclass_and_mapping_existing_results_are_supported():
    candidate=IntactMassCandidate(99.0,1,[1],1,1.0)
    result=compare_sciex_intact_masses(Detection([peak(101.0)]),100,[candidate,{"Reconstructed_Mass":101.25}])
    assert result.details()[0]["Nearest_Existing_Intact_Mass"]==101.25
    assert result.summaries()[0]["Existing_Intact_Mass_Count"]==2

def test_peak_fields_are_preserved_and_formal_flags_are_false():
    result=compare_sciex_intact_masses(Detection([peak(100.0,Possible_Shoulder=True,Broad_Peak_Flag=True,Peak_Area_Raw=12.0)]),100)
    row=result.details()[0]
    assert row["Possible_Shoulder"] is True and row["Broad_Peak_Flag"] is True and row["Peak_Area_Raw"]==12
    assert row["Shadow_Only"] is True
    assert not row["Applied_To_Formal_Score"] and not row["Applied_To_Ranking"] and not row["Applied_To_Candidate_Filtering"]
    assert not row["Molecular_Identity_Assigned"] and not row["Modification_Lookup_Performed"]

def test_output_is_mass_sorted_and_input_is_not_mutated():
    source=[peak(102,peak_id="late"),peak(98,peak_id="early")]
    snapshot=[dict(x) for x in source]
    rows=compare_sciex_intact_masses(Detection(source),100).details()
    assert [x["Peak_ID"] for x in rows]==["early","late"] and source==snapshot

def test_closest_peak_tie_uses_intensity_then_mass():
    detection=Detection([peak(99,2,"low"),peak(101,5,"high")])
    assert compare_sciex_intact_masses(detection,100).summaries()[0]["Closest_Peak_ID"]=="high"

def test_strongest_peak_tie_uses_lower_mass():
    detection=Detection([peak(101,5,"highmass"),peak(99,5,"lowmass")])
    assert compare_sciex_intact_masses(detection,100).summaries()[0]["Strongest_Peak_ID"]=="lowmass"

def test_noncompleted_detector_is_not_eligible_when_called_directly():
    result=compare_sciex_intact_masses(Detection([peak(100)],status="INVALID_AXIS"),100)
    assert result.details()[0]["Comparison_Status"]=="NOT_ELIGIBLE"
    assert result.summaries()[0]["Comparison_Status"]=="NOT_ELIGIBLE"

def test_no_theoretical_mass_returns_explicit_noneligible_rows_and_summary():
    result=compare([100],theory=None)
    assert result.details()[0]["Comparison_Status"]=="NO_THEORETICAL_MASS"
    assert result.details()[0]["Comparison_Eligible"] is False
    assert result.summaries()[0]["Comparison_Status"]=="NO_THEORETICAL_MASS"

@pytest.mark.parametrize("strict,broad,message",[(0,5,"strict"),(-1,5,"strict"),(1,.5,"broad"),(1,float("nan"),"broad")])
def test_invalid_tolerances_fail(strict,broad,message):
    with pytest.raises(ValueError,match=message): compare([100],strict_tolerance_da=strict,broad_tolerance_da=broad)

def config(enabled=True,comparison=True,strict=1,broad=5):
    return RunConfig(sciex_profile={"enabled":enabled,"path":"x","intact_peak_detection":{"enabled":True},"intact_mass_comparison":{"enabled":comparison,"strict_tolerance_da":strict,"broad_tolerance_da":broad}})

def wrapper(status="DETECTION_COMPLETED",peaks=None):
    return {SCIEX_INTACT_OPTIONAL_RESULT_KEY:{"result":Detection([peak(100)] if peaks is None else peaks,status),"source_file":"x"}}

@pytest.mark.parametrize("comparison,status,peaks",[(False,"DETECTION_COMPLETED",None),(True,"DETECTION_COMPLETED_WITH_WARNINGS",None),(True,"INVALID_AXIS",None),(True,"DETECTION_COMPLETED",[])])
def test_routing_skips_disabled_incomplete_or_zero_peak_inputs(comparison,status,peaks):
    assert main_module.build_sciex_intact_mass_comparison_optional_results(config(comparison=comparison),wrapper(status,peaks),100,[],[])=={}

def test_routing_builds_internal_result_and_accepts_missing_existing_results():
    result=main_module.build_sciex_intact_mass_comparison_optional_results(config(),wrapper(),100,[],[])
    assert SCIEX_MASS_COMPARISON_OPTIONAL_RESULT_KEY in result
    assert result[SCIEX_MASS_COMPARISON_OPTIONAL_RESULT_KEY].details()[0]["Existing_Intact_Comparison_Status"]=="NO_EXISTING_INTACT_RESULT"

def test_comparison_exception_isolated_and_detector_wrapper_retained(monkeypatch):
    source=wrapper(); warnings=[]
    monkeypatch.setattr(main_module,"compare_sciex_intact_masses",lambda *_a,**_k: (_ for _ in ()).throw(RuntimeError("boom")))
    assert main_module.build_sciex_intact_mass_comparison_optional_results(config(),source,100,[],warnings)=={}
    assert SCIEX_INTACT_OPTIONAL_RESULT_KEY in source
    assert warnings[-1]["Source"]=="sciex_intact_mass_comparison"

@pytest.mark.parametrize("strict,broad",[(0,5),(1,.5)])
def test_config_rejects_invalid_comparison_tolerances(strict,broad):
    with pytest.raises(ValueError): validate_config(config(strict=strict,broad=broad))

def test_disabled_sciex_config_does_not_validate_new_comparison_or_warn():
    cfg=RunConfig(sciex_profile={"enabled":False,"path":None,"intact_peak_detection":{"enabled":True},"intact_mass_comparison":"future"})
    warnings=[];validate_config(cfg,warnings)
    assert not any("sciex" in x["Message"] for x in warnings)

@pytest.mark.parametrize("level,present",[("standard",False),("audit",True),("full",True)])
def test_sheet_policy_is_standard_hidden_and_audit_full_visible(level,present):
    names=[SCIEX_MASS_COMPARISON_SUMMARY_SHEET,SCIEX_MASS_COMPARISON_DETAIL_SHEET]
    assert (included_sheet_names(names,AuditPolicy.from_level(level))[0]==names) is present


def writer_config():
    return SimpleNamespace(analysis={"mode":"full"},project={"name":"comparison-excel"},input={},organism={},sequence={},experiment={},instrument={},sciex_profile={},reconstruction={"enabled":False},digestion={"enabled":False},alkaline_phosphatase={},fragment_mapping={},modification_search={},peak_filtering={},p1_annotation={},ms2_annotation={},modification_evidence_ranking={},biological_context={},performance={},reporting={"max_excel_rows_per_sheet":1000,"truncate_large_sheets":True})

@pytest.mark.parametrize("level,expected",[("standard",False),("audit",True),("full",True)])
def test_real_excel_writer_routes_comparison_sheets_by_policy(tmp_path,level,expected):
    result=compare([100])
    report=write_excel_report(tmp_path/level,writer_config(),{},[],[],[],optional_results={SCIEX_MASS_COMPARISON_OPTIONAL_RESULT_KEY:result},audit_policy=AuditPolicy.from_level(level))
    from openpyxl import load_workbook
    workbook=load_workbook(report,read_only=True,data_only=True)
    try: names=workbook.sheetnames
    finally: workbook.close()
    assert (SCIEX_MASS_COMPARISON_DETAIL_SHEET in names) is expected
    assert (SCIEX_MASS_COMPARISON_SUMMARY_SHEET in names) is expected
    if expected:
        detail=pd.read_excel(report,sheet_name=SCIEX_MASS_COMPARISON_DETAIL_SHEET,header=2)
        assert len(detail)==1 and list(detail.columns)==DETAIL_COLUMNS

def test_sheet_categories_and_excel_columns_are_fixed():
    assert sheet_category(SCIEX_MASS_COMPARISON_SUMMARY_SHEET)==AUDIT_SUMMARY
    assert sheet_category(SCIEX_MASS_COMPARISON_DETAIL_SHEET)==AUDIT_DETAIL
    sheets=_sciex_mass_comparison_excel_sheets(compare([100]))
    assert list(sheets[SCIEX_MASS_COMPARISON_DETAIL_SHEET].columns)==DETAIL_COLUMNS
    assert list(sheets[SCIEX_MASS_COMPARISON_SUMMARY_SHEET].columns)==SUMMARY_COLUMNS
    assert len(SCIEX_MASS_COMPARISON_SUMMARY_SHEET)<=31
