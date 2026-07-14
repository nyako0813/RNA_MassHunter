from dataclasses import asdict
import pandas as pd
import pytest
from rna_masshunter.models import Fragment, Peak, RunConfig
from rna_masshunter.ms1_mapping import map_fragments_to_ms1_peaks, theoretical_mz_from_mass
from rna_masshunter.ms1_selection_strategy_audit import (
    STRATEGIES, STRATEGY_COLUMNS, DETAIL_COLUMNS, SUMMARY_COLUMNS, TOP_COLUMNS,
    append_selection_diagnostic_columns, append_top_selection_columns,
    build_ms1_selection_strategy_audit, passes_fragment_ms1_filter,
    select_strategy_matches, tier_then_error_sort_key,
)

def cfg(limit=20, report_limit=100000, max_charge=1):
    return RunConfig(instrument={"polarity":"negative"},fragment_mapping={"enabled":True,"polarity":"negative","min_charge":1,"max_charge":max_charge,"mz_tolerance_ppm":10,"max_matches_per_fragment":limit,"use_peak_tiers":True,"include_trace_peaks":True,"min_fragment_length_for_filtered":3,"filtered_peak_tiers":["Major","Minor"],"filtered_confidence":["High","Medium"]},modification_search={"enabled":False},modification_evidence_ranking={"enabled":False},reporting={"max_excel_rows_per_sheet":report_limit})
def fragment(sequence="ACGU"): return Fragment("F1","target",sequence,1,len(sequence),1,len(sequence),"RNase_T1",0,"default",1000.0)
def peaks(n,tiers=None,intensities=None):
    mz=theoretical_mz_from_mass(1000,1,"negative");out=[]
    for i in range(n):
        ppm=((i//2)+1)*(1 if i%2 else -1)*.15
        out.append(Peak(mz*(1+ppm/1e6),intensities[i] if intensities else 1000-i,1+i/1000,f"S{i}",tier=tiers[i] if tiers else "Major"))
    return out
def run(n,tiers=None,report_limit=100000,limit=20):
    c=cfg(limit,report_limit);ctx={};formal=map_fragments_to_ms1_peaks([fragment()],peaks(n,tiers),c,audit_context=ctx)
    audit=build_ms1_selection_strategy_audit(ctx,c,[],[],formal,[],[],{})
    return c,ctx,formal,audit

@pytest.mark.parametrize("n",[0,1,19,20,21,50])
def test_boundaries_and_all_strategies(n):
    c,ctx,formal,audit=run(n);s=audit["summary"]
    assert len(formal)==min(n,20)
    assert {r["Strategy"] for r in audit["strategy_rows"]}==set(STRATEGIES)
    assert s["Current_Selected_Count"]==min(n,20)
    assert s["Unlimited_Selected_Count"]==n
    assert all(r["Applied_To_Formal_Result"] is False for r in audit["strategy_rows"])

def test_filter_first_does_not_fill_with_failing_rows():
    _,_,_,a=run(25,["Trace"]*20+["Major"]*5)
    s=a["summary"]
    assert s["Current_Filter_Passing_Count"]==0
    assert s["Filter_First_Selected_Count"]==5==s["Filter_First_Filter_Passing_Count"]
    assert s["Unlimited_Selected_Count"]==5

def test_filter_passing_zero():
    _,_,_,a=run(25,["Trace"]*25)
    assert a["summary"]["Filter_First_Selected_Count"]==0
    assert a["summary"]["Unlimited_Selected_Count"]==0

def test_filter_passing_over_limit_and_recovery():
    _,_,_,a=run(25,["Trace"]*5+["Major"]*20)
    assert a["summary"]["Filter_First_Filter_Passing_Count"]==20
    assert a["summary"]["Unlimited_Filter_Passing_Count"]==20
    assert a["summary"]["Filter_First_Recovery_Fraction_Of_Unlimited"]==1

def test_tier_then_error_priorities_and_ties():
    c,ctx,_,_=run(21,["Trace"]*20+["Major"])
    ranked=ctx["fragments"][0]["ranked_matches"]
    chosen=select_strategy_matches(ranked,"tier_then_error",20,c.fragment_mapping)
    assert any(m.peak_tier=="Major" for m in chosen)
    keys=[tier_then_error_sort_key(m,c.fragment_mapping) for m in chosen]
    assert keys==sorted(keys)

def test_confidence_priority_precedes_error_and_intensity():
    c=cfg();ctx={};mz=theoretical_mz_from_mass(1000,1,"negative")
    ps=[Peak(mz*(1+9/1e6),100,1.0,f"x{i}",tier="Major") for i in range(20)]
    ps.append(Peak(mz*(1+.1/1e6),1,1.0,"low",tier="Trace"))
    map_fragments_to_ms1_peaks([fragment()],ps,c,audit_context=ctx)
    chosen=select_strategy_matches(ctx["fragments"][0]["ranked_matches"],"tier_then_error",20,c.fragment_mapping)
    assert all(m.peak_tier=="Major" for m in chosen)

def test_charge_and_physical_peak_tie_are_deterministic():
    c=cfg(max_charge=2);ctx={};mz1=theoretical_mz_from_mass(1000,1,"negative");mz2=theoretical_mz_from_mass(1000,2,"negative")
    ps=[Peak(mz1,100,1,"z1",tier="Major"),Peak(mz2,100,1,"z2",tier="Major")]
    map_fragments_to_ms1_peaks([fragment()],ps,c,audit_context=ctx)
    ranked=ctx["fragments"][0]["ranked_matches"]
    first=sorted(ranked,key=lambda m:tier_then_error_sort_key(m,c.fragment_mapping))
    second=sorted(ranked,key=lambda m:tier_then_error_sort_key(m,c.fragment_mapping))
    assert [getattr(x,"_audit_physical_peak_id") for x in first]==[getattr(x,"_audit_physical_peak_id") for x in second]

def test_status_classification_and_usage_columns():
    _,_,_,a=run(21,["Trace"]*20+["Major"])
    ff=[r for r in a["detail_rows"] if r["Strategy"]=="filter_first"]
    assert any(r["Selected_Vs_Current_Status"]=="added_by_strategy" for r in ff)
    assert any(r["Selected_Vs_Current_Status"]=="removed_vs_current" for r in ff)
    assert all("Used_By_Unlimited" in r for r in ff)

def test_current_reproduces_formal_and_capture_does_not_mutate():
    c=cfg();baseline=map_fragments_to_ms1_peaks([fragment()],peaks(25),c);ctx={};captured=map_fragments_to_ms1_peaks([fragment()],peaks(25),c,audit_context=ctx)
    before=[asdict(x) for x in captured];a=build_ms1_selection_strategy_audit(ctx,c,[],[],captured,[],[],{})
    assert [asdict(x) for x in captured]==before==[asdict(x) for x in baseline]
    assert [asdict(x) for x in a["strategy_matches"]["current"]]==before

def test_detail_truncation_and_rerun_are_deterministic():
    _,_,_,a=run(25,report_limit=17);_,_,_,b=run(25,report_limit=17)
    assert a["detail_rows"]==b["detail_rows"]
    assert a["summary"]["Detail_Original_Row_Count"]==100
    assert a["summary"]["Detail_Written_Row_Count"]==17
    assert a["summary"]["Detail_Truncated"]

def test_short_fragment_fails_exact_filter():
    c=cfg();ctx={};formal=map_fragments_to_ms1_peaks([fragment("AC")],peaks(1),c,audit_context=ctx)
    assert not passes_fragment_ms1_filter(formal[0],c.fragment_mapping)

def test_columns_sheets_and_append_only_behavior():
    _,_,_,a=run(21)
    assert all(len(name)<=31 for name in ("MS1_Selection_Strategy","MS1_Selection_Detail","MS1_Selection_Summary"))
    assert all(len(x)==len(set(x)) for x in (STRATEGY_COLUMNS,DETAIL_COLUMNS,SUMMARY_COLUMNS))
    top=pd.DataFrame([{"Review_Rank":1,"Modification_ID":"m1","Parent_Fragment_ID":"F1","Candidate_Positions_In_tRNA":36,"Best_Final_Score":7.0}])
    out=append_top_selection_columns(top,a)
    assert list(out.columns[:len(top.columns)])==list(top.columns)
    assert list(out.columns[-len(TOP_COLUMNS):])==TOP_COLUMNS
    assert out.iloc[0]["Best_Final_Score"]==7
    assert not bool(out.iloc[0]["MS1_Selection_Applied_To_Formal_Result"])
    diag=append_selection_diagnostic_columns([{"Existing":1}],a)[0]
    assert diag["Existing"]==1 and diag["MS1_Selection_Applied_To_Formal_Result"] is False

def test_internal_context_not_part_of_report_rows():
    _,_,_,a=run(1)
    assert all(not any(str(k).startswith("_") for k in row) for row in a["strategy_rows"]+a["detail_rows"]+a["summary_rows"])
    assert a["summary"]["Formal_Change_Ready"] is False
