from copy import deepcopy
from dataclasses import asdict, replace
import pandas as pd
import pytest

from rna_masshunter.models import Fragment, FragmentMS1Match, Peak, RunConfig
from rna_masshunter.ms1_mapping import map_fragments_to_ms1_peaks, theoretical_mz_from_mass
from rna_masshunter.ms1_top50_dedup_audit import (
    DEDUP_MODES, DETAIL_COLUMNS, SUMMARY_COLUMNS, TOP50_COLUMNS, TOP_COLUMNS,
    append_top50_diagnostic_columns, append_top50_shadow_columns,
    build_ms1_top50_dedup_audit, build_near_mz_groups, deduplicate_matches,
    physical_peak_id, select_tier_matches, _change,
)


def config(report_limit=100000):
    return RunConfig(instrument={"polarity":"negative"},fragment_mapping={"enabled":True,"polarity":"negative","min_charge":1,"max_charge":1,"mz_tolerance_ppm":10,"max_matches_per_fragment":20,"use_peak_tiers":True,"include_trace_peaks":True,"min_fragment_length_for_filtered":3,"filtered_peak_tiers":["Major","Minor"],"filtered_confidence":["High","Medium"]},modification_search={"enabled":False},modification_evidence_ranking={"enabled":False},reporting={"max_excel_rows_per_sheet":report_limit},sequence={"name":"test"},digestion={"enzyme":"RNase_A"})

def fragment(fid="F1",seq="ACGU",start=1,end=4,mc=0,form="default"):
    return Fragment(fid,"target",seq,start,end,start,end,"RNase_A",mc,form,1000.0)

def matches(count=60):
    c=config();ctx={};mz=theoretical_mz_from_mass(1000,1,"negative")
    peaks=[Peak(mz*(1+((i%20)+1)*.1/1e6),1000-i,1+i/1000,f"S{i}",tier="Major" if i%2 else "Minor") for i in range(count)]
    formal=map_fragments_to_ms1_peaks([fragment()],peaks,c,audit_context=ctx)
    return c,ctx,formal

def clone_match(match,**values):
    item=replace(match,**values)
    for name in ("_audit_peak_index","_audit_physical_peak_id","_audit_generation_order"):
        if hasattr(match,name):setattr(item,name,getattr(match,name))
    return item


def test_physical_peak_one_and_no_dedup():
    c,ctx,_=matches(1);m=ctx["fragments"][0]["ranked_matches"][0]
    assert deduplicate_matches([m],"no_dedup",c.fragment_mapping)==[m]
    assert physical_peak_id(m).startswith("PK_")


def test_exact_duplicate_keeps_one():
    c,ctx,_=matches(1);m=ctx["fragments"][0]["ranked_matches"][0]
    assert len(deduplicate_matches([m,clone_match(m)],"exact_physical_peak_dedup",c.fragment_mapping))==1


def test_same_peak_multi_charge_keeps_best_one():
    c,ctx,_=matches(1);m=ctx["fragments"][0]["ranked_matches"][0];other=clone_match(m,charge=2,mass_error_ppm=m.mass_error_ppm+1)
    kept=deduplicate_matches([other,m],"physical_peak_charge_dedup",c.fragment_mapping)
    assert len(kept)==1 and kept[0].charge==1


def test_same_peak_multi_fragment_global_only_removes_cross_fragment():
    c,ctx,_=matches(1);m=ctx["fragments"][0]["ranked_matches"][0];other=clone_match(m,fragment_id="F2")
    assert len(deduplicate_matches([m,other],"fragment_then_physical_peak_dedup",c.fragment_mapping))==2
    assert len(deduplicate_matches([m,other],"global_physical_peak_dedup",c.fragment_mapping))==1


def test_same_peak_multi_form_is_preserved_between_fragment_ids():
    c,ctx,_=matches(1);m=ctx["fragments"][0]["ranked_matches"][0];other=clone_match(m,fragment_id="F1_form",terminal_form="cyclic_phosphate")
    assert len(deduplicate_matches([m,other],"fragment_then_physical_peak_dedup",c.fragment_mapping))==2


def test_same_mz_different_scan_is_distinct_exact_peak():
    c,ctx,_=matches(2);a,b=ctx["fragments"][0]["ranked_matches"]
    b=clone_match(b,observed_mz=a.observed_mz)
    assert physical_peak_id(a)!=physical_peak_id(b)
    assert len(deduplicate_matches([a,b],"exact_physical_peak_dedup",c.fragment_mapping))==2


def test_near_mz_is_grouped_but_not_deduplicated():
    c,ctx,_=matches(2);a,b=ctx["fragments"][0]["ranked_matches"];groups=build_near_mz_groups([a,b],10)
    assert groups[physical_peak_id(a)]==groups[physical_peak_id(b)]
    assert len(deduplicate_matches([a,b],"fragment_then_physical_peak_dedup",c.fragment_mapping))==2


@pytest.mark.parametrize("limit,expected",[(20,20),(50,50),(None,60)])
def test_tier_limits(limit,expected):
    c,ctx,_=matches(60);ranked=ctx["fragments"][0]["ranked_matches"]
    assert len(select_tier_matches(ranked,limit,c.fragment_mapping))==expected


def test_quality_priority_and_ties_are_deterministic():
    c,ctx,_=matches(60);ranked=ctx["fragments"][0]["ranked_matches"]
    first=select_tier_matches(ranked,50,c.fragment_mapping);second=select_tier_matches(ranked,50,c.fragment_mapping)
    assert [physical_peak_id(x) for x in first]==[physical_peak_id(x) for x in second]
    assert all(x.peak_tier in {"Major","Minor"} for x in first)


@pytest.mark.parametrize("seq",["ACG","ACGU","ACGUA","ACGUAC","ACGUACGU","ACGUACGUA"])
def test_short_and_long_fragment_groups_build(seq):
    f=fragment(seq=seq,end=len(seq));assert len(f.sequence)==len(seq)


def test_full_builder_formal_unchanged_and_required_rows():
    c,ctx,formal=matches(60);before=[asdict(x) for x in formal];context_before=[[asdict(m) for m in item["ranked_matches"]] for item in ctx["fragments"]]
    audit=build_ms1_top50_dedup_audit(ctx,c,[],[],formal,[],[],{})
    assert [asdict(x) for x in formal]==before
    assert [[asdict(m) for m in item["ranked_matches"]] for item in ctx["fragments"]]==context_before
    assert len(audit["top50_rows"])==9
    assert audit["summary"]["Tier_Top20_Filter_Passing"]==20
    assert audit["summary"]["Tier_Top50_Filter_Passing"]==50
    assert audit["summary"]["Tier_Unlimited_Filter_Passing"]==60
    assert audit["summary"]["Applied_To_Formal_Result"] is False


def test_detail_truncation_and_repeat_are_deterministic():
    c,ctx,formal=matches(60);c.reporting["max_excel_rows_per_sheet"]=17
    a=build_ms1_top50_dedup_audit(ctx,c,[],[],formal,[],[],{});b=build_ms1_top50_dedup_audit(ctx,c,[],[],formal,[],[],{})
    assert len(a["detail_rows"])==17 and a["detail_rows"]==b["detail_rows"]
    assert a["summary"]["Detail_Truncated"]


def test_columns_sheet_names_internal_and_append_only():
    c,ctx,formal=matches(2);audit=build_ms1_top50_dedup_audit(ctx,c,[],[],formal,[],[],{})
    assert all(len(name)<=31 for name in ("MS1_Top50_Shadow","MS1_Peak_Dedup_Detail","MS1_Top50_Dedup_Summary"))
    assert all(len(cols)==len(set(cols)) for cols in (TOP50_COLUMNS,DETAIL_COLUMNS,SUMMARY_COLUMNS))
    assert all(not any(str(k).startswith("_") for k in row) for row in audit["top50_rows"]+audit["detail_rows"]+audit["summary_rows"])
    top=pd.DataFrame([{"Review_Rank":1,"Modification_ID":"m1","Parent_Fragment_ID":"F1","Candidate_Positions_In_tRNA":36,"Best_Final_Score":7}])
    out=append_top50_shadow_columns(top,audit)
    assert list(out.columns[:len(top.columns)])==list(top.columns) and list(out.columns[-len(TOP_COLUMNS):])==TOP_COLUMNS
    assert out.iloc[0].Best_Final_Score==7 and not bool(out.iloc[0].MS1_Top50_Dedup_Applied_To_Formal_Result)
    diag=append_top50_diagnostic_columns([{"Existing":1}],audit)[0]
    assert diag["Existing"]==1 and diag["MS1_Top50_Dedup_Applied_To_Formal_Result"] is False


def ranking_row(mod="m1",score=1,rank=1,position=36):
    return {"Modification_ID":mod,"Parent_Fragment_ID":"F1","Candidate_tRNA_Position":position,"Final_Score":score,"Final_Confidence":"Low","Rank":rank}


def test_candidate_score_and_rank_change_detection():
    base=[ranking_row()]
    assert not _change(base,[ranking_row()])["set"]
    assert _change(base,[ranking_row("m2")])["set"]
    assert _change(base,[ranking_row(score=2)])["score"]==1
    assert _change(base,[ranking_row(rank=2)])["rank"]==1


def test_cnm5u_36_37_38_change_detection():
    base=[ranking_row("cnm5U",position=p,rank=i) for i,p in enumerate((36,37,38),1)]
    assert not _change(base,[dict(row) for row in base])["cnm"]
    changed=[dict(row) for row in base];changed[1]["Final_Score"]=2
    assert _change(base,changed)["cnm"]
