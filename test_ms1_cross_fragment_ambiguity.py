from copy import deepcopy
from dataclasses import asdict, replace
from types import SimpleNamespace

import pandas as pd
import pytest

from rna_masshunter.models import Fragment, FragmentMS1Match, RunConfig
from rna_masshunter.ms1_cross_fragment_ambiguity import (
    AMBIGUITY_COLUMNS, DETAIL_COLUMNS, DIAGNOSTIC_COLUMNS, SUMMARY_COLUMNS, TOP_COLUMNS,
    append_crossfrag_diagnostic_columns, append_crossfrag_top_columns,
    assignment_sort_key, assignment_weights, build_ms1_cross_fragment_ambiguity_audit,
    _candidate_change_rows,
)


def cfg(limit=100000):
    return RunConfig(
        instrument={"polarity":"negative"},
        fragment_mapping={"mz_tolerance_ppm":10,"MS1_Tier_Top50_Limit":50,"min_fragment_length_for_filtered":3,"filtered_peak_tiers":["Major","Minor"],"filtered_confidence":["High","Medium"]},
        modification_search={"enabled":False}, modification_evidence_ranking={"enabled":False},
        reporting={"max_excel_rows_per_sheet":limit}, sequence={"name":"synthetic"}, digestion={"enzyme":"RNase_A"},
    )


def frag(fid="F1", seq="ACG", start=1, end=3, mc=0, form="default", mass=1000.0):
    return Fragment(fid,"target",seq,start,end,start,end,"RNase_A",mc,form,mass)


def match(f, ppm=1.0, tier="Major", confidence="High", intensity=1000.0, pid="PK_SHARED", charge=1):
    m=FragmentMS1Match(fid := f"M_{f.fragment_id}",f.fragment_id,f.target_id,f.sequence,f.start,f.end,f.standard_start,f.standard_end,f.enzyme,f.missed_cleavages,f.terminal_form,f.unmodified_mass,charge,999.0,999.0,ppm/1e6*f.unmodified_mass,ppm,intensity,1.0,"scan=7",tier,confidence)
    setattr(m,"_audit_physical_peak_id",pid); setattr(m,"_audit_peak_index",7); setattr(m,"_audit_generation_order",1)
    return m


def context(*pairs):
    return {"configured_max_matches":20,"fragments":[{"fragment":f,"ranked_matches":list(ms)} for f,ms in pairs]}


def audit_for(matches, fragments=None, limit=100000, top_candidates=None):
    fragments=fragments or [frag(m.fragment_id,m.sequence,m.start,m.end,m.missed_cleavages,m.terminal_form,m.fragment_mass) for m in matches]
    pairs=[]
    for f in fragments: pairs.append((f,[m for m in matches if m.fragment_id==f.fragment_id]))
    top=None
    if top_candidates is not None:
        top={"candidates":{"tier_top50__no_dedup":top_candidates},"rankings":{"tier_top50__no_dedup":[]}}
    return build_ms1_cross_fragment_ambiguity_audit(context(*pairs),cfg(limit),[],[],[],[],{},top50_audit=top)


def test_single_assignment_has_no_ambiguity_group():
    f=frag(); a=audit_for([match(f)],[f]); assert a["ambiguity_rows"]==[] and a["summary"]["Shared_Physical_Peaks"]==0


def test_shared_same_sequence_different_position_and_position_risk():
    f1=frag("F1","ACG",1,3); f2=frag("F2","ACG",5,7)
    a=audit_for([match(f1),match(f2,ppm=2)],[f1,f2]); row=a["ambiguity_rows"][0]
    assert row["Same_Sequence_Multiple_Position"] and row["Unique_Position_Count"]==2
    assert row["Position_Discriminating"] and a["summary"]["Position_Localization_Risk"]=="High"


@pytest.mark.parametrize("f2,overlap,nonoverlap,multilen,multimc,multiform",[
    (frag("F2","CGU",2,4),True,False,False,False,False),
    (frag("F2","UGC",8,10),False,True,False,False,False),
    (frag("F2","ACGU",1,4,mass=1000),True,False,True,False,False),
    (frag("F2","CGU",2,4,mc=1),True,False,False,True,False),
    (frag("F2","CGU",2,4,form="cyclic"),True,False,False,False,True),
])
def test_structural_ambiguity_classes(f2,overlap,nonoverlap,multilen,multimc,multiform):
    f1=frag(); a=audit_for([match(f1),match(f2,ppm=2)],[f1,f2]); row=a["ambiguity_rows"][0]
    assert (row["Overlapping_Fragment_Count"]>0)==overlap
    assert (row["Nonoverlapping_Fragment_Count"]>0)==nonoverlap
    assert row["Multiple_Lengths"]==multilen and row["Multiple_Missed_Cleavages"]==multimc and row["Multiple_Terminal_Forms"]==multiform


def test_one_passing_one_failing_and_multiple_passing():
    f1=frag(); f2=frag("F2","CGU",2,4)
    one=audit_for([match(f1),match(f2,tier="Trace")],[f1,f2]); assert one["ambiguity_rows"][0]["Filter_Passing_Assignment_Count"]==1
    two=audit_for([match(f1),match(f2)],[f1,f2]); assert two["summary"]["Multiple_Passing_Assignment_Peaks"]==1


@pytest.mark.parametrize("change,expected",[
    ({"ppm":0.0},"equal_ppm_tie_break"),
    ({"ppm":8.0},"clear_ppm_winner"),
    ({"tier":"Minor"},"tier_winner"),
    ({"confidence":"Medium"},"confidence_winner"),
])
def test_best_assignment_dominance(change,expected):
    f1=frag(); f2=frag("F2","CGU",2,4); m1=match(f1,ppm=0); m2=match(f2,ppm=change.get("ppm",0),tier=change.get("tier","Major"),confidence=change.get("confidence","High"))
    a=audit_for([m1,m2],[f1,f2]); assert a["ambiguity_rows"][0]["Best_Assignment_Dominance"]==expected


def test_intensity_length_missed_cleavage_and_id_ties_are_deterministic():
    f1=frag("B","ACG",1,3,mc=1); f2=frag("A","ACGU",1,4,mass=1000)
    m1=match(f1,ppm=0,intensity=5); m2=match(f2,ppm=0,intensity=10)
    assert sorted([m1,m2],key=lambda m:assignment_sort_key(m,cfg().fragment_mapping))[0] is m2
    m1=replace(m1,intensity=10); setattr(m1,"_audit_physical_peak_id","PK_SHARED")
    assert sorted([m1,m2],key=lambda m:assignment_sort_key(m,cfg().fragment_mapping))[0] is m2
    f3=frag("A","ACG",1,3); f4=frag("B","ACG",1,3)
    assert assignment_sort_key(match(f3,ppm=0),cfg().fragment_mapping)<assignment_sort_key(match(f4,ppm=0),cfg().fragment_mapping)


def test_all_weighting_modes_and_sums():
    f1=frag(); f2=frag("F2","CGU",2,4); ms=[match(f1,ppm=1),match(f2,ppm=2)]
    w=assignment_weights(ms,cfg().fragment_mapping)
    assert sum(w["full_count"].values())==2 and sum(w["winner_take_all"].values())==1
    assert sum(w["equal_fraction"].values())==pytest.approx(1) and sum(w["quality_weighted_fraction"].values())==pytest.approx(1)
    assert sum(w["fragment_family_fraction"].values())==pytest.approx(1)
    assert w["ambiguity_flag_only"]==w["full_count"]


def test_fragment_family_keeps_nonoverlapping_families_independent():
    f1=frag("F1","ACG",1,3); f2=frag("F2","UGC",8,10)
    w=assignment_weights([match(f1),match(f2)],cfg().fragment_mapping)
    assert sum(w["fragment_family_fraction"].values())==2


def test_candidate_supporting_modification_eligible_and_cnm5u_group():
    f1=frag(); f2=frag("F2","CGU",2,4); ms=[match(f1),match(f2)]
    c={"source_type":"fragment","source_id":"F1","charge":1,"observed_mz":999.0,"modification_id":"cnm5U"}
    a=audit_for(ms,[f1,f2],top_candidates=[c]); row=a["ambiguity_rows"][0]
    assert row["Candidate_Supporting_Assignment_Count"]==1 and row["Modification_Eligible_Assignment_Count"]==1 and row["cnm5U_Relevant"]
    assert a["summary"]["Winner_Take_All_Candidate_Set_Changed"] is True
    assert a["summary"]["Equal_Fraction_Candidate_Set_Changed"] is False


def test_shadow_score_confidence_rank_and_top50_change_detection():
    base=[{"Modification_ID":"m","Parent_Fragment_ID":"F","Candidate_tRNA_Position":36,"Final_Score":1,"Final_Confidence":"Low","Rank":50}]
    same=_candidate_change_rows(base,[dict(base[0])])
    assert same["winner_take_all"]["score"]==0 and same["winner_take_all"]["rank"]==0 and not same["winner_take_all"]["top50"]
    changed=[dict(base[0],Final_Score=2,Final_Confidence="Medium",Rank=51)]
    result=_candidate_change_rows(base,changed)["winner_take_all"]
    assert result["score"]==1 and result["confidence"]==1 and result["rank"]==1 and result["top50"]


def test_formal_inputs_existing_columns_and_repeat_are_unchanged():
    f1=frag(); f2=frag("F2","CGU",2,4); ms=[match(f1),match(f2)]; ctx=context((f1,[ms[0]]),(f2,[ms[1]])); before=deepcopy([[asdict(m) for m in x["ranked_matches"]] for x in ctx["fragments"]])
    a=build_ms1_cross_fragment_ambiguity_audit(ctx,cfg(),[],[],[],[],{}); b=build_ms1_cross_fragment_ambiguity_audit(ctx,cfg(),[],[],[],[],{})
    assert [[asdict(m) for m in x["ranked_matches"]] for x in ctx["fragments"]]==before
    assert a["ambiguity_rows"]==b["ambiguity_rows"] and a["detail_rows"]==b["detail_rows"]
    top=pd.DataFrame([{"Review_Rank":1,"Modification_ID":"m","Parent_Fragment_ID":"F1","Candidate_Positions_In_tRNA":36,"Best_Final_Score":7}]); out=append_crossfrag_top_columns(top,a)
    assert list(out.columns[:len(top.columns)])==list(top.columns) and list(out.columns[-len(TOP_COLUMNS):])==TOP_COLUMNS and out.iloc[0].Best_Final_Score==7
    diag=append_crossfrag_diagnostic_columns([{"Existing":1}],a)[0]; assert diag["Existing"]==1 and diag["MS1_CrossFrag_Applied_To_Formal_Result"] is False


def test_detail_truncation_columns_and_sheet_names():
    fs=[frag(f"F{i}","ACG",i,i+2) for i in range(4)]; a=audit_for([match(f,ppm=i) for i,f in enumerate(fs)],fs,limit=2)
    assert len(a["detail_rows"])==2 and a["summary"]["Detail_Original_Row_Count"]==4 and a["summary"]["Detail_Truncated"]
    assert all(len(x)==len(set(x)) for x in (AMBIGUITY_COLUMNS,DETAIL_COLUMNS,SUMMARY_COLUMNS,TOP_COLUMNS,DIAGNOSTIC_COLUMNS))
    assert all(len(name)<=31 for name in ("MS1_CrossFrag_Ambiguity","MS1_CrossFrag_Detail","MS1_CrossFrag_Summary"))
    assert all(row["Applied_To_Formal_Result"] is False for row in a["ambiguity_rows"]+a["summary_rows"])
