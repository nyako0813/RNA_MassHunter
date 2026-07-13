"""Non-mutating shadow audit for Fragment MS1 match truncation."""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, is_dataclass
from statistics import median
from typing import Any

import pandas as pd

from rna_masshunter.evidence_ranking import build_modification_evidence_ranking
from rna_masshunter.modification_search import search_known_modifications

AUDIT_COLUMNS = [
"Fragment_ID","Fragment_Start","Fragment_End","Fragment_Sequence","Fragment_Length","Missed_Cleavages","Terminal_Form","Theoretical_Neutral_Mass",
"Pre_Truncation_Match_Count","Configured_Max_Matches","Retained_Match_Count","Discarded_Match_Count","Truncation_Applied","Retained_Fraction",
"Retained_Unique_Physical_Peak_Count","Retained_Charge_State_Count","Retained_Min_Abs_Error_PPM","Retained_Median_Abs_Error_PPM","Retained_Max_Intensity","Retained_Median_Intensity",
"Retained_Major_Count","Retained_Minor_Count","Retained_Trace_Count","Retained_High_Confidence_Count","Retained_Medium_Confidence_Count","Retained_Low_Confidence_Count","Retained_Filter_Passing_Count","Retained_Modification_Eligible_Count",
"Discarded_Unique_Physical_Peak_Count","Discarded_Charge_State_Count","Discarded_Min_Abs_Error_PPM","Discarded_Median_Abs_Error_PPM","Discarded_Max_Intensity","Discarded_Median_Intensity",
"Discarded_Major_Count","Discarded_Minor_Count","Discarded_Trace_Count","Discarded_High_Confidence_Count","Discarded_Medium_Confidence_Count","Discarded_Low_Confidence_Count","Discarded_Filter_Passing_Count","Discarded_Modification_Eligible_Count",
"Better_Error_In_Discarded","Higher_Intensity_In_Discarded","Better_Peak_Tier_In_Discarded","Better_Confidence_In_Discarded","Filter_Passing_Discarded_While_Retained_Fails",
"New_Charge_State_Only_In_Discarded","New_Physical_Peak_Only_In_Discarded","Modification_Candidate_Only_In_Discarded","Candidate_Family_Only_In_Discarded",
"Truncation_Risk_Level","Truncation_Risk_Reason","Recommended_Action","Applied_To_Final_Score"]
DETAIL_COLUMNS = [
"Fragment_ID","Pre_Truncation_Rank","Retained_Or_Discarded","Retention_Cutoff","Retention_Sort_Key","Physical_Peak_ID","Peak_Index","Observed_MZ","Intensity","Peak_Tier","Confidence","Charge",
"Observed_Neutral_Mass","Theoretical_Neutral_Mass","Mass_Error_Da","Mass_Error_PPM","Abs_Mass_Error_PPM","RT","Passes_Fragment_MS1_Filter","Eligible_For_Modification_Search",
"Modification_Mass_Shift","Possible_Modification_Count","Possible_Modification_IDs","Candidate_Key_List","Candidate_Only_From_Discarded",
"Better_Than_Worst_Retained_By_Error","Better_Than_Worst_Retained_By_Intensity","Better_Than_Worst_Retained_By_Tier","Better_Than_Worst_Retained_By_Confidence"]
SUMMARY_COLUMNS = [
"Total_Fragment_Count","Fragments_With_Any_Match","Truncated_Fragment_Count","Total_Pre_Truncation_Matches","Total_Retained_Matches","Total_Discarded_Matches",
"Total_Retained_Unique_Physical_Peaks","Total_Discarded_Unique_Physical_Peaks","Fragments_With_Better_Error_Discarded","Fragments_With_Higher_Intensity_Discarded",
"Fragments_With_Better_Tier_Discarded","Fragments_With_Better_Confidence_Discarded","Fragments_With_Filter_Passing_Discarded","Fragments_With_Modification_Eligible_Discarded",
"Discarded_Only_Modification_Candidate_Count","Discarded_Only_Candidate_Family_Count","Baseline_Max_Matches_Per_Fragment","Shadow_Max_Matches_Per_Fragment","Shadow_Truncation_Disabled",
"Baseline_Fragment_MS1_Match_Count","Shadow_Fragment_MS1_Match_Count","Baseline_Filtered_Match_Count","Shadow_Filtered_Match_Count",
"Baseline_Modification_Candidate_Count","Shadow_Modification_Candidate_Count","Added_Shadow_Modification_Candidate_Count","Removed_Shadow_Modification_Candidate_Count",
"Changed_Shadow_Modification_Candidate_Count","Baseline_Ranked_Candidate_Count","Shadow_Ranked_Candidate_Count","Candidate_Key_Set_Changed","Final_Score_Would_Change",
"Final_Confidence_Would_Change","Rank_Would_Change","Top50_Membership_Would_Change","cnm5U_Result_Would_Change",
"Current_Sort_Filter_Passing_Retained","Filter_First_Filter_Passing_Retained","Unique_Physical_Peak_Filter_Passing_Retained",
"Balanced_Per_Charge_Filter_Passing_Retained","Tier_Then_Error_Filter_Passing_Retained","Current_Truncation_Sort_Key","Current_Sort_Tie_Break",
"Detail_Original_Row_Count","Detail_Written_Row_Count","Detail_Truncated","Detail_Truncation_Reason","Shadow_Score_Recalculation_Status","Shadow_Score_Recalculation_Reason",
"Shadow_Additional_Time_Seconds","Shadow_Peak_Tracked_Memory_MiB",
"Overall_Truncation_Risk","Overall_Truncation_Conclusion","Recommended_Max_Matches_Per_Fragment","Recommended_Next_Action","Audit_Mode","Applied_To_Final_Score"]
TOP_COLUMNS = ["MS1_Truncation_Affected","MS1_Truncated_Fragment_Count","MS1_Discarded_Support_Match_Count","MS1_Discarded_Filter_Passing_Count",
"MS1_Discarded_Modification_Support_Count","Shadow_Expanded_MS1_Support_Count","Shadow_Expanded_Final_Score","Shadow_Expanded_Final_Confidence",
"Shadow_Expanded_Rank","Shadow_Rank_Changed","MS1_Truncation_Risk","MS1_Truncation_Recommendation","MS1_Truncation_Applied_To_Final_Score"]
DIAGNOSTIC_COLUMNS = ["MS1_Truncation_Audit_Available","MS1_Truncation_Total_Fragments","MS1_Truncated_Fragment_Count","MS1_Pre_Truncation_Match_Count",
"MS1_Retained_Match_Count","MS1_Discarded_Match_Count","MS1_Discarded_Filter_Passing_Count","MS1_Discarded_Modification_Eligible_Count",
"MS1_Discarded_Only_Candidate_Count","MS1_Shadow_Candidate_Key_Set_Changed","MS1_Shadow_Final_Score_Changed","MS1_Shadow_Rank_Changed",
"MS1_Truncation_Overall_Risk","MS1_Truncation_Recommendation","MS1_Truncation_Applied_To_Final_Score"]

def _raw(x):
    return asdict(x) if is_dataclass(x) else dict(x)

def _tier(x):
    return {"major":3,"minor":2,"trace":1}.get(str(x or "").lower(),0)

def _conf(x):
    return {"high":3,"medium":2,"low":1}.get(str(x or "").lower(),0)

def _pid(m):
    return str(getattr(m,"_audit_physical_peak_id","") or f"PK_{getattr(m,'_audit_peak_index',0):06d}_{m.observed_mz:.8f}")

def _mkey(m):
    return (str(m.fragment_id),int(m.charge),round(float(m.observed_mz),10))

def _ckey(c):
    r=_raw(c)
    return (str(r.get("source_type") or ""),str(r.get("source_id") or ""),int(r.get("charge") or 0),round(float(r.get("observed_mz") or 0),10),str(r.get("modification_id") or ""))

def _rkey(r):
    p=r.get("Candidate_tRNA_Position",r.get("Candidate_Positions_In_tRNA",""))
    if isinstance(p,float) and p.is_integer(): p=int(p)
    return (str(r.get("Modification_ID") or ""),str(r.get("Parent_Fragment_ID") or ""),str(p or ""))

def _passes(m,cfg):
    tiers={str(x).lower() for x in cfg.get("filtered_peak_tiers",["Major","Minor"]) or []}
    confs={str(x).lower() for x in cfg.get("filtered_confidence",["High","Medium"]) or []}
    return len(m.sequence or "")>=int(cfg.get("min_fragment_length_for_filtered",3) or 3) and (not tiers or str(m.peak_tier or "").lower() in tiers) and (not confs or str(m.confidence or "").lower() in confs)

def _stats(ms,eligible,prefix,cfg):
    errs=[abs(float(m.mass_error_ppm)) for m in ms]; ints=[float(m.intensity or 0) for m in ms]
    t=Counter(str(m.peak_tier or "").lower() for m in ms); c=Counter(str(m.confidence or "").lower() for m in ms)
    return {f"{prefix}_Unique_Physical_Peak_Count":len({_pid(m) for m in ms}),f"{prefix}_Charge_State_Count":len({m.charge for m in ms}),
    f"{prefix}_Min_Abs_Error_PPM":min(errs) if errs else "",f"{prefix}_Median_Abs_Error_PPM":median(errs) if errs else "",
    f"{prefix}_Max_Intensity":max(ints) if ints else "",f"{prefix}_Median_Intensity":median(ints) if ints else "",
    f"{prefix}_Major_Count":t["major"],f"{prefix}_Minor_Count":t["minor"],f"{prefix}_Trace_Count":t["trace"],
    f"{prefix}_High_Confidence_Count":c["high"],f"{prefix}_Medium_Confidence_Count":c["medium"],f"{prefix}_Low_Confidence_Count":c["low"],
    f"{prefix}_Filter_Passing_Count":sum(_passes(m,cfg) for m in ms),f"{prefix}_Modification_Eligible_Count":sum(_mkey(m) in eligible for m in ms)}

def _family(c,lookup):
    r=_raw(c); mid=str(r.get("modification_id") or ""); m=lookup.get(mid); raw=getattr(m,"raw",{}) or {}
    return str(raw.get("chemical_group") or raw.get("near_isobaric_group") or getattr(m,"category","") or mid)

def _recommend(risk,filter_rows,new_charge,new_peak):
    if risk=="none": return "no_action_needed"
    if risk=="high": return "filter_before_truncation" if filter_rows else "increase_max_matches_per_fragment"
    if new_peak:return "deduplicate_physical_peaks_before_truncation"
    if new_charge:return "allocate_limit_per_charge_state"
    if filter_rows:return "filter_before_truncation"
    return "retain_current_limit"

def _strategy_counts(ms,n,cfg):
    current=ms[:n]
    ff=sorted(ms,key=lambda m:(not _passes(m,cfg),abs(m.mass_error_ppm),-m.intensity))[:n]
    unique=[];seen=set()
    for m in ms:
        if _pid(m) not in seen:seen.add(_pid(m));unique.append(m)
    pools=defaultdict(list)
    for m in ms:pools[m.charge].append(m)
    balanced=[]
    while len(balanced)<n and any(pools.values()):
        for z in sorted(pools):
            if pools[z] and len(balanced)<n:balanced.append(pools[z].pop(0))
    tier=sorted(ms,key=lambda m:(-_tier(m.peak_tier),abs(m.mass_error_ppm),-m.intensity))[:n]
    return [sum(_passes(m,cfg) for m in x) for x in (current,ff,unique[:n],balanced,tier)]

def build_ms1_truncation_audit(context,config,modifications,intact_results,baseline_matches,baseline_candidates,baseline_ranking,ms2_results,rule_set=None,pathways=None):
    cfg=config.fragment_mapping or {}; cutoff=int(context.get("configured_max_matches") or cfg.get("max_matches_per_fragment") or 20)
    fragments=list(context.get("fragments") or []); all_matches=[m for x in fragments for m in x.get("ranked_matches",[])]
    shadow_candidates=search_known_modifications(all_matches,intact_results,modifications,config,warnings=None)
    shadow_ranking,_=build_modification_evidence_ranking(config,modifications,[x["fragment"] for x in fragments],all_matches,shadow_candidates,ms2_results,rule_set=rule_set,pathways=pathways)
    base_ck={_ckey(c) for c in baseline_candidates}; shadow_ck={_ckey(c) for c in shadow_candidates}
    cidx=defaultdict(list)
    for c in shadow_candidates:
        r=_raw(c)
        if str(r.get("source_type") or "").lower()=="fragment":cidx[(str(r.get("source_id") or ""),int(r.get("charge") or 0),round(float(r.get("observed_mz") or 0),10))].append(c)
    eligible=set(cidx); lookup={str(getattr(m,"id","")):m for m in modifications}; audit=[];detail=[];internal={}; strategies=[0]*5
    for item in fragments:
        f=item["fragment"]; ranked=list(item.get("ranked_matches") or []); kept=ranked[:cutoff];drop=ranked[cutoff:]
        kc=[c for m in kept for c in cidx.get(_mkey(m),[])]; dc=[c for m in drop for c in cidx.get(_mkey(m),[])]
        kf=[m for m in kept if _passes(m,cfg)];df=[m for m in drop if _passes(m,cfg)]
        be=bool(drop and kept and min(abs(m.mass_error_ppm) for m in drop)<max(abs(m.mass_error_ppm) for m in kept))
        bi=bool(drop and kept and max(m.intensity for m in drop)>max(m.intensity for m in kept))
        bt=bool(drop and kept and max(_tier(m.peak_tier) for m in drop)>min(_tier(m.peak_tier) for m in kept))
        bc=bool(drop and kept and max(_conf(m.confidence) for m in drop)>min(_conf(m.confidence) for m in kept))
        nc=bool({m.charge for m in drop}-{m.charge for m in kept});np=bool({_pid(m) for m in drop}-{_pid(m) for m in kept})
        onlyc={_ckey(c) for c in dc}-base_ck;onlyf={_family(c,lookup) for c in dc}-{_family(c,lookup) for c in kc}
        highq=bool(df and kept and (min(abs(m.mass_error_ppm) for m in df)<max(abs(m.mass_error_ppm) for m in kept) or max(_tier(m.peak_tier) for m in df)>min(_tier(m.peak_tier) for m in kept) or max(_conf(m.confidence) for m in df)>min(_conf(m.confidence) for m in kept)))
        risk="none" if not drop else "high" if onlyc or highq else "moderate" if df or dc or nc or np else "low"
        reason={"none":"no truncation","high":"discarded-only candidate or higher-quality filter-passing match","moderate":"discarded evidence changes filter, candidate, charge, or peak coverage","low":"discarded rows are downstream-ineligible"}[risk]
        rec=_recommend(risk,bool(df),nc,np);st=_strategy_counts(ranked,cutoff,cfg);strategies=[a+b for a,b in zip(strategies,st)]
        row={"Fragment_ID":f.fragment_id,"Fragment_Start":f.start,"Fragment_End":f.end,"Fragment_Sequence":f.sequence,"Fragment_Length":len(f.sequence or ""),
        "Missed_Cleavages":f.missed_cleavages,"Terminal_Form":f.terminal_form,"Theoretical_Neutral_Mass":f.unmodified_mass,
        "Pre_Truncation_Match_Count":len(ranked),"Configured_Max_Matches":cutoff,"Retained_Match_Count":len(kept),"Discarded_Match_Count":len(drop),
        "Truncation_Applied":bool(drop),"Retained_Fraction":len(kept)/len(ranked) if ranked else 1.0,**_stats(kept,eligible,"Retained",cfg),**_stats(drop,eligible,"Discarded",cfg),
        "Better_Error_In_Discarded":be,"Higher_Intensity_In_Discarded":bi,"Better_Peak_Tier_In_Discarded":bt,"Better_Confidence_In_Discarded":bc,
        "Filter_Passing_Discarded_While_Retained_Fails":bool(df and not kf),"New_Charge_State_Only_In_Discarded":nc,"New_Physical_Peak_Only_In_Discarded":np,
        "Modification_Candidate_Only_In_Discarded":bool(onlyc),"Candidate_Family_Only_In_Discarded":bool(onlyf),"Truncation_Risk_Level":risk,
        "Truncation_Risk_Reason":reason,"Recommended_Action":rec,"Applied_To_Final_Score":False};audit.append(row)
        we=max((abs(m.mass_error_ppm) for m in kept),default=None);wi=min((m.intensity for m in kept),default=None);wt=min((_tier(m.peak_tier) for m in kept),default=None);wc=min((_conf(m.confidence) for m in kept),default=None)
        for rank,m in enumerate(ranked,1):
            cs=cidx.get(_mkey(m),[]);mids=sorted({str(_raw(c).get("modification_id") or "") for c in cs});shifts=sorted({float(_raw(c).get("modification_mass_shift") or 0) for c in cs});discard=rank>cutoff
            detail.append({"Fragment_ID":f.fragment_id,"Pre_Truncation_Rank":rank,"Retained_Or_Discarded":"discarded" if discard else "retained","Retention_Cutoff":cutoff,
            "Retention_Sort_Key":f"({abs(m.mass_error_ppm):.12g}, {-float(m.intensity):.12g})","Physical_Peak_ID":_pid(m),"Peak_Index":getattr(m,"_audit_peak_index",""),
            "Observed_MZ":m.observed_mz,"Intensity":m.intensity,"Peak_Tier":m.peak_tier,"Confidence":m.confidence,"Charge":m.charge,
            "Observed_Neutral_Mass":m.fragment_mass+m.mass_error_da*abs(m.charge),"Theoretical_Neutral_Mass":m.fragment_mass,"Mass_Error_Da":m.mass_error_da,
            "Mass_Error_PPM":m.mass_error_ppm,"Abs_Mass_Error_PPM":abs(m.mass_error_ppm),"RT":m.rt,"Passes_Fragment_MS1_Filter":_passes(m,cfg),
            "Eligible_For_Modification_Search":bool(cs),"Modification_Mass_Shift":";".join(f"{x:.10g}" for x in shifts),"Possible_Modification_Count":len(cs),
            "Possible_Modification_IDs":";".join(mids),"Candidate_Key_List":";".join("|".join(map(str,_ckey(c))) for c in cs),
            "Candidate_Only_From_Discarded":bool(discard and any(_ckey(c) not in base_ck for c in cs)),
            "Better_Than_Worst_Retained_By_Error":bool(discard and we is not None and abs(m.mass_error_ppm)<we),
            "Better_Than_Worst_Retained_By_Intensity":bool(discard and wi is not None and m.intensity>wi),
            "Better_Than_Worst_Retained_By_Tier":bool(discard and wt is not None and _tier(m.peak_tier)>wt),
            "Better_Than_Worst_Retained_By_Confidence":bool(discard and wc is not None and _conf(m.confidence)>wc)})
        internal[str(f.fragment_id)]={"retained":kept,"discarded":drop,"discarded_filter":df,"discarded_candidates":dc,"risk":risk,"recommendation":rec}
    bm={_rkey(r):r for r in baseline_ranking};sm={_rkey(r):r for r in shadow_ranking};common=set(bm)&set(sm)
    keychg=set(bm)!=set(sm);scorechg=any(bm[k].get("Final_Score")!=sm[k].get("Final_Score") for k in common);confchg=any(bm[k].get("Final_Confidence")!=sm[k].get("Final_Confidence") for k in common);rankchg=any(bm[k].get("Rank")!=sm[k].get("Rank") for k in common)
    b50={k for k,v in bm.items() if int(v.get("Rank") or 999999)<=50};s50={k for k,v in sm.items() if int(v.get("Rank") or 999999)<=50};cnm={"36","37","38"};cnmchg={k for k in bm if k[2] in cnm}!={k for k in sm if k[2] in cnm}
    original=len(detail);limit=int((config.reporting or {}).get("max_excel_rows_per_sheet",100000) or 100000);detail=detail[:limit]
    high=any(r["Truncation_Risk_Level"]=="high" for r in audit);moderate=any(r["Truncation_Risk_Level"]=="moderate" for r in audit);trunc=any(r["Truncation_Applied"] for r in audit)
    overall="high" if keychg or scorechg or confchg or rankchg or b50!=s50 or cnmchg or high else "moderate" if moderate else "low" if trunc else "none"
    rec=_recommend(overall,any(r["Discarded_Filter_Passing_Count"] for r in audit),any(r["New_Charge_State_Only_In_Discarded"] for r in audit),any(r["New_Physical_Peak_Only_In_Discarded"] for r in audit))
    summary={"Total_Fragment_Count":len(fragments),"Fragments_With_Any_Match":sum(bool(r["Pre_Truncation_Match_Count"]) for r in audit),"Truncated_Fragment_Count":sum(r["Truncation_Applied"] for r in audit),
    "Total_Pre_Truncation_Matches":len(all_matches),"Total_Retained_Matches":len(baseline_matches),"Total_Discarded_Matches":len(all_matches)-len(baseline_matches),
    "Total_Retained_Unique_Physical_Peaks":len({_pid(m) for m in baseline_matches}),"Total_Discarded_Unique_Physical_Peaks":len({_pid(m) for x in internal.values() for m in x["discarded"]}),
    "Fragments_With_Better_Error_Discarded":sum(r["Better_Error_In_Discarded"] for r in audit),"Fragments_With_Higher_Intensity_Discarded":sum(r["Higher_Intensity_In_Discarded"] for r in audit),
    "Fragments_With_Better_Tier_Discarded":sum(r["Better_Peak_Tier_In_Discarded"] for r in audit),"Fragments_With_Better_Confidence_Discarded":sum(r["Better_Confidence_In_Discarded"] for r in audit),
    "Fragments_With_Filter_Passing_Discarded":sum(bool(r["Discarded_Filter_Passing_Count"]) for r in audit),"Fragments_With_Modification_Eligible_Discarded":sum(bool(r["Discarded_Modification_Eligible_Count"]) for r in audit),
    "Discarded_Only_Modification_Candidate_Count":len(shadow_ck-base_ck),"Discarded_Only_Candidate_Family_Count":len({_family(c,lookup) for c in shadow_candidates}-{_family(c,lookup) for c in baseline_candidates}),
    "Baseline_Max_Matches_Per_Fragment":cutoff,"Shadow_Max_Matches_Per_Fragment":"","Shadow_Truncation_Disabled":True,"Baseline_Fragment_MS1_Match_Count":len(baseline_matches),"Shadow_Fragment_MS1_Match_Count":len(all_matches),
    "Baseline_Filtered_Match_Count":sum(_passes(m,cfg) for m in baseline_matches),"Shadow_Filtered_Match_Count":sum(_passes(m,cfg) for m in all_matches),
    "Baseline_Modification_Candidate_Count":len(baseline_candidates),"Shadow_Modification_Candidate_Count":len(shadow_candidates),"Added_Shadow_Modification_Candidate_Count":len(shadow_ck-base_ck),
    "Removed_Shadow_Modification_Candidate_Count":len(base_ck-shadow_ck),"Changed_Shadow_Modification_Candidate_Count":len(base_ck^shadow_ck),
    "Baseline_Ranked_Candidate_Count":len(baseline_ranking),"Shadow_Ranked_Candidate_Count":len(shadow_ranking),"Candidate_Key_Set_Changed":keychg,
    "Final_Score_Would_Change":scorechg,"Final_Confidence_Would_Change":confchg,"Rank_Would_Change":rankchg,"Top50_Membership_Would_Change":b50!=s50,"cnm5U_Result_Would_Change":cnmchg,
    "Current_Sort_Filter_Passing_Retained":strategies[0],"Filter_First_Filter_Passing_Retained":strategies[1],"Unique_Physical_Peak_Filter_Passing_Retained":strategies[2],
    "Balanced_Per_Charge_Filter_Passing_Retained":strategies[3],"Tier_Then_Error_Filter_Passing_Retained":strategies[4],
    "Current_Truncation_Sort_Key":"abs(mass_error_ppm) ascending; intensity descending","Current_Sort_Tie_Break":"stable: charge ascending, then eligible peak input order",
    "Detail_Original_Row_Count":original,"Detail_Written_Row_Count":len(detail),"Detail_Truncated":len(detail)<original,"Detail_Truncation_Reason":"max_excel_rows_per_sheet" if len(detail)<original else "",
    "Shadow_Score_Recalculation_Status":"exact_existing_ranking_function","Shadow_Score_Recalculation_Reason":"","Overall_Truncation_Risk":overall,
    "Overall_Truncation_Conclusion":"Expanded shadow changes formal-equivalent candidate or ranking results." if keychg or scorechg or confchg or rankchg or b50!=s50 or cnmchg else "Candidate and ranking results are unchanged, but higher-quality filter-passing support is discarded." if high else "Review-relevant support changes without a ranking change." if overall=="moderate" else "No material downstream effect.",
    "Recommended_Max_Matches_Per_Fragment":cutoff if overall in {"none","low"} else "review_before_change","Recommended_Next_Action":rec,
    "Audit_Mode":"shadow_unlimited_pre_truncation_capture","Applied_To_Final_Score":False}
    return {"audit_rows":[{c:r.get(c,"") for c in AUDIT_COLUMNS} for r in audit],"detail_rows":[{c:r.get(c,"") for c in DETAIL_COLUMNS} for r in detail],
    "summary_rows":[{c:summary.get(c,"") for c in SUMMARY_COLUMNS}],"summary":summary,"shadow_ranking":shadow_ranking,"fragment_internal":internal,"candidate_index":cidx,"baseline_rank_map":bm,"shadow_rank_map":sm}

def append_top_shadow_columns(rows,audit):
    is_frame=isinstance(rows,pd.DataFrame);source=rows.to_dict("records") if is_frame else list(rows or [])
    original_columns=list(rows.columns) if is_frame else list(source[0]) if source else []
    out=[];frags=audit["fragment_internal"];sm=audit["shadow_rank_map"]
    for original in source:
        row=dict(original);fid=str(row.get("Parent_Fragment_ID") or "");d=frags.get(fid,{"retained":[],"discarded":[],"discarded_filter":[],"discarded_candidates":[],"risk":"none","recommendation":"no_action_needed"});shadow=sm.get(_rkey(row),{});mid=str(row.get("Modification_ID") or "")
        vals={"MS1_Truncation_Affected":bool(d["discarded"]),"MS1_Truncated_Fragment_Count":int(bool(d["discarded"])),"MS1_Discarded_Support_Match_Count":len(d["discarded"]),
        "MS1_Discarded_Filter_Passing_Count":len(d["discarded_filter"]),"MS1_Discarded_Modification_Support_Count":sum(str(_raw(c).get("modification_id") or "")==mid for c in d["discarded_candidates"]),
        "Shadow_Expanded_MS1_Support_Count":len(d["retained"])+len(d["discarded"]),"Shadow_Expanded_Final_Score":shadow.get("Final_Score",""),
        "Shadow_Expanded_Final_Confidence":shadow.get("Final_Confidence",""),"Shadow_Expanded_Rank":shadow.get("Rank",""),
        "Shadow_Rank_Changed":bool(shadow and shadow.get("Rank")!=audit["baseline_rank_map"].get(_rkey(row),{}).get("Rank")),
        "MS1_Truncation_Risk":d["risk"],"MS1_Truncation_Recommendation":d["recommendation"],"MS1_Truncation_Applied_To_Final_Score":False}
        row.update(vals);out.append(row)
    return pd.DataFrame(out,columns=original_columns+TOP_COLUMNS) if is_frame else out

def append_diagnostic_shadow_columns(rows,audit):
    s=audit["summary"];vals={"MS1_Truncation_Audit_Available":True,"MS1_Truncation_Total_Fragments":s["Total_Fragment_Count"],"MS1_Truncated_Fragment_Count":s["Truncated_Fragment_Count"],
    "MS1_Pre_Truncation_Match_Count":s["Total_Pre_Truncation_Matches"],"MS1_Retained_Match_Count":s["Total_Retained_Matches"],"MS1_Discarded_Match_Count":s["Total_Discarded_Matches"],
    "MS1_Discarded_Filter_Passing_Count":sum(r["Discarded_Filter_Passing_Count"] for r in audit["audit_rows"]),"MS1_Discarded_Modification_Eligible_Count":sum(r["Discarded_Modification_Eligible_Count"] for r in audit["audit_rows"]),
    "MS1_Discarded_Only_Candidate_Count":s["Discarded_Only_Modification_Candidate_Count"],"MS1_Shadow_Candidate_Key_Set_Changed":s["Candidate_Key_Set_Changed"],
    "MS1_Shadow_Final_Score_Changed":s["Final_Score_Would_Change"],"MS1_Shadow_Rank_Changed":s["Rank_Would_Change"],"MS1_Truncation_Overall_Risk":s["Overall_Truncation_Risk"],
    "MS1_Truncation_Recommendation":s["Recommended_Next_Action"],"MS1_Truncation_Applied_To_Final_Score":False}
    is_frame=isinstance(rows,pd.DataFrame);source=rows.to_dict("records") if is_frame else list(rows or [{}])
    out=[dict(r,**vals) for r in source]
    return pd.DataFrame(out,columns=list(rows.columns)+DIAGNOSTIC_COLUMNS) if is_frame else out
