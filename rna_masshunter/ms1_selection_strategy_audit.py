"""Non-mutating A/B audit of Fragment MS1 selection strategies."""
from __future__ import annotations
from collections import Counter, defaultdict
from dataclasses import asdict, is_dataclass
from statistics import median
from typing import Any
import pandas as pd
from rna_masshunter.evidence_ranking import build_modification_evidence_ranking
from rna_masshunter.modification_search import search_known_modifications
from rna_masshunter.ms1_mapping import retention_sort_key

STRATEGIES=("current","filter_first","tier_then_error","unlimited")
STRATEGY_COLUMNS="""Fragment_ID Fragment_Start Fragment_End Fragment_Sequence Fragment_Length Strategy Pre_Truncation_Match_Count Pre_Filter_Passing_Count Selection_Limit Selected_Match_Count Selected_Filter_Passing_Count Selected_Filter_Failing_Count Selected_Unique_Physical_Peak_Count Selected_Charge_State_Count Selected_Major_Count Selected_Minor_Count Selected_Trace_Count Selected_High_Confidence_Count Selected_Medium_Confidence_Count Selected_Low_Confidence_Count Selected_Min_Abs_Error_PPM Selected_Median_Abs_Error_PPM Selected_Max_Intensity Selected_Median_Intensity Added_Vs_Current_Count Removed_Vs_Current_Count Added_Filter_Passing_Count Removed_Filter_Passing_Count Added_Unique_Physical_Peak_Count Added_Charge_State_Count Added_Modification_Eligible_Count Added_Modification_Candidate_Count Added_Candidate_Family_Count Candidate_Key_Set_Changed Candidate_Support_Changed Shadow_Final_Score_Changed Shadow_Final_Confidence_Changed Shadow_Rank_Changed Top50_Membership_Changed cnm5U_Result_Changed Strategy_Risk Strategy_Recommendation Applied_To_Formal_Result""".split()
DETAIL_COLUMNS="""Fragment_ID Strategy Pre_Truncation_Rank_Current Strategy_Rank Selected Selected_Vs_Current_Status Physical_Peak_ID Peak_Index Observed_MZ Charge Intensity Peak_Tier Confidence Observed_Neutral_Mass Theoretical_Neutral_Mass Mass_Error_Da Mass_Error_PPM Abs_Mass_Error_PPM RT Passes_Fragment_MS1_Filter Eligible_For_Modification_Search Possible_Modification_Count Possible_Modification_IDs Candidate_Key_List Candidate_Family_List Used_By_Current Used_By_Filter_First Used_By_Tier_Then_Error Used_By_Unlimited""".split()
SUMMARY_COLUMNS="""Total_Fragment_Count Fragments_With_Matches Fragments_Affected_By_Limit Baseline_Max_Matches_Per_Fragment Current_Selected_Count Current_Filter_Passing_Count Filter_First_Selected_Count Filter_First_Filter_Passing_Count Tier_Then_Error_Selected_Count Tier_Then_Error_Filter_Passing_Count Unlimited_Selected_Count Unlimited_Filter_Passing_Count Filter_First_Added_Filter_Passing_Count Tier_Then_Error_Added_Filter_Passing_Count Unlimited_Added_Filter_Passing_Count Filter_First_Recovery_Fraction_Of_Unlimited Tier_Then_Error_Recovery_Fraction_Of_Unlimited Current_Modification_Candidate_Count Filter_First_Modification_Candidate_Count Tier_Then_Error_Modification_Candidate_Count Unlimited_Modification_Candidate_Count Filter_First_Added_Candidate_Count Tier_Then_Error_Added_Candidate_Count Unlimited_Added_Candidate_Count Filter_First_Candidate_Key_Set_Changed Tier_Then_Error_Candidate_Key_Set_Changed Unlimited_Candidate_Key_Set_Changed Filter_First_Final_Score_Changed_Candidate_Count Tier_Then_Error_Final_Score_Changed_Candidate_Count Unlimited_Final_Score_Changed_Candidate_Count Filter_First_Final_Confidence_Changed_Candidate_Count Tier_Then_Error_Final_Confidence_Changed_Candidate_Count Unlimited_Final_Confidence_Changed_Candidate_Count Filter_First_Rank_Changed_Candidate_Count Tier_Then_Error_Rank_Changed_Candidate_Count Unlimited_Rank_Changed_Candidate_Count Filter_First_Top50_Membership_Changed Tier_Then_Error_Top50_Membership_Changed Unlimited_Top50_Membership_Changed Filter_First_cnm5U_Result_Changed Tier_Then_Error_cnm5U_Result_Changed Unlimited_cnm5U_Result_Changed Recommended_Formal_Strategy Recommended_Max_Matches Evidence_For_Recommendation Remaining_Risk Formal_Change_Ready Required_Additional_Validation Detail_Original_Row_Count Detail_Written_Row_Count Detail_Truncated Detail_Truncation_Reason Strategy_Audit_Additional_Time_Seconds Strategy_Audit_Peak_Tracked_Memory_MiB Audit_Mode Applied_To_Formal_Result""".split()
TOP_COLUMNS="""MS1_Selection_Strategy_Affected Current_MS1_Filtered_Support Filter_First_MS1_Filtered_Support Tier_Then_Error_MS1_Filtered_Support Unlimited_MS1_Filtered_Support Filter_First_Support_Delta Tier_Then_Error_Support_Delta Unlimited_Support_Delta Filter_First_Shadow_Final_Score Tier_Then_Error_Shadow_Final_Score Unlimited_Shadow_Final_Score Filter_First_Shadow_Rank Tier_Then_Error_Shadow_Rank Unlimited_Shadow_Rank Recommended_MS1_Selection_Strategy MS1_Selection_Applied_To_Formal_Result""".split()
DIAGNOSTIC_COLUMNS="""MS1_Selection_Audit_Available MS1_Current_Filter_Passing_Count MS1_Filter_First_Filter_Passing_Count MS1_Tier_Then_Error_Filter_Passing_Count MS1_Unlimited_Filter_Passing_Count MS1_Filter_First_Recovery_Fraction MS1_Tier_Then_Error_Recovery_Fraction MS1_Filter_First_Candidate_Set_Changed MS1_Tier_Then_Error_Candidate_Set_Changed MS1_Unlimited_Candidate_Set_Changed MS1_Filter_First_Rank_Changed MS1_Tier_Then_Error_Rank_Changed MS1_Unlimited_Rank_Changed MS1_Recommended_Selection_Strategy MS1_Formal_Change_Ready MS1_Selection_Applied_To_Formal_Result""".split()

def _raw(x): return asdict(x) if is_dataclass(x) else dict(x)
def _tier(x): return {"major":3,"minor":2,"trace":1}.get(str(x or "").lower(),0)
def _conf(x): return {"high":3,"medium":2,"low":1}.get(str(x or "").lower(),0)
def _pid(m): return str(getattr(m,"_audit_physical_peak_id","") or f"PK_{getattr(m,'_audit_peak_index',0):06d}_{m.observed_mz:.8f}")
def _mkey(m): return (str(m.fragment_id),_pid(m),int(m.charge),round(float(m.observed_mz),10))
def _ckey(c):
 r=_raw(c);return (str(r.get("source_type") or ""),str(r.get("source_id") or ""),int(r.get("charge") or 0),round(float(r.get("observed_mz") or 0),10),str(r.get("modification_id") or ""))
def _rkey(r):
 p=r.get("Candidate_tRNA_Position",r.get("Candidate_Positions_In_tRNA",""));p=int(p) if isinstance(p,float) and p.is_integer() else p
 return (str(r.get("Modification_ID") or ""),str(r.get("Parent_Fragment_ID") or ""),str(p or ""))
def _family(c,lookup):
 r=_raw(c);mid=str(r.get("modification_id") or "");m=lookup.get(mid);raw=getattr(m,"raw",{}) or {}
 return str(raw.get("chemical_group") or raw.get("near_isobaric_group") or getattr(m,"category","") or mid)

def passes_fragment_ms1_filter(m,cfg):
 try: n=int(cfg.get("min_fragment_length_for_filtered",3) or 3)
 except (TypeError,ValueError): n=3
 if n<1:n=3
 ts=cfg.get("filtered_peak_tiers",["Major","Minor"]);cs=cfg.get("filtered_confidence",["High","Medium"])
 ts=[ts] if isinstance(ts,str) else (ts or []);cs=[cs] if isinstance(cs,str) else (cs or [])
 tiers={str(x).lower() for x in ts};confs={str(x).lower() for x in cs}
 return len(m.sequence or "")>=n and (not tiers or str(m.peak_tier or "").lower() in tiers) and (not confs or str(m.confidence or "").lower() in confs)

def tier_then_error_sort_key(m,cfg):
 return (not passes_fragment_ms1_filter(m,cfg),-_tier(m.peak_tier),-_conf(m.confidence),abs(float(m.mass_error_ppm)),-float(m.intensity or 0),int(m.charge),_pid(m),int(getattr(m,"_audit_generation_order",0) or 0))
def select_strategy_matches(ms,strategy,limit,cfg):
 if strategy=="current":return list(ms[:limit])
 passing=[m for m in ms if passes_fragment_ms1_filter(m,cfg)]
 if strategy=="filter_first":return sorted(passing,key=retention_sort_key)[:limit]
 if strategy=="tier_then_error":return sorted(ms,key=lambda m:tier_then_error_sort_key(m,cfg))[:limit]
 if strategy=="unlimited":return passing
 raise ValueError(f"Unknown MS1 selection strategy: {strategy}")
def _change(base,shadow):
 bm={_rkey(r):r for r in base};sm={_rkey(r):r for r in shadow};common=set(bm)&set(sm)
 top=lambda x:{k for k,v in x.items() if int(v.get("Rank") or 999999)<=50};cnm=lambda x:{k for k in x if k[2] in {"36","37","38"}}
 return {"map":sm,"set":set(bm)!=set(sm),"score":sum(bm[k].get("Final_Score")!=sm[k].get("Final_Score") for k in common),"confidence":sum(bm[k].get("Final_Confidence")!=sm[k].get("Final_Confidence") for k in common),"rank":sum(bm[k].get("Rank")!=sm[k].get("Rank") for k in common),"top50":top(bm)!=top(sm),"cnm":cnm(bm)!=cnm(sm)}
def _stats(ms):
 t=Counter(str(m.peak_tier or "").lower() for m in ms);c=Counter(str(m.confidence or "").lower() for m in ms);e=[abs(m.mass_error_ppm) for m in ms];i=[float(m.intensity or 0) for m in ms]
 return {"Selected_Unique_Physical_Peak_Count":len({_pid(m) for m in ms}),"Selected_Charge_State_Count":len({m.charge for m in ms}),"Selected_Major_Count":t["major"],"Selected_Minor_Count":t["minor"],"Selected_Trace_Count":t["trace"],"Selected_High_Confidence_Count":c["high"],"Selected_Medium_Confidence_Count":c["medium"],"Selected_Low_Confidence_Count":c["low"],"Selected_Min_Abs_Error_PPM":min(e) if e else "","Selected_Median_Abs_Error_PPM":median(e) if e else "","Selected_Max_Intensity":max(i) if i else "","Selected_Median_Intensity":median(i) if i else ""}
def _recovery(a,b,u): return (1.0 if b>=u else 0.0) if u<=a else max(0,min(1,(b-a)/(u-a)))

def build_ms1_selection_strategy_audit(context,config,modifications,intact_results,baseline_matches,baseline_candidates,baseline_ranking,ms2_results,rule_set=None,pathways=None):
 cfg=config.fragment_mapping or {};limit=int(context.get("configured_max_matches") or cfg.get("max_matches_per_fragment") or 20);limit=limit if limit>0 else 20
 requested=cfg.get("MS1_Shadow_Selection_Strategies",STRATEGIES) or STRATEGIES;requested=[requested] if isinstance(requested,str) else requested
 strategies=tuple(s for s in STRATEGIES if s in {str(x).lower() for x in requested})
 for s in ("current","filter_first","unlimited"):
  if s not in strategies:strategies+=(s,)
 items=list(context.get("fragments") or []);selected={s:{} for s in strategies}
 for item in items:
  fid=str(item["fragment"].fragment_id);ranked=list(item.get("ranked_matches") or [])
  for s in strategies:selected[s][fid]=select_strategy_matches(ranked,s,limit,cfg)
 matches={s:[m for group in selected[s].values() for m in group] for s in strategies}
 candidates={"current":list(baseline_candidates)};rankings={"current":list(baseline_ranking)}
 for s in strategies:
  if s=="current":continue
  candidates[s]=search_known_modifications(matches[s],intact_results,modifications,config,warnings=None)
  rankings[s],_=build_modification_evidence_ranking(config,modifications,[x["fragment"] for x in items],matches[s],candidates[s],ms2_results,rule_set=rule_set,pathways=pathways)
 changes={s:_change(baseline_ranking,rankings[s]) for s in strategies};base_rmap={_rkey(r):r for r in baseline_ranking};lookup={str(getattr(m,"id","")):m for m in modifications}
 cidx={};fkeys={};ffamilies={}
 for s in strategies:
  cidx[s]=defaultdict(list);fkeys[s]=defaultdict(set);ffamilies[s]=defaultdict(set)
  for c in candidates[s]:
   r=_raw(c)
   if str(r.get("source_type") or "").lower()!="fragment":continue
   k=(str(r.get("source_id") or ""),int(r.get("charge") or 0),round(float(r.get("observed_mz") or 0),10));cidx[s][k].append(c);fkeys[s][k[0]].add(_ckey(c));ffamilies[s][k[0]].add(_family(c,lookup))
 rows=[];detail=[];internal={}
 for item in items:
  f=item["fragment"];fid=str(f.fragment_id);ranked=list(item.get("ranked_matches") or []);cur=selected["current"][fid];curkeys={_mkey(m) for m in cur};used={s:{_mkey(m) for m in selected[s][fid]} for s in strategies};internal[fid]={"ranked":ranked,"selected":{s:selected[s][fid] for s in strategies}}
  for s in strategies:
   ms=selected[s][fid];keys=used[s];added=keys-curkeys;removed=curkeys-keys;am=[m for m in ms if _mkey(m) in added];rm=[m for m in cur if _mkey(m) in removed];cks=fkeys[s][fid];basecks=fkeys["current"][fid];common={k for k in changes[s]["map"] if k[1]==fid}&{k for k in base_rmap if k[1]==fid}
   score=any(changes[s]["map"][k].get("Final_Score")!=base_rmap[k].get("Final_Score") for k in common);conf=any(changes[s]["map"][k].get("Final_Confidence")!=base_rmap[k].get("Final_Confidence") for k in common);rank=any(changes[s]["map"][k].get("Rank")!=base_rmap[k].get("Rank") for k in common);support=keys!=curkeys;keychg=cks!=basecks
   risk,rec=("baseline","retain_current_for_formal_result") if s=="current" else (("high","additional_validation_required") if keychg or score or conf or rank else (("moderate","filter_before_truncation") if support else ("none","equivalent_to_current")))
   passing=sum(passes_fragment_ms1_filter(m,cfg) for m in ms)
   rows.append({"Fragment_ID":fid,"Fragment_Start":f.start,"Fragment_End":f.end,"Fragment_Sequence":f.sequence,"Fragment_Length":len(f.sequence or ""),"Strategy":s,"Pre_Truncation_Match_Count":len(ranked),"Pre_Filter_Passing_Count":sum(passes_fragment_ms1_filter(m,cfg) for m in ranked),"Selection_Limit":"" if s=="unlimited" else limit,"Selected_Match_Count":len(ms),"Selected_Filter_Passing_Count":passing,"Selected_Filter_Failing_Count":len(ms)-passing,**_stats(ms),"Added_Vs_Current_Count":len(added),"Removed_Vs_Current_Count":len(removed),"Added_Filter_Passing_Count":sum(passes_fragment_ms1_filter(m,cfg) for m in am),"Removed_Filter_Passing_Count":sum(passes_fragment_ms1_filter(m,cfg) for m in rm),"Added_Unique_Physical_Peak_Count":len({_pid(m) for m in am}-{_pid(m) for m in cur}),"Added_Charge_State_Count":len({m.charge for m in am}-{m.charge for m in cur}),"Added_Modification_Eligible_Count":sum(bool(cidx[s].get((fid,int(m.charge),round(float(m.observed_mz),10)))) for m in am),"Added_Modification_Candidate_Count":len(cks-basecks),"Added_Candidate_Family_Count":len(ffamilies[s][fid]-ffamilies["current"][fid]),"Candidate_Key_Set_Changed":keychg,"Candidate_Support_Changed":support,"Shadow_Final_Score_Changed":score,"Shadow_Final_Confidence_Changed":conf,"Shadow_Rank_Changed":rank,"Top50_Membership_Changed":changes[s]["top50"],"cnm5U_Result_Changed":changes[s]["cnm"],"Strategy_Risk":risk,"Strategy_Recommendation":rec,"Applied_To_Formal_Result":False})
  currank={_mkey(m):i for i,m in enumerate(ranked,1)};sranks={s:{_mkey(m):i for i,m in enumerate(selected[s][fid],1)} for s in strategies}
  for s in strategies:
   for m in ranked:
    k=_mkey(m);chosen=k in used[s];status="retained_in_both" if chosen and k in curkeys else "added_by_strategy" if chosen else "removed_vs_current" if k in curkeys else "not_selected";cs=cidx[s].get((fid,int(m.charge),round(float(m.observed_mz),10)),[])
    detail.append({"Fragment_ID":fid,"Strategy":s,"Pre_Truncation_Rank_Current":currank[k],"Strategy_Rank":sranks[s].get(k,""),"Selected":chosen,"Selected_Vs_Current_Status":status,"Physical_Peak_ID":_pid(m),"Peak_Index":getattr(m,"_audit_peak_index",""),"Observed_MZ":m.observed_mz,"Charge":m.charge,"Intensity":m.intensity,"Peak_Tier":m.peak_tier,"Confidence":m.confidence,"Observed_Neutral_Mass":m.fragment_mass+m.mass_error_da*abs(m.charge),"Theoretical_Neutral_Mass":m.fragment_mass,"Mass_Error_Da":m.mass_error_da,"Mass_Error_PPM":m.mass_error_ppm,"Abs_Mass_Error_PPM":abs(m.mass_error_ppm),"RT":m.rt,"Passes_Fragment_MS1_Filter":passes_fragment_ms1_filter(m,cfg),"Eligible_For_Modification_Search":bool(cs),"Possible_Modification_Count":len(cs),"Possible_Modification_IDs":";".join(sorted({str(_raw(c).get('modification_id') or '') for c in cs})),"Candidate_Key_List":";".join(sorted("|".join(map(str,_ckey(c))) for c in cs)),"Candidate_Family_List":";".join(sorted({_family(c,lookup) for c in cs})),"Used_By_Current":k in used.get("current",set()),"Used_By_Filter_First":k in used.get("filter_first",set()),"Used_By_Tier_Then_Error":k in used.get("tier_then_error",set()),"Used_By_Unlimited":k in used.get("unlimited",set())})
 counts={s:{"selected":len(matches[s]),"passing":sum(passes_fragment_ms1_filter(m,cfg) for m in matches[s])} for s in strategies};a=counts["current"]["passing"];u=counts["unlimited"]["passing"];baseck={_ckey(c) for c in candidates["current"]}
 summary={"Total_Fragment_Count":len(items),"Fragments_With_Matches":sum(bool(x.get("ranked_matches")) for x in items),"Fragments_Affected_By_Limit":sum(len(x.get("ranked_matches") or [])>limit for x in items),"Baseline_Max_Matches_Per_Fragment":limit,"Current_Selected_Count":counts["current"]["selected"],"Current_Filter_Passing_Count":a,"Filter_First_Selected_Count":counts["filter_first"]["selected"],"Filter_First_Filter_Passing_Count":counts["filter_first"]["passing"],"Tier_Then_Error_Selected_Count":counts.get("tier_then_error",{}).get("selected",0),"Tier_Then_Error_Filter_Passing_Count":counts.get("tier_then_error",{}).get("passing",0),"Unlimited_Selected_Count":counts["unlimited"]["selected"],"Unlimited_Filter_Passing_Count":u,"Filter_First_Added_Filter_Passing_Count":counts["filter_first"]["passing"]-a,"Tier_Then_Error_Added_Filter_Passing_Count":counts.get("tier_then_error",{}).get("passing",0)-a,"Unlimited_Added_Filter_Passing_Count":u-a,"Filter_First_Recovery_Fraction_Of_Unlimited":_recovery(a,counts["filter_first"]["passing"],u),"Tier_Then_Error_Recovery_Fraction_Of_Unlimited":_recovery(a,counts.get("tier_then_error",{}).get("passing",0),u),"Current_Modification_Candidate_Count":len(candidates["current"]),"Filter_First_Modification_Candidate_Count":len(candidates["filter_first"]),"Tier_Then_Error_Modification_Candidate_Count":len(candidates.get("tier_then_error",[])),"Unlimited_Modification_Candidate_Count":len(candidates["unlimited"])}
 for s,label in (("filter_first","Filter_First"),("tier_then_error","Tier_Then_Error"),("unlimited","Unlimited")):
  if s not in strategies:continue
  ck={_ckey(c) for c in candidates[s]};summary.update({f"{label}_Added_Candidate_Count":len(ck-baseck),f"{label}_Candidate_Key_Set_Changed":ck!=baseck,f"{label}_Final_Score_Changed_Candidate_Count":changes[s]["score"],f"{label}_Final_Confidence_Changed_Candidate_Count":changes[s]["confidence"],f"{label}_Rank_Changed_Candidate_Count":changes[s]["rank"],f"{label}_Top50_Membership_Changed":changes[s]["top50"],f"{label}_cnm5U_Result_Changed":changes[s]["cnm"]})
 summary.update({"Recommended_Formal_Strategy":"filter_first","Recommended_Max_Matches":limit,"Evidence_For_Recommendation":f"filter_first recovers {summary['Filter_First_Added_Filter_Passing_Count']} filter-passing rows ({summary['Filter_First_Recovery_Fraction_Of_Unlimited']:.3f} of unlimited recoverable evidence)","Remaining_Risk":"single-dataset validation; physical-peak duplication and cross-dataset stability require confirmation","Formal_Change_Ready":False,"Required_Additional_Validation":"independent RNase_A/Tko and RNase_T1 datasets; physical-peak duplication policy","Audit_Mode":"shadow_selection_strategy_ab","Applied_To_Formal_Result":False})
 original=len(detail);maxrows=int((config.reporting or {}).get("max_excel_rows_per_sheet",100000) or 100000);detail=detail[:maxrows];summary.update({"Detail_Original_Row_Count":original,"Detail_Written_Row_Count":len(detail),"Detail_Truncated":len(detail)<original,"Detail_Truncation_Reason":"max_excel_rows_per_sheet; deterministic order" if len(detail)<original else ""})
 return {"strategy_rows":[{c:r.get(c,"") for c in STRATEGY_COLUMNS} for r in rows],"detail_rows":[{c:r.get(c,"") for c in DETAIL_COLUMNS} for r in detail],"summary_rows":[{c:summary.get(c,"") for c in SUMMARY_COLUMNS}],"summary":summary,"strategy_matches":matches,"strategy_candidates":candidates,"strategy_rankings":rankings,"ranking_changes":changes,"fragment_internal":internal,"base_ranking_map":base_rmap}

def append_top_selection_columns(rows,audit):
 frame=isinstance(rows,pd.DataFrame);source=rows.to_dict("records") if frame else list(rows or []);cols=list(rows.columns) if frame else (list(source[0]) if source else []);lookup={(r["Fragment_ID"],r["Strategy"]):r for r in audit["strategy_rows"]};out=[]
 for original in source:
  row=dict(original);fid=str(row.get("Parent_Fragment_ID") or "");support={s:lookup.get((fid,s),{}).get("Selected_Filter_Passing_Count",0) for s in STRATEGIES};key=_rkey(row);maps={s:audit["ranking_changes"].get(s,{}).get("map",{}) for s in STRATEGIES}
  row.update({"MS1_Selection_Strategy_Affected":any(support[s]!=support["current"] for s in STRATEGIES[1:]),"Current_MS1_Filtered_Support":support["current"],"Filter_First_MS1_Filtered_Support":support["filter_first"],"Tier_Then_Error_MS1_Filtered_Support":support["tier_then_error"],"Unlimited_MS1_Filtered_Support":support["unlimited"],"Filter_First_Support_Delta":support["filter_first"]-support["current"],"Tier_Then_Error_Support_Delta":support["tier_then_error"]-support["current"],"Unlimited_Support_Delta":support["unlimited"]-support["current"],"Filter_First_Shadow_Final_Score":maps["filter_first"].get(key,{}).get("Final_Score",""),"Tier_Then_Error_Shadow_Final_Score":maps["tier_then_error"].get(key,{}).get("Final_Score",""),"Unlimited_Shadow_Final_Score":maps["unlimited"].get(key,{}).get("Final_Score",""),"Filter_First_Shadow_Rank":maps["filter_first"].get(key,{}).get("Rank",""),"Tier_Then_Error_Shadow_Rank":maps["tier_then_error"].get(key,{}).get("Rank",""),"Unlimited_Shadow_Rank":maps["unlimited"].get(key,{}).get("Rank",""),"Recommended_MS1_Selection_Strategy":audit["summary"]["Recommended_Formal_Strategy"],"MS1_Selection_Applied_To_Formal_Result":False});out.append(row)
 return pd.DataFrame(out,columns=cols+TOP_COLUMNS) if frame else out

def append_selection_diagnostic_columns(rows,audit):
 s=audit["summary"];vals={"MS1_Selection_Audit_Available":True,"MS1_Current_Filter_Passing_Count":s["Current_Filter_Passing_Count"],"MS1_Filter_First_Filter_Passing_Count":s["Filter_First_Filter_Passing_Count"],"MS1_Tier_Then_Error_Filter_Passing_Count":s["Tier_Then_Error_Filter_Passing_Count"],"MS1_Unlimited_Filter_Passing_Count":s["Unlimited_Filter_Passing_Count"],"MS1_Filter_First_Recovery_Fraction":s["Filter_First_Recovery_Fraction_Of_Unlimited"],"MS1_Tier_Then_Error_Recovery_Fraction":s["Tier_Then_Error_Recovery_Fraction_Of_Unlimited"],"MS1_Filter_First_Candidate_Set_Changed":s["Filter_First_Candidate_Key_Set_Changed"],"MS1_Tier_Then_Error_Candidate_Set_Changed":s["Tier_Then_Error_Candidate_Key_Set_Changed"],"MS1_Unlimited_Candidate_Set_Changed":s["Unlimited_Candidate_Key_Set_Changed"],"MS1_Filter_First_Rank_Changed":bool(s["Filter_First_Rank_Changed_Candidate_Count"]),"MS1_Tier_Then_Error_Rank_Changed":bool(s["Tier_Then_Error_Rank_Changed_Candidate_Count"]),"MS1_Unlimited_Rank_Changed":bool(s["Unlimited_Rank_Changed_Candidate_Count"]),"MS1_Recommended_Selection_Strategy":s["Recommended_Formal_Strategy"],"MS1_Formal_Change_Ready":s["Formal_Change_Ready"],"MS1_Selection_Applied_To_Formal_Result":False};frame=isinstance(rows,pd.DataFrame);source=rows.to_dict("records") if frame else list(rows or [{}]);out=[dict(r,**vals) for r in source]
 return pd.DataFrame(out,columns=list(rows.columns)+DIAGNOSTIC_COLUMNS) if frame else out
