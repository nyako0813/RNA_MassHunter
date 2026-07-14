"""Full non-mutating tier-top50 and physical-peak deduplication shadow audit."""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, is_dataclass
from statistics import median
import time
from typing import Any

import pandas as pd

from rna_masshunter.evidence_ranking import build_modification_evidence_ranking
from rna_masshunter.modification_search import search_known_modifications, summarize_known_modification_candidates
from rna_masshunter.ms1_selection_strategy_audit import passes_fragment_ms1_filter, tier_then_error_sort_key
from rna_masshunter.review_dashboard import build_review_dashboard_results

DEDUP_MODES = ("no_dedup", "exact_physical_peak_dedup", "physical_peak_charge_dedup", "fragment_then_physical_peak_dedup", "global_physical_peak_dedup")
TOP50_COLUMNS = """Dataset_ID Enzyme Fragment_ID Fragment_Start Fragment_End Fragment_Sequence Fragment_Length Strategy Selection_Limit Dedup_Mode Pre_Truncation_Match_Count Pre_Filter_Passing_Count Selected_Match_Count Selected_Filter_Passing_Count Unique_Physical_Peak_Count Unique_Charge_Count Duplicate_Row_Count Duplicate_Physical_Peak_Count Multi_Charge_Physical_Peak_Count Multi_Fragment_Shared_Peak_Count Major_Count Minor_Count Trace_Count High_Confidence_Count Medium_Confidence_Count Low_Confidence_Count Min_Abs_Error_PPM Median_Abs_Error_PPM Max_Intensity Median_Intensity Modification_Candidate_Count Candidate_Key_Count Added_Candidate_Count_Vs_Current Removed_Candidate_Count_Vs_Current Candidate_Support_Changed Shadow_Final_Score_Changed Shadow_Final_Confidence_Changed Shadow_Rank_Changed Shadow_Top50_Changed cnm5U_Result_Changed Recommended_Dedup_Mode Applied_To_Formal_Result""".split()
DETAIL_COLUMNS = """Dataset_ID Enzyme Fragment_ID Physical_Peak_ID Spectrum_ID Scan_Index Peak_Index RT Observed_MZ Intensity Charge Observed_Neutral_Mass Theoretical_Neutral_Mass Mass_Error_Da Mass_Error_PPM Abs_Mass_Error_PPM Peak_Tier Confidence Passes_Fragment_MS1_Filter Tier_Top20_Selected Tier_Top50_Selected Tier_Unlimited_Selected Exact_Physical_Duplicate_Group_Size Charge_Assignment_Count_For_Peak Fragment_Assignment_Count_For_Peak Theoretical_Form_Count_For_Peak Same_Fragment_Multi_Charge Cross_Fragment_Shared Near_MZ_Group_ID Kept_No_Dedup Kept_Exact_Peak_Dedup Kept_Charge_Dedup Kept_Fragment_Peak_Dedup Kept_Global_Dedup Dedup_Selection_Reason Possible_Modification_Count Candidate_Key_List Candidate_Family_List""".split()
SUMMARY_COLUMNS = """Dataset_ID Enzyme Sequence_Name Summary_Scope Fragment_Length_Group Fragment_Count Total_Matches Filter_Passing_Matches Unique_Physical_Peaks Exact_Duplicate_Rows Multi_Charge_Peaks Multi_Fragment_Shared_Peaks Near_MZ_Groups Short_Fragment_Match_Count Short_Fragment_Duplicate_Fraction Average_Charge_Assignments_Per_Peak Average_Fragment_Assignments_Per_Peak Tier_Top20_Filter_Passing Tier_Top50_Filter_Passing Tier_Unlimited_Filter_Passing Top50_Recovery_Fraction Top50_Added_Major Top50_Added_Minor Top50_Added_High Top50_Added_Medium No_Dedup_Selected_Count Exact_Peak_Dedup_Selected_Count Charge_Dedup_Selected_Count Fragment_Peak_Dedup_Selected_Count Global_Dedup_Selected_Count No_Dedup_Filter_Passing Exact_Peak_Dedup_Filter_Passing Charge_Dedup_Filter_Passing Fragment_Peak_Dedup_Filter_Passing Global_Dedup_Filter_Passing Top50_Modification_Candidate_Count Top50_Candidate_Key_Set_Changed Top50_Final_Score_Changed_Count Top50_Final_Confidence_Changed_Count Top50_Rank_Changed_Count Top50_Membership_Changed Exact_Dedup_Candidate_Key_Set_Changed Exact_Dedup_Final_Score_Changed_Count Exact_Dedup_Final_Confidence_Changed_Count Exact_Dedup_Rank_Changed_Count Exact_Dedup_Top50_Changed Charge_Dedup_Candidate_Key_Set_Changed Charge_Dedup_Final_Score_Changed_Count Charge_Dedup_Final_Confidence_Changed_Count Charge_Dedup_Rank_Changed_Count Charge_Dedup_Top50_Changed Fragment_Dedup_Candidate_Key_Set_Changed Fragment_Dedup_Final_Score_Changed_Count Fragment_Dedup_Final_Confidence_Changed_Count Fragment_Dedup_Rank_Changed_Count Fragment_Dedup_Top50_Changed Global_Dedup_Candidate_Key_Set_Changed Global_Dedup_Final_Score_Changed_Count Global_Dedup_Final_Confidence_Changed_Count Global_Dedup_Rank_Changed_Count Global_Dedup_Top50_Changed Recommended_Selection_Strategy Recommended_Max_Matches Recommended_Dedup_Mode Physical_Peak_Dedup_Needed Formal_Change_Ready Evidence_For_Recommendation Remaining_Risk Required_Additional_Validation Detail_Original_Row_Count Detail_Written_Row_Count Detail_Truncated Detail_Truncation_Reason Top50_Shadow_Additional_Time_Seconds Dedup_Audit_Additional_Time_Seconds Audit_Peak_Tracked_Memory_MiB Audit_Mode Applied_To_Formal_Result""".split()
TOP_COLUMNS = """Tier_Top50_MS1_Support Tier_Top50_Support_Delta Tier_Top50_Shadow_Final_Score Tier_Top50_Shadow_Final_Confidence Tier_Top50_Shadow_Rank Tier_Top50_Rank_Changed Exact_Peak_Dedup_Support Charge_Dedup_Support Fragment_Peak_Dedup_Support Recommended_MS1_Dedup_Mode MS1_Top50_Dedup_Affected MS1_Top50_Dedup_Applied_To_Formal_Result""".split()
DIAGNOSTIC_COLUMNS = """MS1_Top50_Dedup_Audit_Available MS1_Tier_Top20_Filter_Passing MS1_Tier_Top50_Filter_Passing MS1_Tier_Unlimited_Filter_Passing MS1_Top50_Recovery_Fraction MS1_Exact_Duplicate_Row_Count MS1_Multi_Charge_Peak_Count MS1_Multi_Fragment_Shared_Peak_Count MS1_Exact_Dedup_Filter_Passing MS1_Charge_Dedup_Filter_Passing MS1_Fragment_Peak_Dedup_Filter_Passing MS1_Top50_Candidate_Set_Changed MS1_Top50_Rank_Changed MS1_Recommended_Dedup_Mode MS1_Top50_Formal_Change_Ready MS1_Top50_Dedup_Applied_To_Formal_Result""".split()


def _raw(value: Any) -> dict[str, Any]:
    return asdict(value) if is_dataclass(value) else dict(value)


def physical_peak_id(match: Any) -> str:
    return str(getattr(match, "_audit_physical_peak_id", "") or f"PK_{int(getattr(match, '_audit_peak_index', 0) or 0):06d}_{float(match.observed_mz):.8f}")


def match_key(match: Any) -> tuple[Any, ...]:
    return (str(match.fragment_id), physical_peak_id(match), int(match.charge), round(float(match.observed_mz), 10))


def candidate_key(candidate: Any) -> tuple[Any, ...]:
    row = _raw(candidate)
    return (str(row.get("source_type") or ""), str(row.get("source_id") or ""), int(row.get("charge") or 0), round(float(row.get("observed_mz") or 0), 10), str(row.get("modification_id") or ""))


def ranking_key(row: dict[str, Any]) -> tuple[str, str, str]:
    position = row.get("Candidate_tRNA_Position", row.get("Candidate_Positions_In_tRNA", ""))
    if isinstance(position, float) and position.is_integer():
        position = int(position)
    return (str(row.get("Modification_ID") or ""), str(row.get("Parent_Fragment_ID") or ""), str(position or ""))


def modification_family(candidate: Any, lookup: dict[str, Any]) -> str:
    row = _raw(candidate); mod_id = str(row.get("modification_id") or ""); mod = lookup.get(mod_id); raw = getattr(mod, "raw", {}) or {}
    return str(raw.get("chemical_group") or raw.get("near_isobaric_group") or getattr(mod, "category", "") or mod_id)


def deduplicate_matches(matches: list[Any], mode: str, mapping_config: dict[str, Any]) -> list[Any]:
    """Keep the deterministic best assignment for the requested exact-ID scope."""
    if mode == "no_dedup":
        return list(matches)
    if mode not in DEDUP_MODES:
        raise ValueError(f"Unknown MS1 dedup mode: {mode}")
    ordered = sorted(matches, key=lambda m: (tier_then_error_sort_key(m, mapping_config), str(m.fragment_id)))
    seen: set[Any] = set(); kept = []
    for match in ordered:
        if mode == "global_physical_peak_dedup":
            key = physical_peak_id(match)
        elif mode == "exact_physical_peak_dedup":
            key = (str(match.fragment_id), physical_peak_id(match))
        elif mode in {"physical_peak_charge_dedup", "fragment_then_physical_peak_dedup"}:
            key = (str(match.fragment_id), physical_peak_id(match))
        else:
            key = match_key(match)
        if key in seen:
            continue
        seen.add(key); kept.append(match)
    return kept


def select_tier_matches(matches: list[Any], limit: int | None, mapping_config: dict[str, Any]) -> list[Any]:
    ordered = sorted(matches, key=lambda m: tier_then_error_sort_key(m, mapping_config))
    return ordered if limit is None else ordered[:limit]


def build_near_mz_groups(matches: list[Any], tolerance_ppm: float) -> dict[str, str]:
    representatives: dict[str, Any] = {}
    for match in matches:
        representatives.setdefault(physical_peak_id(match), match)
    ordered = sorted(representatives.items(), key=lambda item: (float(item[1].observed_mz), item[0]))
    result: dict[str, str] = {}; group = 0; previous_mz = None
    for peak_id, match in ordered:
        mz = float(match.observed_mz)
        if previous_mz is None or abs(mz - previous_mz) / previous_mz * 1_000_000 > tolerance_ppm:
            group += 1
        result[peak_id] = f"NMZ_{group:06d}"; previous_mz = mz
    return result


def _change(base: list[dict[str, Any]], shadow: list[dict[str, Any]]) -> dict[str, Any]:
    bm = {ranking_key(row): row for row in base}; sm = {ranking_key(row): row for row in shadow}; common = set(bm) & set(sm)
    top = lambda rows: {key for key, row in rows.items() if int(row.get("Rank") or 999999) <= 50}
    cnm = lambda rows: {(key, rows[key].get("Final_Score"), rows[key].get("Final_Confidence"), rows[key].get("Rank")) for key in rows if key[0] == "cnm5U" and key[2] in {"36", "37", "38"}}
    return {"map": sm, "set": set(bm) != set(sm), "score": sum(bm[k].get("Final_Score") != sm[k].get("Final_Score") for k in common), "confidence": sum(bm[k].get("Final_Confidence") != sm[k].get("Final_Confidence") for k in common), "rank": sum(bm[k].get("Rank") != sm[k].get("Rank") for k in common), "top50": top(bm) != top(sm), "cnm": cnm(bm) != cnm(sm)}


def _quality(matches: list[Any]) -> dict[str, Any]:
    tiers = Counter(str(m.peak_tier or "").lower() for m in matches); confs = Counter(str(m.confidence or "").lower() for m in matches)
    errors = [abs(float(m.mass_error_ppm)) for m in matches]; intensities = [float(m.intensity or 0) for m in matches]
    return {"Major_Count": tiers["major"], "Minor_Count": tiers["minor"], "Trace_Count": tiers["trace"], "High_Confidence_Count": confs["high"], "Medium_Confidence_Count": confs["medium"], "Low_Confidence_Count": confs["low"], "Min_Abs_Error_PPM": min(errors) if errors else "", "Median_Abs_Error_PPM": median(errors) if errors else "", "Max_Intensity": max(intensities) if intensities else "", "Median_Intensity": median(intensities) if intensities else ""}


def _length_group(length: int) -> str:
    if length <= 2: return "1-2_nt"
    if length == 3: return "3_nt"
    if length == 4: return "4_nt"
    if length == 5: return "5_nt"
    if length <= 8: return "6-8_nt"
    return "9plus_nt"


def _review_map(ms2_results: dict[str, Any], ranking: list[dict[str, Any]], config: Any) -> dict[tuple[str, str, str], dict[str, Any]]:
    shadow_results = dict(ms2_results); shadow_results["Modification_Evidence_Ranking"] = ranking
    review = build_review_dashboard_results(shadow_results, config).get("Top_Modification_Candidates", [])
    source = review.to_dict("records") if isinstance(review, pd.DataFrame) else list(review or [])
    return {ranking_key(row): row for row in source}


def build_ms1_top50_dedup_audit(context, config, modifications, intact_results, baseline_matches, baseline_candidates, baseline_ranking, ms2_results, rule_set=None, pathways=None):
    audit_started = time.perf_counter()
    cfg = config.fragment_mapping or {}; items = list(context.get("fragments") or []); all_matches = [m for item in items for m in item.get("ranked_matches", [])]
    dataset_id = str((config.sequence or {}).get("name") or "dataset"); enzyme = str((config.digestion or {}).get("enzyme") or ""); limit = int(cfg.get("MS1_Tier_Top50_Limit", 50) or 50)
    by_fragment = {str(item["fragment"].fragment_id): list(item.get("ranked_matches") or []) for item in items}
    current = {fid: ranked[:int(context.get("configured_max_matches") or cfg.get("max_matches_per_fragment") or 20)] for fid, ranked in by_fragment.items()}
    tier20 = {fid: select_tier_matches(ranked, 20, cfg) for fid, ranked in by_fragment.items()}
    tier50 = {fid: select_tier_matches(ranked, limit, cfg) for fid, ranked in by_fragment.items()}
    unlimited = {fid: select_tier_matches(ranked, None, cfg) for fid, ranked in by_fragment.items()}
    sets: dict[str, dict[str, list[Any]]] = {"current_top20__no_dedup": current, "tier_top20__no_dedup": tier20, "tier_top50__no_dedup": tier50, "tier_unlimited__no_dedup": unlimited}
    for mode in DEDUP_MODES[1:4]:
        sets[f"tier_top50__{mode}"] = {fid: select_tier_matches(deduplicate_matches(ranked, mode, cfg), limit, cfg) for fid, ranked in by_fragment.items()}
    global_dedup = deduplicate_matches(all_matches, "global_physical_peak_dedup", cfg); global_by = defaultdict(list)
    for match in global_dedup: global_by[str(match.fragment_id)].append(match)
    sets["tier_top50__global_physical_peak_dedup"] = {fid: select_tier_matches(global_by.get(fid, []), limit, cfg) for fid in by_fragment}
    sets["tier_unlimited__fragment_then_physical_peak_dedup"] = {fid: deduplicate_matches(ranked, "fragment_then_physical_peak_dedup", cfg) for fid, ranked in by_fragment.items()}
    dedup_elapsed = time.perf_counter() - audit_started
    flat = {name: [m for group in groups.values() for m in group] for name, groups in sets.items()}
    filtered_matches = {name: [m for m in matches if passes_fragment_ms1_filter(m, cfg)] for name, matches in flat.items()}
    fragment_summaries = {
        name: {fid: {"Selected_Match_Count": len(groups[fid]), "Filtered_Match_Count": sum(passes_fragment_ms1_filter(m, cfg) for m in groups[fid]), "Unique_Physical_Peak_Count": len({physical_peak_id(m) for m in groups[fid]})} for fid in groups}
        for name, groups in sets.items()
    }

    candidates = {"current_top20__no_dedup": list(baseline_candidates)}; rankings = {"current_top20__no_dedup": list(baseline_ranking)}; reviews = {"current_top20__no_dedup": _review_map(ms2_results, baseline_ranking, config)}
    fragments = [item["fragment"] for item in items]
    cache: dict[tuple[tuple[Any, ...], ...], tuple[list[Any], list[dict[str, Any]], dict[Any, Any]]] = {}
    for name, matches in flat.items():
        if name == "current_top20__no_dedup": continue
        signature = tuple(sorted(match_key(m) for m in matches))
        if signature in cache:
            candidates[name], rankings[name], reviews[name] = cache[signature]; continue
        cs = search_known_modifications(matches, intact_results, modifications, config, warnings=None)
        rs, _ = build_modification_evidence_ranking(config, modifications, fragments, matches, cs, ms2_results, rule_set=rule_set, pathways=pathways)
        result = (cs, rs, _review_map(ms2_results, rs, config)); cache[signature] = result; candidates[name], rankings[name], reviews[name] = result
    changes = {name: _change(baseline_ranking, rankings[name]) for name in sets}
    base_candidate_keys = {candidate_key(c) for c in baseline_candidates}; lookup = {str(getattr(m, "id", "")): m for m in modifications}
    candidate_by_fragment = {}; families_by_fragment = {}
    for name, cs in candidates.items():
        candidate_by_fragment[name] = defaultdict(set); families_by_fragment[name] = defaultdict(set)
        for candidate in cs:
            row = _raw(candidate); fid = str(row.get("source_id") or "") if str(row.get("source_type") or "").lower() == "fragment" else ""
            if fid: candidate_by_fragment[name][fid].add(candidate_key(candidate)); families_by_fragment[name][fid].add(modification_family(candidate, lookup))

    pid_fragments = defaultdict(set); pid_charges = defaultdict(set); pid_fragment_charges = defaultdict(set); exact_groups = defaultdict(list); pid_forms = defaultdict(set)
    for match in all_matches:
        pid = physical_peak_id(match); fid = str(match.fragment_id); pid_fragments[pid].add(fid); pid_charges[pid].add(int(match.charge)); pid_fragment_charges[(pid, fid)].add(int(match.charge)); exact_groups[(pid, fid, int(match.charge))].append(match); pid_forms[pid].add((int(match.start), int(match.end), int(match.missed_cleavages), str(match.terminal_form)))
    near_groups = build_near_mz_groups(all_matches, float(cfg.get("mz_tolerance_ppm", 10) or 10)); near_counts = Counter(near_groups.values())
    near_duplicate_groups = {group for group, count in near_counts.items() if count > 1}
    mode_name = {"no_dedup":"No_Dedup", "exact_physical_peak_dedup":"Exact_Peak_Dedup", "physical_peak_charge_dedup":"Charge_Dedup", "fragment_then_physical_peak_dedup":"Fragment_Peak_Dedup", "global_physical_peak_dedup":"Global_Dedup"}
    top_rows = []
    combos = [("current_top20",20,"no_dedup"),("tier_top20",20,"no_dedup"),("tier_top50",limit,"no_dedup"),("tier_top50",limit,"exact_physical_peak_dedup"),("tier_top50",limit,"physical_peak_charge_dedup"),("tier_top50",limit,"fragment_then_physical_peak_dedup"),("tier_top50",limit,"global_physical_peak_dedup"),("tier_unlimited",None,"no_dedup"),("tier_unlimited",None,"fragment_then_physical_peak_dedup")]
    for strategy, selection_limit, mode in combos:
        name = f"{strategy}__{mode}"; groups = sets[name]
        all_candidate_keys = {candidate_key(c) for c in candidates[name]}
        for item in items:
            f = item["fragment"]; fid = str(f.fragment_id); pre = by_fragment[fid]; selected = groups[fid]; pids = [physical_peak_id(m) for m in selected]; candidate_keys = candidate_by_fragment[name][fid]; base_fragment_keys = candidate_by_fragment["current_top20__no_dedup"][fid]
            top_rows.append({"Dataset_ID":dataset_id,"Enzyme":enzyme,"Fragment_ID":fid,"Fragment_Start":f.start,"Fragment_End":f.end,"Fragment_Sequence":f.sequence,"Fragment_Length":len(f.sequence or ""),"Strategy":strategy,"Selection_Limit":"" if selection_limit is None else selection_limit,"Dedup_Mode":mode,"Pre_Truncation_Match_Count":len(pre),"Pre_Filter_Passing_Count":sum(passes_fragment_ms1_filter(m,cfg) for m in pre),"Selected_Match_Count":len(selected),"Selected_Filter_Passing_Count":sum(passes_fragment_ms1_filter(m,cfg) for m in selected),"Unique_Physical_Peak_Count":len(set(pids)),"Unique_Charge_Count":len({m.charge for m in selected}),"Duplicate_Row_Count":len(selected)-len({match_key(m) for m in selected}),"Duplicate_Physical_Peak_Count":len(selected)-len(set(pids)),"Multi_Charge_Physical_Peak_Count":sum(len(pid_fragment_charges[(pid,fid)])>1 for pid in set(pids)),"Multi_Fragment_Shared_Peak_Count":sum(len(pid_fragments[pid])>1 for pid in set(pids)),**_quality(selected),"Modification_Candidate_Count":len(candidate_keys),"Candidate_Key_Count":len(candidate_keys),"Added_Candidate_Count_Vs_Current":len(candidate_keys-base_fragment_keys),"Removed_Candidate_Count_Vs_Current":len(base_fragment_keys-candidate_keys),"Candidate_Support_Changed":candidate_keys!=base_fragment_keys,"Shadow_Final_Score_Changed":bool(changes[name]["score"]),"Shadow_Final_Confidence_Changed":bool(changes[name]["confidence"]),"Shadow_Rank_Changed":bool(changes[name]["rank"]),"Shadow_Top50_Changed":changes[name]["top50"],"cnm5U_Result_Changed":changes[name]["cnm"],"Recommended_Dedup_Mode":"fragment_then_physical_peak_dedup","Applied_To_Formal_Result":False})

    selected_keys = {name:{match_key(m) for m in matches} for name,matches in flat.items()}; top50_candidates = candidates["tier_top50__no_dedup"]; top50_cindex = defaultdict(list)
    for candidate in top50_candidates:
        row=_raw(candidate)
        if str(row.get("source_type") or "").lower()=="fragment": top50_cindex[(str(row.get("source_id") or ""),int(row.get("charge") or 0),round(float(row.get("observed_mz") or 0),10))].append(candidate)
    detail=[]
    for match in sorted(all_matches,key=lambda m:(str(m.fragment_id),tier_then_error_sort_key(m,cfg))):
        fid=str(match.fragment_id);pid=physical_peak_id(match);key=match_key(match);cs=top50_cindex.get((fid,int(match.charge),round(float(match.observed_mz),10)),[]);charge_count=len(pid_fragment_charges[(pid,fid)]);fragment_count=len(pid_fragments[pid]);exact_size=len(exact_groups[(pid,fid,int(match.charge))])
        kept_exact=key in selected_keys["tier_top50__exact_physical_peak_dedup"];kept_charge=key in selected_keys["tier_top50__physical_peak_charge_dedup"];kept_fragment=key in selected_keys["tier_top50__fragment_then_physical_peak_dedup"];kept_global=key in selected_keys["tier_top50__global_physical_peak_dedup"]
        reasons=[]
        if charge_count>1:reasons.append("same_fragment_multi_charge")
        if fragment_count>1:reasons.append("cross_fragment_shared")
        if exact_size>1:reasons.append("exact_duplicate")
        if not (kept_exact and kept_charge and kept_fragment):reasons.append("lower_rank_duplicate_assignment")
        if not kept_global and key in selected_keys["tier_top50__no_dedup"]:reasons.append("removed_by_global_cross_fragment_dedup")
        detail.append({"Dataset_ID":dataset_id,"Enzyme":enzyme,"Fragment_ID":fid,"Physical_Peak_ID":pid,"Spectrum_ID":match.scan_id or "","Scan_Index":"","Peak_Index":getattr(match,"_audit_peak_index",""),"RT":match.rt,"Observed_MZ":match.observed_mz,"Intensity":match.intensity,"Charge":match.charge,"Observed_Neutral_Mass":match.fragment_mass+match.mass_error_da*abs(match.charge),"Theoretical_Neutral_Mass":match.fragment_mass,"Mass_Error_Da":match.mass_error_da,"Mass_Error_PPM":match.mass_error_ppm,"Abs_Mass_Error_PPM":abs(match.mass_error_ppm),"Peak_Tier":match.peak_tier,"Confidence":match.confidence,"Passes_Fragment_MS1_Filter":passes_fragment_ms1_filter(match,cfg),"Tier_Top20_Selected":key in selected_keys["tier_top20__no_dedup"],"Tier_Top50_Selected":key in selected_keys["tier_top50__no_dedup"],"Tier_Unlimited_Selected":key in selected_keys["tier_unlimited__no_dedup"],"Exact_Physical_Duplicate_Group_Size":exact_size,"Charge_Assignment_Count_For_Peak":len(pid_charges[pid]),"Fragment_Assignment_Count_For_Peak":fragment_count,"Theoretical_Form_Count_For_Peak":len(pid_forms[pid]),"Same_Fragment_Multi_Charge":charge_count>1,"Cross_Fragment_Shared":fragment_count>1,"Near_MZ_Group_ID":near_groups[pid],"Kept_No_Dedup":key in selected_keys["tier_top50__no_dedup"],"Kept_Exact_Peak_Dedup":kept_exact,"Kept_Charge_Dedup":kept_charge,"Kept_Fragment_Peak_Dedup":kept_fragment,"Kept_Global_Dedup":kept_global,"Dedup_Selection_Reason":";".join(reasons) or "unique_assignment","Possible_Modification_Count":len(cs),"Candidate_Key_List":";".join(sorted("|".join(map(str,candidate_key(c))) for c in cs)),"Candidate_Family_List":";".join(sorted({modification_family(c,lookup) for c in cs}))})

    def count_passing(name): return sum(passes_fragment_ms1_filter(m,cfg) for m in flat[name])
    current_passing=count_passing("current_top20__no_dedup");top20_passing=count_passing("tier_top20__no_dedup");top50_passing=count_passing("tier_top50__no_dedup");unlimited_passing=count_passing("tier_unlimited__no_dedup")
    top20_keys=selected_keys["tier_top20__no_dedup"];added50=[m for m in flat["tier_top50__no_dedup"] if match_key(m) not in top20_keys and passes_fragment_ms1_filter(m,cfg)]
    exact_duplicate_rows=sum(max(0,len(group)-1) for group in exact_groups.values());multi_charge=sum(len(charges)>1 for charges in pid_charges.values());local_multi_charge=sum(len(charges)>1 for charges in pid_fragment_charges.values());multi_fragment=sum(len(fragments)>1 for fragments in pid_fragments.values())
    no_name="tier_top50__no_dedup"; exact_name="tier_top50__exact_physical_peak_dedup";charge_name="tier_top50__physical_peak_charge_dedup";fragment_name="tier_top50__fragment_then_physical_peak_dedup";global_name="tier_top50__global_physical_peak_dedup"
    summary_base={"Dataset_ID":dataset_id,"Enzyme":enzyme,"Sequence_Name":dataset_id,"Total_Matches":len(all_matches),"Filter_Passing_Matches":sum(passes_fragment_ms1_filter(m,cfg) for m in all_matches),"Unique_Physical_Peaks":len(pid_fragments),"Exact_Duplicate_Rows":exact_duplicate_rows,"Multi_Charge_Peaks":multi_charge,"Multi_Fragment_Shared_Peaks":multi_fragment,"Near_MZ_Groups":len(near_duplicate_groups),"Short_Fragment_Match_Count":sum(len(m.sequence or "")<=4 for m in all_matches),"Short_Fragment_Duplicate_Fraction":1-len({(str(m.fragment_id),physical_peak_id(m)) for m in all_matches if len(m.sequence or "")<=4})/max(1,sum(len(m.sequence or "")<=4 for m in all_matches)),"Average_Charge_Assignments_Per_Peak":sum(len(v) for v in pid_charges.values())/max(1,len(pid_charges)),"Average_Fragment_Assignments_Per_Peak":sum(len(v) for v in pid_fragments.values())/max(1,len(pid_fragments)),"Tier_Top20_Filter_Passing":top20_passing,"Tier_Top50_Filter_Passing":top50_passing,"Tier_Unlimited_Filter_Passing":unlimited_passing,"Top50_Recovery_Fraction":1.0 if unlimited_passing<=current_passing else (top50_passing-current_passing)/(unlimited_passing-current_passing),"Top50_Added_Major":sum(str(m.peak_tier).lower()=="major" for m in added50),"Top50_Added_Minor":sum(str(m.peak_tier).lower()=="minor" for m in added50),"Top50_Added_High":sum(str(m.confidence).lower()=="high" for m in added50),"Top50_Added_Medium":sum(str(m.confidence).lower()=="medium" for m in added50),"No_Dedup_Selected_Count":len(flat[no_name]),"Exact_Peak_Dedup_Selected_Count":len(flat[exact_name]),"Charge_Dedup_Selected_Count":len(flat[charge_name]),"Fragment_Peak_Dedup_Selected_Count":len(flat[fragment_name]),"Global_Dedup_Selected_Count":len(flat[global_name]),"No_Dedup_Filter_Passing":count_passing(no_name),"Exact_Peak_Dedup_Filter_Passing":count_passing(exact_name),"Charge_Dedup_Filter_Passing":count_passing(charge_name),"Fragment_Peak_Dedup_Filter_Passing":count_passing(fragment_name),"Global_Dedup_Filter_Passing":count_passing(global_name),"Top50_Modification_Candidate_Count":len(top50_candidates),"Top50_Candidate_Key_Set_Changed":{candidate_key(c) for c in top50_candidates}!=base_candidate_keys,"Top50_Final_Score_Changed_Count":changes[no_name]["score"],"Top50_Final_Confidence_Changed_Count":changes[no_name]["confidence"],"Top50_Rank_Changed_Count":changes[no_name]["rank"],"Top50_Membership_Changed":changes[no_name]["top50"],"Recommended_Selection_Strategy":"tier_then_error","Recommended_Max_Matches":limit,"Recommended_Dedup_Mode":"fragment_then_physical_peak_dedup" if exact_duplicate_rows or local_multi_charge else "no_dedup","Physical_Peak_Dedup_Needed":bool(exact_duplicate_rows or local_multi_charge),"Formal_Change_Ready":False,"Evidence_For_Recommendation":f"tier top50 recovers {top50_passing} of {unlimited_passing} filter-passing rows; exact duplicate rows={exact_duplicate_rows}, same-fragment multi-charge groups={local_multi_charge}, cross-fragment shared peaks={multi_fragment}","Remaining_Risk":"candidate-positive data and the scientific meaning of cross-fragment peak sharing remain unresolved","Required_Additional_Validation":"candidate-positive RNase dataset; replicate performance; review shared-peak fragment ambiguity","Audit_Mode":"shadow_tier_top50_physical_peak_dedup","Applied_To_Formal_Result":False}
    for name,prefix in ((exact_name,"Exact_Dedup"),(charge_name,"Charge_Dedup"),(fragment_name,"Fragment_Dedup"),(global_name,"Global_Dedup")):
        summary_base.update({f"{prefix}_Candidate_Key_Set_Changed":{candidate_key(c) for c in candidates[name]}!=base_candidate_keys,f"{prefix}_Final_Score_Changed_Count":changes[name]["score"],f"{prefix}_Final_Confidence_Changed_Count":changes[name]["confidence"],f"{prefix}_Rank_Changed_Count":changes[name]["rank"],f"{prefix}_Top50_Changed":changes[name]["top50"]})
    summary_rows=[]
    summary_rows.append(dict(summary_base,Summary_Scope="dataset",Fragment_Length_Group="all",Fragment_Count=len(items)))
    for group in ("1-2_nt","3_nt","4_nt","5_nt","6-8_nt","9plus_nt"):
        fids={str(item["fragment"].fragment_id) for item in items if _length_group(len(item["fragment"].sequence or ""))==group};raw=[m for m in all_matches if str(m.fragment_id) in fids];pids={physical_peak_id(m) for m in raw};fragment_peaks={(str(m.fragment_id),physical_peak_id(m)) for m in raw};row={column:"" for column in SUMMARY_COLUMNS};row.update({"Dataset_ID":dataset_id,"Enzyme":enzyme,"Sequence_Name":dataset_id,"Summary_Scope":"fragment_length","Fragment_Length_Group":group,"Fragment_Count":len(fids),"Total_Matches":len(raw),"Filter_Passing_Matches":sum(passes_fragment_ms1_filter(m,cfg) for m in raw),"Unique_Physical_Peaks":len(pids),"Exact_Duplicate_Rows":sum(max(0,len(v)-1) for k,v in exact_groups.items() if k[1] in fids),"Multi_Charge_Peaks":sum(len(v)>1 for (pid,fid),v in pid_fragment_charges.items() if fid in fids),"Multi_Fragment_Shared_Peaks":sum(len(pid_fragments[pid])>1 for pid in pids),"Near_MZ_Groups":len({near_groups[pid] for pid in pids if near_groups[pid] in near_duplicate_groups}),"Short_Fragment_Match_Count":len(raw) if group in {"1-2_nt","3_nt","4_nt"} else 0,"Short_Fragment_Duplicate_Fraction":1-len(fragment_peaks)/max(1,len(raw)),"Average_Charge_Assignments_Per_Peak":sum(len(pid_charges[pid]) for pid in pids)/max(1,len(pids)),"Average_Fragment_Assignments_Per_Peak":sum(len(pid_fragments[pid]) for pid in pids)/max(1,len(pids)),"Tier_Top20_Filter_Passing":sum(passes_fragment_ms1_filter(m,cfg) for fid in fids for m in tier20[fid]),"Tier_Top50_Filter_Passing":sum(passes_fragment_ms1_filter(m,cfg) for fid in fids for m in tier50[fid]),"Tier_Unlimited_Filter_Passing":sum(passes_fragment_ms1_filter(m,cfg) for fid in fids for m in unlimited[fid]),"Fragment_Peak_Dedup_Filter_Passing":sum(passes_fragment_ms1_filter(m,cfg) for fid in fids for m in sets[fragment_name][fid]),"Applied_To_Formal_Result":False});den=row["Tier_Unlimited_Filter_Passing"];row["Top50_Recovery_Fraction"]=row["Tier_Top50_Filter_Passing"]/den if den else 1.0;summary_rows.append(row)
    original=len(detail);max_rows=int((config.reporting or {}).get("max_excel_rows_per_sheet",100000) or 100000);detail=detail[:max_rows];summary_base.update({"Detail_Original_Row_Count":original,"Detail_Written_Row_Count":len(detail),"Detail_Truncated":len(detail)<original,"Detail_Truncation_Reason":"max_excel_rows_per_sheet; deterministic fragment/tier order" if len(detail)<original else ""})
    for row in top_rows:
        row["Recommended_Dedup_Mode"] = summary_base["Recommended_Dedup_Mode"]
    for key in ("Detail_Original_Row_Count","Detail_Written_Row_Count","Detail_Truncated","Detail_Truncation_Reason"):summary_rows[0][key]=summary_base[key]
    total_elapsed = time.perf_counter() - audit_started
    summary_base["Top50_Shadow_Additional_Time_Seconds"] = total_elapsed
    summary_base["Dedup_Audit_Additional_Time_Seconds"] = dedup_elapsed
    summary_rows[0]["Top50_Shadow_Additional_Time_Seconds"] = total_elapsed
    summary_rows[0]["Dedup_Audit_Additional_Time_Seconds"] = dedup_elapsed
    return {"mapping_config":cfg,"base_ranking_map":{ranking_key(row):row for row in baseline_ranking},"top50_rows":[{c:r.get(c,"") for c in TOP50_COLUMNS} for r in top_rows],"detail_rows":[{c:r.get(c,"") for c in DETAIL_COLUMNS} for r in detail],"summary_rows":[{c:r.get(c,"") for c in SUMMARY_COLUMNS} for r in summary_rows],"summary":summary_base,"sets":sets,"flat_matches":flat,"filtered_matches":filtered_matches,"fragment_summaries":fragment_summaries,"candidates":candidates,"candidate_summaries":{name:summarize_known_modification_candidates(value) for name,value in candidates.items()},"candidate_families":{name:sorted({modification_family(c,lookup) for c in value}) for name,value in candidates.items()},"rankings":rankings,"reviews":reviews,"changes":changes,"fragment_internal":by_fragment}


def append_top50_shadow_columns(rows, audit):
    is_frame=isinstance(rows,pd.DataFrame);source=rows.to_dict("records") if is_frame else list(rows or []);columns=list(rows.columns) if is_frame else (list(source[0]) if source else []);out=[]
    summary=audit["summary"];names={"top50":"tier_top50__no_dedup","exact":"tier_top50__exact_physical_peak_dedup","charge":"tier_top50__physical_peak_charge_dedup","fragment":"tier_top50__fragment_then_physical_peak_dedup"}
    support={name:Counter(str(m.fragment_id) for m in audit["flat_matches"][key] if passes_fragment_ms1_filter(m,audit["mapping_config"])) for name,key in names.items()};base=Counter(str(m.fragment_id) for m in audit["flat_matches"]["current_top20__no_dedup"] if passes_fragment_ms1_filter(m,audit["mapping_config"]));rankmap=audit["changes"][names["top50"]]["map"]
    for original in source:
        row=dict(original);fid=str(row.get("Parent_Fragment_ID") or "");key=ranking_key(row);shadow=rankmap.get(key,{});vals={"Tier_Top50_MS1_Support":support["top50"][fid],"Tier_Top50_Support_Delta":support["top50"][fid]-base[fid],"Tier_Top50_Shadow_Final_Score":shadow.get("Final_Score",""),"Tier_Top50_Shadow_Final_Confidence":shadow.get("Final_Confidence",""),"Tier_Top50_Shadow_Rank":shadow.get("Rank",""),"Tier_Top50_Rank_Changed":bool(shadow and shadow.get("Rank")!=audit["base_ranking_map"].get(key,{}).get("Rank")),"Exact_Peak_Dedup_Support":support["exact"][fid],"Charge_Dedup_Support":support["charge"][fid],"Fragment_Peak_Dedup_Support":support["fragment"][fid],"Recommended_MS1_Dedup_Mode":summary["Recommended_Dedup_Mode"],"MS1_Top50_Dedup_Affected":support["top50"][fid]!=base[fid] or support["fragment"][fid]!=support["top50"][fid],"MS1_Top50_Dedup_Applied_To_Formal_Result":False};row.update(vals);out.append(row)
    return pd.DataFrame(out,columns=columns+TOP_COLUMNS) if is_frame else out


def append_top50_diagnostic_columns(rows,audit):
    s=audit["summary"];vals={"MS1_Top50_Dedup_Audit_Available":True,"MS1_Tier_Top20_Filter_Passing":s["Tier_Top20_Filter_Passing"],"MS1_Tier_Top50_Filter_Passing":s["Tier_Top50_Filter_Passing"],"MS1_Tier_Unlimited_Filter_Passing":s["Tier_Unlimited_Filter_Passing"],"MS1_Top50_Recovery_Fraction":s["Top50_Recovery_Fraction"],"MS1_Exact_Duplicate_Row_Count":s["Exact_Duplicate_Rows"],"MS1_Multi_Charge_Peak_Count":s["Multi_Charge_Peaks"],"MS1_Multi_Fragment_Shared_Peak_Count":s["Multi_Fragment_Shared_Peaks"],"MS1_Exact_Dedup_Filter_Passing":s["Exact_Peak_Dedup_Filter_Passing"],"MS1_Charge_Dedup_Filter_Passing":s["Charge_Dedup_Filter_Passing"],"MS1_Fragment_Peak_Dedup_Filter_Passing":s["Fragment_Peak_Dedup_Filter_Passing"],"MS1_Top50_Candidate_Set_Changed":s["Top50_Candidate_Key_Set_Changed"],"MS1_Top50_Rank_Changed":bool(s["Top50_Rank_Changed_Count"]),"MS1_Recommended_Dedup_Mode":s["Recommended_Dedup_Mode"],"MS1_Top50_Formal_Change_Ready":s["Formal_Change_Ready"],"MS1_Top50_Dedup_Applied_To_Formal_Result":False};is_frame=isinstance(rows,pd.DataFrame);source=rows.to_dict("records") if is_frame else list(rows or [{}]);out=[dict(r,**vals) for r in source]
    return pd.DataFrame(out,columns=list(rows.columns)+DIAGNOSTIC_COLUMNS) if is_frame else out
