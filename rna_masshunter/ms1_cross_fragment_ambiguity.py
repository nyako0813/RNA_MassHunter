"""Cross-fragment physical-peak assignment ambiguity shadow audit.

This module never mutates formal matches, candidates, scores, or ranks.  Assignment
weights describe alternative interpretations of tier-top50 MS1 evidence only.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, is_dataclass
from statistics import median
import re
import time
from typing import Any

import pandas as pd

from rna_masshunter.evidence_ranking import build_modification_evidence_ranking
from rna_masshunter.modification_search import search_known_modifications
from rna_masshunter.ms1_selection_strategy_audit import passes_fragment_ms1_filter
from rna_masshunter.ms1_top50_dedup_audit import (
    _change, candidate_key, modification_family, physical_peak_id, ranking_key,
    select_tier_matches,
)

AMBIGUITY_COLUMNS = """Dataset_ID Enzyme Ambiguity_Group_ID Physical_Peak_ID Spectrum_ID Scan_Index Peak_Index RT Observed_MZ Intensity Charge Assignment_Count Filter_Passing_Assignment_Count Unique_Fragment_Count Unique_Fragment_Sequence_Count Unique_Position_Count Unique_Length_Count Unique_Theoretical_Mass_Count Modification_Eligible_Assignment_Count Candidate_Supporting_Assignment_Count Same_Sequence_Multiple_Position Overlapping_Fragment_Count Nonoverlapping_Fragment_Count Multiple_Lengths Multiple_Missed_Cleavages Multiple_Terminal_Forms Near_Isobaric_Assignments Max_Theoretical_Mass_Separation_Da Max_Theoretical_Mass_Separation_PPM Best_Fragment_ID Best_Fragment_Start Best_Fragment_End Best_Fragment_Sequence Best_Fragment_Length Best_Peak_Tier Best_Confidence Best_Abs_Error_PPM Second_Best_Fragment_ID Second_Best_Abs_Error_PPM Best_Second_Error_Delta_PPM Best_Assignment_Dominance Full_Count_Total_Weight Winner_Take_All_Total_Weight Equal_Fraction_Total_Weight Quality_Weighted_Total_Weight Candidate_Key_Count Candidate_Key_List Candidate_Family_List Position_Discriminating cnm5U_Relevant Ambiguity_Severity Recommended_Handling Applied_To_Formal_Result""".split()
DETAIL_COLUMNS = """Dataset_ID Enzyme Ambiguity_Group_ID Physical_Peak_ID Fragment_ID Fragment_Start Fragment_End Fragment_Sequence Fragment_Length Missed_Cleavage_Count Terminal_Form Observed_MZ Charge Intensity Peak_Tier Confidence Theoretical_Neutral_Mass Observed_Neutral_Mass Mass_Error_Da Mass_Error_PPM Abs_Mass_Error_PPM Passes_Fragment_MS1_Filter Assignment_Rank Is_Best_Assignment Is_Second_Best_Assignment Same_Sequence_As_Best Overlaps_Best_Fragment Position_Distance_From_Best Theoretical_Mass_Delta_From_Best Possible_Modification_Count Possible_Modification_IDs Candidate_Key_List Candidate_Family_List Full_Count_Weight Winner_Take_All_Weight Equal_Fraction_Weight Quality_Weighted_Fraction Fragment_Family_Fraction Ambiguity_Severity Assignment_Reason""".split()
SUMMARY_COLUMNS = """Dataset_ID Enzyme Sequence_Name Summary_Scope Fragment_Length_Group Total_Physical_Peaks Shared_Physical_Peaks Shared_Peak_Fraction Total_Assignments_In_Shared_Peaks Passing_Shared_Peaks Multiple_Passing_Assignment_Peaks Same_Sequence_Multi_Position_Peaks Overlapping_Fragment_Peaks Nonoverlapping_Fragment_Peaks Multi_Length_Peaks Multi_Missed_Cleavage_Peaks Multi_Terminal_Form_Peaks Modification_Eligible_Shared_Peaks Candidate_Supporting_Shared_Peaks Physical_Peak_Count Assignment_Count Passing_Assignment_Count Average_Assignments_Per_Peak Best_Second_PPM_Delta_Median Position_Localization_Contribution Low_Count Medium_Count High_Count Critical_Count Full_Count_Effective_Support Winner_Take_All_Effective_Support Equal_Fraction_Effective_Support Quality_Weighted_Effective_Support Fragment_Family_Effective_Support Equal_Fraction_Support_Reduction Winner_Take_All_Candidate_Set_Changed Equal_Fraction_Candidate_Set_Changed Quality_Weighted_Candidate_Set_Changed Fragment_Family_Candidate_Set_Changed Winner_Take_All_Score_Changed_Count Equal_Fraction_Score_Changed_Count Quality_Weighted_Score_Changed_Count Fragment_Family_Score_Changed_Count Winner_Take_All_Confidence_Changed_Count Equal_Fraction_Confidence_Changed_Count Quality_Weighted_Confidence_Changed_Count Fragment_Family_Confidence_Changed_Count Winner_Take_All_Rank_Changed_Count Equal_Fraction_Rank_Changed_Count Quality_Weighted_Rank_Changed_Count Fragment_Family_Rank_Changed_Count Winner_Take_All_Top50_Changed Equal_Fraction_Top50_Changed Quality_Weighted_Top50_Changed Fragment_Family_Top50_Changed Recommended_Ambiguity_Handling Shared_Peak_Double_Counting_Risk Position_Localization_Risk Candidate_Ranking_Risk Formal_Change_Ready Evidence_For_Recommendation Remaining_Risk Required_Additional_Validation Detail_Original_Row_Count Detail_Written_Row_Count Detail_Truncated Detail_Truncation_Reason Ambiguity_Grouping_Time_Seconds Shadow_Weighting_Time_Seconds Audit_Peak_Tracked_Memory_MiB Audit_Mode Applied_To_Formal_Result""".split()
TOP_COLUMNS = """CrossFrag_Shared_Peak_Support CrossFrag_Unique_Peak_Support CrossFrag_Ambiguous_Support_Fraction Winner_Take_All_Shadow_Support Equal_Fraction_Shadow_Support Quality_Weighted_Shadow_Support Winner_Take_All_Shadow_Final_Score Equal_Fraction_Shadow_Final_Score Quality_Weighted_Shadow_Final_Score Winner_Take_All_Shadow_Rank Equal_Fraction_Shadow_Rank Quality_Weighted_Shadow_Rank CrossFrag_Ambiguity_Affected Recommended_CrossFrag_Handling CrossFrag_Ambiguity_Applied_To_Formal_Result""".split()
DIAGNOSTIC_COLUMNS = """MS1_CrossFrag_Audit_Available MS1_Shared_Physical_Peak_Count MS1_Shared_Physical_Peak_Fraction MS1_Multiple_Passing_Assignment_Peaks MS1_High_Ambiguity_Group_Count MS1_Critical_Ambiguity_Group_Count MS1_Candidate_Supporting_Shared_Peaks MS1_Winner_Take_All_Candidate_Set_Changed MS1_Equal_Fraction_Candidate_Set_Changed MS1_Quality_Weighted_Candidate_Set_Changed MS1_Winner_Take_All_Rank_Changed MS1_Equal_Fraction_Rank_Changed MS1_Quality_Weighted_Rank_Changed MS1_Recommended_CrossFrag_Handling MS1_CrossFrag_Formal_Change_Ready MS1_CrossFrag_Applied_To_Formal_Result""".split()
STRATEGIES = ("full_count", "winner_take_all", "equal_fraction", "quality_weighted_fraction", "ambiguity_flag_only", "fragment_family_fraction")


def _raw(value: Any) -> dict[str, Any]:
    return asdict(value) if is_dataclass(value) else dict(value)


def _length_group(length: int) -> str:
    if length <= 2: return "1-2_nt"
    if length == 3: return "3_nt"
    if length == 4: return "4_nt"
    if length == 5: return "5_nt"
    if length <= 8: return "6-8_nt"
    return "9plus_nt"


def _tier(value: Any) -> int:
    return {"major": 3, "minor": 2, "trace": 1}.get(str(value or "").lower(), 0)


def _confidence(value: Any) -> int:
    return {"high": 3, "medium": 2, "low": 1}.get(str(value or "").lower(), 0)


def _overlap(a: Any, b: Any) -> bool:
    return int(a.start) <= int(b.end) and int(b.start) <= int(a.end)


def _scan_index(scan_id: Any) -> Any:
    found = re.search(r"(?:scan=|index=)(\d+)", str(scan_id or ""), re.I)
    return int(found.group(1)) if found else ""


def assignment_sort_key(match: Any, mapping_config: dict[str, Any]) -> tuple[Any, ...]:
    """Deterministic best: pass, tier, confidence, ppm, intensity, length, MC, ID/order."""
    return (
        -int(passes_fragment_ms1_filter(match, mapping_config)),
        -_tier(match.peak_tier), -_confidence(match.confidence),
        abs(float(match.mass_error_ppm)), -float(match.intensity or 0),
        -len(match.sequence or ""), int(match.missed_cleavages or 0),
        str(match.fragment_id), int(getattr(match, "_audit_generation_order", 0) or 0),
    )


def _quality_value(match: Any, tolerance: float) -> float:
    ppm_part = max(0.0, 1.0 - abs(float(match.mass_error_ppm)) / max(tolerance, 1e-12))
    return float(_tier(match.peak_tier) + _confidence(match.confidence)) + ppm_part


def assignment_weights(matches: list[Any], mapping_config: dict[str, Any]) -> dict[str, dict[int, float]]:
    """Return weights keyed by object identity; failing assignments have zero support."""
    passing = [m for m in matches if passes_fragment_ms1_filter(m, mapping_config)]
    result = {name: {id(m): 0.0 for m in matches} for name in STRATEGIES}
    if not passing:
        return result
    ordered = sorted(passing, key=lambda m: assignment_sort_key(m, mapping_config))
    for m in passing:
        result["full_count"][id(m)] = 1.0
        result["ambiguity_flag_only"][id(m)] = 1.0
        result["equal_fraction"][id(m)] = 1.0 / len(passing)
    result["winner_take_all"][id(ordered[0])] = 1.0
    tolerance = float(mapping_config.get("mz_tolerance_ppm", 10) or 10)
    q = {_id: _quality_value(m, tolerance) for m in passing for _id in [id(m)]}
    qsum = sum(q.values())
    for m in passing:
        result["quality_weighted_fraction"][id(m)] = q[id(m)] / qsum if qsum else 1.0 / len(passing)
    # Connected components under same-sequence or coordinate overlap form a family.
    remaining = set(map(id, passing)); by_id = {id(m): m for m in passing}
    families: list[list[Any]] = []
    while remaining:
        seed = min(remaining, key=lambda key: assignment_sort_key(by_id[key], mapping_config))
        component = {seed}; frontier = [seed]; remaining.remove(seed)
        while frontier:
            current = by_id[frontier.pop()]
            related = [key for key in remaining if by_id[key].sequence == current.sequence or _overlap(by_id[key], current)]
            for key in related:
                remaining.remove(key); component.add(key); frontier.append(key)
        families.append([by_id[key] for key in component])
    for family in families:
        for m in family:
            result["fragment_family_fraction"][id(m)] = 1.0 / len(family)
    return result


def _candidate_indexes(candidates: list[Any], modifications: list[Any]) -> tuple[dict[tuple[str, int, float], list[Any]], dict[str, Any]]:
    index: dict[tuple[str, int, float], list[Any]] = defaultdict(list)
    for candidate in candidates or []:
        row = _raw(candidate)
        if str(row.get("source_type") or "").lower() == "fragment":
            index[(str(row.get("source_id") or ""), int(row.get("charge") or 0), round(float(row.get("observed_mz") or 0), 10))].append(candidate)
    return index, {str(getattr(m, "id", "")): m for m in modifications or []}


def _match_candidates(match: Any, index: dict[tuple[str, int, float], list[Any]]) -> list[Any]:
    return index.get((str(match.fragment_id), int(match.charge), round(float(match.observed_mz), 10)), [])


def _severity(matches: list[Any], passing: list[Any], candidate_count: int, same_multi: bool, overlap_count: int, nonoverlap_count: int, delta: float) -> str:
    if candidate_count and len(passing) > 1 and (same_multi or nonoverlap_count or delta <= 0.1): return "Critical"
    if len(passing) > 1 and (same_multi or overlap_count or nonoverlap_count or delta <= 1.0): return "High"
    if len(passing) > 1 or candidate_count or nonoverlap_count: return "Medium"
    return "Low"


def _dominance(best: Any, second: Any | None, mapping_config: dict[str, Any]) -> str:
    if second is None: return "single_passing_assignment"
    bp = passes_fragment_ms1_filter(best, mapping_config); sp = passes_fragment_ms1_filter(second, mapping_config)
    if bp and not sp: return "filter_pass_winner"
    if _tier(best.peak_tier) != _tier(second.peak_tier): return "tier_winner"
    if _confidence(best.confidence) != _confidence(second.confidence): return "confidence_winner"
    delta = abs(abs(float(second.mass_error_ppm)) - abs(float(best.mass_error_ppm)))
    if delta <= 1e-9: return "equal_ppm_tie_break"
    tolerance = float(mapping_config.get("mz_tolerance_ppm", 10) or 10)
    return "clear_ppm_winner" if delta >= tolerance / 3 else "weak_ppm_winner"


def _candidate_change_rows(base_ranking: list[dict[str, Any]], winner_ranking: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    winner = _change(base_ranking, winner_ranking)
    unchanged = _change(base_ranking, base_ranking)
    return {"winner_take_all": winner, "equal_fraction": unchanged, "quality_weighted": unchanged, "fragment_family": unchanged}


def build_ms1_cross_fragment_ambiguity_audit(context, config, modifications, intact_results, baseline_candidates, baseline_ranking, ms2_results, rule_set=None, pathways=None, top50_audit=None):
    started = time.perf_counter(); cfg = config.fragment_mapping or {}
    items = list(context.get("fragments") or [])
    all_matches = [m for item in items for m in item.get("ranked_matches", [])]
    by_fragment = {str(item["fragment"].fragment_id): list(item.get("ranked_matches") or []) for item in items}
    limit = int(cfg.get("MS1_Tier_Top50_Limit", 50) or 50)
    top50 = [m for fid in sorted(by_fragment) for m in select_tier_matches(by_fragment[fid], limit, cfg)]
    top50_keys = {id(m) for m in top50}
    dataset = str((config.sequence or {}).get("name") or "dataset"); enzyme = str((config.digestion or {}).get("enzyme") or "")
    if top50_audit:
        strategy_name = "tier_top50__no_dedup"
        top_candidates = list(top50_audit["candidates"].get(strategy_name, []))
        top_ranking = list(top50_audit["rankings"].get(strategy_name, []))
    else:
        top_candidates = search_known_modifications(top50, intact_results, modifications, config, warnings=None)
        top_ranking, _ = build_modification_evidence_ranking(config, modifications, [item["fragment"] for item in items], top50, top_candidates, ms2_results, rule_set=rule_set, pathways=pathways)
    candidate_index, modification_lookup = _candidate_indexes(top_candidates, modifications)
    groups: dict[str, list[Any]] = defaultdict(list)
    for match in all_matches: groups[physical_peak_id(match)].append(match)
    shared = [(pid, matches) for pid, matches in groups.items() if len({str(m.fragment_id) for m in matches}) > 1]
    shared.sort(key=lambda item: item[0])
    grouping_elapsed = time.perf_counter() - started

    ambiguity_rows: list[dict[str, Any]] = []; details: list[dict[str, Any]] = []
    weights_started = time.perf_counter(); strategy_support = Counter(); fragment_support = {name: Counter() for name in STRATEGIES}
    group_meta: dict[str, dict[str, Any]] = {}
    for group_number, (pid, matches) in enumerate(shared, 1):
        ordered = sorted(matches, key=lambda m: assignment_sort_key(m, cfg)); best = ordered[0]; second = ordered[1] if len(ordered) > 1 else None
        selected = [m for m in matches if id(m) in top50_keys]
        weights = assignment_weights(selected, cfg)
        passing = [m for m in selected if passes_fragment_ms1_filter(m, cfg)]
        sequences = {str(m.sequence or "") for m in matches}; positions = {(int(m.start), int(m.end)) for m in matches}
        same_multi = any(len({(m.start, m.end) for m in matches if m.sequence == seq}) > 1 for seq in sequences)
        pairs = [(a, b) for i, a in enumerate(matches) for b in matches[i + 1:]]
        overlaps = sum(_overlap(a, b) for a, b in pairs); nonoverlaps = len(pairs) - overlaps
        masses = [float(m.fragment_mass) for m in matches]; mass_sep = max(masses) - min(masses) if masses else 0.0
        mass_ppm = mass_sep / min(masses) * 1e6 if masses and min(masses) else 0.0
        near = mass_ppm <= 2 * float(cfg.get("mz_tolerance_ppm", 10) or 10)
        candidate_assignments = {id(m): _match_candidates(m, candidate_index) for m in matches}
        candidates = [c for values in candidate_assignments.values() for c in values]
        candidate_keys = sorted({"|".join(map(str, candidate_key(c))) for c in candidates})
        families = sorted({modification_family(c, modification_lookup) for c in candidates})
        delta = abs(abs(float(second.mass_error_ppm)) - abs(float(best.mass_error_ppm))) if second else 0.0
        severity = _severity(matches, passing, len(candidates), same_multi, overlaps, nonoverlaps, delta)
        gid = f"CFA_{group_number:06d}"
        total_weights = {name: sum(weights[name].values()) for name in STRATEGIES}
        for name, total in total_weights.items(): strategy_support[name] += total
        for m in selected:
            for name in STRATEGIES: fragment_support[name][str(m.fragment_id)] += weights[name].get(id(m), 0.0)
        position_discriminating = len(positions) > 1 and bool(passing)
        cnm = any(str(_raw(c).get("modification_id") or "") == "cnm5U" for c in candidates)
        same_seq_best = lambda m: m.sequence == best.sequence
        group_meta[pid] = {"severity": severity, "passing": len(passing), "same": same_multi, "overlap": overlaps, "nonoverlap": nonoverlaps, "candidate": bool(candidates), "lengths": {_length_group(len(m.sequence or "")) for m in matches}, "delta": delta, "position": position_discriminating, "weights": total_weights, "assignments": len(matches)}
        ambiguity_rows.append({
            "Dataset_ID":dataset,"Enzyme":enzyme,"Ambiguity_Group_ID":gid,"Physical_Peak_ID":pid,"Spectrum_ID":best.scan_id or "","Scan_Index":_scan_index(best.scan_id),"Peak_Index":getattr(best,"_audit_peak_index",""),"RT":best.rt,"Observed_MZ":best.observed_mz,"Intensity":best.intensity,"Charge":best.charge,
            "Assignment_Count":len(matches),"Filter_Passing_Assignment_Count":len(passing),"Unique_Fragment_Count":len({str(m.fragment_id) for m in matches}),"Unique_Fragment_Sequence_Count":len(sequences),"Unique_Position_Count":len(positions),"Unique_Length_Count":len({len(m.sequence or "") for m in matches}),"Unique_Theoretical_Mass_Count":len({round(float(m.fragment_mass),8) for m in matches}),"Modification_Eligible_Assignment_Count":sum(bool(candidate_assignments[id(m)]) for m in matches),"Candidate_Supporting_Assignment_Count":sum(bool(candidate_assignments[id(m)]) for m in passing),
            "Same_Sequence_Multiple_Position":same_multi,"Overlapping_Fragment_Count":overlaps,"Nonoverlapping_Fragment_Count":nonoverlaps,"Multiple_Lengths":len({len(m.sequence or "") for m in matches})>1,"Multiple_Missed_Cleavages":len({int(m.missed_cleavages or 0) for m in matches})>1,"Multiple_Terminal_Forms":len({str(m.terminal_form or "") for m in matches})>1,"Near_Isobaric_Assignments":near,"Max_Theoretical_Mass_Separation_Da":mass_sep,"Max_Theoretical_Mass_Separation_PPM":mass_ppm,
            "Best_Fragment_ID":best.fragment_id,"Best_Fragment_Start":best.start,"Best_Fragment_End":best.end,"Best_Fragment_Sequence":best.sequence,"Best_Fragment_Length":len(best.sequence or ""),"Best_Peak_Tier":best.peak_tier,"Best_Confidence":best.confidence,"Best_Abs_Error_PPM":abs(best.mass_error_ppm),"Second_Best_Fragment_ID":second.fragment_id if second else "","Second_Best_Abs_Error_PPM":abs(second.mass_error_ppm) if second else "","Best_Second_Error_Delta_PPM":delta,"Best_Assignment_Dominance":_dominance(best, second, cfg),
            "Full_Count_Total_Weight":total_weights["full_count"],"Winner_Take_All_Total_Weight":total_weights["winner_take_all"],"Equal_Fraction_Total_Weight":total_weights["equal_fraction"],"Quality_Weighted_Total_Weight":total_weights["quality_weighted_fraction"],"Candidate_Key_Count":len(candidate_keys),"Candidate_Key_List":";".join(candidate_keys),"Candidate_Family_List":";".join(families),"Position_Discriminating":position_discriminating,"cnm5U_Relevant":cnm,"Ambiguity_Severity":severity,"Recommended_Handling":"ambiguity_flag_only_pending_candidate_positive_validation","Applied_To_Formal_Result":False,
        })
        for rank, match in enumerate(ordered, 1):
            cs = candidate_assignments[id(match)]; ckeys = sorted({"|".join(map(str,candidate_key(c))) for c in cs}); mids=sorted({str(_raw(c).get("modification_id") or "") for c in cs}); cf=sorted({modification_family(c,modification_lookup) for c in cs})
            reasons=[]
            if passes_fragment_ms1_filter(match,cfg): reasons.append("passes_filter")
            else: reasons.append("fails_filter")
            if rank==1: reasons.append(_dominance(best, second, cfg))
            if len(matches)>1: reasons.append("shared_physical_peak")
            details.append({
                "Dataset_ID":dataset,"Enzyme":enzyme,"Ambiguity_Group_ID":gid,"Physical_Peak_ID":pid,"Fragment_ID":match.fragment_id,"Fragment_Start":match.start,"Fragment_End":match.end,"Fragment_Sequence":match.sequence,"Fragment_Length":len(match.sequence or ""),"Missed_Cleavage_Count":match.missed_cleavages,"Terminal_Form":match.terminal_form,"Observed_MZ":match.observed_mz,"Charge":match.charge,"Intensity":match.intensity,"Peak_Tier":match.peak_tier,"Confidence":match.confidence,"Theoretical_Neutral_Mass":match.fragment_mass,"Observed_Neutral_Mass":match.fragment_mass+match.mass_error_da*abs(match.charge),"Mass_Error_Da":match.mass_error_da,"Mass_Error_PPM":match.mass_error_ppm,"Abs_Mass_Error_PPM":abs(match.mass_error_ppm),"Passes_Fragment_MS1_Filter":passes_fragment_ms1_filter(match,cfg),"Assignment_Rank":rank,"Is_Best_Assignment":rank==1,"Is_Second_Best_Assignment":rank==2,"Same_Sequence_As_Best":same_seq_best(match),"Overlaps_Best_Fragment":_overlap(match,best),"Position_Distance_From_Best":min(abs(int(match.start)-int(best.start)),abs(int(match.end)-int(best.end))),"Theoretical_Mass_Delta_From_Best":float(match.fragment_mass)-float(best.fragment_mass),"Possible_Modification_Count":len(mids),"Possible_Modification_IDs":";".join(mids),"Candidate_Key_List":";".join(ckeys),"Candidate_Family_List":";".join(cf),"Full_Count_Weight":weights["full_count"].get(id(match),0.0),"Winner_Take_All_Weight":weights["winner_take_all"].get(id(match),0.0),"Equal_Fraction_Weight":weights["equal_fraction"].get(id(match),0.0),"Quality_Weighted_Fraction":weights["quality_weighted_fraction"].get(id(match),0.0),"Fragment_Family_Fraction":weights["fragment_family_fraction"].get(id(match),0.0),"Ambiguity_Severity":severity,"Assignment_Reason":";".join(reasons),
            })
    weighting_elapsed = time.perf_counter() - weights_started

    # Winner changes candidate membership; fractional strategies retain all positive-weight assignments.
    shared_ids = set(group_meta); winner_matches=[]
    for pid, matches in groups.items():
        selected = [m for m in matches if id(m) in top50_keys]
        if pid not in shared_ids: winner_matches.extend(selected); continue
        w = assignment_weights(selected,cfg)["winner_take_all"]
        winner_matches.extend(m for m in selected if w.get(id(m),0)>0)
    winner_candidates = search_known_modifications(winner_matches, intact_results, modifications, config, warnings=None)
    winner_ranking, _ = build_modification_evidence_ranking(config, modifications, [item["fragment"] for item in items], winner_matches, winner_candidates, ms2_results, rule_set=rule_set, pathways=pathways)
    changes = _candidate_change_rows(top_ranking, winner_ranking)
    base_keys={candidate_key(c) for c in top_candidates}; winner_keys={candidate_key(c) for c in winner_candidates}
    candidate_changed={"winner_take_all":winner_keys!=base_keys,"equal_fraction":False,"quality_weighted":False,"fragment_family":False}

    all_passing_support=sum(passes_fragment_ms1_filter(m,cfg) for m in top50)
    unshared_passing=sum(passes_fragment_ms1_filter(m,cfg) for m in top50 if physical_peak_id(m) not in shared_ids)
    effective={name:unshared_passing+strategy_support[name] for name in STRATEGIES}
    severity_counts=Counter(row["Ambiguity_Severity"] for row in ambiguity_rows)
    summary_base={
        "Dataset_ID":dataset,"Enzyme":enzyme,"Sequence_Name":dataset,"Total_Physical_Peaks":len(groups),"Shared_Physical_Peaks":len(shared),"Shared_Peak_Fraction":len(shared)/max(1,len(groups)),"Total_Assignments_In_Shared_Peaks":sum(len(v) for _,v in shared),"Passing_Shared_Peaks":sum(meta["passing"]>0 for meta in group_meta.values()),"Multiple_Passing_Assignment_Peaks":sum(meta["passing"]>1 for meta in group_meta.values()),"Same_Sequence_Multi_Position_Peaks":sum(meta["same"] for meta in group_meta.values()),"Overlapping_Fragment_Peaks":sum(meta["overlap"]>0 for meta in group_meta.values()),"Nonoverlapping_Fragment_Peaks":sum(meta["nonoverlap"]>0 for meta in group_meta.values()),"Multi_Length_Peaks":sum(row["Multiple_Lengths"] for row in ambiguity_rows),"Multi_Missed_Cleavage_Peaks":sum(row["Multiple_Missed_Cleavages"] for row in ambiguity_rows),"Multi_Terminal_Form_Peaks":sum(row["Multiple_Terminal_Forms"] for row in ambiguity_rows),"Modification_Eligible_Shared_Peaks":sum(row["Modification_Eligible_Assignment_Count"]>0 for row in ambiguity_rows),"Candidate_Supporting_Shared_Peaks":sum(row["Candidate_Supporting_Assignment_Count"]>0 for row in ambiguity_rows),"Low_Count":severity_counts["Low"],"Medium_Count":severity_counts["Medium"],"High_Count":severity_counts["High"],"Critical_Count":severity_counts["Critical"],
        "Full_Count_Effective_Support":effective["full_count"],"Winner_Take_All_Effective_Support":effective["winner_take_all"],"Equal_Fraction_Effective_Support":effective["equal_fraction"],"Quality_Weighted_Effective_Support":effective["quality_weighted_fraction"],"Fragment_Family_Effective_Support":effective["fragment_family_fraction"],"Equal_Fraction_Support_Reduction":effective["full_count"]-effective["equal_fraction"],
        "Recommended_Ambiguity_Handling":"ambiguity_flag_only","Shared_Peak_Double_Counting_Risk":"High" if len(shared) else "None","Position_Localization_Risk":"High" if any(m["position"] for m in group_meta.values()) else "Low","Candidate_Ranking_Risk":"Unknown_without_candidate_positive_shared_peaks" if not any(m["candidate"] for m in group_meta.values()) else "Medium","Formal_Change_Ready":False,"Evidence_For_Recommendation":f"{len(shared)} of {len(groups)} physical peaks have cross-fragment assignments; {sum(m['passing']>1 for m in group_meta.values())} have multiple passing assignments","Remaining_Risk":"fractional support is not represented by the formal binary MS1 scoring model; short-fragment and candidate-positive behavior remain unresolved","Required_Additional_Validation":"candidate-positive RNase datasets, replicates, position-localization review, and a defined fractional-score model","Ambiguity_Grouping_Time_Seconds":grouping_elapsed,"Shadow_Weighting_Time_Seconds":weighting_elapsed,"Audit_Mode":"shadow_cross_fragment_assignment_ambiguity","Applied_To_Formal_Result":False,
    }
    for label,key in (("Winner_Take_All","winner_take_all"),("Equal_Fraction","equal_fraction"),("Quality_Weighted","quality_weighted"),("Fragment_Family","fragment_family")):
        change=changes[key]
        summary_base.update({f"{label}_Candidate_Set_Changed":candidate_changed[key],f"{label}_Score_Changed_Count":change["score"],f"{label}_Confidence_Changed_Count":change["confidence"],f"{label}_Rank_Changed_Count":change["rank"],f"{label}_Top50_Changed":change["top50"]})
    summary_rows=[dict(summary_base,Summary_Scope="dataset",Fragment_Length_Group="all",Physical_Peak_Count=len(groups),Assignment_Count=len(all_matches),Passing_Assignment_Count=all_passing_support,Average_Assignments_Per_Peak=len(all_matches)/max(1,len(groups)),Best_Second_PPM_Delta_Median=median([m["delta"] for m in group_meta.values()]) if group_meta else "",Position_Localization_Contribution=sum(m["position"] for m in group_meta.values()))]
    for length_group in ("1-2_nt","3_nt","4_nt","5_nt","6-8_nt","9plus_nt"):
        relevant=[(pid,matches) for pid,matches in groups.items() if any(_length_group(len(m.sequence or ""))==length_group for m in matches)]
        rpids={pid for pid,_ in relevant}; shared_meta=[meta for pid,meta in group_meta.items() if pid in rpids]
        row={c:"" for c in SUMMARY_COLUMNS}; row.update({"Dataset_ID":dataset,"Enzyme":enzyme,"Sequence_Name":dataset,"Summary_Scope":"fragment_length","Fragment_Length_Group":length_group,"Physical_Peak_Count":len(rpids),"Shared_Physical_Peaks":len(shared_meta),"Shared_Peak_Fraction":len(shared_meta)/max(1,len(rpids)),"Assignment_Count":sum(sum(_length_group(len(m.sequence or ""))==length_group for m in ms) for _,ms in relevant),"Passing_Assignment_Count":sum(sum(_length_group(len(m.sequence or ""))==length_group and id(m) in top50_keys and passes_fragment_ms1_filter(m,cfg) for m in ms) for _,ms in relevant),"Multiple_Passing_Assignment_Peaks":sum(meta["passing"]>1 for meta in shared_meta),"Average_Assignments_Per_Peak":sum(meta["assignments"] for meta in shared_meta)/max(1,len(shared_meta)),"Best_Second_PPM_Delta_Median":median([meta["delta"] for meta in shared_meta]) if shared_meta else "","High_Count":sum(meta["severity"]=="High" for meta in shared_meta),"Critical_Count":sum(meta["severity"]=="Critical" for meta in shared_meta),"Equal_Fraction_Support_Reduction":sum(meta["weights"]["full_count"]-meta["weights"]["equal_fraction"] for meta in shared_meta),"Candidate_Supporting_Shared_Peaks":sum(meta["candidate"] for meta in shared_meta),"Position_Localization_Contribution":sum(meta["position"] for meta in shared_meta),"Audit_Mode":"shadow_cross_fragment_assignment_ambiguity","Applied_To_Formal_Result":False}); summary_rows.append(row)
    original=len(details); max_rows=int((config.reporting or {}).get("max_excel_rows_per_sheet",100000) or 100000); details=details[:max_rows]
    summary_base.update({"Detail_Original_Row_Count":original,"Detail_Written_Row_Count":len(details),"Detail_Truncated":len(details)<original,"Detail_Truncation_Reason":"max_excel_rows_per_sheet; deterministic physical-peak and assignment rank order" if len(details)<original else ""})
    for key in ("Detail_Original_Row_Count","Detail_Written_Row_Count","Detail_Truncated","Detail_Truncation_Reason"): summary_rows[0][key]=summary_base[key]
    return {"ambiguity_rows":[{c:r.get(c,"") for c in AMBIGUITY_COLUMNS} for r in ambiguity_rows],"detail_rows":[{c:r.get(c,"") for c in DETAIL_COLUMNS} for r in details],"summary_rows":[{c:r.get(c,"") for c in SUMMARY_COLUMNS} for r in summary_rows],"summary":summary_base,"fragment_support":fragment_support,"effective_support":effective,"top_ranking":top_ranking,"winner_ranking":winner_ranking,"changes":changes,"candidate_sets":{"top50":base_keys,"winner_take_all":winner_keys}}


def append_crossfrag_top_columns(rows, audit):
    is_frame=isinstance(rows,pd.DataFrame); source=rows.to_dict("records") if is_frame else list(rows or []); original_columns=list(rows.columns) if is_frame else (list(source[0]) if source else [])
    supports=audit["fragment_support"]; summary=audit["summary"]; top_map={ranking_key(r):r for r in audit["top_ranking"]}; winner_map={ranking_key(r):r for r in audit["winner_ranking"]}
    out=[]
    for original in source:
        row=dict(original); fid=str(row.get("Parent_Fragment_ID") or ""); key=ranking_key(row); full=supports["full_count"][fid]; unique=supports["winner_take_all"][fid]; winner=winner_map.get(key,{}); unchanged=top_map.get(key,{})
        row.update({"CrossFrag_Shared_Peak_Support":full,"CrossFrag_Unique_Peak_Support":unique,"CrossFrag_Ambiguous_Support_Fraction":0.0 if not full else max(0.0,(full-unique)/full),"Winner_Take_All_Shadow_Support":unique,"Equal_Fraction_Shadow_Support":supports["equal_fraction"][fid],"Quality_Weighted_Shadow_Support":supports["quality_weighted_fraction"][fid],"Winner_Take_All_Shadow_Final_Score":winner.get("Final_Score",""),"Equal_Fraction_Shadow_Final_Score":unchanged.get("Final_Score",""),"Quality_Weighted_Shadow_Final_Score":unchanged.get("Final_Score",""),"Winner_Take_All_Shadow_Rank":winner.get("Rank",""),"Equal_Fraction_Shadow_Rank":unchanged.get("Rank",""),"Quality_Weighted_Shadow_Rank":unchanged.get("Rank",""),"CrossFrag_Ambiguity_Affected":full>0,"Recommended_CrossFrag_Handling":summary["Recommended_Ambiguity_Handling"],"CrossFrag_Ambiguity_Applied_To_Formal_Result":False}); out.append(row)
    return pd.DataFrame(out,columns=original_columns+TOP_COLUMNS) if is_frame else out


def append_crossfrag_diagnostic_columns(rows, audit):
    s=audit["summary"]; vals={"MS1_CrossFrag_Audit_Available":True,"MS1_Shared_Physical_Peak_Count":s["Shared_Physical_Peaks"],"MS1_Shared_Physical_Peak_Fraction":s["Shared_Peak_Fraction"],"MS1_Multiple_Passing_Assignment_Peaks":s["Multiple_Passing_Assignment_Peaks"],"MS1_High_Ambiguity_Group_Count":s["High_Count"],"MS1_Critical_Ambiguity_Group_Count":s["Critical_Count"],"MS1_Candidate_Supporting_Shared_Peaks":s["Candidate_Supporting_Shared_Peaks"],"MS1_Winner_Take_All_Candidate_Set_Changed":s["Winner_Take_All_Candidate_Set_Changed"],"MS1_Equal_Fraction_Candidate_Set_Changed":s["Equal_Fraction_Candidate_Set_Changed"],"MS1_Quality_Weighted_Candidate_Set_Changed":s["Quality_Weighted_Candidate_Set_Changed"],"MS1_Winner_Take_All_Rank_Changed":bool(s["Winner_Take_All_Rank_Changed_Count"]),"MS1_Equal_Fraction_Rank_Changed":bool(s["Equal_Fraction_Rank_Changed_Count"]),"MS1_Quality_Weighted_Rank_Changed":bool(s["Quality_Weighted_Rank_Changed_Count"]),"MS1_Recommended_CrossFrag_Handling":s["Recommended_Ambiguity_Handling"],"MS1_CrossFrag_Formal_Change_Ready":False,"MS1_CrossFrag_Applied_To_Formal_Result":False}
    is_frame=isinstance(rows,pd.DataFrame); source=rows.to_dict("records") if is_frame else list(rows or [{}]); out=[dict(r,**vals) for r in source]
    return pd.DataFrame(out,columns=list(rows.columns)+DIAGNOSTIC_COLUMNS) if is_frame else out
