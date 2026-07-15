"""Shadow-only PT cross-run keys, recurrence statistics, classes, and pair summaries."""
from __future__ import annotations
from collections import Counter, defaultdict
import hashlib
import json
import math
import statistics
from typing import Any
from rna_masshunter.cross_run_manifest import INDEPENDENCE_ORDER, classify_run_independence, strongest_independence

FALSE_FLAGS = {"Applied_To_Formal_Result": False, "Formal_Change_Ready": False, "Formal_Result_Changed": False}

def _canon(value: Any) -> str:
    if value in (None, ""): return "unmodified"
    if isinstance(value, (list, tuple, set)): return ";".join(sorted(_canon(x) for x in value))
    return str(value).strip()

def candidate_key(row: dict[str, Any], *, include_charge: bool = True) -> str:
    fields = {
        "sequence": _canon(row.get("Sequence_ID")), "enzyme": _canon(row.get("Enzyme")),
        "start": int(row.get("Fragment_Start") or 0), "end": int(row.get("Fragment_End") or 0),
        "fragment_type": _canon(row.get("Fragment_Type") or "missed_cleavage"),
        "terminal": _canon(row.get("Terminal_Form")),
        "charge": int(row.get("Charge") or 0) if include_charge else 0,
        "nucleoside": _canon(row.get("Nucleoside_State") or row.get("Shared_Nucleoside_States")),
        "backbone": _canon(row.get("Backbone_State")), "bond": _canon(row.get("Bond_ID")),
        "composition": _canon(row.get("Elemental_Composition")),
        "theoretical_mz": round(float(row.get("Theoretical_mz") or 0), 8) if include_charge else 0,
    }
    digest = hashlib.sha256(json.dumps(fields, sort_keys=True, separators=(",", ":")).encode()).hexdigest()[:20]
    prefix = f"{fields['sequence']}|{fields['enzyme']}|{fields['start']}-{fields['end']}|{fields['backbone']}"
    return f"CRK|{prefix}|{digest}"

def neutral_candidate_key(row: dict[str, Any]) -> str:
    base = dict(row); base["Charge"] = 0; base["Theoretical_mz"] = 0
    return candidate_key(base, include_charge=False).replace("CRK|", "NCK|", 1)

def _num(values):
    return [float(x) for x in values if x not in (None, "") and math.isfinite(float(x))]

def _sd(values): return statistics.stdev(values) if len(values) > 1 else 0.0 if values else ""
def _mad(values):
    if not values: return ""
    med = statistics.median(values); return statistics.median(abs(x-med) for x in values)

def _sign_consistency(values):
    if not values: return "NOT_EVALUABLE"
    signs = {0 if x == 0 else (1 if x > 0 else -1) for x in values}
    return "CONSISTENT" if len(signs) == 1 else "MIXED"

def mass_consistent(values: list[float], tolerance_ppm: float) -> bool:
    return len(values) >= 2 and (_mad(values) <= max(1.0, tolerance_ppm * 0.25)) and (max(values)-min(values) <= tolerance_ppm)

def recurrence_class(rows: list[dict[str, Any]], *, tolerance_ppm: float = 10.0,
    normal_detected_runs: int = 0) -> tuple[str, str]:
    evaluable = [r for r in rows if bool(r.get("Observable"))]
    detected = [r for r in evaluable if bool(r.get("Matched"))]
    pt = bool(rows and rows[0].get("Backbone_State") == "phosphorothioate")
    if len(evaluable) < 2: return "NOT_EVALUABLE", "fewer than two evaluable runs"
    if not pt:
        return ("NO_RECURRENT_SUPPORT", "normal-state candidate summary")
    if normal_detected_runs >= 2 and len(detected) < 2: return "NORMAL_DOMINANT", "normal counterpart recurrent; PT is not recurrent"
    if normal_detected_runs and len(detected) >= 2: return "MIXED_NORMAL_AND_PT", "normal and PT counterparts recur"
    if not detected: return "NO_RECURRENT_SUPPORT", "no PT match in evaluable runs"
    if len(detected) == 1: return "SINGLE_RUN_PT_MS1_SUPPORT", "PT matched in one of multiple evaluable runs"
    ppm = _num(r.get("Mass_Error_ppm") for r in detected)
    if not mass_consistent(ppm, tolerance_ppm): return "RECURRENT_MASS_INCONSISTENT", "PT proximity recurs but ppm errors are inconsistent"
    ambiguous = any(not bool(r.get("Candidate_Specific")) or bool(r.get("Isotope_Ambiguity")) or
        bool(r.get("Charge_Ambiguity")) or int(r.get("Legacy_Competition_Count") or 0) or
        int(r.get("Composite_Competition_Count") or 0) for r in detected)
    if ambiguous: return "RECURRENT_BUT_AMBIGUOUS", "recurrence retains isotope/charge/candidate competition"
    level = strongest_independence([r["_run_metadata"] for r in detected])
    if INDEPENDENCE_ORDER[level] >= INDEPENDENCE_ORDER["INDEPENDENT_DIGESTION"]:
        return "RECURRENT_INDEPENDENT_PT_MS1_SUPPORT", f"candidate-specific, mass-consistent recurrence at {level}"
    return "RECURRENT_PT_MS1_SUPPORT", f"candidate-specific, mass-consistent recurrence at {level}"

def independence_evidence_level(level: str) -> tuple[str, float]:
    mapping = {
        "BIOLOGICAL_REPLICATE": ("BIOLOGICALLY_REPLICATED", 1.0),
        "INDEPENDENT_SAMPLE_PREPARATION": ("PREPARATION_REPLICATED", .85),
        "INDEPENDENT_DIGESTION": ("DIGESTION_REPLICATED", .7),
        "INDEPENDENT_INJECTION": ("INJECTION_REPLICATED", .5),
        "TECHNICAL_REPLICATE": ("TECHNICAL_ONLY", .35),
        "SAME_INJECTION": ("SINGLE_INJECTION", .1),
        "UNKNOWN_INDEPENDENCE": ("UNKNOWN", 0.0),
    }
    return mapping[level]

def aggregate_candidates(detail_rows: list[dict[str, Any]], run_metadata: list[dict[str, Any]], *, tolerance_ppm=10.0):
    run_by_id = {r["run_id"]: r for r in run_metadata}; groups = defaultdict(list)
    for row in detail_rows: groups[row["Cross_Run_Candidate_Key"]].append(row)
    normal_by_pair = Counter()
    for rows in groups.values():
        sample = rows[0]
        if sample.get("Backbone_State") == "normal_phosphate":
            normal_by_pair[(sample.get("Candidate_ID"), sample.get("Charge"))] = sum(bool(r.get("Matched")) for r in rows)
    summaries=[]
    for key, rows in sorted(groups.items()):
        detected=[r for r in rows if r.get("Matched")]; evaluable=[r for r in rows if r.get("Observable")]
        for r in rows: r["_run_metadata"] = run_by_id[r["Run_ID"]]
        sample=rows[0]; normal_count=normal_by_pair[(sample.get("Candidate_ID"),sample.get("Charge"))]
        klass, reason=recurrence_class(rows,tolerance_ppm=tolerance_ppm,normal_detected_runs=normal_count)
        ppm=_num(r.get("Mass_Error_ppm") for r in detected); rts=_num(r.get("RT") for r in detected)
        aligned=_num(r.get("Aligned_RT") for r in detected); intens=_num(r.get("Intensity_Percentile") for r in detected)
        detected_meta=[run_by_id[r["Run_ID"]] for r in detected]; level=strongest_independence(detected_meta)
        evidence_level,weight=independence_evidence_level(level)
        charge_counts=Counter(int(r.get("Charge") or 0) for r in detected); dominant=charge_counts.most_common(1)[0] if charge_counts else ("",0)
        def replicate_count(target):
            if len(detected_meta)<2:return 0
            return sum(any(classify_run_independence(m,o)==target for o in detected_meta if o is not m) for m in detected_meta)
        summary={
            "Cross_Run_Candidate_Key":key,"Neutral_Candidate_Key":sample.get("Neutral_Candidate_Key"),
            "Candidate_ID":sample.get("Candidate_ID"),"Hypothesis_ID":sample.get("Hypothesis_ID"),"Search_Mode":sample.get("Search_Mode"),
            "Sequence_ID":sample.get("Sequence_ID"),"Enzyme":sample.get("Enzyme"),"Bond_ID":sample.get("Bond_ID"),
            "Fragment_Start":sample.get("Fragment_Start"),"Fragment_End":sample.get("Fragment_End"),
            "Nucleoside_State":sample.get("Nucleoside_State"),"Backbone_State":sample.get("Backbone_State"),
            "Charge_or_Neutral_Summary":sample.get("Charge"),"Total_Run_Count":len(run_metadata),
            "Evaluable_Run_Count":len(evaluable),"Detected_Run_Count":len(detected),
            "Detection_Rate":len(detected)/len(evaluable) if evaluable else 0.0,
            "Candidate_Specific_Run_Count":sum(bool(r.get("Candidate_Specific")) for r in detected),
            "Candidate_Specific_Detection_Rate":sum(bool(r.get("Candidate_Specific")) for r in detected)/len(evaluable) if evaluable else 0.0,
            "Ambiguous_Run_Count":sum(not bool(r.get("Candidate_Specific")) for r in detected),
            "Technical_Replicate_Detected_Count":replicate_count("TECHNICAL_REPLICATE"),
            "Independent_Injection_Detected_Count":replicate_count("INDEPENDENT_INJECTION"),
            "Independent_Digestion_Detected_Count":replicate_count("INDEPENDENT_DIGESTION"),
            "Independent_Preparation_Detected_Count":replicate_count("INDEPENDENT_SAMPLE_PREPARATION"),
            "Biological_Replicate_Detected_Count":replicate_count("BIOLOGICAL_REPLICATE"),
            "Mean_Error_ppm":statistics.mean(ppm) if ppm else "","Median_Error_ppm":statistics.median(ppm) if ppm else "",
            "SD_Error_ppm":_sd(ppm),"MAD_Error_ppm":_mad(ppm),"Min_Error_ppm":min(ppm) if ppm else "",
            "Max_Error_ppm":max(ppm) if ppm else "","Error_Sign_Consistency":_sign_consistency(ppm),
            "Observed_Charge_Set":";".join(map(str,sorted(charge_counts))),"Dominant_Charge":dominant[0],
            "Dominant_Charge_Run_Count":dominant[1],"Charge_Consistency_Rate":dominant[1]/len(detected) if detected else 0.0,
            "Charge_Consistency_Status":"CONSISTENT" if detected and dominant[1]==len(detected) else "VARIABLE" if detected else "NOT_EVALUABLE",
            "Mean_RT":statistics.mean(rts) if rts else "","Median_RT":statistics.median(rts) if rts else "","RT_SD":_sd(rts),
            "Aligned_RT_SD":_sd(aligned),"RT_Consistency_Status":("ALIGNED_CONSISTENT" if len(aligned)>=2 and _sd(aligned) <= 0.5 else "ALIGNED_INCONSISTENT" if len(aligned)>=2 else "RAW_CONSISTENT" if len(rts)>=2 and _sd(rts) <= 0.5 else "RAW_INCONSISTENT" if len(rts)>=2 else "unavailable"),
            "Mean_Intensity_Percentile":statistics.mean(intens) if intens else "","Median_Intensity_Percentile":statistics.median(intens) if intens else "",
            "Multi_Scan_Run_Count":sum(r.get("Continuity_Status") in {"Multi_Scan_Continuous","Multi_Scan_Discontinuous"} for r in detected),
            "Isotope_Clear_Run_Count":sum(not bool(r.get("Isotope_Ambiguity")) for r in detected),
            "Isotope_Ambiguous_Run_Count":sum(bool(r.get("Isotope_Ambiguity")) for r in detected),
            "Isotope_Clear_Detection_Rate":sum(not bool(r.get("Isotope_Ambiguity")) for r in detected)/len(detected) if detected else 0.0,
            "Competition_Free_Run_Count":sum(bool(r.get("Candidate_Specific")) for r in detected),
            "Recurrence_Evidence_Class":klass,"Recurrence_Independence_Level":evidence_level,
            "Run_Independence_Class":level,"Recurrence_Weight_Suggestion":weight,"Evidence_Reason":reason,**FALSE_FLAGS}
        summaries.append(summary)
        for r in rows:r.pop("_run_metadata",None)
    return summaries

def aggregate_neutral_candidates(detail_rows, run_metadata, *, tolerance_ppm=10.0):
    groups=defaultdict(list)
    for row in detail_rows: groups[row["Neutral_Candidate_Key"]].append(row)
    out=[]
    for key, rows in sorted(groups.items()):
        by_run=defaultdict(list)
        for r in rows: by_run[r["Run_ID"]].append(r)
        detected_runs=[rid for rid,rr in by_run.items() if any(x.get("Matched") for x in rr)]
        charges=sorted({int(x.get("Charge") or 0) for x in rows if x.get("Matched")})
        out.append({"Neutral_Candidate_Key":key,"Candidate_ID":rows[0].get("Candidate_ID"),
            "Backbone_State":rows[0].get("Backbone_State"),"Total_Run_Count":len(run_metadata),
            "Detected_Run_Count":len(detected_runs),"Observed_Charge_Set":";".join(map(str,charges)),**FALSE_FLAGS})
    return out

def aggregate_pairs(pair_rows: list[dict[str, Any]], run_metadata: list[dict[str, Any]]):
    groups=defaultdict(list)
    for row in pair_rows:groups[row["Pair_Key"]].append(row)
    out=[]; run_by_id={r["run_id"]:r for r in run_metadata}
    for key,rows in sorted(groups.items()):
        counts=Counter(r["Pair_State"] for r in rows); pt_meta=[run_by_id[r["Run_ID"]] for r in rows if r["Pair_State"] in {"PT_ONLY","BOTH_PRESENT"}]
        normal_meta=[run_by_id[r["Run_ID"]] for r in rows if r["Pair_State"] in {"NORMAL_ONLY","BOTH_PRESENT"}]
        if counts["BOTH_PRESENT"] or (counts["PT_ONLY"] and counts["NORMAL_ONLY"]): klass="MIXED_NORMAL_AND_PT"
        elif counts["PT_ONLY"]>=2: klass="PT_RECURRENT"
        elif counts["NORMAL_ONLY"]>=2: klass="NORMAL_DOMINANT"
        elif len(rows)<2:klass="NOT_EVALUABLE"
        else:klass="NO_RECURRENT_SUPPORT"
        out.append({"Pair_Key":key,"Normal_Candidate_Key":rows[0]["Normal_Candidate_Key"],"PT_Candidate_Key":rows[0]["PT_Candidate_Key"],
            "Total_Run_Count":len(run_metadata),"Evaluable_Run_Count":sum(r["Pair_State"]!="NOT_OBSERVABLE" for r in rows),
            "PT_Only_Run_Count":counts["PT_ONLY"],"Normal_Only_Run_Count":counts["NORMAL_ONLY"],"Both_Present_Run_Count":counts["BOTH_PRESENT"],
            "Neither_Run_Count":counts["NEITHER_PRESENT"],"Not_Observable_Run_Count":counts["NOT_OBSERVABLE"],"Ambiguous_Run_Count":counts["AMBIGUOUS"],
            "Independent_PT_Detection_Count":len(pt_meta) if INDEPENDENCE_ORDER[strongest_independence(pt_meta)]>=4 else 0,
            "Independent_Normal_Detection_Count":len(normal_meta) if INDEPENDENCE_ORDER[strongest_independence(normal_meta)]>=4 else 0,
            "Pair_Recurrence_Class":klass,**FALSE_FLAGS})
    return out
