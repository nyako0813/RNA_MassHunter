"""Streaming orchestration for explicit-manifest PT cross-run shadow auditing."""
from __future__ import annotations
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
import bisect
import os
import re
import resource
import statistics
import time
import tracemalloc
from typing import Any, Callable
import numpy as np
from rna_masshunter.cross_run_manifest import INDEPENDENCE_ORDER, load_cross_run_manifest, classify_run_independence, strongest_independence
from rna_masshunter.models import Fragment, Peak
from rna_masshunter.mzml_diagnostics import _rt_minutes
from rna_masshunter.ms2_annotation import extract_ms2_spectra
from rna_masshunter.composite_ms2_matcher import match_composite_ms2
from rna_masshunter.composite_ms2_propagation import generate_composite_theoretical_ions
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.mzml_reader import iter_spectra
from rna_masshunter.peak_filtering import classify_peak_tiers
from rna_masshunter.phosphorothioate_evidence import build_pt_evidence
from rna_masshunter.pt_cross_run_recurrence import (
    FALSE_FLAGS, aggregate_candidates, aggregate_neutral_candidates, aggregate_pairs,
    candidate_key, neutral_candidate_key,
)

@dataclass(frozen=True)
class PTCrossRunAuditResult:
    sheets: dict[str, list[dict[str, Any]]]
    metrics: dict[str, Any]

def _value(item, name, default=""):
    if hasattr(item,name):return getattr(item,name)
    if isinstance(item,dict):
        for key in item:
            if str(key).lower()==name.lower():return item[key]
    return default

def _build_ms2_ions(pairs,config,audit_level):
    ions=[]
    for pair in pairs:
        sequence=pair.spec.sequence
        parent=Fragment(pair.spec.candidate_id+"|parent",pair.spec.sequence_id,sequence[pair.spec.fragment_start-1:pair.spec.fragment_end],
            pair.spec.fragment_start,pair.spec.fragment_end,pair.spec.fragment_start,pair.spec.fragment_end,pair.spec.enzyme,1,
            pair.spec.terminal_form,pair.normal_fragment.neutral_exact_mass)
        ions.extend(generate_composite_theoretical_ions([pair.modified_structure],[parent],sequence,config,audit_level=audit_level))
    return ions

def _ms2_run_audit(ions,spectra,config,audit_level):
    matches=match_composite_ms2(ions,spectra,config,audit_level=audit_level) if ions and spectra else []
    compatible=set();ms2=getattr(config,"ms2_annotation",{}) or {};tolerance=float(ms2.get("precursor_match_tolerance_ppm",20) or 20)
    polarity=str((getattr(config,"instrument",{}) or {}).get("polarity") or "negative").lower()
    for spectrum in spectra or ():
        observed=getattr(spectrum,"precursor_mz",None);charge=getattr(spectrum,"precursor_charge",None)
        if observed in (None,"") or charge in (None,"",0):continue
        for ion in ions:
            theoretical=mz_from_neutral_mass(float(ion["Parent_Neutral_Mass"]),abs(int(charge)),polarity)
            if theoretical and abs((float(observed)-theoretical)/theoretical*1e6)<=tolerance:
                compatible.add(str(ion.get("Candidate_ID") or "").removesuffix("|pt"))
    matched={str(row.get("Candidate_ID") or "").removesuffix("|pt") for row in matches}
    position={str(row.get("Candidate_ID") or "").removesuffix("|pt") for row in matches if row.get("Position_Informative")}
    backbone={str(row.get("Candidate_ID") or "").removesuffix("|pt") for row in matches if row.get("Backbone_Informative")}
    return {"compatible":compatible,"matched":matched,"position":position,"backbone":backbone,"rows":matches}

def _load_run(path: str | Path, reconstruction: dict[str, Any]):
    peaks=[]; ms1=ms2=0; mz_lo=None; mz_hi=None
    rt_min=reconstruction.get("rt_min");rt_max=reconstruction.get("rt_max")
    low=float(reconstruction.get("mz_min",0));high=float(reconstruction.get("mz_max",float("inf")))
    threshold=float(reconstruction.get("intensity_threshold",0))
    for spectrum in iter_spectra(path):
        level=int(spectrum.get("ms level",0)); ms1 += level==1; ms2 += level==2
        if level != 1: continue
        rt=_rt_minutes(spectrum)
        if rt is not None and ((rt_min is not None and rt<float(rt_min)) or (rt_max is not None and rt>float(rt_max))):continue
        mz=np.asarray(spectrum.get("m/z array",[]),dtype=float); inten=np.asarray(spectrum.get("intensity array",[]),dtype=float)
        if mz.size != inten.size:continue
        if mz.size:
            mz_lo=float(mz.min()) if mz_lo is None else min(mz_lo,float(mz.min()))
            mz_hi=float(mz.max()) if mz_hi is None else max(mz_hi,float(mz.max()))
        mask=(mz>=low)&(mz<=high)&(inten>=threshold)
        for m,i in zip(mz[mask],inten[mask],strict=False):peaks.append(Peak(float(m),float(i),rt,str(spectrum.get("id",""))))
    return peaks,{"MS1_Spectrum_Count":ms1,"MS2_Spectrum_Count":ms2,"Acquisition_mz_Min":mz_lo or "","Acquisition_mz_Max":mz_hi or ""}

def _scan_number(value):
    numbers=re.findall(r"\d+",str(value or ""));return int(numbers[-1]) if numbers else None

def _continuity(peaks, theoretical, tolerance):
    width=abs(theoretical)*tolerance/1e6
    near=sorted((p for p in peaks if abs(float(_value(p,"mz",0))-theoretical)<=width),key=lambda p:(_value(p,"rt",0) or 0,_scan_number(_value(p,"scan_id")) or 0))
    if not near:return "Feature_Not_Assessed",0,0,"","",""
    scans=[_scan_number(_value(p,"scan_id")) for p in near];best=1
    if all(x is not None for x in scans):
        current=1
        for a,b in zip(scans,scans[1:]):
            current=current+1 if b==a+1 else 1;best=max(best,current)
    rts=[float(_value(p,"rt")) for p in near if _value(p,"rt") not in (None,"")]
    status="Single_Scan" if len(near)==1 else "Multi_Scan_Continuous" if best>=2 else "Multi_Scan_Discontinuous"
    return status,best,len(near),min(rts) if rts else "",max(rts) if rts else "",(max(rts)-min(rts)) if rts else ""

def _percentile(sorted_intensities, value):
    if value in (None,"") or not sorted_intensities:return ""
    return bisect.bisect_right(sorted_intensities,float(value))/len(sorted_intensities)

def _detail_from_state(state, run, peaks, tolerance):
    matched=bool(state.get("Matched")); observed=state.get("Nearest_Observed_mz") if matched else ""
    intensity=state.get("Intensity") if matched else ""; rt=state.get("RT") if matched else ""
    competition=int(state.get("Competition_Count") or 0); competitors=str(state.get("Competing_Candidate_IDs") or "")
    status,consecutive,total,start,end,span=_continuity(peaks,float(state["Theoretical_mz"]),tolerance) if matched else ("Feature_Not_Assessed",0,0,"","","")
    intensities=sorted(float(_value(p,"intensity",0)) for p in peaks if float(_value(p,"intensity",0))>0)
    row={"Audit_Level":"full","Run_ID":run["run_id"],"Candidate_ID":state.get("Candidate_ID"),
        "Hypothesis_ID":state.get("Candidate_State_ID") if state.get("Candidate_State_ID") in {"H1","H2","H3","H4"} else state.get("Hypothesis_ID"),"Search_Mode":state.get("Search_Mode"),"Sequence_ID":state.get("Sequence_ID"),
        "Enzyme":state.get("Enzyme"),"Bond_ID":state.get("Bond_ID"),"Fragment_ID":f"{state.get('Candidate_ID')}|{state.get('Fragment_Start')}_{state.get('Fragment_End')}",
        "Fragment_Start":state.get("Fragment_Start"),"Fragment_End":state.get("Fragment_End"),"Fragment_Type":"missed_cleavage",
        "Terminal_Form":state.get("Terminal_Form"),"Nucleoside_State":state.get("Shared_Nucleoside_States"),
        "Backbone_State":state.get("Backbone_State"),"Charge":state.get("Charge"),"Elemental_Composition":state.get("Elemental_Composition"),
        "Neutral_Mass":state.get("Neutral_Mass"),"Theoretical_mz":state.get("Theoretical_mz"),"Observable":bool(state.get("Observable")),
        "Matched":matched,"Observed_mz":observed,"Mass_Error_Da":state.get("Mass_Error_Da") if matched else "",
        "Mass_Error_ppm":state.get("Mass_Error_ppm") if matched else "","Intensity":intensity,"Raw_Intensity":intensity,
        "Peak_Tier":state.get("Peak_Tier", ""),"Scan":state.get("Scan") if matched else "","RT":rt,"Observed_RT":rt,
        "Physical_Peak_ID":state.get("Physical_Peak_ID") if matched else "","Consecutive_Scan_Count":consecutive,
        "Total_Nearby_Scan_Count":total,"RT_Feature_Start":start,"RT_Feature_End":end,"RT_Feature_Span":span,
        "Continuity_Status":status,"Run_Total_Eligible_Intensity":sum(intensities),
        "Run_Median_Eligible_Intensity":statistics.median(intensities) if intensities else "",
        "Relative_Intensity_To_Run_Median":float(intensity)/statistics.median(intensities) if intensity not in (None,"") and intensities and statistics.median(intensities)>0 else "",
        "Intensity_Percentile":_percentile(intensities,intensity),
        "Isotope_Ambiguity":matched and competition>0 and "iso" in competitors,
        "Run_Isotope_Status":"AMBIGUOUS" if matched and competition>0 and "iso" in competitors else "CLEAR" if matched else "NOT_DETECTED",
        "Charge_Ambiguity":matched and competition>0 and "|z" in competitors,
        "Legacy_Competition_Count":0,"Composite_Competition_Count":competition,
        "Candidate_Specific":matched and competition==0,"Evidence_Class":"RUN_LEVEL_PT_MS1_MATCH" if matched and state.get("Backbone_State")=="phosphorothioate" else "RUN_LEVEL_MATCH" if matched else "NO_MATCH",
        "RT_Alignment_Status":"unavailable","RT_Alignment_Shift":"","Aligned_RT":"","Aligned_RT_Delta":"","RT_Consistency_Status":"unavailable",**FALSE_FLAGS}
    row["Cross_Run_Candidate_Key"]=candidate_key(row);row["Neutral_Candidate_Key"]=neutral_candidate_key(row)
    return row

def _decoy_rows(state_rows, peaks, run, tolerance):
    """Small nonchemical, charge-preserving m/z shifts; reference-only, never formal."""
    mz_values=sorted(float(_value(p,"mz",0)) for p in peaks)
    seen=set();out=[]
    for state in state_rows:
        core=(state.get("Sequence_ID"),state.get("Enzyme"),state.get("Bond_ID"),state.get("Fragment_Start"),state.get("Fragment_End"),
              state.get("Shared_Nucleoside_States"),state.get("Backbone_State"),state.get("Charge"),state.get("Elemental_Composition"),state.get("Theoretical_mz"))
        if core in seen:continue
        seen.add(core);target=candidate_key({"Sequence_ID":state.get("Sequence_ID"),"Enzyme":state.get("Enzyme"),"Bond_ID":state.get("Bond_ID"),
            "Fragment_Start":state.get("Fragment_Start"),"Fragment_End":state.get("Fragment_End"),"Terminal_Form":state.get("Terminal_Form"),"Charge":state.get("Charge"),
            "Nucleoside_State":state.get("Shared_Nucleoside_States"),"Backbone_State":state.get("Backbone_State"),"Elemental_Composition":state.get("Elemental_Composition"),"Theoretical_mz":state.get("Theoretical_mz")})
        theoretical=float(state.get("Theoretical_mz") or 0)
        for label,shift in (("PLUS_0.050_MZ",0.050),("MINUS_0.047_MZ",-0.047)):
            shifted=theoretical+shift;width=abs(shifted)*tolerance/1e6;left=bisect.bisect_left(mz_values,shifted-width);right=bisect.bisect_right(mz_values,shifted+width)
            out.append({"Audit_Level":"full","Run_ID":run["run_id"],"Target_Candidate_Key":target,"Decoy_Candidate_Key":f"DECOY|{target}|{label}",
                "Decoy_Type":"charge_preserving_shifted_theoretical_mz","Decoy_Shift_mz":shift,"Decoy_Theoretical_mz":shifted,
                "Matched":right>left,"Candidate_Specific":right>left,"Matched_Peak_Count":right-left,**FALSE_FLAGS})
    return out

def _pair_rows(evidence, state_rows, run):
    keys={(s.get("Candidate_ID"),s.get("Backbone_State"),s.get("Charge")):candidate_key({
        "Sequence_ID":s.get("Sequence_ID"),"Enzyme":s.get("Enzyme"),"Fragment_Start":s.get("Fragment_Start"),"Fragment_End":s.get("Fragment_End"),
        "Terminal_Form":s.get("Terminal_Form"),"Charge":s.get("Charge"),"Nucleoside_State":s.get("Shared_Nucleoside_States"),
        "Backbone_State":s.get("Backbone_State"),"Bond_ID":s.get("Bond_ID"),"Elemental_Composition":s.get("Elemental_Composition"),"Theoretical_mz":s.get("Theoretical_mz")}) for s in state_rows}
    out=[]
    for e in evidence:
        n=bool(e.get("Normal_Observed_mz"));p=bool(e.get("Modified_Observed_mz"));amb=e.get("Peak_Ambiguity_Status") in {"same_physical_peak","competing"}
        state="NOT_OBSERVABLE" if not e.get("Observable") else "AMBIGUOUS" if amb else "BOTH_PRESENT" if n and p else "PT_ONLY" if p else "NORMAL_ONLY" if n else "NEITHER_PRESENT"
        pair_key=f"PAIR|{e.get('Sequence_ID')}|{e.get('Enzyme')}|{e.get('Bond_ID')}|{e.get('Fragment_Start')}-{e.get('Fragment_End')}|{e.get('Shared_Nucleoside_States')}|z{e.get('Charge')}"
        out.append({"Run_ID":run["run_id"],"Pair_Key":pair_key,"Normal_Candidate_Key":keys.get((e.get("Candidate_ID"),"normal_phosphate",e.get("Charge")),""),
            "PT_Candidate_Key":keys.get((e.get("Candidate_ID"),"phosphorothioate",e.get("Charge")),""),"Pair_State":state,**FALSE_FLAGS})
    return out

def _align_rt(rows, runs):
    by_key=defaultdict(dict)
    for r in rows:
        if r.get("Matched") and r.get("Candidate_Specific") and r.get("RT") not in (None,""):by_key[r["Cross_Run_Candidate_Key"]][r["Run_ID"]]=float(r["RT"])
    reference=runs[0]["run_id"]; shifts={reference:0.0}
    for run in runs[1:]:
        rid=run["run_id"]
        if run["instrument_method_id"]!=runs[0]["instrument_method_id"] or run["acquisition_batch_id"]!=runs[0]["acquisition_batch_id"]:continue
        deltas=[values[rid]-values[reference] for values in by_key.values() if rid in values and reference in values]
        if len(deltas)>=3:shifts[rid]=statistics.median(deltas)
    for r in rows:
        rid=r["Run_ID"];reference_rt=by_key.get(r["Cross_Run_Candidate_Key"],{}).get(reference,"")
        r["Reference_RT"]=reference_rt
        r["RT_Delta_Raw"]=float(r["RT"])-reference_rt if r.get("RT") not in (None,"") and reference_rt not in (None,"") else ""
        if rid in shifts and r.get("RT") not in (None,""):
            r["RT_Alignment_Status"]="available";r["RT_Alignment_Shift"]=shifts[rid];r["Aligned_RT"]=float(r["RT"])-shifts[rid]
            r["Aligned_RT_Delta"]=r["Aligned_RT"]-reference_rt if reference_rt not in (None,"") else ""
        else:r["RT_Alignment_Status"]="unavailable"

def _deduplicate(rows):
    selected={}
    for r in rows:
        key=(r["Run_ID"],r["Cross_Run_Candidate_Key"])
        old=selected.get(key)
        if old is None or (bool(r.get("Matched")),bool(r.get("Candidate_Specific")),-abs(float(r.get("Mass_Error_ppm") or 1e9))) > (bool(old.get("Matched")),bool(old.get("Candidate_Specific")),-abs(float(old.get("Mass_Error_ppm") or 1e9))):selected[key]=r
    return list(selected.values())

def _targeted(rows,runs,summaries):
    labels={}
    for r in rows:
        if r.get("Search_Mode")=="hypothesis_driven" and r.get("Bond_ID")=="10_11" and r.get("Hypothesis_ID"):
            key=(r["Run_ID"],r["Hypothesis_ID"]);old=labels.get(key)
            if old is None or (bool(r.get("Matched")),-abs(float(r.get("Mass_Error_ppm") or 1e9))) > (bool(old.get("Matched")),-abs(float(old.get("Mass_Error_ppm") or 1e9))):labels[key]=r
    table=[]
    for i,run in enumerate(runs):
        row={"Run_ID":run["run_id"],"Replicate_Type":"REFERENCE" if i==0 else classify_run_independence(runs[0],run)}
        for h in ("H1","H2","H3","H4"):
            item=labels.get((run["run_id"],h),{});row.update({f"{h}_Matched":bool(item.get("Matched")),f"{h}_mz":item.get("Observed_mz",""),f"{h}_ppm":item.get("Mass_Error_ppm",""),f"{h}_Intensity":item.get("Intensity",""),f"{h}_RT":item.get("RT","")})
        row.update(FALSE_FLAGS);table.append(row)
    summary={"Run_ID":"SUMMARY","Replicate_Type":"CROSS_RUN_SUMMARY"}
    for h in ("H1","H2","H3","H4"):
        detected=sum(bool(row.get(f"{h}_Matched")) for row in table);summary[f"{h}_Detection_Count"]=detected;summary[f"{h}_Detection_Rate"]=detected/len(runs) if runs else 0.0
    for h in ("H3","H4"):
        hrows=[r for r in rows if r.get("Hypothesis_ID")==h and r.get("Matched")]
        summary[f"{h}_Candidate_Specific_Run_Count"]=sum(bool(r.get("Candidate_Specific")) for r in hrows)
        matched_meta=[next(x for x in runs if x["run_id"]==r["Run_ID"]) for r in hrows]
        level=strongest_independence(matched_meta)
        summary[f"{h}_Independent_Run_Count"]=len(hrows) if INDEPENDENCE_ORDER[level]>=INDEPENDENCE_ORDER["INDEPENDENT_DIGESTION"] else 0
    targeted_classes=[x["Recurrence_Evidence_Class"] for x in summaries if x.get("Hypothesis_ID") in {"H3","H4"}]
    priority=("RECURRENT_INDEPENDENT_PT_MS1_SUPPORT","RECURRENT_PT_MS1_SUPPORT","RECURRENT_BUT_AMBIGUOUS","RECURRENT_MASS_INCONSISTENT","SINGLE_RUN_PT_MS1_SUPPORT","NO_RECURRENT_SUPPORT","NOT_EVALUABLE")
    summary["Final_Bond_10_11_Recurrence_Class"]=next((x for x in priority if x in targeted_classes),"NOT_EVALUABLE")
    summary.update(FALSE_FLAGS);table.append(summary)
    return table

def build_pt_cross_run_audit(manifest_path, pairs, config, *, audit_level="audit", legacy_matches=(), other_composite_matches=(), peak_loader: Callable|None=None):
    if audit_level=="standard":return PTCrossRunAuditResult({}, {"Manifest_Loaded":False,"Reason":"standard_mode"})
    started=time.perf_counter();tracemalloc.start();manifest=load_cross_run_manifest(manifest_path,require_files=False);runs=list(manifest.runs)
    all_detail=[];all_pairs=[];all_decoys=[];run_rows=[];per_run={};ms2_by_run={};ms2_detail=[]
    tolerance=float((getattr(config,"fragment_mapping",{}) or {}).get("mz_tolerance_ppm",10) or 10)
    ms2_ions=_build_ms2_ions(pairs,config,audit_level)
    for run in runs:
        tick=time.perf_counter()
        run_valid=True;invalid_reason=""
        if peak_loader:
            loaded=peak_loader(run);peaks,diag=loaded if isinstance(loaded,tuple) else (loaded,{})
        elif not Path(run["mzml_path"]).is_file():
            peaks=[];diag={"MS1_Spectrum_Count":0,"MS2_Spectrum_Count":0};run_valid=False;invalid_reason=f"missing_mzml_file:{run['mzml_path']}"
        else:
            try:peaks,diag=_load_run(run["mzml_path"],getattr(config,"reconstruction",{}) or {})
            except Exception as exc:
                peaks=[];diag={"MS1_Spectrum_Count":0,"MS2_Spectrum_Count":0};run_valid=False;invalid_reason=f"mzml_read_error:{type(exc).__name__}:{exc}"
        classify_peak_tiers(peaks,getattr(config,"peak_filtering",{}) or {})
        if peak_loader:
            spectra=diag.pop("_MS2_Spectra",[])
        elif run_valid and int(diag.get("MS2_Spectrum_Count") or 0)>0:
            spectra=extract_ms2_spectra(run["mzml_path"],getattr(config,"ms2_annotation",{}) or {})
        else:spectra=[]
        ms2_by_run[run["run_id"]]=_ms2_run_audit(ms2_ions,spectra,config,audit_level)
        if audit_level=="full":
            ms2_detail.extend(dict(row,Run_ID=run["run_id"],Formal_Result_Changed=False) for row in ms2_by_run[run["run_id"]]["rows"])
        evidence,states=build_pt_evidence(list(pairs),peaks,config,legacy_matches=legacy_matches,other_composite_matches=other_composite_matches,audit_level=audit_level,include_detail=True)
        details=[_detail_from_state(s,run,peaks,tolerance) for s in states]
        all_detail.extend(details);all_pairs.extend(_pair_rows(evidence,states,run));all_decoys.extend(_decoy_rows(states,peaks,run,tolerance));per_run[run["run_id"]]=time.perf_counter()-tick
        intensities=[float(_value(p,"intensity",0)) for p in peaks]
        run_rows.append({"Run_ID":run["run_id"],"mzML_Path":run["mzml_path"],"Sample_ID":run["sample_id"],
            "Biological_Replicate_ID":run["biological_replicate_id"],"Sample_Preparation_ID":run["sample_preparation_id"],"Digestion_ID":run["digestion_id"],
            "Technical_Replicate_ID":run["technical_replicate_id"],"Acquisition_Batch_ID":run["acquisition_batch_id"],"Instrument_Method_ID":run["instrument_method_id"],
            "Condition":run["condition"],"Enzyme":run["enzyme"],"Sequence_ID":run["sequence_id"],"Organism":run["organism"],"Notes":run.get("notes",""),
            "Run_Valid":run_valid,"Invalid_Reason":invalid_reason,"MS1_Spectrum_Count":diag.get("MS1_Spectrum_Count",""),"MS2_Spectrum_Count":diag.get("MS2_Spectrum_Count",""),
            "Eligible_Peak_Count":len(peaks),"Acquisition_mz_Min":diag.get("Acquisition_mz_Min",min((_value(p,"mz") for p in peaks),default="")),
            "Acquisition_mz_Max":diag.get("Acquisition_mz_Max",max((_value(p,"mz") for p in peaks),default="")),
            "Run_Total_Eligible_Intensity":sum(intensities),"Run_Median_Eligible_Intensity":statistics.median(intensities) if intensities else "",**FALSE_FLAGS})
        del peaks,states,evidence,details
    all_detail=_deduplicate(all_detail);_align_rt(all_detail,runs)
    pair_state_index={}
    for item in all_pairs:
        for key_name in ("Normal_Candidate_Key","PT_Candidate_Key"):
            if item.get(key_name):pair_state_index[(item["Run_ID"],item[key_name])]=item["Pair_State"]
    for row in all_detail:row["Pair_State"]=pair_state_index.get((row["Run_ID"],row["Cross_Run_Candidate_Key"]),"")
    aggregation_started=time.perf_counter()
    summary=aggregate_candidates(all_detail,runs,tolerance_ppm=tolerance);neutral=aggregate_neutral_candidates(all_detail,runs,tolerance_ppm=tolerance);pairs_summary=aggregate_pairs(all_pairs,runs)
    decoy_groups=defaultdict(list)
    for row in all_decoys:decoy_groups[row["Decoy_Candidate_Key"]].append(row)
    decoy_recurrent=sum(sum(bool(x["Matched"]) for x in group)>=2 for group in decoy_groups.values())
    decoy_specific_recurrent=sum(sum(bool(x["Matched"] and x["Candidate_Specific"]) for x in group)>=2 for group in decoy_groups.values())
    decoy=[{"Audit_Level":audit_level,"Target_Candidate_Count":len(summary),"Target_Recurrent_Count":sum(str(r["Recurrence_Evidence_Class"]).startswith("RECURRENT") for r in summary),
        "Decoy_Candidate_Count":len(decoy_groups),"Decoy_Recurrent_Count":decoy_recurrent,"Target_Candidate_Specific_Recurrent_Count":sum(r["Recurrence_Evidence_Class"] in {"RECURRENT_PT_MS1_SUPPORT","RECURRENT_INDEPENDENT_PT_MS1_SUPPORT"} for r in summary),
        "Decoy_Candidate_Specific_Recurrent_Count":decoy_specific_recurrent,"Decoy_Status":"reference_only_not_FDR","Decoy_Recurrence_Rate":decoy_recurrent/len(decoy_groups) if decoy_groups else 0.0,**FALSE_FLAGS}]
    current,peak=tracemalloc.get_traced_memory();tracemalloc.stop();metrics={"Run_Count":len(runs),"Total_Wall_Time":time.perf_counter()-started,"Per_Run_Runtime":per_run,
        "Cross_Run_Aggregation_Runtime":time.perf_counter()-aggregation_started,"Maximum_RSS_MiB":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        "Cross_Run_Tracemalloc_Peak_MiB":peak/1024/1024,"Detail_Row_Count":len(all_detail) if audit_level=="full" else 0,"Summary_Row_Count":len(summary)}
    ms2_evaluable=sum(int(row.get("MS2_Spectrum_Count") or 0)>0 for row in run_rows)
    for row in summary:
        cid=str(row.get("Candidate_ID") or "")
        compatible=sum(cid in result["compatible"] for result in ms2_by_run.values());matched=sum(cid in result["matched"] for result in ms2_by_run.values())
        position=sum(cid in result["position"] for result in ms2_by_run.values());backbone=sum(cid in result["backbone"] for result in ms2_by_run.values())
        row.update({"MS2_Evaluable_Run_Count":ms2_evaluable,"MS2_Precursor_Compatible_Run_Count":compatible,"MS2_Matched_Run_Count":matched,
            "Position_Localizing_Run_Count":position,"Backbone_Localizing_Run_Count":backbone,
            "MS2_Cross_Run_Status":"no_ms2_spectra" if not ms2_evaluable else "matched" if matched else "no_candidate_ms2_match"})
    for row in run_rows:row.update(metrics)
    sheets={"PT_Cross_Run_Runs":run_rows,"PT_Cross_Run_Summary":summary,"PT_Cross_Run_Neutral":neutral,"PT_Cross_Run_Pairs":pairs_summary,"PT_Cross_Run_Decoy":decoy,"PT_Cross_Run_Target":_targeted(all_detail,runs,summary)}
    if audit_level=="full":
        sheets["PT_Cross_Run_Detail"]=all_detail
        sheets["PT_Cross_Run_Decoy_Detail"]=all_decoys
        sheets["PT_Cross_Run_MS2_Detail"]=ms2_detail
    return PTCrossRunAuditResult(sheets,metrics)
