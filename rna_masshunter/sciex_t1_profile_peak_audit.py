"""Deterministic shadow peak audit for SCIEX RNase-T1 m/z profiles."""
from __future__ import annotations
from dataclasses import dataclass, replace
from enum import Enum
from hashlib import sha256
from math import ceil
from typing import Sequence
import numpy as np
from scipy.signal import find_peaks, peak_prominences, peak_widths

ALGORITHM_VERSION = "sciex-t1-profile-peak-audit-v1"
ISOTOPE_MASS_DIFFERENCE_DA = 1.003355

class T1PeakQualityClass(str, Enum):
    MAJOR_SHARP="MAJOR_SHARP"; MAJOR_BROAD="MAJOR_BROAD"
    MINOR_SHARP="MINOR_SHARP"; MINOR_BROAD="MINOR_BROAD"
    SHOULDER_OR_OVERLAP="SHOULDER_OR_OVERLAP"
    ISOTOPE_OR_ENVELOPE_COMPONENT="ISOTOPE_OR_ENVELOPE_COMPONENT"
    LOW_SUPPORT="LOW_SUPPORT"

@dataclass(frozen=True)
class T1PeakDetectionParameters:
    minimum_relative_apex_intensity: float = 0.0005
    minimum_relative_integrated_intensity: float = 0.0002
    minimum_relative_prominence: float = 0.0005
    minor_sharp_minimum_relative_prominence: float = 0.001
    major_relative_intensity: float = 0.01
    minimum_peak_separation_mz: float = 0.015
    maximum_peak_width_mz: float = 0.15
    sharp_fwhm_max_mz: float = 0.05
    maximum_selected_peaks_per_profile: int = 300
    isotope_spacing_tolerance_mz: float = 0.01
    shared_valley_fraction: float = 0.60
    duplicate_separation_mz: float = 0.03
    def validate(self):
        for name in ("minimum_relative_apex_intensity","minimum_relative_integrated_intensity",
                     "minimum_relative_prominence","minor_sharp_minimum_relative_prominence",
                     "major_relative_intensity","shared_valley_fraction"):
            if not 0 <= float(getattr(self,name)) <= 1: raise ValueError(f"{name} must be in [0,1]")
        for name in ("minimum_peak_separation_mz","maximum_peak_width_mz","sharp_fwhm_max_mz",
                     "isotope_spacing_tolerance_mz","duplicate_separation_mz"):
            if float(getattr(self,name)) <= 0: raise ValueError(f"{name} must be positive")
        if self.maximum_selected_peaks_per_profile < 1: raise ValueError("selection bound must be positive")

@dataclass(frozen=True, kw_only=True)
class T1Safeguards:
    shadow_analysis_only: bool=True; mass_evidence_only: bool=True
    rna_identity_confirmed: bool=False; target_rna_identity_confirmed_by_mass: bool=False
    fragment_assigned: bool=False; unique_fragment_assigned: bool=False
    modification_assigned: bool=False; modification_composition_assigned: bool=False
    position_assigned: bool=False; structure_assigned: bool=False
    charge_state_confirmed: bool=False; ion_mode_confirmed: bool=False; polarity_confirmed: bool=False
    isotope_assigned: bool=False; co_captured_rna_excluded: bool=False
    background_component_excluded: bool=False; contaminant_excluded: bool=False
    biological_cause_assigned: bool=False
    applied_to_formal_score: bool=False; applied_to_ranking: bool=False
    applied_to_candidate_filtering: bool=False; applied_to_final_consensus: bool=False

@dataclass(frozen=True, kw_only=True)
class T1ProfilePeak(T1Safeguards):
    t1_peak_id: str; source_id: str; measurement_id: str; rna_identity: str
    apex_mz: float; centroid_mz: float | None; apex_intensity: float
    integrated_intensity: float; relative_apex_intensity: float
    relative_integrated_intensity: float; left_boundary_mz: float; right_boundary_mz: float
    peak_width_mz: float; fwhm_mz: float; prominence: float; relative_prominence: float
    sharpness_score: float; nearest_peak_separation_mz: float | None
    shared_valley_fraction: float | None
    possible_shoulder: bool; possible_duplicate: bool; possible_isotope_component: bool
    possible_isotope_charge: int | None; possible_isotope_spacing_error: float | None
    possible_overlapping_envelope: bool; peak_quality_class: T1PeakQualityClass
    selected_as_primary: bool=False

@dataclass(frozen=True)
class T1ProfilePeakAuditResult:
    status: str; source_id: str; measurement_id: str; rna_identity: str
    coordinate_type: str; observed_mass_scale: str; observed_output_species: str
    parameters: T1PeakDetectionParameters; peaks: tuple[T1ProfilePeak,...]
    selected_peaks: tuple[T1ProfilePeak,...]; algorithm_version: str=ALGORITHM_VERSION

_QUALITY_ORDER={T1PeakQualityClass.MAJOR_SHARP:0,T1PeakQualityClass.MAJOR_BROAD:1,
 T1PeakQualityClass.MINOR_SHARP:2,T1PeakQualityClass.MINOR_BROAD:3,
 T1PeakQualityClass.ISOTOPE_OR_ENVELOPE_COMPONENT:4,T1PeakQualityClass.SHOULDER_OR_OVERLAP:5,
 T1PeakQualityClass.LOW_SUPPORT:6}

def _interp(x: np.ndarray, index: float) -> float:
    low=int(np.floor(index)); high=min(low+1,len(x)-1); frac=index-low
    return float(x[low]*(1-frac)+x[high]*frac)

def _isotope_annotations(peaks: list[T1ProfilePeak], tolerance: float):
    ordered=sorted(range(len(peaks)),key=lambda i:(peaks[i].apex_mz,peaks[i].t1_peak_id))
    annotations={i:[] for i in ordered}
    for pos,i in enumerate(ordered):
        for j in ordered[pos+1:]:
            delta=peaks[j].apex_mz-peaks[i].apex_mz
            if delta > ISOTOPE_MASS_DIFFERENCE_DA+tolerance: break
            for charge in range(1,6):
                error=delta-ISOTOPE_MASS_DIFFERENCE_DA/charge
                if abs(error)<=tolerance: annotations[j].append((abs(error),charge,error,i))
    return annotations

def detect_t1_profile_peaks(coordinates: Sequence[float], intensities: Sequence[float], *,
                            source_id: str, measurement_id: str, rna_identity: str,
                            parameters: T1PeakDetectionParameters|None=None) -> T1ProfilePeakAuditResult:
    params=parameters or T1PeakDetectionParameters(); params.validate()
    x=np.asarray(tuple(coordinates),dtype=float); y=np.asarray(tuple(intensities),dtype=float)
    if x.ndim!=1 or y.ndim!=1 or len(x)!=len(y) or len(x)<3: raise ValueError("aligned 1-D profile required")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or np.any(np.diff(x)<=0): raise ValueError("finite increasing m/z axis required")
    y=np.maximum(y,0); base=float(y.max())
    if base<=0: return T1ProfilePeakAuditResult("NO_PEAKS",source_id,measurement_id,rna_identity,"MZ","MZ","CHARGED_ION_UNKNOWN",params,(),())
    spacing=float(np.median(np.diff(x))); distance=max(1,int(ceil(params.minimum_peak_separation_mz/spacing)))
    indices,_=find_peaks(y,distance=distance,prominence=base*params.minimum_relative_prominence)
    if not len(indices): return T1ProfilePeakAuditResult("NO_PEAKS",source_id,measurement_id,rna_identity,"MZ","MZ","CHARGED_ION_UNKNOWN",params,(),())
    prominences,left_bases,right_bases=peak_prominences(y,indices)
    widths,_,left_ips,right_ips=peak_widths(y,indices,rel_height=0.5,
                                            prominence_data=(prominences,left_bases,right_bases))
    total_area=float(np.trapezoid(y,x)); interim=[]
    for serial,(idx,prom,lb,rb,lip,rip) in enumerate(zip(indices,prominences,left_bases,right_bases,left_ips,right_ips),1):
        apex_mz=float(x[idx]); half_window=params.maximum_peak_width_mz/2
        local_left=max(0,int(np.searchsorted(x,apex_mz-half_window,side="left")))
        local_right=min(len(x)-1,int(np.searchsorted(x,apex_mz+half_window,side="right"))-1)
        lb=max(int(lb),local_left); rb=min(int(rb),local_right)
        segment_x=x[lb:rb+1]; segment_y=y[lb:rb+1]
        area=float(np.trapezoid(segment_y,segment_x)); centroid=float(np.average(segment_x,weights=segment_y)) if segment_y.sum()>0 else apex_mz
        apex=float(y[idx]); rel_apex=apex/base; rel_area=area/total_area if total_area>0 else 0
        left_mz=float(x[lb]); right_mz=float(x[rb]); width=right_mz-left_mz
        fwhm=_interp(x,rip)-_interp(x,lip); sharp=prom/max(fwhm,spacing)
        interim.append(T1ProfilePeak(t1_peak_id=f"SCIEX_T1_P{serial:05d}",source_id=source_id,
          measurement_id=measurement_id,rna_identity=rna_identity,apex_mz=float(x[idx]),centroid_mz=centroid,
          apex_intensity=apex,integrated_intensity=area,relative_apex_intensity=rel_apex,
          relative_integrated_intensity=rel_area,left_boundary_mz=left_mz,right_boundary_mz=right_mz,
          peak_width_mz=width,fwhm_mz=fwhm,prominence=float(prom),relative_prominence=float(prom/base),
          sharpness_score=sharp,nearest_peak_separation_mz=None,shared_valley_fraction=None,
          possible_shoulder=False,possible_duplicate=False,possible_isotope_component=False,
          possible_isotope_charge=None,possible_isotope_spacing_error=None,
          possible_overlapping_envelope=False,peak_quality_class=T1PeakQualityClass.LOW_SUPPORT))
    ordered=sorted(range(len(interim)),key=lambda i:interim[i].apex_mz)
    isotope=_isotope_annotations(interim,params.isotope_spacing_tolerance_mz)
    output=[]
    for pos,i in enumerate(ordered):
        peak=interim[i]; neighbors=[]
        if pos: neighbors.append(ordered[pos-1])
        if pos+1<len(ordered): neighbors.append(ordered[pos+1])
        nearest=min((abs(interim[j].apex_mz-peak.apex_mz) for j in neighbors),default=None)
        valley_fraction=None; shoulder=False; duplicate=False
        for j in neighbors:
            lo,hi=sorted((peak.apex_mz,interim[j].apex_mz)); mask=(x>=lo)&(x<=hi)
            valley=float(y[mask].min()) if mask.any() else 0
            frac=valley/max(min(peak.apex_intensity,interim[j].apex_intensity),1e-30)
            valley_fraction=max(valley_fraction or 0,frac)
            separation=hi-lo
            shoulder |= frac>=params.shared_valley_fraction and separation<max(peak.peak_width_mz,interim[j].peak_width_mz)
            duplicate |= frac>=0.8 and separation<params.duplicate_separation_mz
        iso=min(isotope[i]) if isotope[i] else None
        overlap=shoulder or bool(iso)
        sharp=peak.fwhm_mz<=params.sharp_fwhm_max_mz and peak.peak_width_mz<=params.maximum_peak_width_mz
        major=peak.relative_apex_intensity>=params.major_relative_intensity or peak.relative_integrated_intensity>=params.major_relative_intensity
        supported=peak.relative_apex_intensity>=params.minimum_relative_apex_intensity and peak.relative_integrated_intensity>=params.minimum_relative_integrated_intensity
        if duplicate or shoulder: quality=T1PeakQualityClass.SHOULDER_OR_OVERLAP
        elif iso: quality=T1PeakQualityClass.ISOTOPE_OR_ENVELOPE_COMPONENT
        elif major: quality=T1PeakQualityClass.MAJOR_SHARP if sharp else T1PeakQualityClass.MAJOR_BROAD
        elif supported: quality=T1PeakQualityClass.MINOR_SHARP if sharp else T1PeakQualityClass.MINOR_BROAD
        else: quality=T1PeakQualityClass.LOW_SUPPORT
        output.append(replace(peak,nearest_peak_separation_mz=nearest,shared_valley_fraction=valley_fraction,
          possible_shoulder=shoulder,possible_duplicate=duplicate,possible_isotope_component=bool(iso),
          possible_isotope_charge=iso[1] if iso else None,possible_isotope_spacing_error=iso[2] if iso else None,
          possible_overlapping_envelope=overlap,peak_quality_class=quality))
    output=sorted(output,key=lambda p:(p.apex_mz,p.t1_peak_id))
    eligible=[p for p in output if not p.possible_shoulder and not p.possible_duplicate and
      (p.peak_quality_class in {T1PeakQualityClass.MAJOR_SHARP,T1PeakQualityClass.MAJOR_BROAD} or
       (p.peak_quality_class is T1PeakQualityClass.MINOR_SHARP and p.relative_prominence>=params.minor_sharp_minimum_relative_prominence))]
    eligible.sort(key=lambda p:(p.possible_isotope_component,_QUALITY_ORDER[p.peak_quality_class],-p.relative_apex_intensity,-p.relative_integrated_intensity,p.apex_mz,p.t1_peak_id))
    selected_ids={p.t1_peak_id for p in eligible[:params.maximum_selected_peaks_per_profile]}
    output=tuple(replace(p,selected_as_primary=p.t1_peak_id in selected_ids) for p in output)
    return T1ProfilePeakAuditResult("COMPLETED",source_id,measurement_id,rna_identity,"MZ","MZ","CHARGED_ION_UNKNOWN",params,output,tuple(p for p in output if p.selected_as_primary))

def analyze_loaded_t1_profile(loaded, source, *, parameters=None):
    if loaded.profile_source_id!=source.profile_source_id: raise ValueError("loaded/source conflict")
    if tuple(loaded.header)!=("Mass/Charge","Intensity"): raise ValueError("T1 profile header must be Mass/Charge, Intensity")
    return detect_t1_profile_peaks(loaded.coordinates,loaded.intensities,source_id=source.profile_source_id,
      measurement_id=source.measurement_id,rna_identity=source.rna_identity_id,parameters=parameters)
