"""One-to-one shadow comparison for RNase-T1 m/z-profile peaks."""
from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from math import isfinite
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr
from rna_masshunter.sciex_t1_profile_peak_audit import T1ProfilePeak,T1ProfilePeakAuditResult,T1Safeguards

class T1ComparisonLayer(str,Enum):
 ALL_DETECTED_T1_PEAKS="ALL_DETECTED_T1_PEAKS"; SELECTED_T1_PEAKS="SELECTED_T1_PEAKS"; NON_ISOTOPE_SELECTED_T1_PEAKS="NON_ISOTOPE_SELECTED_T1_PEAKS"
class T1MassMatchClass(str,Enum): STRICT="STRICT"; EXPLORATORY="EXPLORATORY"
class T1ShapeSimilarityClass(str,Enum):
 HIGHLY_SIMILAR_T1_PEAK="HIGHLY_SIMILAR_T1_PEAK"; MODERATELY_SIMILAR_T1_PEAK="MODERATELY_SIMILAR_T1_PEAK"
 MZ_MATCH_SHAPE_DIFFERENT="MZ_MATCH_SHAPE_DIFFERENT"; MZ_ONLY_MATCH="MZ_ONLY_MATCH"; ISOTOPE_OR_ENVELOPE_AMBIGUOUS="ISOTOPE_OR_ENVELOPE_AMBIGUOUS"
class T1SelectedClassification(str,Enum):
 COMMON_T1_SELECTED_PEAK="COMMON_T1_SELECTED_PEAK"; UAA_ONLY_T1_SELECTED_PEAK="UAA_ONLY_T1_SELECTED_PEAK"; UAG_ONLY_T1_SELECTED_PEAK="UAG_ONLY_T1_SELECTED_PEAK"
 AMBIGUOUS_T1_CROSS_PROFILE_MATCH="AMBIGUOUS_T1_CROSS_PROFILE_MATCH"; ISOTOPE_OR_ENVELOPE_ONLY_MATCH="ISOTOPE_OR_ENVELOPE_ONLY_MATCH"
class CommonalityConcordanceClass(str,Enum):
 COMMONALITY_HIGH_IN_BOTH="COMMONALITY_HIGH_IN_BOTH"; COMMONALITY_HIGH_ONLY_IN_FULL="COMMONALITY_HIGH_ONLY_IN_FULL"; COMMONALITY_HIGH_ONLY_IN_T1="COMMONALITY_HIGH_ONLY_IN_T1"
 SAMPLE_SPECIFICITY_HIGH_IN_T1="SAMPLE_SPECIFICITY_HIGH_IN_T1"; UNRESOLVED="UNRESOLVED"
@dataclass(frozen=True)
class T1CrossProfileParameters:
 strict_cross_profile_tolerance_mz:float=.01; exploratory_cross_profile_tolerance_mz:float=.02; ambiguity_margin_mz:float=.003
 highly_similar_log_ratio_max:float=.3; moderately_similar_log_ratio_max:float=.7; minimum_shape_fields:int=4
 def validate(self):
  if self.strict_cross_profile_tolerance_mz<=0 or self.exploratory_cross_profile_tolerance_mz<self.strict_cross_profile_tolerance_mz: raise ValueError("invalid tolerance")
@dataclass(frozen=True,kw_only=True)
class T1CrossProfileMatch(T1Safeguards):
 t1_cross_profile_match_id:str; comparison_layer:T1ComparisonLayer; uaa_peak_id:str; uag_peak_id:str
 uaa_apex_mz:float; uag_apex_mz:float; apex_mz_difference:float; uaa_centroid_mz:float|None; uag_centroid_mz:float|None; centroid_mz_difference:float|None
 uaa_relative_apex_intensity:float; uag_relative_apex_intensity:float; uaa_relative_integrated_intensity:float; uag_relative_integrated_intensity:float
 uaa_fwhm:float; uag_fwhm:float; uaa_peak_width:float; uag_peak_width:float; uaa_prominence:float; uag_prominence:float
 uaa_quality_class:str; uag_quality_class:str; uaa_isotope_flag:bool; uag_isotope_flag:bool
 mass_match_class:T1MassMatchClass; shape_similarity_class:T1ShapeSimilarityClass; ambiguous_assignment:bool
@dataclass(frozen=True)
class T1SelectedPeakStatus:
 profile:str; peak_id:str; apex_mz:float; classification:T1SelectedClassification; matched_peak_id:str|None
@dataclass(frozen=True,kw_only=True)
class T1ProfileCommonSummary(T1Safeguards):
 profile:str; selected_peak_count:int; common_selected_count:int; sample_specific_selected_count:int; isotope_only_selected_count:int
 common_apex_intensity_fraction:float; sample_specific_apex_intensity_fraction:float; common_integrated_intensity_fraction:float; sample_specific_integrated_intensity_fraction:float
 hypotheses:tuple[str,...]
@dataclass(frozen=True)
class T1CrossCorrelations:
 method:str; apex_mz:float|None; relative_apex_intensity:float|None; relative_integrated_intensity:float|None; prominence:float|None; fwhm:float|None
@dataclass(frozen=True)
class FullT1CommonalitySummary:
 full_profile_common_component_fraction:float; t1_profile_common_component_fraction:float; full_profile_sample_specific_fraction:float; t1_profile_sample_specific_fraction:float
 full_and_t1_commonality_concordance_class:CommonalityConcordanceClass; direct_mass_connection_performed:bool=False
@dataclass(frozen=True)
class T1CrossProfileAuditResult:
 parameters:T1CrossProfileParameters; all_detected_matches:tuple[T1CrossProfileMatch,...]; selected_matches:tuple[T1CrossProfileMatch,...]
 non_isotope_selected_matches:tuple[T1CrossProfileMatch,...]; statuses:tuple[T1SelectedPeakStatus,...]; uaa_summary:T1ProfileCommonSummary; uag_summary:T1ProfileCommonSummary; correlations:T1CrossCorrelations

def _ordered(values): return tuple(sorted(values,key=lambda p:(p.apex_mz,p.t1_peak_id)))
def _ratio(a,b): return b/a if a and b and a>0 and b>0 else None
def _shape(a,b):
 if a.possible_isotope_component or b.possible_isotope_component:return T1ShapeSimilarityClass.ISOTOPE_OR_ENVELOPE_AMBIGUOUS
 vals=[_ratio(a.fwhm_mz,b.fwhm_mz),_ratio(a.peak_width_mz,b.peak_width_mz),_ratio(a.prominence,b.prominence),_ratio(a.relative_apex_intensity,b.relative_apex_intensity),_ratio(a.relative_integrated_intensity,b.relative_integrated_intensity)]
 ds=[abs(np.log10(v)) for v in vals if v is not None]; cd=abs(a.centroid_mz-b.centroid_mz) if a.centroid_mz is not None and b.centroid_mz is not None else None
 if len(ds)<4 or cd is None:return T1ShapeSimilarityClass.MZ_ONLY_MATCH
 if a.peak_quality_class is b.peak_quality_class and cd<=.01 and max(ds)<=.3:return T1ShapeSimilarityClass.HIGHLY_SIMILAR_T1_PEAK
 if cd<=.02 and sum(d<=.7 for d in ds)>=4:return T1ShapeSimilarityClass.MODERATELY_SIMILAR_T1_PEAK
 return T1ShapeSimilarityClass.MZ_MATCH_SHAPE_DIFFERENT

def _components(a,b,tol):
 nodes=sorted([(p.apex_mz,0,i) for i,p in enumerate(a)]+[(p.apex_mz,1,j) for j,p in enumerate(b)])
 out=[]; current=[]; last=None
 for node in nodes:
  if last is not None and node[0]-last>tol:
   out.append(current);current=[]
  current.append(node);last=node[0]
 if current:out.append(current)
 return out

def _assignment(a,b,tol):
 pairs=[]
 for component in _components(a,b,tol):
  ai=sorted(n[2] for n in component if n[1]==0); bj=sorted(n[2] for n in component if n[1]==1)
  if not ai or not bj:continue
  n,m=len(ai),len(bj);size=n+m;unmatched=tol+1;forbid=unmatched*(size+2);cost=np.full((size,size),forbid)
  for ii,i in enumerate(ai):
   for jj,j in enumerate(bj):
    e=abs(a[i].apex_mz-b[j].apex_mz)
    if e<=tol:cost[ii,jj]=e+(ii*(m+1)+jj)*1e-12
   cost[ii,m+ii]=unmatched
  for jj in range(m):cost[n+jj,jj]=unmatched
  cost[n:,m:]=0
  rows,cols=linear_sum_assignment(cost)
  pairs.extend((ai[ii],bj[jj]) for ii,jj in zip(rows,cols) if ii<n and jj<m and cost[ii,jj]<unmatched)
 return tuple(sorted(pairs,key=lambda q:(a[q[0]].apex_mz,a[q[0]].t1_peak_id,b[q[1]].t1_peak_id)))

def match_t1_cross_profile(uaa,uag,*,layer,parameters=None):
 p=parameters or T1CrossProfileParameters();p.validate();a=_ordered(uaa);b=_ordered(uag);out=[]
 for i,j in _assignment(a,b,p.exploratory_cross_profile_tolerance_mz):
  x,y=a[i],b[j];err=abs(x.apex_mz-y.apex_mz);cd=abs(x.centroid_mz-y.centroid_mz) if x.centroid_mz is not None and y.centroid_mz is not None else None
  alternatives=[abs(x.apex_mz-q.apex_mz) for k,q in enumerate(b) if k!=j and abs(x.apex_mz-q.apex_mz)<=p.exploratory_cross_profile_tolerance_mz]
  ambiguous=bool(alternatives and min(alternatives)<=err+p.ambiguity_margin_mz)
  out.append(T1CrossProfileMatch(t1_cross_profile_match_id=f"T1CROSS__{layer.value}__{x.t1_peak_id}__{y.t1_peak_id}",comparison_layer=layer,
   uaa_peak_id=x.t1_peak_id,uag_peak_id=y.t1_peak_id,uaa_apex_mz=x.apex_mz,uag_apex_mz=y.apex_mz,apex_mz_difference=err,
   uaa_centroid_mz=x.centroid_mz,uag_centroid_mz=y.centroid_mz,centroid_mz_difference=cd,uaa_relative_apex_intensity=x.relative_apex_intensity,uag_relative_apex_intensity=y.relative_apex_intensity,
   uaa_relative_integrated_intensity=x.relative_integrated_intensity,uag_relative_integrated_intensity=y.relative_integrated_intensity,uaa_fwhm=x.fwhm_mz,uag_fwhm=y.fwhm_mz,
   uaa_peak_width=x.peak_width_mz,uag_peak_width=y.peak_width_mz,uaa_prominence=x.prominence,uag_prominence=y.prominence,uaa_quality_class=x.peak_quality_class.value,uag_quality_class=y.peak_quality_class.value,
   uaa_isotope_flag=x.possible_isotope_component,uag_isotope_flag=y.possible_isotope_component,mass_match_class=T1MassMatchClass.STRICT if err<=p.strict_cross_profile_tolerance_mz+1e-12 else T1MassMatchClass.EXPLORATORY,
   shape_similarity_class=_shape(x,y),ambiguous_assignment=ambiguous))
 return tuple(out)
def _fraction(values,flags):
 total=sum(values);common=sum(v for v,f in zip(values,flags) if f);return (common/total,(total-common)/total) if total else (0.,0.)
def _summary(profile,peaks,statuses):
 own=[s for s in statuses if s.profile==profile]; common=[s.classification in {T1SelectedClassification.COMMON_T1_SELECTED_PEAK,T1SelectedClassification.AMBIGUOUS_T1_CROSS_PROFILE_MATCH,T1SelectedClassification.ISOTOPE_OR_ENVELOPE_ONLY_MATCH} for s in own]
 ca,sa=_fraction([p.relative_apex_intensity for p in peaks],common);ci,si=_fraction([p.relative_integrated_intensity for p in peaks],common)
 return T1ProfileCommonSummary(profile=profile,selected_peak_count=len(peaks),common_selected_count=sum(common),sample_specific_selected_count=len(peaks)-sum(common),isotope_only_selected_count=sum(s.classification is T1SelectedClassification.ISOTOPE_OR_ENVELOPE_ONLY_MATCH for s in own),common_apex_intensity_fraction=ca,sample_specific_apex_intensity_fraction=sa,common_integrated_intensity_fraction=ci,sample_specific_integrated_intensity_fraction=si,hypotheses=("COMMON_T1_COMPONENT_POSSIBLE","COMMON_CO_CAPTURED_RNA_FRAGMENT_POSSIBLE","COMMON_BACKGROUND_OR_CONTAMINANT_POSSIBLE","COMMON_IN_SOURCE_FRAGMENT_POSSIBLE","TARGET_SPECIFIC_FRAGMENT_POSSIBLE","SAMPLE_SPECIFIC_CONTAMINANT_POSSIBLE","LOW_ABUNDANCE_SHARED_FRAGMENT_POSSIBLE","IONIZATION_OR_DETECTION_DIFFERENCE_POSSIBLE"))
def _corr(a,b):
 pairs=[(x,y) for x,y in zip(a,b) if x is not None and y is not None]
 if len(pairs)<3 or len({x for x,y in pairs})<2 or len({y for x,y in pairs})<2:return None
 v=float(spearmanr([x for x,y in pairs],[y for x,y in pairs]).statistic);return v if isfinite(v) else None
def audit_t1_cross_profiles(uaa:T1ProfilePeakAuditResult,uag:T1ProfilePeakAuditResult,*,parameters=None):
 p=parameters or T1CrossProfileParameters();a=uaa.selected_peaks;b=uag.selected_peaks
 allm=match_t1_cross_profile(uaa.peaks,uag.peaks,layer=T1ComparisonLayer.ALL_DETECTED_T1_PEAKS,parameters=p)
 sm=match_t1_cross_profile(a,b,layer=T1ComparisonLayer.SELECTED_T1_PEAKS,parameters=p)
 non=match_t1_cross_profile([x for x in a if not x.possible_isotope_component],[x for x in b if not x.possible_isotope_component],layer=T1ComparisonLayer.NON_ISOTOPE_SELECTED_T1_PEAKS,parameters=p)
 am={x.uaa_peak_id:x for x in sm};bm={x.uag_peak_id:x for x in sm};statuses=[]
 for profile,peaks,mapping,other in (("UAA",a,am,"uag_peak_id"),("UAG",b,bm,"uaa_peak_id")):
  for peak in peaks:
   m=mapping.get(peak.t1_peak_id)
   cls=(T1SelectedClassification.AMBIGUOUS_T1_CROSS_PROFILE_MATCH if m and m.ambiguous_assignment else T1SelectedClassification.ISOTOPE_OR_ENVELOPE_ONLY_MATCH if m and (m.uaa_isotope_flag or m.uag_isotope_flag) else T1SelectedClassification.COMMON_T1_SELECTED_PEAK if m else T1SelectedClassification.UAA_ONLY_T1_SELECTED_PEAK if profile=="UAA" else T1SelectedClassification.UAG_ONLY_T1_SELECTED_PEAK)
   statuses.append(T1SelectedPeakStatus(profile,peak.t1_peak_id,peak.apex_mz,cls,getattr(m,other) if m else None))
 statuses=tuple(statuses);corr=T1CrossCorrelations("SPEARMAN_RANK_CORRELATION",_corr([m.uaa_apex_mz for m in sm],[m.uag_apex_mz for m in sm]),_corr([m.uaa_relative_apex_intensity for m in sm],[m.uag_relative_apex_intensity for m in sm]),_corr([m.uaa_relative_integrated_intensity for m in sm],[m.uag_relative_integrated_intensity for m in sm]),_corr([m.uaa_prominence for m in sm],[m.uag_prominence for m in sm]),_corr([m.uaa_fwhm for m in sm],[m.uag_fwhm for m in sm]))
 return T1CrossProfileAuditResult(p,allm,sm,non,statuses,_summary("UAA",a,statuses),_summary("UAG",b,statuses),corr)
def connect_full_t1_commonality(full_cross,t1_cross):
 full=(full_cross.uaa_summary.common_selected_apex_intensity_fraction+full_cross.uag_summary.common_selected_apex_intensity_fraction)/2
 t1=(t1_cross.uaa_summary.common_apex_intensity_fraction+t1_cross.uag_summary.common_apex_intensity_fraction)/2
 cls=CommonalityConcordanceClass.COMMONALITY_HIGH_IN_BOTH if full>=.5 and t1>=.5 else CommonalityConcordanceClass.COMMONALITY_HIGH_ONLY_IN_FULL if full>=.5 else CommonalityConcordanceClass.COMMONALITY_HIGH_ONLY_IN_T1 if t1>=.5 else CommonalityConcordanceClass.SAMPLE_SPECIFICITY_HIGH_IN_T1 if t1<.3 else CommonalityConcordanceClass.UNRESOLVED
 return FullT1CommonalitySummary(full,t1,1-full,1-t1,cls)
