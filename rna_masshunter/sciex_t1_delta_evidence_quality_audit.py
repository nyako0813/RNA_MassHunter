"""Independent evidence-quality shadow audit for T1 chemical-delta relations."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from hashlib import sha256
from typing import Iterable,Mapping,Sequence
from rna_masshunter.sciex_intact_peak_family import DeltaMassDefinition
from rna_masshunter.sciex_t1_profile_peak_audit import T1ProfilePeak,T1PeakQualityClass
from rna_masshunter.sciex_t1_fragment_shadow_match import T1IonMode,TheoreticalDiscriminationClass
from rna_masshunter.sciex_t1_fragment_delta_audit import (
 ChemicalDeltaAuditParameters,ChemicalDeltaReference,ChemicalHypothesisClass,
 ChemicalReferenceCategory,DeltaRelationClass,DeltaSafeguards,ReferenceAmbiguityClass,
 T1ChemicalDeltaRelation,T1ChemicalStateSeries,T1FragmentDeltaAuditResult,T1StateSeriesPattern,
 _ref_matches,
)

class ObservedPolarityStatus(str,Enum):
 POLARITY_CONFIRMED_NEGATIVE="POLARITY_CONFIRMED_NEGATIVE";POLARITY_CONFIRMED_POSITIVE="POLARITY_CONFIRMED_POSITIVE";POLARITY_UNKNOWN="POLARITY_UNKNOWN";POLARITY_CONFLICT="POLARITY_CONFLICT"
class PolaritySupportClass(str,Enum):
 CONFIRMED_COMPATIBLE="CONFIRMED_COMPATIBLE";CONFIRMED_INCOMPATIBLE="CONFIRMED_INCOMPATIBLE";UNKNOWN_POLARITY="UNKNOWN_POLARITY";CONFLICTING_METADATA="CONFLICTING_METADATA"
class ChargeSupportClass(str,Enum):
 MULTI_CHARGE_CONCORDANT_SUPPORT="MULTI_CHARGE_CONCORDANT_SUPPORT";SINGLE_CHARGE_ONLY="SINGLE_CHARGE_ONLY";MULTI_CHARGE_CONFLICTING="MULTI_CHARGE_CONFLICTING";NO_CHARGE_SUPPORT="NO_CHARGE_SUPPORT"
class ApexCentroidConcordanceClass(str,Enum):
 BOTH_STRICT_SAME_REFERENCE="BOTH_STRICT_SAME_REFERENCE";APEX_STRICT_CENTROID_EXPLORATORY_SAME_REFERENCE="APEX_STRICT_CENTROID_EXPLORATORY_SAME_REFERENCE";BOTH_EXPLORATORY_SAME_REFERENCE="BOTH_EXPLORATORY_SAME_REFERENCE";APEX_ONLY_SUPPORT="APEX_ONLY_SUPPORT";CENTROID_ONLY_SUPPORT="CENTROID_ONLY_SUPPORT";REFERENCE_CATEGORY_CONFLICT="REFERENCE_CATEGORY_CONFLICT";NO_CONCORDANT_SUPPORT="NO_CONCORDANT_SUPPORT"
class PeakShapeSupportClass(str,Enum):
 STRONG_PEAK_SHAPE_SUPPORT="STRONG_PEAK_SHAPE_SUPPORT";MODERATE_PEAK_SHAPE_SUPPORT="MODERATE_PEAK_SHAPE_SUPPORT";WEAK_PEAK_SHAPE_SUPPORT="WEAK_PEAK_SHAPE_SUPPORT";EXCLUDED_PEAK_SHAPE="EXCLUDED_PEAK_SHAPE"
class FragmentAmbiguityClass(str,Enum):
 LOW_FRAGMENT_AMBIGUITY="LOW_FRAGMENT_AMBIGUITY";MODERATE_FRAGMENT_AMBIGUITY="MODERATE_FRAGMENT_AMBIGUITY";HIGH_FRAGMENT_AMBIGUITY="HIGH_FRAGMENT_AMBIGUITY";CROSS_RNA_AMBIGUITY="CROSS_RNA_AMBIGUITY"
class DiscriminatorySupportClass(str,Enum):
 TARGET_SPECIFIC_FRAGMENT_SUPPORT="TARGET_SPECIFIC_FRAGMENT_SUPPORT";SHARED_FRAGMENT_SUPPORT="SHARED_FRAGMENT_SUPPORT";SAME_MASS_SEQUENCE_AMBIGUITY="SAME_MASS_SEQUENCE_AMBIGUITY";CROSS_RNA_AMBIGUITY="CROSS_RNA_AMBIGUITY";NO_DISCRIMINATORY_INFORMATION="NO_DISCRIMINATORY_INFORMATION"
class ReferenceResolutionClass(str,Enum):
 WELL_RESOLVED_REFERENCE_CATEGORY="WELL_RESOLVED_REFERENCE_CATEGORY";MARGINALLY_RESOLVED_REFERENCE_CATEGORY="MARGINALLY_RESOLVED_REFERENCE_CATEGORY";UNRESOLVED_REFERENCE_CATEGORY="UNRESOLVED_REFERENCE_CATEGORY";SINGLE_REFERENCE_ONLY="SINGLE_REFERENCE_ONLY"
class MassDefinitionCompatibilityStatus(str,Enum):
 CONFIRMED_COMPATIBLE="CONFIRMED_COMPATIBLE";UNKNOWN_COMPATIBILITY="UNKNOWN_COMPATIBILITY";CONFIRMED_MISMATCH="CONFIRMED_MISMATCH"
class RecurrentSupportClass(str,Enum):
 MULTI_PEAK_MULTI_CHARGE_RECURRENT="MULTI_PEAK_MULTI_CHARGE_RECURRENT";MULTI_PEAK_SINGLE_CHARGE_RECURRENT="MULTI_PEAK_SINGLE_CHARGE_RECURRENT";SINGLE_PEAK_ONLY="SINGLE_PEAK_ONLY";CROSS_SAMPLE_RECURRENT="CROSS_SAMPLE_RECURRENT";NO_RECURRENT_SUPPORT="NO_RECURRENT_SUPPORT"
class CrossSampleSupportStatus(str,Enum):
 SUPPORTED_IN_BOTH_SAMPLES="SUPPORTED_IN_BOTH_SAMPLES";SUPPORTED_ONLY_IN_UAA="SUPPORTED_ONLY_IN_UAA";SUPPORTED_ONLY_IN_UAG="SUPPORTED_ONLY_IN_UAG";CROSS_RNA_SEQUENCE_AMBIGUOUS="CROSS_RNA_SEQUENCE_AMBIGUOUS";NOT_COMPARABLE="NOT_COMPARABLE"
class StateSeriesSupportClass(str,Enum):
 STRONG_LINEAR_SERIES_SUPPORT="STRONG_LINEAR_SERIES_SUPPORT";MODERATE_SERIES_SUPPORT="MODERATE_SERIES_SUPPORT";BRANCHED_AMBIGUOUS_SERIES="BRANCHED_AMBIGUOUS_SERIES";SINGLE_EDGE_ONLY="SINGLE_EDGE_ONLY";NO_SERIES_SUPPORT="NO_SERIES_SUPPORT"
class EvidenceTier(str,Enum):
 TIER_A_HIGH_SUPPORT="TIER_A_HIGH_SUPPORT";TIER_B_MODERATE_SUPPORT="TIER_B_MODERATE_SUPPORT";TIER_C_WEAK_SUPPORT="TIER_C_WEAK_SUPPORT";TIER_D_DIAGNOSTIC_ONLY="TIER_D_DIAGNOSTIC_ONLY";TIER_E_BLOCKED="TIER_E_BLOCKED"
class EvidenceBlockReason(str,Enum):
 INVALID_RELATION="INVALID_RELATION";CONFIRMED_POLARITY_MISMATCH="CONFIRMED_POLARITY_MISMATCH";ISOTOPE_OR_ENVELOPE_COMPONENT="ISOTOPE_OR_ENVELOPE_COMPONENT";SHOULDER_OR_DUPLICATE="SHOULDER_OR_DUPLICATE";MASS_DEFINITION_MISMATCH="MASS_DEFINITION_MISMATCH";CROSS_RNA_AMBIGUITY="CROSS_RNA_AMBIGUITY";MULTIPLE_FRAGMENT_AMBIGUITY="MULTIPLE_FRAGMENT_AMBIGUITY";UNRESOLVED_REFERENCE_CATEGORY="UNRESOLVED_REFERENCE_CATEGORY";APEX_CENTROID_DISCORDANCE="APEX_CENTROID_DISCORDANCE";SINGLE_CHARGE_ONLY="SINGLE_CHARGE_ONLY";UNKNOWN_POLARITY="UNKNOWN_POLARITY";UNKNOWN_MASS_DEFINITION="UNKNOWN_MASS_DEFINITION";NO_RECURRENT_SUPPORT="NO_RECURRENT_SUPPORT";NO_BLOCK="NO_BLOCK"
class PriorityOSClass(str,Enum):
 O_TO_S_PRIORITY_CANDIDATE="O_TO_S_PRIORITY_CANDIDATE";S_TO_O_PRIORITY_CANDIDATE="S_TO_O_PRIORITY_CANDIDATE";O_TO_S_DIAGNOSTIC_ONLY="O_TO_S_DIAGNOSTIC_ONLY";S_TO_O_DIAGNOSTIC_ONLY="S_TO_O_DIAGNOSTIC_ONLY";BLOCKED_O_S_CANDIDATE="BLOCKED_O_S_CANDIDATE"

@dataclass(frozen=True)
class EvidenceQualityParameters:
 well_resolved_margin_da:float=.05;marginally_resolved_margin_da:float=.02
 both_strict_score:int=20;strong_shape_score:int=15;low_ambiguity_score:int=15;resolved_reference_score:int=10;multi_charge_score:int=10;recurrent_score:int=10;target_specific_score:int=10;strong_series_score:int=10
 artifact_penalty:int=25;cross_rna_penalty:int=20;unresolved_reference_penalty:int=20;unknown_polarity_penalty:int=15;unknown_mass_definition_penalty:int=15;branched_series_penalty:int=15;single_charge_penalty:int=10
 tier_b_minimum:int=60;tier_c_minimum:int=25
 def validate(self):
  if not 0<=self.marginally_resolved_margin_da<=self.well_resolved_margin_da:raise ValueError("invalid resolution margins")
  if not 0<=self.tier_c_minimum<=self.tier_b_minimum<=100:raise ValueError("invalid tier thresholds")

@dataclass(frozen=True,kw_only=True)
class QualitySafeguards(DeltaSafeguards):
 evidence_tier_formal:bool=False;evidence_support_score_formal:bool=False
 in_source_fragmentation_excluded:bool=False;different_fragment_explanation_excluded:bool=False

@dataclass(frozen=True,kw_only=True)
class T1DeltaEvidenceQualityRecord(QualitySafeguards):
 evidence_quality_record_id:str;source_id:str;measurement_id:str;rna_identity_candidate:str;t1_delta_relation_id:str;observed_peak_id:str;theoretical_fragment_id:str;fragment_sequence:str;start_position:int;end_position:int;cca_state:str;ion_mode:T1IonMode;charge:int
 chemical_hypothesis_class:ChemicalHypothesisClass;reference_id:str;reference_name:str;reference_mass_definition:DeltaMassDefinition;observed_apex_mz:float;observed_centroid_mz:float|None;theoretical_mz:float;apex_neutral_delta:float;centroid_neutral_delta:float|None;reference_delta_da:float;apex_error_da:float|None;centroid_error_da:float|None
 observed_polarity_status:ObservedPolarityStatus;ion_mode_polarity_compatible:bool|None;polarity_support_class:PolaritySupportClass
 charge_support_count:int;supported_charge_states:tuple[int,...];distinct_observed_peak_count:int;charge_support_class:ChargeSupportClass
 apex_reference_match_class:str;centroid_reference_match_class:str;apex_reference_category:ChemicalReferenceCategory|None;centroid_reference_category:ChemicalReferenceCategory|None;apex_centroid_same_reference:bool;apex_centroid_error_difference:float|None;apex_centroid_concordance_class:ApexCentroidConcordanceClass
 peak_quality_class:T1PeakQualityClass;fwhm:float;prominence:float;sharpness:float;possible_overlapping_envelope:bool;peak_shape_support_class:PeakShapeSupportClass
 distinct_fragment_count:int;distinct_charge_count:int;distinct_cca_state_count:int;distinct_rna_identity_count:int;candidate_reference_count:int;distinct_reference_category_count:int;raw_candidate_count:int;cca_collapsed_candidate_count:int;fragment_identity_collapsed_count:int;fragment_ambiguity_class:FragmentAmbiguityClass
 fragment_discrimination_class:TheoreticalDiscriminationClass|None;uaa_specific_theoretical:bool;uag_specific_theoretical:bool;shared_theoretical:bool;same_mass_different_sequence:bool;cross_rna_identity_ambiguous:bool;discriminatory_support_class:DiscriminatorySupportClass
 reference_category_count:int;reference_category_list:tuple[str,...];best_category_error:float|None;second_best_category_error:float|None;reference_resolution_margin_da:float|None;reference_resolution_class:ReferenceResolutionClass
 observed_mass_definition_status:str;mass_definition_compatibility_status:MassDefinitionCompatibilityStatus
 recurrent_support_group_id:str;recurrent_observed_peak_count:int;recurrent_fragment_identity_count:int;recurrent_charge_count:int;recurrent_source_count:int;recurrent_support_class:RecurrentSupportClass
 cross_sample_support_status:CrossSampleSupportStatus;supported_in_uaa:bool;supported_in_uag:bool;cross_sample_fragment_identity_match:bool;cross_sample_chemical_category_match:bool;common_component_or_shared_fragment_possible:bool
 member_of_state_series:bool;state_series_id:str|None;state_series_member_count:int;state_series_strict_edge_count:int;state_series_pattern:T1StateSeriesPattern|None;state_series_support_class:StateSeriesSupportClass
 alternative_explanation_count:int;alternative_explanation_categories:tuple[str,...];adduct_alternative_possible:bool;water_loss_alternative_possible:bool;oxidation_alternative_possible:bool;known_modification_alternative_possible:bool;different_fragment_alternative_possible:bool
 expected_mz_for_alternative_charge:float|None;observed_supporting_peak_id:str|None;alternative_charge_error_mz:float|None;alternative_charge_support_class:str
 evidence_tier:EvidenceTier;evidence_support_score:int;evidence_block_reason:EvidenceBlockReason;evidence_warning_flags:tuple[str,...];eligible_for_high_quality_shadow_evidence:bool

@dataclass(frozen=True,kw_only=True)
class OSPriorityQualityRecord(QualitySafeguards):
 quality_record:T1DeltaEvidenceQualityRecord;priority_class:PriorityOSClass
@dataclass(frozen=True,kw_only=True)
class T1DeltaEvidenceQualityAuditResult(QualitySafeguards):
 source_id:str;records:tuple[T1DeltaEvidenceQualityRecord,...];os_priority_records:tuple[OSPriorityQualityRecord,...];input_delta_relation_count:int;quality_record_count:int;oxygen_equivalent_pattern_also_observed_in_glu_intact:bool=False;water_equivalent_pattern_also_observed_in_glu_intact:bool=False;direct_chemical_identity_link_assigned:bool=False

_MATCHED={"STRICT","EXPLORATORY"}
_GOOD_CONCORD={ApexCentroidConcordanceClass.BOTH_STRICT_SAME_REFERENCE,ApexCentroidConcordanceClass.APEX_STRICT_CENTROID_EXPLORATORY_SAME_REFERENCE}

def determine_observed_polarity(metadata_values:Iterable[object]=()):
 values=set()
 for raw in metadata_values:
  text=str(raw or "").strip().lower()
  if text in {"negative","neg","-"}:values.add("negative")
  elif text in {"positive","pos","+"}:values.add("positive")
 if len(values)>1:return ObservedPolarityStatus.POLARITY_CONFLICT
 if values=={"negative"}:return ObservedPolarityStatus.POLARITY_CONFIRMED_NEGATIVE
 if values=={"positive"}:return ObservedPolarityStatus.POLARITY_CONFIRMED_POSITIVE
 return ObservedPolarityStatus.POLARITY_UNKNOWN

def _polarity(status,mode):
 if status is ObservedPolarityStatus.POLARITY_UNKNOWN:return None,PolaritySupportClass.UNKNOWN_POLARITY
 if status is ObservedPolarityStatus.POLARITY_CONFLICT:return None,PolaritySupportClass.CONFLICTING_METADATA
 compatible=(status is ObservedPolarityStatus.POLARITY_CONFIRMED_NEGATIVE and mode is T1IonMode.NEGATIVE_DEPROTONATED) or (status is ObservedPolarityStatus.POLARITY_CONFIRMED_POSITIVE and mode is T1IonMode.POSITIVE_PROTONATED)
 return compatible,PolaritySupportClass.CONFIRMED_COMPATIBLE if compatible else PolaritySupportClass.CONFIRMED_INCOMPATIBLE

def _best_match(delta,z,refs,delta_parameters):
 if delta is None:return None
 matches=_ref_matches(delta,z,refs,delta_parameters)
 return matches[0] if matches else None

def _concordance(r,refs,dp):
 a=_best_match(r.apex_neutral_delta,r.charge,refs,dp);c=_best_match(r.centroid_neutral_delta,r.charge,refs,dp)
 ac=a[1] if a else "NO_MATCH";cc=c[1] if c else "NO_MATCH";ar=a[2] if a else None;cr=c[2] if c else None;same=bool(ar and cr and ar.reference_id==cr.reference_id);diff=abs(abs(a[3])-abs(c[3])) if a and c else None
 if same and ac==cc=="STRICT":klass=ApexCentroidConcordanceClass.BOTH_STRICT_SAME_REFERENCE
 elif same and ac=="STRICT" and cc=="EXPLORATORY":klass=ApexCentroidConcordanceClass.APEX_STRICT_CENTROID_EXPLORATORY_SAME_REFERENCE
 elif same and ac==cc=="EXPLORATORY":klass=ApexCentroidConcordanceClass.BOTH_EXPLORATORY_SAME_REFERENCE
 elif a and not c:klass=ApexCentroidConcordanceClass.APEX_ONLY_SUPPORT
 elif c and not a:klass=ApexCentroidConcordanceClass.CENTROID_ONLY_SUPPORT
 elif a and c and ar.reference_category is not cr.reference_category:klass=ApexCentroidConcordanceClass.REFERENCE_CATEGORY_CONFLICT
 else:klass=ApexCentroidConcordanceClass.NO_CONCORDANT_SUPPORT
 return ac,cc,ar,cr,same,diff,klass

def _shape(p):
 if p.possible_isotope_component or p.possible_shoulder or p.possible_duplicate:return PeakShapeSupportClass.EXCLUDED_PEAK_SHAPE
 if p.peak_quality_class is T1PeakQualityClass.MAJOR_SHARP and not p.possible_overlapping_envelope:return PeakShapeSupportClass.STRONG_PEAK_SHAPE_SUPPORT
 if p.peak_quality_class is T1PeakQualityClass.MINOR_SHARP:return PeakShapeSupportClass.MODERATE_PEAK_SHAPE_SUPPORT
 return PeakShapeSupportClass.WEAK_PEAK_SHAPE_SUPPORT

def _resolution(r,refs,dp,qp):
 matches=_ref_matches(r.apex_neutral_delta,r.charge,refs,dp);best={}
 for err,_,ref,_ in matches:best[ref.reference_category]=min(best.get(ref.reference_category,float("inf")),err)
 ordered=sorted(best.items(),key=lambda x:(x[1],x[0].value));cats=tuple(x[0].value for x in ordered)
 if len(ordered)<=1:return cats,(ordered[0][1] if ordered else None),None,None,ReferenceResolutionClass.SINGLE_REFERENCE_ONLY
 margin=ordered[1][1]-ordered[0][1];klass=ReferenceResolutionClass.WELL_RESOLVED_REFERENCE_CATEGORY if margin>=qp.well_resolved_margin_da else ReferenceResolutionClass.MARGINALLY_RESOLVED_REFERENCE_CATEGORY if margin>=qp.marginally_resolved_margin_da else ReferenceResolutionClass.UNRESOLVED_REFERENCE_CATEGORY
 return cats,ordered[0][1],ordered[1][1],margin,klass

def _series_support(series):
 if series is None:return StateSeriesSupportClass.NO_SERIES_SUPPORT
 if series.series_pattern is T1StateSeriesPattern.BRANCHED_SERIES:return StateSeriesSupportClass.BRANCHED_AMBIGUOUS_SERIES
 edges=series.strict_edge_count+series.exploratory_edge_count
 if series.member_count>=3 and series.strict_edge_count>=max(2,series.exploratory_edge_count):return StateSeriesSupportClass.STRONG_LINEAR_SERIES_SUPPORT
 if edges==1:return StateSeriesSupportClass.SINGLE_EDGE_ONLY
 return StateSeriesSupportClass.MODERATE_SERIES_SUPPORT

def _discrimination(r,rows):
 candidates=[x for x in rows if x.rna_identity==r.rna_identity_candidate and x.fragment_sequence==r.fragment_sequence]
 classes={x.discrimination_class for x in candidates};cls=sorted(classes,key=lambda x:x.value)[0] if len(classes)==1 else None
 u=TheoreticalDiscriminationClass.UAA_SPECIFIC_THEORETICAL_FRAGMENT in classes;g=TheoreticalDiscriminationClass.UAG_SPECIFIC_THEORETICAL_FRAGMENT in classes;shared=TheoreticalDiscriminationClass.SEQUENCE_IDENTICAL_FRAGMENT in classes;same=TheoreticalDiscriminationClass.SAME_MASS_DIFFERENT_SEQUENCE in classes
 support=DiscriminatorySupportClass.TARGET_SPECIFIC_FRAGMENT_SUPPORT if u or g else DiscriminatorySupportClass.SAME_MASS_SEQUENCE_AMBIGUITY if same else DiscriminatorySupportClass.SHARED_FRAGMENT_SUPPORT if shared else DiscriminatorySupportClass.NO_DISCRIMINATORY_INFORMATION
 return cls,u,g,shared,same,support

def _sample_flags(rows,r):
 tail=r.theoretical_fragment_id.split(":",1)[-1];support={"UAA":False,"UAG":False};rnas=set()
 for x in rows:
  if x.reference_match_class not in _MATCHED or x.chemical_hypothesis_class is not r.chemical_hypothesis_class:continue
  if x.theoretical_fragment_id.split(":",1)[-1]!=tail:continue
  label="UAA" if "UAA" in x.observed_source_id else "UAG" if "UAG" in x.observed_source_id else None
  if label:support[label]=True
  rnas.add(x.rna_identity_candidate)
 ambiguous=len(rnas)>1
 if ambiguous:status=CrossSampleSupportStatus.CROSS_RNA_SEQUENCE_AMBIGUOUS
 elif all(support.values()):status=CrossSampleSupportStatus.SUPPORTED_IN_BOTH_SAMPLES
 elif support["UAA"]:status=CrossSampleSupportStatus.SUPPORTED_ONLY_IN_UAA
 elif support["UAG"]:status=CrossSampleSupportStatus.SUPPORTED_ONLY_IN_UAG
 else:status=CrossSampleSupportStatus.NOT_COMPARABLE
 return status,support["UAA"],support["UAG"],ambiguous

def _fragment_key(r):return (r.rna_identity_candidate,r.fragment_sequence,r.start_position,r.end_position,r.ion_mode,r.chemical_hypothesis_class)

def audit_t1_delta_evidence_quality(delta_result:T1FragmentDeltaAuditResult,peaks:Sequence[T1ProfilePeak],*,measurement_id:str,polarity_status:ObservedPolarityStatus=ObservedPolarityStatus.POLARITY_UNKNOWN,all_sample_relations:Sequence[T1ChemicalDeltaRelation]=(),cross_profile_matches:Sequence[object]=(),discrimination:Sequence[object]=(),parameters:EvidenceQualityParameters|None=None,glu_summary_flags:Mapping[str,bool]|None=None):
 qp=parameters or EvidenceQualityParameters();qp.validate();dp=delta_result.parameters;refs=delta_result.references;peak_map={p.t1_peak_id:p for p in peaks};relations=tuple(r for r in delta_result.relations if r.reference_match_class in _MATCHED);all_rows=tuple(all_sample_relations) or relations;cross_profile_peak_ids={str(getattr(m,n)) for m in cross_profile_matches for n in ("uaa_peak_id","uag_peak_id") if getattr(m,n,None)}
 by_peak=defaultdict(list);by_group=defaultdict(list);by_recur=defaultdict(list);series_index=defaultdict(list)
 for r in relations:by_peak[r.observed_peak_id].append(r);by_group[_fragment_key(r)].append(r);by_recur[(r.rna_identity_candidate,r.chemical_hypothesis_class,r.theoretical_fragment_id,r.reference_mass_definition)].append(r)
 for s in delta_result.state_series:
  for pid in s.member_observed_peak_ids:series_index[(pid,s.theoretical_fragment_identity,s.charge,s.ion_mode)].append(s)
 records=[]
 for r in sorted(relations,key=lambda x:x.t1_delta_relation_id):
  p=peak_map[r.observed_peak_id];group=by_group[_fragment_key(r)];peak_rows=by_peak[r.observed_peak_id];recur=by_recur[(r.rna_identity_candidate,r.chemical_hypothesis_class,r.theoretical_fragment_id,r.reference_mass_definition)]
  charge_to_peaks=defaultdict(set);charge_to_refs=defaultdict(set)
  for x in group:charge_to_peaks[x.charge].add(x.observed_peak_id);charge_to_refs[x.charge].add(x.reference_id)
  charges=tuple(sorted(charge_to_peaks));distinct_peaks=set().union(*charge_to_peaks.values()) if charge_to_peaks else set()
  if not charges:charge_class=ChargeSupportClass.NO_CHARGE_SUPPORT
  elif len(charges)==1:charge_class=ChargeSupportClass.SINGLE_CHARGE_ONLY
  elif len(distinct_peaks)>=len(charges) and len({ref for values in charge_to_refs.values() for ref in values})==1:charge_class=ChargeSupportClass.MULTI_CHARGE_CONCORDANT_SUPPORT
  else:charge_class=ChargeSupportClass.MULTI_CHARGE_CONFLICTING
  ac,cc,ar,cr,same,diff,concord=_concordance(r,refs,dp);shape=_shape(p);polcomp,polsupport=_polarity(polarity_status,r.ion_mode)
  raw=len(peak_rows);frags={(x.rna_identity_candidate,x.fragment_sequence,x.start_position,x.end_position) for x in peak_rows};rnas={x.rna_identity_candidate for x in peak_rows};allcharges={x.charge for x in peak_rows};cca={v for x in peak_rows for v in x.cca_state.removeprefix("COLLAPSED:").split("|")};categories={x.reference_category for x in peak_rows}
  if len(rnas)>1:ambiguity=FragmentAmbiguityClass.CROSS_RNA_AMBIGUITY
  elif len(frags)>1 or len(allcharges)>1 or len(categories)>1:ambiguity=FragmentAmbiguityClass.HIGH_FRAGMENT_AMBIGUITY
  elif len(cca)>1:ambiguity=FragmentAmbiguityClass.MODERATE_FRAGMENT_AMBIGUITY
  else:ambiguity=FragmentAmbiguityClass.LOW_FRAGMENT_AMBIGUITY
  cats,best,second,margin,resolution=_resolution(r,refs,dp,qp);disc,u,g,shared,same_mass,disc_support=_discrimination(r,discrimination)
  mass_status=MassDefinitionCompatibilityStatus.CONFIRMED_MISMATCH if r.delta_relation_class is DeltaRelationClass.MASS_DEFINITION_MISMATCH_DIAGNOSTIC else MassDefinitionCompatibilityStatus.CONFIRMED_COMPATIBLE if r.mass_definition_compatible else MassDefinitionCompatibilityStatus.UNKNOWN_COMPATIBILITY
  recurrent_peaks={x.observed_peak_id for x in recur};recurrent_charges={x.charge for x in recur};recurrent_sources={x.observed_source_id for x in all_rows if x.chemical_hypothesis_class is r.chemical_hypothesis_class and x.theoretical_fragment_id.split(":",1)[-1]==r.theoretical_fragment_id.split(":",1)[-1]}
  if len(recurrent_sources)>1:rec_class=RecurrentSupportClass.CROSS_SAMPLE_RECURRENT
  elif len(recurrent_peaks)>1 and len(recurrent_charges)>1:rec_class=RecurrentSupportClass.MULTI_PEAK_MULTI_CHARGE_RECURRENT
  elif len(recurrent_peaks)>1:rec_class=RecurrentSupportClass.MULTI_PEAK_SINGLE_CHARGE_RECURRENT
  elif len(recurrent_peaks)==1:rec_class=RecurrentSupportClass.SINGLE_PEAK_ONLY
  else:rec_class=RecurrentSupportClass.NO_RECURRENT_SUPPORT
  cross,inu,ing,crossamb=_sample_flags(all_rows,r);series_rows=series_index.get((r.observed_peak_id,r.theoretical_fragment_id,r.charge,r.ion_mode),[]);series=max(series_rows,key=lambda s:(s.member_count,s.strict_edge_count),default=None);series_support=_series_support(series)
  alternatives=set(cats)-{r.reference_category.value}
  if len(frags)>1:alternatives.add("DIFFERENT_FRAGMENT")
  if len(allcharges)>1:alternatives.add("CHARGE_AMBIGUITY")
  known=any(x.reference_category is ChemicalReferenceCategory.KNOWN_MODIFICATION_DELTA for x in peak_rows)
  alt_rel=[x for x in group if x.charge!=r.charge and x.reference_id==r.reference_id and x.observed_peak_id!=r.observed_peak_id]
  alt=min(alt_rel,key=lambda x:(abs(x.apex_delta_error_da or 0),x.observed_peak_id),default=None);expected=alt.theoretical_mz+r.reference_delta_da/alt.charge if alt else None;alt_error=alt.observed_apex_mz-expected if alt else None
  score=0
  if concord is ApexCentroidConcordanceClass.BOTH_STRICT_SAME_REFERENCE:score+=qp.both_strict_score
  if shape is PeakShapeSupportClass.STRONG_PEAK_SHAPE_SUPPORT:score+=qp.strong_shape_score
  if ambiguity is FragmentAmbiguityClass.LOW_FRAGMENT_AMBIGUITY:score+=qp.low_ambiguity_score
  if resolution in {ReferenceResolutionClass.WELL_RESOLVED_REFERENCE_CATEGORY,ReferenceResolutionClass.SINGLE_REFERENCE_ONLY}:score+=qp.resolved_reference_score
  if charge_class is ChargeSupportClass.MULTI_CHARGE_CONCORDANT_SUPPORT:score+=qp.multi_charge_score
  if rec_class in {RecurrentSupportClass.MULTI_PEAK_MULTI_CHARGE_RECURRENT,RecurrentSupportClass.MULTI_PEAK_SINGLE_CHARGE_RECURRENT,RecurrentSupportClass.CROSS_SAMPLE_RECURRENT}:score+=qp.recurrent_score
  if disc_support is DiscriminatorySupportClass.TARGET_SPECIFIC_FRAGMENT_SUPPORT:score+=qp.target_specific_score
  if series_support is StateSeriesSupportClass.STRONG_LINEAR_SERIES_SUPPORT:score+=qp.strong_series_score
  if shape is PeakShapeSupportClass.EXCLUDED_PEAK_SHAPE:score-=qp.artifact_penalty
  if ambiguity is FragmentAmbiguityClass.CROSS_RNA_AMBIGUITY or crossamb:score-=qp.cross_rna_penalty
  if resolution is ReferenceResolutionClass.UNRESOLVED_REFERENCE_CATEGORY:score-=qp.unresolved_reference_penalty
  if polsupport in {PolaritySupportClass.UNKNOWN_POLARITY,PolaritySupportClass.CONFLICTING_METADATA}:score-=qp.unknown_polarity_penalty
  if mass_status is MassDefinitionCompatibilityStatus.UNKNOWN_COMPATIBILITY:score-=qp.unknown_mass_definition_penalty
  if series_support is StateSeriesSupportClass.BRANCHED_AMBIGUOUS_SERIES:score-=qp.branched_series_penalty
  if charge_class is ChargeSupportClass.SINGLE_CHARGE_ONLY:score-=qp.single_charge_penalty
  score=max(0,min(100,score));warnings=[]
  for condition,name in [(polsupport is PolaritySupportClass.UNKNOWN_POLARITY,"UNKNOWN_POLARITY"),(mass_status is MassDefinitionCompatibilityStatus.UNKNOWN_COMPATIBILITY,"UNKNOWN_MASS_DEFINITION"),(charge_class is ChargeSupportClass.SINGLE_CHARGE_ONLY,"SINGLE_CHARGE_ONLY"),(ambiguity is FragmentAmbiguityClass.CROSS_RNA_AMBIGUITY,"CROSS_RNA_AMBIGUITY"),(ambiguity is FragmentAmbiguityClass.HIGH_FRAGMENT_AMBIGUITY,"MULTIPLE_FRAGMENT_AMBIGUITY"),(resolution is ReferenceResolutionClass.UNRESOLVED_REFERENCE_CATEGORY,"UNRESOLVED_REFERENCE_CATEGORY"),(concord not in _GOOD_CONCORD,"APEX_CENTROID_DISCORDANCE"),(series_support is StateSeriesSupportClass.BRANCHED_AMBIGUOUS_SERIES,"BRANCHED_SERIES")]:
   if condition:warnings.append(name)
  invalid=r.charge<=0 or r.reference_match_class not in _MATCHED
  priorities=[(invalid,EvidenceBlockReason.INVALID_RELATION),(polsupport is PolaritySupportClass.CONFIRMED_INCOMPATIBLE,EvidenceBlockReason.CONFIRMED_POLARITY_MISMATCH),(p.possible_isotope_component or p.possible_overlapping_envelope,EvidenceBlockReason.ISOTOPE_OR_ENVELOPE_COMPONENT),(p.possible_shoulder or p.possible_duplicate,EvidenceBlockReason.SHOULDER_OR_DUPLICATE),(mass_status is MassDefinitionCompatibilityStatus.CONFIRMED_MISMATCH,EvidenceBlockReason.MASS_DEFINITION_MISMATCH),(ambiguity is FragmentAmbiguityClass.CROSS_RNA_AMBIGUITY or crossamb,EvidenceBlockReason.CROSS_RNA_AMBIGUITY),(ambiguity is FragmentAmbiguityClass.HIGH_FRAGMENT_AMBIGUITY,EvidenceBlockReason.MULTIPLE_FRAGMENT_AMBIGUITY),(resolution is ReferenceResolutionClass.UNRESOLVED_REFERENCE_CATEGORY,EvidenceBlockReason.UNRESOLVED_REFERENCE_CATEGORY),(concord not in _GOOD_CONCORD,EvidenceBlockReason.APEX_CENTROID_DISCORDANCE),(charge_class is ChargeSupportClass.SINGLE_CHARGE_ONLY,EvidenceBlockReason.SINGLE_CHARGE_ONLY),(polsupport in {PolaritySupportClass.UNKNOWN_POLARITY,PolaritySupportClass.CONFLICTING_METADATA},EvidenceBlockReason.UNKNOWN_POLARITY),(mass_status is MassDefinitionCompatibilityStatus.UNKNOWN_COMPATIBILITY,EvidenceBlockReason.UNKNOWN_MASS_DEFINITION),(rec_class in {RecurrentSupportClass.SINGLE_PEAK_ONLY,RecurrentSupportClass.NO_RECURRENT_SUPPORT},EvidenceBlockReason.NO_RECURRENT_SUPPORT)]
  block=next((reason for condition,reason in priorities if condition),EvidenceBlockReason.NO_BLOCK)
  hard=block in {EvidenceBlockReason.INVALID_RELATION,EvidenceBlockReason.CONFIRMED_POLARITY_MISMATCH,EvidenceBlockReason.ISOTOPE_OR_ENVELOPE_COMPONENT,EvidenceBlockReason.SHOULDER_OR_DUPLICATE}
  diagnostic=mass_status is MassDefinitionCompatibilityStatus.CONFIRMED_MISMATCH or resolution is ReferenceResolutionClass.UNRESOLVED_REFERENCE_CATEGORY or series_support is StateSeriesSupportClass.BRANCHED_AMBIGUOUS_SERIES
  complete=(polsupport is PolaritySupportClass.CONFIRMED_COMPATIBLE and mass_status is MassDefinitionCompatibilityStatus.CONFIRMED_COMPATIBLE and concord in _GOOD_CONCORD and ambiguity is FragmentAmbiguityClass.LOW_FRAGMENT_AMBIGUITY and shape is PeakShapeSupportClass.STRONG_PEAK_SHAPE_SUPPORT and resolution is not ReferenceResolutionClass.UNRESOLVED_REFERENCE_CATEGORY and (charge_class is ChargeSupportClass.MULTI_CHARGE_CONCORDANT_SUPPORT or series_support is StateSeriesSupportClass.STRONG_LINEAR_SERIES_SUPPORT))
  unknowns=sum((polsupport in {PolaritySupportClass.UNKNOWN_POLARITY,PolaritySupportClass.CONFLICTING_METADATA},mass_status is MassDefinitionCompatibilityStatus.UNKNOWN_COMPATIBILITY))
  if hard:tier=EvidenceTier.TIER_E_BLOCKED
  elif diagnostic:tier=EvidenceTier.TIER_D_DIAGNOSTIC_ONLY
  elif complete:tier=EvidenceTier.TIER_A_HIGH_SUPPORT
  elif score>=qp.tier_b_minimum and unknowns<=1:tier=EvidenceTier.TIER_B_MODERATE_SUPPORT
  else:tier=EvidenceTier.TIER_C_WEAK_SUPPORT
  records.append(T1DeltaEvidenceQualityRecord(evidence_quality_record_id="T1QUALITY__"+sha256(r.t1_delta_relation_id.encode()).hexdigest()[:20].upper(),source_id=r.observed_source_id,measurement_id=measurement_id,rna_identity_candidate=r.rna_identity_candidate,t1_delta_relation_id=r.t1_delta_relation_id,observed_peak_id=r.observed_peak_id,theoretical_fragment_id=r.theoretical_fragment_id,fragment_sequence=r.fragment_sequence,start_position=r.start_position,end_position=r.end_position,cca_state=r.cca_state,ion_mode=r.ion_mode,charge=r.charge,chemical_hypothesis_class=r.chemical_hypothesis_class,reference_id=r.reference_id,reference_name=r.reference_name,reference_mass_definition=r.reference_mass_definition,observed_apex_mz=r.observed_apex_mz,observed_centroid_mz=r.observed_centroid_mz,theoretical_mz=r.theoretical_mz,apex_neutral_delta=r.apex_neutral_delta,centroid_neutral_delta=r.centroid_neutral_delta,reference_delta_da=r.reference_delta_da,apex_error_da=r.apex_delta_error_da,centroid_error_da=r.centroid_delta_error_da,observed_polarity_status=polarity_status,ion_mode_polarity_compatible=polcomp,polarity_support_class=polsupport,charge_support_count=len(charges),supported_charge_states=charges,distinct_observed_peak_count=len(distinct_peaks),charge_support_class=charge_class,apex_reference_match_class=ac,centroid_reference_match_class=cc,apex_reference_category=ar.reference_category if ar else None,centroid_reference_category=cr.reference_category if cr else None,apex_centroid_same_reference=same,apex_centroid_error_difference=diff,apex_centroid_concordance_class=concord,peak_quality_class=p.peak_quality_class,fwhm=p.fwhm_mz,prominence=p.prominence,sharpness=p.sharpness_score,possible_overlapping_envelope=p.possible_overlapping_envelope,peak_shape_support_class=shape,distinct_fragment_count=len(frags),distinct_charge_count=len(allcharges),distinct_cca_state_count=len(cca),distinct_rna_identity_count=len(rnas),candidate_reference_count=r.candidate_reference_count,distinct_reference_category_count=r.distinct_reference_category_count,raw_candidate_count=raw,cca_collapsed_candidate_count=len(frags),fragment_identity_collapsed_count=len(frags),fragment_ambiguity_class=ambiguity,fragment_discrimination_class=disc,uaa_specific_theoretical=u,uag_specific_theoretical=g,shared_theoretical=shared,same_mass_different_sequence=same_mass,cross_rna_identity_ambiguous=crossamb,discriminatory_support_class=DiscriminatorySupportClass.CROSS_RNA_AMBIGUITY if crossamb else disc_support,reference_category_count=len(cats),reference_category_list=cats,best_category_error=best,second_best_category_error=second,reference_resolution_margin_da=margin,reference_resolution_class=resolution,observed_mass_definition_status="UNKNOWN",mass_definition_compatibility_status=mass_status,recurrent_support_group_id="T1RECUR__"+sha256(repr((r.rna_identity_candidate,r.chemical_hypothesis_class.value,r.theoretical_fragment_id,r.reference_mass_definition.value)).encode()).hexdigest()[:16].upper(),recurrent_observed_peak_count=len(recurrent_peaks),recurrent_fragment_identity_count=len({x.theoretical_fragment_id for x in recur}),recurrent_charge_count=len(recurrent_charges),recurrent_source_count=len(recurrent_sources),recurrent_support_class=rec_class,cross_sample_support_status=cross,supported_in_uaa=inu,supported_in_uag=ing,cross_sample_fragment_identity_match=inu and ing,cross_sample_chemical_category_match=inu and ing,common_component_or_shared_fragment_possible=(inu and ing) or r.observed_peak_id in cross_profile_peak_ids,member_of_state_series=series is not None,state_series_id=series.t1_state_series_id if series else None,state_series_member_count=series.member_count if series else 0,state_series_strict_edge_count=series.strict_edge_count if series else 0,state_series_pattern=series.series_pattern if series else None,state_series_support_class=series_support,alternative_explanation_count=len(alternatives),alternative_explanation_categories=tuple(sorted(alternatives)),adduct_alternative_possible=True,water_loss_alternative_possible=r.chemical_hypothesis_class in {ChemicalHypothesisClass.H2O_EQUIVALENT,ChemicalHypothesisClass.H2O_LOSS_EQUIVALENT,ChemicalHypothesisClass.O_TO_S_EQUIVALENT,ChemicalHypothesisClass.S_TO_O_EQUIVALENT},oxidation_alternative_possible=r.chemical_hypothesis_class in {ChemicalHypothesisClass.O_EQUIVALENT,ChemicalHypothesisClass.O_TO_S_EQUIVALENT,ChemicalHypothesisClass.S_TO_O_EQUIVALENT},known_modification_alternative_possible=known,different_fragment_alternative_possible=len(frags)>1,expected_mz_for_alternative_charge=expected,observed_supporting_peak_id=alt.observed_peak_id if alt else None,alternative_charge_error_mz=alt_error,alternative_charge_support_class="CONCORDANT_OBSERVED" if alt else "NOT_OBSERVED",evidence_tier=tier,evidence_support_score=score,evidence_block_reason=block,evidence_warning_flags=tuple(warnings),eligible_for_high_quality_shadow_evidence=tier in {EvidenceTier.TIER_A_HIGH_SUPPORT,EvidenceTier.TIER_B_MODERATE_SUPPORT}))
 records=tuple(sorted(records,key=lambda x:x.evidence_quality_record_id));priority=[]
 for q in records:
  if q.chemical_hypothesis_class not in {ChemicalHypothesisClass.O_TO_S_EQUIVALENT,ChemicalHypothesisClass.S_TO_O_EQUIVALENT}:continue
  if q.evidence_tier is EvidenceTier.TIER_E_BLOCKED:pc=PriorityOSClass.BLOCKED_O_S_CANDIDATE
  elif q.evidence_tier is EvidenceTier.TIER_D_DIAGNOSTIC_ONLY:pc=PriorityOSClass.O_TO_S_DIAGNOSTIC_ONLY if q.chemical_hypothesis_class is ChemicalHypothesisClass.O_TO_S_EQUIVALENT else PriorityOSClass.S_TO_O_DIAGNOSTIC_ONLY
  else:pc=PriorityOSClass.O_TO_S_PRIORITY_CANDIDATE if q.chemical_hypothesis_class is ChemicalHypothesisClass.O_TO_S_EQUIVALENT else PriorityOSClass.S_TO_O_PRIORITY_CANDIDATE
  priority.append(OSPriorityQualityRecord(quality_record=q,priority_class=pc))
 flags=glu_summary_flags or {}
 return T1DeltaEvidenceQualityAuditResult(source_id=delta_result.source_id,records=records,os_priority_records=tuple(priority),input_delta_relation_count=len(delta_result.relations),quality_record_count=len(records),oxygen_equivalent_pattern_also_observed_in_glu_intact=bool(flags.get("oxygen")),water_equivalent_pattern_also_observed_in_glu_intact=bool(flags.get("water")))
