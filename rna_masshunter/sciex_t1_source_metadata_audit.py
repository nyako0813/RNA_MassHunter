"""Read-only source-level metadata shadow audit for SCIEX T1 profiles."""
from __future__ import annotations
from bisect import bisect_left,bisect_right
from collections import Counter,defaultdict
from dataclasses import dataclass
from enum import Enum
from hashlib import sha256
from pathlib import Path
from statistics import mean,median,pstdev
from typing import Iterable,Mapping,Sequence
import numpy as np
from rna_masshunter.sciex_t1_profile_peak_audit import ISOTOPE_MASS_DIFFERENCE_DA,T1ProfilePeak,T1Safeguards
from rna_masshunter.sciex_t1_fragment_shadow_match import T1FragmentMatch,T1FragmentMatchClass,T1IonCandidate,T1IonMode
from rna_masshunter.sciex_t1_fragment_delta_audit import ChemicalHypothesisClass
from rna_masshunter.sciex_t1_delta_evidence_quality_audit import (
 EvidenceTier,MassDefinitionCompatibilityStatus,ObservedPolarityStatus,
 RecurrentSupportClass,T1DeltaEvidenceQualityRecord,
)

class MetadataSourceType(str,Enum):
 TEXT_FILE_HEADER="TEXT_FILE_HEADER";TEXT_FILE_PREAMBLE="TEXT_FILE_PREAMBLE";TEXT_FILE_COMMENTS="TEXT_FILE_COMMENTS";TEXT_FILE_COLUMNS="TEXT_FILE_COLUMNS";FILENAME="FILENAME";PARENT_DIRECTORY_NAME="PARENT_DIRECTORY_NAME";MANIFEST="MANIFEST";SOURCE_REGISTRY="SOURCE_REGISTRY";PROJECT_CONFIG="PROJECT_CONFIG";ADJACENT_METADATA_FILE="ADJACENT_METADATA_FILE";MZML_OR_RAW_SOURCE_REFERENCE="MZML_OR_RAW_SOURCE_REFERENCE";INSTRUMENT_EXPORT_CONVENTION="INSTRUMENT_EXPORT_CONVENTION";EMPIRICAL_PROFILE_STRUCTURE="EMPIRICAL_PROFILE_STRUCTURE";USER_CONFIRMED_METADATA="USER_CONFIRMED_METADATA"
class MetadataTrustLevel(str,Enum):
 AUTHORITATIVE_SOURCE_SPECIFIC="AUTHORITATIVE_SOURCE_SPECIFIC";AUTHORITATIVE_EXPERIMENT_LEVEL="AUTHORITATIVE_EXPERIMENT_LEVEL";PROJECT_DECLARED_SOURCE_SPECIFIC="PROJECT_DECLARED_SOURCE_SPECIFIC";PROJECT_DECLARED_GENERAL="PROJECT_DECLARED_GENERAL";INFERRED_FROM_DATA="INFERRED_FROM_DATA";FILENAME_ONLY="FILENAME_ONLY";UNVERIFIED="UNVERIFIED";NOT_AVAILABLE="NOT_AVAILABLE"
class PolarityStatus(str,Enum):NEGATIVE="NEGATIVE";POSITIVE="POSITIVE";MIXED="MIXED";UNKNOWN="UNKNOWN";CONFLICT="CONFLICT";NOT_APPLICABLE="NOT_APPLICABLE"
class DataRepresentationStatus(str,Enum):
 PROFILE_CONTINUOUS="PROFILE_CONTINUOUS";CENTROID_PEAK_LIST="CENTROID_PEAK_LIST";RESAMPLED_PROFILE="RESAMPLED_PROFILE";DECONVOLUTED_PROFILE="DECONVOLUTED_PROFILE";RECONSTRUCTED_PROFILE="RECONSTRUCTED_PROFILE";UNKNOWN_REPRESENTATION="UNKNOWN_REPRESENTATION";CONFLICTING_REPRESENTATION="CONFLICTING_REPRESENTATION"
class ObservedMassDefinitionStatus(str,Enum):
 MONOISOTOPIC_MZ="MONOISOTOPIC_MZ";AVERAGE_MZ="AVERAGE_MZ";CENTROID_OF_ISOTOPE_ENVELOPE="CENTROID_OF_ISOTOPE_ENVELOPE";APEX_OF_PROFILE_SIGNAL="APEX_OF_PROFILE_SIGNAL";DECONVOLUTED_NEUTRAL_TO_MZ_EXPORT="DECONVOLUTED_NEUTRAL_TO_MZ_EXPORT";UNKNOWN_MASS_DEFINITION="UNKNOWN_MASS_DEFINITION";CONFLICTING_MASS_DEFINITION="CONFLICTING_MASS_DEFINITION"
class ExportTransformationStatus(str,Enum):
 DIRECT_PROFILE_EXPORT="DIRECT_PROFILE_EXPORT";CENTROID_EXPORT="CENTROID_EXPORT";DECONVOLUTED_EXPORT="DECONVOLUTED_EXPORT";RECONSTRUCTED_EXPORT="RECONSTRUCTED_EXPORT";RESAMPLED_EXPORT="RESAMPLED_EXPORT";UNKNOWN_EXPORT_TRANSFORMATION="UNKNOWN_EXPORT_TRANSFORMATION"
class ConflictResolutionStatus(str,Enum):RESOLVED_BY_HIGHER_TRUST_SOURCE="RESOLVED_BY_HIGHER_TRUST_SOURCE";UNRESOLVED_CONFLICT="UNRESOLVED_CONFLICT";LOWER_TRUST_INFERENCE_IGNORED="LOWER_TRUST_INFERENCE_IGNORED";NO_CONFLICT="NO_CONFLICT"
class MetadataCompletenessClass(str,Enum):COMPLETE_FOR_CHEMICAL_ASSIGNMENT="COMPLETE_FOR_CHEMICAL_ASSIGNMENT";SUFFICIENT_FOR_TIER_B_SHADOW_EVIDENCE="SUFFICIENT_FOR_TIER_B_SHADOW_EVIDENCE";SUFFICIENT_FOR_TIER_C_DIAGNOSTIC="SUFFICIENT_FOR_TIER_C_DIAGNOSTIC";INSUFFICIENT_FOR_CHEMICAL_ASSIGNMENT="INSUFFICIENT_FOR_CHEMICAL_ASSIGNMENT";CONFLICTING_METADATA="CONFLICTING_METADATA"
class EnvelopeSupportClass(str,Enum):STRONG_MULTI_MEMBER_ENVELOPE_SUPPORT="STRONG_MULTI_MEMBER_ENVELOPE_SUPPORT";MODERATE_TWO_MEMBER_SUPPORT="MODERATE_TWO_MEMBER_SUPPORT";WEAK_SPACING_ONLY="WEAK_SPACING_ONLY";BRANCHED_OR_OVERLAPPING="BRANCHED_OR_OVERLAPPING"
class ChargePlausibilityClass(str,Enum):SUPPORTED_CHARGE="SUPPORTED_CHARGE";WEAKLY_SUPPORTED_CHARGE="WEAKLY_SUPPORTED_CHARGE";INSUFFICIENT_SUPPORT="INSUFFICIENT_SUPPORT";CONFLICTING_SUPPORT="CONFLICTING_SUPPORT"
class SimulationScenarioID(str,Enum):CURRENT_METADATA="CURRENT_METADATA";NEGATIVE_POLARITY_CONFIRMED_ONLY="NEGATIVE_POLARITY_CONFIRMED_ONLY";POSITIVE_POLARITY_CONFIRMED_ONLY="POSITIVE_POLARITY_CONFIRMED_ONLY";MONOISOTOPIC_MZ_CONFIRMED_ONLY="MONOISOTOPIC_MZ_CONFIRMED_ONLY";NEGATIVE_PLUS_MONO_CONFIRMED="NEGATIVE_PLUS_MONO_CONFIRMED";POSITIVE_PLUS_MONO_CONFIRMED="POSITIVE_PLUS_MONO_CONFIRMED"

_TRUST_RANK={MetadataTrustLevel.AUTHORITATIVE_SOURCE_SPECIFIC:7,MetadataTrustLevel.AUTHORITATIVE_EXPERIMENT_LEVEL:6,MetadataTrustLevel.PROJECT_DECLARED_SOURCE_SPECIFIC:5,MetadataTrustLevel.PROJECT_DECLARED_GENERAL:4,MetadataTrustLevel.INFERRED_FROM_DATA:3,MetadataTrustLevel.FILENAME_ONLY:2,MetadataTrustLevel.UNVERIFIED:1,MetadataTrustLevel.NOT_AVAILABLE:0}
@dataclass(frozen=True)
class MetadataSourceInventoryRecord:
 metadata_source_id:str;metadata_source_type:MetadataSourceType;source_path_or_reference:str;available:bool;read_performed:bool;content_summary:str;trust_level:MetadataTrustLevel;applicable_to_runtime_source:bool
@dataclass(frozen=True)
class MetadataFieldEvidence:
 field_name:str;value:str;source_id:str;trust_level:MetadataTrustLevel;source_specific:bool;confirmed:bool
@dataclass(frozen=True)
class MetadataConflict:
 metadata_conflict_id:str;metadata_field:str;evidence_source_a:str;value_a:str;trust_a:MetadataTrustLevel;evidence_source_b:str;value_b:str;trust_b:MetadataTrustLevel;conflict_resolution_status:ConflictResolutionStatus
@dataclass(frozen=True)
class SpacingDiagnostics:
 row_count:int;mz_minimum:float|None;mz_maximum:float|None;median_spacing:float|None;mean_spacing:float|None;minimum_positive_spacing:float|None;maximum_spacing:float|None;spacing_standard_deviation:float|None;spacing_coefficient_of_variation:float|None;unique_spacing_count:int;most_common_spacing:float|None;zero_or_negative_spacing_count:int;duplicate_coordinate_count:int
@dataclass(frozen=True)
class IntensityDiagnostics:
 minimum_intensity:float|None;maximum_intensity:float|None;median_intensity:float|None;mean_intensity:float|None;zero_intensity_count:int;negative_intensity_count:int;positive_intensity_count:int;non_finite_intensity_count:int;dynamic_range:float|None
@dataclass(frozen=True,kw_only=True)
class MetadataSafeguards(T1Safeguards):
 metadata_audit_only:bool=True;scenario_assumption_only:bool=False;applied_to_runtime_result:bool=False;representation_confirmed:bool=False;export_transformation_confirmed:bool=False;mass_definition_confirmed:bool=False
@dataclass(frozen=True,kw_only=True)
class IsotopeFamilyCandidate(MetadataSafeguards):
 isotope_family_id:str;member_peak_ids:tuple[str,...];member_apex_mzs:tuple[float,...];member_count:int;possible_charge:int;mean_spacing:float;maximum_spacing_error:float;intensity_pattern:str;envelope_support_class:EnvelopeSupportClass;isotope_envelope_confirmed:bool=False
@dataclass(frozen=True,kw_only=True)
class ChargeSupportSummary(MetadataSafeguards):
 charge:int;theoretical_ion_candidate_count:int;direct_fragment_match_count:int;chemical_delta_relation_count:int;multi_charge_concordant_count:int;isotope_family_support_count:int;distinct_observed_peak_count:int;charge_support_score:int;charge_plausibility_class:ChargePlausibilityClass
@dataclass(frozen=True,kw_only=True)
class SourceMetadataAuditRecord(MetadataSafeguards):
 t1_source_metadata_audit_id:str;source_id:str;measurement_id:str;rna_identity:str;runtime_path:str;sha256:str;coordinate_type:str;observed_mass_scale:str;observed_output_species:str;acquisition_polarity:PolarityStatus;exported_ion_polarity:PolarityStatus;polarity_confirmed:bool;polarity_evidence_source:tuple[str,...];polarity_evidence_trust_level:MetadataTrustLevel;polarity_evidence_count:int;polarity_conflict_count:int;polarity_final_status:PolarityStatus;data_representation_status:DataRepresentationStatus;representation_evidence_sources:tuple[str,...];representation_trust_level:MetadataTrustLevel;observed_mass_definition_status:ObservedMassDefinitionStatus;mass_definition_evidence_source:tuple[str,...];mass_definition_trust_level:MetadataTrustLevel;export_transformation_status:ExportTransformationStatus;configured_charge_range:tuple[int,...];empirically_supported_charges:tuple[int,...];weakly_supported_charges:tuple[int,...];unsupported_charges:tuple[int,...];recommended_shadow_charge_range:str;metadata_completeness_class:MetadataCompletenessClass;metadata_conflict_status:ConflictResolutionStatus;source_linked_raw_metadata_available:bool
@dataclass(frozen=True,kw_only=True)
class TierSimulationResult(MetadataSafeguards):
 scenario_id:SimulationScenarioID;assumed_polarity:PolarityStatus;assumed_mass_definition:ObservedMassDefinitionStatus;relation_count:int;tier_a_count:int;tier_b_count:int;tier_c_count:int;tier_d_count:int;tier_e_count:int;surviving_chemical_relation_count:int;surviving_distinct_peak_count:int;category_survival:tuple[tuple[str,int],...];multi_charge_concordant_count:int;average_reference_surviving_count:int;monoisotopic_reference_surviving_count:int;known_modification_diagnostic_count:int;reference_category_conflict_count:int;apex_centroid_concordant_count:int;scenario_assumption_only:bool=True
@dataclass(frozen=True,kw_only=True)
class PriorityCandidateScenarioStatus(MetadataSafeguards):
 observed_peak_id:str;t1_delta_relation_id:str;scenario_id:SimulationScenarioID;current_tier:EvidenceTier;simulated_tier:EvidenceTier;remains_candidate:bool;blocked_by_polarity:bool;blocked_by_mass_definition:bool;score_change:int;scenario_assumption_only:bool=True

@dataclass(frozen=True,kw_only=True)
class T1SourceMetadataAuditResult(MetadataSafeguards):
 source_record:SourceMetadataAuditRecord;inventory:tuple[MetadataSourceInventoryRecord,...];conflicts:tuple[MetadataConflict,...];spacing_diagnostics:SpacingDiagnostics;intensity_diagnostics:IntensityDiagnostics;isotope_families:tuple[IsotopeFamilyCandidate,...];charge_summaries:tuple[ChargeSupportSummary,...];simulations:tuple[TierSimulationResult,...]

@dataclass(frozen=True)
class SourceMetadataAuditParameters:
 isotope_spacing_tolerance_mz:float=.01;strong_isotope_spacing_error_mz:float=.003;supported_charge_score:int=60;weak_charge_score:int=20;configured_charges:tuple[int,...]=(1,2,3,4,5)
 def validate(self):
  if self.isotope_spacing_tolerance_mz<=0 or not 0<self.strong_isotope_spacing_error_mz<=self.isotope_spacing_tolerance_mz:raise ValueError("invalid isotope tolerances")
  if not 0<=self.weak_charge_score<=self.supported_charge_score<=100:raise ValueError("invalid charge thresholds")

def resolve_metadata_field(evidence:Sequence[MetadataFieldEvidence]):
 available=[e for e in evidence if e.value and e.value not in {"UNKNOWN","NOT_AVAILABLE"}]
 if not available:return "UNKNOWN",False,(),ConflictResolutionStatus.NO_CONFLICT
 highest=max((_TRUST_RANK[e.trust_level] for e in available),default=0);top=[e for e in available if _TRUST_RANK[e.trust_level]==highest];values={e.value for e in top};conflicts=[]
 for i,a in enumerate(available):
  for b in available[i+1:]:
   if a.value==b.value:continue
   if _TRUST_RANK[a.trust_level]==_TRUST_RANK[b.trust_level]:status=ConflictResolutionStatus.UNRESOLVED_CONFLICT
   elif max(_TRUST_RANK[a.trust_level],_TRUST_RANK[b.trust_level])==highest:status=ConflictResolutionStatus.RESOLVED_BY_HIGHER_TRUST_SOURCE
   else:status=ConflictResolutionStatus.LOWER_TRUST_INFERENCE_IGNORED
   conflicts.append(MetadataConflict("METACONFLICT__"+sha256(f"{a.field_name}|{a.source_id}|{b.source_id}".encode()).hexdigest()[:16].upper(),a.field_name,a.source_id,a.value,a.trust_level,b.source_id,b.value,b.trust_level,status))
 if len(values)>1:return "CONFLICT",False,tuple(conflicts),ConflictResolutionStatus.UNRESOLVED_CONFLICT
 winner=top[0];confirmed=winner.confirmed and winner.source_specific and winner.trust_level in {MetadataTrustLevel.AUTHORITATIVE_SOURCE_SPECIFIC,MetadataTrustLevel.AUTHORITATIVE_EXPERIMENT_LEVEL,MetadataTrustLevel.PROJECT_DECLARED_SOURCE_SPECIFIC}
 status=ConflictResolutionStatus.RESOLVED_BY_HIGHER_TRUST_SOURCE if conflicts else ConflictResolutionStatus.NO_CONFLICT
 return winner.value,confirmed,tuple(conflicts),status

def calculate_spacing_diagnostics(coordinates:Sequence[float]):
 x=np.asarray(tuple(coordinates),dtype=float);n=len(x)
 if not n:return SpacingDiagnostics(0,None,None,None,None,None,None,None,None,0,None,0,0)
 dif=np.diff(x);positive=dif[dif>0];rounded=np.round(positive,9);counts=Counter(float(v) for v in rounded);common=min(counts,key=lambda v:(-counts[v],v)) if counts else None;avg=float(np.mean(positive)) if len(positive) else None;sd=float(np.std(positive)) if len(positive) else None
 return SpacingDiagnostics(n,float(np.min(x)),float(np.max(x)),float(np.median(positive)) if len(positive) else None,avg,float(np.min(positive)) if len(positive) else None,float(np.max(positive)) if len(positive) else None,sd,sd/avg if avg else None,len(counts),common,int(np.sum(dif<=0)),int(np.sum(dif==0)))
def calculate_intensity_diagnostics(intensities:Sequence[float]):
 y=np.asarray(tuple(intensities),dtype=float);finite=y[np.isfinite(y)]
 if not len(finite):return IntensityDiagnostics(None,None,None,None,0,0,0,len(y),None)
 positive=finite[finite>0];minimum=float(np.min(finite));maximum=float(np.max(finite));small=float(np.min(positive)) if len(positive) else None
 return IntensityDiagnostics(minimum,maximum,float(np.median(finite)),float(np.mean(finite)),int(np.sum(finite==0)),int(np.sum(finite<0)),int(np.sum(finite>0)),int(np.sum(~np.isfinite(y))),maximum/small if small and small>0 else None)

def build_metadata_source_inventory(runtime_path:Path,loaded,source,measurement,*,project_config_path:Path|None=None,user_metadata:Mapping[str,str]|None=None):
 path=Path(runtime_path);items=[]
 def add(kind,ref,available,read,summary,trust,applicable=True):items.append(MetadataSourceInventoryRecord("METASOURCE__"+sha256(f"{source.profile_source_id}|{kind.value}|{ref}".encode()).hexdigest()[:16].upper(),kind,str(ref),available,read,summary,trust,applicable))
 add(MetadataSourceType.TEXT_FILE_HEADER,path,True,True,";".join(loaded.header),MetadataTrustLevel.AUTHORITATIVE_SOURCE_SPECIFIC)
 add(MetadataSourceType.TEXT_FILE_PREAMBLE,path,False,True,"HEADER_IMMEDIATELY_FOLLOWED_BY_NUMERIC_ROWS",MetadataTrustLevel.NOT_AVAILABLE)
 add(MetadataSourceType.TEXT_FILE_COMMENTS,path,False,True,"NO_COMMENT_LINES_OBSERVED_AT_FILE_START",MetadataTrustLevel.NOT_AVAILABLE)
 add(MetadataSourceType.TEXT_FILE_COLUMNS,path,True,True,",".join(loaded.header),MetadataTrustLevel.AUTHORITATIVE_SOURCE_SPECIFIC)
 add(MetadataSourceType.FILENAME,path.name,True,True,"FILENAME_HAS_SAMPLE_LABEL_ONLY",MetadataTrustLevel.FILENAME_ONLY)
 add(MetadataSourceType.PARENT_DIRECTORY_NAME,path.parent.name,True,True,"RUNTIME_CACHE_DIRECTORY",MetadataTrustLevel.FILENAME_ONLY,False)
 add(MetadataSourceType.MANIFEST,"data/sciex_sample_manifest.yaml",True,True,f"measurement={measurement.measurement_id};source_reference={measurement.source_file_name}",MetadataTrustLevel.AUTHORITATIVE_EXPERIMENT_LEVEL)
 add(MetadataSourceType.SOURCE_REGISTRY,"sciex_reconstructed_profile_registry",True,True,f"type={source.profile_type.value};coordinate={source.mass_column}",MetadataTrustLevel.PROJECT_DECLARED_SOURCE_SPECIFIC)
 config_available=bool(project_config_path and project_config_path.is_file());config_text=project_config_path.read_text(encoding="utf-8",errors="replace") if config_available else "";config_summary="GENERAL_POLARITY_DECLARATION_PRESENT_NOT_SOURCE_SPECIFIC" if "polarity:" in config_text else "GENERAL_PROJECT_CONFIGURATION_NOT_SOURCE_SPECIFIC";add(MetadataSourceType.PROJECT_CONFIG,project_config_path or "NOT_PROVIDED",config_available,config_available,config_summary,MetadataTrustLevel.PROJECT_DECLARED_GENERAL,False)
 patterns={".txt",".csv",".tsv",".xml",".mzml",".wiff",".scan",".json",".yaml",".yml",".log",".method",".report"};adj=[p for p in path.parent.iterdir() if p.is_file() and p!=path and (p.suffix.lower() in patterns or p.name.lower().endswith(".wiff.scan"))];add(MetadataSourceType.ADJACENT_METADATA_FILE,path.parent,bool(adj),True,";".join(f"{p.name}:{p.stat().st_size}" for p in sorted(adj)) or "NO_METADATA_CANDIDATE",MetadataTrustLevel.UNVERIFIED,False)
 raw_name=str(getattr(measurement,"source_file_name","") or "");raw_hint=str(getattr(measurement,"source_file_path_hint","") or "");raw_available=bool(raw_hint and Path(raw_hint).is_file());add(MetadataSourceType.MZML_OR_RAW_SOURCE_REFERENCE,raw_hint or raw_name or "NOT_AVAILABLE",raw_available,False,"REFERENCE_PRESENT_BUT_RUNTIME_PATH_UNAVAILABLE" if raw_name and not raw_available else "AVAILABLE" if raw_available else "NOT_AVAILABLE",MetadataTrustLevel.AUTHORITATIVE_EXPERIMENT_LEVEL)
 add(MetadataSourceType.INSTRUMENT_EXPORT_CONVENTION,"REPOSITORY_DOCUMENTATION",False,True,"NO_SOURCE_SPECIFIC_EXPORT_CONVENTION",MetadataTrustLevel.NOT_AVAILABLE)
 add(MetadataSourceType.EMPIRICAL_PROFILE_STRUCTURE,path,True,True,"NUMERIC_MZ_INTENSITY_STRUCTURE_DIAGNOSTIC_ONLY",MetadataTrustLevel.INFERRED_FROM_DATA)
 user=user_metadata or {};add(MetadataSourceType.USER_CONFIRMED_METADATA,"RUNTIME_ARGUMENT",bool(user),bool(user),";".join(f"{k}={v}" for k,v in sorted(user.items())) or "NOT_PROVIDED",MetadataTrustLevel.AUTHORITATIVE_SOURCE_SPECIFIC)
 return tuple(items)

def build_isotope_families(peaks:Sequence[T1ProfilePeak],parameters:SourceMetadataAuditParameters|None=None):
 p=parameters or SourceMetadataAuditParameters();p.validate();ordered=tuple(sorted(peaks,key=lambda x:(x.apex_mz,x.t1_peak_id)));masses=[x.apex_mz for x in ordered];by_charge=defaultdict(list)
 for child in ordered:
  z=child.possible_isotope_charge
  if not child.possible_isotope_component or not z:continue
  expected=ISOTOPE_MASS_DIFFERENCE_DA/z;target=child.apex_mz-expected;lo=bisect_left(masses,target-p.isotope_spacing_tolerance_mz);hi=bisect_right(masses,target+p.isotope_spacing_tolerance_mz);parents=[x for x in ordered[lo:hi] if x.t1_peak_id!=child.t1_peak_id]
  if parents:
   parent=min(parents,key=lambda x:(abs((child.apex_mz-x.apex_mz)-expected),x.t1_peak_id));by_charge[z].append((parent.t1_peak_id,child.t1_peak_id,abs((child.apex_mz-parent.apex_mz)-expected)))
 peakmap={x.t1_peak_id:x for x in ordered};out=[]
 for z,edges in sorted(by_charge.items()):
  adj=defaultdict(set);errors={}
  for a,b,e in edges:adj[a].add(b);adj[b].add(a);errors[tuple(sorted((a,b)))]=e
  unseen=set(adj)
  while unseen:
   start=min(unseen,key=lambda i:(peakmap[i].apex_mz,i));stack=[start];members=set()
   while stack:
    cur=stack.pop()
    if cur in members:continue
    members.add(cur);unseen.discard(cur);stack.extend(adj[cur]-members)
   ids=tuple(sorted(members,key=lambda i:(peakmap[i].apex_mz,i)));component_errors=[e for pair,e in errors.items() if pair[0] in members and pair[1] in members];degrees=[len(adj[i]&members) for i in members];spacings=[peakmap[ids[i]].apex_mz-peakmap[ids[i-1]].apex_mz for i in range(1,len(ids))];maxerr=max(component_errors,default=0.)
   if max(degrees,default=0)>2:klass=EnvelopeSupportClass.BRANCHED_OR_OVERLAPPING
   elif len(ids)>=3 and maxerr<=p.strong_isotope_spacing_error_mz:klass=EnvelopeSupportClass.STRONG_MULTI_MEMBER_ENVELOPE_SUPPORT
   elif len(ids)==2 and maxerr<=p.strong_isotope_spacing_error_mz:klass=EnvelopeSupportClass.MODERATE_TWO_MEMBER_SUPPORT
   else:klass=EnvelopeSupportClass.WEAK_SPACING_ONLY
   ints=[peakmap[i].apex_intensity for i in ids];pattern="MONOTONIC_DECREASING" if all(a>=b for a,b in zip(ints,ints[1:])) else "NON_MONOTONIC"
   out.append(IsotopeFamilyCandidate(isotope_family_id="T1ISOFAMILY__"+sha256(f"{z}|{'|'.join(ids)}".encode()).hexdigest()[:20].upper(),member_peak_ids=ids,member_apex_mzs=tuple(peakmap[i].apex_mz for i in ids),member_count=len(ids),possible_charge=z,mean_spacing=mean(spacings) if spacings else 0.,maximum_spacing_error=maxerr,intensity_pattern=pattern,envelope_support_class=klass))
 return tuple(sorted(out,key=lambda x:(x.possible_charge,x.member_apex_mzs,x.isotope_family_id)))

def build_charge_support_summaries(ions:Sequence[T1IonCandidate],fragment_matches:Sequence[T1FragmentMatch],quality_records:Sequence[T1DeltaEvidenceQualityRecord],families:Sequence[IsotopeFamilyCandidate],*,polarity:PolarityStatus=PolarityStatus.UNKNOWN,parameters:SourceMetadataAuditParameters|None=None):
 p=parameters or SourceMetadataAuditParameters();p.validate();out=[]
 for z in p.configured_charges:
  im=[x for x in ions if x.charge==z];fm=[x for x in fragment_matches if x.charge==z and x.match_class is not T1FragmentMatchClass.NO_MATCH];qr=[x for x in quality_records if x.charge==z];ff=[x for x in families if x.possible_charge==z];peaks={x.observed_peak_id for x in qr}|{x.observed_peak_id for x in fm};multi=sum(x.charge_support_class.value=="MULTI_CHARGE_CONCORDANT_SUPPORT" for x in qr);recur=sum(x.recurrent_support_class in {RecurrentSupportClass.MULTI_PEAK_MULTI_CHARGE_RECURRENT,RecurrentSupportClass.MULTI_PEAK_SINGLE_CHARGE_RECURRENT,RecurrentSupportClass.CROSS_SAMPLE_RECURRENT} for x in qr);score=0
  if any(x.envelope_support_class is EnvelopeSupportClass.STRONG_MULTI_MEMBER_ENVELOPE_SUPPORT for x in ff):score+=30
  if len(fm)>=2:score+=20
  if multi:score+=20
  if recur:score+=15
  if len(peaks)>=5:score+=10
  if len(peaks)<=1 and (fm or qr) and not multi:score-=20
  if qr and all(x.cross_rna_identity_ambiguous for x in qr):score-=20
  compatible=[x for x in qr if (polarity is PolarityStatus.NEGATIVE and x.ion_mode is T1IonMode.NEGATIVE_DEPROTONATED) or (polarity is PolarityStatus.POSITIVE and x.ion_mode is T1IonMode.POSITIVE_PROTONATED)]
  conflict=polarity in {PolarityStatus.NEGATIVE,PolarityStatus.POSITIVE} and qr and not compatible
  if conflict:score-=20
  score=max(0,min(100,score));klass=ChargePlausibilityClass.CONFLICTING_SUPPORT if conflict else ChargePlausibilityClass.SUPPORTED_CHARGE if score>=p.supported_charge_score else ChargePlausibilityClass.WEAKLY_SUPPORTED_CHARGE if score>=p.weak_charge_score else ChargePlausibilityClass.INSUFFICIENT_SUPPORT
  out.append(ChargeSupportSummary(charge=z,theoretical_ion_candidate_count=len(im),direct_fragment_match_count=len(fm),chemical_delta_relation_count=len(qr),multi_charge_concordant_count=multi,isotope_family_support_count=len(ff),distinct_observed_peak_count=len(peaks),charge_support_score=score,charge_plausibility_class=klass))
 return tuple(out)

def simulate_quality_tiers(records:Sequence[T1DeltaEvidenceQualityRecord]):
 scenarios=((SimulationScenarioID.CURRENT_METADATA,PolarityStatus.UNKNOWN,ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION),(SimulationScenarioID.NEGATIVE_POLARITY_CONFIRMED_ONLY,PolarityStatus.NEGATIVE,ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION),(SimulationScenarioID.POSITIVE_POLARITY_CONFIRMED_ONLY,PolarityStatus.POSITIVE,ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION),(SimulationScenarioID.MONOISOTOPIC_MZ_CONFIRMED_ONLY,PolarityStatus.UNKNOWN,ObservedMassDefinitionStatus.MONOISOTOPIC_MZ),(SimulationScenarioID.NEGATIVE_PLUS_MONO_CONFIRMED,PolarityStatus.NEGATIVE,ObservedMassDefinitionStatus.MONOISOTOPIC_MZ),(SimulationScenarioID.POSITIVE_PLUS_MONO_CONFIRMED,PolarityStatus.POSITIVE,ObservedMassDefinitionStatus.MONOISOTOPIC_MZ));out=[]
 for sid,pol,mass in scenarios:
  tiers=[];survivors=[]
  for r in records:
   if sid is SimulationScenarioID.CURRENT_METADATA:tiers.append(r.evidence_tier);survivors.append(r);continue
   pol_ok=pol is PolarityStatus.UNKNOWN or (pol is PolarityStatus.NEGATIVE and r.ion_mode is T1IonMode.NEGATIVE_DEPROTONATED) or (pol is PolarityStatus.POSITIVE and r.ion_mode is T1IonMode.POSITIVE_PROTONATED)
   mass_ok=mass is ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION or r.reference_mass_definition.value=="MONOISOTOPIC_DELTA"
   if not pol_ok or not mass_ok:tiers.append(EvidenceTier.TIER_E_BLOCKED);continue
   survivors.append(r)
   if r.evidence_tier is EvidenceTier.TIER_D_DIAGNOSTIC_ONLY:tiers.append(EvidenceTier.TIER_D_DIAGNOSTIC_ONLY);continue
   confirmed=sum((pol is not PolarityStatus.UNKNOWN,mass is not ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION));score=min(100,r.evidence_support_score+15*confirmed)
   if confirmed==2 and score>=60:tiers.append(EvidenceTier.TIER_A_HIGH_SUPPORT)
   elif confirmed>=1 and score>=40:tiers.append(EvidenceTier.TIER_B_MODERATE_SUPPORT)
   else:tiers.append(EvidenceTier.TIER_C_WEAK_SUPPORT)
  counts=Counter(tiers);cats=Counter(r.chemical_hypothesis_class.value for r in survivors);avg=sum(r.reference_mass_definition.value=="AVERAGE_DELTA" for r in survivors);mono=sum(r.reference_mass_definition.value=="MONOISOTOPIC_DELTA" for r in survivors);known=sum(r.chemical_hypothesis_class is ChemicalHypothesisClass.KNOWN_MODIFICATION_EQUIVALENT for r in survivors);conflicts=sum(r.apex_centroid_concordance_class.value=="REFERENCE_CATEGORY_CONFLICT" for r in survivors);concord=sum(r.apex_centroid_concordance_class.value in {"BOTH_STRICT_SAME_REFERENCE","APEX_STRICT_CENTROID_EXPLORATORY_SAME_REFERENCE"} for r in survivors)
  out.append(TierSimulationResult(scenario_id=sid,assumed_polarity=pol,assumed_mass_definition=mass,relation_count=len(records),tier_a_count=counts[EvidenceTier.TIER_A_HIGH_SUPPORT],tier_b_count=counts[EvidenceTier.TIER_B_MODERATE_SUPPORT],tier_c_count=counts[EvidenceTier.TIER_C_WEAK_SUPPORT],tier_d_count=counts[EvidenceTier.TIER_D_DIAGNOSTIC_ONLY],tier_e_count=counts[EvidenceTier.TIER_E_BLOCKED],surviving_chemical_relation_count=sum(r.chemical_hypothesis_class is not ChemicalHypothesisClass.KNOWN_MODIFICATION_EQUIVALENT for r in survivors),surviving_distinct_peak_count=len({r.observed_peak_id for r in survivors}),category_survival=tuple(sorted(cats.items())),multi_charge_concordant_count=sum(r.charge_support_class.value=="MULTI_CHARGE_CONCORDANT_SUPPORT" for r in survivors),average_reference_surviving_count=avg,monoisotopic_reference_surviving_count=mono,known_modification_diagnostic_count=known,reference_category_conflict_count=conflicts,apex_centroid_concordant_count=concord))
 return tuple(out)

def simulate_priority_candidates(records:Sequence[T1DeltaEvidenceQualityRecord],observed_peak_ids:Iterable[str]):
 targets=set(observed_peak_ids);scenario_map={x.scenario_id:(x.assumed_polarity,x.assumed_mass_definition) for x in simulate_quality_tiers(())};out=[]
 for r in records:
  if r.observed_peak_id not in targets:continue
  for sid,(pol,mass) in scenario_map.items():
   pol_block=pol is not PolarityStatus.UNKNOWN and not ((pol is PolarityStatus.NEGATIVE and r.ion_mode is T1IonMode.NEGATIVE_DEPROTONATED) or (pol is PolarityStatus.POSITIVE and r.ion_mode is T1IonMode.POSITIVE_PROTONATED));mass_block=mass is not ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION and r.reference_mass_definition.value!="MONOISOTOPIC_DELTA"
   if sid is SimulationScenarioID.CURRENT_METADATA:tier=r.evidence_tier;change=0
   elif pol_block or mass_block:tier=EvidenceTier.TIER_E_BLOCKED;change=0
   elif r.evidence_tier is EvidenceTier.TIER_D_DIAGNOSTIC_ONLY:tier=r.evidence_tier;change=15*sum((pol is not PolarityStatus.UNKNOWN,mass is not ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION))
   else:
    confirmed=sum((pol is not PolarityStatus.UNKNOWN,mass is not ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION));change=15*confirmed;score=min(100,r.evidence_support_score+change);tier=EvidenceTier.TIER_A_HIGH_SUPPORT if confirmed==2 and score>=60 else EvidenceTier.TIER_B_MODERATE_SUPPORT if confirmed>=1 and score>=40 else EvidenceTier.TIER_C_WEAK_SUPPORT
   out.append(PriorityCandidateScenarioStatus(observed_peak_id=r.observed_peak_id,t1_delta_relation_id=r.t1_delta_relation_id,scenario_id=sid,current_tier=r.evidence_tier,simulated_tier=tier,remains_candidate=tier is not EvidenceTier.TIER_E_BLOCKED,blocked_by_polarity=pol_block,blocked_by_mass_definition=mass_block,score_change=change))
 return tuple(sorted(out,key=lambda x:(x.observed_peak_id,x.t1_delta_relation_id,x.scenario_id.value)))

def audit_t1_source_metadata(loaded,source,measurement,peak_result,ions,fragment_matches,quality_records,*,project_config_path:Path|None=None,user_metadata:Mapping[str,str]|None=None,parameters:SourceMetadataAuditParameters|None=None):
 p=parameters or SourceMetadataAuditParameters();p.validate();inventory=build_metadata_source_inventory(Path(loaded.runtime_path),loaded,source,measurement,project_config_path=project_config_path,user_metadata=user_metadata);spacing=calculate_spacing_diagnostics(loaded.coordinates);intensity=calculate_intensity_diagnostics(loaded.intensities);families=build_isotope_families(peak_result.peaks,p);user=user_metadata or {}
 def resolve_user_field(name):
  rows=[MetadataFieldEvidence(name,str(user[name]).upper(),"USER_CONFIRMED_METADATA",MetadataTrustLevel.AUTHORITATIVE_SOURCE_SPECIFIC,True,True)] if name in user else []
  return resolve_metadata_field(rows),rows
 (acq_value,acq_confirmed,acq_conflicts,acq_status),acq_e=resolve_user_field("acquisition_polarity");(exp_value,exp_confirmed,exp_conflicts,exp_status),exp_e=resolve_user_field("exported_ion_polarity")
 acquisition=PolarityStatus(acq_value) if acq_value in PolarityStatus._value2member_map_ else PolarityStatus.UNKNOWN;exported=PolarityStatus(exp_value) if exp_value in PolarityStatus._value2member_map_ else PolarityStatus.UNKNOWN;conflicts=acq_conflicts+exp_conflicts
 if acquisition not in {PolarityStatus.UNKNOWN,PolarityStatus.NOT_APPLICABLE} and exported not in {PolarityStatus.UNKNOWN,PolarityStatus.NOT_APPLICABLE} and acquisition is not exported:
  final_pol=PolarityStatus.CONFLICT;pol_confirmed=False;confstatus=ConflictResolutionStatus.UNRESOLVED_CONFLICT
 else:
  final_pol=exported if exported is not PolarityStatus.UNKNOWN else acquisition;pol_confirmed=acq_confirmed and exp_confirmed and acquisition is exported;confstatus=ConflictResolutionStatus.UNRESOLVED_CONFLICT if acq_status is ConflictResolutionStatus.UNRESOLVED_CONFLICT or exp_status is ConflictResolutionStatus.UNRESOLVED_CONFLICT else ConflictResolutionStatus.NO_CONFLICT
 rep_value=str(user.get("data_representation",DataRepresentationStatus.UNKNOWN_REPRESENTATION.value)).upper();representation=DataRepresentationStatus(rep_value) if rep_value in DataRepresentationStatus._value2member_map_ else DataRepresentationStatus.UNKNOWN_REPRESENTATION;rep_confirmed="data_representation" in user and representation not in {DataRepresentationStatus.UNKNOWN_REPRESENTATION,DataRepresentationStatus.CONFLICTING_REPRESENTATION}
 mass_value=str(user.get("observed_mass_definition",ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION.value)).upper();mass=ObservedMassDefinitionStatus(mass_value) if mass_value in ObservedMassDefinitionStatus._value2member_map_ else ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION;mass_confirmed="observed_mass_definition" in user and mass not in {ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION,ObservedMassDefinitionStatus.CONFLICTING_MASS_DEFINITION}
 export_value=str(user.get("export_transformation",ExportTransformationStatus.UNKNOWN_EXPORT_TRANSFORMATION.value)).upper();export=ExportTransformationStatus(export_value) if export_value in ExportTransformationStatus._value2member_map_ else ExportTransformationStatus.UNKNOWN_EXPORT_TRANSFORMATION;export_confirmed="export_transformation" in user and export is not ExportTransformationStatus.UNKNOWN_EXPORT_TRANSFORMATION
 charges=build_charge_support_summaries(ions,fragment_matches,quality_records,families,polarity=final_pol,parameters=p);supported=tuple(x.charge for x in charges if x.charge_plausibility_class is ChargePlausibilityClass.SUPPORTED_CHARGE);weak=tuple(x.charge for x in charges if x.charge_plausibility_class is ChargePlausibilityClass.WEAKLY_SUPPORTED_CHARGE);unsupported=tuple(x.charge for x in charges if x.charge_plausibility_class is ChargePlausibilityClass.INSUFFICIENT_SUPPORT);raw=next(x for x in inventory if x.metadata_source_type is MetadataSourceType.MZML_OR_RAW_SOURCE_REFERENCE)
 if confstatus is ConflictResolutionStatus.UNRESOLVED_CONFLICT:completeness=MetadataCompletenessClass.CONFLICTING_METADATA
 elif pol_confirmed and mass_confirmed and rep_confirmed and export_confirmed:completeness=MetadataCompletenessClass.COMPLETE_FOR_CHEMICAL_ASSIGNMENT
 elif (pol_confirmed or mass_confirmed) and rep_confirmed:completeness=MetadataCompletenessClass.SUFFICIENT_FOR_TIER_B_SHADOW_EVIDENCE
 else:completeness=MetadataCompletenessClass.INSUFFICIENT_FOR_CHEMICAL_ASSIGNMENT
 recommend="KEEP_CONFIGURED_RANGE" if completeness in {MetadataCompletenessClass.INSUFFICIENT_FOR_CHEMICAL_ASSIGNMENT,MetadataCompletenessClass.CONFLICTING_METADATA} else ",".join(map(str,supported or weak));pol_e=acq_e+exp_e
 record=SourceMetadataAuditRecord(t1_source_metadata_audit_id="T1SOURCEAUDIT__"+sha256(f"{source.profile_source_id}|{loaded.source_sha256}".encode()).hexdigest()[:20].upper(),source_id=source.profile_source_id,measurement_id=source.measurement_id,rna_identity=source.rna_identity_id,runtime_path=str(loaded.runtime_path),sha256=loaded.source_sha256,coordinate_type="MZ",observed_mass_scale="MZ",observed_output_species="CHARGED_ION_UNKNOWN",acquisition_polarity=acquisition,exported_ion_polarity=exported,polarity_confirmed=pol_confirmed,polarity_evidence_source=tuple(e.source_id for e in pol_e),polarity_evidence_trust_level=max((e.trust_level for e in pol_e),key=lambda x:_TRUST_RANK[x],default=MetadataTrustLevel.NOT_AVAILABLE),polarity_evidence_count=len(pol_e),polarity_conflict_count=len(conflicts)+(1 if final_pol is PolarityStatus.CONFLICT else 0),polarity_final_status=final_pol,data_representation_status=representation,representation_confirmed=rep_confirmed,representation_evidence_sources=("USER_CONFIRMED_METADATA",) if rep_confirmed else ("EMPIRICAL_PROFILE_STRUCTURE",),representation_trust_level=MetadataTrustLevel.AUTHORITATIVE_SOURCE_SPECIFIC if rep_confirmed else MetadataTrustLevel.INFERRED_FROM_DATA,observed_mass_definition_status=mass,mass_definition_confirmed=mass_confirmed,mass_definition_evidence_source=("USER_CONFIRMED_METADATA",) if mass_confirmed else (),mass_definition_trust_level=MetadataTrustLevel.AUTHORITATIVE_SOURCE_SPECIFIC if mass_confirmed else MetadataTrustLevel.NOT_AVAILABLE,export_transformation_status=export,export_transformation_confirmed=export_confirmed,configured_charge_range=p.configured_charges,empirically_supported_charges=supported,weakly_supported_charges=weak,unsupported_charges=unsupported,recommended_shadow_charge_range=recommend,metadata_completeness_class=completeness,metadata_conflict_status=confstatus,source_linked_raw_metadata_available=raw.available)
 return T1SourceMetadataAuditResult(source_record=record,inventory=inventory,conflicts=conflicts,spacing_diagnostics=spacing,intensity_diagnostics=intensity,isotope_families=families,charge_summaries=charges,simulations=simulate_quality_tiers(quality_records))
