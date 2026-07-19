"""Bounded shadow audit of chemical neutral deltas in SCIEX RNase-T1 profiles."""
from __future__ import annotations
from bisect import bisect_left,bisect_right
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from hashlib import sha256
from math import isfinite
from statistics import median
from typing import Iterable
from rna_masshunter.intact_rna_average_mass import AVERAGE_ATOMIC_MASSES,calculate_average_neutral_mass_from_composition
from rna_masshunter.masses import MONOISOTOPIC_ATOMIC_MASSES
from rna_masshunter.sciex_intact_peak_family import DeltaMassDefinition
from rna_masshunter.sciex_t1_fragment_shadow_match import T1IonCandidate,T1IonMode
from rna_masshunter.sciex_t1_profile_peak_audit import T1ProfilePeak,T1Safeguards

class ChemicalReferenceCategory(str,Enum):
 OXYGEN_ADDITION_EQUIVALENT="OXYGEN_ADDITION_EQUIVALENT"; WATER_ADDITION_EQUIVALENT="WATER_ADDITION_EQUIVALENT"; WATER_LOSS_EQUIVALENT="WATER_LOSS_EQUIVALENT"; SULFUR_ADDITION_EQUIVALENT="SULFUR_ADDITION_EQUIVALENT"; O_TO_S_SUBSTITUTION_EQUIVALENT="O_TO_S_SUBSTITUTION_EQUIVALENT"; S_TO_O_SUBSTITUTION_EQUIVALENT="S_TO_O_SUBSTITUTION_EQUIVALENT"; KNOWN_MODIFICATION_DELTA="KNOWN_MODIFICATION_DELTA"; OTHER_DIAGNOSTIC="OTHER_DIAGNOSTIC"
class DeltaRelationClass(str,Enum):
 CHEMICAL_DELTA_STRICT="CHEMICAL_DELTA_STRICT"; CHEMICAL_DELTA_EXPLORATORY="CHEMICAL_DELTA_EXPLORATORY"; KNOWN_MODIFICATION_DELTA_STRICT="KNOWN_MODIFICATION_DELTA_STRICT"; KNOWN_MODIFICATION_DELTA_EXPLORATORY="KNOWN_MODIFICATION_DELTA_EXPLORATORY"; MASS_DEFINITION_MISMATCH_DIAGNOSTIC="MASS_DEFINITION_MISMATCH_DIAGNOSTIC"; NO_REFERENCE_MATCH="NO_REFERENCE_MATCH"
class ChemicalHypothesisClass(str,Enum):
 O_EQUIVALENT="O_EQUIVALENT"; H2O_EQUIVALENT="H2O_EQUIVALENT"; H2O_LOSS_EQUIVALENT="H2O_LOSS_EQUIVALENT"; S_EQUIVALENT="S_EQUIVALENT"; O_TO_S_EQUIVALENT="O_TO_S_EQUIVALENT"; S_TO_O_EQUIVALENT="S_TO_O_EQUIVALENT"; KNOWN_MODIFICATION_EQUIVALENT="KNOWN_MODIFICATION_EQUIVALENT"; OTHER="OTHER"
class ReferenceAmbiguityClass(str,Enum):
 UNIQUE_REFERENCE_CATEGORY="UNIQUE_REFERENCE_CATEGORY"; MULTIPLE_REFERENCE_SAME_CATEGORY="MULTIPLE_REFERENCE_SAME_CATEGORY"; MULTIPLE_REFERENCE_CATEGORIES="MULTIPLE_REFERENCE_CATEGORIES"; MASS_DEFINITION_AMBIGUITY="MASS_DEFINITION_AMBIGUITY"; NO_REFERENCE="NO_REFERENCE"
class StateEdgeClass(str,Enum): STRICT="STRICT"; EXPLORATORY="EXPLORATORY"
class T1StateSeriesPattern(str,Enum):
 SINGLE_O_STEP="SINGLE_O_STEP"; MULTIPLE_O_STEPS="MULTIPLE_O_STEPS"; SINGLE_H2O_STEP="SINGLE_H2O_STEP"; MULTIPLE_H2O_STEPS="MULTIPLE_H2O_STEPS"; SINGLE_O_TO_S_STEP="SINGLE_O_TO_S_STEP"; MIXED_O_H2O_SERIES="MIXED_O_H2O_SERIES"; MIXED_O_S_SERIES="MIXED_O_S_SERIES"; BRANCHED_SERIES="BRANCHED_SERIES"; UNRESOLVED_SERIES="UNRESOLVED_SERIES"
@dataclass(frozen=True)
class ChemicalDeltaAuditParameters:
 strict_tolerance_mz:float=.01; exploratory_tolerance_mz:float=.02; strict_reference_delta_tolerance_da:float=.05; exploratory_reference_delta_tolerance_da:float=.10
 minimum_neutral_delta_da:float=-250.; maximum_neutral_delta_da:float=300.; maximum_theoretical_candidates_per_observed_peak_per_charge:int=50; maximum_reported_delta_relations_per_peak:int=20
 def validate(self):
  if self.strict_tolerance_mz<=0 or self.exploratory_tolerance_mz<self.strict_tolerance_mz:raise ValueError("invalid m/z tolerances")
  if self.strict_reference_delta_tolerance_da<=0 or self.exploratory_reference_delta_tolerance_da<self.strict_reference_delta_tolerance_da:raise ValueError("invalid reference tolerances")
  if self.minimum_neutral_delta_da>=self.maximum_neutral_delta_da:raise ValueError("invalid neutral delta range")
  if self.maximum_theoretical_candidates_per_observed_peak_per_charge<1 or self.maximum_reported_delta_relations_per_peak<1:raise ValueError("candidate bounds must be positive")
@dataclass(frozen=True)
class ChemicalDeltaReference:
 reference_id:str; reference_name:str; reference_category:ChemicalReferenceCategory; signed_delta_da:float; reference_mass_definition:DeltaMassDefinition; reference_provenance:str; chemical_interpretation:str; diagnostic_only:bool=True; comparison_role:str="CHEMICAL_DELTA_DIAGNOSTIC"
@dataclass(frozen=True,kw_only=True)
class DeltaSafeguards(T1Safeguards):
 modification_position_assigned:bool=False; exact_residue_assigned:bool=False; oxygen_addition_assigned:bool=False; hydration_assigned:bool=False; dehydration_assigned:bool=False
 sulfur_addition_assigned:bool=False; o_to_s_substitution_assigned:bool=False; s_to_o_substitution_assigned:bool=False; sulfur_atom_assigned:bool=False
 thioamide_assigned:bool=False; thioamide_position_assigned:bool=False; thioamide_oxidation_state_assigned:bool=False
 oxidation_assigned:bool=False; oxidation_number_assigned:bool=False; oxidation_pathway_assigned:bool=False; observed_mass_definition_confirmed:bool=False
 reaction_direction_assigned:bool=False; precursor_product_assigned:bool=False; in_source_water_loss_excluded:bool=False; adduct_explanation_excluded:bool=False
@dataclass(frozen=True,kw_only=True)
class T1ChemicalDeltaRelation(DeltaSafeguards):
 t1_delta_relation_id:str; observed_source_id:str; observed_peak_id:str; observed_apex_mz:float; observed_centroid_mz:float|None; rna_identity_candidate:str
 theoretical_fragment_id:str; fragment_sequence:str; start_position:int; end_position:int; cca_state:str; ion_mode:T1IonMode; charge:int; theoretical_mz:float
 apex_mz_delta:float; centroid_mz_delta:float|None; apex_neutral_delta:float; centroid_neutral_delta:float|None
 reference_id:str; reference_name:str; reference_category:ChemicalReferenceCategory; reference_delta_da:float; reference_mass_definition:DeltaMassDefinition; comparison_role:str
 apex_delta_error_da:float|None; centroid_delta_error_da:float|None; delta_relation_class:DeltaRelationClass; reference_match_class:str; chemical_hypothesis_class:ChemicalHypothesisClass
 observed_theoretical_mass_definition_compatibility:str; mass_definition_compatible:bool; ion_convention_compatible:bool; charge_convention_compatible:bool; polarity_compatible:bool
 possible_isotope_component:bool; possible_shoulder:bool; possible_duplicate:bool; eligible_for_neutral_delta_audit:bool; eligible_for_chemical_delta_evidence:bool
 candidate_reference_count:int; distinct_reference_category_count:int; best_reference_id:str; best_absolute_error_da:float|None; reference_ambiguity_class:ReferenceAmbiguityClass
 thioamide_hypothesis_possible:bool; hydration_or_dehydration_hypothesis_possible:bool
@dataclass(frozen=True,kw_only=True)
class T1ChemicalStateEdge(DeltaSafeguards):
 t1_chemical_state_edge_id:str; lower_observed_peak_id:str; higher_observed_peak_id:str; shared_theoretical_fragment_id:str; shared_charge:int; shared_ion_mode:T1IonMode
 observed_neutral_interpeak_delta:float; reference_name:str; reference_category:ChemicalReferenceCategory; reference_delta:float; error_da:float; edge_class:StateEdgeClass; eligible_for_state_series:bool
@dataclass(frozen=True,kw_only=True)
class T1ChemicalStateSeries(DeltaSafeguards):
 t1_state_series_id:str; rna_identity_candidate:str; theoretical_fragment_identity:str; member_observed_peak_ids:tuple[str,...]; member_apex_mzs:tuple[float,...]; member_centroid_mzs:tuple[float|None,...]; member_count:int; charge:int; ion_mode:T1IonMode
 o_edge_count:int; h2o_edge_count:int; s_edge_count:int; o_to_s_edge_count:int; s_to_o_edge_count:int; strict_edge_count:int; exploratory_edge_count:int; series_pattern:T1StateSeriesPattern; mass_span_neutral_da:float
 sequential_oxygen_equivalent_series_detected:bool; oxidation_state_series_possible:bool
@dataclass(frozen=True,kw_only=True)
class CandidateExplosionSummary(DeltaSafeguards):
 raw_theoretical_combination_count:int; post_mz_window_candidate_count:int; post_neutral_delta_range_candidate_count:int; post_per_peak_cap_candidate_count:int; final_reported_relation_count:int; maximum_candidates_for_one_observed_peak:int; median_candidates_per_observed_peak:float
@dataclass(frozen=True,kw_only=True)
class T1FragmentDeltaAuditResult(DeltaSafeguards):
 source_id:str; parameters:ChemicalDeltaAuditParameters; references:tuple[ChemicalDeltaReference,...]; primary_peak_count:int; relations:tuple[T1ChemicalDeltaRelation,...]; state_edges:tuple[T1ChemicalStateEdge,...]; state_series:tuple[T1ChemicalStateSeries,...]; candidate_explosion:CandidateExplosionSummary
@dataclass(frozen=True,kw_only=True)
class CrossSampleCategoryComparison(DeltaSafeguards):
 reference_category:ChemicalReferenceCategory; uaa_strict_peak_count:int; uaa_exploratory_peak_count:int; uag_strict_peak_count:int; uag_exploratory_peak_count:int; detected_in_both:bool; detected_only_in_uaa:bool; detected_only_in_uag:bool; uaa_fraction_of_primary_selected_peaks:float; uag_fraction_of_primary_selected_peaks:float; uaa_fraction_of_all_eligible_relations:float; uag_fraction_of_all_eligible_relations:float
@dataclass(frozen=True,kw_only=True)
class CrossSampleFragmentDeltaComparison(DeltaSafeguards):
 fragment_identity:str; chemical_hypothesis_class:ChemicalHypothesisClass; uaa_observed_support:int; uag_observed_support:int; shared_delta_hypothesis:bool; uaa_only_delta_hypothesis:bool; uag_only_delta_hypothesis:bool; cross_rna_sequence_identity_ambiguous:bool

def build_chemical_delta_reference_registry(modifications:Iterable[object]=()):
 ah,ao,ass=AVERAGE_ATOMIC_MASSES["H"],AVERAGE_ATOMIC_MASSES["O"],AVERAGE_ATOMIC_MASSES["S"]
 mh,mo,ms=MONOISOTOPIC_ATOMIC_MASSES["H"],MONOISOTOPIC_ATOMIC_MASSES["O"],MONOISOTOPIC_ATOMIC_MASSES["S"]
 aw=calculate_average_neutral_mass_from_composition({"H":2,"O":1});mw=2*mh+mo
 refs=[]
 def add(prefix,name,cat,value,definition,provenance,interpretation,comparison_role="CHEMICAL_DELTA_DIAGNOSTIC"):refs.append(ChemicalDeltaReference(prefix,name,cat,float(value),definition,provenance,interpretation,comparison_role=comparison_role))
 for suffix,definition,prov,o,w,s in (("AVERAGE",DeltaMassDefinition.AVERAGE_DELTA,"AVERAGE_ATOMIC_MASSES",ao,aw,ass),("MONOISOTOPIC",DeltaMassDefinition.MONOISOTOPIC_DELTA,"MONOISOTOPIC_ATOMIC_MASSES",mo,mw,ms)):
  add(f"O_ADDITION_{suffix}",f"O_ADDITION_{suffix}",ChemicalReferenceCategory.OXYGEN_ADDITION_EQUIVALENT,o,definition,prov,"OXYGEN_ADDITION_EQUIVALENT")
  add(f"H2O_ADDITION_{suffix}",f"H2O_ADDITION_{suffix}",ChemicalReferenceCategory.WATER_ADDITION_EQUIVALENT,w,definition,prov,"WATER_ADDITION_EQUIVALENT")
  add(f"H2O_LOSS_{suffix}",f"H2O_LOSS_{suffix}",ChemicalReferenceCategory.WATER_LOSS_EQUIVALENT,-w,definition,prov,"WATER_LOSS_EQUIVALENT")
  add(f"S_ADDITION_{suffix}",f"S_ADDITION_{suffix}",ChemicalReferenceCategory.SULFUR_ADDITION_EQUIVALENT,s,definition,prov,"SULFUR_ADDITION_EQUIVALENT")
  add(f"O_TO_S_{suffix}",f"O_TO_S_{suffix}",ChemicalReferenceCategory.O_TO_S_SUBSTITUTION_EQUIVALENT,s-o,definition,prov,"O_TO_S_SUBSTITUTION_EQUIVALENT")
  add(f"S_TO_O_{suffix}",f"S_TO_O_{suffix}",ChemicalReferenceCategory.S_TO_O_SUBSTITUTION_EQUIVALENT,o-s,definition,prov,"S_TO_O_SUBSTITUTION_EQUIVALENT")
 for mod in modifications:
  value=getattr(mod,"mass_shift_from_unmodified",float("nan"));name=str(getattr(mod,"id","") or getattr(mod,"symbol","") or "")
  if name and isfinite(float(value)) and float(value)!=0:add("KNOWN_MOD__"+sha256(name.encode()).hexdigest()[:16].upper(),name,ChemicalReferenceCategory.KNOWN_MODIFICATION_DELTA,float(value),DeltaMassDefinition.MONOISOTOPIC_DELTA,"data/modifications.yaml:Mass(mono)","KNOWN_MODIFICATION_MASS_SHIFT","MASS_DEFINITION_MISMATCH_DIAGNOSTIC")
 return tuple(refs)
def convert_mz_delta_to_neutral_delta(observed_mz,theoretical_mz,charge,*,observed_ion_mode=None,theoretical_ion_mode=None):
 z=abs(int(charge))
 if z==0:raise ValueError("charge cannot be zero")
 compatible=observed_ion_mode is None or observed_ion_mode==theoretical_ion_mode
 delta=float(observed_mz)-float(theoretical_mz)
 return delta,delta*z,compatible
_HYP={ChemicalReferenceCategory.OXYGEN_ADDITION_EQUIVALENT:ChemicalHypothesisClass.O_EQUIVALENT,ChemicalReferenceCategory.WATER_ADDITION_EQUIVALENT:ChemicalHypothesisClass.H2O_EQUIVALENT,ChemicalReferenceCategory.WATER_LOSS_EQUIVALENT:ChemicalHypothesisClass.H2O_LOSS_EQUIVALENT,ChemicalReferenceCategory.SULFUR_ADDITION_EQUIVALENT:ChemicalHypothesisClass.S_EQUIVALENT,ChemicalReferenceCategory.O_TO_S_SUBSTITUTION_EQUIVALENT:ChemicalHypothesisClass.O_TO_S_EQUIVALENT,ChemicalReferenceCategory.S_TO_O_SUBSTITUTION_EQUIVALENT:ChemicalHypothesisClass.S_TO_O_EQUIVALENT,ChemicalReferenceCategory.KNOWN_MODIFICATION_DELTA:ChemicalHypothesisClass.KNOWN_MODIFICATION_EQUIVALENT,ChemicalReferenceCategory.OTHER_DIAGNOSTIC:ChemicalHypothesisClass.OTHER}
def _ref_matches(delta,z,refs,p):
 strict=max(p.strict_reference_delta_tolerance_da,p.strict_tolerance_mz*z);expl=max(p.exploratory_reference_delta_tolerance_da,p.exploratory_tolerance_mz*z);out=[]
 for r in refs:
  e=delta-r.signed_delta_da
  if abs(e)<=expl:out.append((abs(e),"STRICT" if abs(e)<=strict else "EXPLORATORY",r,e))
 return sorted(out,key=lambda x:(x[0],x[2].reference_category.value,x[2].reference_id))
def _ambiguity(matches):
 if not matches:return ReferenceAmbiguityClass.NO_REFERENCE
 cats={m[2].reference_category for m in matches};defs={m[2].reference_mass_definition for m in matches}
 if len(cats)>1:return ReferenceAmbiguityClass.MULTIPLE_REFERENCE_CATEGORIES
 if len(defs)>1:return ReferenceAmbiguityClass.MASS_DEFINITION_AMBIGUITY
 if len(matches)>1:return ReferenceAmbiguityClass.MULTIPLE_REFERENCE_SAME_CATEGORY
 return ReferenceAmbiguityClass.UNIQUE_REFERENCE_CATEGORY
def _collapse_ions(ions):
 groups={}
 for c in ions:
  key=(c.rna_identity,c.fragment_sequence,c.start_position,c.end_position,c.ion_mode,c.charge,round(c.theoretical_mz,9))
  groups.setdefault(key,[]).append(c)
 out=[]
 for key,items in groups.items():
  c=min(items,key=lambda x:x.ion_candidate_id);states="|".join(sorted({x.cca_state for x in items}));fid=f"{c.rna_identity}:{c.start_position}-{c.end_position}:{c.fragment_sequence}"
  out.append((c,states,fid))
 return tuple(sorted(out,key=lambda x:(x[0].theoretical_mz,x[2],x[1])))

def audit_t1_fragment_chemical_deltas(peaks,ion_candidates,*,source_id,references=None,modifications=(),parameters=None):
 p=parameters or ChemicalDeltaAuditParameters();p.validate();refs=tuple(references or build_chemical_delta_reference_registry(modifications));primary=tuple(sorted((x for x in peaks if x.selected_as_primary and not x.possible_isotope_component and not x.possible_shoulder and not x.possible_duplicate),key=lambda x:(x.apex_mz,x.t1_peak_id)))
 collapsed=_collapse_ions(ion_candidates);by_charge={z:[] for z in range(1,6)}
 for item in collapsed:by_charge.setdefault(item[0].charge,[]).append(item)
 masses={z:[x[0].theoretical_mz for x in items] for z,items in by_charge.items()};raw=len(primary)*len(ion_candidates);window_count=range_count=cap_count=0;relations=[];per_peak=[]
 for peak in primary:
  candidates=[]
  for z,items in by_charge.items():
   if not items:continue
   low=peak.apex_mz-p.maximum_neutral_delta_da/z;high=peak.apex_mz-p.minimum_neutral_delta_da/z;lo=bisect_left(masses[z],low);hi=bisect_right(masses[z],high);found=items[lo:hi];window_count+=len(found);valid=[]
   for item in found:
    c,states,fid=item;_,nd,_=convert_mz_delta_to_neutral_delta(peak.apex_mz,c.theoretical_mz,z); 
    if p.minimum_neutral_delta_da<=nd<=p.maximum_neutral_delta_da:valid.append((item,nd,_ref_matches(nd,z,refs,p)))
   range_count+=len(valid);valid.sort(key=lambda x:((x[2][0][0] if x[2] else float("inf")),abs(x[1]),x[0][2],x[0][0].ion_candidate_id));valid=valid[:p.maximum_theoretical_candidates_per_observed_peak_per_charge];cap_count+=len(valid);candidates.extend(valid)
  rows=[]
  for (c,states,fid),nd,matches in candidates:
   mz_delta=peak.apex_mz-c.theoretical_mz;cent_mz=peak.centroid_mz-c.theoretical_mz if peak.centroid_mz is not None else None;cent_nd=cent_mz*abs(c.charge) if cent_mz is not None else None;ambiguity=_ambiguity(matches)
   if matches:_,tier,ref,error=matches[0];cent_error=cent_nd-ref.signed_delta_da if cent_nd is not None else None;known=ref.reference_category is ChemicalReferenceCategory.KNOWN_MODIFICATION_DELTA;klass=DeltaRelationClass.MASS_DEFINITION_MISMATCH_DIAGNOSTIC if known else DeltaRelationClass.CHEMICAL_DELTA_STRICT if tier=="STRICT" else DeltaRelationClass.CHEMICAL_DELTA_EXPLORATORY
   else:tier="NO_MATCH";ref=ChemicalDeltaReference("NO_REFERENCE","NO_REFERENCE",ChemicalReferenceCategory.OTHER_DIAGNOSTIC,0,DeltaMassDefinition.UNKNOWN,"NONE","NO_REFERENCE");error=cent_error=None;klass=DeltaRelationClass.NO_REFERENCE_MATCH
   eligible=bool(matches and not known);hyp=_HYP[ref.reference_category]
   rows.append(T1ChemicalDeltaRelation(t1_delta_relation_id="T1DELTA__"+sha256(f"{peak.t1_peak_id}|{fid}|{c.ion_mode.value}|{c.charge}|{ref.reference_id}".encode()).hexdigest()[:20].upper(),observed_source_id=source_id,observed_peak_id=peak.t1_peak_id,observed_apex_mz=peak.apex_mz,observed_centroid_mz=peak.centroid_mz,rna_identity_candidate=c.rna_identity,theoretical_fragment_id=fid,fragment_sequence=c.fragment_sequence,start_position=c.start_position,end_position=c.end_position,cca_state="COLLAPSED:"+states,ion_mode=c.ion_mode,charge=c.charge,theoretical_mz=c.theoretical_mz,apex_mz_delta=mz_delta,centroid_mz_delta=cent_mz,apex_neutral_delta=nd,centroid_neutral_delta=cent_nd,reference_id=ref.reference_id,reference_name=ref.reference_name,reference_category=ref.reference_category,reference_delta_da=ref.signed_delta_da,reference_mass_definition=ref.reference_mass_definition,comparison_role=ref.comparison_role,apex_delta_error_da=error,centroid_delta_error_da=cent_error,delta_relation_class=klass,reference_match_class=tier,chemical_hypothesis_class=hyp,observed_theoretical_mass_definition_compatibility="UNKNOWN",mass_definition_compatible=False,ion_convention_compatible=True,charge_convention_compatible=True,polarity_compatible=False,possible_isotope_component=peak.possible_isotope_component,possible_shoulder=peak.possible_shoulder,possible_duplicate=peak.possible_duplicate,eligible_for_neutral_delta_audit=True,eligible_for_chemical_delta_evidence=eligible,candidate_reference_count=len(matches),distinct_reference_category_count=len({x[2].reference_category for x in matches}),best_reference_id=ref.reference_id,best_absolute_error_da=abs(error) if error is not None else None,reference_ambiguity_class=ambiguity,thioamide_hypothesis_possible=hyp in {ChemicalHypothesisClass.O_TO_S_EQUIVALENT,ChemicalHypothesisClass.S_TO_O_EQUIVALENT,ChemicalHypothesisClass.S_EQUIVALENT},hydration_or_dehydration_hypothesis_possible=hyp in {ChemicalHypothesisClass.H2O_EQUIVALENT,ChemicalHypothesisClass.H2O_LOSS_EQUIVALENT}))
  rows.sort(key=lambda r:(r.delta_relation_class is DeltaRelationClass.NO_REFERENCE_MATCH,r.best_absolute_error_da if r.best_absolute_error_da is not None else float("inf"),r.theoretical_fragment_id,r.ion_mode.value,r.charge));kept=rows[:p.maximum_reported_delta_relations_per_peak];relations.extend(kept);per_peak.append(len(candidates))
 edges=_build_edges(primary,tuple(relations),refs,p);series=_build_series(primary,edges,relations);summary=CandidateExplosionSummary(raw_theoretical_combination_count=raw,post_mz_window_candidate_count=window_count,post_neutral_delta_range_candidate_count=range_count,post_per_peak_cap_candidate_count=cap_count,final_reported_relation_count=len(relations),maximum_candidates_for_one_observed_peak=max(per_peak,default=0),median_candidates_per_observed_peak=median(per_peak) if per_peak else 0.)
 return T1FragmentDeltaAuditResult(source_id=source_id,parameters=p,references=refs,primary_peak_count=len(primary),relations=tuple(relations),state_edges=edges,state_series=series,candidate_explosion=summary)

def _build_edges(peaks,relations,refs,p):
 by_peak={x.t1_peak_id:x for x in peaks};groups=defaultdict(set)
 for r in relations:
  if r.eligible_for_neutral_delta_audit:groups[(r.rna_identity_candidate,r.theoretical_fragment_id,r.charge,r.ion_mode)].add(r.observed_peak_id)
 chemical=[r for r in refs if r.reference_category is not ChemicalReferenceCategory.KNOWN_MODIFICATION_DELTA and r.signed_delta_da>0]
 edges=[];seen=set()
 for (rna,fid,z,mode),ids in sorted(groups.items(),key=lambda x:str(x[0])):
  ordered=sorted(ids,key=lambda i:(by_peak[i].apex_mz,i))
  for i,left_id in enumerate(ordered):
   for right_id in ordered[i+1:]:
    left,right=by_peak[left_id],by_peak[right_id];delta=(right.apex_mz-left.apex_mz)*z;matches=_ref_matches(delta,z,chemical,p)
    if not matches:continue
    _,tier,ref,error=matches[0];key=(left_id,right_id,fid,z,mode,ref.reference_category)
    if key in seen:continue
    seen.add(key);edges.append(T1ChemicalStateEdge(t1_chemical_state_edge_id="T1STATEEDGE__"+sha256("|".join(map(str,key)).encode()).hexdigest()[:20].upper(),lower_observed_peak_id=left_id,higher_observed_peak_id=right_id,shared_theoretical_fragment_id=fid,shared_charge=z,shared_ion_mode=mode,observed_neutral_interpeak_delta=delta,reference_name=ref.reference_name,reference_category=ref.reference_category,reference_delta=ref.signed_delta_da,error_da=error,edge_class=StateEdgeClass(tier),eligible_for_state_series=True))
 return tuple(sorted(edges,key=lambda e:(e.shared_theoretical_fragment_id,e.shared_charge,e.shared_ion_mode.value,by_peak[e.lower_observed_peak_id].apex_mz,by_peak[e.higher_observed_peak_id].apex_mz,e.reference_name)))
def _series_pattern(edges,degrees):
 cats={e.reference_category for e in edges}
 if any(v>2 for v in degrees.values()):return T1StateSeriesPattern.BRANCHED_SERIES
 o=ChemicalReferenceCategory.OXYGEN_ADDITION_EQUIVALENT in cats;w=bool(cats&{ChemicalReferenceCategory.WATER_ADDITION_EQUIVALENT,ChemicalReferenceCategory.WATER_LOSS_EQUIVALENT});s=bool(cats&{ChemicalReferenceCategory.SULFUR_ADDITION_EQUIVALENT,ChemicalReferenceCategory.O_TO_S_SUBSTITUTION_EQUIVALENT,ChemicalReferenceCategory.S_TO_O_SUBSTITUTION_EQUIVALENT})
 if o and w:return T1StateSeriesPattern.MIXED_O_H2O_SERIES
 if o and s:return T1StateSeriesPattern.MIXED_O_S_SERIES
 if len(edges)==1 and o:return T1StateSeriesPattern.SINGLE_O_STEP
 if len(edges)>1 and o:return T1StateSeriesPattern.MULTIPLE_O_STEPS
 if len(edges)==1 and w:return T1StateSeriesPattern.SINGLE_H2O_STEP
 if len(edges)>1 and w:return T1StateSeriesPattern.MULTIPLE_H2O_STEPS
 if len(edges)==1 and ChemicalReferenceCategory.O_TO_S_SUBSTITUTION_EQUIVALENT in cats:return T1StateSeriesPattern.SINGLE_O_TO_S_STEP
 return T1StateSeriesPattern.UNRESOLVED_SERIES
def _build_series(peaks,edges,relations):
 by_peak={x.t1_peak_id:x for x in peaks};rna_by_fragment={r.theoretical_fragment_id:r.rna_identity_candidate for r in relations};groups=defaultdict(list)
 for e in edges:groups[(e.shared_theoretical_fragment_id,e.shared_charge,e.shared_ion_mode)].append(e)
 out=[]
 for (fid,z,mode),items in groups.items():
  adjacency=defaultdict(set)
  for e in items:adjacency[e.lower_observed_peak_id].add(e.higher_observed_peak_id);adjacency[e.higher_observed_peak_id].add(e.lower_observed_peak_id)
  unseen=set(adjacency)
  while unseen:
   start=min(unseen,key=lambda i:(by_peak[i].apex_mz,i));stack=[start];members=set()
   while stack:
    cur=stack.pop()
    if cur in members:continue
    members.add(cur);unseen.discard(cur);stack.extend(adjacency[cur]-members)
   ordered=tuple(sorted(members,key=lambda i:(by_peak[i].apex_mz,i)));component=tuple(e for e in items if e.lower_observed_peak_id in members and e.higher_observed_peak_id in members);degrees={i:0 for i in members}
   for e in component:degrees[e.lower_observed_peak_id]+=1;degrees[e.higher_observed_peak_id]+=1
   count=lambda cat:sum(e.reference_category is cat for e in component);o=count(ChemicalReferenceCategory.OXYGEN_ADDITION_EQUIVALENT);h=count(ChemicalReferenceCategory.WATER_ADDITION_EQUIVALENT)+count(ChemicalReferenceCategory.WATER_LOSS_EQUIVALENT);s=count(ChemicalReferenceCategory.SULFUR_ADDITION_EQUIVALENT);os=count(ChemicalReferenceCategory.O_TO_S_SUBSTITUTION_EQUIVALENT);so=count(ChemicalReferenceCategory.S_TO_O_SUBSTITUTION_EQUIVALENT);span=(by_peak[ordered[-1]].apex_mz-by_peak[ordered[0]].apex_mz)*z
   out.append(T1ChemicalStateSeries(t1_state_series_id="T1STATESERIES__"+sha256(f"{fid}|{z}|{mode.value}|{'|'.join(ordered)}".encode()).hexdigest()[:20].upper(),rna_identity_candidate=rna_by_fragment.get(fid,"UNKNOWN"),theoretical_fragment_identity=fid,member_observed_peak_ids=ordered,member_apex_mzs=tuple(by_peak[i].apex_mz for i in ordered),member_centroid_mzs=tuple(by_peak[i].centroid_mz for i in ordered),member_count=len(ordered),charge=z,ion_mode=mode,o_edge_count=o,h2o_edge_count=h,s_edge_count=s,o_to_s_edge_count=os,s_to_o_edge_count=so,strict_edge_count=sum(e.edge_class is StateEdgeClass.STRICT for e in component),exploratory_edge_count=sum(e.edge_class is StateEdgeClass.EXPLORATORY for e in component),series_pattern=_series_pattern(component,degrees),mass_span_neutral_da=span,sequential_oxygen_equivalent_series_detected=o>1,oxidation_state_series_possible=o>1))
 return tuple(sorted(out,key=lambda s:(-s.member_count,-s.strict_edge_count,s.theoretical_fragment_identity,s.charge,s.ion_mode.value,s.t1_state_series_id)))
def compare_cross_sample_categories(uaa,uag):
 cats=tuple(ChemicalReferenceCategory);out=[]
 def counts(result,cat):
  rows=[r for r in result.relations if r.reference_category is cat and r.eligible_for_chemical_delta_evidence];strict={r.observed_peak_id for r in rows if r.reference_match_class=="STRICT"};expl={r.observed_peak_id for r in rows if r.reference_match_class=="EXPLORATORY"};return rows,strict,expl
 for cat in cats:
  ur,us,ue=counts(uaa,cat);gr,gs,ge=counts(uag,cat);ud=bool(us or ue);gd=bool(gs or ge)
  out.append(CrossSampleCategoryComparison(reference_category=cat,uaa_strict_peak_count=len(us),uaa_exploratory_peak_count=len(ue),uag_strict_peak_count=len(gs),uag_exploratory_peak_count=len(ge),detected_in_both=ud and gd,detected_only_in_uaa=ud and not gd,detected_only_in_uag=gd and not ud,uaa_fraction_of_primary_selected_peaks=(len(us|ue)/uaa.primary_peak_count if uaa.primary_peak_count else 0),uag_fraction_of_primary_selected_peaks=(len(gs|ge)/uag.primary_peak_count if uag.primary_peak_count else 0),uaa_fraction_of_all_eligible_relations=(len(ur)/max(1,sum(r.eligible_for_chemical_delta_evidence for r in uaa.relations))),uag_fraction_of_all_eligible_relations=(len(gr)/max(1,sum(r.eligible_for_chemical_delta_evidence for r in uag.relations)))))
 return tuple(out)
def compare_cross_sample_fragment_deltas(uaa,uag):
 def support(result):
  d=defaultdict(set)
  for r in result.relations:
   if r.eligible_for_chemical_delta_evidence:d[(r.theoretical_fragment_id,r.chemical_hypothesis_class)].add(r.observed_peak_id)
  return d
 u,g=support(uaa),support(uag);out=[];keys=set(u)|set(g);identities=defaultdict(set)
 for fid,hyp in keys:
  prefix,tail=fid.split(":",1);identities[(tail,hyp)].add(prefix)
 for key in sorted(keys,key=lambda x:(x[0],x[1].value)):
  us,gs=len(u.get(key,set())),len(g.get(key,set()));fid,hyp=key;tail=fid.split(":",1)[1]
  out.append(CrossSampleFragmentDeltaComparison(fragment_identity=fid,chemical_hypothesis_class=hyp,uaa_observed_support=us,uag_observed_support=gs,shared_delta_hypothesis=bool(us and gs),uaa_only_delta_hypothesis=bool(us and not gs),uag_only_delta_hypothesis=bool(gs and not us),cross_rna_sequence_identity_ambiguous=len(identities[(tail,hyp)])>1))
 return tuple(out)
