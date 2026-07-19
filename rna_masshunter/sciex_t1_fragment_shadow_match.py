"""RNase-T1 theoretical-fragment and observed m/z shadow matching."""
from __future__ import annotations
from bisect import bisect_left,bisect_right
from dataclasses import dataclass
from enum import Enum
from statistics import median
from rna_masshunter.cca_tail_state import generate_cca_tail_variants
from rna_masshunter.digestion import digest_sequence
from rna_masshunter.intact_rna_average_mass import calculate_intact_rna_average_mass
from rna_masshunter.intact_rna_mass import IntactRnaMassParameters
from rna_masshunter.masses import PROTON_MASS
from rna_masshunter.models import RunConfig
from rna_masshunter.ms1_mapping import theoretical_mz_from_mass
from rna_masshunter.sciex_sample_manifest import get_rna_identity
from rna_masshunter.sciex_t1_profile_peak_audit import T1ProfilePeak,T1Safeguards

class T1IonMode(str,Enum): NEGATIVE_DEPROTONATED="NEGATIVE_DEPROTONATED"; POSITIVE_PROTONATED="POSITIVE_PROTONATED"; UNKNOWN_POLARITY_DIAGNOSTIC="UNKNOWN_POLARITY_DIAGNOSTIC"
class T1FragmentMatchClass(str,Enum): STRICT="STRICT"; EXPLORATORY="EXPLORATORY"; NO_MATCH="NO_MATCH"
class T1FragmentAmbiguityClass(str,Enum):
 UNIQUE_FRAGMENT_CANDIDATE="UNIQUE_FRAGMENT_CANDIDATE"; MULTIPLE_CHARGE_AMBIGUITY="MULTIPLE_CHARGE_AMBIGUITY"; MULTIPLE_FRAGMENT_AMBIGUITY="MULTIPLE_FRAGMENT_AMBIGUITY"; MULTIPLE_CCA_STATE_AMBIGUITY="MULTIPLE_CCA_STATE_AMBIGUITY"; CROSS_RNA_IDENTITY_AMBIGUITY="CROSS_RNA_IDENTITY_AMBIGUITY"; MASS_DEFINITION_AMBIGUITY="MASS_DEFINITION_AMBIGUITY"; NO_MATCH="NO_MATCH"
class TheoreticalDiscriminationClass(str,Enum):
 SEQUENCE_IDENTICAL_FRAGMENT="SEQUENCE_IDENTICAL_FRAGMENT"; UAA_SPECIFIC_THEORETICAL_FRAGMENT="UAA_SPECIFIC_THEORETICAL_FRAGMENT"; UAG_SPECIFIC_THEORETICAL_FRAGMENT="UAG_SPECIFIC_THEORETICAL_FRAGMENT"; SAME_MASS_DIFFERENT_SEQUENCE="SAME_MASS_DIFFERENT_SEQUENCE"; SHARED_MASS_AMBIGUOUS_FRAGMENT="SHARED_MASS_AMBIGUOUS_FRAGMENT"
class IdentitySupportLevel(str,Enum): NO_DISCRIMINATORY_SUPPORT="NO_DISCRIMINATORY_SUPPORT"; WEAK_DISCRIMINATORY_SUPPORT="WEAK_DISCRIMINATORY_SUPPORT"; MODERATE_DISCRIMINATORY_SUPPORT="MODERATE_DISCRIMINATORY_SUPPORT"; STRONG_DISCRIMINATORY_SUPPORT="STRONG_DISCRIMINATORY_SUPPORT"
@dataclass(frozen=True)
class T1FragmentMatchParameters:
 strict_tolerance_mz:float=.01; exploratory_tolerance_mz:float=.02; charges:tuple[int,...]=(1,2,3,4,5)
 def validate(self):
  if self.strict_tolerance_mz<=0 or self.exploratory_tolerance_mz<self.strict_tolerance_mz:raise ValueError("invalid tolerance")
  if not self.charges or any(z<1 for z in self.charges):raise ValueError("positive charge magnitudes required")
@dataclass(frozen=True)
class TheoreticalT1Fragment:
 theoretical_t1_fragment_id:str; rna_identity:str; sequence_candidate_id:str; cca_state:str; five_prime_state:str; three_prime_state:str
 fragment_index:int; start_position:int; end_position:int; fragment_sequence:str; cleavage_context:str; neutral_monoisotopic_mass:float; neutral_average_mass:float
 native_modifications_expected:bool=True; unmodified_fragment_reference_only:bool=True; modified_fragment_composition_not_enumerated:bool=True
@dataclass(frozen=True)
class T1IonCandidate:
 ion_candidate_id:str; theoretical_t1_fragment_id:str; rna_identity:str; sequence_candidate_id:str; cca_state:str; fragment_sequence:str; start_position:int; end_position:int
 ion_mode:T1IonMode; charge:int; theoretical_mz:float; neutral_mass_definition:str; ion_mass_definition:str; proton_count:int
@dataclass(frozen=True,kw_only=True)
class T1FragmentMatch(T1Safeguards):
 t1_fragment_match_id:str; observed_peak_id:str; observed_apex_mz:float; observed_centroid_mz:float|None; observed_t1_mz_mass_definition:str
 theoretical_t1_fragment_id:str; rna_identity:str; sequence_candidate_id:str; cca_state:str; fragment_sequence:str; start_position:int|None; end_position:int|None
 ion_mode:T1IonMode|None; charge:int|None; theoretical_mz:float|None; apex_error_mz:float|None; centroid_error_mz:float|None; absolute_apex_error_mz:float|None; absolute_centroid_error_mz:float|None
 match_class:T1FragmentMatchClass; mass_definition_compatible:bool; polarity_compatible:bool; eligible_for_target_evidence:bool
 observed_peak_candidate_count:int; distinct_fragment_count:int; distinct_charge_count:int; distinct_cca_state_count:int; distinct_rna_identity_count:int
 unique_fragment_assignment:bool; unique_rna_identity_assignment:bool; ambiguity_class:T1FragmentAmbiguityClass
 native_modifications_expected:bool=True; unmodified_fragment_reference_only:bool=True; modified_fragment_composition_not_enumerated:bool=True
@dataclass(frozen=True)
class TheoreticalFragmentDiscrimination:
 fragment_id:str; rna_identity:str; fragment_sequence:str; neutral_monoisotopic_mass:float; discrimination_class:TheoreticalDiscriminationClass
@dataclass(frozen=True,kw_only=True)
class T1IdentityEvidenceSummary(T1Safeguards):
 rna_identity:str; theoretical_specific_fragment_count:int; observed_strict_support_count:int; observed_exploratory_support_count:int; unique_discriminatory_support_count:int; ambiguous_discriminatory_support_count:int; observed_shared_fragment_evidence:int; observed_cross_rna_ambiguous_evidence:int; identity_support_level:IdentitySupportLevel
 native_modifications_expected:bool=True; unmodified_fragment_reference_only:bool=True; modified_fragment_composition_not_enumerated:bool=True

def generate_theoretical_t1_fragments(manifest,rna_identity_id,base_masses,*,candidate_states=None):
 identity=get_rna_identity(manifest,rna_identity_id);variants=generate_cca_tail_variants(identity.sequence,identity.registered_sequence_cca_mode,candidate_states)
 config=RunConfig(digestion={"enabled":True,"enzyme":"RNase_T1","digestion_mode":"specific","missed_cleavages":0,"min_length":1,"include_terminal_forms":False},alkaline_phosphatase={"enabled":False})
 output=[]
 for variant in variants:
  cid=f"{rna_identity_id}__CCA_{variant.candidate_cca_tail_state.value}"
  seq=variant.complete_candidate_sequence;fragments=digest_sequence(cid,seq,{i:i for i in range(1,len(seq)+1)},config,base_masses)
  for index,f in enumerate(fragments,1):
   avg=calculate_intact_rna_average_mass(f.sequence, parameters=IntactRnaMassParameters()).average_neutral_molecular_mass_m
   context=("RNASE_T1_AFTER_G" if f.sequence.endswith("G") else "THREE_PRIME_TERMINAL_FRAGMENT")
   output.append(TheoreticalT1Fragment(f"{cid}__T1_{index:02d}_{f.start}_{f.end}",rna_identity_id,cid,variant.candidate_cca_tail_state.value,"DIGEST_TERMINUS_UNKNOWN","DIGEST_TERMINUS_UNKNOWN",index,f.start,f.end,f.sequence,context,f.unmodified_mass,avg))
 return tuple(output)

def generate_t1_ion_candidates(fragments,*,parameters=None):
 p=parameters or T1FragmentMatchParameters();p.validate();out=[]
 for f in fragments:
  for mode in T1IonMode:
   for z in p.charges:
    if mode is T1IonMode.NEGATIVE_DEPROTONATED:mz=theoretical_mz_from_mass(f.neutral_monoisotopic_mass,z,"negative");protons=-z;definition="MONOISOTOPIC_NEGATIVE_ION_MZ"
    elif mode is T1IonMode.POSITIVE_PROTONATED:mz=theoretical_mz_from_mass(f.neutral_monoisotopic_mass,z,"positive");protons=z;definition="MONOISOTOPIC_POSITIVE_ION_MZ"
    else:mz=f.neutral_monoisotopic_mass/z;protons=0;definition="UNKNOWN_POLARITY_DIAGNOSTIC_MZ"
    out.append(T1IonCandidate(f"{f.theoretical_t1_fragment_id}__{mode.value}__Z{z}",f.theoretical_t1_fragment_id,f.rna_identity,f.sequence_candidate_id,f.cca_state,f.fragment_sequence,f.start_position,f.end_position,mode,z,mz,"MONOISOTOPIC_NEUTRAL",definition,protons))
 return tuple(out)
def _ambiguity(candidates):
 if not candidates:return T1FragmentAmbiguityClass.NO_MATCH
 if len({c.rna_identity for c in candidates})>1:return T1FragmentAmbiguityClass.CROSS_RNA_IDENTITY_AMBIGUITY
 if len({c.theoretical_t1_fragment_id for c in candidates})>1:return T1FragmentAmbiguityClass.MULTIPLE_FRAGMENT_AMBIGUITY
 if len({c.cca_state for c in candidates})>1:return T1FragmentAmbiguityClass.MULTIPLE_CCA_STATE_AMBIGUITY
 if len({c.charge for c in candidates})>1:return T1FragmentAmbiguityClass.MULTIPLE_CHARGE_AMBIGUITY
 if len({c.ion_mass_definition for c in candidates})>1:return T1FragmentAmbiguityClass.MASS_DEFINITION_AMBIGUITY
 return T1FragmentAmbiguityClass.UNIQUE_FRAGMENT_CANDIDATE

def match_observed_t1_fragments(peaks,ion_candidates,*,parameters=None):
 p=parameters or T1FragmentMatchParameters();p.validate();ions=tuple(sorted(ion_candidates,key=lambda c:(c.theoretical_mz,c.ion_candidate_id)));masses=[c.theoretical_mz for c in ions];out=[]
 for peak in sorted(peaks,key=lambda x:(x.apex_mz,x.t1_peak_id)):
  lo=bisect_left(masses,peak.apex_mz-p.exploratory_tolerance_mz);hi=bisect_right(masses,peak.apex_mz+p.exploratory_tolerance_mz);candidates=ions[lo:hi];klass=_ambiguity(candidates)
  counts=(len(candidates),len({c.theoretical_t1_fragment_id for c in candidates}),len({c.charge for c in candidates}),len({c.cca_state for c in candidates}),len({c.rna_identity for c in candidates}))
  if not candidates:
   out.append(T1FragmentMatch(t1_fragment_match_id=f"T1MATCH__{peak.t1_peak_id}__NO_MATCH",observed_peak_id=peak.t1_peak_id,observed_apex_mz=peak.apex_mz,observed_centroid_mz=peak.centroid_mz,observed_t1_mz_mass_definition="UNKNOWN",theoretical_t1_fragment_id="",rna_identity="",sequence_candidate_id="",cca_state="",fragment_sequence="",start_position=None,end_position=None,ion_mode=None,charge=None,theoretical_mz=None,apex_error_mz=None,centroid_error_mz=None,absolute_apex_error_mz=None,absolute_centroid_error_mz=None,match_class=T1FragmentMatchClass.NO_MATCH,mass_definition_compatible=False,polarity_compatible=False,eligible_for_target_evidence=False,observed_peak_candidate_count=0,distinct_fragment_count=0,distinct_charge_count=0,distinct_cca_state_count=0,distinct_rna_identity_count=0,unique_fragment_assignment=False,unique_rna_identity_assignment=False,ambiguity_class=klass));continue
  unique_fragment=counts[1]==1;unique_rna=counts[4]==1
  for c in candidates:
   ae=peak.apex_mz-c.theoretical_mz;ce=peak.centroid_mz-c.theoretical_mz if peak.centroid_mz is not None else None;mc=T1FragmentMatchClass.STRICT if abs(ae)<=p.strict_tolerance_mz else T1FragmentMatchClass.EXPLORATORY
   eligible=mc is T1FragmentMatchClass.STRICT and unique_fragment and unique_rna and not peak.possible_isotope_component
   out.append(T1FragmentMatch(t1_fragment_match_id=f"T1MATCH__{peak.t1_peak_id}__{c.ion_candidate_id}",observed_peak_id=peak.t1_peak_id,observed_apex_mz=peak.apex_mz,observed_centroid_mz=peak.centroid_mz,observed_t1_mz_mass_definition="UNKNOWN",theoretical_t1_fragment_id=c.theoretical_t1_fragment_id,rna_identity=c.rna_identity,sequence_candidate_id=c.sequence_candidate_id,cca_state=c.cca_state,fragment_sequence=c.fragment_sequence,start_position=c.start_position,end_position=c.end_position,ion_mode=c.ion_mode,charge=c.charge,theoretical_mz=c.theoretical_mz,apex_error_mz=ae,centroid_error_mz=ce,absolute_apex_error_mz=abs(ae),absolute_centroid_error_mz=abs(ce) if ce is not None else None,match_class=mc,mass_definition_compatible=False,polarity_compatible=False,eligible_for_target_evidence=eligible,observed_peak_candidate_count=counts[0],distinct_fragment_count=counts[1],distinct_charge_count=counts[2],distinct_cca_state_count=counts[3],distinct_rna_identity_count=counts[4],unique_fragment_assignment=unique_fragment,unique_rna_identity_assignment=unique_rna,ambiguity_class=klass))
 return tuple(out)
def classify_theoretical_fragments(uaa,uag,tolerance=.01):
 out=[]
 for own,other,label in ((uaa,uag,"UAA"),(uag,uaa,"UAG")):
  for f in own:
   identical=any(x.fragment_sequence==f.fragment_sequence for x in other);same_mass=[x for x in other if abs(x.neutral_monoisotopic_mass-f.neutral_monoisotopic_mass)<=tolerance]
   cls=TheoreticalDiscriminationClass.SEQUENCE_IDENTICAL_FRAGMENT if identical else TheoreticalDiscriminationClass.SAME_MASS_DIFFERENT_SEQUENCE if any(x.fragment_sequence!=f.fragment_sequence for x in same_mass) else TheoreticalDiscriminationClass.SHARED_MASS_AMBIGUOUS_FRAGMENT if same_mass else TheoreticalDiscriminationClass.UAA_SPECIFIC_THEORETICAL_FRAGMENT if label=="UAA" else TheoreticalDiscriminationClass.UAG_SPECIFIC_THEORETICAL_FRAGMENT
   out.append(TheoreticalFragmentDiscrimination(f.theoretical_t1_fragment_id,f.rna_identity,f.fragment_sequence,f.neutral_monoisotopic_mass,cls))
 return tuple(out)
def summarize_identity_evidence(rna_identity,discrimination,matches):
 specific_class=TheoreticalDiscriminationClass.UAA_SPECIFIC_THEORETICAL_FRAGMENT if "UAA" in rna_identity else TheoreticalDiscriminationClass.UAG_SPECIFIC_THEORETICAL_FRAGMENT
 ids={x.fragment_id for x in discrimination if x.rna_identity==rna_identity and x.discrimination_class is specific_class};rows=[m for m in matches if m.rna_identity==rna_identity and m.theoretical_t1_fragment_id in ids and m.match_class is not T1FragmentMatchClass.NO_MATCH]
 strict={m.observed_peak_id for m in rows if m.match_class is T1FragmentMatchClass.STRICT};expl={m.observed_peak_id for m in rows if m.match_class is T1FragmentMatchClass.EXPLORATORY};unique={m.observed_peak_id for m in rows if m.eligible_for_target_evidence};amb={m.observed_peak_id for m in rows if not m.eligible_for_target_evidence}
 n=len(unique);level=IdentitySupportLevel.STRONG_DISCRIMINATORY_SUPPORT if n>=3 else IdentitySupportLevel.MODERATE_DISCRIMINATORY_SUPPORT if n==2 else IdentitySupportLevel.WEAK_DISCRIMINATORY_SUPPORT if n==1 or strict or expl else IdentitySupportLevel.NO_DISCRIMINATORY_SUPPORT
 return T1IdentityEvidenceSummary(rna_identity=rna_identity,theoretical_specific_fragment_count=len(ids),observed_strict_support_count=len(strict),observed_exploratory_support_count=len(expl),unique_discriminatory_support_count=len(unique),ambiguous_discriminatory_support_count=len(amb),observed_shared_fragment_evidence=len({m.observed_peak_id for m in matches if m.rna_identity==rna_identity and m.theoretical_t1_fragment_id not in ids and m.match_class is not T1FragmentMatchClass.NO_MATCH}),observed_cross_rna_ambiguous_evidence=len({m.observed_peak_id for m in matches if m.ambiguity_class is T1FragmentAmbiguityClass.CROSS_RNA_IDENTITY_AMBIGUITY}),identity_support_level=level)
