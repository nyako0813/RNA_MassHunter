from pathlib import Path
from unittest.mock import patch
import pytest
from rna_masshunter.masses import load_base_masses,PROTON_MASS
from rna_masshunter.sciex_sample_manifest import load_sciex_sample_manifest
from rna_masshunter.sciex_t1_profile_peak_audit import T1ProfilePeak,T1PeakQualityClass
from rna_masshunter.sciex_t1_fragment_shadow_match import *
ROOT=Path(__file__).parent
@pytest.fixture(scope='module')
def manifest():return load_sciex_sample_manifest(ROOT/'data/sciex_sample_manifest.yaml')
@pytest.fixture(scope='module')
def base():return load_base_masses(ROOT/'data/base_masses.yaml')
def peak(mz=100,iso=False):
 return T1ProfilePeak(t1_peak_id='P',source_id='S',measurement_id='M',rna_identity='R',apex_mz=mz,centroid_mz=mz+.002,apex_intensity=100,integrated_intensity=10,relative_apex_intensity=1,relative_integrated_intensity=1,left_boundary_mz=mz-.03,right_boundary_mz=mz+.03,peak_width_mz=.06,fwhm_mz=.02,prominence=90,relative_prominence=.9,sharpness_score=4500,nearest_peak_separation_mz=None,shared_valley_fraction=0,possible_shoulder=False,possible_duplicate=False,possible_isotope_component=iso,possible_isotope_charge=1 if iso else None,possible_isotope_spacing_error=0 if iso else None,possible_overlapping_envelope=iso,peak_quality_class=T1PeakQualityClass.MAJOR_SHARP,selected_as_primary=True)
def fragment(fid='F',rna='TRNA_LEU_UAA',seq='AG',cca='CCA',mass=1000):return TheoreticalT1Fragment(fid,rna,rna+'__'+cca,cca,'DIGEST_TERMINUS_UNKNOWN','DIGEST_TERMINUS_UNKNOWN',1,1,len(seq),seq,'RNASE_T1_AFTER_G',mass,mass+1)
def ion(fid='F',rna='TRNA_LEU_UAA',mz=100,charge=1,cca='CCA',seq='AG',mode=T1IonMode.NEGATIVE_DEPROTONATED):return T1IonCandidate(fid+mode.value+str(charge)+rna,fid,rna,rna+'__'+cca,cca,seq,1,len(seq),mode,charge,mz,'MONOISOTOPIC_NEUTRAL','MONOISOTOPIC_ION_MZ',-charge)
def test_existing_digest_logic_is_reused(manifest,base):
 with patch('rna_masshunter.sciex_t1_fragment_shadow_match.digest_sequence',wraps=digest_sequence) as wrapped:
  out=generate_theoretical_t1_fragments(manifest,'TRNA_LEU_UAA',base,candidate_states=['CCA']);assert out and wrapped.called
@pytest.mark.parametrize('rna', ['TRNA_LEU_UAA','TRNA_LEU_UAG'])
def test_registered_theoretical_fragments_and_fields(manifest,base,rna):
 out=generate_theoretical_t1_fragments(manifest,rna,base);assert out;assert len({f.sequence_candidate_id for f in out})==4;assert {f.cca_state for f in out}=={'NONE','C','CC','CCA'}
 assert all(f.neutral_monoisotopic_mass>0 and f.neutral_average_mass>0 and f.start_position<=f.end_position for f in out)
@pytest.mark.parametrize('charge',[1,2,3,4,5])
def test_charge_generation(charge):
 ions=generate_t1_ion_candidates([fragment()],parameters=T1FragmentMatchParameters(charges=(charge,)));assert len(ions)==3 and {i.charge for i in ions}=={charge};assert {i.ion_mode for i in ions}==set(T1IonMode)
def test_positive_negative_unknown_modes_and_project_proton_mass():
 ions={i.ion_mode:i for i in generate_t1_ion_candidates([fragment(mass=1000)],parameters=T1FragmentMatchParameters(charges=(2,)))}
 assert ions[T1IonMode.NEGATIVE_DEPROTONATED].theoretical_mz==pytest.approx((1000-2*PROTON_MASS)/2)
 assert ions[T1IonMode.POSITIVE_PROTONATED].theoretical_mz==pytest.approx((1000+2*PROTON_MASS)/2)
 assert ions[T1IonMode.UNKNOWN_POLARITY_DIAGNOSTIC].theoretical_mz==pytest.approx(500)
@pytest.mark.parametrize('offset,klass',[(.005,T1FragmentMatchClass.STRICT),(.015,T1FragmentMatchClass.EXPLORATORY),(.03,T1FragmentMatchClass.NO_MATCH)])
def test_match_classes_apex_primary_centroid_secondary(offset,klass):
 rows=match_observed_t1_fragments([peak(100+offset)],[ion(mz=100)]);assert rows[0].match_class is klass
 if klass is not T1FragmentMatchClass.NO_MATCH:assert rows[0].apex_error_mz==pytest.approx(offset);assert rows[0].centroid_error_mz==pytest.approx(offset+.002)
def test_multiple_charge_ambiguity():
 rows=match_observed_t1_fragments([peak()],[ion(charge=1),ion(charge=2)]);assert all(r.ambiguity_class is T1FragmentAmbiguityClass.MULTIPLE_CHARGE_AMBIGUITY for r in rows);assert rows[0].distinct_charge_count==2
def test_multiple_fragment_and_cca_ambiguity():
 rows=match_observed_t1_fragments([peak()],[ion(fid='F1'),ion(fid='F2',cca='CC')]);assert all(r.ambiguity_class is T1FragmentAmbiguityClass.MULTIPLE_FRAGMENT_AMBIGUITY for r in rows);assert rows[0].distinct_fragment_count==2 and rows[0].distinct_cca_state_count==2
def test_cross_rna_ambiguity():
 rows=match_observed_t1_fragments([peak()],[ion(rna='TRNA_LEU_UAA'),ion(rna='TRNA_LEU_UAG')]);assert all(r.ambiguity_class is T1FragmentAmbiguityClass.CROSS_RNA_IDENTITY_AMBIGUITY for r in rows);assert rows[0].distinct_rna_identity_count==2
def test_unique_candidate_is_computational_only_not_assignment():
 row=match_observed_t1_fragments([peak()],[ion()])[0];assert row.ambiguity_class is T1FragmentAmbiguityClass.UNIQUE_FRAGMENT_CANDIDATE;assert row.unique_fragment_assignment is True;assert row.fragment_assigned is False and row.unique_fragment_assigned is False;assert row.observed_t1_mz_mass_definition=='UNKNOWN';assert row.mass_definition_compatible is False and row.polarity_compatible is False
@pytest.mark.parametrize('iso,eligible',[(False,True),(True,False)])
def test_isotope_excluded_from_unique_target_evidence(iso,eligible):assert match_observed_t1_fragments([peak(iso=iso)],[ion()])[0].eligible_for_target_evidence is eligible
def test_discrimination_classes_same_mass_and_specific():
 u=[fragment('U1',seq='AG',mass=100),fragment('U2',seq='AA',mass=200)];g=[fragment('G1','TRNA_LEU_UAG',seq='AG',mass=100),fragment('G2','TRNA_LEU_UAG',seq='CC',mass=200)]
 classes={x.discrimination_class for x in classify_theoretical_fragments(u,g)};assert TheoreticalDiscriminationClass.SEQUENCE_IDENTICAL_FRAGMENT in classes;assert TheoreticalDiscriminationClass.SAME_MASS_DIFFERENT_SEQUENCE in classes
def test_identity_support_and_safeguards():
 d=(TheoreticalFragmentDiscrimination('F','TRNA_LEU_UAA','AG',100,TheoreticalDiscriminationClass.UAA_SPECIFIC_THEORETICAL_FRAGMENT),);m=match_observed_t1_fragments([peak()],[ion()]);s=summarize_identity_evidence('TRNA_LEU_UAA',d,m);assert s.identity_support_level is IdentitySupportLevel.WEAK_DISCRIMINATORY_SUPPORT
 assert s.target_rna_identity_confirmed_by_mass is False and s.native_modifications_expected is True and s.modified_fragment_composition_not_enumerated is True
 for n in ('modification_assigned','position_assigned','structure_assigned','applied_to_formal_score','applied_to_ranking','applied_to_candidate_filtering','applied_to_final_consensus'):assert getattr(s,n) is False

@pytest.mark.parametrize('flag',['charge_state_confirmed','ion_mode_confirmed','polarity_confirmed'])
def test_fragment_relation_ion_certainty_flags_false(flag):
 assert getattr(match_observed_t1_fragments([peak()],[ion()])[0],flag) is False

@pytest.mark.parametrize('flag',['co_captured_rna_excluded','background_component_excluded','contaminant_excluded'])
def test_fragment_relation_exclusion_flags_false(flag):
 assert getattr(match_observed_t1_fragments([peak()],[ion()])[0],flag) is False
