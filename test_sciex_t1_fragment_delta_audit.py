from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import pytest
from rna_masshunter.modifications import load_modifications
from rna_masshunter.sciex_intact_peak_family import DeltaMassDefinition
from rna_masshunter.sciex_t1_fragment_shadow_match import T1IonCandidate,T1IonMode
from rna_masshunter.sciex_t1_profile_peak_audit import T1ProfilePeak,T1PeakQualityClass
from rna_masshunter.sciex_t1_fragment_delta_audit import *
ROOT=Path(__file__).parent

def peak(pid,mz,*,iso=False,shoulder=False,duplicate=False):
 return T1ProfilePeak(t1_peak_id=pid,source_id='S',measurement_id='M',rna_identity='R',apex_mz=mz,centroid_mz=mz+.002,apex_intensity=100,integrated_intensity=10,relative_apex_intensity=1,relative_integrated_intensity=1,left_boundary_mz=mz-.03,right_boundary_mz=mz+.03,peak_width_mz=.06,fwhm_mz=.02,prominence=90,relative_prominence=.9,sharpness_score=4500,nearest_peak_separation_mz=None,shared_valley_fraction=0,possible_shoulder=shoulder,possible_duplicate=duplicate,possible_isotope_component=iso,possible_isotope_charge=1 if iso else None,possible_isotope_spacing_error=0 if iso else None,possible_overlapping_envelope=iso,peak_quality_class=T1PeakQualityClass.MAJOR_SHARP,selected_as_primary=True)
def ion(fid='F',mz=100,charge=1,mode=T1IonMode.NEGATIVE_DEPROTONATED,rna='R',seq='AG',start=1,end=2,cca='CCA'):
 return T1IonCandidate(f'{fid}_{cca}_{mode.value}_{charge}',fid,rna,rna+'_'+cca,cca,seq,start,end,mode,charge,mz,'MONOISOTOPIC_NEUTRAL','MONOISOTOPIC_ION_MZ',-charge)
def core_refs():return build_chemical_delta_reference_registry()
def audit(peaks,ions,**kwargs):return audit_t1_fragment_chemical_deltas(peaks,ions,source_id='S',references=core_refs(),parameters=ChemicalDeltaAuditParameters(**kwargs))
def relation(delta,charge=1,centroid=True):
 p=peak('P',100+delta/charge)
 if not centroid:p=replace(p,centroid_mz=None)
 return audit([p],[ion(mz=100,charge=charge)]).relations[0]

def refs():return {r.reference_id:r for r in core_refs()}
@pytest.mark.parametrize('rid,value,definition',[('O_ADDITION_AVERAGE',15.9994,DeltaMassDefinition.AVERAGE_DELTA),('O_ADDITION_MONOISOTOPIC',15.99491461957,DeltaMassDefinition.MONOISOTOPIC_DELTA),('H2O_ADDITION_AVERAGE',18.01528,DeltaMassDefinition.AVERAGE_DELTA),('H2O_ADDITION_MONOISOTOPIC',18.01056468403,DeltaMassDefinition.MONOISOTOPIC_DELTA),('S_ADDITION_AVERAGE',32.065,DeltaMassDefinition.AVERAGE_DELTA),('S_ADDITION_MONOISOTOPIC',31.9720711744,DeltaMassDefinition.MONOISOTOPIC_DELTA),('O_TO_S_AVERAGE',16.0656,DeltaMassDefinition.AVERAGE_DELTA),('O_TO_S_MONOISOTOPIC',15.97715655483,DeltaMassDefinition.MONOISOTOPIC_DELTA),('S_TO_O_AVERAGE',-16.0656,DeltaMassDefinition.AVERAGE_DELTA)])
def test_reference_values_mass_definition_and_provenance(rid,value,definition):
 r=refs()[rid];assert r.signed_delta_da==pytest.approx(value);assert r.reference_mass_definition is definition;assert r.reference_provenance

def test_modifications_yaml_classified_mono():
 registry=build_chemical_delta_reference_registry(load_modifications(ROOT/'data/modifications.yaml'));mods=[r for r in registry if r.reference_category is ChemicalReferenceCategory.KNOWN_MODIFICATION_DELTA];assert mods and all(r.reference_mass_definition is DeltaMassDefinition.MONOISOTOPIC_DELTA and 'Mass(mono)' in r.reference_provenance for r in mods)
@pytest.mark.parametrize('charge,delta,expected',[(1,.1,.1),(2,.1,.2),(3,.1,.3),(1,-.1,-.1),(2,-.1,-.2)])
def test_signed_neutral_delta_conversion(charge,delta,expected):
 mz,neutral,compatible=convert_mz_delta_to_neutral_delta(100+delta,100,charge);assert mz==pytest.approx(delta);assert neutral==pytest.approx(expected);assert compatible

def test_ion_convention_mismatch():assert convert_mz_delta_to_neutral_delta(101,100,1,observed_ion_mode=T1IonMode.POSITIVE_PROTONATED,theoretical_ion_mode=T1IonMode.NEGATIVE_DEPROTONATED)[2] is False
def test_zero_charge_rejected():
 with pytest.raises(ValueError):convert_mz_delta_to_neutral_delta(101,100,0)
def test_conversion_deterministic():assert convert_mz_delta_to_neutral_delta(101,100,2)==convert_mz_delta_to_neutral_delta(101,100,2)
@pytest.mark.parametrize('delta,hypothesis',[ (15.9994,ChemicalHypothesisClass.O_EQUIVALENT),(18.01528,ChemicalHypothesisClass.H2O_EQUIVALENT),(-18.01528,ChemicalHypothesisClass.H2O_LOSS_EQUIVALENT),(32.065,ChemicalHypothesisClass.S_EQUIVALENT),(16.0656,ChemicalHypothesisClass.O_TO_S_EQUIVALENT),(-16.0656,ChemicalHypothesisClass.S_TO_O_EQUIVALENT)])
def test_chemical_strict_matches(delta,hypothesis):
 r=relation(delta);assert r.delta_relation_class is DeltaRelationClass.CHEMICAL_DELTA_STRICT;assert r.chemical_hypothesis_class is hypothesis;assert r.apex_neutral_delta==pytest.approx(delta)
def test_o_exploratory():
 ref=refs()['O_ADDITION_AVERAGE'];r=audit_t1_fragment_chemical_deltas([peak('P',116.08)],[ion()],source_id='S',references=(ref,)).relations[0];assert r.reference_match_class=='EXPLORATORY' and r.chemical_hypothesis_class is ChemicalHypothesisClass.O_EQUIVALENT
def test_known_modification_mismatch_diagnostic():
 mod=SimpleNamespace(id='MODX',symbol=None,mass_shift_from_unmodified=10.0);r=audit_t1_fragment_chemical_deltas([peak('P',110)],[ion()],source_id='S',modifications=[mod]).relations[0];assert r.delta_relation_class is DeltaRelationClass.MASS_DEFINITION_MISMATCH_DIAGNOSTIC;assert not r.eligible_for_chemical_delta_evidence

def test_multiple_reference_mass_definition_ambiguity():
 r=relation(15.9994);assert r.candidate_reference_count>=2;assert r.reference_ambiguity_class is ReferenceAmbiguityClass.MULTIPLE_REFERENCE_CATEGORIES
def test_no_reference_match():assert relation(123.456).delta_relation_class is DeltaRelationClass.NO_REFERENCE_MATCH
def test_apex_primary_centroid_secondary():
 r=relation(15.9994);assert r.apex_delta_error_da==pytest.approx(0,abs=.01);assert r.centroid_delta_error_da!=r.apex_delta_error_da
def test_charge_scaled_tolerance():
 ref=refs()['O_ADDITION_AVERAGE'];r=audit_t1_fragment_chemical_deltas([peak('P',101.608)],[ion(mz=100,charge=10)],source_id='S',references=(ref,)).relations[0];assert r.reference_match_class=='STRICT'
def test_candidate_and_report_bounds():
 ions=[ion(str(i),100+i*.01,charge=1) for i in range(100)];r=audit([peak('P',116)],ions,maximum_theoretical_candidates_per_observed_peak_per_charge=5,maximum_reported_delta_relations_per_peak=3);assert r.candidate_explosion.post_per_peak_cap_candidate_count<=5;assert len(r.relations)<=3

def test_same_fragment_o_edge_and_no_reverse_duplicate():
 r=audit([peak('A',100),peak('B',115.9994)],[ion()]);assert r.state_edges;assert all(e.lower_observed_peak_id=='A' and e.higher_observed_peak_id=='B' for e in r.state_edges)
@pytest.mark.parametrize('delta,category',[(18.01528,ChemicalReferenceCategory.WATER_ADDITION_EQUIVALENT),(16.0656,ChemicalReferenceCategory.O_TO_S_SUBSTITUTION_EQUIVALENT)])
def test_same_fragment_other_edges(delta,category):assert any(e.reference_category is category for e in audit([peak('A',100),peak('B',100+delta)],[ion()]).state_edges)
def test_different_fragment_never_combined_in_edge_identity():
 r=audit([peak('A',100),peak('B',116)],[ion('F1'),ion('F2')]);assert all(e.shared_theoretical_fragment_id in {'R:1-2:AG'} for e in r.state_edges)
def test_different_charge_and_mode_are_separate_groups():
 r=audit([peak('A',100),peak('B',116)],[ion(charge=1),ion(charge=2,mode=T1IonMode.POSITIVE_PROTONATED)]);assert all(e.shared_charge in {1,2} and isinstance(e.shared_ion_mode,T1IonMode) for e in r.state_edges)
@pytest.mark.parametrize('flag',['iso','shoulder','duplicate'])
def test_artifact_peak_excluded_from_edges(flag):
 kwargs={flag:True};r=audit([peak('A',100),peak('B',116,**kwargs)],[ion()]);assert not r.state_edges

def test_connected_mixed_o_h2o_series():
 r=audit([peak('A',100),peak('B',115.9994),peak('C',134.01468)],[ion()]);s=r.state_series[0];assert s.member_count==3 and s.series_pattern is T1StateSeriesPattern.MIXED_O_H2O_SERIES

def test_connected_mixed_o_s_series_and_deterministic_id():
 args=([peak('A',100),peak('B',115.9994),peak('C',132.065)],[ion()]);a=audit(*args);b=audit(*args);assert a==b;assert any(s.series_pattern is T1StateSeriesPattern.MIXED_O_S_SERIES for s in a.state_series);assert a.state_series[0].t1_state_series_id==b.state_series[0].t1_state_series_id

def test_input_nonmutation():
 peaks=[peak('A',100),peak('B',116)];ions=[ion()];before=(repr(peaks),repr(ions));audit(peaks,ions);assert before==(repr(peaks),repr(ions))
@pytest.mark.parametrize('flag',['thioamide_assigned','oxidation_assigned','hydration_assigned','modification_position_assigned','structure_assigned','charge_state_confirmed','polarity_confirmed','applied_to_formal_score','applied_to_ranking','applied_to_candidate_filtering','applied_to_final_consensus'])
def test_relation_safeguards_false(flag):assert getattr(relation(16.0656),flag) is False

def test_all_references_are_diagnostic_only():assert all(r.diagnostic_only for r in core_refs())
def test_signed_o_to_s_and_s_to_o_are_opposites():assert refs()['O_TO_S_AVERAGE'].signed_delta_da==-refs()['S_TO_O_AVERAGE'].signed_delta_da
def test_relation_ion_charge_compatibility_fields():
 r=relation(15.9994);assert r.ion_convention_compatible and r.charge_convention_compatible and r.eligible_for_neutral_delta_audit;assert r.observed_theoretical_mass_definition_compatibility=='UNKNOWN'
def test_cross_sample_category_comparison_both_and_only():
 u=audit([peak('U',116)],[ion()]);g=audit([peak('G',116)],[ion()]);rows=compare_cross_sample_categories(u,g);o=next(x for x in rows if x.reference_category is ChemicalReferenceCategory.OXYGEN_ADDITION_EQUIVALENT);assert o.detected_in_both and not o.detected_only_in_uaa

def test_cross_sample_fragment_delta_comparison():
 u=audit([peak('U',116)],[ion()]);g=audit([peak('G',116)],[ion()]);rows=compare_cross_sample_fragment_deltas(u,g);assert any(x.shared_delta_hypothesis for x in rows)
def test_series_and_edge_false_certainty_flags():
 r=audit([peak('A',100),peak('B',116)],[ion()]);records=(*r.state_edges,*r.state_series);assert records
 for x in records:assert x.thioamide_assigned is False and x.oxidation_assigned is False and x.applied_to_formal_score is False


def test_known_modification_comparison_role_is_retained():
 mod=SimpleNamespace(id='MODX',symbol=None,mass_shift_from_unmodified=10.0);r=audit_t1_fragment_chemical_deltas([peak('P',110)],[ion()],source_id='S',modifications=[mod]).relations[0];assert r.comparison_role=='MASS_DEFINITION_MISMATCH_DIAGNOSTIC'

def test_summary_false_certainty_flags():
 r=audit([peak('P',116)],[ion()]);records=(r,r.candidate_explosion,*compare_cross_sample_categories(r,r),*compare_cross_sample_fragment_deltas(r,r))
 for record in records:
  assert record.shadow_analysis_only is True and record.mass_evidence_only is True
  assert record.thioamide_assigned is False and record.modification_position_assigned is False and record.applied_to_formal_score is False

def test_cross_rna_sequence_identity_ambiguity_is_retained():
 u=audit([peak('U',116)],[ion(rna='TRNA_LEU_UAA')]);g=audit([peak('G',116)],[ion(rna='TRNA_LEU_UAG')]);combined=compare_cross_sample_fragment_deltas(u,g);assert combined and all(x.cross_rna_sequence_identity_ambiguous for x in combined)
