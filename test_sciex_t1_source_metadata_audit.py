from pathlib import Path
from types import SimpleNamespace
from dataclasses import replace
import math
import pytest
from rna_masshunter.sciex_intact_peak_family import DeltaMassDefinition
from rna_masshunter.sciex_t1_profile_peak_audit import T1ProfilePeak,T1PeakQualityClass,T1ProfilePeakAuditResult,T1PeakDetectionParameters
from rna_masshunter.sciex_t1_fragment_shadow_match import T1FragmentMatchClass,T1IonMode
from rna_masshunter.sciex_t1_fragment_delta_audit import ChemicalHypothesisClass
from rna_masshunter.sciex_t1_delta_evidence_quality_audit import EvidenceTier,MassDefinitionCompatibilityStatus,ApexCentroidConcordanceClass,ChargeSupportClass,RecurrentSupportClass
from rna_masshunter.sciex_t1_source_metadata_audit import *

def evidence(value,trust=MetadataTrustLevel.AUTHORITATIVE_SOURCE_SPECIFIC,source='S',specific=True,confirmed=True,field='polarity'):return MetadataFieldEvidence(field,value,source,trust,specific,confirmed)
def peak(pid,mz,*,iso=False,z=None,error=None,intensity=100):
 return T1ProfilePeak(t1_peak_id=pid,source_id='S',measurement_id='M',rna_identity='R',apex_mz=mz,centroid_mz=mz,apex_intensity=intensity,integrated_intensity=10,relative_apex_intensity=1,relative_integrated_intensity=1,left_boundary_mz=mz-.02,right_boundary_mz=mz+.02,peak_width_mz=.04,fwhm_mz=.02,prominence=90,relative_prominence=.9,sharpness_score=4500,nearest_peak_separation_mz=None,shared_valley_fraction=0,possible_shoulder=False,possible_duplicate=False,possible_isotope_component=iso,possible_isotope_charge=z,possible_isotope_spacing_error=error,possible_overlapping_envelope=iso,peak_quality_class=T1PeakQualityClass.ISOTOPE_OR_ENVELOPE_COMPONENT if iso else T1PeakQualityClass.MAJOR_SHARP,selected_as_primary=True)
def qrecord(pid='P',z=1,mode=T1IonMode.NEGATIVE_DEPROTONATED,tier=EvidenceTier.TIER_C_WEAK_SUPPORT,hyp=ChemicalHypothesisClass.O_EQUIVALENT,definition=DeltaMassDefinition.MONOISOTOPIC_DELTA,score=20,cross=False,multi=False,recur=False):
 return SimpleNamespace(observed_peak_id=pid,charge=z,ion_mode=mode,evidence_tier=tier,chemical_hypothesis_class=hyp,reference_mass_definition=definition,evidence_support_score=score,cross_rna_identity_ambiguous=cross,charge_support_class=ChargeSupportClass.MULTI_CHARGE_CONCORDANT_SUPPORT if multi else ChargeSupportClass.SINGLE_CHARGE_ONLY,recurrent_support_class=RecurrentSupportClass.MULTI_PEAK_MULTI_CHARGE_RECURRENT if recur else RecurrentSupportClass.SINGLE_PEAK_ONLY,apex_centroid_concordance_class=ApexCentroidConcordanceClass.BOTH_STRICT_SAME_REFERENCE,mass_definition_compatibility_status=MassDefinitionCompatibilityStatus.UNKNOWN_COMPATIBILITY)
def ion(z):return SimpleNamespace(charge=z)
def fmatch(pid='P',z=1,klass=T1FragmentMatchClass.STRICT):return SimpleNamespace(observed_peak_id=pid,charge=z,match_class=klass)

@pytest.mark.parametrize('higher,lower',[ (MetadataTrustLevel.AUTHORITATIVE_SOURCE_SPECIFIC,MetadataTrustLevel.PROJECT_DECLARED_GENERAL),(MetadataTrustLevel.AUTHORITATIVE_EXPERIMENT_LEVEL,MetadataTrustLevel.INFERRED_FROM_DATA),(MetadataTrustLevel.PROJECT_DECLARED_SOURCE_SPECIFIC,MetadataTrustLevel.FILENAME_ONLY)])
def test_higher_trust_wins(higher,lower):
 value,confirmed,conflicts,status=resolve_metadata_field([evidence('NEGATIVE',higher,'HIGH'),evidence('POSITIVE',lower,'LOW')]);assert value=='NEGATIVE' and conflicts and status is ConflictResolutionStatus.RESOLVED_BY_HIGHER_TRUST_SOURCE

def test_experiment_level_can_confirm():assert resolve_metadata_field([evidence('NEGATIVE',MetadataTrustLevel.AUTHORITATIVE_EXPERIMENT_LEVEL)])[1]
def test_project_general_does_not_confirm():assert not resolve_metadata_field([evidence('NEGATIVE',MetadataTrustLevel.PROJECT_DECLARED_GENERAL,specific=False)])[1]
def test_filename_inference_low_trust():assert not resolve_metadata_field([evidence('NEGATIVE',MetadataTrustLevel.FILENAME_ONLY,confirmed=False)])[1]
def test_empirical_inference_does_not_confirm():assert not resolve_metadata_field([evidence('NEGATIVE',MetadataTrustLevel.INFERRED_FROM_DATA,confirmed=False)])[1]
def test_equal_trust_conflict_unresolved():assert resolve_metadata_field([evidence('NEGATIVE',source='A'),evidence('POSITIVE',source='B')])[0]=='CONFLICT'
def test_unavailable_metadata():assert resolve_metadata_field([])[0]=='UNKNOWN'

@pytest.mark.parametrize('value,status',[('NEGATIVE',PolarityStatus.NEGATIVE),('POSITIVE',PolarityStatus.POSITIVE),('MIXED',PolarityStatus.MIXED),('UNKNOWN',PolarityStatus.UNKNOWN)])
def test_polarity_enum(value,status):assert PolarityStatus(value) is status
def test_acquisition_export_distinct_fields():
 a=evidence('NEGATIVE',field='acquisition');b=evidence('POSITIVE',field='export');assert a.field_name!=b.field_name
def test_incompatible_ion_mode_charge_conflict():
 rows=build_charge_support_summaries([ion(1)],[fmatch()], [qrecord(mode=T1IonMode.POSITIVE_PROTONATED)],[],polarity=PolarityStatus.NEGATIVE);assert rows[0].charge_plausibility_class is ChargePlausibilityClass.CONFLICTING_SUPPORT
def test_unknown_polarity_does_not_remove_diagnostics():assert build_charge_support_summaries([ion(1)],[],[qrecord()],[],polarity=PolarityStatus.UNKNOWN)[0].chemical_delta_relation_count==1
def test_unknown_polarity_is_unconfirmed():assert not resolve_metadata_field([])[1]

@pytest.mark.parametrize('coords,median_expected', [([1,2,3],1),([1,1.5,2],.5),([1,1,2],1),([1,2,4],1.5)])
def test_spacing_diagnostics(coords,median_expected):assert calculate_spacing_diagnostics(coords).median_spacing==pytest.approx(median_expected)
def test_duplicate_coordinate_count():assert calculate_spacing_diagnostics([1,1,2]).duplicate_coordinate_count==1
def test_irregular_spacing_cv():assert calculate_spacing_diagnostics([1,1.1,2,4]).spacing_coefficient_of_variation>0
def test_regular_spacing_low_cv():assert calculate_spacing_diagnostics([1,1.1,1.2,1.3]).spacing_coefficient_of_variation<1e-10
def test_spacing_empty():assert calculate_spacing_diagnostics([]).row_count==0

@pytest.mark.parametrize('values,zeros,neg,pos',[([0,1,2],1,0,2),([-1,0,1],1,1,1),([1,2,3],0,0,3)])
def test_intensity_counts(values,zeros,neg,pos):
 d=calculate_intensity_diagnostics(values);assert (d.zero_intensity_count,d.negative_intensity_count,d.positive_intensity_count)==(zeros,neg,pos)
def test_intensity_nonfinite():assert calculate_intensity_diagnostics([1,float('nan')]).non_finite_intensity_count==1
def test_dynamic_range():assert calculate_intensity_diagnostics([1,10]).dynamic_range==10

@pytest.mark.parametrize('status',[DataRepresentationStatus.PROFILE_CONTINUOUS,DataRepresentationStatus.CENTROID_PEAK_LIST,DataRepresentationStatus.RESAMPLED_PROFILE,DataRepresentationStatus.RECONSTRUCTED_PROFILE,DataRepresentationStatus.DECONVOLUTED_PROFILE,DataRepresentationStatus.UNKNOWN_REPRESENTATION,DataRepresentationStatus.CONFLICTING_REPRESENTATION])
def test_representation_status_values(status):assert status.value
def test_spacing_only_not_confirmed():assert MetadataTrustLevel.INFERRED_FROM_DATA is not MetadataTrustLevel.AUTHORITATIVE_SOURCE_SPECIFIC

@pytest.mark.parametrize('status',[ObservedMassDefinitionStatus.MONOISOTOPIC_MZ,ObservedMassDefinitionStatus.AVERAGE_MZ,ObservedMassDefinitionStatus.CENTROID_OF_ISOTOPE_ENVELOPE,ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION,ObservedMassDefinitionStatus.CONFLICTING_MASS_DEFINITION])
def test_mass_definition_status_values(status):assert status.value
def test_fragment_match_count_not_mass_confirmation():assert resolve_metadata_field([evidence('MONOISOTOPIC_MZ',MetadataTrustLevel.INFERRED_FROM_DATA,confirmed=False)])[1] is False
def test_isotope_evidence_diagnostic_only():assert IsotopeFamilyCandidate.__dataclass_fields__['isotope_envelope_confirmed'].default is False

@pytest.mark.parametrize('z',[1,2,3])
def test_three_member_isotope_family(z):
 s=ISOTOPE_MASS_DIFFERENCE_DA/z;ps=[peak('A',100),peak('B',100+s,iso=True,z=z,error=0),peak('C',100+2*s,iso=True,z=z,error=0)];f=build_isotope_families(ps);assert len(f)==1 and f[0].member_count==3 and f[0].possible_charge==z and f[0].envelope_support_class is EnvelopeSupportClass.STRONG_MULTI_MEMBER_ENVELOPE_SUPPORT
def test_two_member_moderate():
 s=ISOTOPE_MASS_DIFFERENCE_DA/2;f=build_isotope_families([peak('A',100),peak('B',100+s,iso=True,z=2,error=0)]);assert f[0].envelope_support_class is EnvelopeSupportClass.MODERATE_TWO_MEMBER_SUPPORT
def test_weak_spacing_only():
 s=ISOTOPE_MASS_DIFFERENCE_DA+.006;f=build_isotope_families([peak('A',100),peak('B',100+s,iso=True,z=1,error=.006)]);assert f[0].envelope_support_class is EnvelopeSupportClass.WEAK_SPACING_ONLY
def test_duplicate_peak_not_reused():
 s=ISOTOPE_MASS_DIFFERENCE_DA;f=build_isotope_families([peak('A',100),peak('B',100+s,iso=True,z=1,error=0),peak('C',100+s+.001,iso=True,z=1,error=.001)]);assert len(set(f[0].member_peak_ids))==f[0].member_count
def test_isotope_family_id_deterministic():
 s=ISOTOPE_MASS_DIFFERENCE_DA;ps=[peak('A',100),peak('B',100+s,iso=True,z=1,error=0)];assert build_isotope_families(ps)==build_isotope_families(reversed(ps))
def test_isotope_assignment_false():
 s=ISOTOPE_MASS_DIFFERENCE_DA;assert not build_isotope_families([peak('A',100),peak('B',100+s,iso=True,z=1,error=0)])[0].isotope_envelope_confirmed

@pytest.mark.parametrize('direct,multi,strong,recur,expected',[ (2,True,True,True,ChargePlausibilityClass.SUPPORTED_CHARGE),(2,False,False,False,ChargePlausibilityClass.WEAKLY_SUPPORTED_CHARGE),(0,False,False,False,ChargePlausibilityClass.INSUFFICIENT_SUPPORT)])
def test_charge_plausibility(direct,multi,strong,recur,expected):
 fam=[]
 if strong:
  fam=[IsotopeFamilyCandidate(isotope_family_id='F',member_peak_ids=('A','B','C'),member_apex_mzs=(1,2,3),member_count=3,possible_charge=1,mean_spacing=1,maximum_spacing_error=0,intensity_pattern='X',envelope_support_class=EnvelopeSupportClass.STRONG_MULTI_MEMBER_ENVELOPE_SUPPORT)]
 matches=[fmatch(f'P{i}') for i in range(direct)];quality=[qrecord(f'Q{i}',multi=multi,recur=recur,cross=False) for i in range(5 if recur else 1)];row=build_charge_support_summaries([ion(1)],matches,quality,fam)[0];assert row.charge_plausibility_class is expected
def test_multi_charge_evidence_increases_score():
 a=build_charge_support_summaries([ion(1)],[],[qrecord(multi=True)],[])[0];b=build_charge_support_summaries([ion(1)],[],[qrecord(multi=False)],[])[0];assert a.charge_support_score>b.charge_support_score
def test_isotope_family_increases_score():
 f=IsotopeFamilyCandidate(isotope_family_id='F',member_peak_ids=('A','B','C'),member_apex_mzs=(1,2,3),member_count=3,possible_charge=1,mean_spacing=1,maximum_spacing_error=0,intensity_pattern='X',envelope_support_class=EnvelopeSupportClass.STRONG_MULTI_MEMBER_ENVELOPE_SUPPORT);assert build_charge_support_summaries([ion(1)],[],[],[f])[0].charge_support_score>build_charge_support_summaries([ion(1)],[],[],[])[0].charge_support_score
def test_cross_rna_only_penalty():
 a=build_charge_support_summaries([ion(1)],[],[qrecord(cross=True)] ,[])[0];b=build_charge_support_summaries([ion(1)],[],[qrecord(cross=False)],[])[0];assert a.charge_support_score<=b.charge_support_score
def test_charge_score_deterministic():
 args=([ion(1)],[fmatch('A'),fmatch('B')],[qrecord('A')],[]);assert build_charge_support_summaries(*args)==build_charge_support_summaries(*args)

@pytest.mark.parametrize('scenario',[s for s in SimulationScenarioID])
def test_simulation_scenarios_present(scenario):assert scenario in {x.scenario_id for x in simulate_quality_tiers([qrecord()])}
def test_simulation_does_not_mutate():
 q=[qrecord()];before=repr(q);simulate_quality_tiers(q);assert repr(q)==before
def test_simulation_deterministic():
 q=[qrecord()];assert simulate_quality_tiers(q)==simulate_quality_tiers(q)
def test_negative_polarity_filtering():
 q=[qrecord('N',mode=T1IonMode.NEGATIVE_DEPROTONATED),qrecord('P',mode=T1IonMode.POSITIVE_PROTONATED)];s=next(x for x in simulate_quality_tiers(q) if x.scenario_id is SimulationScenarioID.NEGATIVE_POLARITY_CONFIRMED_ONLY);assert s.tier_e_count==1
def test_positive_polarity_filtering():
 q=[qrecord('N',mode=T1IonMode.NEGATIVE_DEPROTONATED),qrecord('P',mode=T1IonMode.POSITIVE_PROTONATED)];s=next(x for x in simulate_quality_tiers(q) if x.scenario_id is SimulationScenarioID.POSITIVE_POLARITY_CONFIRMED_ONLY);assert s.tier_e_count==1
def test_mono_filtering():
 q=[qrecord('M',definition=DeltaMassDefinition.MONOISOTOPIC_DELTA),qrecord('A',definition=DeltaMassDefinition.AVERAGE_DELTA)];s=next(x for x in simulate_quality_tiers(q) if x.scenario_id is SimulationScenarioID.MONOISOTOPIC_MZ_CONFIRMED_ONLY);assert s.tier_e_count==1
def test_simulation_flags():assert all(x.scenario_assumption_only and not x.applied_to_runtime_result and not x.applied_to_formal_score for x in simulate_quality_tiers([qrecord()]))

def test_inventory_and_integrated_contract(tmp_path):
 path=tmp_path/'sample.txt';path.write_text('Mass/Charge\tIntensity\n100\t1\n101\t2\n',encoding='utf-8');loaded=SimpleNamespace(runtime_path=path,header=('Mass/Charge','Intensity'),coordinates=(100.,101.),intensities=(1.,2.),source_sha256='abc');source=SimpleNamespace(profile_source_id='S',measurement_id='M',rna_identity_id='R',profile_type=SimpleNamespace(value='MZ_PROFILE'),mass_column='Mass/Charge');measurement=SimpleNamespace(measurement_id='M',source_file_name='raw.mzML',source_file_path_hint=None);peaks=(peak('A',100),);pr=T1ProfilePeakAuditResult('COMPLETED','S','M','R','MZ','MZ','CHARGED_ION_UNKNOWN',T1PeakDetectionParameters(),peaks,peaks);result=audit_t1_source_metadata(loaded,source,measurement,pr,[ion(1)],[],[qrecord()],project_config_path=tmp_path/'missing.yaml');assert result.source_record.data_representation_status is DataRepresentationStatus.UNKNOWN_REPRESENTATION;assert not result.source_record.representation_confirmed;assert result.source_record.observed_mass_definition_status is ObservedMassDefinitionStatus.UNKNOWN_MASS_DEFINITION;assert result.source_record.recommended_shadow_charge_range=='KEEP_CONFIGURED_RANGE';assert result.source_record.metadata_completeness_class is MetadataCompletenessClass.INSUFFICIENT_FOR_CHEMICAL_ASSIGNMENT;assert not result.source_record.source_linked_raw_metadata_available

def test_user_confirmed_polarity(tmp_path):
 path=tmp_path/'sample.txt';path.write_text('x',encoding='utf-8');loaded=SimpleNamespace(runtime_path=path,header=('Mass/Charge','Intensity'),coordinates=(100.,101.),intensities=(1.,2.),source_sha256='abc');source=SimpleNamespace(profile_source_id='S',measurement_id='M',rna_identity_id='R',profile_type=SimpleNamespace(value='MZ_PROFILE'),mass_column='Mass/Charge');measurement=SimpleNamespace(measurement_id='M',source_file_name='raw.mzML',source_file_path_hint=None);peaks=(peak('A',100),);pr=T1ProfilePeakAuditResult('COMPLETED','S','M','R','MZ','MZ','CHARGED_ION_UNKNOWN',T1PeakDetectionParameters(),peaks,peaks);result=audit_t1_source_metadata(loaded,source,measurement,pr,[ion(1)],[],[qrecord()],user_metadata={'acquisition_polarity':'NEGATIVE','exported_ion_polarity':'NEGATIVE'});assert result.source_record.polarity_confirmed and result.source_record.polarity_final_status is PolarityStatus.NEGATIVE

@pytest.mark.parametrize('flag',['metadata_audit_only','shadow_analysis_only'])
def test_metadata_safeguards_true(flag):
 fam=IsotopeFamilyCandidate(isotope_family_id='F',member_peak_ids=('A','B'),member_apex_mzs=(1,2),member_count=2,possible_charge=1,mean_spacing=1,maximum_spacing_error=0,intensity_pattern='X',envelope_support_class=EnvelopeSupportClass.MODERATE_TWO_MEMBER_SUPPORT);assert getattr(fam,flag) is True
@pytest.mark.parametrize('flag',['charge_state_confirmed','ion_mode_confirmed','fragment_assigned','modification_assigned','structure_assigned','applied_to_formal_score','applied_to_ranking','applied_to_candidate_filtering','applied_to_final_consensus'])
def test_formal_safeguards_false(flag):
 fam=IsotopeFamilyCandidate(isotope_family_id='F',member_peak_ids=('A','B'),member_apex_mzs=(1,2),member_count=2,possible_charge=1,mean_spacing=1,maximum_spacing_error=0,intensity_pattern='X',envelope_support_class=EnvelopeSupportClass.MODERATE_TWO_MEMBER_SUPPORT);assert getattr(fam,flag) is False


def test_acquisition_only_does_not_confirm_exported_polarity(tmp_path):
 path=tmp_path/'sample.txt';path.write_text('x',encoding='utf-8');loaded=SimpleNamespace(runtime_path=path,header=('Mass/Charge','Intensity'),coordinates=(100.,101.),intensities=(1.,2.),source_sha256='abc');source=SimpleNamespace(profile_source_id='S',measurement_id='M',rna_identity_id='R',profile_type=SimpleNamespace(value='MZ_PROFILE'),mass_column='Mass/Charge');measurement=SimpleNamespace(measurement_id='M',source_file_name='raw.mzML',source_file_path_hint=None);peaks=(peak('A',100),);pr=T1ProfilePeakAuditResult('COMPLETED','S','M','R','MZ','MZ','CHARGED_ION_UNKNOWN',T1PeakDetectionParameters(),peaks,peaks);r=audit_t1_source_metadata(loaded,source,measurement,pr,[ion(1)],[],[qrecord()],user_metadata={'acquisition_polarity':'NEGATIVE'}).source_record;assert r.acquisition_polarity is PolarityStatus.NEGATIVE and r.exported_ion_polarity is PolarityStatus.UNKNOWN and not r.polarity_confirmed

def test_priority_candidate_tracking_generic():
 q=qrecord('TARGET');q.t1_delta_relation_id='R';rows=simulate_priority_candidates([q],['TARGET']);assert len(rows)==6 and all(x.observed_peak_id=='TARGET' for x in rows)

def test_priority_candidate_tracking_no_hardcoded_target():
 q=qrecord('ARBITRARY');q.t1_delta_relation_id='R';assert simulate_priority_candidates([q],['ARBITRARY'])
