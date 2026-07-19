from dataclasses import replace
from types import SimpleNamespace
import pytest
from rna_masshunter.sciex_intact_peak_family import DeltaMassDefinition
from rna_masshunter.sciex_t1_profile_peak_audit import T1ProfilePeak,T1PeakQualityClass
from rna_masshunter.sciex_t1_fragment_shadow_match import T1IonCandidate,T1IonMode,TheoreticalFragmentDiscrimination,TheoreticalDiscriminationClass
from rna_masshunter.sciex_t1_fragment_delta_audit import *
from rna_masshunter.sciex_t1_delta_evidence_quality_audit import *
from rna_masshunter.sciex_t1_delta_evidence_quality_audit import _concordance

def peak(pid='P',mz=116,centroid=None,quality=T1PeakQualityClass.MAJOR_SHARP,iso=False,shoulder=False,duplicate=False,overlap=False):
 return T1ProfilePeak(t1_peak_id=pid,source_id='LEU_UAA_WT_T1_MZ',measurement_id='M',rna_identity='TRNA_LEU_UAA',apex_mz=mz,centroid_mz=mz if centroid is None else centroid,apex_intensity=100,integrated_intensity=10,relative_apex_intensity=1,relative_integrated_intensity=1,left_boundary_mz=mz-.03,right_boundary_mz=mz+.03,peak_width_mz=.06,fwhm_mz=.02,prominence=90,relative_prominence=.9,sharpness_score=4500,nearest_peak_separation_mz=None,shared_valley_fraction=0,possible_shoulder=shoulder,possible_duplicate=duplicate,possible_isotope_component=iso,possible_isotope_charge=1 if iso else None,possible_isotope_spacing_error=0 if iso else None,possible_overlapping_envelope=overlap or iso,peak_quality_class=quality,selected_as_primary=True)
def ion(fid='F',mz=100,charge=1,mode=T1IonMode.NEGATIVE_DEPROTONATED,rna='TRNA_LEU_UAA',seq='AG',cca='CCA'):
 return T1IonCandidate(f'{fid}_{rna}_{cca}_{mode.value}_{charge}',fid,rna,rna+'_'+cca,cca,seq,1,len(seq),mode,charge,mz,'MONOISOTOPIC_NEUTRAL','MONOISOTOPIC_ION_MZ',-charge)
def one_ref(rid='O_ADDITION_AVERAGE'):
 return next(r for r in build_chemical_delta_reference_registry() if r.reference_id==rid)
def delta(peaks=None,ions=None,refs=None,mods=(),source='LEU_UAA_WT_T1_MZ'):
 return audit_t1_fragment_chemical_deltas(peaks or [peak()],ions or [ion()],source_id=source,references=refs,modifications=mods)
def quality(d=None,peaks=None,status=ObservedPolarityStatus.POLARITY_UNKNOWN,all_rows=(),discrimination=(),params=None):
 d=d or delta(refs=(one_ref(),));peaks=peaks or [peak()]
 return audit_t1_delta_evidence_quality(d,peaks,measurement_id='M',polarity_status=status,all_sample_relations=all_rows,discrimination=discrimination,parameters=params)
def record(**kwargs):return quality(**kwargs).records[0]

@pytest.mark.parametrize('values,expected',[(('negative',),ObservedPolarityStatus.POLARITY_CONFIRMED_NEGATIVE),(('positive',),ObservedPolarityStatus.POLARITY_CONFIRMED_POSITIVE),((),ObservedPolarityStatus.POLARITY_UNKNOWN),(('negative','positive'),ObservedPolarityStatus.POLARITY_CONFLICT)])
def test_polarity_metadata(values,expected):assert determine_observed_polarity(values) is expected
@pytest.mark.parametrize('status,mode,support',[ (ObservedPolarityStatus.POLARITY_CONFIRMED_NEGATIVE,T1IonMode.NEGATIVE_DEPROTONATED,PolaritySupportClass.CONFIRMED_COMPATIBLE),(ObservedPolarityStatus.POLARITY_CONFIRMED_NEGATIVE,T1IonMode.POSITIVE_PROTONATED,PolaritySupportClass.CONFIRMED_INCOMPATIBLE),(ObservedPolarityStatus.POLARITY_UNKNOWN,T1IonMode.NEGATIVE_DEPROTONATED,PolaritySupportClass.UNKNOWN_POLARITY),(ObservedPolarityStatus.POLARITY_CONFLICT,T1IonMode.NEGATIVE_DEPROTONATED,PolaritySupportClass.CONFLICTING_METADATA)])
def test_polarity_support(status,mode,support):
 d=delta(ions=[ion(mode=mode)],refs=(one_ref(),));assert record(d=d,status=status).polarity_support_class is support

def test_mass_definition_unknown_and_tier_a_blocked():
 r=record(status=ObservedPolarityStatus.POLARITY_CONFIRMED_NEGATIVE);assert r.mass_definition_compatibility_status is MassDefinitionCompatibilityStatus.UNKNOWN_COMPATIBILITY;assert r.evidence_tier is not EvidenceTier.TIER_A_HIGH_SUPPORT

def test_mass_definition_mismatch_known_mod():
 mod=SimpleNamespace(id='MODX',symbol=None,mass_shift_from_unmodified=16.0);d=delta(mods=[mod]);r=record(d=d);assert r.mass_definition_compatibility_status is MassDefinitionCompatibilityStatus.CONFIRMED_MISMATCH and r.evidence_tier is EvidenceTier.TIER_D_DIAGNOSTIC_ONLY

def test_mass_definition_confirmed_compatible_synthetic():
 d=delta(refs=(one_ref(),));d=replace(d,relations=tuple(replace(r,mass_definition_compatible=True) for r in d.relations));assert record(d=d).mass_definition_compatibility_status is MassDefinitionCompatibilityStatus.CONFIRMED_COMPATIBLE

def test_tier_a_blocked_by_unknown_polarity():
 d=delta(refs=(one_ref(),));d=replace(d,relations=tuple(replace(r,mass_definition_compatible=True) for r in d.relations));assert record(d=d).evidence_tier is not EvidenceTier.TIER_A_HIGH_SUPPORT

def multi_charge(conflict=False,same_peak=False):
 ps=[peak('P1',115.9994),peak('P1' if same_peak else 'P2',107.9997)]
 ions=[ion('F1',100,1),ion('F1',100,2)]
 d=delta(peaks=ps,ions=ions,refs=(one_ref(),));
 if conflict:d=replace(d,relations=tuple(replace(r,reference_id=r.reference_id+str(r.charge)) for r in d.relations))
 return d,ps
@pytest.mark.parametrize('n',[2,3])
def test_multi_charge_concordant(n):
 ps=[peak(f'P{z}',100+15.9994/z) for z in range(1,n+1)];ions=[ion('F',100,z) for z in range(1,n+1)];d=delta(peaks=ps,ions=ions,refs=(one_ref(),));assert all(r.charge_support_class is ChargeSupportClass.MULTI_CHARGE_CONCORDANT_SUPPORT for r in quality(d,ps).records)
def test_same_peak_duplicate_charge_not_independent():
 d,ps=multi_charge(same_peak=True);assert all(r.charge_support_class is not ChargeSupportClass.MULTI_CHARGE_CONCORDANT_SUPPORT for r in quality(d,ps).records)
def test_single_charge_only():assert record().charge_support_class is ChargeSupportClass.SINGLE_CHARGE_ONLY
def test_multi_charge_conflicting():
 d,ps=multi_charge(conflict=True);assert all(r.charge_support_class is ChargeSupportClass.MULTI_CHARGE_CONFLICTING for r in quality(d,ps).records)
def test_cca_duplicate_collapsed():
 d=delta(ions=[ion(cca='C'),ion(cca='CCA')],refs=(one_ref(),));r=record(d=d);assert r.cca_collapsed_candidate_count==1 and r.fragment_identity_collapsed_count==1
def test_charge_support_deterministic():
 d,ps=multi_charge();assert quality(d,ps)==quality(d,ps)

@pytest.mark.parametrize('centroid,klass',[(115.9994,ApexCentroidConcordanceClass.BOTH_STRICT_SAME_REFERENCE),(116.0694,ApexCentroidConcordanceClass.APEX_STRICT_CENTROID_EXPLORATORY_SAME_REFERENCE),(116.0794,ApexCentroidConcordanceClass.APEX_STRICT_CENTROID_EXPLORATORY_SAME_REFERENCE),(116.2,ApexCentroidConcordanceClass.APEX_ONLY_SUPPORT)])
def test_apex_centroid_classes(centroid,klass):
 p=peak(centroid=centroid);assert record(d=delta([p],refs=(one_ref(),)),peaks=[p]).apex_centroid_concordance_class is klass
def test_both_exploratory():
 p=peak(mz=116.0794,centroid=116.0794);assert record(d=delta([p],refs=(one_ref(),)),peaks=[p]).apex_centroid_concordance_class is ApexCentroidConcordanceClass.BOTH_EXPLORATORY_SAME_REFERENCE
def test_centroid_only_support():
 p=peak(mz=116.2,centroid=115.9994);d=delta([p],refs=(one_ref(),));assert _concordance(d.relations[0],d.references,d.parameters)[-1] is ApexCentroidConcordanceClass.CENTROID_ONLY_SUPPORT

def test_reference_category_conflict():
 refs=(one_ref('O_ADDITION_AVERAGE'),one_ref('O_TO_S_MONOISOTOPIC'));p=peak(mz=115.9994,centroid=115.97715655483);r=record(d=delta([p],refs=refs),peaks=[p]);assert r.apex_centroid_concordance_class is ApexCentroidConcordanceClass.REFERENCE_CATEGORY_CONFLICT
def test_no_concordance_with_missing_centroid():
 p=replace(peak(),centroid_mz=None);r=record(d=delta([p],refs=(one_ref(),)),peaks=[p]);assert r.apex_centroid_concordance_class is ApexCentroidConcordanceClass.APEX_ONLY_SUPPORT
@pytest.mark.parametrize('quality_class,expected',[(T1PeakQualityClass.MAJOR_SHARP,PeakShapeSupportClass.STRONG_PEAK_SHAPE_SUPPORT),(T1PeakQualityClass.MINOR_SHARP,PeakShapeSupportClass.MODERATE_PEAK_SHAPE_SUPPORT),(T1PeakQualityClass.MAJOR_BROAD,PeakShapeSupportClass.WEAK_PEAK_SHAPE_SUPPORT)])
def test_peak_shape(quality_class,expected):
 p=peak(quality=quality_class);assert record(d=delta([p],refs=(one_ref(),)),peaks=[p]).peak_shape_support_class is expected

def test_low_fragment_ambiguity():assert record().fragment_ambiguity_class is FragmentAmbiguityClass.LOW_FRAGMENT_AMBIGUITY
def test_cca_only_moderate_ambiguity():
 d=delta(ions=[ion(cca='C'),ion(cca='CCA')],refs=(one_ref(),));assert record(d=d).fragment_ambiguity_class is FragmentAmbiguityClass.MODERATE_FRAGMENT_AMBIGUITY
def test_multiple_fragment_high_ambiguity():
 d=delta(ions=[ion('F1'),ion('F2',seq='CG')],refs=(one_ref(),));assert record(d=d).fragment_ambiguity_class is FragmentAmbiguityClass.HIGH_FRAGMENT_AMBIGUITY
def test_cross_rna_ambiguity():
 d=delta(ions=[ion(rna='TRNA_LEU_UAA'),ion(rna='TRNA_LEU_UAG')],refs=(one_ref(),));assert record(d=d).fragment_ambiguity_class is FragmentAmbiguityClass.CROSS_RNA_AMBIGUITY

def test_single_reference_resolution():assert record().reference_resolution_class is ReferenceResolutionClass.SINGLE_REFERENCE_ONLY
@pytest.mark.parametrize('margin,expected',[(.06,ReferenceResolutionClass.WELL_RESOLVED_REFERENCE_CATEGORY),(.03,ReferenceResolutionClass.MARGINALLY_RESOLVED_REFERENCE_CATEGORY),(.01,ReferenceResolutionClass.UNRESOLVED_REFERENCE_CATEGORY)])
def test_resolution_thresholds(margin,expected):
 base=one_ref();other=replace(base,reference_id='OTHER',reference_category=ChemicalReferenceCategory.OTHER_DIAGNOSTIC,signed_delta_da=base.signed_delta_da+margin);r=record(d=delta(refs=(base,other)));assert r.reference_resolution_class is expected

def test_target_discrimination():
 disc=[TheoreticalFragmentDiscrimination('F','TRNA_LEU_UAA','AG',100,TheoreticalDiscriminationClass.UAA_SPECIFIC_THEORETICAL_FRAGMENT)];assert record(discrimination=disc).discriminatory_support_class is DiscriminatorySupportClass.TARGET_SPECIFIC_FRAGMENT_SUPPORT
@pytest.mark.parametrize('cls,expected',[(TheoreticalDiscriminationClass.SEQUENCE_IDENTICAL_FRAGMENT,DiscriminatorySupportClass.SHARED_FRAGMENT_SUPPORT),(TheoreticalDiscriminationClass.SAME_MASS_DIFFERENT_SEQUENCE,DiscriminatorySupportClass.SAME_MASS_SEQUENCE_AMBIGUITY)])
def test_other_discrimination(cls,expected):
 disc=[TheoreticalFragmentDiscrimination('F','TRNA_LEU_UAA','AG',100,cls)];assert record(discrimination=disc).discriminatory_support_class is expected

@pytest.mark.parametrize('series,expected',[(None,StateSeriesSupportClass.NO_SERIES_SUPPORT)])
def test_no_series(series,expected):assert record().state_series_support_class is expected
def series(pattern=T1StateSeriesPattern.SINGLE_O_STEP,members=('P','Q','R'),strict=2,expl=0):
 return T1ChemicalStateSeries(t1_state_series_id='S',rna_identity_candidate='TRNA_LEU_UAA',theoretical_fragment_identity='TRNA_LEU_UAA:1-2:AG',member_observed_peak_ids=members,member_apex_mzs=tuple(range(len(members))),member_centroid_mzs=tuple(range(len(members))),member_count=len(members),charge=1,ion_mode=T1IonMode.NEGATIVE_DEPROTONATED,o_edge_count=1,h2o_edge_count=0,s_edge_count=0,o_to_s_edge_count=0,s_to_o_edge_count=0,strict_edge_count=strict,exploratory_edge_count=expl,series_pattern=pattern,mass_span_neutral_da=16,sequential_oxygen_equivalent_series_detected=False,oxidation_state_series_possible=False)
@pytest.mark.parametrize('s,expected',[(series(),StateSeriesSupportClass.STRONG_LINEAR_SERIES_SUPPORT),(series(members=('P','Q'),strict=1),StateSeriesSupportClass.SINGLE_EDGE_ONLY),(series(strict=1,expl=2),StateSeriesSupportClass.MODERATE_SERIES_SUPPORT),(series(T1StateSeriesPattern.BRANCHED_SERIES),StateSeriesSupportClass.BRANCHED_AMBIGUOUS_SERIES)])
def test_series_support(s,expected):
 d=delta(refs=(one_ref(),));d=replace(d,state_series=(s,));assert record(d=d).state_series_support_class is expected

def test_recurrent_multi_peak_single_charge():
 ps=[peak('P1'),peak('P2',116.001)];d=delta(ps,refs=(one_ref(),));assert all(r.recurrent_support_class is RecurrentSupportClass.MULTI_PEAK_SINGLE_CHARGE_RECURRENT for r in quality(d,ps).records)
def test_recurrent_multi_peak_multi_charge():
 d,ps=multi_charge();assert all(r.recurrent_support_class is RecurrentSupportClass.MULTI_PEAK_MULTI_CHARGE_RECURRENT for r in quality(d,ps).records)
def test_single_peak_recurrent():assert record().recurrent_support_class is RecurrentSupportClass.SINGLE_PEAK_ONLY
def test_cross_sample_recurrent():
 u=delta(refs=(one_ref(),));g=delta(refs=(one_ref(),),source='LEU_UAG_WT_T1_MZ');r=record(d=u,all_rows=u.relations+g.relations);assert r.recurrent_support_class is RecurrentSupportClass.CROSS_SAMPLE_RECURRENT

def test_tier_a_complete_synthetic():
 d,ps=multi_charge();d=replace(d,relations=tuple(replace(r,mass_definition_compatible=True) for r in d.relations));assert all(r.evidence_tier is EvidenceTier.TIER_A_HIGH_SUPPORT for r in quality(d,ps,status=ObservedPolarityStatus.POLARITY_CONFIRMED_NEGATIVE).records)
def test_tier_b_one_unknown():
 d,ps=multi_charge();d=replace(d,relations=tuple(replace(r,mass_definition_compatible=True) for r in d.relations));q=quality(d,ps,status=ObservedPolarityStatus.POLARITY_UNKNOWN);assert all(r.evidence_tier is EvidenceTier.TIER_B_MODERATE_SUPPORT for r in q.records)
def test_tier_c_weak():assert record().evidence_tier is EvidenceTier.TIER_C_WEAK_SUPPORT
def test_tier_d_diagnostic():
 base=one_ref();other=replace(base,reference_id='X',reference_category=ChemicalReferenceCategory.OTHER_DIAGNOSTIC,signed_delta_da=base.signed_delta_da+.005);assert record(d=delta(refs=(base,other))).evidence_tier is EvidenceTier.TIER_D_DIAGNOSTIC_ONLY
@pytest.mark.parametrize('kind',["isotope","polarity"])
def test_tier_e_blocked(kind):
 d=delta(refs=(one_ref(),));p=peak()
 if kind=='isotope':p=replace(p,possible_isotope_component=True,possible_overlapping_envelope=True)
 status=ObservedPolarityStatus.POLARITY_CONFIRMED_POSITIVE if kind=='polarity' else ObservedPolarityStatus.POLARITY_UNKNOWN
 assert record(d=d,peaks=[p],status=status).evidence_tier is EvidenceTier.TIER_E_BLOCKED
@pytest.mark.parametrize('penalty,expected',[(1000,0),(-1000,100)])
def test_score_clipping(penalty,expected):
 params=EvidenceQualityParameters(unknown_polarity_penalty=penalty,unknown_mass_definition_penalty=penalty,single_charge_penalty=penalty) if penalty>0 else EvidenceQualityParameters(both_strict_score=1000,strong_shape_score=1000,low_ambiguity_score=1000,resolved_reference_score=1000,unknown_polarity_penalty=0,unknown_mass_definition_penalty=0,single_charge_penalty=0);assert record(params=params).evidence_support_score==expected

def test_block_reason_priority():
 d=delta(refs=(one_ref(),));p=replace(peak(),possible_isotope_component=True,possible_overlapping_envelope=True);r=record(d=d,peaks=[p],status=ObservedPolarityStatus.POLARITY_CONFIRMED_POSITIVE);assert r.evidence_block_reason is EvidenceBlockReason.CONFIRMED_POLARITY_MISMATCH
def test_repeated_complete_match():assert quality()==quality()
def test_input_nonmutation():
 d=delta(refs=(one_ref(),));ps=[peak()];before=(repr(d),repr(ps));quality(d,ps);assert before==(repr(d),repr(ps))

@pytest.mark.parametrize('rid,hyp,priority',[('O_TO_S_AVERAGE',ChemicalHypothesisClass.O_TO_S_EQUIVALENT,PriorityOSClass.O_TO_S_PRIORITY_CANDIDATE),('S_TO_O_AVERAGE',ChemicalHypothesisClass.S_TO_O_EQUIVALENT,PriorityOSClass.S_TO_O_PRIORITY_CANDIDATE)])
def test_os_priority(rid,hyp,priority):
 ref=one_ref(rid);p=peak(mz=100+ref.signed_delta_da);q=quality(delta([p],refs=(ref,)),[p]);assert q.os_priority_records[0].quality_record.chemical_hypothesis_class is hyp;assert q.os_priority_records[0].priority_class is priority

def test_os_diagnostic_only():
 ref=one_ref('O_TO_S_AVERAGE');other=replace(ref,reference_id='X',reference_category=ChemicalReferenceCategory.OTHER_DIAGNOSTIC,signed_delta_da=ref.signed_delta_da+.001);p=peak(mz=100+ref.signed_delta_da);assert quality(delta([p],refs=(ref,other)),[p]).os_priority_records[0].priority_class is PriorityOSClass.O_TO_S_DIAGNOSTIC_ONLY

def test_blocked_isotope_os():
 ref=one_ref('O_TO_S_AVERAGE');p=peak(mz=100+ref.signed_delta_da);d=delta([p],refs=(ref,));p=replace(p,possible_isotope_component=True,possible_overlapping_envelope=True);assert quality(d,[p]).os_priority_records[0].priority_class is PriorityOSClass.BLOCKED_O_S_CANDIDATE

def test_cross_rna_ambiguous_os():
 ref=one_ref('O_TO_S_AVERAGE');p=peak(mz=100+ref.signed_delta_da);d=delta([p],[ion(rna='TRNA_LEU_UAA'),ion(rna='TRNA_LEU_UAG')],refs=(ref,));assert all(r.cross_rna_identity_ambiguous for r in quality(d,[p]).records)
def test_alternative_charge_support():
 ref=one_ref('O_TO_S_AVERAGE');ps=[peak('P1',100+ref.signed_delta_da),peak('P2',100+ref.signed_delta_da/2)];d=delta(ps,[ion('F',100,1),ion('F',100,2)],refs=(ref,));assert all(r.alternative_charge_support_class=='CONCORDANT_OBSERVED' for r in quality(d,ps).records)
def test_no_alternative_charge_support():assert record().alternative_charge_support_class=='NOT_OBSERVED'
@pytest.mark.parametrize('flag',['thioamide_assigned','sulfur_atom_assigned','o_to_s_substitution_assigned','s_to_o_substitution_assigned','modification_position_assigned','structure_assigned','charge_state_confirmed','polarity_confirmed','applied_to_formal_score','applied_to_ranking','applied_to_candidate_filtering','applied_to_final_consensus'])
def test_safeguards(flag):assert getattr(record(),flag) is False

def test_glu_summary_flags_no_direct_link():
 q=audit_t1_delta_evidence_quality(delta(refs=(one_ref(),)),[peak()],measurement_id='M',glu_summary_flags={'oxygen':True,'water':True});assert q.oxygen_equivalent_pattern_also_observed_in_glu_intact and q.water_equivalent_pattern_also_observed_in_glu_intact and not q.direct_chemical_identity_link_assigned


def test_existing_cross_profile_match_marks_common_component_possible():
 d=delta(refs=(one_ref(),));match=SimpleNamespace(uaa_peak_id='P',uag_peak_id='G');q=audit_t1_delta_evidence_quality(d,[peak()],measurement_id='M',cross_profile_matches=[match]);assert q.records[0].common_component_or_shared_fragment_possible is True
