from dataclasses import replace
import pytest
from rna_masshunter.sciex_t1_profile_peak_audit import *
from rna_masshunter.sciex_t1_cross_profile_audit import *

def peak(prefix,i,mz,rel=.5,iso=False,fwhm=.02,width=.06,prom=50,cent=.001,quality=T1PeakQualityClass.MAJOR_SHARP):
 return T1ProfilePeak(t1_peak_id=f'{prefix}{i}',source_id=prefix,measurement_id=prefix,rna_identity=prefix,apex_mz=mz,centroid_mz=mz+cent,apex_intensity=rel*100,integrated_intensity=rel*20,relative_apex_intensity=rel,relative_integrated_intensity=rel,left_boundary_mz=mz-width/2,right_boundary_mz=mz+width/2,peak_width_mz=width,fwhm_mz=fwhm,prominence=prom,relative_prominence=.1,sharpness_score=prom/fwhm,nearest_peak_separation_mz=None,shared_valley_fraction=0,possible_shoulder=False,possible_duplicate=False,possible_isotope_component=iso,possible_isotope_charge=2 if iso else None,possible_isotope_spacing_error=0 if iso else None,possible_overlapping_envelope=iso,peak_quality_class=quality,selected_as_primary=True)
def result(prefix,peaks):
 p=T1PeakDetectionParameters();v=tuple(peaks);return T1ProfilePeakAuditResult('COMPLETED',prefix,prefix,prefix,'MZ','MZ','CHARGED_ION_UNKNOWN',p,v,v)
def match(a,b,p=None):return match_t1_cross_profile(a,b,layer=T1ComparisonLayer.SELECTED_T1_PEAKS,parameters=p)
@pytest.mark.parametrize('d,klass',[(.01,T1MassMatchClass.STRICT),(.015,T1MassMatchClass.EXPLORATORY)])
def test_strict_exploratory(d,klass):assert match([peak('A',1,100)],[peak('B',1,100+d)])[0].mass_match_class is klass
def test_maximum_cardinality_minimum_error_unique_and_deterministic():
 a=[peak('A',1,100),peak('A',2,100.018)];b=[peak('B',1,100.001),peak('B',2,100.019)]
 x=match(a,b);y=match(reversed(a),reversed(b));assert x==y and len(x)==2;assert len({m.uaa_peak_id for m in x})==2 and sum(m.apex_mz_difference for m in x)==pytest.approx(.002)
def test_unmatched_and_classifications_and_fractions():
 r=audit_t1_cross_profiles(result('A',[peak('A',1,100,.75),peak('A',2,200,.25)]),result('B',[peak('B',1,100,.6),peak('B',2,300,.4)]))
 classes={s.classification for s in r.statuses};assert T1SelectedClassification.COMMON_T1_SELECTED_PEAK in classes;assert T1SelectedClassification.UAA_ONLY_T1_SELECTED_PEAK in classes;assert T1SelectedClassification.UAG_ONLY_T1_SELECTED_PEAK in classes
 assert r.uaa_summary.common_apex_intensity_fraction==pytest.approx(.75);assert r.uag_summary.sample_specific_apex_intensity_fraction==pytest.approx(.4)
def test_isotope_only_match():
 r=audit_t1_cross_profiles(result('A',[peak('A',1,100,iso=True)]),result('B',[peak('B',1,100)]));assert all(s.classification is T1SelectedClassification.ISOTOPE_OR_ENVELOPE_ONLY_MATCH for s in r.statuses);assert r.selected_matches[0].shape_similarity_class is T1ShapeSimilarityClass.ISOTOPE_OR_ENVELOPE_AMBIGUOUS
@pytest.mark.parametrize('kwargs,klass',[({},T1ShapeSimilarityClass.HIGHLY_SIMILAR_T1_PEAK),({'rel':.2,'fwhm':.03,'width':.08},T1ShapeSimilarityClass.MODERATELY_SIMILAR_T1_PEAK),({'rel':.001,'fwhm':.3,'width':.6,'prom':1,'quality':T1PeakQualityClass.MAJOR_BROAD},T1ShapeSimilarityClass.MZ_MATCH_SHAPE_DIFFERENT)])
def test_shape_classes(kwargs,klass):assert match([peak('A',1,100)],[peak('B',1,100.001,**kwargs)])[0].shape_similarity_class is klass
def test_match_fields_and_formal_flags():
 m=match([peak('A',1,100,.5)],[peak('B',1,100.005,.8,fwhm=.04,width=.08,prom=80)])[0]
 assert m.centroid_mz_difference==pytest.approx(.005);assert m.uaa_fwhm==.02 and m.uag_fwhm==.04;assert m.uaa_prominence==50 and m.uag_prominence==80;assert m.uaa_relative_apex_intensity==.5
 for n in ('rna_identity_confirmed','fragment_assigned','position_assigned','structure_assigned','applied_to_formal_score','applied_to_final_consensus'):assert getattr(m,n) is False
def test_correlations_and_three_layers():
 a=[peak('A',i,100*i,i/4) for i in range(1,4)];b=[peak('B',i,100*i,i/4) for i in range(1,4)];r=audit_t1_cross_profiles(result('A',a),result('B',b));assert r.correlations.apex_mz==pytest.approx(1);assert len(r.all_detected_matches)==len(r.selected_matches)==len(r.non_isotope_selected_matches)==3

@pytest.mark.parametrize('field',['centroid_mz_difference','uaa_fwhm','uag_fwhm','uaa_prominence','uag_prominence','uaa_relative_integrated_intensity'])
def test_required_cross_profile_field_is_populated(field):
 m=match([peak('A',1,100)],[peak('B',1,100.001)])[0];assert getattr(m,field) is not None

@pytest.mark.parametrize('flag',['applied_to_formal_score','applied_to_ranking','applied_to_candidate_filtering','applied_to_final_consensus'])
def test_cross_profile_formal_nonpropagation(flag):
 m=match([peak('A',1,100)],[peak('B',1,100)])[0];assert getattr(m,flag) is False
