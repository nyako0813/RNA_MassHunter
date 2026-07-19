from dataclasses import replace
import numpy as np
import pytest
from rna_masshunter.sciex_t1_profile_peak_audit import *

def profile(centers=(100.0,),amps=(100.0,),sigma=.012):
 x=np.arange(98,104,.002);y=np.ones_like(x)*.01
 for c,a in zip(centers,amps):y+=a*np.exp(-.5*((x-c)/sigma)**2)
 return tuple(x),tuple(y)
def detect(centers=(100.,),amps=(100.,),**kwargs):
 x,y=profile(centers,amps,kwargs.pop('sigma',.012));p=T1PeakDetectionParameters(minimum_relative_apex_intensity=1e-5,minimum_relative_integrated_intensity=1e-5,minimum_relative_prominence=1e-4,minimum_peak_separation_mz=.01,**kwargs)
 return detect_t1_profile_peaks(x,y,source_id='S',measurement_id='M',rna_identity='R',parameters=p)
def test_mz_model_and_peak_metrics():
 r=detect();p=r.peaks[0]
 assert (r.coordinate_type,r.observed_mass_scale,r.observed_output_species)==('MZ','MZ','CHARGED_ION_UNKNOWN')
 assert p.apex_mz==pytest.approx(100,abs=.002);assert p.centroid_mz==pytest.approx(100,abs=.002)
 assert p.integrated_intensity>0 and p.fwhm_mz>0 and p.prominence>0 and p.sharpness_score>0
 assert p.relative_apex_intensity==pytest.approx(1,rel=1e-3);assert p.relative_integrated_intensity>0
@pytest.mark.parametrize('charge',[1,2,3,4,5])
def test_isotope_spacing_diagnostic(charge):
 spacing=ISOTOPE_MASS_DIFFERENCE_DA/charge;r=detect((100,100+spacing),(100,50),isotope_spacing_tolerance_mz=.004)
 flagged=[p for p in r.peaks if p.possible_isotope_component]
 assert flagged and any(p.possible_isotope_charge==charge for p in flagged)
 assert all(p.isotope_assigned is False for p in flagged)
def test_quality_and_bounded_selection():
 centers=tuple(99+i*.08 for i in range(40));x=np.arange(98,104,.002);y=np.ones_like(x)*.01
 for i,c in enumerate(centers):y+=(100-i)*np.exp(-.5*((x-c)/.008)**2)
 params=T1PeakDetectionParameters(minimum_relative_apex_intensity=1e-5,minimum_relative_integrated_intensity=1e-6,minimum_relative_prominence=1e-4,minimum_peak_separation_mz=.01,maximum_selected_peaks_per_profile=10,major_relative_intensity=.01)
 r=detect_t1_profile_peaks(tuple(x),tuple(y),source_id='S',measurement_id='M',rna_identity='R',parameters=params)
 assert len(r.selected_peaks)<=10;assert any(p.peak_quality_class in {T1PeakQualityClass.MAJOR_SHARP,T1PeakQualityClass.MAJOR_BROAD} for p in r.peaks)
def test_shoulder_duplicate_and_overlap_flags():
 r=detect((100,100.02),(100,90),sigma=.006,shared_valley_fraction=.1,duplicate_separation_mz=.04)
 assert any(p.possible_shoulder for p in r.peaks) or any(p.possible_duplicate for p in r.peaks)
def test_determinism_and_nonmutation():
 x,y=profile((100,101),(100,50));before=(x,y);a=detect_t1_profile_peaks(x,y,source_id='S',measurement_id='M',rna_identity='R');b=detect_t1_profile_peaks(x,y,source_id='S',measurement_id='M',rna_identity='R')
 assert a==b and (x,y)==before
@pytest.mark.parametrize('bad',[((),()),((1,2),(1,)),((2,1),(1,2))])
def test_invalid_input_rejected(bad):
 with pytest.raises(ValueError):detect_t1_profile_peaks(*bad,source_id='S',measurement_id='M',rna_identity='R')
def test_all_false_certainty_and_formal_flags():
 p=detect().peaks[0]
 for name in ('rna_identity_confirmed','fragment_assigned','modification_assigned','position_assigned','structure_assigned','charge_state_confirmed','ion_mode_confirmed','polarity_confirmed','applied_to_formal_score','applied_to_ranking','applied_to_candidate_filtering','applied_to_final_consensus'):assert getattr(p,name) is False

@pytest.mark.parametrize('field',['apex_mz','centroid_mz','integrated_intensity','fwhm_mz','prominence','sharpness_score','relative_apex_intensity'])
def test_required_peak_metric_is_populated(field):
 assert getattr(detect().peaks[0],field) is not None

def test_loaded_profile_header_validation():
 from types import SimpleNamespace
 loaded=SimpleNamespace(profile_source_id='S',header=('wrong','Intensity'),coordinates=(1,2,3),intensities=(0,1,0));source=SimpleNamespace(profile_source_id='S',measurement_id='M',rna_identity_id='R')
 with pytest.raises(ValueError,match='header'):analyze_loaded_t1_profile(loaded,source)

def test_loaded_profile_source_identity_validation():
 from types import SimpleNamespace
 loaded=SimpleNamespace(profile_source_id='X',header=('Mass/Charge','Intensity'),coordinates=(1,2,3),intensities=(0,1,0));source=SimpleNamespace(profile_source_id='S',measurement_id='M',rna_identity_id='R')
 with pytest.raises(ValueError,match='conflict'):analyze_loaded_t1_profile(loaded,source)
