from copy import deepcopy
from types import SimpleNamespace

import pytest

from rna_masshunter.models import Peak
from rna_masshunter.p1_sap_dinucleotide_candidates import generate_dinucleotide_candidates
from rna_masshunter.p1_sap_dinucleotide_feature_audit import (
    _competition_rows, _map_isotopes, audit_dinucleotide_features,
    build_ms2_provenance, classify_mass_accuracy, match_raw_groups,
)
from rna_masshunter.p1_sap_dinucleotide_interpretation import (
    build_p1_sap_dinucleotide_audit, classify_features, interpret_groups,
)
from rna_masshunter.p1_sap_feature_quality import _expected_isotope_abundance


def cfg(*, charges=(1,), polarity="positive", search=(100, 1000), tolerance=10, targets=None):
    return SimpleNamespace(
        p1_sap_dinucleotide={
            "candidate_generation": {"max_modifications_per_side": 0, "max_composite_states_per_position": 16, "max_candidate_count": 1000, "include_normal_phosphate": True, "include_phosphorothioate": True, "charges": list(charges), "polarity": polarity},
            "search": {"mz_min": search[0], "mz_max": search[1], "tolerance_ppm": tolerance},
            "mass_accuracy": {"strong_ppm": 2, "moderate_ppm": 5, "search_ppm": tolerance},
            "feature_quality": {"min_spectrum_count": 2, "min_profile_point_count": 2, "max_rt_gap_min": 0.08, "background_window_rt_min": 0.5, "background_mz_tolerance_ppm": 10},
            "isotope": {"enabled": True, "tolerance_ppm": 20, "require_same_scan": True},
            "ms2_provenance": {"enabled": True}, "targets": targets or [],
        }, sequence={}, organism={}, instrument={}, p1_annotation={},
    )


def one_group(sequence="AG", linkage="NORMAL_PHOSPHATE", **kwargs):
    config = cfg(**kwargs); result = generate_dinucleotide_candidates(sequence, ".", config=config)
    return deepcopy(next(row for row in result.candidates if row["Linkage_State"] == linkage and row["Charge"] == kwargs.get("charges", (1,))[0])), config


def profile_feature_peaks(group, scans=3, start=1.0, step=0.02, points_per_scan=1):
    target = float(group["Theoretical_mz"]); rows = []
    for scan in range(scans):
        for point in range(points_per_scan):
            rows.append(Peak(target + (point-(points_per_scan-1)/2)*target*0.2e-6, 1000+scan*100+point*10, start+scan*step, f"scan{scan}"))
    return rows


def test_binary_raw_matching_arbitrary_theoretical_mz():
    group, config = one_group(); target = group["Theoretical_mz"]
    raw = [{"mz": target, "intensity": 10, "rt": 1.0, "scan_id": "a"}, {"mz": target*1.0001, "intensity": 10, "rt": 1.0, "scan_id": "a"}]
    rows, matches = match_raw_groups([group], raw, 10)
    assert rows[0]["Raw_Profile_Point_Count"] == 1
    assert len(matches[group["Dinucleotide_Group_ID"]]) == 1


@pytest.mark.parametrize("charge,polarity", [(1, "positive"), (2, "positive"), (1, "negative"), (2, "negative")])
def test_raw_matching_charge_and_polarity_are_candidate_driven(charge, polarity):
    group, config = one_group(charges=(charge,), polarity=polarity)
    audit = audit_dinucleotide_features([group], profile_feature_peaks(group, scans=2), config)
    assert audit.raw_group_rows[0]["Raw_Profile_Point_Count"] == 2


def test_search_range_outside_skips_execution():
    group, config = one_group(search=(100, 200))
    rows, matches = match_raw_groups([group], [{"mz": group["Theoretical_mz"], "intensity": 1, "rt": 1, "scan_id": "s"}], 10)
    assert rows[0]["Search_Executed"] is False and not matches[group["Dinucleotide_Group_ID"]]


def test_same_scan_profile_points_form_one_spectrum_peak():
    group, config = one_group(tolerance=20)
    audit = audit_dinucleotide_features([group], profile_feature_peaks(group, scans=1, points_per_scan=3), config)
    assert len(audit.spectrum_peaks) == 1
    assert audit.spectrum_peaks[0]["Local_Profile_Point_Count"] == 3
    assert audit.spectrum_peaks[0]["Spectrum_ID"] == "scan0"


def test_multiple_scans_form_one_chromatographic_feature():
    group, config = one_group()
    audit = audit_dinucleotide_features([group], profile_feature_peaks(group, scans=4), config)
    assert len(audit.features) == 1
    assert audit.features[0]["Unique_Spectrum_Count"] == 4
    assert audit.features[0]["RT_Start"] < audit.features[0]["RT_End"]


def test_rt_gap_splits_and_retains_multiple_features():
    group, config = one_group()
    peaks = profile_feature_peaks(group, scans=2, start=1.0) + profile_feature_peaks(group, scans=2, start=2.0)
    audit = audit_dinucleotide_features([group], peaks, config)
    assert len(audit.features) == 2
    assert {row["RT_Apex"] < 1.5 for row in audit.features} == {True, False}


@pytest.mark.parametrize("error,expected", [(1.9, "WITHIN_STRONG_TOLERANCE"), (4.9, "WITHIN_MODERATE_TOLERANCE"), (9.9, "WITHIN_SEARCH_TOLERANCE"), (10.1, "OUTSIDE_SEARCH_TOLERANCE"), (None, "NOT_EVALUABLE")])
def test_generic_mass_accuracy_thresholds(error, expected):
    assert classify_mass_accuracy(error, 2, 5, 10) == expected


def test_apex_and_centroid_mass_accuracy_are_separate():
    group, config = one_group(tolerance=20); target = group["Theoretical_mz"]
    peaks = [Peak(target*(1+1e-6), 1000, 1.0, "a"), Peak(target*(1+7e-6), 100, 1.02, "b")]
    audit = audit_dinucleotide_features([group], peaks, config)
    feature = audit.features[0]
    assert feature["Mass_Accuracy_Class_Apex"] == "WITHIN_STRONG_TOLERANCE"
    assert feature["Mass_Accuracy_Class_Centroid"] in {"WITHIN_STRONG_TOLERANCE", "WITHIN_MODERATE_TOLERANCE"}


def test_single_point_background_classification():
    group, config = one_group(); audit = audit_dinucleotide_features([group], profile_feature_peaks(group, scans=1), config)
    assert audit.features[0]["Background_Status"] == "SINGLE_POINT_EVENT"


def test_persistent_background_trace_is_rejected_by_generic_interpretation():
    group, config = one_group(); peaks = profile_feature_peaks(group, scans=25, start=1.0, step=0.05)
    audit = audit_dinucleotide_features([group], peaks, config)
    classify_features(audit, [group], config); interpret_groups([group], audit)
    assert audit.features[0]["Background_Status"] == "PERSISTENT_BACKGROUND_TRACE"
    assert audit.features[0]["Feature_Quality_Status"] == "BACKGROUND_TRACE_REJECTED"


def test_same_scan_CHNOSP_isotope_envelope_and_sulfur_zero():
    group, config = one_group(); target = float(group["Theoretical_mz"]); composition = {"C": 20, "H": 25, "N": 10, "O": 11, "P": 1}
    expected = _expected_isotope_abundance(composition); peaks = []
    for scan in range(2):
        rt = 1+scan*.02; mono = 10000
        peaks.extend([Peak(target, mono, rt, f"s{scan}"), Peak(target+1.00335483507, mono*expected["M+1"], rt, f"s{scan}"), Peak(target+2*1.00335483507, mono*expected["M+2"], rt, f"s{scan}")])
    audit = audit_dinucleotide_features([group], peaks, config)
    isotope = audit.isotopes[0]
    assert isotope["Apex_Spectrum_ID"] in {"s0", "s1"}
    assert isotope["Envelope_Status"] == "ENVELOPE_COMPATIBLE"
    assert isotope["Sulfur_Count"] == 0 and isotope["Sulfur_Envelope_Status"] == "NOT_APPLICABLE"


def test_sulfur_one_is_assessed_from_complete_composition():
    group, config = one_group(linkage="PHOSPHOROTHIOATE")
    audit = audit_dinucleotide_features([group], profile_feature_peaks(group, scans=2), config)
    assert audit.isotopes[0]["Sulfur_Count"] == 1
    assert audit.isotopes[0]["Sulfur_Isotope_Contribution_Assessed"] is True


def test_sulfur_more_than_one_is_generic():
    group, config = one_group(linkage="PHOSPHOROTHIOATE"); group["Elemental_Composition"] = group["Elemental_Composition"].replace("S1", "S2"); group["Final_Elemental_Composition"] = group["Elemental_Composition"]
    audit = audit_dinucleotide_features([group], profile_feature_peaks(group, scans=2), config)
    assert audit.isotopes[0]["Sulfur_Count"] == 2


def test_isotope_requires_same_scan_not_cross_scan():
    group, config = one_group(); target = group["Theoretical_mz"]
    peaks = [Peak(target, 1000, 1.0, "mono"), Peak(target+1.00335483507, 100, 1.01, "plus1"), Peak(target+2.00670967014, 50, 1.02, "plus2"), Peak(target, 900, 1.02, "mono2")]
    audit = audit_dinucleotide_features([group], peaks, config)
    isotope = audit.isotopes[0]
    assert isotope["Envelope_Point_Count"] == 1
    assert isotope["Envelope_Status"] == "ENVELOPE_INCOMPATIBLE"


def test_isotope_peak_overlap_marks_confounded():
    group, config = one_group(); other = deepcopy(group); other["Dinucleotide_Group_ID"] = other["Chemical_State_ID"] = "OTHER"; other["Theoretical_mz"] = group["Theoretical_mz"] + 1.00335483507
    feature = {"Dinucleotide_Feature_ID": "F", "Dinucleotide_Group_ID": group["Dinucleotide_Group_ID"]}
    legacy = [{"Feature_ID": "F", "Candidate_ID": group["Dinucleotide_Group_ID"], "Envelope_Status": "ENVELOPE_COMPATIBLE", "Envelope_Assessed": True, "Sulfur_Count": 0}]
    mapped = _map_isotopes(legacy, [feature], {group["Dinucleotide_Group_ID"]: group, "OTHER": other}, {**__import__('rna_masshunter.p1_sap_dinucleotide_candidates', fromlist=['dinucleotide_settings']).dinucleotide_settings(config)})
    assert mapped[0]["Isotope_Peak_Shared"] is True
    assert mapped[0]["Envelope_Status"] == "ENVELOPE_CONFOUNDED"


def test_shared_physical_feature_competition_distinguishes_composition_and_linkage():
    normal, _config = one_group(); pt = deepcopy(normal); pt["Dinucleotide_Group_ID"] = pt["Chemical_State_ID"] = "PT"; pt["Linkage_State"] = "PHOSPHOROTHIOATE"; pt["Final_Elemental_Composition"] = "C1S1"
    features = [{"Dinucleotide_Feature_ID": "F1", "Dinucleotide_Group_ID": normal["Dinucleotide_Group_ID"], "Physical_Feature_ID": "PF"}, {"Dinucleotide_Feature_ID": "F2", "Dinucleotide_Group_ID": "PT", "Physical_Feature_ID": "PF"}]
    rows = _competition_rows(features, {normal["Dinucleotide_Group_ID"]: normal, "PT": pt})
    assert all(row["Competing_Dinucleotide_Group_Count"] == 1 for row in rows)
    assert all("NORMAL_PHOSPHATE_VS_PT_COMPETITION" in row["Competition_Types"] for row in rows)
    assert all(row["Candidate_Specific"] is False and row["Linkage_Specific"] is False for row in rows)


def test_same_composition_structural_isomers_are_not_structure_specific():
    group, _config = one_group(); group["Structural_Assignment_Count"] = 3
    features = [{"Dinucleotide_Feature_ID": "F", "Dinucleotide_Group_ID": group["Dinucleotide_Group_ID"], "Physical_Feature_ID": "PF"}]
    row = _competition_rows(features, {group["Dinucleotide_Group_ID"]: group})[0]
    assert row["Structure_Specific"] is False


def test_generic_ms2_provenance_keeps_low_product_mz(monkeypatch):
    import rna_masshunter.p1_sap_dinucleotide_feature_audit as module
    group, _config = one_group(); feature = {"Dinucleotide_Feature_ID": "F", "Dinucleotide_Group_ID": group["Dinucleotide_Group_ID"], "RT_Start": .9, "RT_End": 1.1}
    spectrum = {"ms level": 2, "id": "ms2", "m/z array": [50.0, 499.9, 500.0, 700.0], "scanList": {"scan": [{"scan start time": 1.0}]}, "precursorList": {"precursor": [{"selectedIonList": {"selectedIon": [{"selected ion m/z": group["Theoretical_mz"]}]}, "isolationWindow": {"isolation window target m/z": group["Theoretical_mz"], "isolation window lower offset": 1.0, "isolation window upper offset": 1.0}}]}}
    monkeypatch.setattr(module, "iter_spectra", lambda _path: iter([spectrum]))
    row = build_ms2_provenance("unused", [feature], {group["Dinucleotide_Group_ID"]: group}, 10)[0]
    assert row["MS2_Product_Min_mz"] == "50.0"
    assert row["MS2_Product_Count_Below_500"] == 2 and row["MS2_Product_Count_At_Or_Above_500"] == 2
    assert row["MS2_Model_Applicable"] is False


def test_generic_feature_classification_qualified_but_structurally_ambiguous():
    group, config = one_group(); group["Structural_Assignment_Count"] = 2
    audit = audit_dinucleotide_features([group], profile_feature_peaks(group, scans=3), config)
    classify_features(audit, [group], config); interpret_groups([group], audit)
    assert audit.features[0]["Feature_Quality_Status"] in {"QUALIFIED_BUT_STRUCTURALLY_AMBIGUOUS", "ISOTOPE_INCOMPATIBLE"}
    assert group["Source_Bond_Resolution_Status"] == "SOURCE_BOND_UNRESOLVED"
    assert group["Sequence_Position_Localized"] is False


def test_targets_do_not_change_generic_groups_or_features():
    base = build_p1_sap_dinucleotide_audit(".", "AG", [], cfg(targets=[]), audit_level="full")
    mz = base["generated"].candidates[0]["Theoretical_mz"]
    targeted = build_p1_sap_dinucleotide_audit(".", "AG", [], cfg(targets=[{"label": "arbitrary", "theoretical_mz": mz, "tolerance_ppm": 1}]), audit_level="full")
    normalize = lambda rows: [{key: value for key, value in row.items() if not key.endswith("Runtime") and "RSS" not in key and "Tracemalloc" not in key} for row in rows]
    assert normalize(base["generated"].candidates) == normalize(targeted["generated"].candidates)
    assert base["features"] == targeted["features"]
    assert targeted["sheets"]["P1_SAP_Dinuc_Targets"][0]["Target_Label"] == "arbitrary"
