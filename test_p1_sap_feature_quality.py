from types import SimpleNamespace

import pytest

from rna_masshunter.audit_policy import AuditPolicy, included_sheet_names
from rna_masshunter.p1_sap_feature_quality import (
    FORMAL_FALSE,
    ISOTOPE_SPACING,
    P1_SAP_FEATURE_QUALITY_COLUMNS,
    P1_SAP_ISOTOPE_AUDIT_COLUMNS,
    P1_SAP_QUALITY_SUMMARY_COLUMNS,
    P1_SAP_REFINED_FEATURE_COLUMNS,
    P1_SAP_SPECTRUM_PEAK_COLUMNS,
    _expected_isotope_abundance,
    _parse_composition,
    assess_isotope_envelope,
    build_p1_sap_feature_quality,
    build_spectrum_level_peaks,
)


def cfg(**overrides):
    return SimpleNamespace(p1_annotation={"mz_tolerance_ppm": 20.0, **overrides})


def candidate(cid="c1", mz=100.0, charge=1, composition=None, family="PHOSPHOROTHIOATE"):
    return {
        "Chemical_State_ID": cid,
        "Chemical_Family": family,
        "Charge": charge,
        "Theoretical_mz": mz,
        "Elemental_Composition": composition,
    }


def point(index, mz, intensity, rt, scan):
    return {"index": index, "mz": mz, "intensity": intensity, "rt": rt, "scan_id": scan}


def feature(
    fid="F1", cid="c1", point_ids=None, start=0.4, end=0.6, apex=0.5,
    charge=1, physical=None, family="PHOSPHOROTHIOATE", mz=100.0,
):
    point_ids = list(point_ids or [])
    return {
        "Feature_ID": fid,
        "Physical_Feature_ID": physical or f"PF_{fid}",
        "Chemical_State_ID": cid,
        "Chemical_Family": family,
        "Charge": charge,
        "RT_Start": start,
        "RT_End": end,
        "RT_Apex": apex,
        "RT_Span": end - start,
        "Spectrum_Count": 999,
        "Profile_Point_Count": 999,
        "Apex_mz": mz,
        "mz_Centroid": mz,
        "Mass_Error_ppm_at_Apex": 0.0,
        "Mass_Error_ppm_at_Centroid": 0.0,
        "mz_SD": 0.0,
        "Apex_Intensity": 100.0,
        "Integrated_Intensity": 140.0,
        "_point_ids": point_ids,
    }


def localized_points(mz=100.0):
    return [
        point(0, mz, 20.0, 0.4, "s1"),
        point(1, mz, 100.0, 0.5, "s2"),
        point(2, mz, 20.0, 0.6, "s3"),
    ]


def run(candidates, features, raw, **overrides):
    return build_p1_sap_feature_quality(candidates, features, raw, cfg(**overrides), 20.0)


def test_same_scan_profile_points_form_one_local_peak_with_centroid_and_integrated_intensity():
    raw = [point(0, 100.0, 10.0, 0.5, "s1"), point(1, 100.001, 30.0, 0.5, "s1")]
    rows = build_spectrum_level_peaks([feature(point_ids=[0, 1])], raw, 20.0)
    assert len(rows) == 1
    assert rows[0]["Local_Profile_Point_Count"] == 2
    assert rows[0]["Local_Integrated_Intensity"] == 40.0
    assert rows[0]["Local_Centroid_mz"] == pytest.approx(100.00075)


def test_spectrum_count_is_unique_scan_count_not_legacy_metadata():
    raw = localized_points()
    result = run([candidate()], [feature(point_ids=[0, 1, 2])], raw)
    assert result["quality_rows"][0]["Spectrum_Count"] == 3
    assert result["quality_rows"][0]["Profile_Point_Count"] == 3


def test_local_peak_mz_tolerance_splits_points():
    raw = [point(0, 100.0, 10.0, 0.5, "s1"), point(1, 100.01, 20.0, 0.5, "s1")]
    rows = build_spectrum_level_peaks([feature(point_ids=[0, 1])], raw, 20.0)
    assert len(rows) == 2


def test_refined_rt_gap_splits_at_configured_maximum():
    raw = [
        point(0, 100.0, 20.0, 0.0, "s1"),
        point(1, 100.0, 100.0, 0.05, "s2"),
        point(2, 100.0, 20.0, 0.20, "s3"),
    ]
    f = feature(point_ids=[0, 1, 2], start=0.0, end=0.2, apex=0.05)
    result = run([candidate()], [f], raw, max_refined_rt_gap_min=0.08)
    assert len(result["refined_features"]) == 2


def test_refined_features_are_split_by_charge():
    raw = localized_points()
    candidates = [candidate("z1", charge=1), candidate("z2", charge=2)]
    features = [
        feature("Fz1", "z1", [0, 1, 2], charge=1, physical="PF"),
        feature("Fz2", "z2", [0, 1, 2], charge=2, physical="PF"),
    ]
    result = run(candidates, features, raw)
    assert {row["Charge"] for row in result["refined_features"]} == {1, 2}


def test_background_uses_only_candidate_mz_neighborhood():
    raw = [
        point(0, 100.0, 2.0, 0.2, "b1"),
        point(1, 200.0, 10000.0, 0.2, "b1"),
        *[dict(row) for row in localized_points()],
        point(5, 100.0, 2.0, 0.8, "b2"),
        point(6, 200.0, 10000.0, 0.8, "b2"),
    ]
    for index, row in enumerate(raw):
        row["index"] = index
    result = run([candidate()], [feature(point_ids=[2, 3, 4])], raw)
    quality = result["quality_rows"][0]
    assert quality["Baseline_Intensity"] == 2.0
    assert quality["Apex_Local_Contrast"] == 50.0


@pytest.mark.parametrize("charge", [1, 2, 3])
def test_isotope_spacing_is_charge_adjusted(charge):
    composition, _ = _parse_composition("C10H10S")
    expected = _expected_isotope_abundance(composition)
    mz = 300.0
    spacing = ISOTOPE_SPACING / charge
    raw = [
        point(0, mz, 1000.0, 0.5, "s1"),
        point(1, mz + spacing, 1000.0 * expected["M+1"], 0.5, "s1"),
        point(2, mz + 2 * spacing, 1000.0 * expected["M+2"], 0.5, "s1"),
    ]
    isotope = assess_isotope_envelope(feature(mz=mz, point_ids=[0], charge=charge), raw, candidate(mz=mz, charge=charge, composition="C10H10S"), 20.0)
    assert isotope["Envelope_Point_Count"] == 3
    assert isotope["Observed_Mz_Mplus1"] == pytest.approx(mz + spacing)


def test_m_mplus1_mplus2_are_taken_from_one_scan():
    composition, _ = _parse_composition("C10H10S")
    expected = _expected_isotope_abundance(composition)
    raw = [
        point(0, 100.0, 1000.0, 0.5, "apex"),
        point(1, 100.0 + ISOTOPE_SPACING, expected["M+1"] * 1000, 0.5, "apex"),
        point(2, 100.0 + 2 * ISOTOPE_SPACING, expected["M+2"] * 1000, 0.5, "apex"),
        point(3, 100.0, 500.0, 0.55, "other"),
    ]
    isotope = assess_isotope_envelope(feature(point_ids=[0]), raw, candidate(composition="C10H10S"), 20.0)
    assert isotope["Envelope_Spectrum_ID"] == "apex"
    assert isotope["Envelope_Status"] == "ENVELOPE_COMPATIBLE"


def test_isotope_peaks_from_different_scans_are_not_one_envelope():
    raw = [
        point(0, 100.0, 1000.0, 0.45, "s1"),
        point(1, 100.0 + ISOTOPE_SPACING, 100.0, 0.50, "s2"),
        point(2, 100.0 + 2 * ISOTOPE_SPACING, 50.0, 0.55, "s3"),
    ]
    isotope = assess_isotope_envelope(feature(point_ids=[0]), raw, candidate(composition="C10H10S"), 20.0)
    assert isotope["Envelope_Spectrum_ID"] == "s1"
    assert isotope["Envelope_Point_Count"] == 1
    assert isotope["Envelope_Status"] == "ENVELOPE_INCOMPATIBLE"


@pytest.mark.parametrize("composition,status", [(None, "NOT_ASSESSED"), ("", "NOT_ASSESSED"), ("MODEL_NOT_DEFINED", "NOT_ASSESSED")])
def test_undefined_composition_is_not_assessed(composition, status):
    isotope = assess_isotope_envelope(feature(), localized_points(), candidate(composition=composition), 20.0)
    assert isotope["Envelope_Assessed"] is False
    assert isotope["Envelope_Status"] == status
    assert isotope["Sulfur_Envelope_Compatible"] == "NOT_ASSESSED"


@pytest.mark.parametrize("composition", ["O-1", "C10H10O-1", "S-1"])
def test_delta_composition_is_not_used_as_positive_composition(composition):
    isotope = assess_isotope_envelope(feature(), localized_points(), candidate(composition=composition), 20.0)
    assert isotope["Envelope_Assessed"] is False
    assert isotope["Envelope_Status"] == "MODEL_NOT_APPLICABLE"


def test_nominal_distribution_includes_cross_terms_and_sulfur_mplus2():
    carbon, _ = _parse_composition("C2")
    carbon_one, _ = _parse_composition("C1")
    no_sulfur, _ = _parse_composition("C10H10")
    sulfur, _ = _parse_composition("C10H10S")
    c2 = _expected_isotope_abundance(carbon)
    c1 = _expected_isotope_abundance(carbon_one)
    assert c2["M+2"] > 0.0
    assert c2["M+2"] == pytest.approx(c1["M+1"] ** 2)
    assert _expected_isotope_abundance(sulfur)["M+2"] > _expected_isotope_abundance(no_sulfur)["M+2"]


def test_shared_isotope_peak_and_isomer_are_explicit():
    composition, _ = _parse_composition("C10H10S")
    expected = _expected_isotope_abundance(composition)
    raw = []
    for scan_number, (rt, scale) in enumerate(((0.4, 0.2), (0.5, 1.0), (0.6, 0.2))):
        base = len(raw)
        raw.extend([
            point(base, 100.0, 1000.0 * scale, rt, f"s{scan_number}"),
            point(base + 1, 100.0 + ISOTOPE_SPACING, expected["M+1"] * 1000 * scale, rt, f"s{scan_number}"),
            point(base + 2, 100.0 + 2 * ISOTOPE_SPACING, expected["M+2"] * 1000 * scale, rt, f"s{scan_number}"),
        ])
    ids = [0, 3, 6]
    candidates = [candidate("a", composition="C10H10S"), candidate("b", composition="C10H10S")]
    features = [
        feature("Fa", "a", ids, physical="PF"),
        feature("Fb", "b", ids, physical="PF"),
    ]
    result = run(candidates, features, raw)
    assert all(row["Isotope_Peak_Shared_With_Other_Candidate"] for row in result["isotope_rows"])
    assert all(row["Isomer_Isotope_Indistinguishable"] for row in result["isotope_rows"])
    assert all(row["Feature_Quality_Status"] == "COMPETITION_UNRESOLVED" for row in result["quality_rows"])

def test_competition_requires_another_candidate_sharing_physical_or_local_peak():
    raw = []
    ids_100 = []
    ids_200 = []
    for scan_number, rt in enumerate((0.4, 0.5, 0.6)):
        ids_100.append(len(raw))
        raw.append(point(len(raw), 100.0, 100.0, rt, f"s{scan_number}"))
        ids_200.append(len(raw))
        raw.append(point(len(raw), 200.0, 100.0, rt, f"s{scan_number}"))
    candidates = [candidate("a"), candidate("b", mz=200.0)]
    features = [
        feature("Fa", "a", ids_100, physical="PFa"),
        feature("Fb", "b", ids_200, physical="PFb", mz=200.0),
    ]
    result = run(candidates, features, raw)
    assert result["quality_rows"][0]["Competition_Count"] == 0
    features[1]["Physical_Feature_ID"] = "PFa"
    result = run(candidates, features, raw)
    assert all(row["Competition_Count"] == 1 for row in result["quality_rows"])

def test_not_assessed_isotope_does_not_hard_exclude_but_incompatible_does():
    raw = localized_points()
    not_assessed = run([candidate(composition=None)], [feature(point_ids=[0, 1, 2])], raw)
    assert not_assessed["quality_rows"][0]["Feature_Quality_Status"] == "QUALIFIED_CHROMATOGRAPHIC_FEATURE"
    assert not_assessed["quality_rows"][0]["Isotope_Envelope_Component"] == 0.5

    incompatible = run([candidate(composition="C10H10")], [feature(point_ids=[0, 1, 2])], raw)
    assert incompatible["isotope_rows"][0]["Envelope_Status"] == "ENVELOPE_INCOMPATIBLE"
    assert incompatible["quality_rows"][0]["Feature_Quality_Status"] == "ISOTOPE_INCOMPATIBLE"
    assert incompatible["quality_rows"][0]["Isotope_Envelope_Component"] == 0.0


def test_background_thresholds_are_configurable():
    raw = localized_points()
    default = run([candidate()], [feature(point_ids=[0, 1, 2])], raw)
    assert default["quality_rows"][0]["Feature_Quality_Status"] == "QUALIFIED_CHROMATOGRAPHIC_FEATURE"
    overridden = run(
        [candidate()], [feature(point_ids=[0, 1, 2])], raw,
        background_rt_coverage_fraction=0.0,
        local_maximum_count_threshold=0,
    )
    assert overridden["quality_rows"][0]["Feature_Quality_Status"] == "BACKGROUND_TRACE_REJECTED"


def test_refined_quality_status_propagates_from_legacy_feature():
    raw = localized_points()
    result = run([candidate()], [feature(point_ids=[0, 1, 2])], raw)
    assert result["refined_features"]
    assert {row["Feature_Quality_Status"] for row in result["refined_features"]} == {"QUALIFIED_CHROMATOGRAPHIC_FEATURE"}


def test_single_point_feature_is_rejected_and_pt_count_is_zero():
    raw = [point(0, 100.0, 100.0, 0.5, "s1")]
    result = run([candidate()], [feature(point_ids=[0])], raw)
    assert result["quality_rows"][0]["Feature_Quality_Status"] == "SINGLE_POINT_REJECTED"
    assert result["summary_row"]["Qualified_PT_Feature_Count"] == 0
    assert result["summary_row"]["PT_Final_Interpretation"] == "NO_QUALIFIED_PT_FEATURE"


def test_feature_quality_sheet_is_audit_summary_only():
    names = ["P1_SAP_Feature_Quality"]
    standard, unknown_standard = included_sheet_names(names, AuditPolicy.from_level("standard"))
    audit, unknown_audit = included_sheet_names(names, AuditPolicy.from_level("audit"))
    full, unknown_full = included_sheet_names(names, AuditPolicy.from_level("full"))
    assert standard == []
    assert audit == names
    assert full == names
    assert unknown_standard == unknown_audit == unknown_full == []


def test_all_schemas_are_unique_and_formal_flags_stay_false():
    schemas = [
        P1_SAP_SPECTRUM_PEAK_COLUMNS,
        P1_SAP_REFINED_FEATURE_COLUMNS,
        P1_SAP_FEATURE_QUALITY_COLUMNS,
        P1_SAP_ISOTOPE_AUDIT_COLUMNS,
        P1_SAP_QUALITY_SUMMARY_COLUMNS,
    ]
    assert all(len(columns) == len(set(columns)) for columns in schemas)
    result = run([candidate()], [feature(point_ids=[0, 1, 2])], localized_points())
    rows = (
        result["spectrum_peaks"]
        + result["refined_features"]
        + result["quality_rows"]
        + result["isotope_rows"]
        + [result["summary_row"]]
    )
    for row in rows:
        assert all(row[key] is value for key, value in FORMAL_FALSE.items())
