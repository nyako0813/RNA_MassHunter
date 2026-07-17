from dataclasses import replace

import numpy as np
import pytest

from rna_masshunter.sciex_intact_peak_detection import (
    FORMAL_FALSE,
    SciexIntactPeakDetectionParameters,
    _centroid_and_areas,
    detect_sciex_intact_peaks,
)


BASE_PARAMETERS = SciexIntactPeakDetectionParameters(
    baseline_window_da=30.0,
    noise_window_da=10.0,
    prominence_window_da=40.0,
    absolute_height_floor=0.1,
    absolute_prominence_floor=0.1,
)


def run(masses, intensities, parameters=BASE_PARAMETERS, **metadata):
    return detect_sciex_intact_peaks(
        masses,
        intensities,
        profile_type=metadata.get("profile_type", "NEUTRAL_MASS_PROFILE"),
        input_status=metadata.get("input_status", "SUPPORTED_INPUT"),
        eligible_for_neutral_mass_analysis=metadata.get("eligible", True),
        parameters=parameters,
    )


def axis(step=0.1, stop=100.0):
    return np.arange(0.0, stop + step / 2.0, step)


def gaussian(x, center, amplitude=20.0, sigma=1.0):
    return amplitude * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def test_ineligible_mz_profile_is_skipped_without_conversion():
    result = run([500, 501, 502, 503, 504], [1, 2, 3, 2, 1], profile_type="MZ_PROFILE", eligible=False)
    assert result.diagnostics["Detection_Status"] == "SKIPPED_INELIGIBLE_PROFILE"
    assert result.peaks == ()


def test_unsupported_input_is_skipped():
    result = run([1, 2, 3, 4, 5], [1, 2, 3, 2, 1], input_status="UNSUPPORTED_PROFILE_TYPE", eligible=False)
    assert result.diagnostics["Detection_Status"] == "SKIPPED_INELIGIBLE_PROFILE"


def test_length_mismatch_returns_invalid_axis():
    result = run([1, 2, 3, 4, 5], [1, 2, 3])
    assert result.diagnostics["Detection_Status"] == "INVALID_AXIS"


def test_nonfinite_mass_returns_invalid_axis():
    result = run([1, 2, np.nan, 4, 5], [1, 2, 3, 2, 1])
    assert result.diagnostics["Detection_Status"] == "INVALID_AXIS"
    assert result.diagnostics["Nonfinite_Value_Count"] == 1


def test_nonfinite_intensity_returns_invalid_intensity():
    result = run([1, 2, 3, 4, 5], [1, 2, np.inf, 2, 1])
    assert result.diagnostics["Detection_Status"] == "INVALID_INTENSITY"


def test_negative_intensity_is_invalid_and_not_clamped():
    result = run([1, 2, 3, 4, 5], [1, 2, -1, 2, 1])
    assert result.diagnostics["Detection_Status"] == "INVALID_INTENSITY"
    assert result.diagnostics["Negative_Intensity_Count"] == 1


def test_nonmonotonic_mass_returns_invalid_axis():
    result = run([1, 2, 4, 3, 5], [1, 2, 3, 2, 1])
    assert result.diagnostics["Detection_Status"] == "INVALID_AXIS"


def test_duplicate_mass_is_reported():
    result = run([1, 2, 2, 3, 4], [1, 2, 3, 2, 1])
    assert result.diagnostics["Detection_Status"] == "INVALID_AXIS"
    assert result.diagnostics["Duplicate_Mass_Count"] == 1


def test_insufficient_points_returns_diagnostic_result():
    result = run([1, 2, 3, 4], [1, 2, 1, 0])
    assert result.diagnostics["Detection_Status"] == "INSUFFICIENT_POINTS"


def test_flat_baseline_has_no_peaks():
    x = axis()
    result = run(x, np.full_like(x, 7.0))
    assert result.diagnostics["Detected_Sensitive_Peak_Count"] == 0


def test_single_gaussian_peak_has_expected_apex_and_width():
    x = axis()
    result = run(x, 5 + gaussian(x, 50, sigma=1.2))
    assert len(result.peaks) == 1
    peak = result.peaks[0]
    assert peak["Apex_Mass"] == pytest.approx(50.0, abs=0.1)
    assert peak["FWHM_Da"] == pytest.approx(2.355 * 1.2, rel=0.08)


def test_two_well_separated_peaks_are_detected():
    x = axis()
    y = 4 + gaussian(x, 30, sigma=0.8) + gaussian(x, 70, amplitude=15, sigma=1.0)
    result = run(x, y)
    assert [peak["Apex_Mass"] for peak in result.peaks] == pytest.approx([30, 70], abs=0.1)


def test_nearby_peaks_survive_one_da_minimum_distance():
    x = axis()
    y = 3 + gaussian(x, 49, sigma=0.45) + gaussian(x, 51, sigma=0.45)
    result = run(x, y)
    assert len(result.peaks) == 2
    assert result.diagnostics["Suppressed_By_Distance_Count"] == 0


def test_peaks_less_than_five_da_apart_are_not_overmerged():
    x = axis()
    y = 3 + gaussian(x, 48.5, sigma=0.45) + gaussian(x, 51.5, amplitude=18, sigma=0.45)
    result = run(x, y)
    assert len(result.peaks) == 2
    assert result.peaks[1]["Apex_Mass"] - result.peaks[0]["Apex_Mass"] == pytest.approx(3.0, abs=0.2)


def test_mass_dependent_baseline_is_removed_without_losing_peak():
    x = axis()
    y = 5 + 0.05 * x + gaussian(x, 55, sigma=1.0)
    result = run(x, y)
    assert any(peak["Apex_Mass"] == pytest.approx(55, abs=0.1) for peak in result.peaks)
    assert result.diagnostics["Baseline_Max"] > result.diagnostics["Baseline_Min"]


def test_signed_corrected_signal_can_contain_negative_values():
    x = axis()
    y = 5 + gaussian(x, 50, sigma=1.0)
    y[100] = 4.0
    result = run(x, y)
    assert min(result.signed_baseline_corrected_intensity) < 0
    assert min(result.nonnegative_corrected_quantification_weights) == 0
    assert result.diagnostics["Detection_Status"].startswith("DETECTION_COMPLETED")


def test_local_noise_profile_tracks_different_noise_regions():
    rng = np.random.default_rng(7)
    x = axis(stop=200)
    y = 10 + np.r_[rng.normal(0, 0.1, len(x) // 2), rng.normal(0, 1.0, len(x) - len(x) // 2)]
    y -= min(y) - 1
    result = run(x, y, replace(BASE_PARAMETERS, absolute_height_floor=10, absolute_prominence_floor=10))
    assert result.diagnostics["Estimated_Noise_Local_Max"] > 4 * result.diagnostics["Estimated_Noise_Local_Min"]


def test_exact_plateau_uses_center_and_reports_extent():
    x = axis()
    y = np.ones_like(x)
    y[490:511] = 20
    result = run(x, y)
    peak = result.peaks[0]
    assert peak["Plateau_Size_Points"] == 21
    assert peak["Apex_Index"] == 500
    assert peak["Plateau_Center_Mass"] == pytest.approx(50.0)


def test_even_length_plateau_center_mass_is_midpoint():
    x = axis()
    y = np.ones_like(x)
    y[490:510] = 20
    result = run(x, y)
    peak = result.peaks[0]
    assert peak["Plateau_Size_Points"] == 20
    assert peak["Apex_Index"] == 499
    assert peak["Plateau_Center_Mass"] == pytest.approx((x[490] + x[509]) / 2)


def test_shallow_valley_is_diagnostic_not_identity_assignment():
    x = axis()
    y = 2 + gaussian(x, 48, sigma=1.3) + gaussian(x, 52, amplitude=18, sigma=1.3)
    result = run(x, y, replace(BASE_PARAMETERS, minimum_width_da=0.5))
    assert len(result.peaks) == 2
    assert any(peak["Shallow_Valley_Neighbor_Flag"] for peak in result.peaks)
    assert all(peak["Molecular_Identity_Assigned"] is False for peak in result.peaks)


def test_broad_peak_is_flagged_but_not_removed():
    x = axis(stop=160)
    y = 3 + sum(gaussian(x, center, sigma=0.8) for center in [20, 45, 70, 95]) + gaussian(x, 130, sigma=8)
    params = replace(BASE_PARAMETERS, broad_peak_quantile=0.8, severe_broad_peak_width_da=50)
    result = run(x, y, params)
    broad = [peak for peak in result.peaks if peak["Broad_Peak_Flag"]]
    assert broad
    assert any(abs(peak["Apex_Mass"] - 130) < 0.2 for peak in broad)


def test_severe_broad_peak_is_flagged_and_retained():
    x = axis(stop=160)
    y = 3 + gaussian(x, 80, sigma=10)
    result = run(x, y, replace(BASE_PARAMETERS, baseline_window_da=100, prominence_window_da=150, severe_broad_peak_width_da=20))
    assert len(result.peaks) == 1
    assert result.peaks[0]["Severe_Broad_Peak_Flag"] is True


def test_edge_peak_marks_incomplete_quantification():
    x = axis()
    y = 2 + gaussian(x, 0.4, sigma=1.2)
    result = run(x, y, replace(BASE_PARAMETERS, height_sigma_multiplier=0, prominence_sigma_multiplier=0, strict_prominence_sigma_multiplier=0, minimum_width_da=0.2))
    assert result.peaks
    peak = result.peaks[0]
    assert peak["Edge_Peak_Flag"] is True
    assert peak["Peak_Area_Complete"] is False
    assert peak["Centroid_Complete"] is False


def test_shared_valley_boundaries_are_identical_and_areas_do_not_overlap():
    x = axis()
    y = 3 + gaussian(x, 45, sigma=1.2) + gaussian(x, 55, amplitude=18, sigma=1.2)
    result = run(x, y)
    left, right = result.peaks
    assert left["Right_Boundary_Index"] == right["Left_Boundary_Index"]
    assert left["Right_Boundary_Mass"] == right["Left_Boundary_Mass"]


def test_sensitive_boundaries_are_not_recomputed_for_strict_subset():
    x = axis()
    y = 3 + gaussian(x, 45, amplitude=20, sigma=1) + gaussian(x, 55, amplitude=4, sigma=1)
    params = replace(
        BASE_PARAMETERS,
        absolute_prominence_floor=0.1,
        prominence_sigma_multiplier=1,
        strict_prominence_sigma_multiplier=1,
    )
    result = run(x, y, params)
    assert len(result.peaks) == 2
    assert all(peak["Boundary_Peak_Set_Tier"] == "SENSITIVE" for peak in result.peaks)
    assert all(peak["Boundary_Recomputed_After_Filtering"] is False for peak in result.peaks)


def test_raw_and_corrected_areas_are_separate():
    x = axis()
    result = run(x, 10 + gaussian(x, 50, sigma=1))
    peak = result.peaks[0]
    assert peak["Peak_Area_Raw"] > peak["Peak_Area_Baseline_Corrected"] > 0


def test_area_uses_raw_not_smoothed_signal():
    x = axis()
    raw = 7 + gaussian(x, 50, sigma=0.7)
    result = run(x, raw)
    peak = result.peaks[0]
    left, right = peak["Left_Boundary_Index"], peak["Right_Boundary_Index"]
    expected = np.trapezoid(raw[left:right + 1], x[left:right + 1])
    assert peak["Peak_Area_Raw"] == pytest.approx(expected)


def test_centroid_zero_weight_falls_back_to_apex():
    masses = np.array([1.0, 2.0, 4.0])
    centroid, raw_area, corrected_area, fallback = _centroid_and_areas(
        masses, np.array([5.0, 5.0, 5.0]), np.zeros(3), 0, 2, 2.0,
    )
    assert centroid == 2.0 and corrected_area == 0 and raw_area > 0 and fallback is True


def test_nonuniform_axis_uses_actual_mass_for_area():
    x = np.r_[np.arange(0, 45, 0.2), np.arange(45, 55, 0.1), np.arange(55, 100.1, 0.3)]
    y = 5 + gaussian(x, 50, sigma=1)
    result = run(x, y)
    peak = max(result.peaks, key=lambda item: item["Apex_Intensity_Raw"])
    left, right = peak["Left_Boundary_Index"], peak["Right_Boundary_Index"]
    assert result.diagnostics["Mass_Axis_Uniform"] is False
    assert peak["Peak_Area_Raw"] == pytest.approx(np.trapezoid(y[left:right + 1], x[left:right + 1]))


def test_parameter_provenance_distinguishes_explicit_and_derived():
    x = axis()
    result = run(x, 5 + gaussian(x, 50))
    rows = result.provenance_rows()
    assert any(row["Parameter_Name"] == "baseline_quantile" and row["Parameter_Source"] == "explicit" for row in rows)
    assert any(row["Parameter_Name"] == "Baseline_Window_Points" and row["Parameter_Source"] == "derived" for row in rows)


def test_formal_nonpropagation_flags_are_always_false():
    x = axis()
    result = run(x, 5 + gaussian(x, 50))
    for name, expected in FORMAL_FALSE.items():
        assert result.diagnostics[name] is expected
        assert all(peak[name] is expected for peak in result.peaks)


def test_result_is_deterministic():
    rng = np.random.default_rng(42)
    x = axis()
    y = 5 + gaussian(x, 50) + rng.normal(0, 0.05, len(x))
    first = run(x, y)
    second = run(x, y)
    assert first.diagnostics_row() == second.diagnostics_row()
    assert first.peak_rows() == second.peak_rows()
    assert first.parameter_provenance == second.parameter_provenance


def test_input_arrays_are_not_modified():
    x = axis()
    y = 5 + gaussian(x, 50)
    x_before = x.copy()
    y_before = y.copy()
    run(x, y)
    np.testing.assert_array_equal(x, x_before)
    np.testing.assert_array_equal(y, y_before)


@pytest.mark.parametrize(
    "parameters",
    [
        replace(BASE_PARAMETERS, baseline_quantile=1.1),
        replace(BASE_PARAMETERS, smoothing_window_points=4),
        replace(BASE_PARAMETERS, smoothing_polyorder=5),
    ],
)
def test_invalid_parameters_are_rejected(parameters):
    with pytest.raises(ValueError):
        run(axis(), np.ones_like(axis()), parameters)
