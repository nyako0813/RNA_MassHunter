"""Shadow-only peak detection for SCIEX deconvoluted neutral-mass profiles."""
from __future__ import annotations

from dataclasses import dataclass, fields
from math import ceil, sqrt
from types import MappingProxyType
from typing import Any, Iterable, Mapping, Sequence
import warnings as python_warnings

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_prominences, peak_widths, savgol_filter

# numpy.trapezoid was introduced in numpy 2.0, replacing the older numpy.trapz name.
# Support both so this module works on numpy < 2.0 (where only trapz exists) and
# numpy >= 2.0 (where trapz is deprecated in favor of trapezoid).
_trapezoid = getattr(np, "trapezoid", None) or np.trapz

NEUTRAL_MASS_PROFILE = "NEUTRAL_MASS_PROFILE"
SUPPORTED_INPUT = "SUPPORTED_INPUT"
ALGORITHM_VERSION = "sciex-intact-neutral-mass-v1"
FORMAL_FALSE = {
    "SCIEX_Intact_Peak_Detection_Applied_To_Formal_Score": False,
    "SCIEX_Intact_Peak_Detection_Applied_To_Ranking": False,
    "SCIEX_Intact_Peak_Detection_Applied_To_Candidate_Filtering": False,
}


@dataclass(frozen=True)
class SciexIntactPeakDetectionParameters:
    baseline_quantile: float = 0.10
    baseline_window_da: float = 500.0
    smoothing_enabled: bool = True
    smoothing_window_points: int = 5
    smoothing_polyorder: int = 2
    noise_window_da: float = 250.0
    height_sigma_multiplier: float = 3.0
    prominence_sigma_multiplier: float = 5.0
    strict_prominence_sigma_multiplier: float = 8.0
    minimum_width_da: float = 1.0
    minimum_distance_da: float = 1.0
    broad_peak_quantile: float = 0.95
    severe_broad_peak_width_da: float = 20.0
    positive_residual_quantile: float = 0.05
    absolute_height_floor: float = 0.0
    absolute_prominence_floor: float = 0.0
    prominence_window_da: float = 1000.0
    shoulder_valley_ratio_threshold: float = 0.50
    shoulder_max_separation_da: float = 10.0
    shoulder_width_to_separation_ratio: float = 2.0 / 3.0
    minimum_points: int = 5

    def validate(self) -> None:
        if not 0 <= self.baseline_quantile <= 1:
            raise ValueError("baseline_quantile must be between 0 and 1")
        if not 0 <= self.positive_residual_quantile <= 1:
            raise ValueError("positive_residual_quantile must be between 0 and 1")
        if not 0 < self.broad_peak_quantile <= 1:
            raise ValueError("broad_peak_quantile must be in (0, 1]")
        for name in (
            "baseline_window_da", "noise_window_da", "minimum_width_da",
            "minimum_distance_da", "severe_broad_peak_width_da", "prominence_window_da",
        ):
            if float(getattr(self, name)) <= 0:
                raise ValueError(f"{name} must be positive")
        for name in (
            "height_sigma_multiplier", "prominence_sigma_multiplier",
            "strict_prominence_sigma_multiplier", "absolute_height_floor",
            "absolute_prominence_floor", "shoulder_valley_ratio_threshold",
            "shoulder_max_separation_da", "shoulder_width_to_separation_ratio",
        ):
            if float(getattr(self, name)) < 0:
                raise ValueError(f"{name} must be nonnegative")
        if self.smoothing_window_points < 1 or self.smoothing_window_points % 2 == 0:
            raise ValueError("smoothing_window_points must be a positive odd integer")
        if self.smoothing_polyorder < 0 or self.smoothing_polyorder >= self.smoothing_window_points:
            raise ValueError("smoothing_polyorder must be nonnegative and smaller than the smoothing window")
        if self.minimum_points < 3:
            raise ValueError("minimum_points must be at least 3")


@dataclass(frozen=True)
class SciexIntactParameterProvenance:
    parameter_name: str
    parameter_value: Any
    parameter_unit: str
    parameter_source: str
    parameter_reason: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "Parameter_Name": self.parameter_name,
            "Parameter_Value": self.parameter_value,
            "Parameter_Unit": self.parameter_unit,
            "Parameter_Source": self.parameter_source,
            "Parameter_Reason": self.parameter_reason,
        }


@dataclass(frozen=True)
class _FrozenRecord:
    values: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "values", MappingProxyType(dict(self.values)))

    def __getitem__(self, key: str) -> Any:
        return self.values[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self.values.get(key, default)

    def to_dict(self) -> dict[str, Any]:
        return dict(self.values)


@dataclass(frozen=True)
class SciexIntactPeakDiagnostics(_FrozenRecord):
    pass


@dataclass(frozen=True)
class SciexIntactDetectedPeak(_FrozenRecord):
    pass


@dataclass(frozen=True)
class SciexIntactPeakDetectionResult:
    parameters: SciexIntactPeakDetectionParameters
    diagnostics: SciexIntactPeakDiagnostics
    peaks: tuple[SciexIntactDetectedPeak, ...]
    parameter_provenance: tuple[SciexIntactParameterProvenance, ...]
    raw_intensity: tuple[float, ...] = ()
    baseline: tuple[float, ...] = ()
    signed_baseline_corrected_intensity: tuple[float, ...] = ()
    smoothed_detection_signal: tuple[float, ...] = ()
    nonnegative_corrected_quantification_weights: tuple[float, ...] = ()
    warnings: tuple[str, ...] = ()

    def diagnostics_row(self) -> dict[str, Any]:
        return self.diagnostics.to_dict()

    def peak_rows(self) -> list[dict[str, Any]]:
        return [peak.to_dict() for peak in self.peaks]

    def provenance_rows(self) -> list[dict[str, Any]]:
        return [item.to_dict() for item in self.parameter_provenance]


def _as_tuple(values: np.ndarray) -> tuple[float, ...]:
    return tuple(float(value) for value in values)


def _odd_window_from_da(window_da: float, step: float, point_count: int) -> int:
    value = max(1, int(round(window_da / step)) + 1)
    if value % 2 == 0:
        value += 1
    maximum = point_count if point_count % 2 else point_count - 1
    return max(1, min(value, maximum))


def _uniform_axis(steps: np.ndarray, median_step: float) -> bool:
    tolerance = max(1e-12, abs(median_step) * 1e-6)
    return bool(np.all(np.abs(steps - median_step) <= tolerance))


def _rolling_quantile_baseline(
    masses: np.ndarray,
    intensities: np.ndarray,
    quantile: float,
    window_da: float,
    uniform: bool,
    window_points: int,
) -> tuple[np.ndarray, str]:
    if uniform:
        baseline = pd.Series(intensities).rolling(
            window_points, center=True, min_periods=1,
        ).quantile(quantile).to_numpy(dtype=float)
        return baseline, "centered_truncated_point_window"
    half = window_da / 2.0
    baseline = np.empty_like(intensities)
    for index, mass in enumerate(masses):
        left = int(np.searchsorted(masses, mass - half, side="left"))
        right = int(np.searchsorted(masses, mass + half, side="right"))
        baseline[index] = float(np.quantile(intensities[left:right], quantile))
    return baseline, "centered_truncated_mass_window"


def _scaled_diff_mad(values: np.ndarray, *, nonzero_only: bool = False) -> float:
    diffs = np.diff(values)
    diffs = diffs[np.isfinite(diffs)]
    finite_count = diffs.size
    if nonzero_only:
        diffs = diffs[diffs != 0]
        if diffs.size < max(5, int(0.10 * finite_count)):
            return 0.0
    if diffs.size == 0:
        return 0.0
    center = float(np.median(diffs))
    return float(1.4826 * np.median(np.abs(diffs - center)) / sqrt(2.0))


def _local_noise_profile(
    masses: np.ndarray,
    corrected: np.ndarray,
    window_da: float,
    positive_fallback: float,
    absolute_fallback: float,
) -> tuple[np.ndarray, np.ndarray, float, list[str]]:
    global_sigma = _scaled_diff_mad(corrected)
    global_nonzero_sigma = _scaled_diff_mad(corrected, nonzero_only=True)
    median_step = float(np.median(np.diff(masses)))
    approximate_points = max(3, int(round(window_da / median_step)) + 1)
    anchor_stride = max(1, approximate_points // 10)
    anchors = np.unique(np.r_[0, np.arange(0, len(masses), anchor_stride), len(masses) - 1]).astype(int)
    anchor_values: list[float] = []
    anchor_methods: list[str] = []
    half = window_da / 2.0
    fallbacks: list[str] = []
    for index in anchors:
        left = int(np.searchsorted(masses, masses[index] - half, side="left"))
        right = int(np.searchsorted(masses, masses[index] + half, side="right"))
        local = corrected[left:right]
        sigma = _scaled_diff_mad(local)
        method = "local_diff_residual_mad"
        if not np.isfinite(sigma) or sigma <= 0:
            sigma = _scaled_diff_mad(local, nonzero_only=True)
            method = "local_nonzero_diff_residual_mad"
        if not np.isfinite(sigma) or sigma <= 0:
            sigma = global_sigma if global_sigma > 0 else global_nonzero_sigma
            method = "global_diff_residual_mad"
        if not np.isfinite(sigma) or sigma <= 0:
            sigma = positive_fallback
            method = "positive_residual_lower_quantile"
        if not np.isfinite(sigma) or sigma <= 0:
            sigma = absolute_fallback
            method = "absolute_configured_floor"
        if method != "local_diff_residual_mad":
            fallbacks.append(method)
        anchor_values.append(float(max(0.0, sigma)))
        anchor_methods.append(method)
    profile = np.interp(np.arange(len(masses)), anchors, np.asarray(anchor_values, dtype=float))
    method_profile = np.asarray([
        anchor_methods[int(np.argmin(np.abs(anchors - index)))] for index in range(len(masses))
    ], dtype=object)
    return profile, method_profile, float(global_sigma), sorted(set(fallbacks))


def _interpolated_mass(masses: np.ndarray, fractional_index: float) -> float:
    return float(np.interp(fractional_index, np.arange(len(masses), dtype=float), masses))


def _crossing_index(signal: np.ndarray, peak: int, direction: int, level: float, limit: int) -> float:
    index = peak
    while index != limit and signal[index] > level:
        index += direction
    previous = index - direction
    if index == previous or signal[previous] == signal[index]:
        return float(index)
    fraction = (signal[previous] - level) / (signal[previous] - signal[index])
    return float(previous + fraction * direction)


def _fwhm(
    signal: np.ndarray, masses: np.ndarray, peak: int, left_limit: int, right_limit: int,
) -> tuple[float, float, float, float]:
    if signal[peak] <= 0:
        return 0.0, 0.0, float(peak), float(peak)
    level = float(signal[peak]) / 2.0
    left = _crossing_index(signal, peak, -1, level, left_limit)
    right = _crossing_index(signal, peak, 1, level, right_limit)
    left_mass = _interpolated_mass(masses, left)
    right_mass = _interpolated_mass(masses, right)
    return float(right - left), float(right_mass - left_mass), left, right


def _suppress_by_mass_distance(
    candidates: list[dict[str, Any]], masses: np.ndarray, signal: np.ndarray, minimum_distance_da: float,
) -> tuple[list[dict[str, Any]], int]:
    kept: list[dict[str, Any]] = []
    for candidate in sorted(candidates, key=lambda row: (signal[row["index"]], row["prominence"]), reverse=True):
        if all(abs(masses[candidate["index"]] - masses[item["index"]]) >= minimum_distance_da for item in kept):
            kept.append(candidate)
    return sorted(kept, key=lambda row: row["index"]), len(candidates) - len(kept)


def _middle_minimum(signal: np.ndarray, left_peak: int, right_peak: int) -> int:
    region = signal[left_peak:right_peak + 1]
    minimum = float(np.min(region))
    indices = np.flatnonzero(region == minimum)
    return int(left_peak + indices[len(indices) // 2])


def _edge_boundary(
    corrected: np.ndarray, peak: int, direction: int, local_noise: float,
) -> tuple[int, str]:
    edge = 0 if direction < 0 else len(corrected) - 1
    index = peak
    while index != edge:
        index += direction
        if corrected[index] <= local_noise:
            return index, "noise_threshold_crossing"
    index = peak
    while index != edge:
        index += direction
        if corrected[index] <= 0:
            return index, "baseline_crossing"
    return edge, "profile_edge"


def _centroid_and_areas(
    masses: np.ndarray,
    raw: np.ndarray,
    quantification_weights: np.ndarray,
    left: int,
    right: int,
    apex_mass: float,
) -> tuple[float, float, float, bool]:
    local_mass = masses[left:right + 1]
    local_raw = raw[left:right + 1]
    local_weights = quantification_weights[left:right + 1]
    raw_area = float(_trapezoid(local_raw, local_mass)) if len(local_mass) > 1 else 0.0
    corrected_area = float(_trapezoid(local_weights, local_mass)) if len(local_mass) > 1 else 0.0
    if corrected_area > 0:
        centroid = float(_trapezoid(local_mass * local_weights, local_mass) / corrected_area)
        return centroid, raw_area, corrected_area, False
    return float(apex_mass), raw_area, corrected_area, True


def _provenance(
    parameters: SciexIntactPeakDetectionParameters,
    explicit: bool,
    derived: Mapping[str, tuple[Any, str, str] | tuple[Any, str, str, str]],
) -> tuple[SciexIntactParameterProvenance, ...]:
    output = []
    source = "explicit" if explicit else "default"
    units = {
        "baseline_window_da": "Da", "noise_window_da": "Da", "minimum_width_da": "Da",
        "minimum_distance_da": "Da", "severe_broad_peak_width_da": "Da",
        "prominence_window_da": "Da", "absolute_height_floor": "intensity",
        "absolute_prominence_floor": "intensity", "smoothing_window_points": "points",
    }
    for field in fields(parameters):
        output.append(SciexIntactParameterProvenance(
            field.name, getattr(parameters, field.name), units.get(field.name, "dimensionless"), source,
            "caller-supplied parameter" if explicit else "empirical SCIEX neutral-mass default",
        ))
    for name, item in derived.items():
        value, unit, reason = item[:3]
        derived_source = item[3] if len(item) == 4 else "derived"
        output.append(SciexIntactParameterProvenance(name, value, unit, derived_source, reason))
    return tuple(output)


def _empty_result(
    parameters: SciexIntactPeakDetectionParameters,
    profile_type: str,
    input_status: str,
    eligible: bool,
    detection_status: str,
    warnings: Iterable[str],
    diagnostics_extra: Mapping[str, Any] | None = None,
    *,
    explicit_parameters: bool,
) -> SciexIntactPeakDetectionResult:
    warning_list = tuple(warnings)
    values = {
        "Profile_Type": profile_type,
        "Input_Status": input_status,
        "Eligible_For_Neutral_Mass_Analysis": eligible,
        "Input_Validation_Status": "SUPPORTED" if detection_status == "SKIPPED_INELIGIBLE_PROFILE" else detection_status,
        "Detection_Status": detection_status,
        "Detection_Method": "SCIPY_SIGNAL_WITH_RNA_MASSHUNTER_QUANTIFICATION",
        "Algorithm_Version": ALGORITHM_VERSION,
        "Detected_Sensitive_Peak_Count": 0,
        "Detected_Strict_Peak_Count": 0,
        "Warning_Count": len(warning_list),
        "Automatic_Parameter_Fallbacks": "",
        **FORMAL_FALSE,
    }
    values.update(diagnostics_extra or {})
    return SciexIntactPeakDetectionResult(
        parameters, SciexIntactPeakDiagnostics(values), (),
        _provenance(parameters, explicit_parameters, {}), warnings=warning_list,
    )


def detect_sciex_intact_peaks(
    masses: Sequence[float] | np.ndarray,
    intensities: Sequence[float] | np.ndarray,
    *,
    profile_type: str,
    input_status: str,
    eligible_for_neutral_mass_analysis: bool,
    parameters: SciexIntactPeakDetectionParameters | None = None,
) -> SciexIntactPeakDetectionResult:
    """Detect peaks without assigning molecular identity or altering formal results."""
    explicit_parameters = parameters is not None
    params = parameters or SciexIntactPeakDetectionParameters()
    params.validate()
    eligible = bool(eligible_for_neutral_mass_analysis)
    if profile_type != NEUTRAL_MASS_PROFILE or input_status != SUPPORTED_INPUT or not eligible:
        return _empty_result(
            params, profile_type, input_status, eligible, "SKIPPED_INELIGIBLE_PROFILE",
            ["input_not_eligible_for_neutral_mass_peak_detection"],
            explicit_parameters=explicit_parameters,
        )
    try:
        mass = np.asarray(masses, dtype=float)
        raw = np.asarray(intensities, dtype=float)
    except (TypeError, ValueError):
        return _empty_result(
            params, profile_type, input_status, eligible, "INVALID_AXIS",
            ["input_arrays_are_not_numeric"], explicit_parameters=explicit_parameters,
        )
    if mass.ndim != 1 or raw.ndim != 1 or len(mass) != len(raw):
        return _empty_result(
            params, profile_type, input_status, eligible, "INVALID_AXIS",
            ["mass_and_intensity_must_be_one_dimensional_and_equal_length"],
            {"Parsed_Row_Count": min(mass.size, raw.size)}, explicit_parameters=explicit_parameters,
        )
    mass = mass.copy()
    raw = raw.copy()
    n = len(mass)
    missing_count = int(np.isnan(mass).sum() + np.isnan(raw).sum())
    nonfinite_mass = int((~np.isfinite(mass)).sum())
    nonfinite_intensity = int((~np.isfinite(raw)).sum())
    negative_count = int((raw < 0).sum())
    zero_count = int((raw == 0).sum())
    if n < params.minimum_points:
        return _empty_result(
            params, profile_type, input_status, eligible, "INSUFFICIENT_POINTS",
            ["profile_has_too_few_points"],
            {"Parsed_Row_Count": n, "Missing_Value_Count": missing_count,
             "Nonfinite_Value_Count": nonfinite_mass + nonfinite_intensity,
             "Negative_Intensity_Count": negative_count, "Zero_Intensity_Count": zero_count},
            explicit_parameters=explicit_parameters,
        )
    if nonfinite_mass:
        return _empty_result(
            params, profile_type, input_status, eligible, "INVALID_AXIS", ["mass_axis_contains_nonfinite_values"],
            {"Parsed_Row_Count": n, "Missing_Value_Count": missing_count,
             "Nonfinite_Value_Count": nonfinite_mass + nonfinite_intensity},
            explicit_parameters=explicit_parameters,
        )
    if nonfinite_intensity or negative_count:
        return _empty_result(
            params, profile_type, input_status, eligible, "INVALID_INTENSITY",
            ["intensity_contains_nonfinite_values" if nonfinite_intensity else "negative_intensity_is_not_supported"],
            {"Parsed_Row_Count": n, "Missing_Value_Count": missing_count,
             "Nonfinite_Value_Count": nonfinite_mass + nonfinite_intensity,
             "Negative_Intensity_Count": negative_count, "Zero_Intensity_Count": zero_count},
            explicit_parameters=explicit_parameters,
        )
    steps = np.diff(mass)
    duplicate_count = int(n - len(np.unique(mass)))
    if duplicate_count or np.any(steps <= 0):
        return _empty_result(
            params, profile_type, input_status, eligible, "INVALID_AXIS",
            ["mass_axis_must_be_strictly_increasing_without_duplicates"],
            {"Parsed_Row_Count": n, "Duplicate_Mass_Count": duplicate_count,
             "Mass_Axis_Strictly_Increasing": False}, explicit_parameters=explicit_parameters,
        )

    step_min = float(np.min(steps))
    step_median = float(np.median(steps))
    step_max = float(np.max(steps))
    uniform = _uniform_axis(steps, step_median)
    baseline_points = _odd_window_from_da(params.baseline_window_da, step_median, n)
    noise_points = _odd_window_from_da(params.noise_window_da, step_median, n)
    prominence_points = _odd_window_from_da(params.prominence_window_da, step_median, n)
    warning_list: list[str] = []
    fallback_list: list[str] = []
    if not uniform:
        warning_list.append("nonuniform_mass_axis_uses_mass_based_baseline_and_interpolated_widths")
    baseline, baseline_edge_mode = _rolling_quantile_baseline(
        mass, raw, params.baseline_quantile, params.baseline_window_da, uniform, baseline_points,
    )
    corrected = raw - baseline
    quantification_weights = np.maximum(corrected, 0.0)
    positive = corrected[corrected > 0]
    positive_floor = float(np.quantile(positive, params.positive_residual_quantile)) if positive.size else 0.0
    absolute_noise_fallback = max(params.absolute_height_floor, params.absolute_prominence_floor)
    # A sparse positive set is usually an isolated synthetic/real peak, not a noise population.
    # In that case using its lower quantile as sigma would make the peak reject itself.
    positive_noise_fallback = (
        positive_floor if positive.size >= max(10, int(0.10 * n)) else absolute_noise_fallback
    )
    local_noise, noise_methods, global_noise, noise_fallbacks = _local_noise_profile(
        mass, corrected, params.noise_window_da, positive_noise_fallback, absolute_noise_fallback,
    )
    fallback_list.extend(noise_fallbacks)

    smoothing_method = "NONE"
    smoothing_window = 1
    smoothing_polyorder: int | str = ""
    detection_signal = corrected.copy()
    if params.smoothing_enabled:
        smoothing_window = params.smoothing_window_points
        if smoothing_window > n:
            smoothing_window = n if n % 2 else n - 1
            fallback_list.append("smoothing_window_reduced_for_input_length")
        if smoothing_window > params.smoothing_polyorder and smoothing_window >= 3:
            detection_signal = savgol_filter(
                corrected, smoothing_window, params.smoothing_polyorder, mode="interp",
            )
            smoothing_method = "SAVITZKY_GOLAY"
            smoothing_polyorder = params.smoothing_polyorder
        else:
            smoothing_window = 1
            fallback_list.append("smoothing_disabled_for_input_length")

    # Savitzky-Golay can create symmetric overshoots inside an exact flat top.
    # Preserve raw plateaus in the detection-only signal so SciPy returns their midpoint once.
    raw_plateau_indices, raw_plateau_properties = find_peaks(raw, plateau_size=(2, None))
    for left, right in zip(
        raw_plateau_properties.get("left_edges", []),
        raw_plateau_properties.get("right_edges", []),
        strict=False,
    ):
        left_index = int(left)
        right_index = int(right)
        detection_signal[left_index:right_index + 1] = float(
            np.max(detection_signal[left_index:right_index + 1])
        )
    if raw_plateau_indices.size and smoothing_method == "SAVITZKY_GOLAY":
        smoothing_method = "SAVITZKY_GOLAY_EXACT_PLATEAU_PRESERVING"

    candidate_indices, plateau_properties = find_peaks(detection_signal, plateau_size=(1, None))
    if candidate_indices.size:
        with python_warnings.catch_warnings():
            python_warnings.simplefilter("ignore")
            prominences, left_bases, right_bases = peak_prominences(
                detection_signal, candidate_indices, wlen=prominence_points,
            )
            half_widths, _, left_ips, right_ips = peak_widths(
                detection_signal, candidate_indices, rel_height=0.5,
                prominence_data=(prominences, left_bases, right_bases),
            )
    else:
        prominences = np.array([], dtype=float)
        left_bases = right_bases = np.array([], dtype=int)
        half_widths = left_ips = right_ips = np.array([], dtype=float)

    height_thresholds = np.maximum.reduce([
        params.height_sigma_multiplier * local_noise,
        np.full(n, positive_floor),
        np.full(n, params.absolute_height_floor),
    ])
    prominence_thresholds = np.maximum(
        params.prominence_sigma_multiplier * local_noise, params.absolute_prominence_floor,
    )
    strict_thresholds = np.maximum(
        params.strict_prominence_sigma_multiplier * local_noise, params.absolute_prominence_floor,
    )
    rejected_height = rejected_prominence = rejected_width = 0
    sensitive_candidates: list[dict[str, Any]] = []
    for serial, peak_index in enumerate(candidate_indices):
        index = int(peak_index)
        left_mass = _interpolated_mass(mass, float(left_ips[serial]))
        right_mass = _interpolated_mass(mass, float(right_ips[serial]))
        width_da = float(right_mass - left_mass)
        height_pass = bool(detection_signal[index] >= height_thresholds[index])
        prominence_pass = bool(prominences[serial] >= prominence_thresholds[index])
        width_pass = bool(width_da >= params.minimum_width_da)
        if not height_pass:
            rejected_height += 1
            continue
        if not prominence_pass:
            rejected_prominence += 1
            continue
        if not width_pass:
            rejected_width += 1
            continue
        fwhm_points, fwhm_da, fwhm_left, fwhm_right = _fwhm(
            detection_signal, mass, index, int(left_bases[serial]), int(right_bases[serial]),
        )
        sensitive_candidates.append({
            "index": index,
            "prominence": float(prominences[serial]),
            "left_base": int(left_bases[serial]), "right_base": int(right_bases[serial]),
            "half_width_points": float(half_widths[serial]), "half_width_da": width_da,
            "half_left_ip": float(left_ips[serial]), "half_right_ip": float(right_ips[serial]),
            "fwhm_points": fwhm_points, "fwhm_da": fwhm_da,
            "fwhm_left_ip": fwhm_left, "fwhm_right_ip": fwhm_right,
            "plateau_left": int(plateau_properties["left_edges"][serial]),
            "plateau_right": int(plateau_properties["right_edges"][serial]),
            "strict": bool(prominences[serial] >= strict_thresholds[index]),
        })
    sensitive_candidates, suppressed_by_distance = _suppress_by_mass_distance(
        sensitive_candidates, mass, detection_signal, params.minimum_distance_da,
    )

    pair_valleys: dict[tuple[int, int], int] = {}
    for left_candidate, right_candidate in zip(sensitive_candidates, sensitive_candidates[1:]):
        pair_valleys[(left_candidate["index"], right_candidate["index"])] = _middle_minimum(
            detection_signal, left_candidate["index"], right_candidate["index"],
        )
    widths_da = np.asarray([item["half_width_da"] for item in sensitive_candidates], dtype=float)
    if len(widths_da) >= 5:
        broad_threshold = float(np.quantile(widths_da, params.broad_peak_quantile))
        broad_source = "automatic"
    else:
        broad_threshold = float(params.severe_broad_peak_width_da)
        broad_source = "fallback"
        fallback_list.append("broad_width_quantile_insufficient_peak_count")

    peak_dicts: list[dict[str, Any]] = []
    centroid_fallback_count = 0
    for position, candidate in enumerate(sensitive_candidates):
        index = candidate["index"]
        peak_id = f"SCIEX_INT_P{position + 1:05d}"
        left_neighbor = sensitive_candidates[position - 1] if position else None
        right_neighbor = sensitive_candidates[position + 1] if position + 1 < len(sensitive_candidates) else None
        fallback_parts: list[str] = []
        if left_neighbor is not None:
            left_boundary = pair_valleys[(left_neighbor["index"], index)]
        else:
            left_boundary, method = _edge_boundary(corrected, index, -1, float(local_noise[index]))
            fallback_parts.append(f"left:{method}")
        if right_neighbor is not None:
            right_boundary = pair_valleys[(index, right_neighbor["index"])]
        else:
            right_boundary, method = _edge_boundary(corrected, index, 1, float(local_noise[index]))
            fallback_parts.append(f"right:{method}")
        edge_peak = left_boundary == 0 or right_boundary == n - 1
        centroid, raw_area, corrected_area, centroid_fallback = _centroid_and_areas(
            mass, raw, quantification_weights, left_boundary, right_boundary, float(mass[index]),
        )
        if centroid_fallback:
            centroid_fallback_count += 1
            warning_list.append(f"{peak_id}:centroid_zero_weight_fallback_to_apex")
        plateau_left = candidate["plateau_left"]
        plateau_right = candidate["plateau_right"]
        plateau_center_mass = float((mass[plateau_left] + mass[plateau_right]) / 2.0)
        values = {
            "Peak_ID": peak_id,
            "Apex_Index": index,
            "Apex_Mass": float(mass[index]),
            "Apex_Intensity_Raw": float(raw[index]),
            "Local_Baseline_At_Apex": float(baseline[index]),
            "Apex_Intensity_Baseline_Corrected": float(corrected[index]),
            "Detection_Signal_Apex": float(detection_signal[index]),
            "Local_Noise_Sigma": float(local_noise[index]),
            "Noise_Estimation_Method_Local": str(noise_methods[index]),
            "Noise_Fallback_Used": str(noise_methods[index]) != "local_diff_residual_mad",
            "Height_Threshold_Local": float(height_thresholds[index]),
            "Prominence": candidate["prominence"],
            "Prominence_Threshold_Local": float(prominence_thresholds[index]),
            "Strict_Prominence_Threshold_Local": float(strict_thresholds[index]),
            "Detection_Tier": "STRICT" if candidate["strict"] else "SENSITIVE",
            "Sensitive_Threshold_Passed": True,
            "Strict_Threshold_Passed": candidate["strict"],
            "Molecular_Identity_Assigned": False,
            "Prominence_Base_Left_Index": candidate["left_base"],
            "Prominence_Base_Right_Index": candidate["right_base"],
            "Prominence_Base_Left_Mass": float(mass[candidate["left_base"]]),
            "Prominence_Base_Right_Mass": float(mass[candidate["right_base"]]),
            "Half_Prominence_Width_Points": candidate["half_width_points"],
            "Half_Prominence_Width_Da": candidate["half_width_da"],
            "Half_Prominence_Left_IP": candidate["half_left_ip"],
            "Half_Prominence_Right_IP": candidate["half_right_ip"],
            "FWHM_Points": candidate["fwhm_points"],
            "FWHM_Da": candidate["fwhm_da"],
            "Left_Boundary_Index": left_boundary,
            "Right_Boundary_Index": right_boundary,
            "Left_Boundary_Mass": float(mass[left_boundary]),
            "Right_Boundary_Mass": float(mass[right_boundary]),
            "Boundary_Width_Da": float(mass[right_boundary] - mass[left_boundary]),
            "Boundary_Method": "SENSITIVE_ACCEPTED_PEAK_SHARED_VALLEY",
            "Boundary_Fallback_Used": ";".join(fallback_parts),
            "Boundary_Left_Neighbor_Peak_ID": f"SCIEX_INT_P{position:05d}" if left_neighbor else "",
            "Boundary_Right_Neighbor_Peak_ID": f"SCIEX_INT_P{position + 2:05d}" if right_neighbor else "",
            "Boundary_Peak_Set_Tier": "SENSITIVE",
            "Boundary_Recomputed_After_Filtering": False,
            "Centroid_Mass": centroid,
            "Centroid_Minus_Apex_Da": float(centroid - mass[index]),
            "Centroid_Fallback_Used": centroid_fallback,
            "Peak_Area_Raw": raw_area,
            "Peak_Area_Baseline_Corrected": corrected_area,
            "Area_Unit": "intensity_x_Da",
            "Plateau_Start_Index": plateau_left,
            "Plateau_End_Index": plateau_right,
            "Plateau_Size_Points": plateau_right - plateau_left + 1,
            "Plateau_Width_Da": float(mass[plateau_right] - mass[plateau_left]),
            "Plateau_Center_Mass": plateau_center_mass,
            "Broad_Peak_Flag": candidate["half_width_da"] >= broad_threshold,
            "Severe_Broad_Peak_Flag": candidate["half_width_da"] >= params.severe_broad_peak_width_da,
            "Edge_Peak_Flag": edge_peak,
            "Peak_Area_Complete": not edge_peak,
            "Centroid_Complete": not edge_peak and not centroid_fallback,
            **FORMAL_FALSE,
        }
        peak_dicts.append(values)

    shallow_count = shoulder_count = 0
    for position, values in enumerate(peak_dicts):
        neighbors: list[tuple[int, int]] = []
        if position:
            neighbors.append((position - 1, pair_valleys[(sensitive_candidates[position - 1]["index"], sensitive_candidates[position]["index"])]))
        if position + 1 < len(peak_dicts):
            neighbors.append((position + 1, pair_valleys[(sensitive_candidates[position]["index"], sensitive_candidates[position + 1]["index"])]))
        if not neighbors:
            values.update({
                "Neighbor_Peak_ID": "", "Neighbor_Separation_Da": "", "Shared_Valley_Index": "",
                "Shared_Valley_Mass": "", "Valley_To_Smaller_Apex_Ratio": "",
                "Width_To_Separation_Ratio": "", "Shallow_Valley_Neighbor_Flag": False,
                "Shoulder_Diagnostic_Reason": "no_neighbor", "Possible_Shoulder": False,
            })
            continue
        neighbor_position, valley = min(
            neighbors, key=lambda item: abs(values["Apex_Mass"] - peak_dicts[item[0]]["Apex_Mass"]),
        )
        neighbor = peak_dicts[neighbor_position]
        separation = abs(float(values["Apex_Mass"]) - float(neighbor["Apex_Mass"]))
        smaller_apex = min(float(values["Detection_Signal_Apex"]), float(neighbor["Detection_Signal_Apex"]))
        valley_ratio = float(detection_signal[valley] / smaller_apex) if smaller_apex > 0 else float("nan")
        width_ratio = float(values["Half_Prominence_Width_Da"] / separation) if separation > 0 else float("inf")
        shallow = bool(np.isfinite(valley_ratio) and valley_ratio >= params.shoulder_valley_ratio_threshold)
        possible = shallow and (
            separation <= params.shoulder_max_separation_da
            or width_ratio >= params.shoulder_width_to_separation_ratio
        )
        if shallow:
            shallow_count += 1
        if possible:
            shoulder_count += 1
        values.update({
            "Neighbor_Peak_ID": neighbor["Peak_ID"],
            "Neighbor_Separation_Da": separation,
            "Shared_Valley_Index": valley,
            "Shared_Valley_Mass": float(mass[valley]),
            "Valley_To_Smaller_Apex_Ratio": valley_ratio,
            "Width_To_Separation_Ratio": width_ratio,
            "Shallow_Valley_Neighbor_Flag": shallow,
            "Shoulder_Diagnostic_Reason": "shallow_valley_close_or_width_overlap" if possible else "criteria_not_met",
            "Possible_Shoulder": possible,
        })

    fallback_list = sorted(set(fallback_list))
    if fallback_list:
        warning_list.extend(f"automatic_parameter_fallback:{item}" for item in fallback_list)
    if not uniform:
        warning_list.append("minimum_distance_enforced_in_mass_coordinates")
    warning_list = list(dict.fromkeys(warning_list))
    detection_status = "DETECTION_COMPLETED_WITH_WARNINGS" if warning_list else "DETECTION_COMPLETED"
    strict_count = sum(bool(row["Strict_Threshold_Passed"]) for row in peak_dicts)
    edge_count = sum(bool(row["Edge_Peak_Flag"]) for row in peak_dicts)
    broad_count = sum(bool(row["Broad_Peak_Flag"]) for row in peak_dicts)
    severe_broad_count = sum(bool(row["Severe_Broad_Peak_Flag"]) for row in peak_dicts)
    derived = {
        "Baseline_Window_Points": (baseline_points, "points", "derived from Da and median mass step; odd window"),
        "Noise_Window_Points": (noise_points, "points", "derived from Da and median mass step; odd window"),
        "Prominence_Window_Points": (prominence_points, "points", "derived from Da and median mass step; odd window"),
        "Minimum_Distance_Points": (max(1, int(ceil(params.minimum_distance_da / step_median))), "points", "reporting equivalent; enforcement uses mass coordinates"),
        "Minimum_Width_Points": (params.minimum_width_da / step_median, "points", "reporting equivalent; filtering uses interpolated mass width"),
        "Positive_Residual_Quantile_Value": (positive_floor, "intensity", "automatic lower positive residual floor", "automatic"),
        "Estimated_Noise_Global": (global_noise, "intensity", "global diff-residual MAD fallback", "automatic"),
        "Broad_Peak_Threshold_Da": (broad_threshold, "Da", f"{broad_source} accepted-peak width threshold", broad_source),
    }
    diagnostics_values = {
        "Profile_Type": profile_type,
        "Input_Status": input_status,
        "Eligible_For_Neutral_Mass_Analysis": eligible,
        "Input_Validation_Status": "SUPPORTED",
        "Detection_Status": detection_status,
        "Detection_Method": "SCIPY_SIGNAL_WITH_RNA_MASSHUNTER_QUANTIFICATION",
        "Algorithm_Version": ALGORITHM_VERSION,
        "Parsed_Row_Count": n,
        "Mass_Min_Da": float(mass[0]), "Mass_Max_Da": float(mass[-1]),
        "Mass_Step_Min_Da": step_min, "Mass_Step_Median_Da": step_median, "Mass_Step_Max_Da": step_max,
        "Mass_Axis_Strictly_Increasing": True, "Mass_Axis_Uniform": uniform,
        "Duplicate_Mass_Count": 0, "Missing_Value_Count": missing_count,
        "Nonfinite_Value_Count": 0, "Negative_Intensity_Count": 0, "Zero_Intensity_Count": zero_count,
        "Baseline_Method": "ROLLING_QUANTILE", "Baseline_Quantile": params.baseline_quantile,
        "Baseline_Window_Points": baseline_points, "Baseline_Window_Da": params.baseline_window_da,
        "Baseline_Edge_Mode": baseline_edge_mode,
        "Baseline_Min": float(np.min(baseline)), "Baseline_Median": float(np.median(baseline)), "Baseline_Max": float(np.max(baseline)),
        "Baseline_Negative_Residual_Fraction": float(np.mean(corrected < 0)),
        "Noise_Estimation_Method": "LOCAL_DIFF_RESIDUAL_MAD_INTERPOLATED",
        "Noise_Window_Points": noise_points, "Noise_Window_Da": params.noise_window_da,
        "Estimated_Noise_Global": global_noise,
        "Estimated_Noise_Local_Min": float(np.min(local_noise)),
        "Estimated_Noise_Local_Median": float(np.median(local_noise)),
        "Estimated_Noise_Local_Max": float(np.max(local_noise)),
        "Height_Threshold_Method": "MAX_HEIGHT_SIGMA_POSITIVE_Q_ABSOLUTE_FLOOR",
        "Prominence_Threshold_Method": "MAX_PROMINENCE_SIGMA_ABSOLUTE_FLOOR",
        "Strict_Prominence_Threshold_Method": "MAX_STRICT_SIGMA_ABSOLUTE_FLOOR",
        "Positive_Residual_Quantile_Value": positive_floor,
        "Smoothing_Method": smoothing_method, "Smoothing_Window_Points": smoothing_window,
        "Smoothing_Window_Da": float((smoothing_window - 1) * step_median), "Smoothing_Polyorder": smoothing_polyorder,
        "Minimum_Distance_Points": max(1, int(ceil(params.minimum_distance_da / step_median))),
        "Minimum_Distance_Da": params.minimum_distance_da,
        "Minimum_Width_Points": params.minimum_width_da / step_median,
        "Minimum_Width_Da": params.minimum_width_da,
        "Boundary_Method": "SENSITIVE_ACCEPTED_PEAK_SHARED_VALLEY",
        "Boundary_Peak_Set_Tier": "SENSITIVE", "Boundary_Recomputed_After_Filtering": False,
        "Centroid_Method": "POSITIVE_BASELINE_CORRECTED_TRAPEZOID_WEIGHTED",
        "Area_Method": "RAW_AND_POSITIVE_BASELINE_CORRECTED_TRAPEZOID",
        "Detected_Sensitive_Peak_Count": len(peak_dicts), "Detected_Strict_Peak_Count": strict_count,
        "Rejected_Height_Count": rejected_height, "Rejected_Prominence_Count": rejected_prominence,
        "Rejected_Width_Count": rejected_width, "Suppressed_By_Distance_Count": suppressed_by_distance,
        "Shallow_Valley_Neighbor_Count": shallow_count, "Possible_Shoulder_Count": shoulder_count,
        "Broad_Peak_Width_Threshold_Da": broad_threshold, "Broad_Peak_Threshold_Source": broad_source,
        "Broad_Peak_Count": broad_count, "Severe_Broad_Peak_Count": severe_broad_count,
        "Edge_Peak_Count": edge_count, "Centroid_Fallback_Count": centroid_fallback_count,
        "Warning_Count": len(warning_list), "Automatic_Parameter_Fallbacks": ";".join(fallback_list),
        **FORMAL_FALSE,
    }
    return SciexIntactPeakDetectionResult(
        params,
        SciexIntactPeakDiagnostics(diagnostics_values),
        tuple(SciexIntactDetectedPeak(values) for values in peak_dicts),
        _provenance(params, explicit_parameters, derived),
        _as_tuple(raw), _as_tuple(baseline), _as_tuple(corrected), _as_tuple(detection_signal),
        _as_tuple(quantification_weights), tuple(warning_list),
    )
