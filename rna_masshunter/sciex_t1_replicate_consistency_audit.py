"""Independent MS1 technical-replicate consistency shadow audit for SCIEX mzML."""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import asdict, dataclass, field, replace
from enum import Enum
from math import log2
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Mapping, Sequence
import re

import numpy as np
from pyteomics import mzml
from scipy.signal import find_peaks, peak_prominences

from rna_masshunter.sciex_mzml_source_metadata_audit import (
    MzMLSourceMetadataRecord,
    PolarityStatus,
    RepresentationStatus,
)
from rna_masshunter.sciex_t1_profile_peak_audit import (
    T1PeakDetectionParameters,
    T1PeakQualityClass,
    detect_t1_profile_peaks,
)

OPTIONAL_RESULT_KEY = "sciex_t1_replicate_consistency_audit"
ALGORITHM_VERSION = "sciex-t1-replicate-consistency-audit-v1"

_BLOCK_ORDER = (
    "INPUT_FILE_NOT_FOUND", "INPUT_FILE_UNREADABLE", "SOURCE_METADATA_RECORD_MISSING",
    "USER_MANIFEST_CONTEXT_MISSING", "REPLICATE_CONTEXT_MISMATCH",
    "POLARITY_MISMATCH_BETWEEN_REPLICATES", "MIXED_POLARITY_INPUT",
    "MISSING_POLARITY_METADATA", "REPRESENTATION_MISMATCH_BETWEEN_REPLICATES",
    "MIXED_REPRESENTATION_INPUT", "MISSING_REPRESENTATION_METADATA",
    "MISSING_MS_LEVEL_METADATA", "NON_MS1_SPECTRUM_EXCLUDED", "NO_MS1_SPECTRA",
    "INSUFFICIENT_MS1_SPECTRA", "PROFILE_EXTRACTION_FAILED", "PEAK_DETECTION_FAILED",
    "NO_DETECTED_PEAKS_RUN_1", "NO_DETECTED_PEAKS_RUN_2", "NO_MATCHED_PEAKS",
    "INSUFFICIENT_MATCHES_FOR_DRIFT", "AMBIGUOUS_PEAK_MATCH", "LOW_SCAN_RECURRENCE",
    "LOW_PROMINENCE", "INVALID_FWHM", "EXTREME_INTENSITY_VARIATION",
    "MISSING_INTENSITY_NORMALIZATION", "MISSING_PEAK_SHAPE_METRICS",
)
_LOW_QUALITY_CLASSES = {T1PeakQualityClass.LOW_SUPPORT.value, T1PeakQualityClass.SHOULDER_OR_OVERLAP.value}


def _ordered_blocks(values: Iterable[str]) -> tuple[str, ...]:
    found = set(values)
    return tuple(x for x in _BLOCK_ORDER if x in found) + tuple(sorted(found - set(_BLOCK_ORDER)))


def _safe_label(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_").upper() or "RUN"


def _ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator is None or denominator <= 0:
        return None
    return numerator / denominator


def _median(values: Iterable[float | int | None]) -> float | None:
    usable = [float(x) for x in values if x is not None and np.isfinite(float(x))]
    return float(median(usable)) if usable else None


class MatchAmbiguityStatus(str, Enum):
    UNAMBIGUOUS_ONE_TO_ONE = "UNAMBIGUOUS_ONE_TO_ONE"
    AMBIGUOUS_LEFT_TO_MULTIPLE = "AMBIGUOUS_LEFT_TO_MULTIPLE"
    AMBIGUOUS_RIGHT_TO_MULTIPLE = "AMBIGUOUS_RIGHT_TO_MULTIPLE"
    AMBIGUOUS_MANY_TO_MANY = "AMBIGUOUS_MANY_TO_MANY"
    UNMATCHED = "UNMATCHED"


class ReplicateConsistencyStatus(str, Enum):
    REPRODUCED_HIGH_CONFIDENCE = "REPRODUCED_HIGH_CONFIDENCE"
    REPRODUCED_SUPPORTIVE = "REPRODUCED_SUPPORTIVE"
    REPRODUCED_WITH_INTENSITY_VARIATION = "REPRODUCED_WITH_INTENSITY_VARIATION"
    REPRODUCED_WITH_SHAPE_VARIATION = "REPRODUCED_WITH_SHAPE_VARIATION"
    REPRODUCED_AMBIGUOUS_MATCH = "REPRODUCED_AMBIGUOUS_MATCH"
    RUN_1_ONLY = "RUN_1_ONLY"
    RUN_2_ONLY = "RUN_2_ONLY"
    LOW_RECURRENCE_BOTH_RUNS = "LOW_RECURRENCE_BOTH_RUNS"
    INSUFFICIENT_PEAK_QUALITY = "INSUFFICIENT_PEAK_QUALITY"
    UNRESOLVED = "UNRESOLVED"


class SystematicMZShiftStatus(str, Enum):
    NO_MEANINGFUL_SHIFT = "NO_MEANINGFUL_SHIFT"
    SMALL_SYSTEMATIC_SHIFT = "SMALL_SYSTEMATIC_SHIFT"
    POSSIBLE_SYSTEMATIC_SHIFT = "POSSIBLE_SYSTEMATIC_SHIFT"
    INSUFFICIENT_MATCHES = "INSUFFICIENT_MATCHES"


class OverallReplicateStatus(str, Enum):
    STRONG_TECHNICAL_REPRODUCIBILITY = "STRONG_TECHNICAL_REPRODUCIBILITY"
    MODERATE_TECHNICAL_REPRODUCIBILITY = "MODERATE_TECHNICAL_REPRODUCIBILITY"
    LIMITED_TECHNICAL_REPRODUCIBILITY = "LIMITED_TECHNICAL_REPRODUCIBILITY"
    SUBSTANTIAL_RUN_SPECIFIC_SIGNAL = "SUBSTANTIAL_RUN_SPECIFIC_SIGNAL"
    INSUFFICIENT_COMPARABLE_PEAKS = "INSUFFICIENT_COMPARABLE_PEAKS"
    BLOCKED_BY_INPUT_QUALITY = "BLOCKED_BY_INPUT_QUALITY"


@dataclass(frozen=True)
class ReplicateAuditParameters:
    mz_grid_step: float = 0.005
    minimum_valid_scan_points: int = 3
    minimum_ms1_spectra: int = 3
    scan_peak_minimum_relative_prominence: float = 0.0005
    maximum_scan_peaks: int = 2000
    absolute_tolerance_da: float = 0.01
    ppm_tolerance: float = 10.0
    strict_absolute_tolerance_da: float = 0.005
    strict_ppm_tolerance: float = 5.0
    ambiguity_margin_da: float = 0.002
    high_recurrence_fraction: float = 0.02
    supportive_recurrence_fraction: float = 0.005
    minimum_relative_prominence: float = 0.001
    high_rank_difference: int = 10
    extreme_log2_intensity_ratio: float = 2.0
    minimum_fwhm_ratio: float = 0.5
    maximum_fwhm_ratio: float = 2.0
    minimum_drift_matches: int = 3
    no_shift_delta_da: float = 0.002
    small_shift_delta_da: float = 0.01
    strong_jaccard: float = 0.60
    moderate_jaccard: float = 0.30

    def validate(self) -> None:
        positive = (
            "mz_grid_step", "minimum_valid_scan_points", "minimum_ms1_spectra",
            "maximum_scan_peaks", "absolute_tolerance_da", "ppm_tolerance",
            "strict_absolute_tolerance_da", "strict_ppm_tolerance", "minimum_drift_matches",
        )
        if any(float(getattr(self, name)) <= 0 for name in positive):
            raise ValueError("positive replicate-audit parameters required")
        fractions = (
            "scan_peak_minimum_relative_prominence", "high_recurrence_fraction",
            "supportive_recurrence_fraction", "minimum_relative_prominence",
            "minimum_fwhm_ratio", "strong_jaccard", "moderate_jaccard",
        )
        if any(not 0 <= float(getattr(self, name)) <= 1 for name in fractions):
            raise ValueError("fraction parameters must be in [0,1]")
        if self.supportive_recurrence_fraction > self.high_recurrence_fraction:
            raise ValueError("supportive recurrence cannot exceed high recurrence")
        if self.moderate_jaccard > self.strong_jaccard:
            raise ValueError("moderate Jaccard cannot exceed strong Jaccard")


@dataclass(frozen=True, kw_only=True)
class ReplicateRunPeak:
    run_label: str
    peak_id: str
    apex_mz: float
    centroid_mz: float | None
    neutral_or_observed_space: str = "OBSERVED_MZ"
    raw_apex_intensity: float
    normalized_apex_intensity: float
    raw_integrated_intensity: float
    normalized_integrated_intensity: float
    relative_intensity: float
    intensity_rank: int
    prominence: float | None
    relative_prominence: float | None
    fwhm: float | None
    left_bound_mz: float
    right_bound_mz: float
    supporting_ms1_scan_count: int
    total_ms1_scan_count: int
    scan_recurrence_fraction: float
    first_supporting_scan_time: float | None
    last_supporting_scan_time: float | None
    detection_status: str
    detection_block_reasons: tuple[str, ...]
    formal_propagation: bool = False
    chemical_identity_assigned: bool = False
    charge_state_confirmed: bool = False


@dataclass(frozen=True)
class ReplicateRunPeakProfile:
    run_label: str
    input_path: str
    status: str
    aggregation_method: str
    ms1_spectra_total: int
    ms1_spectra_used: int
    ms1_spectra_excluded: int
    ms2_spectra_excluded: int
    missing_ms_level_spectra: int
    mz_grid_method: str
    intensity_normalization_method: str
    baseline_method: str
    smoothing_method: str
    peak_detection_method: str
    detected_peak_count: int
    peaks: tuple[ReplicateRunPeak, ...]
    polarity_status: str
    representation_status: str
    block_reasons: tuple[str, ...]
    technical_replicate_only: bool = True
    biological_replicate_claim: bool = False
    chemical_identity_assigned: bool = False
    formal_propagation: bool = False
    comparison_mz_grid: np.ndarray | None = field(default=None, repr=False, compare=False)
    comparison_raw_profile: np.ndarray | None = field(default=None, repr=False, compare=False)
    comparison_normalized_profile: np.ndarray | None = field(default=None, repr=False, compare=False)


@dataclass(frozen=True, kw_only=True)
class ReplicatePeakMatch:
    replicate_pair_id: str
    run_1_label: str
    run_2_label: str
    run_1_peak_id: str | None
    run_2_peak_id: str | None
    run_1_apex_mz: float | None
    run_2_apex_mz: float | None
    delta_mz: float | None
    absolute_delta_mz: float | None
    ppm_error: float | None
    run_1_centroid_mz: float | None
    run_2_centroid_mz: float | None
    centroid_delta_mz: float | None
    run_1_normalized_intensity: float | None
    run_2_normalized_intensity: float | None
    normalized_intensity_ratio: float | None
    log2_normalized_intensity_ratio: float | None
    run_1_intensity_rank: int | None
    run_2_intensity_rank: int | None
    intensity_rank_difference: int | None
    run_1_prominence: float | None
    run_2_prominence: float | None
    prominence_ratio: float | None
    run_1_fwhm: float | None
    run_2_fwhm: float | None
    fwhm_ratio: float | None
    run_1_scan_recurrence_fraction: float | None
    run_2_scan_recurrence_fraction: float | None
    candidate_match_count_left: int
    candidate_match_count_right: int
    match_ambiguity_status: MatchAmbiguityStatus
    alternative_match_peak_ids: tuple[str, ...]
    replicate_consistency_status: ReplicateConsistencyStatus
    replicate_consistency_confidence: str
    consistency_block_reasons: tuple[str, ...]
    formal_propagation: bool = False
    technical_replicate_only: bool = True
    biological_replicate_claim: bool = False
    chemical_identity_assigned: bool = False


@dataclass(frozen=True)
class MZDriftSummary:
    median_delta_mz: float | None
    median_ppm_error: float | None
    mad_delta_mz: float | None
    mad_ppm_error: float | None
    matched_peak_count_for_drift: int
    systematic_mz_shift_status: SystematicMZShiftStatus
    drift_adjustment_applied: bool = False


@dataclass(frozen=True, kw_only=True)
class ReplicateConsistencySummary:
    run_1_label: str
    run_2_label: str
    run_1_detected_peak_count: int
    run_2_detected_peak_count: int
    matched_peak_pair_count: int
    run_1_only_peak_count: int
    run_2_only_peak_count: int
    ambiguous_match_count: int
    high_confidence_reproduced_count: int
    supportive_reproduced_count: int
    low_quality_reproduced_count: int
    union_peak_count: int
    intersection_peak_count: int
    jaccard_index: float
    run_1_overlap_fraction: float
    run_2_overlap_fraction: float
    median_absolute_delta_mz: float | None
    median_ppm_error: float | None
    median_intensity_rank_difference: float | None
    median_log2_normalized_intensity_ratio: float | None
    drift: MZDriftSummary
    replicate_consistency_overall_status: OverallReplicateStatus
    overall_confidence: str
    overall_block_reasons: tuple[str, ...]
    technical_replicate_only: bool = True
    biological_replicate_claim: bool = False
    chemical_identity_assigned: bool = False
    formal_propagation: bool = False


@dataclass(frozen=True)
class ReplicateConsistencyAuditResult:
    parameters: ReplicateAuditParameters
    run_profiles: tuple[ReplicateRunPeakProfile, ...]
    matches: tuple[ReplicatePeakMatch, ...]
    summaries: tuple[ReplicateConsistencySummary, ...]
    algorithm_version: str = ALGORITHM_VERSION
    formal_propagation: bool = False


def _rt_minutes(spectrum: Mapping[str, Any]) -> float | None:
    scans = spectrum.get("scanList", {}).get("scan", [])
    if not scans:
        return None
    scan = scans[0]
    value = scan.get("scan start time")
    if value is None:
        return None
    unit = str(scan.get("unitName", "")).lower()
    return float(value) / 60.0 if "second" in unit else float(value)


def _decode_array(value: Any) -> np.ndarray:
    if hasattr(value, "decode") and not isinstance(value, (bytes, bytearray, str)):
        value = value.decode()
    return np.asarray(value if value is not None else (), dtype=float)


def _scan_peak_mzs(mz_values: np.ndarray, intensities: np.ndarray, parameters: ReplicateAuditParameters) -> np.ndarray:
    if len(mz_values) < parameters.minimum_valid_scan_points:
        return np.array([], dtype=float)
    positive = np.maximum(intensities, 0)
    base = float(positive.max())
    if base <= 0:
        return np.array([], dtype=float)
    diffs = np.diff(mz_values)
    positive_diffs = diffs[diffs > 0]
    spacing = float(np.median(positive_diffs)) if len(positive_diffs) else parameters.mz_grid_step
    gap = max(0.05, spacing * 5)
    cuts = np.flatnonzero(diffs > gap) + 1
    starts = np.r_[0, cuts]
    ends = np.r_[cuts, len(mz_values)]
    candidates: list[tuple[float, float]] = []
    for start, end in zip(starts, ends, strict=False):
        if end - start < 3:
            continue
        segment = positive[start:end]
        distance = max(1, int(np.ceil(0.015 / max(spacing, 1e-12))))
        indices, _ = find_peaks(segment, distance=distance, prominence=base * parameters.scan_peak_minimum_relative_prominence)
        if not len(indices):
            continue
        prominences = peak_prominences(segment, indices)[0]
        candidates.extend((float(prom), float(mz_values[start + index])) for index, prom in zip(indices, prominences, strict=False))
    candidates.sort(key=lambda item: (-item[0], item[1]))
    return np.asarray(sorted(mz for _, mz in candidates[:parameters.maximum_scan_peaks]), dtype=float)


def _add_binned_scan(
    mz_values: np.ndarray,
    intensities: np.ndarray,
    parameters: ReplicateAuditParameters,
    state: dict[str, Any],
) -> None:
    bins = np.rint(mz_values / parameters.mz_grid_step).astype(np.int64)
    order = np.argsort(bins, kind="stable")
    bins = bins[order]
    values = np.maximum(intensities[order], 0)
    unique, starts = np.unique(bins, return_index=True)
    maxima = np.maximum.reduceat(values, starts)
    low, high = int(unique[0]), int(unique[-1])
    if state["origin"] is None:
        state["origin"] = low
        size = high - low + 1
        state["raw"] = np.zeros(size, dtype=float)
        state["normalized"] = np.zeros(size, dtype=float)
    elif low < state["origin"] or high >= state["origin"] + len(state["raw"]):
        new_origin = min(low, state["origin"])
        new_high = max(high, state["origin"] + len(state["raw"]) - 1)
        size = new_high - new_origin + 1
        raw = np.zeros(size, dtype=float)
        normalized = np.zeros(size, dtype=float)
        offset = state["origin"] - new_origin
        raw[offset:offset + len(state["raw"])] = state["raw"]
        normalized[offset:offset + len(state["normalized"])] = state["normalized"]
        state.update(origin=new_origin, raw=raw, normalized=normalized)
    index = unique - state["origin"]
    state["raw"][index] += maxima
    base = float(maxima.max())
    if base > 0:
        state["normalized"][index] += maxima / base


def _metadata_status(metadata_record: MzMLSourceMetadataRecord | None) -> tuple[str, str, list[str]]:
    if metadata_record is None:
        return "UNKNOWN", "UNKNOWN", ["SOURCE_METADATA_RECORD_MISSING"]
    blocks: list[str] = []
    if metadata_record.context_source != "USER_PROVIDED_RUNTIME_MANIFEST":
        blocks.append("USER_MANIFEST_CONTEXT_MISSING")
    polarity = metadata_record.polarity_status.value
    representation = metadata_record.representation_status.value
    if metadata_record.polarity_status is PolarityStatus.MIXED_POLARITY:
        blocks.append("MIXED_POLARITY_INPUT")
    elif metadata_record.polarity_status is PolarityStatus.NOT_RECORDED:
        blocks.append("MISSING_POLARITY_METADATA")
    if metadata_record.representation_status is RepresentationStatus.MIXED_REPRESENTATION:
        blocks.append("MIXED_REPRESENTATION_INPUT")
    elif metadata_record.representation_status is RepresentationStatus.NOT_RECORDED:
        blocks.append("MISSING_REPRESENTATION_METADATA")
    return polarity, representation, blocks


def build_ms1_peak_profile_from_spectra(
    spectra: Iterable[Mapping[str, Any]], *, run_label: str,
    input_path: str = "IN_MEMORY", metadata_record: MzMLSourceMetadataRecord | None = None,
    parameters: ReplicateAuditParameters | None = None,
    peak_detection_parameters: T1PeakDetectionParameters | None = None,
) -> ReplicateRunPeakProfile:
    params = parameters or ReplicateAuditParameters()
    params.validate()
    polarity, representation, blocks = _metadata_status(metadata_record)
    state: dict[str, Any] = {"origin": None, "raw": None, "normalized": None}
    scan_observations: list[tuple[float | None, np.ndarray]] = []
    ms1_total = ms1_used = ms1_excluded = ms2_excluded = missing_level = 0
    for spectrum in spectra:
        raw_level = spectrum.get("ms level")
        try:
            level = int(raw_level)
        except (TypeError, ValueError):
            missing_level += 1
            continue
        if level != 1:
            ms2_excluded += 1
            continue
        ms1_total += 1
        try:
            mz_values = _decode_array(spectrum.get("m/z array"))
            intensities = _decode_array(spectrum.get("intensity array"))
        except Exception:
            ms1_excluded += 1
            continue
        valid = (
            len(mz_values) == len(intensities) and len(mz_values) >= params.minimum_valid_scan_points
            and np.all(np.isfinite(mz_values)) and np.all(np.isfinite(intensities))
            and np.all(np.diff(mz_values) > 0)
        )
        if not valid:
            ms1_excluded += 1
            continue
        ms1_used += 1
        _add_binned_scan(mz_values, intensities, params, state)
        scan_observations.append((_rt_minutes(spectrum), _scan_peak_mzs(mz_values, intensities, params)))
    if missing_level:
        blocks.append("MISSING_MS_LEVEL_METADATA")
    if ms2_excluded:
        blocks.append("NON_MS1_SPECTRUM_EXCLUDED")
    if not ms1_total:
        blocks.append("NO_MS1_SPECTRA")
    if 0 < ms1_used < params.minimum_ms1_spectra:
        blocks.append("INSUFFICIENT_MS1_SPECTRA")
    peaks: tuple[ReplicateRunPeak, ...] = ()
    status = "COMPLETED"
    if ms1_used == 0 or state["origin"] is None:
        status = "BLOCKED"
        blocks.append("PROFILE_EXTRACTION_FAILED")
    else:
        grid = (state["origin"] + np.arange(len(state["normalized"]))) * params.mz_grid_step
        raw_mean = state["raw"] / ms1_used
        normalized_mean = state["normalized"] / ms1_used
        smoothed = np.convolve(normalized_mean, np.array([0.25, 0.5, 0.25]), mode="same")
        try:
            detected = detect_t1_profile_peaks(
                grid, smoothed, source_id=run_label, measurement_id=run_label,
                rna_identity=metadata_record.rna_identity if metadata_record else "UNKNOWN",
                parameters=peak_detection_parameters,
            )
        except Exception:
            detected = None
            status = "BLOCKED"
            blocks.append("PEAK_DETECTION_FAILED")
        if detected is not None:
            interim: list[ReplicateRunPeak] = []
            for serial, peak in enumerate(detected.peaks, 1):
                tolerance = max(params.absolute_tolerance_da, peak.apex_mz * params.ppm_tolerance * 1e-6)
                supporting_times: list[float] = []
                support = 0
                for rt, observed in scan_observations:
                    index = bisect_left(observed, peak.apex_mz - tolerance)
                    if index < len(observed) and observed[index] <= peak.apex_mz + tolerance:
                        support += 1
                        if rt is not None:
                            supporting_times.append(float(rt))
                apex_index = int(np.clip(round(peak.apex_mz / params.mz_grid_step) - state["origin"], 0, len(raw_mean) - 1))
                left_index = int(np.clip(round(peak.left_boundary_mz / params.mz_grid_step) - state["origin"], 0, len(raw_mean) - 1))
                right_index = int(np.clip(round(peak.right_boundary_mz / params.mz_grid_step) - state["origin"], left_index, len(raw_mean) - 1))
                raw_area = float(np.trapezoid(raw_mean[left_index:right_index + 1], grid[left_index:right_index + 1])) if right_index > left_index else 0.0
                peak_blocks: list[str] = []
                recurrence = support / ms1_used
                if recurrence < params.supportive_recurrence_fraction:
                    peak_blocks.append("LOW_SCAN_RECURRENCE")
                if peak.relative_prominence < params.minimum_relative_prominence:
                    peak_blocks.append("LOW_PROMINENCE")
                if peak.fwhm_mz is None or peak.fwhm_mz <= 0:
                    peak_blocks.append("INVALID_FWHM")
                interim.append(ReplicateRunPeak(
                    run_label=run_label, peak_id=f"T1REPLPEAK__{_safe_label(run_label)}__{serial:05d}",
                    apex_mz=peak.apex_mz, centroid_mz=peak.centroid_mz,
                    raw_apex_intensity=float(raw_mean[apex_index]), normalized_apex_intensity=peak.relative_apex_intensity,
                    raw_integrated_intensity=raw_area, normalized_integrated_intensity=peak.relative_integrated_intensity,
                    relative_intensity=peak.relative_apex_intensity, intensity_rank=0,
                    prominence=peak.prominence, relative_prominence=peak.relative_prominence,
                    fwhm=peak.fwhm_mz, left_bound_mz=peak.left_boundary_mz, right_bound_mz=peak.right_boundary_mz,
                    supporting_ms1_scan_count=support, total_ms1_scan_count=ms1_used,
                    scan_recurrence_fraction=recurrence,
                    first_supporting_scan_time=min(supporting_times) if supporting_times else None,
                    last_supporting_scan_time=max(supporting_times) if supporting_times else None,
                    detection_status=peak.peak_quality_class.value,
                    detection_block_reasons=_ordered_blocks(peak_blocks),
                ))
            ranked = sorted(interim, key=lambda item: (-item.normalized_apex_intensity, item.apex_mz, item.peak_id))
            rank = {peak.peak_id: index for index, peak in enumerate(ranked, 1)}
            peaks = tuple(replace(peak, intensity_rank=rank[peak.peak_id]) for peak in sorted(interim, key=lambda item: (item.apex_mz, item.peak_id)))
    return ReplicateRunPeakProfile(
        run_label=run_label, input_path=input_path, status=status,
        aggregation_method="BASE_PEAK_NORMALIZED_SCAN_BIN_MEAN",
        ms1_spectra_total=ms1_total, ms1_spectra_used=ms1_used, ms1_spectra_excluded=ms1_excluded,
        ms2_spectra_excluded=ms2_excluded, missing_ms_level_spectra=missing_level,
        mz_grid_method=f"ROUND_TO_FIXED_GRID_{params.mz_grid_step:.6f}_MZ",
        intensity_normalization_method="PER_SCAN_BASE_PEAK_THEN_RUN_MEAN",
        baseline_method="ZERO_FLOOR_NO_BASELINE_SUBTRACTION",
        smoothing_method="THREE_POINT_TRIANGULAR_0.25_0.50_0.25",
        peak_detection_method="detect_t1_profile_peaks",
        detected_peak_count=len(peaks), peaks=peaks, polarity_status=polarity,
        representation_status=representation, block_reasons=_ordered_blocks(blocks),
        comparison_mz_grid=grid.copy() if ms1_used and state["origin"] is not None else None,
        comparison_raw_profile=raw_mean.copy() if ms1_used and state["origin"] is not None else None,
        comparison_normalized_profile=smoothed.copy() if ms1_used and state["origin"] is not None else None,
    )


def build_ms1_peak_profile_for_run(
    mzml_path: Path, *, metadata_record: MzMLSourceMetadataRecord | None = None,
    detection_config: Mapping[str, Any] | None = None,
) -> ReplicateRunPeakProfile:
    path = Path(mzml_path)
    config = dict(detection_config or {})
    params = config.pop("parameters", None) or ReplicateAuditParameters(**config.pop("replicate_parameters", {}))
    peak_params = config.pop("peak_detection_parameters", None)
    if config:
        raise ValueError(f"unsupported detection_config keys: {sorted(config)}")
    run_label = metadata_record.technical_run_label if metadata_record and metadata_record.technical_run_label != "UNKNOWN" else path.stem
    if not path.exists():
        return ReplicateRunPeakProfile(run_label, str(path), "BLOCKED", "NOT_APPLICABLE", 0, 0, 0, 0, 0,
            "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", 0, (),
            "UNKNOWN", "UNKNOWN", ("INPUT_FILE_NOT_FOUND",))
    try:
        with mzml.MzML(str(path), decode_binary=False) as reader:
            return build_ms1_peak_profile_from_spectra(
                reader, run_label=run_label, input_path=str(path), metadata_record=metadata_record,
                parameters=params, peak_detection_parameters=peak_params,
            )
    except (OSError, IOError):
        return ReplicateRunPeakProfile(run_label, str(path), "BLOCKED", "NOT_APPLICABLE", 0, 0, 0, 0, 0,
            "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", 0, (),
            "UNKNOWN", "UNKNOWN", ("INPUT_FILE_UNREADABLE",))
    except Exception:
        return ReplicateRunPeakProfile(run_label, str(path), "BLOCKED", "NOT_APPLICABLE", 0, 0, 0, 0, 0,
            "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE", 0, (),
            "UNKNOWN", "UNKNOWN", ("PROFILE_EXTRACTION_FAILED",))


def _candidate_edges(left: Sequence[ReplicateRunPeak], right: Sequence[ReplicateRunPeak], parameters: ReplicateAuditParameters):
    right_mz = [peak.apex_mz for peak in right]
    edges: list[tuple[float, float, float, str, str, int, int]] = []
    left_candidates: dict[int, list[int]] = {index: [] for index in range(len(left))}
    right_candidates: dict[int, list[int]] = {index: [] for index in range(len(right))}
    maximum = max(parameters.absolute_tolerance_da, (max(right_mz, default=0.0) + parameters.absolute_tolerance_da) * parameters.ppm_tolerance * 1e-6)
    for i, peak in enumerate(left):
        lo = bisect_left(right_mz, peak.apex_mz - maximum)
        hi = bisect_right(right_mz, peak.apex_mz + maximum)
        for j in range(lo, hi):
            delta = right[j].apex_mz - peak.apex_mz
            tolerance = max(parameters.absolute_tolerance_da, peak.apex_mz * parameters.ppm_tolerance * 1e-6)
            if abs(delta) <= tolerance + 1e-12:
                ppm = abs(delta) / peak.apex_mz * 1e6
                edges.append((abs(delta), ppm, -min(peak.relative_intensity, right[j].relative_intensity), peak.peak_id, right[j].peak_id, i, j))
                left_candidates[i].append(j)
                right_candidates[j].append(i)
    return sorted(edges), left_candidates, right_candidates


def _ambiguity(left_count: int, right_count: int) -> MatchAmbiguityStatus:
    if left_count > 1 and right_count > 1:
        return MatchAmbiguityStatus.AMBIGUOUS_MANY_TO_MANY
    if left_count > 1:
        return MatchAmbiguityStatus.AMBIGUOUS_LEFT_TO_MULTIPLE
    if right_count > 1:
        return MatchAmbiguityStatus.AMBIGUOUS_RIGHT_TO_MULTIPLE
    return MatchAmbiguityStatus.UNAMBIGUOUS_ONE_TO_ONE


def _classify_match(left: ReplicateRunPeak, right: ReplicateRunPeak, ambiguity: MatchAmbiguityStatus, parameters: ReplicateAuditParameters):
    blocks: list[str] = []
    if ambiguity is not MatchAmbiguityStatus.UNAMBIGUOUS_ONE_TO_ONE:
        blocks.append("AMBIGUOUS_PEAK_MATCH")
        return ReplicateConsistencyStatus.REPRODUCED_AMBIGUOUS_MATCH, "LOW", blocks
    if left.detection_status in _LOW_QUALITY_CLASSES or right.detection_status in _LOW_QUALITY_CLASSES:
        return ReplicateConsistencyStatus.INSUFFICIENT_PEAK_QUALITY, "LOW", blocks
    if left.scan_recurrence_fraction < parameters.supportive_recurrence_fraction and right.scan_recurrence_fraction < parameters.supportive_recurrence_fraction:
        blocks.append("LOW_SCAN_RECURRENCE")
        return ReplicateConsistencyStatus.LOW_RECURRENCE_BOTH_RUNS, "LOW", blocks
    intensity_ratio = _ratio(right.normalized_apex_intensity, left.normalized_apex_intensity)
    log_ratio = abs(log2(intensity_ratio)) if intensity_ratio and intensity_ratio > 0 else None
    if log_ratio is None:
        blocks.append("MISSING_INTENSITY_NORMALIZATION")
    elif log_ratio > parameters.extreme_log2_intensity_ratio:
        blocks.append("EXTREME_INTENSITY_VARIATION")
        return ReplicateConsistencyStatus.REPRODUCED_WITH_INTENSITY_VARIATION, "MEDIUM", blocks
    if left.relative_prominence is None or right.relative_prominence is None:
        blocks.append("MISSING_PEAK_SHAPE_METRICS")
        return ReplicateConsistencyStatus.REPRODUCED_WITH_SHAPE_VARIATION, "LOW", blocks
    low_prominence = left.relative_prominence < parameters.minimum_relative_prominence or right.relative_prominence < parameters.minimum_relative_prominence
    if low_prominence:
        blocks.append("LOW_PROMINENCE")
    fwhm_ratio = _ratio(right.fwhm, left.fwhm)
    if fwhm_ratio is None:
        blocks.append("MISSING_PEAK_SHAPE_METRICS")
        return ReplicateConsistencyStatus.REPRODUCED_WITH_SHAPE_VARIATION, "LOW", blocks
    if not parameters.minimum_fwhm_ratio <= fwhm_ratio <= parameters.maximum_fwhm_ratio:
        blocks.append("INVALID_FWHM")
        return ReplicateConsistencyStatus.REPRODUCED_WITH_SHAPE_VARIATION, "MEDIUM", blocks
    delta = abs(right.apex_mz - left.apex_mz)
    strict = delta <= max(parameters.strict_absolute_tolerance_da, left.apex_mz * parameters.strict_ppm_tolerance * 1e-6) + 1e-12
    high = (
        strict and left.scan_recurrence_fraction >= parameters.high_recurrence_fraction
        and right.scan_recurrence_fraction >= parameters.high_recurrence_fraction
        and not low_prominence
        and abs(left.intensity_rank - right.intensity_rank) <= parameters.high_rank_difference
    )
    return (ReplicateConsistencyStatus.REPRODUCED_HIGH_CONFIDENCE, "HIGH", blocks) if high else (ReplicateConsistencyStatus.REPRODUCED_SUPPORTIVE, "MEDIUM", blocks)


def _match_record(left: ReplicateRunPeak | None, right: ReplicateRunPeak | None, *, run_1_label: str, run_2_label: str,
                  left_count: int, right_count: int, ambiguity: MatchAmbiguityStatus,
                  alternatives: tuple[str, ...], parameters: ReplicateAuditParameters) -> ReplicatePeakMatch:
    if left is None:
        status, confidence, blocks = ReplicateConsistencyStatus.RUN_2_ONLY, "LOW", []
    elif right is None:
        status, confidence, blocks = ReplicateConsistencyStatus.RUN_1_ONLY, "LOW", []
    else:
        status, confidence, blocks = _classify_match(left, right, ambiguity, parameters)
    delta = right.apex_mz - left.apex_mz if left and right else None
    intensity_ratio = _ratio(right.normalized_apex_intensity, left.normalized_apex_intensity) if left and right else None
    fwhm_ratio = _ratio(right.fwhm, left.fwhm) if left and right else None
    return ReplicatePeakMatch(
        replicate_pair_id=f"T1REPLPAIR__{left.peak_id if left else 'NONE'}__{right.peak_id if right else 'NONE'}",
        run_1_label=run_1_label, run_2_label=run_2_label,
        run_1_peak_id=left.peak_id if left else None, run_2_peak_id=right.peak_id if right else None,
        run_1_apex_mz=left.apex_mz if left else None, run_2_apex_mz=right.apex_mz if right else None,
        delta_mz=delta, absolute_delta_mz=abs(delta) if delta is not None else None,
        ppm_error=delta / left.apex_mz * 1e6 if delta is not None and left else None,
        run_1_centroid_mz=left.centroid_mz if left else None, run_2_centroid_mz=right.centroid_mz if right else None,
        centroid_delta_mz=right.centroid_mz - left.centroid_mz if left and right and left.centroid_mz is not None and right.centroid_mz is not None else None,
        run_1_normalized_intensity=left.normalized_apex_intensity if left else None,
        run_2_normalized_intensity=right.normalized_apex_intensity if right else None,
        normalized_intensity_ratio=intensity_ratio,
        log2_normalized_intensity_ratio=log2(intensity_ratio) if intensity_ratio and intensity_ratio > 0 else None,
        run_1_intensity_rank=left.intensity_rank if left else None, run_2_intensity_rank=right.intensity_rank if right else None,
        intensity_rank_difference=abs(left.intensity_rank - right.intensity_rank) if left and right else None,
        run_1_prominence=left.prominence if left else None, run_2_prominence=right.prominence if right else None,
        prominence_ratio=_ratio(right.prominence, left.prominence) if left and right else None,
        run_1_fwhm=left.fwhm if left else None, run_2_fwhm=right.fwhm if right else None, fwhm_ratio=fwhm_ratio,
        run_1_scan_recurrence_fraction=left.scan_recurrence_fraction if left else None,
        run_2_scan_recurrence_fraction=right.scan_recurrence_fraction if right else None,
        candidate_match_count_left=left_count, candidate_match_count_right=right_count,
        match_ambiguity_status=ambiguity, alternative_match_peak_ids=alternatives,
        replicate_consistency_status=status, replicate_consistency_confidence=confidence,
        consistency_block_reasons=_ordered_blocks(blocks),
    )


def match_replicate_peaks(
    left: ReplicateRunPeakProfile, right: ReplicateRunPeakProfile, *,
    absolute_tolerance_da: float | None = None, ppm_tolerance: float | None = None,
    parameters: ReplicateAuditParameters | None = None,
) -> tuple[ReplicatePeakMatch, ...]:
    params = parameters or ReplicateAuditParameters()
    if absolute_tolerance_da is not None:
        params = replace(params, absolute_tolerance_da=float(absolute_tolerance_da))
    if ppm_tolerance is not None:
        params = replace(params, ppm_tolerance=float(ppm_tolerance))
    params.validate()
    if (right.run_label, right.input_path) < (left.run_label, left.input_path):
        left, right = right, left
    a = tuple(sorted(left.peaks, key=lambda peak: (peak.apex_mz, peak.peak_id)))
    b = tuple(sorted(right.peaks, key=lambda peak: (peak.apex_mz, peak.peak_id)))
    edges, left_candidates, right_candidates = _candidate_edges(a, b, params)
    used_left: set[int] = set()
    used_right: set[int] = set()
    assigned: list[tuple[int, int]] = []
    for _, _, _, _, _, i, j in edges:
        if i not in used_left and j not in used_right:
            used_left.add(i)
            used_right.add(j)
            assigned.append((i, j))
    output: list[ReplicatePeakMatch] = []
    for i, j in sorted(assigned, key=lambda pair: (a[pair[0]].apex_mz, a[pair[0]].peak_id, b[pair[1]].peak_id)):
        ambiguity = _ambiguity(len(left_candidates[i]), len(right_candidates[j]))
        alternatives = tuple(sorted(
            {b[k].peak_id for k in left_candidates[i] if k != j}
            | {a[k].peak_id for k in right_candidates[j] if k != i}
        ))
        output.append(_match_record(a[i], b[j], run_1_label=left.run_label, run_2_label=right.run_label,
            left_count=len(left_candidates[i]), right_count=len(right_candidates[j]), ambiguity=ambiguity,
            alternatives=alternatives, parameters=params))
    for i, peak in enumerate(a):
        if i not in used_left:
            output.append(_match_record(peak, None, run_1_label=left.run_label, run_2_label=right.run_label,
                left_count=len(left_candidates[i]), right_count=0, ambiguity=MatchAmbiguityStatus.UNMATCHED,
                alternatives=tuple(sorted(b[k].peak_id for k in left_candidates[i])), parameters=params))
    for j, peak in enumerate(b):
        if j not in used_right:
            output.append(_match_record(None, peak, run_1_label=left.run_label, run_2_label=right.run_label,
                left_count=0, right_count=len(right_candidates[j]), ambiguity=MatchAmbiguityStatus.UNMATCHED,
                alternatives=tuple(sorted(a[k].peak_id for k in right_candidates[j])), parameters=params))
    return tuple(sorted(output, key=lambda item: (
        item.run_1_apex_mz if item.run_1_apex_mz is not None else float("inf"),
        item.run_2_apex_mz if item.run_2_apex_mz is not None else float("inf"), item.replicate_pair_id,
    )))


def _drift_summary(matches: Sequence[ReplicatePeakMatch], parameters: ReplicateAuditParameters) -> MZDriftSummary:
    usable = [match for match in matches if match.delta_mz is not None and match.ppm_error is not None
              and match.match_ambiguity_status is MatchAmbiguityStatus.UNAMBIGUOUS_ONE_TO_ONE
              and match.replicate_consistency_status is not ReplicateConsistencyStatus.INSUFFICIENT_PEAK_QUALITY]
    if len(usable) < parameters.minimum_drift_matches:
        return MZDriftSummary(None, None, None, None, len(usable), SystematicMZShiftStatus.INSUFFICIENT_MATCHES)
    delta = np.asarray([match.delta_mz for match in usable], dtype=float)
    ppm = np.asarray([match.ppm_error for match in usable], dtype=float)
    med_delta = float(np.median(delta))
    med_ppm = float(np.median(ppm))
    mad_delta = float(np.median(np.abs(delta - med_delta)))
    mad_ppm = float(np.median(np.abs(ppm - med_ppm)))
    absolute = abs(med_delta)
    status = SystematicMZShiftStatus.NO_MEANINGFUL_SHIFT if absolute <= parameters.no_shift_delta_da else SystematicMZShiftStatus.SMALL_SYSTEMATIC_SHIFT if absolute <= parameters.small_shift_delta_da else SystematicMZShiftStatus.POSSIBLE_SYSTEMATIC_SHIFT
    return MZDriftSummary(med_delta, med_ppm, mad_delta, mad_ppm, len(usable), status)


def summarize_replicate_consistency(
    left: ReplicateRunPeakProfile, right: ReplicateRunPeakProfile,
    matches: Sequence[ReplicatePeakMatch], *, parameters: ReplicateAuditParameters | None = None,
) -> ReplicateConsistencySummary:
    params = parameters or ReplicateAuditParameters()
    matched = [match for match in matches if match.run_1_peak_id and match.run_2_peak_id]
    run_1_only = [match for match in matches if match.run_1_peak_id and not match.run_2_peak_id]
    run_2_only = [match for match in matches if match.run_2_peak_id and not match.run_1_peak_id]
    ambiguous = [match for match in matched if match.match_ambiguity_status is not MatchAmbiguityStatus.UNAMBIGUOUS_ONE_TO_ONE]
    high = [match for match in matched if match.replicate_consistency_status is ReplicateConsistencyStatus.REPRODUCED_HIGH_CONFIDENCE]
    supportive_statuses = {ReplicateConsistencyStatus.REPRODUCED_SUPPORTIVE}
    supportive = [match for match in matched if match.replicate_consistency_status in supportive_statuses]
    low = [match for match in matched if match not in high and match not in supportive]
    union = left.detected_peak_count + right.detected_peak_count - len(matched)
    jaccard = len(matched) / union if union else 0.0
    overlap_left = len(matched) / left.detected_peak_count if left.detected_peak_count else 0.0
    overlap_right = len(matched) / right.detected_peak_count if right.detected_peak_count else 0.0
    blocks: list[str] = []
    if not left.peaks:
        blocks.append("NO_DETECTED_PEAKS_RUN_1")
    if not right.peaks:
        blocks.append("NO_DETECTED_PEAKS_RUN_2")
    if not matched:
        blocks.append("NO_MATCHED_PEAKS")
    if ambiguous:
        blocks.append("AMBIGUOUS_PEAK_MATCH")
    drift = _drift_summary(matches, params)
    if drift.systematic_mz_shift_status is SystematicMZShiftStatus.INSUFFICIENT_MATCHES:
        blocks.append("INSUFFICIENT_MATCHES_FOR_DRIFT")
    input_blockers = set(left.block_reasons + right.block_reasons) & {
        "INPUT_FILE_NOT_FOUND", "INPUT_FILE_UNREADABLE", "NO_MS1_SPECTRA", "PROFILE_EXTRACTION_FAILED",
        "PEAK_DETECTION_FAILED", "POLARITY_MISMATCH_BETWEEN_REPLICATES",
        "REPRESENTATION_MISMATCH_BETWEEN_REPLICATES", "REPLICATE_CONTEXT_MISMATCH",
    }
    if input_blockers:
        overall, confidence = OverallReplicateStatus.BLOCKED_BY_INPUT_QUALITY, "BLOCKED"
    elif not matched:
        overall, confidence = OverallReplicateStatus.INSUFFICIENT_COMPARABLE_PEAKS, "LOW"
    elif jaccard >= params.strong_jaccard and len(high) >= len(matched) / 2:
        overall, confidence = OverallReplicateStatus.STRONG_TECHNICAL_REPRODUCIBILITY, "HIGH"
    elif jaccard >= params.moderate_jaccard:
        overall, confidence = OverallReplicateStatus.MODERATE_TECHNICAL_REPRODUCIBILITY, "MEDIUM"
    elif min(overlap_left, overlap_right) < 0.20:
        overall, confidence = OverallReplicateStatus.SUBSTANTIAL_RUN_SPECIFIC_SIGNAL, "LOW"
    else:
        overall, confidence = OverallReplicateStatus.LIMITED_TECHNICAL_REPRODUCIBILITY, "LOW"
    return ReplicateConsistencySummary(
        run_1_label=left.run_label, run_2_label=right.run_label,
        run_1_detected_peak_count=left.detected_peak_count, run_2_detected_peak_count=right.detected_peak_count,
        matched_peak_pair_count=len(matched), run_1_only_peak_count=len(run_1_only), run_2_only_peak_count=len(run_2_only),
        ambiguous_match_count=len(ambiguous), high_confidence_reproduced_count=len(high),
        supportive_reproduced_count=len(supportive), low_quality_reproduced_count=len(low),
        union_peak_count=union, intersection_peak_count=len(matched), jaccard_index=jaccard,
        run_1_overlap_fraction=overlap_left, run_2_overlap_fraction=overlap_right,
        median_absolute_delta_mz=_median(match.absolute_delta_mz for match in matched),
        median_ppm_error=_median(match.ppm_error for match in matched),
        median_intensity_rank_difference=_median(match.intensity_rank_difference for match in matched),
        median_log2_normalized_intensity_ratio=_median(match.log2_normalized_intensity_ratio for match in matched),
        drift=drift, replicate_consistency_overall_status=overall, overall_confidence=confidence,
        overall_block_reasons=_ordered_blocks(blocks + list(input_blockers)),
    )


def _validate_pair_context(left: MzMLSourceMetadataRecord | None, right: MzMLSourceMetadataRecord | None) -> tuple[str, ...]:
    blocks: list[str] = []
    if left is None or right is None:
        return ("SOURCE_METADATA_RECORD_MISSING",)
    if left.context_source != "USER_PROVIDED_RUNTIME_MANIFEST" or right.context_source != "USER_PROVIDED_RUNTIME_MANIFEST":
        blocks.append("USER_MANIFEST_CONTEXT_MISSING")
    if (left.rna_identity, left.digest_type) != (right.rna_identity, right.digest_type):
        blocks.append("REPLICATE_CONTEXT_MISMATCH")
    if left.polarity_status is not right.polarity_status:
        blocks.append("POLARITY_MISMATCH_BETWEEN_REPLICATES")
    if left.representation_status is not right.representation_status:
        blocks.append("REPRESENTATION_MISMATCH_BETWEEN_REPLICATES")
    return _ordered_blocks(blocks)


def audit_sciex_t1_replicate_consistency(
    replicate_paths: Sequence[Path], *,
    source_metadata_records: Sequence[MzMLSourceMetadataRecord] | None = None,
    runtime_context: Mapping[str, Any] | None = None,
    parameters: ReplicateAuditParameters | None = None,
    peak_detection_parameters: T1PeakDetectionParameters | None = None,
) -> ReplicateConsistencyAuditResult:
    params = parameters or ReplicateAuditParameters()
    params.validate()
    paths = tuple(Path(path) for path in replicate_paths)
    metadata_by_path: dict[str, MzMLSourceMetadataRecord] = {}
    for record in source_metadata_records or ():
        metadata_by_path[record.input_path] = record
        metadata_by_path[record.file_name] = record
    context = dict(runtime_context or {})
    profiles: list[ReplicateRunPeakProfile] = []
    records: list[MzMLSourceMetadataRecord | None] = []
    for path in paths:
        record = metadata_by_path.get(str(path)) or metadata_by_path.get(path.name)
        records.append(record)
        label = context.get(path.name) or context.get(str(path))
        if label and record is not None:
            record = replace(record, technical_run_label=str(label))
        profile = build_ms1_peak_profile_for_run(path, metadata_record=record, detection_config={
            "parameters": params, "peak_detection_parameters": peak_detection_parameters,
        })
        profiles.append(profile)
    ordered = tuple(sorted(profiles, key=lambda item: (item.run_label, item.input_path)))
    record_by_label = {profile.run_label: records[index] for index, profile in enumerate(profiles)}
    matches: list[ReplicatePeakMatch] = []
    summaries: list[ReplicateConsistencySummary] = []
    for left_index in range(len(ordered)):
        for right_index in range(left_index + 1, len(ordered)):
            left, right = ordered[left_index], ordered[right_index]
            pair_blocks = _validate_pair_context(record_by_label.get(left.run_label), record_by_label.get(right.run_label))
            if pair_blocks:
                left = replace(left, block_reasons=_ordered_blocks(left.block_reasons + pair_blocks))
                right = replace(right, block_reasons=_ordered_blocks(right.block_reasons + pair_blocks))
            pair_matches = match_replicate_peaks(left, right, parameters=params)
            matches.extend(pair_matches)
            summaries.append(summarize_replicate_consistency(left, right, pair_matches, parameters=params))
    return ReplicateConsistencyAuditResult(params, ordered, tuple(matches), tuple(summaries))


def _scalarize(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return ";".join(str(_scalarize(item)) for item in value)
    if isinstance(value, dict):
        return {key: _scalarize(item) for key, item in value.items()}
    return value


def audit_optional_result(result: ReplicateConsistencyAuditResult) -> dict[str, Any]:
    run_peak_records = [_scalarize(asdict(peak)) for profile in result.run_profiles for peak in profile.peaks]
    match_records = [_scalarize(asdict(match)) for match in result.matches]
    summary_records = []
    for summary in result.summaries:
        row = asdict(summary)
        drift = row.pop("drift")
        row.update({f"drift_{key}": value for key, value in drift.items()})
        summary_records.append(_scalarize(row))
    return {"run_peak_records": run_peak_records, "match_records": match_records, "summary_records": summary_records}
