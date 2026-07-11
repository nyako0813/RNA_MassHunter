from statistics import pstdev
from typing import Any

from rna_masshunter.masses import neutral_mass_from_mz
from rna_masshunter.models import IntactMassCandidate, PeakTierResult


def _confidence(charge_state_count: int, min_charge_states: int) -> str:
    if charge_state_count >= min_charge_states + 2:
        return "High"
    if charge_state_count >= min_charge_states:
        return "Medium"
    return "Low"


DEFAULT_INTACT_QC_CONFIG = {
    "min_charge_states_for_reliable": 3,
    "min_charge_states_for_review": 2,
    "require_contiguous_charge_states": True,
    "max_neutral_mass_sd_da": 0.5,
    "max_neutral_mass_range_da": 1.5,
    "max_mass_error_ppm": 20,
    "max_envelope_internal_error_ppm": 20,
    "min_relative_intensity_percent": 0.5,
    "min_relative_envelope_intensity_percent_for_reliable": 1.0,
    "min_relative_envelope_intensity_percent_for_review": 0.1,
    "max_competing_envelopes": 3,
    "comparison_ready_statuses": ["Reliable", "Review"],
    "max_rt_range_min_for_reliable": 0.15,
    "max_rt_range_min_for_review": 0.30,
    "allow_trace_only_reliable": False,
    "search_mode": "untargeted",
    "reference_masses": [],
    "reference_mass_tolerance_ppm": 20,
    "neutral_mass_range": {"enabled": True, "min_da": 20000, "max_da": 30000},
    "target_review_mass_range": {"enabled": False, "min_da": None, "max_da": None},
}

QC_COLUMNS = [
    "Cluster_ID",
    "Reconstructed_Mass",
    "Observed_Mass",
    "In_Neutral_Mass_Search_Range",
    "Neutral_Mass_Search_Min_Da",
    "Neutral_Mass_Search_Max_Da",
    "Neutral_Mass_Range_Status",
    "In_Target_Review_Mass_Range",
    "Target_Review_Mass_Range_Status",
    "Target_Review_Priority",
    "Envelope_QC_Eligible",
    "Intact_Review_Eligible",
    "Intact_Strict_Eligible",
    "Intact_Envelope_QC_Score",
    "Intact_Envelope_QC_Rank",
    "Strict_Eligible_Rank",
    "Review_Eligible_Rank",
    "Dominant_Intact_Envelope_Flag",
    "Reconstruction_Status",
    "Reconstruction_Confidence",
    "Comparison_Ready_Strict",
    "Comparison_Ready_Review",
    "Comparison_Ready",
    "Comparison_Readiness_Reason",
    "Total_Supporting_Intensity",
    "Mean_Supporting_Intensity",
    "Max_Supporting_Intensity",
    "Relative_Envelope_Intensity_Percent",
    "Relative_Overall_Envelope_Intensity_Percent",
    "Relative_In_Range_Raw_Intensity_Percent",
    "Relative_Intact_Eligible_Intensity_Percent",
    "Supporting_Peak_Classes",
    "Trace_Only_Envelope",
    "Num_Supporting_Charge_States",
    "Charge_State_Range",
    "Charge_State_Continuity",
    "RT_Min",
    "RT_Max",
    "RT_Mean",
    "RT_Range_Min",
    "Max_RT_Difference_Min",
    "RT_Consistency",
    "Neutral_Mass_SD",
    "Neutral_Mass_Range",
    "Envelope_Internal_Error_ppm",
    "Max_Mass_Error_ppm",
    "Unmodified_Theory_Delta_Da",
    "Unmodified_Theory_Delta_ppm",
    "Best_Reference_Label",
    "Best_Reference_Mass_Da",
    "Reference_Mass_Error_Da",
    "Reference_Mass_Error_ppm",
    "Reference_Mass_Matched",
    "Competing_Envelope_Count",
    "Limiting_Factors",
    "Severe_Limiting_Factors",
    "Num_Limiting_Factors",
    "Primary_Limiting_Factor",
]

DIAGNOSTIC_COLUMNS = [
    "Total_Reconstruction_Candidates",
    "Reliable_Count",
    "Review_Count",
    "Insufficient_Count",
    "Failed_Count",
    "Envelope_QC_Eligible_Count",
    "Intact_Strict_Eligible_Count",
    "Intact_Review_Eligible_Count",
    "Comparison_Ready_Strict_Count",
    "Comparison_Ready_Review_Count",
    "Comparison_Ready_Count",
    "Trace_Only_Envelope_Count",
    "Noncontiguous_Envelope_Count",
    "RT_Inconsistent_Count",
    "Internal_Mass_Error_Count",
    "Theory_Near_Match_Count",
    "Reference_Match_Count",
    "Dominant_Envelope_Mass",
    "Dominant_Envelope_Intensity",
    "Dominant_Envelope_Status",
    "Dominant_Envelope_Comparison_Ready",
    "Dominant_Envelope_Overall_Mass",
    "Dominant_Envelope_Overall_Intensity",
    "Dominant_Envelope_In_Mass_Range_Mass",
    "Dominant_Envelope_In_Mass_Range_Intensity",
    "Dominant_Envelope_In_Mass_Range_Status",
    "Dominant_Envelope_In_Mass_Range_Comparison_Ready",
    "Dominant_Envelope_In_Search_Range_Raw_Mass",
    "Dominant_Envelope_In_Search_Range_Raw_Intensity",
    "Dominant_Intact_Strict_Envelope_Mass",
    "Dominant_Intact_Strict_Envelope_Intensity",
    "Dominant_Intact_Strict_QC_Score",
    "Dominant_Intact_Review_Envelope_Mass",
    "Dominant_Intact_Review_Envelope_Intensity",
    "Dominant_Intact_Review_QC_Score",
    "Dominant_Intact_Eligible_Envelope_Mass",
    "Dominant_Intact_Eligible_Envelope_Intensity",
    "Dominant_Intact_Eligible_QC_Score",
    "Dominant_Intact_Eligible_Reference_Label",
    "Failure_Reason_Counts",
    "Reconstruction_Enabled",
    "Neutral_Mass_Search_Min_Da",
    "Neutral_Mass_Search_Max_Da",
    "Total_Candidates_Before_Mass_Range_Filter",
    "Total_Candidates_In_Mass_Range",
    "Total_Candidates_Outside_Mass_Range",
    "Target_Review_Mass_Range_Settings",
    "Target_Review_Candidate_Count",
    "Search_Mode",
    "Intensity_Normalization_Method",
    "RT_Tolerance_Settings",
    "Reference_Masses_Used",
    "Min_Charge_States_For_Reliable",
    "Min_Charge_States_For_Review",
    "Require_Contiguous_Charge_States",
    "Max_Neutral_Mass_SD_Da",
    "Max_Neutral_Mass_Range_Da",
    "Max_Envelope_Internal_Error_ppm",
    "Min_Relative_Intensity_Percent",
    "Min_Relative_Envelope_Intensity_Percent_For_Reliable",
    "Min_Relative_Envelope_Intensity_Percent_For_Review",
    "Max_Competing_Envelopes",
    "Comparison_Ready_Statuses",
    "Notes",
]

SEVERE_LIMITING_FACTORS = {
    "reconstruction_disabled",
    "no_charge_state_candidates",
    "insufficient_charge_states",
    "mass_spread_too_large",
    "internal_mass_error_too_large",
    "rt_inconsistent",
    "insufficient_intensity_support",
    "multiple_competing_envelopes",
}


def _as_bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _as_reference_masses(value: Any) -> list[dict[str, Any]]:
    references = []
    if not value:
        return references
    raw_items = value if isinstance(value, list) else [value]
    for index, item in enumerate(raw_items, start=1):
        if isinstance(item, dict):
            mass = item.get("mass_da") or item.get("mass") or item.get("Mass_Da")
            label = item.get("label") or item.get("name") or f"reference_{index}"
        else:
            mass = item
            label = f"reference_{index}"
        try:
            references.append({"label": str(label), "mass_da": float(mass)})
        except (TypeError, ValueError):
            continue
    return references


def _optional_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _qc_config(reconstruction_config: dict[str, Any]) -> dict[str, Any]:
    raw = reconstruction_config.get("intact_reconstruction") or reconstruction_config.get("qc") or {}
    merged = {**DEFAULT_INTACT_QC_CONFIG, **raw}
    merged["min_charge_states_for_reliable"] = int(merged.get("min_charge_states_for_reliable") or 3)
    merged["min_charge_states_for_review"] = int(merged.get("min_charge_states_for_review") or 2)
    merged["require_contiguous_charge_states"] = _as_bool(merged.get("require_contiguous_charge_states"), True)
    merged["max_neutral_mass_sd_da"] = float(merged.get("max_neutral_mass_sd_da") or 0.5)
    merged["max_neutral_mass_range_da"] = float(merged.get("max_neutral_mass_range_da") or 1.5)
    merged["max_mass_error_ppm"] = float(merged.get("max_mass_error_ppm") or 20)
    merged["max_envelope_internal_error_ppm"] = float(
        merged.get("max_envelope_internal_error_ppm") or merged["max_mass_error_ppm"]
    )
    merged["min_relative_intensity_percent"] = float(merged.get("min_relative_intensity_percent") or 0.5)
    merged["min_relative_envelope_intensity_percent_for_reliable"] = float(
        merged.get("min_relative_envelope_intensity_percent_for_reliable") or 1.0
    )
    merged["min_relative_envelope_intensity_percent_for_review"] = float(
        merged.get("min_relative_envelope_intensity_percent_for_review") or 0.1
    )
    merged["max_competing_envelopes"] = int(merged.get("max_competing_envelopes") or 3)
    merged["max_rt_range_min_for_reliable"] = float(merged.get("max_rt_range_min_for_reliable") or 0.15)
    merged["max_rt_range_min_for_review"] = float(merged.get("max_rt_range_min_for_review") or 0.30)
    merged["allow_trace_only_reliable"] = _as_bool(merged.get("allow_trace_only_reliable"), False)
    search_mode = str(merged.get("search_mode") or "untargeted")
    merged["search_mode"] = search_mode if search_mode in {"untargeted", "theoretical_targeted", "reference_targeted"} else "untargeted"
    statuses = merged.get("comparison_ready_statuses") or ["Reliable", "Review"]
    if isinstance(statuses, str):
        statuses = [item.strip() for item in statuses.split(",") if item.strip()]
    merged["comparison_ready_statuses"] = list(statuses)
    merged["reference_masses"] = _as_reference_masses(merged.get("reference_masses"))
    merged["reference_mass_tolerance_ppm"] = float(merged.get("reference_mass_tolerance_ppm") or 20)
    neutral_range = merged.get("neutral_mass_range") or {}
    if not isinstance(neutral_range, dict):
        neutral_range = {}
    merged["neutral_mass_range"] = {
        "enabled": _as_bool(neutral_range.get("enabled"), True),
        "min_da": float(neutral_range.get("min_da", 20000) if neutral_range.get("min_da", None) is not None else 20000),
        "max_da": float(neutral_range.get("max_da", 30000) if neutral_range.get("max_da", None) is not None else 30000),
    }
    if merged["neutral_mass_range"]["min_da"] > merged["neutral_mass_range"]["max_da"]:
        merged["neutral_mass_range"]["min_da"], merged["neutral_mass_range"]["max_da"] = (
            merged["neutral_mass_range"]["max_da"],
            merged["neutral_mass_range"]["min_da"],
        )
    target_range = merged.get("target_review_mass_range") or {}
    if not isinstance(target_range, dict):
        target_range = {}
    merged["target_review_mass_range"] = {
        "enabled": _as_bool(target_range.get("enabled"), False),
        "min_da": _optional_float(target_range.get("min_da")),
        "max_da": _optional_float(target_range.get("max_da")),
    }
    target_min = merged["target_review_mass_range"]["min_da"]
    target_max = merged["target_review_mass_range"]["max_da"]
    if target_min is not None and target_max is not None and target_min > target_max:
        merged["target_review_mass_range"]["min_da"], merged["target_review_mass_range"]["max_da"] = (
            target_max,
            target_min,
        )
    return merged


def _charge_state_range(charges: list[int]) -> str:
    if not charges:
        return ""
    if len(charges) == 1:
        return str(charges[0])
    return f"{min(charges)}-{max(charges)}"


def _charge_continuity(charges: list[int]) -> str:
    if not charges:
        return "missing"
    expected = set(range(min(charges), max(charges) + 1))
    return "contiguous" if set(charges) == expected else "non_contiguous"


def _primary_factor(factors: list[str]) -> str:
    priority = [
        "reconstruction_disabled",
        "no_charge_state_candidates",
        "insufficient_charge_states",
        "internal_mass_error_too_large",
        "mass_spread_too_large",
        "rt_inconsistent",
        "trace_only_envelope",
        "insufficient_intensity_support",
        "non_contiguous_charge_states",
        "multiple_competing_envelopes",
        "rt_not_available",
    ]
    for item in priority:
        if item in factors:
            return item
    return factors[0] if factors else ""


def _ppm(delta: float | None, reference: float | None) -> float | None:
    if delta is None or not reference:
        return None
    return delta / reference * 1_000_000


def _max_abs_ppm(values: list[float], center: float | None) -> float | None:
    if not values or not center:
        return None
    return max(abs(value - center) / center * 1_000_000 for value in values)


def _rt_metrics(cluster_peaks: list[dict[str, Any]], qc_config: dict[str, Any]) -> tuple[float | None, float | None, float | None, float | None, str]:
    rts = []
    for row in cluster_peaks:
        value = row.get("RT")
        if value is None or value == "":
            continue
        try:
            rts.append(float(value))
        except (TypeError, ValueError):
            continue
    if not rts:
        return None, None, None, None, "not_available"
    rt_min = min(rts)
    rt_max = max(rts)
    rt_mean = sum(rts) / len(rts)
    rt_range = rt_max - rt_min
    if rt_range <= qc_config["max_rt_range_min_for_reliable"]:
        consistency = "consistent"
    elif rt_range <= qc_config["max_rt_range_min_for_review"]:
        consistency = "review"
    else:
        consistency = "inconsistent"
    return rt_min, rt_max, rt_mean, rt_range, consistency


def _neutral_mass_range_status(observed_mass: float, qc_config: dict[str, Any]) -> tuple[bool, float, float, str]:
    neutral_range = qc_config["neutral_mass_range"]
    min_da = float(neutral_range["min_da"])
    max_da = float(neutral_range["max_da"])
    if not neutral_range.get("enabled", True):
        return True, min_da, max_da, "not_applied"
    in_range = min_da <= observed_mass <= max_da
    return in_range, min_da, max_da, "in_range" if in_range else "outside_range"


def _target_review_range_status(observed_mass: float, qc_config: dict[str, Any]) -> tuple[bool, str, str]:
    target_range = qc_config["target_review_mass_range"]
    if not target_range.get("enabled", False):
        return False, "not_configured", "not_configured"
    if target_range.get("min_da") is None or target_range.get("max_da") is None:
        return False, "not_configured", "not_configured"
    min_da = float(target_range["min_da"])
    max_da = float(target_range["max_da"])
    in_range = min_da <= observed_mass <= max_da
    status = "in_range" if in_range else "outside_range"
    priority = "target_review" if in_range else "outside_target_review"
    return in_range, status, priority


def _target_review_settings(qc_config: dict[str, Any]) -> str:
    target_range = qc_config["target_review_mass_range"]
    if not target_range.get("enabled", False):
        return "disabled"
    if target_range.get("min_da") is None or target_range.get("max_da") is None:
        return "not_configured"
    return f"enabled:{target_range['min_da']}-{target_range['max_da']} Da"


def _rt_rank_score(value: str) -> int:
    return {"consistent": 2, "review": 1}.get(str(value or ""), 0)


def _safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _small_metric_score(value: Any, limit: float, points: float) -> float:
    if value is None or value == "":
        return points
    numeric = _safe_float(value, limit)
    if limit <= 0:
        return points
    return max(0.0, points * (1.0 - min(numeric / limit, 1.0)))


def _intact_qc_score(row: dict[str, Any], qc_config: dict[str, Any]) -> float:
    score = 0.0
    if row.get("Intact_Strict_Eligible"):
        score += 30.0
    elif row.get("Intact_Review_Eligible"):
        score += 22.0
    elif row.get("Envelope_QC_Eligible"):
        score += 12.0
    if row.get("Charge_State_Continuity") == "contiguous":
        score += 12.0
    score += min(_safe_float(row.get("Num_Supporting_Charge_States")), 6.0) * 4.0
    score += {"consistent": 12.0, "review": 6.0}.get(str(row.get("RT_Consistency") or ""), 0.0)
    score += _small_metric_score(row.get("Envelope_Internal_Error_ppm"), qc_config["max_envelope_internal_error_ppm"], 10.0)
    score += _small_metric_score(row.get("Neutral_Mass_SD"), qc_config["max_neutral_mass_sd_da"], 8.0)
    score += _small_metric_score(row.get("Neutral_Mass_Range"), qc_config["max_neutral_mass_range_da"], 8.0)
    score += min(_safe_float(row.get("Relative_Overall_Envelope_Intensity_Percent")), 100.0) / 100.0 * 8.0
    return round(score, 3)


def _dominant_sort_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        bool(row.get("Intact_Strict_Eligible")),
        bool(row.get("Intact_Review_Eligible")),
        row.get("Charge_State_Continuity") == "contiguous",
        _safe_float(row.get("Num_Supporting_Charge_States")),
        _rt_rank_score(str(row.get("RT_Consistency") or "")),
        -_safe_float(row.get("Envelope_Internal_Error_ppm"), 1_000_000.0),
        -_safe_float(row.get("Neutral_Mass_SD"), 1_000_000.0),
        -_safe_float(row.get("Neutral_Mass_Range"), 1_000_000.0),
        _safe_float(row.get("Total_Supporting_Intensity")),
    )


def _dominant_row(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return max(rows, key=_dominant_sort_key, default={})


def _reference_match(observed_mass: float, qc_config: dict[str, Any]) -> tuple[str, float | None, float | None, float | None, bool]:
    references = qc_config.get("reference_masses") or []
    if not references:
        return "not_configured", None, None, None, False
    best = None
    for reference in references:
        mass = float(reference["mass_da"])
        error_da = observed_mass - mass
        error_ppm = _ppm(error_da, mass)
        score = abs(error_ppm) if error_ppm is not None else float("inf")
        if best is None or score < best[0]:
            best = (score, str(reference["label"]), mass, error_da, error_ppm)
    if best is None:
        return "not_configured", None, None, None, False
    matched = best[4] is not None and abs(best[4]) <= qc_config["reference_mass_tolerance_ppm"]
    return best[1], best[2], best[3], best[4], matched


def _class_summary(cluster_peaks: list[dict[str, Any]]) -> tuple[str, bool]:
    classes = []
    for row in cluster_peaks:
        value = str(row.get("Peak_Tier") or "").strip()
        if value and value not in classes:
            classes.append(value)
    trace_only = bool(classes) and all(value.lower() == "trace" for value in classes)
    return "; ".join(classes), trace_only


def build_intact_reconstruction_qc(
    candidates: list[IntactMassCandidate],
    charge_state_peaks: list[dict[str, Any]],
    reconstruction_config: dict[str, Any],
    reconstruction_enabled: bool = True,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    qc_config = _qc_config(reconstruction_config or {})
    peaks_by_cluster: dict[str, list[dict[str, Any]]] = {}
    for row in charge_state_peaks or []:
        peaks_by_cluster.setdefault(str(row.get("Cluster_ID") or ""), []).append(row)

    max_intensity = max((float(getattr(candidate, "total_intensity", 0.0) or 0.0) for candidate in candidates), default=0.0)
    qc_rows: list[dict[str, Any]] = []
    for candidate in candidates:
        cluster_id = candidate.cluster_id or ""
        cluster_peaks = peaks_by_cluster.get(cluster_id, [])
        charges = sorted({int(charge) for charge in candidate.charge_states})
        neutral_masses = [float(row.get("Neutral_Mass")) for row in cluster_peaks if row.get("Neutral_Mass") is not None]
        if not neutral_masses and candidate.observed_mass is not None:
            neutral_masses = [float(candidate.observed_mass)]
        reconstructed_mass = float(candidate.observed_mass)
        neutral_sd = pstdev(neutral_masses) if len(neutral_masses) > 1 else 0.0 if neutral_masses else None
        neutral_range = (max(neutral_masses) - min(neutral_masses)) if neutral_masses else None
        envelope_internal_error_ppm = _max_abs_ppm(neutral_masses, reconstructed_mass)
        continuity = _charge_continuity(charges)
        rt_min, rt_max, rt_mean, rt_range, rt_consistency = _rt_metrics(cluster_peaks, qc_config)
        peak_classes, trace_only = _class_summary(cluster_peaks)
        intensities = []
        for row in cluster_peaks:
            try:
                intensities.append(float(row.get("Intensity") or 0.0))
            except (TypeError, ValueError):
                continue
        total_intensity = float(candidate.total_intensity or sum(intensities) or 0.0)
        mean_intensity = sum(intensities) / len(intensities) if intensities else total_intensity / max(len(charges), 1)
        max_supporting_intensity = max(intensities) if intensities else total_intensity
        relative_intensity = (total_intensity / max_intensity * 100.0) if max_intensity else 0.0
        competing = sum(
            1
            for other in candidates
            if other is not candidate
            and abs(float(other.observed_mass) - reconstructed_mass) <= qc_config["max_neutral_mass_range_da"]
        )
        unmodified_delta_da = candidate.mass_error_da
        unmodified_delta_ppm = candidate.mass_error_ppm
        reference_label, reference_mass, reference_error_da, reference_error_ppm, reference_matched = _reference_match(reconstructed_mass, qc_config)
        in_mass_range, neutral_search_min, neutral_search_max, neutral_range_status = _neutral_mass_range_status(reconstructed_mass, qc_config)
        in_target_range, target_range_status, target_review_priority = _target_review_range_status(reconstructed_mass, qc_config)

        factors: list[str] = []
        if len(charges) < qc_config["min_charge_states_for_review"]:
            factors.append("insufficient_charge_states")
        if qc_config["require_contiguous_charge_states"] and continuity == "non_contiguous":
            factors.append("non_contiguous_charge_states")
        if neutral_sd is not None and neutral_sd > qc_config["max_neutral_mass_sd_da"]:
            factors.append("mass_spread_too_large")
        if neutral_range is not None and neutral_range > qc_config["max_neutral_mass_range_da"]:
            factors.append("mass_spread_too_large")
        if envelope_internal_error_ppm is not None and envelope_internal_error_ppm > qc_config["max_envelope_internal_error_ppm"]:
            factors.append("internal_mass_error_too_large")
        if rt_consistency == "inconsistent":
            factors.append("rt_inconsistent")
        elif rt_consistency == "not_available":
            factors.append("rt_not_available")
        if relative_intensity < qc_config["min_relative_envelope_intensity_percent_for_review"]:
            factors.append("insufficient_intensity_support")
        if trace_only:
            factors.append("trace_only_envelope")
        if competing > qc_config["max_competing_envelopes"]:
            factors.append("multiple_competing_envelopes")
        if not in_mass_range:
            factors.append("outside_neutral_mass_search_range")
        factors = list(dict.fromkeys(factors))
        severe_factors = [factor for factor in factors if factor in SEVERE_LIMITING_FACTORS]

        basic_internal_ok = (
            (neutral_sd is None or neutral_sd <= qc_config["max_neutral_mass_sd_da"])
            and (neutral_range is None or neutral_range <= qc_config["max_neutral_mass_range_da"])
            and (envelope_internal_error_ppm is None or envelope_internal_error_ppm <= qc_config["max_envelope_internal_error_ppm"])
        )
        contiguous_ok = continuity == "contiguous" or not qc_config["require_contiguous_charge_states"]
        reliable_intensity_ok = relative_intensity >= qc_config["min_relative_envelope_intensity_percent_for_reliable"]
        review_intensity_ok = relative_intensity >= qc_config["min_relative_envelope_intensity_percent_for_review"]
        reliable_rt_ok = rt_consistency == "consistent"
        review_rt_ok = rt_consistency in {"consistent", "review"}
        trace_ok_for_reliable = not trace_only or qc_config["allow_trace_only_reliable"]

        if not reconstruction_enabled:
            status = "Failed"
            factors.insert(0, "reconstruction_disabled")
            severe_factors.insert(0, "reconstruction_disabled")
        elif len(charges) < qc_config["min_charge_states_for_review"]:
            status = "Insufficient"
        elif (
            len(charges) >= qc_config["min_charge_states_for_reliable"]
            and contiguous_ok
            and basic_internal_ok
            and reliable_rt_ok
            and reliable_intensity_ok
            and trace_ok_for_reliable
            and competing <= qc_config["max_competing_envelopes"]
        ):
            status = "Reliable"
        else:
            status = "Review"
        factors = list(dict.fromkeys(factors))
        severe_factors = list(dict.fromkeys(severe_factors))
        confidence = {"Reliable": "High", "Review": "Medium", "Insufficient": "Low", "Failed": "None"}.get(status, "Low")
        envelope_qc_eligible = (
            in_mass_range
            and reconstruction_enabled
            and len(charges) >= qc_config["min_charge_states_for_review"]
            and basic_internal_ok
            and review_rt_ok
            and not severe_factors
        )
        intact_strict_eligible = (
            envelope_qc_eligible
            and status == "Reliable"
            and contiguous_ok
            and reliable_rt_ok
            and reliable_intensity_ok
            and trace_ok_for_reliable
        )
        intact_review_eligible = (
            envelope_qc_eligible
            and not intact_strict_eligible
            and "Review" in qc_config["comparison_ready_statuses"]
            and review_intensity_ok
        )
        comparison_ready_strict = intact_strict_eligible
        comparison_ready_review = intact_review_eligible
        comparison_ready = comparison_ready_strict or comparison_ready_review
        readiness_reason = "strict" if comparison_ready_strict else "review" if comparison_ready_review else _primary_factor(factors) or "not_ready"
        primary_factor = _primary_factor(factors)

        candidate.reconstruction_status = status
        candidate.reconstruction_confidence = confidence
        candidate.num_supporting_charge_states = len(charges)
        candidate.charge_state_range = _charge_state_range(charges)
        candidate.charge_state_continuity = continuity
        candidate.neutral_mass_sd = neutral_sd
        candidate.neutral_mass_range = neutral_range
        candidate.envelope_internal_error_ppm = envelope_internal_error_ppm
        candidate.max_mass_error_ppm = envelope_internal_error_ppm
        candidate.unmodified_theory_delta_da = unmodified_delta_da
        candidate.unmodified_theory_delta_ppm = unmodified_delta_ppm
        candidate.best_reference_label = reference_label
        candidate.best_reference_mass_da = reference_mass
        candidate.reference_mass_error_da = reference_error_da
        candidate.reference_mass_error_ppm = reference_error_ppm
        candidate.reference_mass_matched = reference_matched
        candidate.in_neutral_mass_search_range = in_mass_range
        candidate.neutral_mass_search_min_da = neutral_search_min
        candidate.neutral_mass_search_max_da = neutral_search_max
        candidate.neutral_mass_range_status = neutral_range_status
        candidate.in_target_review_mass_range = in_target_range
        candidate.target_review_mass_range_status = target_range_status
        candidate.target_review_priority = target_review_priority
        candidate.envelope_qc_eligible = envelope_qc_eligible
        candidate.intact_review_eligible = intact_review_eligible
        candidate.intact_strict_eligible = intact_strict_eligible
        candidate.rt_min = rt_min
        candidate.rt_max = rt_max
        candidate.rt_mean = rt_mean
        candidate.rt_range_min = rt_range
        candidate.max_rt_difference_min = rt_range
        candidate.rt_consistency = rt_consistency
        candidate.total_supporting_intensity = total_intensity
        candidate.mean_supporting_intensity = mean_intensity
        candidate.max_supporting_intensity = max_supporting_intensity
        candidate.relative_envelope_intensity_percent = relative_intensity
        candidate.supporting_peak_classes = peak_classes
        candidate.trace_only_envelope = trace_only
        candidate.competing_envelope_count = competing
        candidate.limiting_factors = "; ".join(factors)
        candidate.severe_limiting_factors = "; ".join(severe_factors)
        candidate.num_limiting_factors = len(factors)
        candidate.primary_limiting_factor = primary_factor
        candidate.comparison_ready_strict = comparison_ready_strict
        candidate.comparison_ready_review = comparison_ready_review
        candidate.comparison_ready = comparison_ready
        candidate.comparison_readiness_reason = readiness_reason

        qc_rows.append({
            "Cluster_ID": cluster_id,
            "Reconstructed_Mass": reconstructed_mass,
            "Observed_Mass": reconstructed_mass,
            "In_Neutral_Mass_Search_Range": in_mass_range,
            "Neutral_Mass_Search_Min_Da": neutral_search_min,
            "Neutral_Mass_Search_Max_Da": neutral_search_max,
            "Neutral_Mass_Range_Status": neutral_range_status,
            "In_Target_Review_Mass_Range": in_target_range,
            "Target_Review_Mass_Range_Status": target_range_status,
            "Target_Review_Priority": target_review_priority,
            "Envelope_QC_Eligible": envelope_qc_eligible,
            "Intact_Review_Eligible": intact_review_eligible,
            "Intact_Strict_Eligible": intact_strict_eligible,
            "Intact_Envelope_QC_Score": None,
            "Intact_Envelope_QC_Rank": None,
            "Strict_Eligible_Rank": None,
            "Review_Eligible_Rank": None,
            "Dominant_Intact_Envelope_Flag": False,
            "Reconstruction_Status": status,
            "Reconstruction_Confidence": confidence,
            "Comparison_Ready_Strict": comparison_ready_strict,
            "Comparison_Ready_Review": comparison_ready_review,
            "Comparison_Ready": comparison_ready,
            "Comparison_Readiness_Reason": readiness_reason,
            "Total_Supporting_Intensity": total_intensity,
            "Mean_Supporting_Intensity": mean_intensity,
            "Max_Supporting_Intensity": max_supporting_intensity,
            "Relative_Envelope_Intensity_Percent": relative_intensity,
            "Relative_Overall_Envelope_Intensity_Percent": relative_intensity,
            "Relative_In_Range_Raw_Intensity_Percent": None,
            "Relative_Intact_Eligible_Intensity_Percent": None,
            "Supporting_Peak_Classes": peak_classes,
            "Trace_Only_Envelope": trace_only,
            "Num_Supporting_Charge_States": len(charges),
            "Charge_State_Range": candidate.charge_state_range,
            "Charge_State_Continuity": continuity,
            "RT_Min": rt_min,
            "RT_Max": rt_max,
            "RT_Mean": rt_mean,
            "RT_Range_Min": rt_range,
            "Max_RT_Difference_Min": rt_range,
            "RT_Consistency": rt_consistency,
            "Neutral_Mass_SD": neutral_sd,
            "Neutral_Mass_Range": neutral_range,
            "Envelope_Internal_Error_ppm": envelope_internal_error_ppm,
            "Max_Mass_Error_ppm": envelope_internal_error_ppm,
            "Unmodified_Theory_Delta_Da": unmodified_delta_da,
            "Unmodified_Theory_Delta_ppm": unmodified_delta_ppm,
            "Best_Reference_Label": reference_label,
            "Best_Reference_Mass_Da": reference_mass,
            "Reference_Mass_Error_Da": reference_error_da,
            "Reference_Mass_Error_ppm": reference_error_ppm,
            "Reference_Mass_Matched": reference_matched,
            "Competing_Envelope_Count": competing,
            "Limiting_Factors": candidate.limiting_factors,
            "Severe_Limiting_Factors": candidate.severe_limiting_factors,
            "Num_Limiting_Factors": len(factors),
            "Primary_Limiting_Factor": primary_factor,
        })

    max_in_range_intensity = max((_safe_float(row.get("Total_Supporting_Intensity")) for row in qc_rows if row.get("In_Neutral_Mass_Search_Range")), default=0.0)
    eligible_rows = [row for row in qc_rows if row.get("Intact_Strict_Eligible") or row.get("Intact_Review_Eligible")]
    max_eligible_intensity = max((_safe_float(row.get("Total_Supporting_Intensity")) for row in eligible_rows), default=0.0)
    for row in qc_rows:
        total = _safe_float(row.get("Total_Supporting_Intensity"))
        row["Relative_In_Range_Raw_Intensity_Percent"] = (total / max_in_range_intensity * 100.0) if max_in_range_intensity and row.get("In_Neutral_Mass_Search_Range") else 0.0
        row["Relative_Intact_Eligible_Intensity_Percent"] = (total / max_eligible_intensity * 100.0) if max_eligible_intensity and (row.get("Intact_Strict_Eligible") or row.get("Intact_Review_Eligible")) else 0.0
        row["Intact_Envelope_QC_Score"] = _intact_qc_score(row, qc_config)

    ranked_rows = sorted(qc_rows, key=_dominant_sort_key, reverse=True)
    for rank, row in enumerate(ranked_rows, start=1):
        row["Intact_Envelope_QC_Rank"] = rank
    strict_rows = sorted([row for row in qc_rows if row.get("Intact_Strict_Eligible")], key=_dominant_sort_key, reverse=True)
    for rank, row in enumerate(strict_rows, start=1):
        row["Strict_Eligible_Rank"] = rank
    review_rows = sorted([row for row in qc_rows if row.get("Intact_Review_Eligible")], key=_dominant_sort_key, reverse=True)
    for rank, row in enumerate(review_rows, start=1):
        row["Review_Eligible_Rank"] = rank
    dominant_eligible = strict_rows[0] if strict_rows else review_rows[0] if review_rows else None
    if dominant_eligible is not None:
        dominant_eligible["Dominant_Intact_Envelope_Flag"] = True

    candidates_by_cluster = {candidate.cluster_id or "": candidate for candidate in candidates}
    for row in qc_rows:
        candidate = candidates_by_cluster.get(str(row.get("Cluster_ID") or ""))
        if candidate is None:
            continue
        candidate.relative_overall_envelope_intensity_percent = row["Relative_Overall_Envelope_Intensity_Percent"]
        candidate.relative_in_range_raw_intensity_percent = row["Relative_In_Range_Raw_Intensity_Percent"]
        candidate.relative_intact_eligible_intensity_percent = row["Relative_Intact_Eligible_Intensity_Percent"]
        candidate.intact_envelope_qc_score = row["Intact_Envelope_QC_Score"]
        candidate.intact_envelope_qc_rank = row["Intact_Envelope_QC_Rank"]
        candidate.strict_eligible_rank = row["Strict_Eligible_Rank"]
        candidate.review_eligible_rank = row["Review_Eligible_Rank"]
        candidate.dominant_intact_envelope_flag = row["Dominant_Intact_Envelope_Flag"]


    diagnostic_rows = build_intact_reconstruction_diagnostics(qc_rows, reconstruction_config, reconstruction_enabled)
    return qc_rows, diagnostic_rows


def build_intact_reconstruction_diagnostics(
    qc_rows: list[dict[str, Any]],
    reconstruction_config: dict[str, Any],
    reconstruction_enabled: bool = True,
) -> list[dict[str, Any]]:
    qc_config = _qc_config(reconstruction_config or {})
    status_counts = {status: 0 for status in ["Reliable", "Review", "Insufficient", "Failed"]}
    reason_counts: dict[str, int] = {}
    for row in qc_rows:
        status = str(row.get("Reconstruction_Status") or "")
        if status in status_counts:
            status_counts[status] += 1
        for reason in str(row.get("Limiting_Factors") or row.get("Primary_Limiting_Factor") or "").split(";"):
            reason = reason.strip()
            if reason:
                reason_counts[reason] = reason_counts.get(reason, 0) + 1
    if not reconstruction_enabled:
        reason_counts["reconstruction_disabled"] = max(1, reason_counts.get("reconstruction_disabled", 0))
    elif not qc_rows:
        reason_counts["no_charge_state_candidates"] = 1
    reason_summary = "; ".join(f"{key}:{value}" for key, value in sorted(reason_counts.items()))
    in_range_rows = [row for row in qc_rows if row.get("In_Neutral_Mass_Search_Range")]
    outside_rows = [row for row in qc_rows if not row.get("In_Neutral_Mass_Search_Range")]
    dominant = max(qc_rows, key=lambda row: _safe_float(row.get("Total_Supporting_Intensity")), default={})
    dominant_in_range = max(in_range_rows, key=lambda row: _safe_float(row.get("Total_Supporting_Intensity")), default={})
    strict_rows = [row for row in qc_rows if row.get("Intact_Strict_Eligible")]
    review_rows = [row for row in qc_rows if row.get("Intact_Review_Eligible")]
    dominant_strict = _dominant_row(strict_rows)
    dominant_review = _dominant_row(review_rows)
    dominant_eligible = dominant_strict or dominant_review
    target_review_rows = [row for row in qc_rows if row.get("In_Target_Review_Mass_Range")]
    references = qc_config.get("reference_masses") or []
    reference_summary = "; ".join(f"{item['label']}={item['mass_da']}" for item in references) or "not_configured"
    rt_settings = (
        f"reliable<={qc_config['max_rt_range_min_for_reliable']} min; "
        f"review<={qc_config['max_rt_range_min_for_review']} min"
    )
    return [{
        "Total_Reconstruction_Candidates": len(qc_rows),
        "Reliable_Count": status_counts["Reliable"],
        "Review_Count": status_counts["Review"],
        "Insufficient_Count": status_counts["Insufficient"],
        "Failed_Count": status_counts["Failed"],
        "Envelope_QC_Eligible_Count": sum(1 for row in qc_rows if row.get("Envelope_QC_Eligible")),
        "Intact_Strict_Eligible_Count": len(strict_rows),
        "Intact_Review_Eligible_Count": len(review_rows),
        "Comparison_Ready_Strict_Count": sum(1 for row in qc_rows if row.get("Comparison_Ready_Strict")),
        "Comparison_Ready_Review_Count": sum(1 for row in qc_rows if row.get("Comparison_Ready_Review")),
        "Comparison_Ready_Count": sum(1 for row in qc_rows if row.get("Comparison_Ready")),
        "Trace_Only_Envelope_Count": sum(1 for row in qc_rows if row.get("Trace_Only_Envelope")),
        "Noncontiguous_Envelope_Count": sum(1 for row in qc_rows if row.get("Charge_State_Continuity") == "non_contiguous"),
        "RT_Inconsistent_Count": sum(1 for row in qc_rows if row.get("RT_Consistency") == "inconsistent"),
        "Internal_Mass_Error_Count": sum(1 for row in qc_rows if float(row.get("Envelope_Internal_Error_ppm") or 0.0) > qc_config["max_envelope_internal_error_ppm"]),
        "Theory_Near_Match_Count": sum(1 for row in qc_rows if row.get("Unmodified_Theory_Delta_ppm") is not None and abs(float(row.get("Unmodified_Theory_Delta_ppm") or 0.0)) <= qc_config["max_mass_error_ppm"]),
        "Reference_Match_Count": sum(1 for row in qc_rows if row.get("Reference_Mass_Matched")),
        "Dominant_Envelope_Mass": dominant.get("Reconstructed_Mass"),
        "Dominant_Envelope_Intensity": dominant.get("Total_Supporting_Intensity"),
        "Dominant_Envelope_Status": dominant.get("Reconstruction_Status"),
        "Dominant_Envelope_Comparison_Ready": dominant.get("Comparison_Ready"),
        "Dominant_Envelope_Overall_Mass": dominant.get("Reconstructed_Mass"),
        "Dominant_Envelope_Overall_Intensity": dominant.get("Total_Supporting_Intensity"),
        "Dominant_Envelope_In_Mass_Range_Mass": dominant_in_range.get("Reconstructed_Mass"),
        "Dominant_Envelope_In_Mass_Range_Intensity": dominant_in_range.get("Total_Supporting_Intensity"),
        "Dominant_Envelope_In_Mass_Range_Status": dominant_in_range.get("Reconstruction_Status"),
        "Dominant_Envelope_In_Mass_Range_Comparison_Ready": dominant_in_range.get("Comparison_Ready"),
        "Dominant_Envelope_In_Search_Range_Raw_Mass": dominant_in_range.get("Reconstructed_Mass"),
        "Dominant_Envelope_In_Search_Range_Raw_Intensity": dominant_in_range.get("Total_Supporting_Intensity"),
        "Dominant_Intact_Strict_Envelope_Mass": dominant_strict.get("Reconstructed_Mass"),
        "Dominant_Intact_Strict_Envelope_Intensity": dominant_strict.get("Total_Supporting_Intensity"),
        "Dominant_Intact_Strict_QC_Score": dominant_strict.get("Intact_Envelope_QC_Score"),
        "Dominant_Intact_Review_Envelope_Mass": dominant_review.get("Reconstructed_Mass"),
        "Dominant_Intact_Review_Envelope_Intensity": dominant_review.get("Total_Supporting_Intensity"),
        "Dominant_Intact_Review_QC_Score": dominant_review.get("Intact_Envelope_QC_Score"),
        "Dominant_Intact_Eligible_Envelope_Mass": dominant_eligible.get("Reconstructed_Mass"),
        "Dominant_Intact_Eligible_Envelope_Intensity": dominant_eligible.get("Total_Supporting_Intensity"),
        "Dominant_Intact_Eligible_QC_Score": dominant_eligible.get("Intact_Envelope_QC_Score"),
        "Dominant_Intact_Eligible_Reference_Label": dominant_eligible.get("Best_Reference_Label"),
        "Failure_Reason_Counts": reason_summary,
        "Reconstruction_Enabled": reconstruction_enabled,
        "Neutral_Mass_Search_Min_Da": qc_config["neutral_mass_range"]["min_da"],
        "Neutral_Mass_Search_Max_Da": qc_config["neutral_mass_range"]["max_da"],
        "Total_Candidates_Before_Mass_Range_Filter": len(qc_rows),
        "Total_Candidates_In_Mass_Range": len(in_range_rows),
        "Total_Candidates_Outside_Mass_Range": len(outside_rows),
        "Target_Review_Mass_Range_Settings": _target_review_settings(qc_config),
        "Target_Review_Candidate_Count": len(target_review_rows),
        "Search_Mode": qc_config["search_mode"],
        "Intensity_Normalization_Method": "overall, in-range raw, and intact-eligible relative intensity are reported separately",
        "RT_Tolerance_Settings": rt_settings,
        "Reference_Masses_Used": reference_summary,
        "Min_Charge_States_For_Reliable": qc_config["min_charge_states_for_reliable"],
        "Min_Charge_States_For_Review": qc_config["min_charge_states_for_review"],
        "Require_Contiguous_Charge_States": qc_config["require_contiguous_charge_states"],
        "Max_Neutral_Mass_SD_Da": qc_config["max_neutral_mass_sd_da"],
        "Max_Neutral_Mass_Range_Da": qc_config["max_neutral_mass_range_da"],
        "Max_Envelope_Internal_Error_ppm": qc_config["max_envelope_internal_error_ppm"],
        "Min_Relative_Intensity_Percent": qc_config["min_relative_intensity_percent"],
        "Min_Relative_Envelope_Intensity_Percent_For_Reliable": qc_config["min_relative_envelope_intensity_percent_for_reliable"],
        "Min_Relative_Envelope_Intensity_Percent_For_Review": qc_config["min_relative_envelope_intensity_percent_for_review"],
        "Max_Competing_Envelopes": qc_config["max_competing_envelopes"],
        "Comparison_Ready_Statuses": "; ".join(map(str, qc_config["comparison_ready_statuses"])),
        "Notes": "Reliable emphasizes charge-envelope internal quality, RT consistency, and signal support. Comparison_Ready requires intact eligibility, not only neutral mass range membership. Neutral mass search range is the absolute intact reconstruction range; default is 20000-30000 Da. Target review mass range is optional prioritization only. Reference mass matches do not confirm modification identity. Raw in-range dominant is intensity-only; Dominant_Intact_Eligible uses QC eligibility and ranking.",
    }]


def reconstruct_intact_masses(
    tier_result: PeakTierResult,
    reconstruction_config: dict[str, Any],
    instrument_config: dict[str, Any],
    theoretical_mass: float | None = None,
    warnings: list[dict[str, Any]] | None = None,
) -> tuple[list[IntactMassCandidate], list[dict[str, Any]]]:
    min_charge = int(reconstruction_config.get("min_charge", 5))
    max_charge = int(reconstruction_config.get("max_charge", 40))
    min_charge_states = int(reconstruction_config.get("min_charge_states", 3))
    tolerance = float(reconstruction_config.get("mass_cluster_tolerance_da", 1.0))
    polarity = instrument_config.get("polarity", "negative")

    observations = []
    for peak in tier_result.usable_peaks:
        for charge in range(min_charge, max_charge + 1):
            neutral_mass = neutral_mass_from_mz(peak.mz, charge, polarity)
            observations.append({"peak": peak, "charge": charge, "neutral_mass": neutral_mass})

    observations.sort(key=lambda row: row["neutral_mass"])
    clusters: list[list[dict[str, Any]]] = []
    for observation in observations:
        if not clusters:
            clusters.append([observation])
            continue
        current_mean = sum(row["neutral_mass"] for row in clusters[-1]) / len(clusters[-1])
        if abs(observation["neutral_mass"] - current_mean) <= tolerance:
            clusters[-1].append(observation)
        else:
            clusters.append([observation])

    candidates: list[IntactMassCandidate] = []
    charge_state_peaks: list[dict[str, Any]] = []
    for index, cluster in enumerate(clusters, start=1):
        cluster_id = f"C{index:04d}"
        observed_mass = sum(row["neutral_mass"] for row in cluster) / len(cluster)
        charges = sorted({int(row["charge"]) for row in cluster})
        total_intensity = sum(float(row["peak"].intensity) for row in cluster)
        mass_error_da = observed_mass - theoretical_mass if theoretical_mass is not None else None
        mass_error_ppm = (mass_error_da / theoretical_mass * 1_000_000) if theoretical_mass else None
        candidate = IntactMassCandidate(
            observed_mass=observed_mass,
            charge_state_count=len(charges),
            charge_states=charges,
            supporting_peak_count=len(cluster),
            total_intensity=total_intensity,
            theoretical_mass=theoretical_mass,
            mass_error_da=mass_error_da,
            mass_error_ppm=mass_error_ppm,
            confidence=_confidence(len(charges), min_charge_states),
            cluster_id=cluster_id,
        )
        candidates.append(candidate)
        for row in cluster:
            peak = row["peak"]
            charge_state_peaks.append(
                {
                    "Cluster_ID": cluster_id,
                    "mz": peak.mz,
                    "Intensity": peak.intensity,
                    "RT": peak.rt,
                    "Scan_ID": peak.scan_id,
                    "Charge": row["charge"],
                    "Neutral_Mass": row["neutral_mass"],
                    "Peak_Tier": peak.tier,
                }
            )

    candidates.sort(key=lambda item: (item.charge_state_count < min_charge_states, -item.charge_state_count, -item.total_intensity))
    build_intact_reconstruction_qc(candidates, charge_state_peaks, reconstruction_config, reconstruction_enabled=True)
    return candidates, charge_state_peaks
