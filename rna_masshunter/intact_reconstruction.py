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
    "min_relative_intensity_percent": 0.5,
    "max_competing_envelopes": 3,
    "comparison_ready_statuses": ["Reliable", "Review"],
}

QC_COLUMNS = [
    "Cluster_ID",
    "Observed_Mass",
    "Reconstruction_Status",
    "Reconstruction_Confidence",
    "Num_Supporting_Charge_States",
    "Charge_State_Range",
    "Charge_State_Continuity",
    "Neutral_Mass_SD",
    "Neutral_Mass_Range",
    "Max_Mass_Error_ppm",
    "Total_Supporting_Intensity",
    "Competing_Envelope_Count",
    "Primary_Limiting_Factor",
    "Comparison_Ready",
]

DIAGNOSTIC_COLUMNS = [
    "Total_Reconstruction_Candidates",
    "Reliable_Count",
    "Review_Count",
    "Insufficient_Count",
    "Failed_Count",
    "Comparison_Ready_Count",
    "Failure_Reason_Counts",
    "Reconstruction_Enabled",
    "Min_Charge_States_For_Reliable",
    "Min_Charge_States_For_Review",
    "Require_Contiguous_Charge_States",
    "Max_Neutral_Mass_SD_Da",
    "Max_Neutral_Mass_Range_Da",
    "Max_Mass_Error_ppm",
    "Min_Relative_Intensity_Percent",
    "Max_Competing_Envelopes",
    "Comparison_Ready_Statuses",
    "Notes",
]


def _as_bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _qc_config(reconstruction_config: dict[str, Any]) -> dict[str, Any]:
    raw = reconstruction_config.get("intact_reconstruction") or reconstruction_config.get("qc") or {}
    merged = {**DEFAULT_INTACT_QC_CONFIG, **raw}
    merged["min_charge_states_for_reliable"] = int(merged.get("min_charge_states_for_reliable") or 3)
    merged["min_charge_states_for_review"] = int(merged.get("min_charge_states_for_review") or 2)
    merged["require_contiguous_charge_states"] = _as_bool(merged.get("require_contiguous_charge_states"), True)
    merged["max_neutral_mass_sd_da"] = float(merged.get("max_neutral_mass_sd_da") or 0.5)
    merged["max_neutral_mass_range_da"] = float(merged.get("max_neutral_mass_range_da") or 1.5)
    merged["max_mass_error_ppm"] = float(merged.get("max_mass_error_ppm") or 20)
    merged["min_relative_intensity_percent"] = float(merged.get("min_relative_intensity_percent") or 0.5)
    merged["max_competing_envelopes"] = int(merged.get("max_competing_envelopes") or 3)
    statuses = merged.get("comparison_ready_statuses") or ["Reliable", "Review"]
    if isinstance(statuses, str):
        statuses = [item.strip() for item in statuses.split(",") if item.strip()]
    merged["comparison_ready_statuses"] = list(statuses)
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
        "non_contiguous_charge_states",
        "mass_spread_too_large",
        "mass_error_too_large",
        "insufficient_intensity_support",
        "multiple_competing_envelopes",
    ]
    for item in priority:
        if item in factors:
            return item
    return factors[0] if factors else ""


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
        neutral_sd = pstdev(neutral_masses) if len(neutral_masses) > 1 else 0.0 if neutral_masses else None
        neutral_range = (max(neutral_masses) - min(neutral_masses)) if neutral_masses else None
        continuity = _charge_continuity(charges)
        competing = 0
        for other in candidates:
            if other is candidate:
                continue
            if abs(float(other.observed_mass) - float(candidate.observed_mass)) <= qc_config["max_neutral_mass_range_da"]:
                competing += 1

        factors: list[str] = []
        if len(charges) < qc_config["min_charge_states_for_review"]:
            factors.append("insufficient_charge_states")
        if qc_config["require_contiguous_charge_states"] and continuity == "non_contiguous":
            factors.append("non_contiguous_charge_states")
        if (neutral_sd is not None and neutral_sd > qc_config["max_neutral_mass_sd_da"]) or (
            neutral_range is not None and neutral_range > qc_config["max_neutral_mass_range_da"]
        ):
            factors.append("mass_spread_too_large")
        max_mass_error_ppm = abs(float(candidate.mass_error_ppm)) if candidate.mass_error_ppm is not None else None
        if max_mass_error_ppm is not None and max_mass_error_ppm > qc_config["max_mass_error_ppm"]:
            factors.append("mass_error_too_large")
        relative_intensity = (float(candidate.total_intensity or 0.0) / max_intensity * 100.0) if max_intensity else 0.0
        if max_intensity and relative_intensity < qc_config["min_relative_intensity_percent"]:
            factors.append("insufficient_intensity_support")
        if competing > qc_config["max_competing_envelopes"]:
            factors.append("multiple_competing_envelopes")

        if not reconstruction_enabled:
            status = "Failed"
            factors.insert(0, "reconstruction_disabled")
        elif len(charges) < qc_config["min_charge_states_for_review"]:
            status = "Insufficient"
        elif len(charges) >= qc_config["min_charge_states_for_reliable"] and not factors:
            status = "Reliable"
        else:
            status = "Review"
        confidence = {"Reliable": "High", "Review": "Medium", "Insufficient": "Low", "Failed": "None"}.get(status, "Low")
        comparison_ready = status in qc_config["comparison_ready_statuses"]
        primary_factor = _primary_factor(factors)

        candidate.reconstruction_status = status
        candidate.reconstruction_confidence = confidence
        candidate.num_supporting_charge_states = len(charges)
        candidate.charge_state_range = _charge_state_range(charges)
        candidate.charge_state_continuity = continuity
        candidate.neutral_mass_sd = neutral_sd
        candidate.neutral_mass_range = neutral_range
        candidate.max_mass_error_ppm = max_mass_error_ppm
        candidate.total_supporting_intensity = float(candidate.total_intensity or 0.0)
        candidate.competing_envelope_count = competing
        candidate.primary_limiting_factor = primary_factor
        candidate.comparison_ready = comparison_ready

        qc_rows.append({
            "Cluster_ID": cluster_id,
            "Observed_Mass": candidate.observed_mass,
            "Reconstruction_Status": status,
            "Reconstruction_Confidence": confidence,
            "Num_Supporting_Charge_States": len(charges),
            "Charge_State_Range": candidate.charge_state_range,
            "Charge_State_Continuity": continuity,
            "Neutral_Mass_SD": neutral_sd,
            "Neutral_Mass_Range": neutral_range,
            "Max_Mass_Error_ppm": max_mass_error_ppm,
            "Total_Supporting_Intensity": candidate.total_supporting_intensity,
            "Competing_Envelope_Count": competing,
            "Primary_Limiting_Factor": primary_factor,
            "Comparison_Ready": comparison_ready,
        })

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
        reason = str(row.get("Primary_Limiting_Factor") or "")
        if reason:
            reason_counts[reason] = reason_counts.get(reason, 0) + 1
    if not reconstruction_enabled:
        reason_counts["reconstruction_disabled"] = max(1, reason_counts.get("reconstruction_disabled", 0))
    elif not qc_rows:
        reason_counts["no_charge_state_candidates"] = 1
    reason_summary = "; ".join(f"{key}:{value}" for key, value in sorted(reason_counts.items()))
    return [{
        "Total_Reconstruction_Candidates": len(qc_rows),
        "Reliable_Count": status_counts["Reliable"],
        "Review_Count": status_counts["Review"],
        "Insufficient_Count": status_counts["Insufficient"],
        "Failed_Count": status_counts["Failed"],
        "Comparison_Ready_Count": sum(1 for row in qc_rows if row.get("Comparison_Ready")),
        "Failure_Reason_Counts": reason_summary,
        "Reconstruction_Enabled": reconstruction_enabled,
        "Min_Charge_States_For_Reliable": qc_config["min_charge_states_for_reliable"],
        "Min_Charge_States_For_Review": qc_config["min_charge_states_for_review"],
        "Require_Contiguous_Charge_States": qc_config["require_contiguous_charge_states"],
        "Max_Neutral_Mass_SD_Da": qc_config["max_neutral_mass_sd_da"],
        "Max_Neutral_Mass_Range_Da": qc_config["max_neutral_mass_range_da"],
        "Max_Mass_Error_ppm": qc_config["max_mass_error_ppm"],
        "Min_Relative_Intensity_Percent": qc_config["min_relative_intensity_percent"],
        "Max_Competing_Envelopes": qc_config["max_competing_envelopes"],
        "Comparison_Ready_Statuses": "; ".join(map(str, qc_config["comparison_ready_statuses"])),
        "Notes": "Comparison_Ready indicates intact-mass quality for cross-condition comparison; it does not confirm modification identity.",
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
