from typing import Any

from rna_masshunter.masses import neutral_mass_from_mz
from rna_masshunter.models import IntactMassCandidate, PeakTierResult


def _confidence(charge_state_count: int, min_charge_states: int) -> str:
    if charge_state_count >= min_charge_states + 2:
        return "High"
    if charge_state_count >= min_charge_states:
        return "Medium"
    return "Low"


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
    return candidates, charge_state_peaks
