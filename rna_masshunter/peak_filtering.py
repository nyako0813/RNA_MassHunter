from typing import Any

from rna_masshunter.models import Peak, PeakTierResult


def classify_peak_tiers(peaks: list[Peak], peak_filtering_config: dict[str, Any], warnings: list[dict[str, Any]] | None = None) -> PeakTierResult:
    major = float(peak_filtering_config.get("major_intensity_threshold", 25000))
    minor = float(peak_filtering_config.get("minor_intensity_threshold", 5000))
    trace = float(peak_filtering_config.get("trace_intensity_threshold", 1000))
    result = PeakTierResult()
    for peak in peaks:
        if peak.intensity >= major:
            peak.tier = "Major"
            result.major.append(peak)
        elif peak.intensity >= minor:
            peak.tier = "Minor"
            result.minor.append(peak)
        elif peak.intensity >= trace:
            peak.tier = "Trace"
            result.trace.append(peak)
        else:
            peak.tier = "Below reporting threshold"
            result.below_threshold.append(peak)
    return result
