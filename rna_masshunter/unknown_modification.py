"""Mass-shift-only search for modifications not present in the curated
modification catalog. Tries simple elemental deltas (e.g. +O, +S) directly
against raw MS1 peaks via binary search, independent of any known
modification ID or prior unmodified-only match.
"""
from __future__ import annotations

from math import isfinite
from typing import Any

from rna_masshunter.mass_shift_ms1_search import SortedPeakIndex, build_sorted_peak_index, find_peaks_near_mz
from rna_masshunter.masses import elemental_delta_mass, mz_from_neutral_mass, neutral_mass_from_mz
from rna_masshunter.models import (
    CompoundModificationCandidate, IntactMassCandidate, Modification, RunConfig, UnknownModificationCandidate,
)
from rna_masshunter.modification_search import (
    _as_bool, _as_float, _as_positive_int, _base_compatible, _candidate_policy_allows_mass_search,
    _modification_name, _normalize_values, _ppm_error, _priority_score, _should_skip_isobaric_shift,
    _target_base_label,
)
from rna_masshunter.ms1_mapping import _eligible_peaks
from rna_masshunter.warnings_manager import add_warning



def _elements_label(elements: dict[str, int]) -> str:
    parts = []
    for element, count in sorted(elements.items()):
        sign = "+" if count >= 0 else "-"
        parts.append(f"{sign}{element}{abs(count) if abs(count) != 1 else ''}")
    return "".join(parts)


def _fragment_unknown_candidates(
    theoretical_fragments: list[Any],
    peak_index: SortedPeakIndex,
    deltas: list[dict[str, Any]],
    config: RunConfig,
    warnings: list[dict[str, Any]] | None,
) -> list[UnknownModificationCandidate]:
    search_config = config.unknown_modification_search or {}
    tolerance_ppm = _as_float(search_config.get("mz_tolerance_ppm"), 10.0)
    max_candidates = _as_positive_int(search_config.get("max_candidates_per_match"), 10)
    allowed_tiers = _normalize_values(search_config.get("min_peak_tier", ["Major", "Minor"]))
    allowed_confidence = _normalize_values(search_config.get("min_confidence", ["High", "Medium"]))
    mapping_config = config.fragment_mapping or {}
    polarity = str(mapping_config.get("polarity", "auto") or "auto").lower()
    if polarity == "auto":
        polarity = str((config.instrument or {}).get("polarity", "negative") or "negative").lower()
    min_charge = _as_positive_int(mapping_config.get("min_charge"), 1)
    max_charge = _as_positive_int(mapping_config.get("max_charge"), 8)

    candidates: list[UnknownModificationCandidate] = []
    for fragment in theoretical_fragments:
        per_fragment: list[UnknownModificationCandidate] = []
        for delta in deltas:
            elements = delta.get("elements", {})
            try:
                shift = elemental_delta_mass(elements)
            except ValueError:
                continue
            if not isfinite(shift):
                continue
            modified_mass = fragment.unmodified_mass + shift
            if modified_mass <= 0:
                continue
            for charge in range(min_charge, max_charge + 1):
                theoretical_mz = mz_from_neutral_mass(modified_mass, charge, polarity)
                for peak_match in find_peaks_near_mz(peak_index, theoretical_mz, tolerance_ppm):
                    tier = getattr(peak_match.peak, "tier", None)
                    if allowed_tiers and str(tier or "").lower() not in allowed_tiers:
                        continue
                    if allowed_confidence and peak_match.confidence.lower() not in allowed_confidence:
                        continue
                    observed_mass = neutral_mass_from_mz(peak_match.observed_mz, charge, polarity)
                    unmodified_mass = fragment.unmodified_mass
                    mass_error_modified_ppm = _ppm_error(observed_mass, modified_mass)
                    priority = _priority_score(peak_match.confidence, tier, mass_error_modified_ppm, tolerance_ppm)
                    per_fragment.append(
                        UnknownModificationCandidate(
                            candidate_id="",
                            source_type="fragment",
                            source_id=fragment.fragment_id,
                            target_id=fragment.target_id,
                            sequence=fragment.sequence,
                            start=fragment.start,
                            end=fragment.end,
                            standard_start=fragment.standard_start,
                            standard_end=fragment.standard_end,
                            observed_mz=peak_match.observed_mz,
                            theoretical_mz=theoretical_mz,
                            observed_mass=observed_mass,
                            unmodified_mass=unmodified_mass,
                            mass_error_unmodified_da=observed_mass - unmodified_mass,
                            mass_error_unmodified_ppm=_ppm_error(observed_mass, unmodified_mass),
                            delta_label=delta.get("label", _elements_label(elements)),
                            delta_elements=_elements_label(elements),
                            delta_mass_shift=shift,
                            modified_mass=modified_mass,
                            mass_error_modified_da=observed_mass - modified_mass,
                            mass_error_modified_ppm=mass_error_modified_ppm,
                            charge=charge,
                            intensity=float(getattr(peak_match.peak, "intensity", 0.0) or 0.0),
                            rt=getattr(peak_match.peak, "rt", None),
                            peak_tier=tier,
                            confidence=peak_match.confidence,
                            priority_score=priority,
                            notes="mass-shift-only candidate; not in curated modification catalog",
                        )
                    )
        per_fragment.sort(key=lambda item: (-item.priority_score, abs(item.mass_error_modified_ppm), -item.intensity))
        candidates.extend(per_fragment[:max_candidates])
        if len(per_fragment) > max_candidates and warnings is not None:
            add_warning(
                warnings, "WARNING", "unknown_modification_search",
                "Unknown modification candidates were truncated by max_candidates_per_match.",
                {"source_id": fragment.fragment_id, "before": len(per_fragment), "after": max_candidates},
            )
    return candidates


def _intact_unknown_candidates(
    intact_results: list[IntactMassCandidate],
    deltas: list[dict[str, Any]],
    config: RunConfig,
) -> list[UnknownModificationCandidate]:
    search_config = config.unknown_modification_search or {}
    tolerance_ppm = _as_float(search_config.get("mz_tolerance_ppm"), 10.0)
    max_candidates = _as_positive_int(search_config.get("max_candidates_per_match"), 10)
    sequence = (config.sequence or {}).get("sequence", "") or ""
    target_id = (config.sequence or {}).get("name", "target_RNA") or "target_RNA"

    candidates: list[UnknownModificationCandidate] = []
    for result in intact_results:
        theoretical = result.theoretical_mass
        if theoretical is None:
            continue
        observed = result.observed_mass
        per_result: list[UnknownModificationCandidate] = []
        for delta in deltas:
            elements = delta.get("elements", {})
            try:
                shift = elemental_delta_mass(elements)
            except ValueError:
                continue
            if not isfinite(shift):
                continue
            modified_mass = theoretical + shift
            error_da = observed - modified_mass
            error_ppm = _ppm_error(observed, modified_mass)
            if abs(error_ppm) > tolerance_ppm:
                continue
            priority = _priority_score(result.confidence, None, error_ppm, tolerance_ppm)
            per_result.append(
                UnknownModificationCandidate(
                    candidate_id="",
                    source_type="intact",
                    source_id=result.cluster_id or "intact_mass",
                    target_id=target_id,
                    sequence=sequence,
                    start=1 if sequence else None,
                    end=len(sequence) if sequence else None,
                    standard_start=None,
                    standard_end=None,
                    observed_mz=None,
                    theoretical_mz=None,
                    observed_mass=observed,
                    unmodified_mass=theoretical,
                    mass_error_unmodified_da=observed - theoretical,
                    mass_error_unmodified_ppm=_ppm_error(observed, theoretical),
                    delta_label=delta.get("label", _elements_label(elements)),
                    delta_elements=_elements_label(elements),
                    delta_mass_shift=shift,
                    modified_mass=modified_mass,
                    mass_error_modified_da=error_da,
                    mass_error_modified_ppm=error_ppm,
                    charge=None,
                    intensity=result.total_intensity,
                    rt=None,
                    peak_tier=None,
                    confidence=result.confidence,
                    priority_score=priority,
                    notes="mass-shift-only candidate; not in curated modification catalog",
                )
            )
        per_result.sort(key=lambda item: (-item.priority_score, abs(item.mass_error_modified_ppm), -item.intensity))
        candidates.extend(per_result[:max_candidates])
    return candidates


def generate_unknown_modification_candidates(
    theoretical_fragments: list[Any],
    peaks: list[Any],
    intact_results: list[IntactMassCandidate],
    config: RunConfig,
    warnings: list[dict[str, Any]] | None = None,
) -> list[UnknownModificationCandidate]:
    search_config = config.unknown_modification_search or {}
    if not _as_bool(search_config.get("enabled"), True):
        return []
    deltas = search_config.get("candidate_deltas", [])
    if not deltas:
        return []
    source = search_config.get("source", {}) if isinstance(search_config.get("source", {}), dict) else {}
    use_fragments = _as_bool(source.get("use_fragments"), True)
    use_intact = _as_bool(source.get("use_intact"), False)

    candidates: list[UnknownModificationCandidate] = []
    if use_fragments:
        if theoretical_fragments and peaks:
            eligible = _eligible_peaks(peaks, config.fragment_mapping or {})
            peak_index = build_sorted_peak_index(eligible)
            candidates.extend(_fragment_unknown_candidates(theoretical_fragments, peak_index, deltas, config, warnings))
        elif warnings is not None:
            add_warning(warnings, "INFO", "unknown_modification_search", "Unknown modification search was enabled but theoretical fragments or MS1 peaks are empty.")
    if use_intact:
        candidates.extend(_intact_unknown_candidates(intact_results, deltas, config))

    candidates.sort(key=lambda item: (-item.priority_score, abs(item.mass_error_modified_ppm), -item.intensity))
    for index, candidate in enumerate(candidates, start=1):
        candidate.candidate_id = f"UMOD_{index:06d}"
    return candidates


def summarize_unknown_modification_candidates(candidates: list[UnknownModificationCandidate]) -> list[dict[str, Any]]:
    grouped: dict[str, list[UnknownModificationCandidate]] = {}
    for candidate in candidates:
        grouped.setdefault(candidate.delta_label, []).append(candidate)
    rows = []
    for delta_label, items in grouped.items():
        best = sorted(items, key=lambda item: (-item.priority_score, abs(item.mass_error_modified_ppm), -item.intensity))[0]
        rows.append(
            {
                "Delta_Label": delta_label,
                "Delta_Elements": best.delta_elements,
                "Delta_Mass_Shift": best.delta_mass_shift,
                "Candidate_Count": len(items),
                "Best_Source_ID": best.source_id,
                "Best_Sequence": best.sequence,
                "Best_Mass_Error_Modified_ppm": best.mass_error_modified_ppm,
                "Best_Intensity": best.intensity,
                "Best_Peak_Tier": best.peak_tier,
                "Best_Confidence": best.confidence,
                "Best_Priority_Score": best.priority_score,
            }
        )
    return sorted(rows, key=lambda row: (-row["Best_Priority_Score"], row["Delta_Label"]))

def _fragment_compound_candidates(
    theoretical_fragments: list[Any],
    peak_index: SortedPeakIndex,
    modifications: list[Modification],
    deltas: list[dict[str, Any]],
    config: RunConfig,
    warnings: list[dict[str, Any]] | None,
) -> list[CompoundModificationCandidate]:
    """Search unmodified_mass + known_modification_shift + extra_elemental_delta,
    e.g. a known modification (ncm5s2U) plus an additional +S/+O on top of it."""
    search_config = config.unknown_modification_search or {}
    tolerance_ppm = _as_float(search_config.get("mz_tolerance_ppm"), 10.0)
    max_candidates = _as_positive_int(search_config.get("max_candidates_per_match"), 10)
    allowed_tiers = _normalize_values(search_config.get("min_peak_tier", ["Major", "Minor"]))
    allowed_confidence = _normalize_values(search_config.get("min_confidence", ["High", "Medium"]))
    require_base = _as_bool((config.modification_search or {}).get("require_base_compatibility"), True)
    mapping_config = config.fragment_mapping or {}
    polarity = str(mapping_config.get("polarity", "auto") or "auto").lower()
    if polarity == "auto":
        polarity = str((config.instrument or {}).get("polarity", "negative") or "negative").lower()
    min_charge = _as_positive_int(mapping_config.get("min_charge"), 1)
    max_charge = _as_positive_int(mapping_config.get("max_charge"), 8)

    candidates: list[CompoundModificationCandidate] = []
    for fragment in theoretical_fragments:
        per_fragment: list[CompoundModificationCandidate] = []
        for modification in modifications:
            if not _candidate_policy_allows_mass_search(modification):
                continue
            known_shift = modification.mass_shift_from_unmodified
            if not isfinite(known_shift):
                continue
            if require_base and not _base_compatible(fragment.sequence, modification):
                continue
            for delta in deltas:
                elements = delta.get("elements", {})
                try:
                    delta_shift = elemental_delta_mass(elements)
                except ValueError:
                    continue
                if not isfinite(delta_shift):
                    continue
                combined_shift = known_shift + delta_shift
                if _should_skip_isobaric_shift(combined_shift, search_config):
                    continue
                modified_mass = fragment.unmodified_mass + combined_shift
                if modified_mass <= 0:
                    continue
                for charge in range(min_charge, max_charge + 1):
                    theoretical_mz = mz_from_neutral_mass(modified_mass, charge, polarity)
                    for peak_match in find_peaks_near_mz(peak_index, theoretical_mz, tolerance_ppm):
                        tier = getattr(peak_match.peak, "tier", None)
                        if allowed_tiers and str(tier or "").lower() not in allowed_tiers:
                            continue
                        if allowed_confidence and peak_match.confidence.lower() not in allowed_confidence:
                            continue
                        observed_mass = neutral_mass_from_mz(peak_match.observed_mz, charge, polarity)
                        unmodified_mass = fragment.unmodified_mass
                        mass_error_modified_ppm = _ppm_error(observed_mass, modified_mass)
                        priority = _priority_score(peak_match.confidence, tier, mass_error_modified_ppm, tolerance_ppm)
                        per_fragment.append(
                            CompoundModificationCandidate(
                                candidate_id="",
                                source_type="fragment",
                                source_id=fragment.fragment_id,
                                target_id=fragment.target_id,
                                sequence=fragment.sequence,
                                start=fragment.start,
                                end=fragment.end,
                                standard_start=fragment.standard_start,
                                standard_end=fragment.standard_end,
                                observed_mz=peak_match.observed_mz,
                                theoretical_mz=theoretical_mz,
                                observed_mass=observed_mass,
                                unmodified_mass=unmodified_mass,
                                mass_error_unmodified_da=observed_mass - unmodified_mass,
                                mass_error_unmodified_ppm=_ppm_error(observed_mass, unmodified_mass),
                                modification_id=modification.id,
                                modification_symbol=modification.symbol,
                                modification_name=_modification_name(modification),
                                target_base=_target_base_label(modification),
                                modification_mass_shift=known_shift,
                                delta_label=delta.get("label", _elements_label(elements)),
                                delta_elements=_elements_label(elements),
                                delta_mass_shift=delta_shift,
                                combined_mass_shift=combined_shift,
                                modified_mass=modified_mass,
                                mass_error_modified_da=observed_mass - modified_mass,
                                mass_error_modified_ppm=mass_error_modified_ppm,
                                charge=charge,
                                intensity=float(getattr(peak_match.peak, "intensity", 0.0) or 0.0),
                                rt=getattr(peak_match.peak, "rt", None),
                                peak_tier=tier,
                                confidence=peak_match.confidence,
                                priority_score=priority,
                                notes=f"known modification ({modification.symbol or modification.id}) plus additional mass shift; not in curated catalog as a combined entry",
                            )
                        )
        per_fragment.sort(key=lambda item: (-item.priority_score, abs(item.mass_error_modified_ppm), -item.intensity))
        candidates.extend(per_fragment[:max_candidates])
        if len(per_fragment) > max_candidates and warnings is not None:
            add_warning(
                warnings, "WARNING", "compound_modification_search",
                "Compound modification candidates were truncated by max_candidates_per_match.",
                {"source_id": fragment.fragment_id, "before": len(per_fragment), "after": max_candidates},
            )
    return candidates


def generate_compound_modification_candidates(
    theoretical_fragments: list[Any],
    peaks: list[Any],
    modifications: list[Modification],
    config: RunConfig,
    warnings: list[dict[str, Any]] | None = None,
) -> list[CompoundModificationCandidate]:
    search_config = config.unknown_modification_search or {}
    if not _as_bool(search_config.get("enabled"), True):
        return []
    if not _as_bool(search_config.get("include_known_modification_composites"), True):
        return []
    deltas = search_config.get("candidate_deltas", [])
    if not deltas or not modifications:
        return []
    if not (theoretical_fragments and peaks):
        if warnings is not None:
            add_warning(warnings, "INFO", "compound_modification_search", "Compound modification search was enabled but theoretical fragments or MS1 peaks are empty.")
        return []

    eligible = _eligible_peaks(peaks, config.fragment_mapping or {})
    peak_index = build_sorted_peak_index(eligible)
    candidates = _fragment_compound_candidates(theoretical_fragments, peak_index, modifications, deltas, config, warnings)

    candidates.sort(key=lambda item: (-item.priority_score, abs(item.mass_error_modified_ppm), -item.intensity))
    for index, candidate in enumerate(candidates, start=1):
        candidate.candidate_id = f"CMOD_{index:06d}"
    return candidates


def summarize_compound_modification_candidates(candidates: list[CompoundModificationCandidate]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[CompoundModificationCandidate]] = {}
    for candidate in candidates:
        key = (candidate.modification_id, candidate.delta_label)
        grouped.setdefault(key, []).append(candidate)
    rows = []
    for (mod_id, delta_label), items in grouped.items():
        best = sorted(items, key=lambda item: (-item.priority_score, abs(item.mass_error_modified_ppm), -item.intensity))[0]
        rows.append(
            {
                "Modification_ID": mod_id,
                "Delta_Label": delta_label,
                "Combined_Mass_Shift": best.combined_mass_shift,
                "Candidate_Count": len(items),
                "Best_Source_ID": best.source_id,
                "Best_Sequence": best.sequence,
                "Best_Mass_Error_Modified_ppm": best.mass_error_modified_ppm,
                "Best_Intensity": best.intensity,
                "Best_Peak_Tier": best.peak_tier,
                "Best_Confidence": best.confidence,
                "Best_Priority_Score": best.priority_score,
            }
        )
    return sorted(rows, key=lambda row: (-row["Best_Priority_Score"], row["Modification_ID"]))