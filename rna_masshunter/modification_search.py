from dataclasses import asdict, is_dataclass
from math import isfinite
from typing import Any

from rna_masshunter.masses import neutral_mass_from_mz
from rna_masshunter.models import FragmentMS1Match, IntactMassCandidate, KnownModificationCandidate, Modification, RunConfig
from rna_masshunter.warnings_manager import add_warning


def _as_bool(value: Any, default: bool = True) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _as_positive_int(value: Any, default: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed > 0 else default


def _as_float(value: Any, default: float) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return parsed if isfinite(parsed) else default


def _normalize_values(values: Any) -> set[str]:
    if values is None:
        return set()
    if isinstance(values, str):
        return {values.lower()}
    return {str(value).lower() for value in values}


def _ppm_error(observed_mass: float, theoretical_mass: float) -> float:
    if theoretical_mass == 0:
        return 0.0
    return (observed_mass - theoretical_mass) / theoretical_mass * 1_000_000


def _confidence_score(confidence: str | None) -> int:
    return {"high": 3, "medium": 2, "low": 1}.get(str(confidence or "").lower(), 0)


def _peak_tier_score(peak_tier: str | None) -> int:
    return {"major": 3, "minor": 2, "trace": 1}.get(str(peak_tier or "").lower(), 0)


def _modification_name(modification: Modification) -> str:
    for key in ("name", "modification_name", "full_name", "description"):
        value = modification.raw.get(key) if isinstance(modification.raw, dict) else None
        if value:
            return str(value)
    return modification.id or modification.symbol or "Unknown modification"


def _target_bases(modification: Modification) -> list[str]:
    bases = modification.target_bases or []
    normalized = []
    for base in bases:
        value = str(base).upper().replace("T", "U")
        if value and value not in {"UNKNOWN", "N", "ANY"}:
            normalized.append(value)
    return normalized


def _target_base_label(modification: Modification) -> str:
    bases = _target_bases(modification)
    return ",".join(bases) if bases else "any"


def _base_compatible(sequence: str, modification: Modification) -> bool:
    bases = _target_bases(modification)
    if not bases:
        return True
    sequence = (sequence or "").upper().replace("T", "U")
    return any(base in sequence for base in bases)


def _position_overlap(standard_start: int | None, standard_end: int | None, prioritized_positions: list[int]) -> tuple[float, bool]:
    if standard_start is None or standard_end is None:
        return 0.0, False
    start = min(standard_start, standard_end)
    end = max(standard_start, standard_end)
    overlap = any(start <= position <= end for position in prioritized_positions)
    wobble_overlap = start <= 34 <= end
    return (1.0 if overlap else 0.0), wobble_overlap


def _priority_score(
    confidence: str,
    peak_tier: str | None,
    position_overlap_score: float,
    mass_error_modified_ppm: float,
    tolerance_ppm: float,
) -> float:
    tolerance = tolerance_ppm if tolerance_ppm > 0 else 1.0
    return (
        _confidence_score(confidence)
        + _peak_tier_score(peak_tier)
        + position_overlap_score
        - abs(mass_error_modified_ppm) / tolerance
    )


def _candidate_rows(candidates: list[KnownModificationCandidate]) -> list[dict[str, Any]]:
    rows = []
    for item in candidates:
        row = asdict(item) if is_dataclass(item) else dict(item)
        warnings = row.get("warnings", [])
        if isinstance(warnings, list):
            row["warnings"] = "; ".join(map(str, warnings))
        rows.append(row)
    return rows


def known_modification_candidate_rows(candidates: list[KnownModificationCandidate]) -> list[dict[str, Any]]:
    return _candidate_rows(candidates)


def summarize_known_modification_candidates(candidates: list[KnownModificationCandidate]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str, str, bool], list[KnownModificationCandidate]] = {}
    for candidate in candidates:
        target_base = candidate.target_base
        key = (
            candidate.modification_id,
            candidate.modification_name,
            candidate.modification_symbol or "",
            target_base,
            candidate.wobble_overlap,
        )
        grouped.setdefault(key, []).append(candidate)

    rows = []
    for (mod_id, mod_name, symbol, target_base, wobble_overlap), items in grouped.items():
        best = sorted(items, key=lambda item: (-item.priority_score, abs(item.mass_error_modified_ppm), -item.intensity))[0]
        rows.append(
            {
                "Modification_ID": mod_id,
                "Modification_Name": mod_name,
                "Symbol": symbol,
                "Target_Base": target_base,
                "Wobble_Overlap": wobble_overlap,
                "Candidate_Count": len(items),
                "Wobble_Overlap_Count": sum(1 for item in items if item.wobble_overlap),
                "Best_Source_ID": best.source_id,
                "Best_Sequence": best.sequence,
                "Best_Standard_Start": best.standard_start,
                "Best_Standard_End": best.standard_end,
                "Best_Mass_Error_Modified_ppm": best.mass_error_modified_ppm,
                "Best_Intensity": best.intensity,
                "Best_Peak_Tier": best.peak_tier,
                "Best_Confidence": best.confidence,
                "Best_Priority_Score": best.priority_score,
            }
        )
    return sorted(rows, key=lambda row: (-row["Best_Priority_Score"], row["Modification_ID"], str(row["Wobble_Overlap"])))


def _fragment_candidates(
    fragment_matches: list[FragmentMS1Match],
    modifications: list[Modification],
    config: RunConfig,
    warnings: list[dict[str, Any]] | None,
) -> list[KnownModificationCandidate]:
    search_config = config.modification_search or {}
    tolerance_ppm = _as_float(search_config.get("mz_tolerance_ppm"), 10.0)
    max_candidates = _as_positive_int(search_config.get("max_candidates_per_match"), 10)
    prioritized = [int(value) for value in search_config.get("prioritize_standard_positions", [34]) or []]
    require_base = _as_bool(search_config.get("require_base_compatibility"), True)
    allowed_tiers = _normalize_values(search_config.get("min_peak_tier", ["Major", "Minor"]))
    allowed_confidence = _normalize_values(search_config.get("min_confidence", ["High", "Medium"]))
    polarity = str((config.fragment_mapping or {}).get("polarity", "auto") or "auto").lower()
    if polarity == "auto":
        polarity = str((config.instrument or {}).get("polarity", "negative") or "negative").lower()

    candidates: list[KnownModificationCandidate] = []
    for match in fragment_matches:
        if allowed_tiers and str(match.peak_tier or "").lower() not in allowed_tiers:
            continue
        if allowed_confidence and str(match.confidence or "").lower() not in allowed_confidence:
            continue
        observed_mass = neutral_mass_from_mz(match.observed_mz, match.charge, polarity)
        unmodified_mass = match.fragment_mass
        mass_error_unmodified_da = observed_mass - unmodified_mass
        mass_error_unmodified_ppm = _ppm_error(observed_mass, unmodified_mass)
        position_score, wobble_overlap = _position_overlap(match.standard_start, match.standard_end, prioritized)
        per_match: list[KnownModificationCandidate] = []
        for modification in modifications:
            shift = modification.mass_shift_from_unmodified
            if not isfinite(shift):
                continue
            if require_base and not _base_compatible(match.sequence, modification):
                continue
            modified_mass = unmodified_mass + shift
            mass_error_modified_da = observed_mass - modified_mass
            mass_error_modified_ppm = _ppm_error(observed_mass, modified_mass)
            if abs(mass_error_modified_ppm) > tolerance_ppm:
                continue
            priority = _priority_score(match.confidence, match.peak_tier, position_score, mass_error_modified_ppm, tolerance_ppm)
            per_match.append(
                KnownModificationCandidate(
                    candidate_id="",
                    source_type="fragment",
                    source_id=match.fragment_id,
                    target_id=match.target_id,
                    sequence=match.sequence,
                    start=match.start,
                    end=match.end,
                    standard_start=match.standard_start,
                    standard_end=match.standard_end,
                    observed_mz=match.observed_mz,
                    theoretical_mz=match.theoretical_mz,
                    observed_mass=observed_mass,
                    unmodified_mass=unmodified_mass,
                    mass_error_unmodified_da=mass_error_unmodified_da,
                    mass_error_unmodified_ppm=mass_error_unmodified_ppm,
                    modification_id=modification.id,
                    modification_symbol=modification.symbol,
                    modification_name=_modification_name(modification),
                    target_base=_target_base_label(modification),
                    modification_mass_shift=shift,
                    modified_mass=modified_mass,
                    mass_error_modified_da=mass_error_modified_da,
                    mass_error_modified_ppm=mass_error_modified_ppm,
                    charge=match.charge,
                    intensity=match.intensity,
                    rt=match.rt,
                    peak_tier=match.peak_tier,
                    confidence=match.confidence,
                    position_overlap_score=position_score,
                    wobble_overlap=wobble_overlap,
                    priority_score=priority,
                    notes="base compatible" if _base_compatible(match.sequence, modification) else "base compatibility not required",
                )
            )
        per_match.sort(key=lambda item: (-item.priority_score, abs(item.mass_error_modified_ppm), -item.intensity))
        candidates.extend(per_match[:max_candidates])
        if len(per_match) > max_candidates and warnings is not None:
            add_warning(
                warnings,
                "WARNING",
                "modification_search",
                "Known modification candidates were truncated by max_candidates_per_match.",
                {"source_id": match.fragment_id, "before": len(per_match), "after": max_candidates},
            )
    return candidates


def _intact_candidates(
    intact_results: list[IntactMassCandidate],
    modifications: list[Modification],
    config: RunConfig,
) -> list[KnownModificationCandidate]:
    search_config = config.modification_search or {}
    tolerance_ppm = _as_float(search_config.get("mz_tolerance_ppm"), 10.0)
    max_candidates = _as_positive_int(search_config.get("max_candidates_per_match"), 10)
    candidates: list[KnownModificationCandidate] = []
    sequence = (config.sequence or {}).get("sequence", "") or ""
    target_id = (config.sequence or {}).get("name", "target_RNA") or "target_RNA"
    for result in intact_results:
        theoretical = result.theoretical_mass
        if theoretical is None:
            continue
        observed = result.observed_mass
        per_result: list[KnownModificationCandidate] = []
        for modification in modifications:
            shift = modification.mass_shift_from_unmodified
            if not isfinite(shift):
                continue
            modified_mass = theoretical + shift
            error_da = observed - modified_mass
            error_ppm = _ppm_error(observed, modified_mass)
            if abs(error_ppm) > tolerance_ppm:
                continue
            priority = _priority_score(result.confidence, None, 0.0, error_ppm, tolerance_ppm)
            per_result.append(
                KnownModificationCandidate(
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
                    modification_id=modification.id,
                    modification_symbol=modification.symbol,
                    modification_name=_modification_name(modification),
                    target_base=_target_base_label(modification),
                    modification_mass_shift=shift,
                    modified_mass=modified_mass,
                    mass_error_modified_da=error_da,
                    mass_error_modified_ppm=error_ppm,
                    charge=None,
                    intensity=result.total_intensity,
                    rt=None,
                    peak_tier=None,
                    confidence=result.confidence,
                    position_overlap_score=0.0,
                    wobble_overlap=False,
                    priority_score=priority,
                    notes="intact mode preliminary candidate",
                )
            )
        per_result.sort(key=lambda item: (-item.priority_score, abs(item.mass_error_modified_ppm), -item.intensity))
        candidates.extend(per_result[:max_candidates])
    return candidates


def search_known_modifications(
    fragment_ms1_matches: list[FragmentMS1Match],
    intact_results: list[IntactMassCandidate],
    modifications: list[Modification],
    config: RunConfig,
    warnings: list[dict[str, Any]] | None = None,
) -> list[KnownModificationCandidate]:
    search_config = config.modification_search or {}
    if not _as_bool(search_config.get("enabled"), True):
        return []

    source = search_config.get("source", {}) if isinstance(search_config.get("source", {}), dict) else {}
    use_fragments = _as_bool(source.get("use_fragments"), True)
    use_intact = _as_bool(source.get("use_intact"), False)
    candidates: list[KnownModificationCandidate] = []
    if use_fragments:
        if fragment_ms1_matches:
            candidates.extend(_fragment_candidates(fragment_ms1_matches, modifications, config, warnings))
        elif warnings is not None:
            add_warning(warnings, "INFO", "modification_search", "Fragment modification search was enabled but Fragment_MS1_matches is empty.")
    if use_intact:
        candidates.extend(_intact_candidates(intact_results, modifications, config))

    candidates.sort(key=lambda item: (-item.priority_score, abs(item.mass_error_modified_ppm), -item.intensity))
    for index, candidate in enumerate(candidates, start=1):
        candidate.candidate_id = f"KMOD_{index:06d}"
    return candidates
