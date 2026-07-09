from dataclasses import asdict, is_dataclass
from typing import Any

from rna_masshunter.enzymes import normalize_enzyme_name
from rna_masshunter.masses import (
    DEFAULT_PHOSPHATE_MASS,
    PROTON_MASS,
    calculate_unmodified_rna_mass,
    neutral_mass_from_mz,
)
from rna_masshunter.models import Modification, PeakTierResult

BASES = ("A", "C", "G", "U")

P1_THEORETICAL_COLUMNS = [
    "Structure_ID",
    "Base",
    "Sequence",
    "Length",
    "Composition",
    "Terminal_Form",
    "Phosphate_State",
    "Modification_ID",
    "Modification_Name",
    "Theoretical_Mass",
    "Theoretical_mz_by_charge",
    "Charge",
    "Comment",
]

P1_ANNOTATION_COLUMNS = [
    "Observed_mz",
    "Observed_neutral_mass",
    "Charge",
    "Intensity",
    "RT",
    "Best_Matched_Structure_ID",
    "Best_Matched_Structure",
    "Best_Matched_Base",
    "Best_Modification_ID",
    "Best_Modification_Name",
    "Best_Theoretical_Mass",
    "Mass_Error_Da",
    "Mass_Error_ppm",
    "Match_Status",
    "Confidence",
    "Alternative_Candidates",
    "Unmatched_Reason",
    "Comment",
]

P1_UNMATCHED_COLUMNS = [
    "Observed_mz",
    "Observed_neutral_mass",
    "Charge",
    "Intensity",
    "RT",
    "Nearest_Theoretical_Structure",
    "Nearest_Theoretical_Mass",
    "Nearest_Mass_Difference_Da",
    "Nearest_Mass_Difference_ppm",
    "Unmatched_Reason",
    "Possible_Interpretation",
    "Comment",
]

P1_SUMMARY_COLUMNS = [
    "Total_Observed_Peaks",
    "Matched_Peaks",
    "Unmatched_Peaks",
    "Multiple_Candidate_Peaks",
    "Major_Matched_Structures",
    "Major_Unmatched_Masses",
    "Notes",
]


def is_p1_enabled(config: Any) -> bool:
    digestion = getattr(config, "digestion", {}) or {}
    p1_config = getattr(config, "p1_annotation", {}) or {}
    return normalize_enzyme_name(digestion.get("enzyme", "")) == "Nuclease_P1" and _as_bool(p1_config.get("enabled"), True)


def build_p1_optional_results(
    config: Any,
    tier_result: PeakTierResult,
    base_masses: dict[str, Any],
    modifications: list[Modification] | None = None,
) -> dict[str, list[dict[str, Any]]]:
    structures = build_p1_theoretical_structures(config, base_masses, modifications)
    peak_results = annotate_p1_peaks(config, tier_result, structures)
    return {
        "P1_Summary": build_p1_summary(peak_results["annotations"], peak_results["unmatched"]),
        "P1_Theoretical_Structures": structures,
        "P1_Peak_Annotations": peak_results["annotations"],
        "P1_Unmatched_Peaks": peak_results["unmatched"],
    }


def build_p1_theoretical_structures(
    config: Any,
    base_masses: dict[str, Any],
    modifications: list[Modification] | None = None,
) -> list[dict[str, Any]]:
    p1_config = getattr(config, "p1_annotation", {}) or {}
    digestion = getattr(config, "digestion", {}) or {}
    include_modified = _as_bool(p1_config.get("include_modified_monomers"), True)
    include_phosphate = _as_bool(p1_config.get("include_phosphate_forms"), True)
    max_length = _optional_positive_int(digestion.get("max_length"), 1) or 1
    charges = _charges(config)
    polarity = getattr(config, "instrument", {}).get("polarity", "negative")

    terminal_forms = ["default"]
    if include_phosphate:
        terminal_forms.append("residual_phosphate")
        ap_config = getattr(config, "alkaline_phosphatase", {}) or {}
        if _as_bool(ap_config.get("enabled"), False) and not _as_bool(ap_config.get("assume_complete"), False):
            if _as_bool(ap_config.get("allow_cyclic_phosphate"), True):
                terminal_forms.append("cyclic_phosphate")

    rows: list[dict[str, Any]] = []
    for sequence in _candidate_sequences(config, max_length):
        for terminal_form in terminal_forms:
            _add_structure_row(rows, sequence, terminal_form, None, charges, polarity, base_masses)
        if include_modified and len(sequence) == 1:
            for modification in modifications or []:
                if sequence in modification.target_bases:
                    _add_structure_row(rows, sequence, "default", modification, charges, polarity, base_masses)
                    if include_phosphate:
                        _add_structure_row(rows, sequence, "residual_phosphate", modification, charges, polarity, base_masses)
    return rows


def annotate_p1_peaks(
    config: Any,
    tier_result: PeakTierResult,
    structures: list[dict[str, Any]],
) -> dict[str, list[dict[str, Any]]]:
    p1_config = getattr(config, "p1_annotation", {}) or {}
    instrument = getattr(config, "instrument", {}) or {}
    tolerance_ppm = float(p1_config.get("mz_tolerance_ppm") or instrument.get("ms1_tolerance_ppm", 10) or 10)
    min_intensity = float(p1_config.get("min_intensity", 0) or 0)
    include_unmatched = _as_bool(p1_config.get("include_unmatched_peaks"), True)
    include_below_threshold = _as_bool(p1_config.get("include_below_threshold_peaks"), False)
    polarity = instrument.get("polarity", "negative")
    charges = _charges(config)

    annotation_rows: list[dict[str, Any]] = []
    unmatched_rows: list[dict[str, Any]] = []
    for peak in _observed_peaks(tier_result, include_below_threshold):
        raw = asdict(peak) if is_dataclass(peak) else dict(peak)
        mz = float(raw.get("mz") or 0)
        intensity = float(raw.get("intensity") or 0)
        best = _best_peak_match(mz, charges, polarity, structures, tolerance_ppm)

        if best["matched"] is None:
            charge = best["charge"] or charges[0]
            observed_mass = neutral_mass_from_mz(mz, charge, polarity)
            nearest = best["nearest"]
            nearest_mass = float(nearest["Theoretical_Mass"]) if nearest else None
            reason = "low_intensity" if intensity < min_intensity else "outside_tolerance"
            status = "low_intensity" if intensity < min_intensity else "possible_adduct_or_unknown"
            annotation_rows.append({
                "Observed_mz": mz,
                "Observed_neutral_mass": observed_mass,
                "Charge": charge,
                "Intensity": intensity,
                "RT": raw.get("rt"),
                "Best_Matched_Structure_ID": "",
                "Best_Matched_Structure": "",
                "Best_Matched_Base": "",
                "Best_Modification_ID": "",
                "Best_Modification_Name": "",
                "Best_Theoretical_Mass": "",
                "Mass_Error_Da": "",
                "Mass_Error_ppm": "",
                "Match_Status": status,
                "Confidence": "Unmatched",
                "Alternative_Candidates": "",
                "Unmatched_Reason": reason,
                "Comment": "No P1 theoretical structure was within tolerance.",
            })
            if include_unmatched:
                unmatched_rows.append({
                    "Observed_mz": mz,
                    "Observed_neutral_mass": observed_mass,
                    "Charge": charge,
                    "Intensity": intensity,
                    "RT": raw.get("rt"),
                    "Nearest_Theoretical_Structure": _structure_label(nearest) if nearest else "",
                    "Nearest_Theoretical_Mass": nearest_mass or "",
                    "Nearest_Mass_Difference_Da": best["nearest_error"] if best["nearest_error"] is not None else "",
                    "Nearest_Mass_Difference_ppm": best["nearest_ppm"] if best["nearest_ppm"] is not None else "",
                    "Unmatched_Reason": reason,
                    "Possible_Interpretation": _possible_interpretation(best["nearest_error"], observed_mass, nearest_mass, intensity, min_intensity),
                    "Comment": "Retained for unknown modification/adduct/phosphate-state review.",
                })
            continue

        row, error_da, error_ppm, observed_mass, charge = best["matched"]
        alternatives = best["alternatives"]
        status = "multiple_candidates" if alternatives else "matched"
        confidence = "High" if abs(error_ppm) <= tolerance_ppm / 2 else "Medium"
        annotation_rows.append({
            "Observed_mz": mz,
            "Observed_neutral_mass": observed_mass,
            "Charge": charge,
            "Intensity": intensity,
            "RT": raw.get("rt"),
            "Best_Matched_Structure_ID": row.get("Structure_ID"),
            "Best_Matched_Structure": _structure_label(row),
            "Best_Matched_Base": row.get("Base"),
            "Best_Modification_ID": row.get("Modification_ID"),
            "Best_Modification_Name": row.get("Modification_Name"),
            "Best_Theoretical_Mass": row.get("Theoretical_Mass"),
            "Mass_Error_Da": error_da,
            "Mass_Error_ppm": error_ppm,
            "Match_Status": status,
            "Confidence": confidence,
            "Alternative_Candidates": "; ".join(_structure_label(item[0]) for item in alternatives[:5]),
            "Unmatched_Reason": "",
            "Comment": "",
        })

    return {"annotations": annotation_rows, "unmatched": unmatched_rows}


def build_p1_summary(annotation_rows: list[dict[str, Any]], unmatched_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    matched = [row for row in annotation_rows if row.get("Match_Status") == "matched"]
    multiple = [row for row in annotation_rows if row.get("Match_Status") == "multiple_candidates"]
    major_matched = sorted(matched + multiple, key=lambda row: float(row.get("Intensity") or 0), reverse=True)[:10]
    major_unmatched = sorted(unmatched_rows, key=lambda row: float(row.get("Intensity") or 0), reverse=True)[:10]
    return [{
        "Total_Observed_Peaks": len(annotation_rows),
        "Matched_Peaks": len(matched),
        "Unmatched_Peaks": len(unmatched_rows),
        "Multiple_Candidate_Peaks": len(multiple),
        "Major_Matched_Structures": "; ".join(str(row.get("Best_Matched_Structure")) for row in major_matched if row.get("Best_Matched_Structure")),
        "Major_Unmatched_Masses": "; ".join(str(round(float(row.get("Observed_neutral_mass") or 0), 4)) for row in major_unmatched),
        "Notes": "Unmatched peaks are retained for unknown modification, adduct, phosphate-state, incomplete digestion, and background review.",
    }]


def _add_structure_row(
    rows: list[dict[str, Any]],
    sequence: str,
    terminal_form: str,
    modification: Modification | None,
    charges: list[int],
    polarity: str,
    base_masses: dict[str, Any],
) -> None:
    mass = calculate_unmodified_rna_mass(sequence, base_masses, terminal_form=terminal_form)
    if mass is None:
        return
    if modification is not None:
        mass += modification.mass_shift_from_unmodified
    if mass != mass:
        return
    for charge in charges:
        mod_id = modification.id if modification else ""
        rows.append({
            "Structure_ID": f"P1_{sequence}_{terminal_form}_{mod_id or 'unmodified'}_z{charge}",
            "Base": sequence if len(sequence) == 1 else "mixed",
            "Sequence": sequence,
            "Length": len(sequence),
            "Composition": _composition(sequence),
            "Terminal_Form": terminal_form,
            "Phosphate_State": _phosphate_state(terminal_form),
            "Modification_ID": mod_id,
            "Modification_Name": _modification_name(modification),
            "Theoretical_Mass": float(mass),
            "Theoretical_mz_by_charge": _mz_from_neutral_mass(float(mass), charge, polarity),
            "Charge": charge,
            "Comment": "modified monomer candidate" if modification else _phosphate_state(terminal_form),
        })


def _best_peak_match(mz: float, charges: list[int], polarity: str, structures: list[dict[str, Any]], tolerance_ppm: float) -> dict[str, Any]:
    best_within = None
    best_alternatives = []
    nearest = None
    nearest_error = None
    nearest_ppm = None
    nearest_charge = None
    for charge in charges:
        observed_mass = neutral_mass_from_mz(mz, charge, polarity)
        same_charge = [row for row in structures if int(row.get("Charge") or 1) == charge] or structures
        ranked = sorted(same_charge, key=lambda row: abs(observed_mass - float(row.get("Theoretical_Mass") or 0)))
        if not ranked:
            continue
        current = ranked[0]
        error_da = observed_mass - float(current["Theoretical_Mass"])
        error_ppm = error_da / float(current["Theoretical_Mass"]) * 1_000_000
        if nearest_error is None or abs(error_da) < abs(nearest_error):
            nearest = current
            nearest_error = error_da
            nearest_ppm = error_ppm
            nearest_charge = charge

        within = []
        for row in ranked:
            theoretical_mass = float(row["Theoretical_Mass"])
            row_error_da = observed_mass - theoretical_mass
            row_error_ppm = row_error_da / theoretical_mass * 1_000_000
            if abs(row_error_ppm) <= tolerance_ppm:
                within.append((row, row_error_da, row_error_ppm, observed_mass, charge))
        if within and (best_within is None or abs(within[0][1]) < abs(best_within[1])):
            best_within = within[0]
            best_alternatives = within[1:]
    return {
        "matched": best_within,
        "alternatives": best_alternatives,
        "nearest": nearest,
        "nearest_error": nearest_error,
        "nearest_ppm": nearest_ppm,
        "charge": nearest_charge,
    }


def _candidate_sequences(config: Any, max_length: int) -> list[str]:
    sequences = list(BASES)
    source_sequence = (getattr(config, "sequence", {}) or {}).get("sequence", "").upper().replace("T", "U")
    for length in range(2, max_length + 1):
        for start in range(0, max(0, len(source_sequence) - length + 1)):
            candidate = source_sequence[start:start + length]
            if set(candidate) <= set(BASES) and candidate not in sequences:
                sequences.append(candidate)
    return sequences


def _observed_peaks(tier_result: PeakTierResult, include_below_threshold: bool) -> list[Any]:
    peaks = list(tier_result.usable_peaks)
    if include_below_threshold:
        peaks.extend(tier_result.below_threshold)
    return peaks


def _charges(config: Any) -> list[int]:
    p1_config = getattr(config, "p1_annotation", {}) or {}
    raw = p1_config.get("charge_states") or p1_config.get("charges") or [1]
    if isinstance(raw, int):
        raw = [raw]
    values = raw if isinstance(raw, list) else [raw]
    charges = []
    for value in values:
        try:
            charge = abs(int(value))
        except (TypeError, ValueError):
            continue
        if charge > 0 and charge not in charges:
            charges.append(charge)
    return charges or [1]


def _as_bool(value: Any, default: bool = True) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _optional_positive_int(value: Any, default: int | None = None) -> int | None:
    if value in (None, ""):
        return default
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed > 0 else default


def _mz_from_neutral_mass(neutral_mass: float, charge: int, polarity: str) -> float:
    z = abs(int(charge))
    if str(polarity).lower() == "negative":
        return (float(neutral_mass) - z * PROTON_MASS) / z
    return (float(neutral_mass) + z * PROTON_MASS) / z


def _phosphate_state(terminal_form: str) -> str:
    return {
        "default": "dephosphorylated form",
        "dephosphorylated": "dephosphorylated form",
        "residual_phosphate": "residual phosphate form",
        "cyclic_phosphate": "cyclic phosphate form",
    }.get(terminal_form, terminal_form)


def _composition(sequence: str) -> str:
    return "; ".join(f"{base}:{sequence.count(base)}" for base in BASES if sequence.count(base))


def _modification_name(modification: Modification | None) -> str:
    if modification is None:
        return ""
    return str(modification.raw.get("name") or modification.symbol or modification.id)


def _structure_label(row: dict[str, Any] | None) -> str:
    if not row:
        return ""
    prefix = "p" if row.get("Terminal_Form") == "residual_phosphate" else ""
    mod_name = row.get("Modification_Name")
    if mod_name:
        return f"{prefix}modified {row.get('Base')} ({mod_name})"
    return f"{prefix}{row.get('Sequence')}"


def _possible_interpretation(
    diff_da: float | None,
    observed_mass: float,
    nearest_mass: float | None,
    intensity: float,
    min_intensity: float,
) -> str:
    if intensity < min_intensity:
        return "noise/background"
    if diff_da is None:
        return "unknown"
    abs_diff = abs(diff_da)
    if abs(abs_diff - DEFAULT_PHOSPHATE_MASS) <= 1.0:
        return "phosphate state mismatch"
    if abs(abs_diff - 21.9819) <= 0.5 or abs(abs_diff - 37.9559) <= 0.5:
        return "sodium/potassium adduct candidate"
    if nearest_mass is not None and observed_mass > nearest_mass + 100:
        return "incomplete digestion product"
    return "outside mass tolerance" if abs_diff else "unknown"
