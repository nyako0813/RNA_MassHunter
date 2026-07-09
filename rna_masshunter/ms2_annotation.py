from typing import Any

import numpy as np

from rna_masshunter.masses import calculate_unmodified_rna_mass
from rna_masshunter.models import Fragment, MS2IonMatch, MS2SpectrumInfo, TheoreticalMS2Ion
from rna_masshunter.ms1_mapping import ppm_error, theoretical_mz_from_mass
from rna_masshunter.mzml_diagnostics import _rt_minutes
from rna_masshunter.mzml_reader import iter_spectra
from rna_masshunter.warnings_manager import add_warning

MS2_SUMMARY_COLUMNS = [
    "Total_MS2_Spectra",
    "Annotated_Spectra",
    "Total_MS2_Peaks",
    "Matched_MS2_Peaks",
    "Unmatched_MS2_Peaks",
    "Total_Theoretical_Ions",
    "Best_Matched_Parent_Fragments",
    "Notes",
]

MS2_SPECTRA_COLUMNS = [
    "Spectrum_ID",
    "Scan_Index",
    "RT",
    "Precursor_mz",
    "Precursor_Charge",
    "Precursor_Intensity",
    "Num_Peaks",
    "Base_Peak_mz",
    "Base_Peak_Intensity",
    "Total_Ion_Current",
    "Num_Matched_Peaks",
    "Num_Unmatched_Peaks",
    "Annotation_Status",
    "Comment",
]

MS2_THEORETICAL_ION_COLUMNS = [
    "Ion_ID",
    "Parent_Fragment_ID",
    "Parent_Sequence",
    "Ion_Type",
    "Ion_Sequence",
    "Ion_Start",
    "Ion_End",
    "Charge",
    "Theoretical_Mass",
    "Theoretical_mz",
    "Modification_ID",
    "Modification_Name",
    "Comment",
]

MS2_ION_MATCH_COLUMNS = [
    "Spectrum_ID",
    "Scan_Index",
    "RT",
    "Precursor_mz",
    "Precursor_Charge",
    "Observed_mz",
    "Observed_Intensity",
    "Best_Ion_ID",
    "Best_Ion_Type",
    "Best_Ion_Sequence",
    "Parent_Fragment_ID",
    "Parent_Fragment_Sequence",
    "Ion_Charge",
    "Theoretical_mz",
    "Mass_Error_Da",
    "Mass_Error_ppm",
    "Match_Status",
    "Confidence",
    "Alternative_Candidates",
    "Comment",
]

MS2_UNMATCHED_COLUMNS = [
    "Spectrum_ID",
    "Scan_Index",
    "RT",
    "Precursor_mz",
    "Precursor_Charge",
    "Observed_mz",
    "Observed_Intensity",
    "Nearest_Ion_ID",
    "Nearest_Ion_Type",
    "Nearest_Ion_Sequence",
    "Nearest_Theoretical_mz",
    "Nearest_Error_Da",
    "Nearest_Error_ppm",
    "Possible_Interpretation",
    "Comment",
]


def annotate_ms2(
    mzml_path: str | None,
    theoretical_fragments: list[Fragment],
    config: Any,
    base_masses: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
) -> dict[str, list[dict[str, Any]]]:
    ms2_config = getattr(config, "ms2_annotation", {}) or {}
    if not _as_bool(ms2_config.get("enabled"), True):
        return {}

    spectra = extract_ms2_spectra(mzml_path, ms2_config, warnings) if mzml_path else []
    ions = generate_theoretical_ms2_ions(theoretical_fragments, config, base_masses, warnings)
    matches, unmatched, spectrum_rows = match_ms2_spectra(spectra, ions, config)

    if not spectra:
        spectrum_rows = [{
            "Spectrum_ID": "",
            "Scan_Index": "",
            "RT": "",
            "Precursor_mz": "",
            "Precursor_Charge": "",
            "Precursor_Intensity": "",
            "Num_Peaks": 0,
            "Base_Peak_mz": "",
            "Base_Peak_Intensity": "",
            "Total_Ion_Current": 0,
            "Num_Matched_Peaks": 0,
            "Num_Unmatched_Peaks": 0,
            "Annotation_Status": "no_ms2_spectra",
            "Comment": "MS2 annotation skipped because no MS2 spectra were available.",
        }]
    elif not ions:
        for row in spectrum_rows:
            row["Annotation_Status"] = "no_theoretical_ions" if row["Num_Peaks"] else "no_peaks"
            row["Comment"] = "MS2 annotation skipped because no theoretical ions were generated."

    return {
        "MS2_Summary": build_ms2_summary(spectra, ions, matches, unmatched),
        "MS2_Spectra": spectrum_rows,
        "MS2_Theoretical_Ions": [theoretical_ion_row(ion) for ion in ions],
        "MS2_Ion_Matches": [ms2_match_row(match) for match in matches] + [_unmatched_match_row(row) for row in unmatched],
        "MS2_Unmatched_Peaks": unmatched,
    }


def extract_ms2_spectra(
    mzml_path: str,
    ms2_config: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
) -> list[MS2SpectrumInfo]:
    spectra: list[MS2SpectrumInfo] = []
    min_intensity = float(ms2_config.get("min_peak_intensity", 0) or 0)
    max_peaks = _optional_positive_int(ms2_config.get("max_peaks_per_spectrum"), 500)

    try:
        iterator = iter_spectra(mzml_path)
        for scan_index, spectrum in enumerate(iterator, start=1):
            try:
                ms_level = int(spectrum.get("ms level", 0) or 0)
            except (TypeError, ValueError):
                ms_level = 0
            if ms_level != 2:
                continue

            mz_array = np.asarray(spectrum.get("m/z array", []), dtype=float)
            intensity_array = np.asarray(spectrum.get("intensity array", []), dtype=float)
            if mz_array.size != intensity_array.size:
                if warnings is not None:
                    add_warning(warnings, "WARNING", "ms2_annotation", "MS2 spectrum m/z and intensity arrays had different lengths.", spectrum.get("id"))
                continue
            if max_peaks and mz_array.size > max_peaks:
                order = np.argsort(intensity_array)[-max_peaks:]
                order = order[np.argsort(mz_array[order])]
                mz_array = mz_array[order]
                intensity_array = intensity_array[order]

            precursor = _precursor_info(spectrum)
            base_peak_index = int(np.argmax(intensity_array)) if intensity_array.size else None
            spectra.append(MS2SpectrumInfo(
                spectrum_id=str(spectrum.get("id") or f"scan_{scan_index}"),
                scan_index=scan_index,
                rt=_safe_float(_rt_minutes(spectrum)),
                precursor_mz=precursor.get("mz"),
                precursor_charge=precursor.get("charge"),
                precursor_intensity=precursor.get("intensity"),
                num_peaks=int(mz_array.size),
                base_peak_mz=float(mz_array[base_peak_index]) if base_peak_index is not None else None,
                base_peak_intensity=float(intensity_array[base_peak_index]) if base_peak_index is not None else None,
                total_ion_current=float(np.sum(intensity_array)) if intensity_array.size else 0.0,
                peaks=[(float(mz), float(intensity)) for mz, intensity in zip(mz_array, intensity_array, strict=False)],
            ))
    except Exception as exc:  # mzML metadata varies enough that extraction should never abort the whole run.
        if warnings is not None:
            add_warning(warnings, "WARNING", "ms2_annotation", "MS2 spectrum extraction failed; annotation skipped.", str(exc))
    return spectra


def generate_theoretical_ms2_ions(
    theoretical_fragments: list[Fragment],
    config: Any,
    base_masses: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
) -> list[TheoreticalMS2Ion]:
    ms2_config = getattr(config, "ms2_annotation", {}) or {}
    if not _as_bool(ms2_config.get("use_theoretical_fragments"), True):
        return []
    min_length = max(1, int(ms2_config.get("min_ion_length", 1) or 1))
    max_length = _optional_positive_int(ms2_config.get("max_ion_length"), None)
    polarity = str(getattr(config, "instrument", {}).get("polarity", "negative") or "negative").lower()
    charges = _ion_charges(ms2_config)

    ions: list[TheoreticalMS2Ion] = []
    for fragment in theoretical_fragments or []:
        sequence = (fragment.sequence or "").upper().replace("T", "U")
        if len(sequence) < 2:
            continue
        max_cut = len(sequence) - 1
        for cut in range(1, max_cut + 1):
            candidates = [
                ("c", sequence[:cut], 1, cut),
                ("y", sequence[cut:], cut + 1, len(sequence)),
            ]
            for ion_type, ion_sequence, ion_start, ion_end in candidates:
                ion_length = len(ion_sequence)
                if ion_length < min_length:
                    continue
                if max_length is not None and ion_length > max_length:
                    continue
                mass = calculate_unmodified_rna_mass(ion_sequence, base_masses, warnings=warnings, terminal_form="default")
                if mass is None:
                    continue
                for charge in charges:
                    ions.append(TheoreticalMS2Ion(
                        ion_id=f"MS2ION_{len(ions) + 1:06d}",
                        parent_fragment_id=fragment.fragment_id,
                        parent_sequence=sequence,
                        ion_type=ion_type,
                        ion_sequence=ion_sequence,
                        ion_start=ion_start,
                        ion_end=ion_end,
                        charge=charge,
                        theoretical_mass=float(mass),
                        theoretical_mz=theoretical_mz_from_mass(float(mass), charge, polarity),
                        modification_id="",
                        modification_name="",
                        comment="unmodified c/y series ion",
                    ))
    return ions


def match_ms2_spectra(
    spectra: list[MS2SpectrumInfo],
    ions: list[TheoreticalMS2Ion],
    config: Any,
) -> tuple[list[MS2IonMatch], list[dict[str, Any]], list[dict[str, Any]]]:
    ms2_config = getattr(config, "ms2_annotation", {}) or {}
    tolerance_ppm = float(ms2_config.get("mz_tolerance_ppm", 20) or 20)
    min_intensity = float(ms2_config.get("min_peak_intensity", 0) or 0)

    matches: list[MS2IonMatch] = []
    unmatched_rows: list[dict[str, Any]] = []
    spectrum_rows: list[dict[str, Any]] = []

    for spectrum in spectra:
        spectrum_match_count = 0
        spectrum_unmatched_count = 0
        for observed_mz, observed_intensity in spectrum.peaks:
            if observed_intensity < min_intensity:
                spectrum_unmatched_count += 1
                unmatched_rows.append(_unmatched_row(spectrum, observed_mz, observed_intensity, None, None, None, "noise/background", "Peak was below min_peak_intensity."))
                continue
            best, alternatives, nearest = _best_ion_match(observed_mz, ions, tolerance_ppm)
            if best is None:
                spectrum_unmatched_count += 1
                nearest_ion, nearest_error_da, nearest_error_ppm = nearest
                unmatched_rows.append(_unmatched_row(
                    spectrum,
                    observed_mz,
                    observed_intensity,
                    nearest_ion,
                    nearest_error_da,
                    nearest_error_ppm,
                    _possible_interpretation(nearest_error_ppm, observed_intensity, min_intensity),
                    "No theoretical c/y ion was within tolerance.",
                ))
                continue

            ion, error_da, error_ppm = best
            spectrum_match_count += 1
            status = "multiple_candidates" if alternatives else "matched"
            matches.append(MS2IonMatch(
                spectrum_id=spectrum.spectrum_id,
                scan_index=spectrum.scan_index,
                rt=spectrum.rt,
                precursor_mz=spectrum.precursor_mz,
                precursor_charge=spectrum.precursor_charge,
                observed_mz=observed_mz,
                observed_intensity=observed_intensity,
                best_ion_id=ion.ion_id,
                best_ion_type=ion.ion_type,
                best_ion_sequence=ion.ion_sequence,
                parent_fragment_id=ion.parent_fragment_id,
                parent_fragment_sequence=ion.parent_sequence,
                ion_charge=ion.charge,
                theoretical_mz=ion.theoretical_mz,
                mass_error_da=error_da,
                mass_error_ppm=error_ppm,
                match_status=status,
                confidence=_confidence(error_ppm, tolerance_ppm, observed_intensity, spectrum.base_peak_intensity, status),
                alternative_candidates="; ".join(candidate[0].ion_id for candidate in alternatives[:5]),
                comment="",
            ))

        status = "annotated" if spectrum_match_count else ("no_theoretical_ions" if not ions else "no_peaks" if not spectrum.peaks else "skipped")
        spectrum_rows.append(spectrum_row(spectrum, spectrum_match_count, spectrum_unmatched_count, status))

    return matches, unmatched_rows, spectrum_rows


def build_ms2_summary(
    spectra: list[MS2SpectrumInfo],
    ions: list[TheoreticalMS2Ion],
    matches: list[MS2IonMatch],
    unmatched_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    parent_counts: dict[str, int] = {}
    for match in matches:
        parent_counts[match.parent_fragment_id] = parent_counts.get(match.parent_fragment_id, 0) + 1
    best_parents = sorted(parent_counts.items(), key=lambda item: (-item[1], item[0]))[:10]
    if not spectra:
        notes = "MS2 spectra: 0; MS2 annotation skipped."
    elif not ions:
        notes = "MS2 spectra were found, but no theoretical ions were generated."
    else:
        notes = "MS2 c/y ion annotation against theoretical digestion fragments."
    return [{
        "Total_MS2_Spectra": len(spectra),
        "Annotated_Spectra": len({match.spectrum_id for match in matches}),
        "Total_MS2_Peaks": sum(spectrum.num_peaks for spectrum in spectra),
        "Matched_MS2_Peaks": len(matches),
        "Unmatched_MS2_Peaks": len(unmatched_rows),
        "Total_Theoretical_Ions": len(ions),
        "Best_Matched_Parent_Fragments": "; ".join(f"{fragment_id} ({count})" for fragment_id, count in best_parents),
        "Notes": notes,
    }]


def spectrum_row(spectrum: MS2SpectrumInfo, matched: int, unmatched: int, status: str) -> dict[str, Any]:
    return {
        "Spectrum_ID": spectrum.spectrum_id,
        "Scan_Index": spectrum.scan_index,
        "RT": spectrum.rt,
        "Precursor_mz": spectrum.precursor_mz,
        "Precursor_Charge": spectrum.precursor_charge,
        "Precursor_Intensity": spectrum.precursor_intensity,
        "Num_Peaks": spectrum.num_peaks,
        "Base_Peak_mz": spectrum.base_peak_mz,
        "Base_Peak_Intensity": spectrum.base_peak_intensity,
        "Total_Ion_Current": spectrum.total_ion_current,
        "Num_Matched_Peaks": matched,
        "Num_Unmatched_Peaks": unmatched,
        "Annotation_Status": status,
        "Comment": "" if status == "annotated" else "MS2 annotation skipped or produced no matches.",
    }


def theoretical_ion_row(ion: TheoreticalMS2Ion) -> dict[str, Any]:
    return {
        "Ion_ID": ion.ion_id,
        "Parent_Fragment_ID": ion.parent_fragment_id,
        "Parent_Sequence": ion.parent_sequence,
        "Ion_Type": ion.ion_type,
        "Ion_Sequence": ion.ion_sequence,
        "Ion_Start": ion.ion_start,
        "Ion_End": ion.ion_end,
        "Charge": ion.charge,
        "Theoretical_Mass": ion.theoretical_mass,
        "Theoretical_mz": ion.theoretical_mz,
        "Modification_ID": ion.modification_id,
        "Modification_Name": ion.modification_name,
        "Comment": ion.comment,
    }


def ms2_match_row(match: MS2IonMatch) -> dict[str, Any]:
    return {
        "Spectrum_ID": match.spectrum_id,
        "Scan_Index": match.scan_index,
        "RT": match.rt,
        "Precursor_mz": match.precursor_mz,
        "Precursor_Charge": match.precursor_charge,
        "Observed_mz": match.observed_mz,
        "Observed_Intensity": match.observed_intensity,
        "Best_Ion_ID": match.best_ion_id,
        "Best_Ion_Type": match.best_ion_type,
        "Best_Ion_Sequence": match.best_ion_sequence,
        "Parent_Fragment_ID": match.parent_fragment_id,
        "Parent_Fragment_Sequence": match.parent_fragment_sequence,
        "Ion_Charge": match.ion_charge,
        "Theoretical_mz": match.theoretical_mz,
        "Mass_Error_Da": match.mass_error_da,
        "Mass_Error_ppm": match.mass_error_ppm,
        "Match_Status": match.match_status,
        "Confidence": match.confidence,
        "Alternative_Candidates": match.alternative_candidates,
        "Comment": match.comment,
    }


def _unmatched_match_row(row: dict[str, Any]) -> dict[str, Any]:
    interpretation = str(row.get("Possible_Interpretation") or "")
    status = "low_intensity" if interpretation == "noise/background" else "outside_tolerance"
    return {
        "Spectrum_ID": row.get("Spectrum_ID"),
        "Scan_Index": row.get("Scan_Index"),
        "RT": row.get("RT"),
        "Precursor_mz": row.get("Precursor_mz"),
        "Precursor_Charge": row.get("Precursor_Charge"),
        "Observed_mz": row.get("Observed_mz"),
        "Observed_Intensity": row.get("Observed_Intensity"),
        "Best_Ion_ID": "",
        "Best_Ion_Type": "",
        "Best_Ion_Sequence": "",
        "Parent_Fragment_ID": "",
        "Parent_Fragment_Sequence": "",
        "Ion_Charge": "",
        "Theoretical_mz": "",
        "Mass_Error_Da": "",
        "Mass_Error_ppm": "",
        "Match_Status": status,
        "Confidence": "Unmatched",
        "Alternative_Candidates": "",
        "Comment": row.get("Comment"),
    }


def _best_ion_match(observed_mz: float, ions: list[TheoreticalMS2Ion], tolerance_ppm: float):
    if not ions:
        return None, [], (None, None, None)
    ranked = sorted(
        ((ion, observed_mz - ion.theoretical_mz, ppm_error(observed_mz, ion.theoretical_mz)) for ion in ions),
        key=lambda item: abs(item[2]),
    )
    nearest = ranked[0]
    within = [item for item in ranked if abs(item[2]) <= tolerance_ppm]
    if not within:
        return None, [], nearest
    return within[0], within[1:], nearest


def _unmatched_row(
    spectrum: MS2SpectrumInfo,
    observed_mz: float,
    observed_intensity: float,
    nearest_ion: TheoreticalMS2Ion | None,
    nearest_error_da: float | None,
    nearest_error_ppm: float | None,
    interpretation: str,
    comment: str,
) -> dict[str, Any]:
    return {
        "Spectrum_ID": spectrum.spectrum_id,
        "Scan_Index": spectrum.scan_index,
        "RT": spectrum.rt,
        "Precursor_mz": spectrum.precursor_mz,
        "Precursor_Charge": spectrum.precursor_charge,
        "Observed_mz": observed_mz,
        "Observed_Intensity": observed_intensity,
        "Nearest_Ion_ID": nearest_ion.ion_id if nearest_ion else "",
        "Nearest_Ion_Type": nearest_ion.ion_type if nearest_ion else "",
        "Nearest_Ion_Sequence": nearest_ion.ion_sequence if nearest_ion else "",
        "Nearest_Theoretical_mz": nearest_ion.theoretical_mz if nearest_ion else "",
        "Nearest_Error_Da": nearest_error_da if nearest_error_da is not None else "",
        "Nearest_Error_ppm": nearest_error_ppm if nearest_error_ppm is not None else "",
        "Possible_Interpretation": interpretation,
        "Comment": comment,
    }


def _precursor_info(spectrum: dict[str, Any]) -> dict[str, Any]:
    info = {"mz": None, "charge": None, "intensity": None}
    precursor_list = spectrum.get("precursorList", {}).get("precursor", [])
    if not precursor_list:
        return info
    selected = precursor_list[0].get("selectedIonList", {}).get("selectedIon", [])
    if not selected:
        return info
    ion = selected[0]
    info["mz"] = _safe_float(ion.get("selected ion m/z") or ion.get("isolation window target m/z"))
    charge = ion.get("charge state")
    try:
        info["charge"] = int(charge) if charge is not None else None
    except (TypeError, ValueError):
        info["charge"] = None
    info["intensity"] = _safe_float(ion.get("peak intensity") or ion.get("ion intensity"))
    return info


def _ion_charges(ms2_config: dict[str, Any]) -> list[int]:
    raw = ms2_config.get("charge_states") or [1]
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


def _confidence(error_ppm: float, tolerance_ppm: float, intensity: float, base_intensity: float | None, status: str) -> str:
    if status == "multiple_candidates":
        return "Low"
    relative = intensity / base_intensity if base_intensity else 1.0
    if abs(error_ppm) <= tolerance_ppm / 3 and relative >= 0.1:
        return "High"
    if abs(error_ppm) <= tolerance_ppm:
        return "Medium"
    return "Low"


def _possible_interpretation(nearest_error_ppm: float | None, intensity: float, min_intensity: float) -> str:
    if intensity < min_intensity:
        return "noise/background"
    if nearest_error_ppm is None:
        return "unknown"
    if abs(nearest_error_ppm) <= 100:
        return "modified fragment ion candidate"
    return "outside tolerance"


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


def _safe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
