from typing import Any

import numpy as np

from rna_masshunter.masses import calculate_unmodified_rna_mass
from rna_masshunter.modified_precursor import find_modified_parent_candidates
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
    "Spectra_With_Precursor_Parent_Match",
    "Spectra_Without_Precursor_Parent_Match",
    "Spectra_With_Unmodified_Precursor_Match",
    "Spectra_With_Modified_Precursor_Match",
    "Spectra_Rescued_By_Modified_Precursor",
    "Total_Modified_Parent_Candidates",
    "Top_Modification_Candidates_For_Precursors",
    "Total_Parent_Candidates",
    "Total_Matched_Ion_Rows",
    "Total_Unmatched_Output_Rows",
    "Strong_Evidence_Fragments",
    "Moderate_Evidence_Fragments",
    "Weak_Evidence_Fragments",
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
    "Num_Parent_Candidates",
    "Num_Unmodified_Parent_Candidates",
    "Num_Modified_Parent_Candidates",
    "Precursor_Annotation_Status",
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
    "Ion_Length",
    "Informative_Ion",
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
    "Ion_Length",
    "Informative_Ion",
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

MS2_PARENT_CANDIDATE_COLUMNS = [
    "Spectrum_ID",
    "Scan_Index",
    "RT",
    "Precursor_mz",
    "Precursor_Charge",
    "Candidate_Parent_Fragment_ID",
    "Candidate_Parent_Sequence",
    "Candidate_Parent_Start",
    "Candidate_Parent_End",
    "Candidate_Type",
    "Modification_ID",
    "Modification_Name",
    "Modification_Target_Base",
    "Modification_Mass_Shift",
    "Parent_Unmodified_Mass",
    "Parent_Modified_Mass",
    "Parent_Charge",
    "Parent_Theoretical_mz",
    "Precursor_Error_Da",
    "Precursor_Error_ppm",
    "Parent_Match_Status",
    "Modified_Precursor_Rescue",
    "Parent_Candidate_Rank",
    "Comment",
]

MS2_FRAGMENT_EVIDENCE_COLUMNS = [
    "Spectrum_ID",
    "RT",
    "Precursor_mz",
    "Parent_Fragment_ID",
    "Parent_Sequence",
    "Parent_Start",
    "Parent_End",
    "Parent_Charge",
    "Candidate_Type",
    "Modification_ID",
    "Modification_Name",
    "Modification_Mass_Shift",
    "Parent_Modified_Mass",
    "Modified_Precursor_Rescue",
    "Match_Approximation",
    "Num_Matched_Ions",
    "Num_Informative_Ions",
    "Num_c_Ions",
    "Num_y_Ions",
    "Best_Mass_Error_ppm",
    "Median_Abs_Error_ppm",
    "Total_Matched_Intensity",
    "Sequence_Coverage",
    "Evidence_Score",
    "Evidence_Level",
    "Notes",
]

MS2_MODIFIED_PRECURSOR_COLUMNS = [
    "Spectrum_ID", "Scan_Index", "RT", "Precursor_mz", "Precursor_Charge",
    "Parent_Fragment_ID", "Parent_Sequence", "Parent_Start", "Parent_End",
    "Modification_ID", "Modification_Name", "Modification_Target_Base",
    "Modification_Mass_Shift", "Parent_Unmodified_Mass", "Parent_Modified_Mass",
    "Parent_Charge", "Parent_Theoretical_mz", "Precursor_Error_Da",
    "Precursor_Error_ppm", "Modified_Precursor_Rescue", "Candidate_Rank", "Comment",
]


def annotate_ms2(
    mzml_path: str | None,
    theoretical_fragments: list[Fragment],
    config: Any,
    base_masses: dict[str, Any],
    modifications: list[Any] | None = None,
    warnings: list[dict[str, Any]] | None = None,
) -> dict[str, list[dict[str, Any]]]:
    ms2_config = getattr(config, "ms2_annotation", {}) or {}
    if not _as_bool(ms2_config.get("enabled"), True):
        return {}

    spectra = extract_ms2_spectra(mzml_path, ms2_config, warnings) if mzml_path else []
    ions = generate_theoretical_ms2_ions(theoretical_fragments, config, base_masses, warnings)
    matches, unmatched, spectrum_rows, parent_rows = match_ms2_spectra(spectra, ions, theoretical_fragments, config, modifications or [])
    evidence_rows = build_fragment_evidence(matches, theoretical_fragments, config, parent_rows)

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

    max_match_rows = _optional_positive_int(ms2_config.get("max_ms2_match_rows"), 100000)
    output_unmatched = _as_bool(ms2_config.get("output_unmatched_peaks"), True)
    max_unmatched = _optional_positive_int(ms2_config.get("max_unmatched_peaks"), 50000)
    match_rows = [ms2_match_row(match, config) for match in matches]
    if max_match_rows is not None:
        match_rows = sorted(match_rows, key=lambda row: (abs(float(row.get("Mass_Error_ppm") or 0)), -float(row.get("Observed_Intensity") or 0)))[:max_match_rows]
    unmatched_rows = _prioritize_unmatched(unmatched, max_unmatched) if output_unmatched else []

    results = {
        "MS2_Summary": build_ms2_summary(spectra, ions, matches, unmatched_rows, parent_rows, evidence_rows),
        "MS2_Spectra": spectrum_rows,
        "MS2_Parent_Candidates": parent_rows,
        "MS2_Modified_Precursor_Candidates": modified_precursor_rows(parent_rows),
        "MS2_Theoretical_Ions": [theoretical_ion_row(ion, config) for ion in ions],
        "MS2_Ion_Matches": match_rows,
        "MS2_Unmatched_Peaks": unmatched_rows,
        "MS2_Fragment_Evidence": evidence_rows,
    }
    if _as_bool(ms2_config.get("output_all_peak_annotations"), False):
        results["MS2_Peak_Annotations"] = match_rows + [_unmatched_match_row(row) for row in unmatched_rows]
    return results


def extract_ms2_spectra(
    mzml_path: str,
    ms2_config: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
) -> list[MS2SpectrumInfo]:
    spectra: list[MS2SpectrumInfo] = []
    min_intensity = float(ms2_config.get("min_peak_intensity", 10) or 0)
    min_relative_percent = float(ms2_config.get("min_relative_intensity_percent", 1.0) or 0)
    max_peaks = _optional_positive_int(ms2_config.get("max_peaks_per_spectrum"), 500)

    try:
        for scan_index, spectrum in enumerate(iter_spectra(mzml_path), start=1):
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

            raw_base_peak_index = int(np.argmax(intensity_array)) if intensity_array.size else None
            base_peak_mz = float(mz_array[raw_base_peak_index]) if raw_base_peak_index is not None else None
            base_peak_intensity = float(intensity_array[raw_base_peak_index]) if raw_base_peak_index is not None else None
            total_ion_current = float(np.sum(intensity_array)) if intensity_array.size else 0.0

            if intensity_array.size:
                mask = intensity_array >= min_intensity
                if base_peak_intensity and min_relative_percent > 0:
                    mask = mask & (intensity_array >= base_peak_intensity * min_relative_percent / 100.0)
                mz_array = mz_array[mask]
                intensity_array = intensity_array[mask]

            if max_peaks and mz_array.size > max_peaks:
                order = np.argsort(intensity_array)[::-1][:max_peaks]
                order = order[np.argsort(mz_array[order])]
                mz_array = mz_array[order]
                intensity_array = intensity_array[order]

            precursor = _precursor_info(spectrum)
            spectra.append(MS2SpectrumInfo(
                spectrum_id=str(spectrum.get("id") or f"scan_{scan_index}"),
                scan_index=scan_index,
                rt=_safe_float(_rt_minutes(spectrum)),
                precursor_mz=precursor.get("mz"),
                precursor_charge=precursor.get("charge"),
                precursor_intensity=precursor.get("intensity"),
                num_peaks=int(mz_array.size),
                base_peak_mz=base_peak_mz,
                base_peak_intensity=base_peak_intensity,
                total_ion_current=total_ion_current,
                peaks=[(float(mz), float(intensity)) for mz, intensity in zip(mz_array, intensity_array, strict=False)],
            ))
    except Exception as exc:
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
        for cut in range(1, len(sequence)):
            candidates = [("c", sequence[:cut], 1, cut), ("y", sequence[cut:], cut + 1, len(sequence))]
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
    theoretical_fragments: list[Fragment],
    config: Any,
    modifications: list[Any] | None = None,
) -> tuple[list[MS2IonMatch], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    ms2_config = getattr(config, "ms2_annotation", {}) or {}
    tolerance_ppm = float(ms2_config.get("mz_tolerance_ppm", 20) or 20)
    constrain_by_precursor = _as_bool(ms2_config.get("constrain_by_precursor"), True)
    fallback_to_all = _as_bool(ms2_config.get("fallback_to_all_ions_if_no_precursor_match"), False)

    matches: list[MS2IonMatch] = []
    unmatched_rows: list[dict[str, Any]] = []
    spectrum_rows: list[dict[str, Any]] = []
    parent_rows: list[dict[str, Any]] = []

    for spectrum in spectra:
        candidates = find_parent_candidates(spectrum, theoretical_fragments, config, modifications or [])
        parent_rows.extend(parent_candidate_rows(spectrum, candidates))
        if constrain_by_precursor:
            if candidates:
                candidate_parent_ids = {row["fragment"].fragment_id for row in candidates}
                scoped_ions = [ion for ion in ions if ion.parent_fragment_id in candidate_parent_ids]
                initial_status = "annotated"
            elif fallback_to_all:
                scoped_ions = list(ions)
                initial_status = "fallback_all_ions"
            else:
                scoped_ions = []
                initial_status = "no_precursor_parent_match"
        else:
            scoped_ions = list(ions)
            initial_status = "annotated"

        spectrum_match_count = 0
        spectrum_unmatched_count = 0
        for observed_mz, observed_intensity in spectrum.peaks:
            best, alternatives, nearest = _best_ion_match(observed_mz, scoped_ions, tolerance_ppm)
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
                    _possible_interpretation(nearest_error_ppm),
                    "No precursor-scoped theoretical c/y ion was within tolerance.",
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
                comment="precursor constrained" if candidates else "fallback all ions" if fallback_to_all else "",
            ))

        if not spectrum.peaks:
            status = "no_peaks"
        elif not scoped_ions and initial_status == "no_precursor_parent_match":
            status = "no_precursor_parent_match"
        elif not scoped_ions:
            status = "no_theoretical_ions"
        elif spectrum_match_count:
            status = "annotated"
        else:
            status = initial_status if initial_status != "annotated" else "skipped"
        spectrum_rows.append(spectrum_row(spectrum, spectrum_match_count, spectrum_unmatched_count, status, candidates))

    return matches, unmatched_rows, spectrum_rows, parent_rows


def find_parent_candidates(
    spectrum: MS2SpectrumInfo,
    theoretical_fragments: list[Fragment],
    config: Any,
    modifications: list[Any] | None = None,
) -> list[dict[str, Any]]:
    if spectrum.precursor_mz is None:
        return []
    ms2_config = getattr(config, "ms2_annotation", {}) or {}
    tolerance_ppm = float(ms2_config.get("precursor_match_tolerance_ppm", 20) or 20)
    polarity = str(getattr(config, "instrument", {}).get("polarity", "negative") or "negative").lower()
    charges = [abs(int(spectrum.precursor_charge))] if spectrum.precursor_charge else _ion_charges(ms2_config)
    candidates: list[dict[str, Any]] = []
    for fragment in theoretical_fragments or []:
        for charge in charges:
            theoretical_mz = theoretical_mz_from_mass(fragment.unmodified_mass, charge, polarity)
            error_ppm = ppm_error(float(spectrum.precursor_mz), theoretical_mz)
            if abs(error_ppm) <= tolerance_ppm:
                candidates.append({
                    "fragment": fragment,
                    "charge": charge,
                    "theoretical_mz": theoretical_mz,
                    "error_da": float(spectrum.precursor_mz) - theoretical_mz,
                    "error_ppm": error_ppm,
                    "candidate_type": "unmodified",
                    "modification_id": "", "modification_name": "", "modification_target_base": "",
                    "modification_mass_shift": 0.0, "unmodified_mass": float(fragment.unmodified_mass),
                    "modified_mass": float(fragment.unmodified_mass), "comment": "",
                })
    unmodified_found = bool(candidates)
    candidates.extend(find_modified_parent_candidates(spectrum, theoretical_fragments, modifications or [], config))
    candidates = sorted(candidates, key=lambda row: abs(row["error_ppm"]))
    limit = _optional_positive_int(ms2_config.get("modified_precursor_max_candidates_per_spectrum"), 20)
    if limit is not None:
        candidates = candidates[:limit]
    for rank, candidate in enumerate(candidates, start=1):
        candidate["rank"] = rank
        candidate["modified_rescue"] = candidate["candidate_type"] == "modified" and not unmodified_found
    return candidates


def parent_candidate_rows(spectrum: MS2SpectrumInfo, candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not candidates:
        return [{
            "Spectrum_ID": spectrum.spectrum_id,
            "Scan_Index": spectrum.scan_index,
            "RT": spectrum.rt,
            "Precursor_mz": spectrum.precursor_mz,
            "Precursor_Charge": spectrum.precursor_charge,
            "Candidate_Parent_Fragment_ID": "",
            "Candidate_Parent_Sequence": "",
            "Candidate_Parent_Start": "",
            "Candidate_Parent_End": "",
            "Candidate_Type": "", "Modification_ID": "", "Modification_Name": "",
            "Modification_Target_Base": "", "Modification_Mass_Shift": "",
            "Parent_Unmodified_Mass": "", "Parent_Modified_Mass": "",
            "Parent_Charge": "",
            "Parent_Theoretical_mz": "",
            "Precursor_Error_Da": "",
            "Precursor_Error_ppm": "",
            "Parent_Match_Status": "no_precursor_parent_match",
            "Modified_Precursor_Rescue": False, "Parent_Candidate_Rank": "",
            "Comment": "No theoretical fragment matched the precursor m/z within tolerance.",
        }]
    rows = []
    types = {candidate["candidate_type"] for candidate in candidates}
    match_status = "matched_both" if len(types) > 1 else "matched_modified" if "modified" in types else "matched_unmodified"
    for candidate in candidates:
        fragment = candidate["fragment"]
        rows.append({
            "Spectrum_ID": spectrum.spectrum_id,
            "Scan_Index": spectrum.scan_index,
            "RT": spectrum.rt,
            "Precursor_mz": spectrum.precursor_mz,
            "Precursor_Charge": spectrum.precursor_charge,
            "Candidate_Parent_Fragment_ID": fragment.fragment_id,
            "Candidate_Parent_Sequence": fragment.sequence,
            "Candidate_Parent_Start": fragment.start,
            "Candidate_Parent_End": fragment.end,
            "Candidate_Type": candidate["candidate_type"],
            "Modification_ID": candidate["modification_id"],
            "Modification_Name": candidate["modification_name"],
            "Modification_Target_Base": candidate["modification_target_base"],
            "Modification_Mass_Shift": candidate["modification_mass_shift"],
            "Parent_Unmodified_Mass": candidate["unmodified_mass"],
            "Parent_Modified_Mass": candidate["modified_mass"],
            "Parent_Charge": candidate["charge"],
            "Parent_Theoretical_mz": candidate["theoretical_mz"],
            "Precursor_Error_Da": candidate["error_da"],
            "Precursor_Error_ppm": candidate["error_ppm"],
            "Parent_Match_Status": match_status,
            "Modified_Precursor_Rescue": candidate["modified_rescue"],
            "Parent_Candidate_Rank": candidate["rank"],
            "Comment": candidate.get("comment", ""),
        })
    return rows


def modified_precursor_rows(parent_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for row in parent_rows:
        if row.get("Candidate_Type") != "modified":
            continue
        rows.append({
            "Spectrum_ID": row.get("Spectrum_ID"), "Scan_Index": row.get("Scan_Index"),
            "RT": row.get("RT"), "Precursor_mz": row.get("Precursor_mz"),
            "Precursor_Charge": row.get("Precursor_Charge"),
            "Parent_Fragment_ID": row.get("Candidate_Parent_Fragment_ID"),
            "Parent_Sequence": row.get("Candidate_Parent_Sequence"),
            "Parent_Start": row.get("Candidate_Parent_Start"), "Parent_End": row.get("Candidate_Parent_End"),
            "Modification_ID": row.get("Modification_ID"), "Modification_Name": row.get("Modification_Name"),
            "Modification_Target_Base": row.get("Modification_Target_Base"),
            "Modification_Mass_Shift": row.get("Modification_Mass_Shift"),
            "Parent_Unmodified_Mass": row.get("Parent_Unmodified_Mass"),
            "Parent_Modified_Mass": row.get("Parent_Modified_Mass"), "Parent_Charge": row.get("Parent_Charge"),
            "Parent_Theoretical_mz": row.get("Parent_Theoretical_mz"),
            "Precursor_Error_Da": row.get("Precursor_Error_Da"), "Precursor_Error_ppm": row.get("Precursor_Error_ppm"),
            "Modified_Precursor_Rescue": row.get("Modified_Precursor_Rescue"),
            "Candidate_Rank": row.get("Parent_Candidate_Rank"),
            "Comment": row.get("Comment") or "modified precursor candidate; fragment ions are unmodified MVP-5.2 approximation",
        })
    return rows


def _top_precursor_modifications(rows: list[dict[str, Any]]) -> str:
    counts: dict[tuple[str, str], int] = {}
    for row in rows:
        if row.get("Candidate_Type") != "modified":
            continue
        key = (str(row.get("Modification_ID") or ""), str(row.get("Modification_Name") or ""))
        counts[key] = counts.get(key, 0) + 1
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:10]
    return "; ".join(f"{mod_id}/{name}/{count}" for (mod_id, name), count in ranked)


def build_fragment_evidence(
    matches: list[MS2IonMatch], theoretical_fragments: list[Fragment], config: Any,
    parent_rows: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    fragment_lookup = {fragment.fragment_id: fragment for fragment in theoretical_fragments or []}
    grouped: dict[tuple[str, str], list[MS2IonMatch]] = {}
    for match in matches:
        grouped.setdefault((match.spectrum_id, match.parent_fragment_id), []).append(match)

    rows = []
    candidate_lookup: dict[tuple[str, str], dict[str, Any]] = {}
    for candidate in parent_rows or []:
        key = (str(candidate.get("Spectrum_ID") or ""), str(candidate.get("Candidate_Parent_Fragment_ID") or ""))
        if key[1] and (key not in candidate_lookup or int(candidate.get("Parent_Candidate_Rank") or 999999) < int(candidate_lookup[key].get("Parent_Candidate_Rank") or 999999)):
            candidate_lookup[key] = candidate
    for (spectrum_id, parent_id), group in grouped.items():
        first = group[0]
        fragment = fragment_lookup.get(parent_id)
        candidate = candidate_lookup.get((spectrum_id, parent_id), {})
        modified = candidate.get("Candidate_Type") == "modified"
        informative = [match for match in group if _is_informative_ion(match.best_ion_sequence, config)]
        ion_types = {match.best_ion_type for match in informative}
        abs_errors = [abs(match.mass_error_ppm) for match in group]
        coverage = _sequence_coverage(group, len(first.parent_fragment_sequence or ""))
        score = len(informative) * 2.0 + (len(group) - len(informative)) * 0.25
        if {"c", "y"} <= ion_types:
            score += 1.0
        level = _evidence_level(group, informative, ion_types)
        rows.append({
            "Spectrum_ID": spectrum_id,
            "RT": first.rt,
            "Precursor_mz": first.precursor_mz,
            "Parent_Fragment_ID": parent_id,
            "Parent_Sequence": first.parent_fragment_sequence,
            "Parent_Start": fragment.start if fragment else "",
            "Parent_End": fragment.end if fragment else "",
            "Parent_Charge": first.precursor_charge or "",
            "Candidate_Type": candidate.get("Candidate_Type", "unmodified"),
            "Modification_ID": candidate.get("Modification_ID", ""),
            "Modification_Name": candidate.get("Modification_Name", ""),
            "Modification_Mass_Shift": candidate.get("Modification_Mass_Shift", 0.0),
            "Parent_Modified_Mass": candidate.get("Parent_Modified_Mass", fragment.unmodified_mass if fragment else ""),
            "Modified_Precursor_Rescue": candidate.get("Modified_Precursor_Rescue", False),
            "Match_Approximation": "modified_parent_unmodified_ions_approximation" if modified else "unmodified_parent_ions",
            "Num_Matched_Ions": len(group),
            "Num_Informative_Ions": len(informative),
            "Num_c_Ions": sum(1 for match in group if match.best_ion_type == "c"),
            "Num_y_Ions": sum(1 for match in group if match.best_ion_type == "y"),
            "Best_Mass_Error_ppm": min(abs_errors) if abs_errors else "",
            "Median_Abs_Error_ppm": float(np.median(abs_errors)) if abs_errors else "",
            "Total_Matched_Intensity": sum(match.observed_intensity for match in group),
            "Sequence_Coverage": coverage,
            "Evidence_Score": score,
            "Evidence_Level": level,
            "Notes": ("modified precursor candidate; fragment ions are unmodified MVP-5.2 approximation" if modified else "") or ("1 nt ion dominated" if level == "Weak" and not informative else ""),
        })
    return sorted(rows, key=lambda row: (-float(row["Evidence_Score"]), row["Spectrum_ID"], row["Parent_Fragment_ID"]))


def build_ms2_summary(
    spectra: list[MS2SpectrumInfo],
    ions: list[TheoreticalMS2Ion],
    matches: list[MS2IonMatch],
    unmatched_rows: list[dict[str, Any]],
    parent_rows: list[dict[str, Any]],
    evidence_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    parent_counts: dict[str, int] = {}
    for match in matches:
        parent_counts[match.parent_fragment_id] = parent_counts.get(match.parent_fragment_id, 0) + 1
    best_parents = sorted(parent_counts.items(), key=lambda item: (-item[1], item[0]))[:10]
    matched_rows = [row for row in parent_rows if str(row.get("Parent_Match_Status", "")).startswith("matched_")]
    spectra_with_parent = {row["Spectrum_ID"] for row in matched_rows}
    spectra_without_parent = {row["Spectrum_ID"] for row in parent_rows if row.get("Parent_Match_Status") == "no_precursor_parent_match"}
    if not spectra:
        notes = "MS2 spectra: 0; MS2 annotation skipped."
    elif not ions:
        notes = "MS2 spectra were found, but no theoretical ions were generated."
    else:
        notes = "Precursor-constrained MS2 c/y ion annotation against theoretical digestion fragments."
    return [{
        "Total_MS2_Spectra": len(spectra),
        "Annotated_Spectra": len({match.spectrum_id for match in matches}),
        "Total_MS2_Peaks": sum(spectrum.num_peaks for spectrum in spectra),
        "Matched_MS2_Peaks": len(matches),
        "Unmatched_MS2_Peaks": len(unmatched_rows),
        "Total_Theoretical_Ions": len(ions),
        "Spectra_With_Precursor_Parent_Match": len(spectra_with_parent),
        "Spectra_Without_Precursor_Parent_Match": len(spectra_without_parent),
        "Spectra_With_Unmodified_Precursor_Match": len({row["Spectrum_ID"] for row in matched_rows if row.get("Candidate_Type") == "unmodified"}),
        "Spectra_With_Modified_Precursor_Match": len({row["Spectrum_ID"] for row in matched_rows if row.get("Candidate_Type") == "modified"}),
        "Spectra_Rescued_By_Modified_Precursor": len({row["Spectrum_ID"] for row in matched_rows if row.get("Modified_Precursor_Rescue")}),
        "Total_Modified_Parent_Candidates": sum(1 for row in matched_rows if row.get("Candidate_Type") == "modified"),
        "Top_Modification_Candidates_For_Precursors": _top_precursor_modifications(matched_rows),
        "Total_Parent_Candidates": len(matched_rows),
        "Total_Matched_Ion_Rows": len(matches),
        "Total_Unmatched_Output_Rows": len(unmatched_rows),
        "Strong_Evidence_Fragments": sum(1 for row in evidence_rows if row.get("Evidence_Level") == "Strong"),
        "Moderate_Evidence_Fragments": sum(1 for row in evidence_rows if row.get("Evidence_Level") == "Moderate"),
        "Weak_Evidence_Fragments": sum(1 for row in evidence_rows if row.get("Evidence_Level") == "Weak"),
        "Best_Matched_Parent_Fragments": "; ".join(f"{fragment_id} ({count})" for fragment_id, count in best_parents),
        "Notes": notes,
    }]


def spectrum_row(spectrum: MS2SpectrumInfo, matched: int, unmatched: int, status: str, candidates: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    candidates = candidates or []
    unmodified = sum(1 for item in candidates if item.get("candidate_type") == "unmodified")
    modified = sum(1 for item in candidates if item.get("candidate_type") == "modified")
    precursor_status = "no_precursor_mz" if spectrum.precursor_mz is None else "unmodified_and_modified_parent_match" if unmodified and modified else "unmodified_parent_match" if unmodified else "modified_parent_match" if modified else "no_precursor_parent_match"
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
        "Num_Parent_Candidates": len(candidates),
        "Num_Unmodified_Parent_Candidates": unmodified,
        "Num_Modified_Parent_Candidates": modified,
        "Precursor_Annotation_Status": precursor_status,
        "Comment": "" if status == "annotated" else "MS2 annotation skipped or produced no matches.",
    }


def theoretical_ion_row(ion: TheoreticalMS2Ion, config: Any) -> dict[str, Any]:
    return {
        "Ion_ID": ion.ion_id,
        "Parent_Fragment_ID": ion.parent_fragment_id,
        "Parent_Sequence": ion.parent_sequence,
        "Ion_Type": ion.ion_type,
        "Ion_Sequence": ion.ion_sequence,
        "Ion_Start": ion.ion_start,
        "Ion_End": ion.ion_end,
        "Ion_Length": len(ion.ion_sequence or ""),
        "Informative_Ion": _is_informative_ion(ion.ion_sequence, config),
        "Charge": ion.charge,
        "Theoretical_Mass": ion.theoretical_mass,
        "Theoretical_mz": ion.theoretical_mz,
        "Modification_ID": ion.modification_id,
        "Modification_Name": ion.modification_name,
        "Comment": ion.comment,
    }


def ms2_match_row(match: MS2IonMatch, config: Any) -> dict[str, Any]:
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
        "Ion_Length": len(match.best_ion_sequence or ""),
        "Informative_Ion": _is_informative_ion(match.best_ion_sequence, config),
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
        "Ion_Length": "",
        "Informative_Ion": False,
        "Parent_Fragment_ID": "",
        "Parent_Fragment_Sequence": "",
        "Ion_Charge": "",
        "Theoretical_mz": "",
        "Mass_Error_Da": "",
        "Mass_Error_ppm": "",
        "Match_Status": "outside_tolerance",
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


def _possible_interpretation(nearest_error_ppm: float | None) -> str:
    if nearest_error_ppm is None:
        return "unknown"
    if abs(nearest_error_ppm) <= 100:
        return "modified fragment ion candidate"
    return "outside tolerance"


def _is_informative_ion(sequence: str, config: Any) -> bool:
    ms2_config = getattr(config, "ms2_annotation", {}) or {}
    min_length = int(ms2_config.get("min_ion_length_for_evidence", 2) or 2)
    return len(sequence or "") >= min_length


def _evidence_level(group: list[MS2IonMatch], informative: list[MS2IonMatch], ion_types: set[str]) -> str:
    if len(informative) >= 3 and {"c", "y"} <= ion_types:
        return "Strong"
    if len(informative) >= 2:
        return "Moderate"
    if group:
        return "Weak"
    return "None"


def _sequence_coverage(group: list[MS2IonMatch], parent_length: int) -> float:
    if parent_length <= 0:
        return 0.0
    covered: set[int] = set()
    for match in group:
        ion_length = len(match.best_ion_sequence or "")
        if match.best_ion_type == "c":
            covered.update(range(1, ion_length + 1))
        elif match.best_ion_type == "y":
            covered.update(range(max(1, parent_length - ion_length + 1), parent_length + 1))
    return round(len(covered) / parent_length, 4)


def _prioritize_unmatched(rows: list[dict[str, Any]], max_rows: int | None) -> list[dict[str, Any]]:
    sorted_rows = sorted(
        rows,
        key=lambda row: (
            abs(float(row.get("Nearest_Error_ppm") or 1e9)),
            -float(row.get("Observed_Intensity") or 0),
        ),
    )
    return sorted_rows[:max_rows] if max_rows is not None else sorted_rows


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
