"""Modified c/y ion generation and localization evidence for MVP-5.3."""

from typing import Any

import numpy as np

from rna_masshunter.masses import calculate_unmodified_rna_mass
from rna_masshunter.ms1_mapping import ppm_error, theoretical_mz_from_mass


def _as_bool(value: Any, default: bool) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _positive_int(value: Any, default: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed > 0 else default


def _target_positions(sequence: str, target_bases: Any, require_target: bool, limit: int) -> list[int]:
    if isinstance(target_bases, str):
        bases = [item.strip().upper().replace("T", "U") for item in target_bases.split(",") if item.strip()]
    else:
        bases = [str(item).upper().replace("T", "U") for item in (target_bases or [])]
    known = [base for base in bases if base not in {"", "ANY", "UNKNOWN", "N", "*"}]
    if not known:
        return [] if require_target else list(range(1, min(len(sequence), limit) + 1))
    return [index for index, base in enumerate(sequence, start=1) if base in known][:limit]


def generate_modified_theoretical_ions(
    parent_rows: list[dict[str, Any]],
    config: Any,
    base_masses: dict[str, Any],
) -> list[dict[str, Any]]:
    ms2 = getattr(config, "ms2_annotation", {}) or {}
    if not _as_bool(ms2.get("include_modified_fragment_ions"), True):
        return []
    require_target = _as_bool(ms2.get("modified_fragment_require_target_base"), True)
    include_counterparts = _as_bool(ms2.get("modified_fragment_include_unmodified_counterparts"), True)
    position_limit = _positive_int(ms2.get("modified_fragment_max_positions_per_candidate"), 20)
    min_length = _positive_int(ms2.get("modified_fragment_min_ion_length"), 1)
    informative_length = _positive_int(ms2.get("modified_fragment_min_ion_length_for_localization"), 2)
    max_rows = _positive_int(ms2.get("modified_fragment_max_rows"), 100000)
    polarity = str(getattr(config, "instrument", {}).get("polarity", "negative") or "negative").lower()
    rows: list[dict[str, Any]] = []
    seen_parents: set[tuple[Any, ...]] = set()
    for parent in parent_rows:
        if parent.get("Candidate_Type") != "modified":
            continue
        sequence = str(parent.get("Candidate_Parent_Sequence") or "").upper().replace("T", "U")
        parent_key = (parent.get("Spectrum_ID"), parent.get("Candidate_Parent_Fragment_ID"), parent.get("Modification_ID"), parent.get("Parent_Charge"))
        if parent_key in seen_parents:
            continue
        seen_parents.add(parent_key)
        positions = _target_positions(sequence, parent.get("Modification_Target_Base"), require_target, position_limit)
        try:
            shift = float(parent.get("Modification_Mass_Shift") or 0.0)
            charge = abs(int(parent.get("Parent_Charge") or 1))
        except (TypeError, ValueError):
            continue
        for position in positions:
            for cut in range(1, len(sequence)):
                for ion_type, ion_sequence, ion_start, ion_end in (
                    ("c", sequence[:cut], 1, cut), ("y", sequence[cut:], cut + 1, len(sequence))
                ):
                    if len(ion_sequence) < min_length:
                        continue
                    contains = ion_start <= position <= ion_end
                    if not contains and not include_counterparts:
                        continue
                    unmodified_mass = calculate_unmodified_rna_mass(ion_sequence, base_masses, terminal_form="default")
                    if unmodified_mass is None:
                        continue
                    applied_shift = shift if contains else 0.0
                    mass = float(unmodified_mass) + applied_shift
                    rows.append({
                        "Ion_ID": f"M53ION_{len(rows) + 1:08d}",
                        "Spectrum_ID": parent.get("Spectrum_ID"),
                        "Parent_Fragment_ID": parent.get("Candidate_Parent_Fragment_ID"),
                        "Parent_Sequence": sequence, "Candidate_Type": "modified",
                        "Modification_ID": parent.get("Modification_ID"),
                        "Modification_Name": parent.get("Modification_Name"),
                        "Modification_Target_Base": parent.get("Modification_Target_Base"),
                        "Modification_Mass_Shift": shift,
                        "Candidate_Modification_Position_In_Parent": position,
                        "Candidate_Modification_Base": sequence[position - 1],
                        "Ion_Type": ion_type, "Ion_Sequence": ion_sequence,
                        "Ion_Start": ion_start, "Ion_End": ion_end, "Ion_Length": len(ion_sequence),
                        "Informative_Ion": len(ion_sequence) >= informative_length,
                        "Ion_Contains_Modification": contains,
                        "Modification_Mass_Shift_Applied": applied_shift,
                        "Charge": charge, "Theoretical_Mass": mass,
                        "Theoretical_mz": theoretical_mz_from_mass(mass, charge, polarity),
                        "Parent_Start": parent.get("Candidate_Parent_Start"),
                        "Parent_End": parent.get("Candidate_Parent_End"),
                        "Comment": "modified c/y ion" if contains else "unmodified counterpart for position candidate",
                    })
                    if len(rows) >= max_rows:
                        return rows
    return rows


def match_modified_ions(spectra: list[Any], ions: list[dict[str, Any]], config: Any) -> list[dict[str, Any]]:
    ms2 = getattr(config, "ms2_annotation", {}) or {}
    tolerance = float(ms2.get("mz_tolerance_ppm", 20) or 20)
    max_rows = _positive_int(ms2.get("modified_fragment_max_rows"), 100000)
    spectrum_lookup = {spectrum.spectrum_id: spectrum for spectrum in spectra}
    grouped: dict[tuple[str, str, str, int], list[dict[str, Any]]] = {}
    for ion in ions:
        key = (str(ion["Spectrum_ID"]), str(ion["Parent_Fragment_ID"]), str(ion["Modification_ID"]), int(ion["Candidate_Modification_Position_In_Parent"]))
        grouped.setdefault(key, []).append(ion)
    rows: list[dict[str, Any]] = []
    for key, position_ions in grouped.items():
        spectrum = spectrum_lookup.get(key[0])
        if spectrum is None:
            continue
        for observed_mz, intensity in spectrum.peaks:
            candidates = []
            for ion in position_ions:
                error_ppm = ppm_error(observed_mz, float(ion["Theoretical_mz"]))
                if abs(error_ppm) <= tolerance:
                    candidates.append((abs(error_ppm), error_ppm, ion))
            if not candidates:
                continue
            candidates.sort(key=lambda item: (item[0], not item[2]["Ion_Contains_Modification"], item[2]["Ion_ID"]))
            _, error_ppm, ion = candidates[0]
            status = "multiple_candidates" if len(candidates) > 1 else "matched_modified_ion" if ion["Ion_Contains_Modification"] else "matched_unmodified_counterpart"
            relative = float(intensity) / float(spectrum.base_peak_intensity or intensity or 1.0)
            confidence = "High" if abs(error_ppm) <= tolerance * 0.25 and relative >= 0.05 else "Medium" if abs(error_ppm) <= tolerance * 0.5 else "Low"
            rows.append({
                "Spectrum_ID": spectrum.spectrum_id, "Scan_Index": spectrum.scan_index, "RT": spectrum.rt,
                "Precursor_mz": spectrum.precursor_mz, "Precursor_Charge": spectrum.precursor_charge,
                "Parent_Fragment_ID": ion["Parent_Fragment_ID"], "Parent_Sequence": ion["Parent_Sequence"],
                "Candidate_Type": "modified", "Modification_ID": ion["Modification_ID"],
                "Modification_Name": ion["Modification_Name"],
                "Candidate_Modification_Position_In_Parent": ion["Candidate_Modification_Position_In_Parent"],
                "Candidate_Modification_Base": ion["Candidate_Modification_Base"],
                "Observed_mz": observed_mz, "Observed_Intensity": intensity, "Ion_ID": ion["Ion_ID"],
                "Ion_Type": ion["Ion_Type"], "Ion_Sequence": ion["Ion_Sequence"], "Ion_Start": ion["Ion_Start"],
                "Ion_End": ion["Ion_End"], "Ion_Length": ion["Ion_Length"], "Informative_Ion": ion["Informative_Ion"],
                "Ion_Contains_Modification": ion["Ion_Contains_Modification"],
                "Modification_Mass_Shift_Applied": ion["Modification_Mass_Shift_Applied"],
                "Theoretical_mz": ion["Theoretical_mz"], "Mass_Error_Da": observed_mz - ion["Theoretical_mz"],
                "Mass_Error_ppm": error_ppm, "Match_Status": status, "Confidence": confidence,
                "Comment": ion["Comment"],
            })
            if len(rows) >= max_rows:
                return annotate_position_discrimination(rows, ions, config)
    return annotate_position_discrimination(rows, ions, config)


def annotate_position_discrimination(
    matches: list[dict[str, Any]], ions: list[dict[str, Any]], config: Any,
) -> list[dict[str, Any]]:
    ms2 = getattr(config, "ms2_annotation", {}) or {}
    enabled = _as_bool(ms2.get("annotate_position_discriminating_ions"), True)
    positions_by_group: dict[tuple[str, str, str], set[int]] = {}
    starts_by_group: dict[tuple[str, str, str], Any] = {}
    for ion in ions:
        key = (str(ion.get("Spectrum_ID") or ""), str(ion.get("Parent_Fragment_ID") or ""), str(ion.get("Modification_ID") or ""))
        positions_by_group.setdefault(key, set()).add(int(ion["Candidate_Modification_Position_In_Parent"]))
        starts_by_group[key] = ion.get("Parent_Start")
    for row in matches:
        key = (str(row.get("Spectrum_ID") or ""), str(row.get("Parent_Fragment_ID") or ""), str(row.get("Modification_ID") or ""))
        positions = sorted(positions_by_group.get(key, set()))
        current = int(row.get("Candidate_Modification_Position_In_Parent") or 0)
        covered = [position for position in positions if int(row.get("Ion_Start") or 0) <= position <= int(row.get("Ion_End") or 0)]
        also = [position for position in covered if position != current]
        informative = bool(row.get("Informative_Ion"))
        contains = bool(row.get("Ion_Contains_Modification"))
        discriminating = enabled and len(positions) > 1 and contains and informative and covered == [current]
        parent_start = starts_by_group.get(key)
        try:
            trna_current = int(parent_start) + current - 1
            also_trna = [int(parent_start) + position - 1 for position in also]
        except (TypeError, ValueError):
            trna_current = current if discriminating else ""
            also_trna = also
        if len(positions) <= 1:
            reason = "single-candidate-position"
        elif not informative:
            reason = "low-information-1nt-ion"
        elif len(covered) > 1:
            reason = "ion-covers-multiple-candidate-positions"
        elif not discriminating:
            reason = "ion-does-not-separate-candidate-positions"
        else:
            reason = ""
        row["Position_Discriminating_Ion"] = discriminating
        row["Discriminates_Position"] = trna_current if discriminating else ""
        row["Also_Explains_Positions"] = ";".join(map(str, also_trna))
        row["Non_Discriminating_Reason"] = reason
    return matches


def build_localization_evidence(ions: list[dict[str, Any]], matches: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys: dict[tuple[str, str, str, int], dict[str, Any]] = {}
    for ion in ions:
        key = (str(ion["Spectrum_ID"]), str(ion["Parent_Fragment_ID"]), str(ion["Modification_ID"]), int(ion["Candidate_Modification_Position_In_Parent"]))
        keys.setdefault(key, ion)
    grouped: dict[tuple[str, str, str, int], list[dict[str, Any]]] = {key: [] for key in keys}
    for match in matches:
        key = (str(match["Spectrum_ID"]), str(match["Parent_Fragment_ID"]), str(match["Modification_ID"]), int(match["Candidate_Modification_Position_In_Parent"]))
        grouped.setdefault(key, []).append(match)
    rows = []
    for key, group in grouped.items():
        ion = keys[key]
        modified = [row for row in group if row["Ion_Contains_Modification"]]
        counterparts = [row for row in group if not row["Ion_Contains_Modification"]]
        informative = [row for row in modified if row["Informative_Ion"]]
        discriminating = [row for row in modified if row.get("Position_Discriminating_Ion")]
        informative_discriminating = [row for row in discriminating if row.get("Informative_Ion")]
        non_discriminating = [row for row in modified if not row.get("Position_Discriminating_Ion")]
        types = {row["Ion_Type"] for row in informative}
        errors = [abs(float(row["Mass_Error_ppm"])) for row in modified]
        if len(informative) >= 3 and {"c", "y"} <= types:
            level = "Strong"
        elif len(informative) >= 2:
            level = "Moderate"
        elif modified:
            level = "Weak"
        else:
            level = "None"
        discrimination_types = {row["Ion_Type"] for row in informative_discriminating}
        discrimination_cuts = {(row.get("Ion_Type"), row.get("Ion_Start"), row.get("Ion_End")) for row in informative_discriminating}
        if len(informative_discriminating) >= 2 and ({"c", "y"} <= discrimination_types or len(discrimination_cuts) >= 2):
            discrimination_level = "Strong"
        elif informative_discriminating:
            discrimination_level = "Moderate"
        elif modified:
            discrimination_level = "Weak"
        else:
            discrimination_level = "None"
        score = len(informative) * 2.0 + (len(modified) - len(informative)) * 0.25 + (1.0 if {"c", "y"} <= types else 0.0) - len(counterparts) * 0.1
        parent_start = ion.get("Parent_Start")
        try:
            trna_position = int(parent_start) + key[3] - 1
        except (TypeError, ValueError):
            trna_position = ""
        rows.append({
            "Spectrum_ID": key[0], "RT": next((row["RT"] for row in group), ""),
            "Precursor_mz": next((row["Precursor_mz"] for row in group), ""),
            "Parent_Fragment_ID": key[1], "Parent_Sequence": ion["Parent_Sequence"],
            "Parent_Start": parent_start, "Parent_End": ion.get("Parent_End"),
            "Modification_ID": key[2], "Modification_Name": ion["Modification_Name"],
            "Candidate_Modification_Position_In_Parent": key[3],
            "Candidate_Modification_Position_In_tRNA": trna_position,
            "Candidate_Modification_Base": ion["Candidate_Modification_Base"],
            "Num_Modified_Ion_Matches": len(modified), "Num_Unmodified_Counterpart_Matches": len(counterparts),
            "Num_Informative_Modified_Ion_Matches": len(informative),
            "Num_c_Modified_Ions": sum(row["Ion_Type"] == "c" for row in modified),
            "Num_y_Modified_Ions": sum(row["Ion_Type"] == "y" for row in modified),
            "Num_Position_Discriminating_Modified_Ions": len(discriminating),
            "Num_Informative_Position_Discriminating_Modified_Ions": len(informative_discriminating),
            "Num_Non_Discriminating_Modified_Ions": len(non_discriminating),
            "Has_Position_Discriminating_Evidence": bool(discriminating),
            "Position_Discrimination_Level": discrimination_level,
            "Best_Modified_Ion_Error_ppm": min(errors) if errors else "",
            "Median_Abs_Modified_Ion_Error_ppm": float(np.median(errors)) if errors else "",
            "Total_Modified_Ion_Intensity": sum(float(row["Observed_Intensity"]) for row in modified),
            "Localization_Score": score, "Localization_Level": level,
            "Localization_Interpretation": "modification-supported-on-position" if level in {"Strong", "Moderate"} else "insufficient-localization-evidence" if modified else "no-modified-ion-support",
            "Notes": "1 nt ions are low-information localization evidence" if modified and not informative else "",
        })
    ambiguity_groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for row in rows:
        ambiguity_groups.setdefault((row["Spectrum_ID"], row["Parent_Fragment_ID"], row["Modification_ID"]), []).append(row)
    for group in ambiguity_groups.values():
        supported = [row for row in group if row["Localization_Level"] in {"Strong", "Moderate"}]
        if len(supported) > 1:
            for row in supported:
                row["Localization_Interpretation"] = "ambiguous-multiple-positions"
        elif len(group) > 1 and not any(row["Has_Position_Discriminating_Evidence"] for row in group) and any(row["Num_Modified_Ion_Matches"] for row in group):
            for row in group:
                if row["Num_Modified_Ion_Matches"]:
                    row["Localization_Interpretation"] = "position-ambiguous-non-discriminating-ions"
    return sorted(rows, key=lambda row: (row["Spectrum_ID"], row["Parent_Fragment_ID"], -float(row["Localization_Score"]), row["Candidate_Modification_Position_In_Parent"]))
