"""Match complete-structure shadow ions to precursor-compatible MS2 spectra."""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha1
from typing import Any

from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.ms1_mapping import ppm_error
from rna_masshunter.ms2_identity_evidence import physical_observed_peak_key


def _assignment_sort_key(item: tuple[float, dict[str, Any], float]) -> tuple[Any, ...]:
    absolute_error, ion, _error = item
    return (
        absolute_error,
        str(ion.get("Candidate_ID") or ""),
        str(ion.get("Ion_ID") or ""),
        str(ion.get("Complete_Structure_ID") or ""),
        str(ion.get("Parent_Fragment_ID") or ""),
        float(ion.get("Theoretical_mz") or 0.0),
        str(ion.get("Included_Modified_Positions") or ""),
        str(ion.get("Included_Backbone_Bonds") or ""),
    )


def _intensity_state(value: Any) -> str:
    try:
        intensity = float(value)
    except (TypeError, ValueError):
        return "missing"
    if intensity > 0:
        return "positive"
    if intensity == 0:
        return "zero"
    return "negative"


def _joined(values: set[str]) -> str:
    return ";".join(sorted(value for value in values if value))


def _match_id(
    physical_key: str, candidate_id: str, structure_id: str, ion_id: str,
) -> str:
    source = "|".join((physical_key, candidate_id, structure_id, ion_id))
    return f"CMPMATCH_{sha1(source.encode('utf-8')).hexdigest()[:16].upper()}"


def _deduplicate_candidates(
    physical_key: str,
    candidates: list[tuple[float, dict[str, Any], float]],
) -> list[tuple[float, dict[str, Any], float]]:
    unique: dict[tuple[str, str, str, str], tuple[float, dict[str, Any], float]] = {}
    for item in sorted(candidates, key=_assignment_sort_key):
        ion = item[1]
        key = (
            physical_key,
            str(ion.get("Candidate_ID") or ""),
            str(ion.get("Complete_Structure_ID") or ""),
            str(ion.get("Ion_ID") or ""),
        )
        unique.setdefault(key, item)
    return sorted(unique.values(), key=_assignment_sort_key)


def match_composite_ms2(
    ions: list[dict[str, Any]], spectra: list[Any], config: Any,
    *, audit_level: str = "full", return_competition: bool = False,
) -> list[dict[str, Any]] | tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Return legacy-compatible best rows and optionally all tolerance assignments."""
    ms2 = getattr(config, "ms2_annotation", {}) or {}
    tolerance = float(ms2.get("mz_tolerance_ppm", 20) or 20)
    precursor_tolerance = float(ms2.get("precursor_match_tolerance_ppm", 20) or 20)
    constrain = bool(ms2.get("constrain_by_precursor", True))
    fallback = bool(ms2.get("fallback_to_all_ions_if_no_precursor_match", False))
    polarity = str((getattr(config, "instrument", {}) or {}).get("polarity") or "negative")
    rows: list[dict[str, Any]] = []
    competition_rows: list[dict[str, Any]] = []
    peak_candidate_index: dict[tuple[str, int], list[int]] = defaultdict(list)

    for spectrum in spectra or ():
        spectrum_id = str(getattr(spectrum, "spectrum_id", "") or "")
        scan_index = getattr(spectrum, "scan_index", "")
        rt = getattr(spectrum, "rt", None)
        precursor_mz = getattr(spectrum, "precursor_mz", None)
        precursor_charge = getattr(spectrum, "precursor_charge", None)
        eligible = list(ions)
        if constrain:
            if precursor_mz in (None, "") or precursor_charge in (None, "", 0):
                continue
            z = abs(int(precursor_charge))
            eligible = [
                ion for ion in ions
                if abs(ppm_error(
                    float(precursor_mz),
                    mz_from_neutral_mass(float(ion["Parent_Neutral_Mass"]), z, polarity),
                )) <= precursor_tolerance
            ]
            if not eligible:
                if fallback:
                    eligible = list(ions)
                else:
                    continue

        for peak_index, (observed, intensity) in enumerate(
            getattr(spectrum, "peaks", ())
        ):
            candidates: list[tuple[float, dict[str, Any], float]] = []
            for ion in eligible:
                error = ppm_error(float(observed), float(ion["Theoretical_mz"]))
                if abs(error) <= tolerance:
                    candidates.append((abs(error), ion, error))
            if not candidates:
                continue

            physical_key = physical_observed_peak_key({
                "Spectrum_ID": spectrum_id,
                "Observed_mz": observed,
                "Observed_Intensity": intensity,
                "RT": rt,
            })
            candidates = _deduplicate_candidates(physical_key, candidates)
            if not candidates:
                continue
            best_absolute_error, best_ion, best_error = candidates[0]
            candidate_ids = {
                str(item[1].get("Candidate_ID") or "") for item in candidates
            }
            structure_ids = {
                str(item[1].get("Complete_Structure_ID") or "") for item in candidates
            }
            ion_ids = {str(item[1].get("Ion_ID") or "") for item in candidates}
            position_signatures = {
                str(item[1].get("Included_Modified_Positions") or "")
                for item in candidates
            }
            backbone_signatures = {
                str(item[1].get("Included_Backbone_Bonds") or "")
                for item in candidates
            }
            second_error = candidates[1][0] if len(candidates) > 1 else None
            margin = (
                second_error - best_absolute_error
                if second_error is not None else None
            )
            candidate_discriminating = len(candidate_ids) == 1

            legacy_row = {
                "Candidate_ID": best_ion["Candidate_ID"],
                "Complete_Structure_ID": best_ion["Complete_Structure_ID"],
                "Spectrum_ID": spectrum_id,
                "Precursor_mz": precursor_mz,
                "Precursor_Charge": precursor_charge,
                "Ion_Series": best_ion["Ion_Series"],
                "Ion_Number": best_ion["Ion_Number"],
                "Cleavage_Position": best_ion["Cleavage_Position"],
                "Included_Positions": best_ion["Included_Positions"],
                "Included_Modified_Positions": best_ion["Included_Modified_Positions"],
                "Included_Backbone_Bonds": best_ion["Included_Backbone_Bonds"],
                "Theoretical_Neutral_Mass": best_ion["Theoretical_Neutral_Mass"],
                "Theoretical_mz": best_ion["Theoretical_mz"],
                "Observed_mz": observed,
                "Mass_Error_Da": float(observed) - float(best_ion["Theoretical_mz"]),
                "Mass_Error_ppm": best_error,
                "Observed_Intensity": intensity,
                "Position_Informative": best_ion["Position_Informative"],
                "Backbone_Informative": best_ion["Backbone_Informative"],
                "Candidate_Discriminating": candidate_discriminating,
                "Isomer_Discriminating": False,
                "Legacy_Competition_Class": "OBSERVATION_NONDISCRIMINATING",
                "Audit_Level": audit_level,
                "Applied_To_Formal_Result": False,
                "Formal_Change_Ready": False,
            }
            key = (spectrum_id, peak_index)
            peak_candidate_index[key].append(len(rows))
            rows.append(legacy_row)

            for rank, (_absolute_error, ion, error) in enumerate(candidates, 1):
                candidate_id = str(ion.get("Candidate_ID") or "")
                structure_id = str(ion.get("Complete_Structure_ID") or "")
                ion_id = str(ion.get("Ion_ID") or "")
                competition_rows.append({
                    "Composite_Match_ID": _match_id(
                        physical_key, candidate_id, structure_id, ion_id,
                    ),
                    "Physical_Observed_Peak_Key": physical_key,
                    "Observed_Peak_Index": peak_index,
                    "Raw_Peak_Index": "",
                    "Raw_Peak_Index_Missing_Reason": (
                        "annotation_input_does_not_retain_raw_peak_index"
                    ),
                    "Scan_Index": scan_index,
                    "Spectrum_ID": spectrum_id,
                    "RT": rt if rt is not None else "",
                    "Observed_mz": observed,
                    "Observed_Intensity": intensity,
                    "Observed_Intensity_State": _intensity_state(intensity),
                    "Ion_ID": ion_id,
                    "Candidate_ID": candidate_id,
                    "Complete_Structure_ID": structure_id,
                    "Parent_Fragment_ID": ion.get("Parent_Fragment_ID", ""),
                    "Ion_Series": ion.get("Ion_Series", ""),
                    "Ion_Number": ion.get("Ion_Number", ""),
                    "Cleavage_Position": ion.get("Cleavage_Position", ""),
                    "Charge": ion.get("Charge", ""),
                    "Theoretical_mz": ion.get("Theoretical_mz", ""),
                    "Mass_Error_Da": float(observed) - float(ion["Theoretical_mz"]),
                    "Mass_Error_ppm": error,
                    "Assignment_Rank": rank,
                    "Best_Assignment": rank == 1,
                    "Within_Tolerance_Assignment_Count": len(candidates),
                    "Competing_Candidate_Count": len(candidate_ids - {candidate_id}),
                    "Competing_Candidate_IDs": _joined(candidate_ids - {candidate_id}),
                    "Competing_Complete_Structure_Count": len(
                        structure_ids - {structure_id}
                    ),
                    "Competing_Complete_Structure_IDs": _joined(
                        structure_ids - {structure_id}
                    ),
                    "Competing_Theoretical_Ion_Count": len(ion_ids - {ion_id}),
                    "Competing_Ion_IDs": _joined(ion_ids - {ion_id}),
                    "Best_Error_ppm": best_absolute_error,
                    "Second_Best_Error_ppm": (
                        second_error if second_error is not None else ""
                    ),
                    "Best_vs_Second_Error_Margin_ppm": (
                        margin if margin is not None else ""
                    ),
                    "Candidate_Specific": len(candidate_ids) == 1,
                    "Complete_Structure_Specific": len(structure_ids) == 1,
                    "Theoretical_Ion_Specific": len(ion_ids) == 1,
                    "Position_Specific": (
                        len(position_signatures) == 1
                        and bool(str(ion.get("Included_Modified_Positions") or ""))
                    ),
                    "Backbone_Bond_Specific": (
                        len(backbone_signatures) == 1
                        and bool(str(ion.get("Included_Backbone_Bonds") or ""))
                    ),
                    "Included_Positions": ion.get("Included_Positions", ""),
                    "Included_Modified_Positions": ion.get(
                        "Included_Modified_Positions", ""
                    ),
                    "Included_Backbone_Bonds": ion.get(
                        "Included_Backbone_Bonds", ""
                    ),
                    "Position_Informative": ion.get("Position_Informative", False),
                    "Backbone_Informative": ion.get("Backbone_Informative", False),
                    "Candidate_Discriminating": candidate_discriminating,
                    "Isomer_Discriminating": False,
                    "Legacy_Competition_Class": "OBSERVATION_NONDISCRIMINATING",
                    "Audit_Level": audit_level,
                    "Applied_To_Formal_Result": False,
                    "Formal_Change_Ready": False,
                    "Formal_Result_Changed": False,
                })

    for indexes in peak_candidate_index.values():
        ids = {rows[index]["Candidate_ID"] for index in indexes}
        if len(ids) == 1:
            for index in indexes:
                rows[index]["Legacy_Competition_Class"] = "MS2_DISCRIMINATED"

    if return_competition:
        competition_rows.sort(key=lambda row: (
            str(row["Physical_Observed_Peak_Key"]),
            int(row["Assignment_Rank"]),
            str(row["Candidate_ID"]),
            str(row["Ion_ID"]),
        ))
        return rows, competition_rows
    return rows
