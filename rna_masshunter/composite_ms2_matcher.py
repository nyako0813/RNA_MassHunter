"""Match complete-structure shadow ions to precursor-compatible MS2 spectra."""
from __future__ import annotations
from collections import defaultdict
from typing import Any
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.ms1_mapping import ppm_error

def match_composite_ms2(ions: list[dict[str, Any]], spectra: list[Any], config: Any,
    *, audit_level: str = "full") -> list[dict[str, Any]]:
    ms2 = getattr(config, "ms2_annotation", {}) or {}
    tolerance = float(ms2.get("mz_tolerance_ppm", 20) or 20)
    precursor_tolerance = float(ms2.get("precursor_match_tolerance_ppm", 20) or 20)
    constrain = bool(ms2.get("constrain_by_precursor", True))
    fallback = bool(ms2.get("fallback_to_all_ions_if_no_precursor_match", False))
    polarity = str((getattr(config, "instrument", {}) or {}).get("polarity") or "negative")
    rows: list[dict[str, Any]] = []
    peak_candidate_index: dict[tuple, list[int]] = defaultdict(list)
    for spectrum in spectra or ():
        spectrum_id = getattr(spectrum, "spectrum_id", "")
        precursor_mz = getattr(spectrum, "precursor_mz", None)
        precursor_charge = getattr(spectrum, "precursor_charge", None)
        eligible = list(ions)
        if constrain:
            if precursor_mz in (None, "") or precursor_charge in (None, "", 0):
                continue
            z = abs(int(precursor_charge))
            eligible = [
                ion for ion in ions
                if abs(ppm_error(float(precursor_mz), mz_from_neutral_mass(
                    float(ion["Parent_Neutral_Mass"]), z, polarity
                ))) <= precursor_tolerance
            ]
            if not eligible:
                if fallback:
                    eligible = list(ions)
                else:
                    continue
        for peak_index, (observed, intensity) in enumerate(getattr(spectrum, "peaks", ())):
            candidates = []
            for ion in eligible:
                error = ppm_error(float(observed), float(ion["Theoretical_mz"]))
                if abs(error) <= tolerance:
                    candidates.append((abs(error), ion, error))
            if not candidates:
                continue
            candidates.sort(key=lambda item: (item[0], item[1]["Candidate_ID"], item[1]["Ion_ID"]))
            _, ion, error = candidates[0]
            row = {
                "Candidate_ID": ion["Candidate_ID"], "Complete_Structure_ID": ion["Complete_Structure_ID"],
                "Spectrum_ID": spectrum_id, "Precursor_mz": precursor_mz,
                "Precursor_Charge": precursor_charge,
                "Ion_Series": ion["Ion_Series"], "Ion_Number": ion["Ion_Number"],
                "Cleavage_Position": ion["Cleavage_Position"], "Included_Positions": ion["Included_Positions"],
                "Included_Modified_Positions": ion["Included_Modified_Positions"],
                "Included_Backbone_Bonds": ion["Included_Backbone_Bonds"],
                "Theoretical_Neutral_Mass": ion["Theoretical_Neutral_Mass"],
                "Theoretical_mz": ion["Theoretical_mz"], "Observed_mz": observed,
                "Mass_Error_Da": float(observed)-float(ion["Theoretical_mz"]), "Mass_Error_ppm": error,
                "Observed_Intensity": intensity, "Position_Informative": ion["Position_Informative"],
                "Backbone_Informative": ion["Backbone_Informative"],
                "Candidate_Discriminating": len({item[1]["Candidate_ID"] for item in candidates}) == 1,
                "Isomer_Discriminating": False, "Legacy_Competition_Class": "OBSERVATION_NONDISCRIMINATING",
                "Audit_Level": audit_level, "Applied_To_Formal_Result": False, "Formal_Change_Ready": False,
            }
            key = (spectrum_id, peak_index)
            peak_candidate_index[key].append(len(rows)); rows.append(row)
    for indexes in peak_candidate_index.values():
        ids = {rows[i]["Candidate_ID"] for i in indexes}
        if len(ids) == 1:
            for i in indexes:
                rows[i]["Legacy_Competition_Class"] = "MS2_DISCRIMINATED"
    return rows