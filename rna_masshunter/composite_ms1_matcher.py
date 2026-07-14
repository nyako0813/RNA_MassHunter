"""Shadow-only MS1 matching for complete structure fragments."""
from __future__ import annotations
from collections import defaultdict
from typing import Any
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.ms1_mapping import _eligible_peaks, ppm_error

def _peak_value(peak: Any, name: str, tuple_index: int, default: Any = "") -> Any:
    if hasattr(peak, name):
        return getattr(peak, name)
    if isinstance(peak, (tuple, list)) and len(peak) > tuple_index:
        return peak[tuple_index]
    if isinstance(peak, dict):
        return peak.get(name, default)
    return default

def _physical_key(mz: Any, scan: Any, rt: Any) -> tuple:
    return (str(scan or ""), "" if rt in ("", None) else round(float(rt), 8), round(float(mz), 8))

def _mz_rt_key(mz: Any, rt: Any) -> tuple:
    return ("" if rt in ("", None) else round(float(rt), 8), round(float(mz), 8))

def match_composite_fragments_to_peaks(fragments: list[Any], peaks: list[Any], config: Any,
    *, legacy_matches: list[Any] = (), isomer_groups: dict[str, str] | None = None,
    audit_level: str = "full") -> list[dict[str, Any]]:
    mapping = getattr(config, "fragment_mapping", {}) or {}
    polarity = str(mapping.get("polarity") or "auto").lower()
    if polarity == "auto":
        polarity = str((getattr(config, "instrument", {}) or {}).get("polarity") or "negative").lower()
    min_charge = int(mapping.get("min_charge", 1) or 1)
    max_charge = int(mapping.get("max_charge", 8) or 8)
    tolerance = float(mapping.get("mz_tolerance_ppm", 10) or 10)
    mz_min = (getattr(config, "reconstruction", {}) or {}).get("mz_min")
    mz_max = (getattr(config, "reconstruction", {}) or {}).get("mz_max")
    eligible_peaks = _eligible_peaks(list(peaks or ()), mapping)
    legacy_keys = set()
    legacy_partial_keys = set()
    for match in legacy_matches or ():
        mz = getattr(match, "observed_mz", None)
        if mz is None and isinstance(match, dict):
            mz = match.get("observed_mz", match.get("Observed_mz", match.get("Observed_MZ")))
        scan = getattr(match, "scan_id", None)
        rt = getattr(match, "rt", None)
        if isinstance(match, dict):
            scan = match.get("scan_id", match.get("Scan_ID", scan))
            rt = match.get("rt", match.get("RT", rt))
        if mz not in (None, ""):
            if scan not in (None, ""):
                legacy_keys.add(_physical_key(mz, scan, rt))
            else:
                legacy_partial_keys.add(_mz_rt_key(mz, rt))
    rows: list[dict[str, Any]] = []
    matched_groups: dict[tuple, list[int]] = defaultdict(list)
    for fragment in fragments:
        for charge in range(min_charge, max_charge + 1):
            theoretical = mz_from_neutral_mass(fragment.neutral_exact_mass, charge, polarity)
            reason = ""
            if mz_min is not None and theoretical < float(mz_min):
                reason = "theoretical_mz_below_acquisition_range"
            elif mz_max is not None and theoretical > float(mz_max):
                reason = "theoretical_mz_above_acquisition_range"
            matches = []
            if not reason:
                for peak in eligible_peaks:
                    observed = float(_peak_value(peak, "mz", 0, 0) or 0)
                    error = ppm_error(observed, theoretical)
                    if abs(error) <= tolerance:
                        matches.append((abs(error), -float(_peak_value(peak, "intensity", 1, 0) or 0), peak, error))
            common = {
                "Candidate_ID": fragment.candidate_id, "Complete_Structure_ID": fragment.complete_structure_id,
                "Fragment_ID": fragment.fragment_id, "Fragment_Type": fragment.fragment_type,
                "Start_Position": fragment.start, "End_Position": fragment.end,
                "Included_Modified_Positions": ";".join(map(str, fragment.included_modified_positions)),
                "Included_Backbone_Bonds": ";".join(fragment.included_backbone_bonds),
                "Neutral_Exact_Mass": fragment.neutral_exact_mass, "Charge": charge,
                "Theoretical_mz": theoretical, "Audit_Level": audit_level,
                "Applied_To_Formal_Result": False, "Formal_Change_Ready": False,
            }
            if not matches:
                rows.append(dict(common, Observed_mz="", Mass_Error_Da="", Mass_Error_ppm="",
                    Observed_Intensity="", Observed_Scan="", Observed_RT="",
                    Match_Status="not_observable" if reason else "no_observation",
                    Support_Class="not_observable" if reason else "no_observation",
                    Not_Observable_Reason=reason, Legacy_Competition_Class="",
                    Is_Isomeric=bool((isomer_groups or {}).get(fragment.candidate_id))))
                continue
            matches.sort(key=lambda item: (item[0], item[1]))
            _, _, peak, error = matches[0]
            observed = float(_peak_value(peak, "mz", 0, 0) or 0)
            scan = _peak_value(peak, "scan_id", 3, "")
            rt = _peak_value(peak, "rt", 2, "")
            key = _physical_key(observed, scan, rt)
            candidate_specific = bool(fragment.included_modified_positions or fragment.included_backbone_modifications)
            support = ("observation_nondiscriminating" if not candidate_specific else
                       "isomeric_unresolved" if (isomer_groups or {}).get(fragment.candidate_id) else
                       "shared_with_legacy" if key in legacy_keys or _mz_rt_key(observed, rt) in legacy_partial_keys else
                       "unique_composite_support")
            row = dict(common, Observed_mz=observed, Mass_Error_Da=observed-theoretical,
                Mass_Error_ppm=error, Observed_Intensity=float(_peak_value(peak, "intensity", 1, 0) or 0),
                Observed_Scan=scan, Observed_RT=rt, Match_Status="matched",
                Support_Class=support, Not_Observable_Reason="",
                Legacy_Competition_Class=("OBSERVATION_NONDISCRIMINATING" if support == "observation_nondiscriminating" else
                    "BOTH_EQUIVALENT" if support == "shared_with_legacy" else "COMPOSITE_ONLY"),
                Is_Isomeric=bool((isomer_groups or {}).get(fragment.candidate_id)))
            matched_groups[key].append(len(rows)); rows.append(row)
    for indexes in matched_groups.values():
        candidate_ids = {rows[i]["Candidate_ID"] for i in indexes}
        if len(candidate_ids) > 1:
            for i in indexes:
                rows[i]["Support_Class"] = "shared_with_other_composite"
    return rows