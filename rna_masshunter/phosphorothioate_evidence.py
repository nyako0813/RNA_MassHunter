"""Indexed MS1 matching, competition audit, and evidence classes for PT pairs."""
from __future__ import annotations
from bisect import bisect_left, bisect_right
from collections import defaultdict
from typing import Any
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.ms1_mapping import _eligible_peaks, ppm_error

FALSE_FLAGS = {"Applied_To_Formal_Result": False, "Formal_Change_Ready": False, "Formal_Result_Changed": False}


def _value(item: Any, name: str, default: Any = "") -> Any:
    if hasattr(item, name): return getattr(item, name)
    if isinstance(item, dict):
        variants = (name, name.lower(), name.title(), name.replace("_", " "))
        for key in variants:
            if key in item: return item[key]
        for key in item:
            if str(key).lower() == name.lower(): return item[key]
    return default


def physical_peak_id(mz: Any, scan: Any, rt: Any, spectrum_id: Any = "") -> str:
    rt_text = "" if rt in (None, "") else f"{float(rt):.8f}"
    return f"PK|{spectrum_id or ''}|{scan or ''}|{rt_text}|{float(mz):.8f}"


def _external_keys(rows: list[Any] | tuple[Any, ...]) -> tuple[set[str], set[tuple[float, float]]]:
    full = set(); partial = set()
    for row in rows or ():
        mz = _value(row, "observed_mz", _value(row, "Observed_mz", "")); rt = _value(row, "rt", _value(row, "Observed_RT", ""))
        scan = _value(row, "scan_id", _value(row, "Observed_Scan", ""))
        if mz in (None, ""): continue
        if scan not in (None, ""): full.add(physical_peak_id(mz, scan, rt))
        else: partial.add((round(float(mz), 8), round(float(rt), 8) if rt not in (None, "") else 0.0))
    return full, partial


def _peak_record(peak: Any) -> dict[str, Any]:
    mz = float(_value(peak, "mz", 0) or 0); intensity = float(_value(peak, "intensity", 0) or 0)
    scan = _value(peak, "scan_id", ""); rt = _value(peak, "rt", "")
    return {"mz": mz, "intensity": intensity, "scan": scan, "rt": rt,
        "physical_id": physical_peak_id(mz, scan, rt), "tier": _value(peak, "tier", "")}


def _nearest(sorted_peaks: list[dict[str, Any]], mz_values: list[float], theoretical: float) -> dict[str, Any] | None:
    if not sorted_peaks: return None
    index = bisect_left(mz_values, theoretical); choices = []
    if index < len(sorted_peaks): choices.append(sorted_peaks[index])
    if index: choices.append(sorted_peaks[index - 1])
    return min(choices, key=lambda x: (abs(x["mz"] - theoretical), -x["intensity"]))


def _matches(sorted_peaks: list[dict[str, Any]], mz_values: list[float], theoretical: float,
    tolerance_ppm: float) -> list[dict[str, Any]]:
    width = abs(theoretical) * tolerance_ppm / 1_000_000
    return sorted_peaks[bisect_left(mz_values, theoretical - width):bisect_right(mz_values, theoretical + width)]


def build_pt_evidence(pairs: list[Any], peaks: list[Any], config: Any, *, legacy_matches: list[Any] = (),
    other_composite_matches: list[Any] = (), audit_level: str = "full", include_detail: bool = True,
    include_compact_states: bool = False) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    mapping = getattr(config, "fragment_mapping", {}) or {}; reconstruction = getattr(config, "reconstruction", {}) or {}
    polarity = str(mapping.get("polarity") or "auto").lower()
    if polarity == "auto": polarity = str((getattr(config, "instrument", {}) or {}).get("polarity") or "negative").lower()
    min_charge = int(mapping.get("min_charge", 1) or 1); max_charge = int(mapping.get("max_charge", 8) or 8)
    tolerance = float(mapping.get("mz_tolerance_ppm", 10) or 10); mz_min = reconstruction.get("mz_min"); mz_max = reconstruction.get("mz_max")
    eligible = [_peak_record(p) for p in _eligible_peaks(list(peaks or ()), mapping)]
    eligible.sort(key=lambda x: x["mz"]); mz_values = [x["mz"] for x in eligible]
    observations = []
    for pair_index, pair in enumerate(pairs):
        has_shared = bool(pair.spec.position_states)
        hypothesis_label = "" if pair.spec.search_mode != "hypothesis_driven" else ("H2" if has_shared else "H1")
        modified_label = "" if pair.spec.search_mode != "hypothesis_driven" else ("H4" if has_shared else "H3")
        for backbone, fragment, state_label in (("normal_phosphate", pair.normal_fragment, hypothesis_label),
                ("phosphorothioate", pair.modified_fragment, modified_label)):
            state_id = state_label or f"{pair.spec.candidate_id}|{backbone}"
            for charge in range(min_charge, max_charge + 1):
                theoretical = mz_from_neutral_mass(fragment.neutral_exact_mass, charge, polarity)
                observable = ((mz_min is None or theoretical >= float(mz_min))
                    and (mz_max is None or theoretical <= float(mz_max)) and bool(eligible))
                candidates = _matches(eligible, mz_values, theoretical, tolerance) if observable else []
                candidates.sort(key=lambda x: (abs(ppm_error(x["mz"], theoretical)), -x["intensity"], x["physical_id"]))
                best = candidates[0] if candidates else None; nearest = _nearest(eligible, mz_values, theoretical) if observable else None
                observations.append({"pair_index": pair_index, "pair": pair, "Candidate_State_ID": state_id,
                    "Backbone_State": backbone, "fragment": fragment, "Charge": charge,
                    "Theoretical_mz": theoretical, "Observable": observable, "Matched": bool(best),
                    "Explanation_ID": "|".join((pair.spec.bond_id, f"{fragment.start}-{fragment.end}",
                        ";".join(str(x.position) + ":" + "+".join(x.transform_ids) for x in pair.spec.position_states) or "unmodified",
                        fragment.elemental_composition_canonical, backbone, f"z{charge}")),
                    "match": best, "nearest": nearest, "Mass_Error_Da": (best["mz"] - theoretical) if best else "",
                    "Mass_Error_ppm": ppm_error(best["mz"], theoretical) if best else ""})
    explanations: dict[str, set[str]] = defaultdict(set)
    for obs in observations:
        theoretical = float(obs["Theoretical_mz"])
        for isotope in (-2, -1, 0, 1, 2):
            shifted = theoretical + isotope * 1.00335483507 / abs(int(obs["Charge"]))
            for peak in _matches(eligible, mz_values, shifted, tolerance):
                explanations[peak["physical_id"]].add(f"{obs['Explanation_ID']}|iso{isotope:+d}")
    legacy_full, legacy_partial = _external_keys(list(legacy_matches or ()) + list(other_composite_matches or ()))
    for obs in observations:
        peak = obs["match"]; own = f"{obs['Explanation_ID']}|iso+0"
        competitors = set()
        if peak:
            competitors = explanations[peak["physical_id"]] - {own}
            partial = (round(peak["mz"], 8), round(float(peak["rt"]), 8) if peak["rt"] not in (None, "") else 0.0)
            if peak["physical_id"] in legacy_full or partial in legacy_partial: competitors.add("external_legacy_or_composite")
        obs["competitors"] = sorted(competitors); obs["Competition_Count"] = len(competitors)
    state_rows = []
    for obs in observations:
        pair = obs["pair"]; fragment = obs["fragment"]; nearest = obs["nearest"]; match = obs["match"]
        if not include_detail:
            if include_compact_states:
                state_rows.append({"Candidate_ID":pair.spec.candidate_id,"Hypothesis_ID":pair.spec.hypothesis_id,"Search_Mode":pair.spec.search_mode,
                    "Sequence_ID":pair.spec.sequence_id,"Enzyme":pair.spec.enzyme,"Bond_ID":pair.spec.bond_id,
                    "Shared_Nucleoside_States":";".join(f"{x.position}:{'+'.join(x.transform_ids) or 'unmodified'}" for x in pair.spec.position_states) or "unmodified",
                    "Shared_Modified_Positions":";".join(str(x.position) for x in pair.spec.position_states),"Backbone_State":obs["Backbone_State"],
                    "Fragment_Start":fragment.start,"Fragment_End":fragment.end,"Terminal_Form":fragment.terminal_form,"Charge":obs["Charge"],
                    "Theoretical_mz":obs["Theoretical_mz"],"Observable":obs["Observable"],"Matched":obs["Matched"],"Competition_Count":obs["Competition_Count"],**FALSE_FLAGS})
            continue
        state_rows.append({"Audit_Level": audit_level, "Candidate_ID": pair.spec.candidate_id,
            "Hypothesis_ID": pair.spec.hypothesis_id, "Search_Mode": pair.spec.search_mode,
            "Candidate_State_ID": obs["Candidate_State_ID"], "Sequence_ID": pair.spec.sequence_id,
            "Enzyme": pair.spec.enzyme, "Bond_ID": pair.spec.bond_id,
            "Shared_Nucleoside_States": ";".join(f"{x.position}:{'+'.join(x.transform_ids) or 'unmodified'}" for x in pair.spec.position_states) or "unmodified",
            "Shared_Modified_Positions": ";".join(str(x.position) for x in pair.spec.position_states),
            "Applied_Transformations": ";".join(t for x in pair.spec.position_states for t in x.transform_ids),
            "Backbone_State": obs["Backbone_State"], "Fragment_Start": fragment.start, "Fragment_End": fragment.end,
            "Terminal_Form": fragment.terminal_form, "Charge": obs["Charge"],
            "Elemental_Composition": fragment.elemental_composition_canonical, "Neutral_Mass": fragment.neutral_exact_mass,
            "Theoretical_mz": obs["Theoretical_mz"], "Nearest_Observed_mz": nearest["mz"] if nearest else "",
            "Mass_Error_Da": ((match or nearest)["mz"] - obs["Theoretical_mz"]) if (match or nearest) else "",
            "Mass_Error_ppm": ppm_error((match or nearest)["mz"], obs["Theoretical_mz"]) if (match or nearest) else "",
            "Intensity": (match or nearest)["intensity"] if (match or nearest) else "",
            "Scan": (match or nearest)["scan"] if (match or nearest) else "",
            "RT": (match or nearest)["rt"] if (match or nearest) else "",
            "Physical_Peak_ID": (match or nearest)["physical_id"] if (match or nearest) else "",
            "Competition_Count": obs["Competition_Count"], "Competing_Candidate_IDs": ";".join(obs["competitors"]),
            "Observable": obs["Observable"], "Matched": obs["Matched"], **FALSE_FLAGS})
    by_pair_charge = defaultdict(dict)
    for obs in observations: by_pair_charge[(obs["pair_index"], obs["Charge"])][obs["Backbone_State"]] = obs
    evidence_rows = []
    for (pair_index, charge), group in sorted(by_pair_charge.items()):
        pair = pairs[pair_index]; normal = group["normal_phosphate"]; modified = group["phosphorothioate"]
        npeak = normal["match"]; mpeak = modified["match"]; same = bool(npeak and mpeak and npeak["physical_id"] == mpeak["physical_id"])
        consistent = (pair.composition_delta == pair.expected_backbone_delta and abs(pair.delta_consistency_error) <= 1e-9)
        ambiguous = bool((npeak and normal["Competition_Count"]) or (mpeak and modified["Competition_Count"]) or same)
        candidate_specific = bool(mpeak and not modified["Competition_Count"] and not same)
        if not consistent: evidence_class = "MASS_SHIFT_INCONSISTENT"; reason = "normal/PT exact composition or mass delta is inconsistent"
        elif not normal["Observable"] and not modified["Observable"]: evidence_class = "NOT_OBSERVABLE"; reason = "both counterpart m/z values are outside acquisition/filter evaluation"
        elif ambiguous and (npeak or mpeak): evidence_class = "AMBIGUOUS_PEAK_ASSIGNMENT"; reason = "matched physical peak has isotope/charge/fragment/legacy competition"
        elif npeak and mpeak: evidence_class = "BOTH_PRESENT"; reason = "normal and PT counterpart peaks are both present"
        elif mpeak and pair.cleavage_candidate.is_normal_cleavage_site and pair.block_status == "blocked" and candidate_specific:
            evidence_class = "PT_CANDIDATE_SPECIFIC_MS1_SUPPORT"
            reason = ("candidate-specific PT mass match with exact O-to-S delta, but no observed normal/PT "
                "peak pair, bond-localizing MS2, or control evidence")
        elif mpeak: evidence_class = "PT_ONLY_SUPPORT"; reason = "PT counterpart matched but strong paired criteria were not all satisfied"
        elif npeak: evidence_class = "NORMAL_ONLY_SUPPORT"; reason = "normal missed-cleavage counterpart matched without PT counterpart"
        else: evidence_class = "NEITHER_PRESENT"; reason = "both counterparts were observable but neither matched"
        neutral_delta = pair.modified_fragment.neutral_exact_mass - pair.normal_fragment.neutral_exact_mass
        theoretical_delta = modified["Theoretical_mz"] - normal["Theoretical_mz"]
        peak_status = "same_physical_peak" if same else "competing" if ambiguous else "distinct_or_unmatched"
        evidence_row = {"Audit_Level": audit_level, "Candidate_ID": pair.spec.candidate_id,
            "Hypothesis_ID": pair.spec.hypothesis_id, "Search_Mode": pair.spec.search_mode,
            "Sequence_ID": pair.spec.sequence_id, "Enzyme": pair.spec.enzyme, "Bond_ID": pair.spec.bond_id,
            "Left_Position": pair.cleavage_candidate.left_position, "Right_Position": pair.cleavage_candidate.right_position,
            "Left_Base": pair.cleavage_candidate.left_base, "Right_Base": pair.cleavage_candidate.right_base,
            "Is_Normal_Cleavage_Site": pair.cleavage_candidate.is_normal_cleavage_site,
            "Fragment_ID": f"{pair.spec.candidate_id}|{pair.spec.fragment_start}_{pair.spec.fragment_end}",
            "Fragment_Start": pair.spec.fragment_start, "Fragment_End": pair.spec.fragment_end,
            "Terminal_Form": pair.spec.terminal_form, "Charge": charge,
            "Shared_Nucleoside_States": ";".join(f"{x.position}:{'+'.join(x.transform_ids)}" for x in pair.spec.position_states) or "unmodified",
            "Shared_Modified_Positions": ";".join(str(x.position) for x in pair.spec.position_states),
            "Applied_Transformations": ";".join(t for x in pair.spec.position_states for t in x.transform_ids),
            "Normal_Backbone_State": "normal_phosphate", "Modified_Backbone_State": "phosphorothioate",
            "Normal_Composition": pair.normal_fragment.elemental_composition_canonical,
            "Modified_Composition": pair.modified_fragment.elemental_composition_canonical,
            "Composition_Delta": pair.composition_delta.canonical_string(),
            "Shared_Modification_Composition_Delta": pair.shared_modification_composition.canonical_string(),
            "Backbone_Composition_Delta": pair.expected_backbone_delta.canonical_string(),
            "Normal_Neutral_Mass": pair.normal_fragment.neutral_exact_mass,
            "Modified_Neutral_Mass": pair.modified_fragment.neutral_exact_mass,
            "Neutral_Mass_Delta": neutral_delta, "Expected_O_to_S_Delta": pair.expected_backbone_delta.exact_mass,
            "Delta_Consistency_Error": pair.delta_consistency_error,
            "Normal_Theoretical_mz": normal["Theoretical_mz"], "Modified_Theoretical_mz": modified["Theoretical_mz"],
            "Theoretical_mz_Delta": theoretical_delta,
            "Normal_Observed_mz": npeak["mz"] if npeak else "", "Modified_Observed_mz": mpeak["mz"] if mpeak else "",
            "Normal_Mass_Error_Da": normal["Mass_Error_Da"], "Modified_Mass_Error_Da": modified["Mass_Error_Da"],
            "Normal_Mass_Error_ppm": normal["Mass_Error_ppm"], "Modified_Mass_Error_ppm": modified["Mass_Error_ppm"],
            "Normal_Intensity": npeak["intensity"] if npeak else "", "Modified_Intensity": mpeak["intensity"] if mpeak else "",
            "Normal_Scan": npeak["scan"] if npeak else "", "Modified_Scan": mpeak["scan"] if mpeak else "",
            "Normal_RT": npeak["rt"] if npeak else "", "Modified_RT": mpeak["rt"] if mpeak else "",
            "Normal_Physical_Peak_ID": npeak["physical_id"] if npeak else "",
            "Modified_Physical_Peak_ID": mpeak["physical_id"] if mpeak else "", "Same_Physical_Peak": same,
            "Normal_Competition_Count": normal["Competition_Count"], "Modified_Competition_Count": modified["Competition_Count"],
            "Normal_Competing_Candidate_IDs": ";".join(normal["competitors"]),
            "Modified_Competing_Candidate_IDs": ";".join(modified["competitors"]),
            "Peak_Ambiguity_Status": peak_status,
            "Normal_Cleavage_Mechanism": "stochastic_missed_cleavage",
            "Modified_Cleavage_Mechanism": "phosphorothioate_blocked" if pair.block_status == "blocked" else "phosphorothioate_blocking_unknown",
            "Nucleoside_Blocking_Status": "unknown",
            "Cleavage_Status": "phosphorothioate_blocked" if pair.block_status == "blocked" else "unknown_blocked",
            "Evidence_Class": evidence_class, "Evidence_Reason": reason,
            "Mechanism_Discriminating": evidence_class == "PT_STRONG_PAIRED_SUPPORT",
            "Position_Localizing": False, "Candidate_Specific": candidate_specific,
            "Observable": normal["Observable"] or modified["Observable"], **FALSE_FLAGS}
        if include_detail:
            evidence_rows.append(evidence_row)
        else:
            summary_keys = ("Search_Mode", "Candidate_ID", "Evidence_Class", "Observable",
                "Candidate_Specific", "Modified_Physical_Peak_ID",
                "Normal_Competition_Count", "Modified_Competition_Count")
            evidence_rows.append({key: evidence_row.get(key) for key in summary_keys})
    return evidence_rows, state_rows
