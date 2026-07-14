"""Observation matching for phosphorothioate-blocked cleavage fragments."""
from __future__ import annotations
from typing import Any
from rna_masshunter.cleavage_constraints import evaluate_cleavage
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.ms1_mapping import ppm_error
from rna_masshunter.structure_fragment import extract_fragment_from_structure

def match_blocked_cleavage_fragments(structures: list[Any], sequence: str, peaks: list[Any],
    config: Any, backbone_transform: Any, base_masses: dict, formal_fragments: list[Any] = (),
    *, audit_level: str = "full") -> list[dict[str, Any]]:
    mapping = getattr(config, "fragment_mapping", {}) or {}
    polarity = str(mapping.get("polarity") or "auto")
    if polarity == "auto": polarity = str((getattr(config, "instrument", {}) or {}).get("polarity") or "negative")
    min_charge = int(mapping.get("min_charge", 1) or 1); max_charge = int(mapping.get("max_charge", 8) or 8)
    tolerance = float(mapping.get("mz_tolerance_ppm", 10) or 10)
    rows = []
    for structure in structures:
        for enzyme in ("RNase_T1", "Nuclease_P1"):
            result = evaluate_cleavage(sequence, enzyme, structure.bond_states, backbone_transform, base_masses)
            for shadow in result.fragments:
                if shadow.cleavage_origin != "phosphorothioate_blocked":
                    continue
                fragment = extract_fragment_from_structure(structure, sequence, shadow.start, shadow.end,
                    fragment_id=shadow.fragment_id, fragment_type=f"{enzyme}_blocked", terminal_form=shadow.terminal_form)
                alternative = any(int(item.start) == shadow.start and int(item.end) == shadow.end for item in formal_fragments)
                best = None
                for charge in range(min_charge, max_charge + 1):
                    theoretical = mz_from_neutral_mass(fragment.neutral_exact_mass, charge, polarity)
                    for peak in peaks or ():
                        observed = float(getattr(peak, "mz", peak[0] if isinstance(peak, (tuple, list)) else 0))
                        error = ppm_error(observed, theoretical)
                        if abs(error) <= tolerance:
                            key = (abs(error), -float(getattr(peak, "intensity", peak[1] if isinstance(peak, (tuple, list)) else 0) or 0))
                            if best is None or key < best[0]:
                                best = (key, charge, theoretical, observed, error, peak)
                blocked_bond = shadow.blocked_cleavage_bond_ids[0] if shadow.blocked_cleavage_bond_ids else ""
                blocked_position = shadow.blocked_cleavage_positions[0] if shadow.blocked_cleavage_positions else ""
                row = {
                    "Candidate_ID": structure.candidate_id, "Enzyme": enzyme, "Fragment_ID": shadow.fragment_id,
                    "Start_Position": shadow.start, "End_Position": shadow.end, "Blocked_Bond_ID": blocked_bond,
                    "Blocked_Cleavage_Position": blocked_position, "Cleavage_Status": "phosphorothioate_blocked",
                    "Blocked_Cleavage_Reason": ";".join(shadow.blocked_cleavage_reasons),
                    "Contains_Phosphorothioate": shadow.contains_phosphorothioate,
                    "Backbone_Modification_Count": shadow.backbone_modification_count,
                    "Backbone_Modification_Positions": ";".join(item.split("_")[0] for item in shadow.blocked_cleavage_bond_ids),
                    "Backbone_Mass_Delta": shadow.backbone_mass_delta,
                    "Theoretical_mz": best[2] if best else "",
                    "Observed_mz": best[3] if best else "", "Mass_Error_ppm": best[4] if best else "",
                    "Observed_Intensity": (getattr(best[5], "intensity", best[5][1] if isinstance(best[5], (tuple, list)) else 0) if best else ""),
                    "Alternative_Stochastic_Fragment_Exists": alternative,
                    "Mechanism_Discriminating": bool(best) and not alternative,
                    "Audit_Level": audit_level, "Applied_To_Formal_Result": False, "Formal_Change_Ready": False,
                }
                rows.append(row)
    return rows