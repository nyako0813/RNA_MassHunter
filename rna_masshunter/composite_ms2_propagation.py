"""Propagate complete position and bond state into c/y shadow ions."""
from __future__ import annotations
from typing import Any
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.structure_fragment import extract_fragment_from_structure

def generate_composite_theoretical_ions(structures: list[Any], fragments: list[Any],
    sequence: str, config: Any, *, audit_level: str = "full") -> list[dict[str, Any]]:
    ms2 = getattr(config, "ms2_annotation", {}) or {}
    polarity = str((getattr(config, "instrument", {}) or {}).get("polarity") or "negative")
    min_length = int(ms2.get("modified_fragment_min_ion_length", 1) or 1)
    informative_length = int(ms2.get("modified_fragment_min_ion_length_for_localization", 2) or 2)
    charge = abs(int(ms2.get("default_charge", 1) or 1))
    max_rows = int(ms2.get("modified_fragment_max_rows", 100000) or 100000)
    rows: list[dict[str, Any]] = []
    for structure in structures:
        for parent in fragments:
            if not any(int(parent.start) <= p <= int(parent.end) for p in structure.position_states) and not any(
                int(parent.start) <= b.left_position and b.right_position <= int(parent.end)
                for b in structure.bond_states.values()):
                continue
            parent_sequence = str(parent.sequence)
            parent_state_fragment = extract_fragment_from_structure(
                structure, sequence, int(parent.start), int(parent.end),
                fragment_id=str(parent.fragment_id), fragment_type=str(parent.enzyme),
                terminal_form=str(parent.terminal_form),
            )
            for cut in range(1, len(parent_sequence)):
                for ion_series, ion_start, ion_end in (("c", 1, cut), ("y", cut + 1, len(parent_sequence))):
                    ion_sequence = parent_sequence[ion_start-1:ion_end]
                    if len(ion_sequence) < min_length:
                        continue
                    absolute_start = int(parent.start) + ion_start - 1
                    absolute_end = int(parent.start) + ion_end - 1
                    state_fragment = extract_fragment_from_structure(structure, sequence, absolute_start, absolute_end,
                        fragment_id=f"{parent.fragment_id}_{ion_series}{len(ion_sequence)}",
                        fragment_type=f"MS2_{ion_series}", terminal_form="default")
                    modified = bool(state_fragment.included_modified_positions)
                    backbone = bool(state_fragment.included_backbone_modifications)
                    rows.append({
                        "Candidate_ID": structure.candidate_id,
                        "Complete_Structure_ID": state_fragment.complete_structure_id,
                        "Parent_Fragment_ID": parent.fragment_id, "Parent_Sequence": parent_sequence,
                        "Parent_Start": parent.start, "Parent_End": parent.end,
                        "Parent_Neutral_Mass": parent_state_fragment.neutral_exact_mass,
                        "Ion_ID": f"CMPION_{len(rows)+1:08d}", "Ion_Series": ion_series,
                        "Ion_Number": len(ion_sequence), "Cleavage_Position": int(parent.start)+cut-1,
                        "Ion_Sequence": ion_sequence, "Included_Positions": ";".join(map(str, state_fragment.included_positions)),
                        "Included_Modified_Positions": ";".join(map(str, state_fragment.included_modified_positions)),
                        "Included_Backbone_Bonds": ";".join(state_fragment.included_backbone_bonds),
                        "Theoretical_Neutral_Mass": state_fragment.neutral_exact_mass,
                        "Charge": charge, "Theoretical_mz": mz_from_neutral_mass(state_fragment.neutral_exact_mass, charge, polarity),
                        "Position_Informative": modified and len(ion_sequence) >= informative_length,
                        "Backbone_Informative": backbone and len(ion_sequence) >= informative_length,
                        "Audit_Level": audit_level, "Applied_To_Formal_Result": False, "Formal_Change_Ready": False,
                    })
                    if len(rows) >= max_rows:
                        return rows
    return rows