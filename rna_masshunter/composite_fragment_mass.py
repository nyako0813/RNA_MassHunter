"""Exact mass helpers for complete structure fragments."""
from __future__ import annotations
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.structure_fragment import StructureFragment

def fragment_mass_row(fragment: StructureFragment, audit_level: str = "full") -> dict:
    return {
        "Candidate_ID": fragment.candidate_id,
        "Complete_Structure_ID": fragment.complete_structure_id,
        "Fragment_ID": fragment.fragment_id,
        "Fragment_Type": fragment.fragment_type,
        "Start_Position": fragment.start, "End_Position": fragment.end,
        "Included_Positions": ";".join(map(str, fragment.included_positions)),
        "Included_Modified_Positions": ";".join(map(str, fragment.included_modified_positions)),
        "Included_Backbone_Bonds": ";".join(fragment.included_backbone_bonds),
        "Included_Backbone_Modifications": ";".join(fragment.included_backbone_modifications),
        "Terminal_Form": fragment.terminal_form,
        "Fragment_Elemental_Composition": fragment.elemental_composition_canonical,
        "Neutral_Exact_Mass": fragment.neutral_exact_mass,
        "Audit_Level": audit_level, "Applied_To_Formal_Result": False, "Formal_Change_Ready": False,
    }

def theoretical_mz_from_mass(neutral_mass: float, charge: int, polarity: str = "negative") -> float:
    return mz_from_neutral_mass(neutral_mass, charge, polarity)