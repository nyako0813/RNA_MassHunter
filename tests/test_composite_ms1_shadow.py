from pathlib import Path
from rna_masshunter.backbone_state import load_backbone_transformations
from rna_masshunter.modification_constraints import load_transformations
from rna_masshunter.sample_structure_schema import load_sample_structure_hypotheses
from rna_masshunter.structure_fragment import build_complete_structure_state, extract_fragment_from_structure

ROOT = Path(__file__).resolve().parents[1]
SEQ = "GCUCCGGUAGUGUAGUCCGGCCAAUCAUUCCGGCCUUUCGAGCCGAAGACUCGGGUUCGAAUCCCGGCCGGAGCACCA"
TARGET_IDENTITY = {"name": "Mac_tRNA-Glu-UUC", "organism": "Methanosarcina acetivorans", "rule_set": "methanosarcina_acetivorans"}

def _structure():
    transforms = load_transformations(ROOT / "data/modification_transforms_v2.yaml")
    loaded = load_sample_structure_hypotheses(
        ROOT / "data/sample_structure_hypotheses.yaml", sequence=SEQ, transformations=transforms,
        backbone_bond_ids={f"{i}_{i+1}" for i in range(1, len(SEQ))}, target_identity=TARGET_IDENTITY,
    )
    state, error = build_complete_structure_state(
        next(x for x in loaded.hypotheses if x.hypothesis_id == "U37_side_chain_plus_PT_37_38"), SEQ, transforms, ROOT / "data/nucleoside_slots.yaml",
        load_backbone_transformations(ROOT / "data/backbone_modifications.yaml")[0],
    )
    assert error is None
    return loaded, state

def test_sample_schema_loads_valid_complete_hypothesis():
    loaded, state = _structure()
    assert loaded.enabled
    assert len(loaded.hypotheses) == 3
    assert not loaded.invalid_rows
    assert state.position_states[37].applied_transform_ids == ("s2U", "cnm5U", "side_chain_thioamide_oxo1")
    assert state.bond_states["37_38"].state == "phosphorothioate"

def test_fragment_complete_composition_applies_only_included_states_once():
    _, state = _structure()
    inside = extract_fragment_from_structure(state, SEQ, 35, 40)
    outside = extract_fragment_from_structure(state, SEQ, 1, 4)
    expected_delta = state.position_states[37].exact_mass_delta + state.bond_states["37_38"].exact_mass_delta
    from rna_masshunter.structure_fragment import unmodified_fragment_composition
    assert abs(inside.neutral_exact_mass - unmodified_fragment_composition(SEQ[34:40]).exact_mass - expected_delta) < 1e-9
    assert outside.neutral_exact_mass == unmodified_fragment_composition(SEQ[:4]).exact_mass
    assert inside.included_modified_positions == (37,)
    assert inside.included_backbone_modifications == ("37_38",)

def test_sample_schema_rejects_wrong_target_identity():
    transforms = load_transformations(ROOT / "data/modification_transforms_v2.yaml")
    loaded = load_sample_structure_hypotheses(
        ROOT / "data/sample_structure_hypotheses.yaml", sequence=SEQ, transformations=transforms,
        backbone_bond_ids={f"{i}_{i+1}" for i in range(1, len(SEQ))},
        target_identity={**TARGET_IDENTITY, "name": "different-target"},
    )
    assert not loaded.hypotheses
    assert loaded.invalid_rows[0]["Invalid_Reason"] == "target_identity_mismatch"
