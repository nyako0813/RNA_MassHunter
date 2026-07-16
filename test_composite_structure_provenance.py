from copy import deepcopy
from dataclasses import replace
from pathlib import Path

from rna_masshunter.audit_policy import AuditPolicy, included_sheet_names
from rna_masshunter.backbone_state import load_backbone_transformations
from rna_masshunter.composite_structure_provenance import build_composite_structure_provenance
from rna_masshunter.models import Modification
from rna_masshunter.modification_constraints import load_transformations
from rna_masshunter.sample_structure_schema import load_sample_structure_hypotheses
from rna_masshunter.structure_fragment import CompleteStructureState, build_complete_structure_state, complete_structure_id

ROOT = Path(__file__).resolve().parent
SEQ = "GCUCCGGUAGUGUAGUCCGGCCAAUCAUUCCGGCCUUUCGAGCCGAAGACUCGGGUUCGAAUCCCGGCCGGAGCACCA"
TARGET = {"name": "Mac_tRNA-Glu-UUC", "organism": "Methanosarcina acetivorans", "rule_set": "methanosarcina_acetivorans"}


def structures():
    transforms = load_transformations(ROOT / "data/modification_transforms_v2.yaml")
    loaded = load_sample_structure_hypotheses(
        ROOT / "data/sample_structure_hypotheses.yaml", sequence=SEQ,
        transformations=transforms,
        backbone_bond_ids={f"{i}_{i+1}" for i in range(1, len(SEQ))},
        target_identity=TARGET,
    )
    backbone = load_backbone_transformations(ROOT / "data/backbone_modifications.yaml")[0]
    output = []
    for hypothesis in loaded.hypotheses:
        state, error = build_complete_structure_state(
            hypothesis, SEQ, transforms, ROOT / "data/nucleoside_slots.yaml", backbone,
        )
        assert error is None
        output.append(state)
    return output


def test_position_transform_explicit_mass_only_and_canonical_are_separate():
    state = structures()[0]
    position, nucleoside = next(iter(state.position_states.items()))
    explicit = {
        legacy for transform_id in nucleoside.applied_transform_ids
        for legacy in nucleoside.transform_map[transform_id].legacy_ids
    }
    explicit_id = sorted(explicit)[0] if explicit else "EXPLICIT_TEST"
    modifications = [
        Modification(explicit_id, explicit_id, nucleoside.exact_mass_delta, "x", [nucleoside.parent_base]),
        Modification("MASS_ONLY", "MASS_ONLY", nucleoside.exact_mass_delta, "x", [nucleoside.parent_base]),
    ]
    row = build_composite_structure_provenance([state], modifications).position_rows[0]
    assert row["Applied_Transform_IDs"] == ";".join(sorted(nucleoside.applied_transform_ids))
    if explicit:
        assert explicit_id in row["Explicit_Legacy_Modification_IDs"].split(";")
        assert explicit_id not in row["Mass_Equivalent_Modification_IDs"].split(";")
    assert "MASS_ONLY" in row["Mass_Equivalent_Modification_IDs"].split(";")
    assert row["Canonical_Structure_ID"] == nucleoside.canonical_structure_id
    assert row["Complete_Structure_ID"] == complete_structure_id(state)
    assert row["Composite_Position"] == position


def test_isomer_information_is_derived_from_states():
    original = structures()[0]
    position, state = next(iter(original.position_states.items()))
    changed_slots = tuple(sorted((*state.slot_states, ("synthetic_isomer_marker", "alternate"))))
    alternate_state = replace(state, slot_states=changed_slots)
    alternate = CompleteStructureState(
        "ALT", {position: alternate_state}, original.bond_states,
        original.terminal_five_prime, original.terminal_three_prime,
    )
    first = CompleteStructureState(
        "ORIGINAL", {position: state}, original.bond_states,
        original.terminal_five_prime, original.terminal_three_prime,
    )
    rows = build_composite_structure_provenance([first, alternate]).position_rows
    assert len(rows) == 2
    assert all(row["Is_Isomeric"] is True for row in rows)
    assert len({row["Isomer_Group_ID"] for row in rows}) == 1
    assert len({row["Canonical_Structure_ID"] for row in rows}) == 2


def test_bond_transform_is_retained_and_not_mixed_with_nucleoside_rows():
    state = next(item for item in structures() if item.bond_states)
    result = build_composite_structure_provenance([state])
    bond = next(iter(state.bond_states.values()))
    row = result.bond_rows[0]
    assert row["Applied_Backbone_Transform_IDs"] == ";".join(sorted(bond.applied_transform_ids))
    assert row["Bond_ID"] == bond.bond_id
    assert "Bond_ID" not in result.position_rows[0]
    assert "Composite_Position" not in row


def test_empty_input_and_input_order_determinism():
    empty = build_composite_structure_provenance([])
    assert empty.position_rows == [] and empty.bond_rows == []
    states = structures()
    assert build_composite_structure_provenance(states) == build_composite_structure_provenance(list(reversed(states)))


def test_all_formal_flags_are_false_and_inputs_unchanged():
    states = structures()
    before = deepcopy(states)
    result = build_composite_structure_provenance(states)
    for rows in (result.position_rows, result.bond_rows):
        for row in rows:
            assert row["Applied_To_Formal_Result"] is False
            assert row["Formal_Change_Ready"] is False
            assert row["Formal_Result_Changed"] is False
    assert states == before


def test_sheet_inclusion():
    names = ["Composite_Structure_Position_Ma", "Composite_Structure_Bond_Map"]
    standard, _ = included_sheet_names(names, AuditPolicy.from_level("standard"))
    audit, _ = included_sheet_names(names, AuditPolicy.from_level("audit"))
    full, _ = included_sheet_names(names, AuditPolicy.from_level("full"))
    assert standard == [] and audit == [] and full == names
