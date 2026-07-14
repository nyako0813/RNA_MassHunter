from dataclasses import replace
from pathlib import Path
import pytest
from rna_masshunter.modification_constraints import load_transformations, validate_transformation
from rna_masshunter.nucleoside_state import initial_nucleoside_state
ROOT=Path(__file__).parent; TS=load_transformations(ROOT/'data/modification_transforms_v2.yaml'); T={x.id:x for x in TS}
def state(base='U'): return initial_nucleoside_state(base,37,ROOT/'data/nucleoside_slots.yaml',TS)
@pytest.mark.parametrize('base,slots', [('A',('A_N1','A_N6')),('G',('G_O6','G_N7')),('C',('C_O2','C_N4')),('U',('U_O2','U_O4'))])
def test_initial_states_all_bases(base,slots):
    s=state(base); assert s.parent_base==base and all(x in s.slot_state_dict for x in slots)
def test_immutable_apply():
    before=state(); after,r=before.apply(T['s2U']); assert r.valid and before.slot_state_dict['U_O2']=='oxygen' and after.slot_state_dict['U_O2']=='sulfur'
def test_canonical_structure_id_and_equivalent_order():
    a,_=state('A').apply(T['m6A']); a,_=a.apply(T['Am']); b,_=state('A').apply(T['Am']); b,_=b.apply(T['m6A']); assert a.canonical_structure_id==b.canonical_structure_id
def test_wrong_parent_base(): assert validate_transformation(state('A'),T['s2U']).reason_code=='wrong_parent_base'
def test_from_state_mismatch():
    s,_=state().apply(T['s2U']); assert validate_transformation(s,T['s2U']).reason_code=='duplicate_component'
def test_missing_requirement_thioamide(): assert validate_transformation(state(),T['side_chain_thioamide']).reason_code=='from_state_mismatch'
def test_parent_derived_double_count():
    s,_=state().apply(T['s2U']); assert validate_transformation(s,T['ncm5s2U']).reason_code in {'parent_child_double_count','superseded_component'}
def test_same_slot_conflict():
    s,_=state().apply(T['ncm5U']); assert validate_transformation(s,T['cnm5U']).reason_code=='from_state_mismatch'
def test_oxidation_without_side_chain(): assert validate_transformation(state(),T['side_chain_thioamide_oxo1']).reason_code=='impossible_oxidation_state'
def test_position_compatibility_is_separate():
    tr=replace(T['s2U'],allowed_positions=(1,)); result=validate_transformation(state(),tr); assert result.valid and result.chemically_valid and not result.position_compatible and result.reason_code=='position_disallowed'
def test_pathway_compatibility_is_separate():
    tr=replace(T['s2U'],pathway_tags=('path-a',)); result=validate_transformation(state(),tr,pathway_context='path-b'); assert result.valid and not result.pathway_compatible and result.reason_code=='pathway_disallowed'
def test_schema_mass_is_calculated(): assert T['s2U'].exact_mass_delta==pytest.approx(T['s2U'].composition_delta.exact_mass)

def test_forbidden_state_reason():
    tr=replace(T['s2U'],forbids=(('U_O4',('oxygen',)),)); assert validate_transformation(state(),tr).reason_code=='forbidden_state'

def test_organism_compatibility_is_separate():
    tr=replace(T['s2U'],allowed_organisms=('org-a',)); r=validate_transformation(state(),tr,organism_context='org-b'); assert r.valid and not r.pathway_compatible
