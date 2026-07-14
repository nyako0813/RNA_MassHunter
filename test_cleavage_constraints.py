from pathlib import Path
import pytest
from rna_masshunter.backbone_state import load_backbone_transformations, normal_bond
from rna_masshunter.cleavage_constraints import evaluate_cleavage
from rna_masshunter.masses import load_base_masses
ROOT=Path(__file__).parent; PS=load_backbone_transformations(ROOT/'data/backbone_modifications.yaml')[0]; MASSES=load_base_masses(ROOT/'data/base_masses.yaml')
def bond(left):
    b=normal_bond(left,left+1).apply(PS); return {b.bond_id:b}
def test_t1_normal_g_cleavage(): assert 2 in evaluate_cleavage('AGU','RNase_T1',{},PS).allowed_sites
def test_t1_phosphorothioate_blocked():
    r=evaluate_cleavage('AGU','RNase_T1',bond(2),PS,MASSES); assert r.blocked_sites==(2,) and r.fragments[0].sequence=='AGU'
def test_p1_normal_cleavage(): assert evaluate_cleavage('ACG','Nuclease_P1',{},PS).allowed_sites==(1,2)
def test_p1_dinucleotide_residual_and_mass_once():
    r=evaluate_cleavage('ACG','Nuclease_P1',bond(1),PS,MASSES); f=r.fragments[0]; assert f.sequence=='AC' and f.backbone_modification_count==1 and f.backbone_mass_delta==pytest.approx(PS.exact_mass_delta)
def test_rnase_a_rule_potential_does_not_force_block():
    r=evaluate_cleavage('CU','RNase_A',bond(1),PS); assert r.evaluations[0].status=='potentially_blocked' and r.allowed_sites==(1,)
def test_rnase_t2_rule_potential(): assert evaluate_cleavage('AU','RNase_T2',bond(1),PS).evaluations[0].status=='potentially_blocked'
def test_stochastic_missed_origin():
    r=evaluate_cleavage('ACG','Nuclease_P1',{},PS,stochastic_missed_sites=(1,)); assert r.fragments[0].cleavage_origin=='stochastic_missed' and r.fragments[0].stochastic_missed_cleavage_count==1
def test_blocked_origin_distinct(): assert evaluate_cleavage('ACG','Nuclease_P1',bond(1),PS).fragments[0].cleavage_origin=='phosphorothioate_blocked'
def test_mixed_origin():
    bonds=bond(1); r=evaluate_cleavage('ACG','Nuclease_P1',bonds,PS,stochastic_missed_sites=(2,)); assert r.fragments[0].cleavage_origin=='mixed_stochastic_and_blocked'
def test_adjacent_blocked_bonds_make_oligomer_and_mass_twice():
    bonds={**bond(1),**bond(2)}; r=evaluate_cleavage('ACGU','Nuclease_P1',bonds,PS,MASSES); f=r.fragments[0]; assert f.sequence=='ACG' and f.backbone_modification_count==2 and f.backbone_mass_delta==pytest.approx(2*PS.exact_mass_delta)
def test_fragment_boundaries():
    r=evaluate_cleavage('ACGU','Nuclease_P1',bond(2),PS); assert [(x.start,x.end) for x in r.fragments]==[(1,1),(2,3),(4,4)]
def test_terminal_form_preserved(): assert evaluate_cleavage('AC','Nuclease_P1',bond(1),PS,MASSES,terminal_form='default').fragments[0].terminal_form=='default'

def test_unknown_blocking_rule():
    unknown=type(PS)(PS.id,PS.name,PS.from_state,PS.to_state,PS.composition_delta,PS.stereochemistry,PS.evidence_status,(),PS.notes)
    assert evaluate_cleavage('AC','Nuclease_P1',bond(1),unknown).evaluations[0].status=='unknown'
