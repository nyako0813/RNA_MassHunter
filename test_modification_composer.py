from pathlib import Path
import pytest
from rna_masshunter.modification_composer import compose_modifications
from rna_masshunter.modification_constraints import load_transformations
ROOT=Path(__file__).parent; TRANS=load_transformations(ROOT/'data/modification_transforms_v2.yaml')
def compose(base='U',max_components=3): return compose_modifications(base,37,TRANS,ROOT/'data/nucleoside_slots.yaml',max_components=max_components)
def ids(result): return {tuple(t.id for t in c.transforms) for c in result.valid_candidates}
@pytest.mark.parametrize('wanted',[('s2U',),('ncm5U',),('cnm5U',),('ncm5s2U',),('cnm5s2U',),('Um','ncm5s2U')])
def test_valid_u_candidates(wanted): assert tuple(sorted(wanted)) in ids(compose())
@pytest.mark.parametrize('base,wanted',[('C',('Cm','m5C')),('G',('Gm','m7G')),('A',('Am','m6A'))])
def test_valid_other_base_composites(base,wanted): assert tuple(sorted(wanted)) in ids(compose(base))
def test_side_chain_thioamide_and_c2_thio():
    assert any(c.state.slot_state_dict['U_O2']=='sulfur' and c.state.slot_state_dict['U_C5_side_chain_carbonyl']=='thioamide_sulfur' for c in compose().valid_candidates)
def test_hypothetical_oxidation_example():
    assert any(c.state.slot_state_dict['U_O2']=='sulfur' and c.state.slot_state_dict['U_C5_side_chain_carbonyl']=='oxidized_sulfur_1' for c in compose().valid_candidates)
@pytest.mark.parametrize('pair',[('cnm5U','ncm5U'),('m5U','cm5U'),('s2U','ncm5s2U'),('ncm5U','ncm5s2U'),('side_chain_thioamide','side_chain_thioamide_oxo1')])
def test_invalid_pairs(pair):
    attempts={x.transform_ids:x.result.reason_code for x in compose().invalid_attempts}; assert tuple(sorted(pair)) in attempts
def test_max_components_enforced(): assert max(c.component_count for c in compose(max_components=2).valid_candidates)<=2
def test_deterministic_order():
    a=compose(); b=compose(); assert [c.state.canonical_structure_id for c in a.valid_candidates]==[c.state.canonical_structure_id for c in b.valid_candidates]
def test_equivalent_candidate_dedup():
    r=compose(); keys=[(c.state.slot_states,c.state.elemental_composition_delta) for c in r.valid_candidates]; assert len(keys)==len(set(keys))
def test_isomer_grouping(): assert compose('G').isomer_groups
def test_valid_invalid_separation():
    r=compose(); assert r.valid_candidates and r.invalid_attempts and all(x.result.valid is False for x in r.invalid_attempts)
@pytest.mark.parametrize('base','AGCU')
def test_all_bases_supported(base): assert compose(base).valid_candidates
