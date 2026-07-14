from pathlib import Path
import pytest
from rna_masshunter.backbone_state import load_backbone_transformations, normal_bond
ROOT=Path(__file__).parent; PS=load_backbone_transformations(ROOT/'data/backbone_modifications.yaml')[0]
def test_normal_bond(): assert normal_bond(36,37).state=='normal_phosphate'
def test_bond_id(): assert normal_bond(36,37).bond_id=='36_37'
def test_phosphorothioate_bond(): assert normal_bond(36,37).apply(PS).state=='phosphorothioate'
def test_o_to_s_mass_delta(): assert normal_bond(1,2).apply(PS).exact_mass_delta==pytest.approx(15.97715658043)
def test_stereochemistry_unknown(): assert normal_bond(1,2).apply(PS).stereochemistry=='unknown'
def test_duplicate_phosphorothioate_prevented():
    with pytest.raises(ValueError,match='duplicate_component'): normal_bond(1,2).apply(PS).apply(PS)
def test_nonadjacent_bond_rejected():
    with pytest.raises(ValueError): normal_bond(1,3)
