from pathlib import Path
from types import SimpleNamespace

from rna_masshunter.base_loss_ions import (
    generate_base_loss_ions,
    generate_precursor_ions_for_base_loss,
    match_base_loss_ions,
)
from rna_masshunter.base_loss_masses import build_base_loss_masses
from rna_masshunter.masses import load_base_masses
from rna_masshunter.models import Fragment, MS2SpectrumInfo
from rna_masshunter.ms2_annotation import generate_theoretical_ms2_ions


def _base_masses():
    return load_base_masses(Path(__file__).parent / "data" / "base_masses.yaml")


def _config(**overrides):
    ms2_annotation = {
        "use_theoretical_fragments": True,
        "min_ion_length": 1,
        "include_base_loss": True,
        "base_loss_min_ion_length": 2,
        "base_loss_max_rows": 50000,
        "base_loss_include_precursor": True,
        "mz_tolerance_ppm": 20,
    }
    ms2_annotation.update(overrides)
    return SimpleNamespace(instrument={"polarity": "negative"}, ms2_annotation=ms2_annotation)


def _fragment(sequence="ACGU", terminal_form="default"):
    return Fragment("F1", "target", sequence, 1, len(sequence), 1, len(sequence), "RNase_T1", 0, terminal_form, 1000.0)


def _table():
    return build_base_loss_masses([], _base_masses())


def test_base_loss_variants_generated_for_each_backbone_and_precursor_ion():
    base_masses = _base_masses()
    config = _config()
    table = _table()
    ions = generate_theoretical_ms2_ions([_fragment("ACGU")], config, base_masses)
    precursor_ions = generate_precursor_ions_for_base_loss([_fragment("ACGU")], config, base_masses)
    assert {ion.ion_type for ion in precursor_ions} == {"precursor"}

    base_loss_ions = generate_base_loss_ions(ions + precursor_ions, config, table)
    assert base_loss_ions

    # Spot-check one d ion (prefix "AC") loses each of its two bases correctly.
    d_ac = next(ion for ion in ions if ion.ion_type == "d" and ion.ion_sequence == "AC")
    variants = [row for row in base_loss_ions if row.ion_type == "d" and row.ion_sequence == "AC"]
    assert {row.base_loss_base for row in variants} == {"A", "C"}
    for row in variants:
        expected = d_ac.theoretical_mass - table[row.base_loss_base]
        assert abs(row.theoretical_mass - expected) < 1e-9
        assert row.base_loss_position == (1 if row.base_loss_base == "A" else 2)

    # Precursor ion (full length "ACGU") also gets base-loss variants.
    precursor_variants = [row for row in base_loss_ions if row.ion_type == "precursor"]
    assert {row.base_loss_position for row in precursor_variants} == {1, 2, 3, 4}


def test_min_ion_length_excludes_short_ions():
    base_masses = _base_masses()
    config = _config()
    table = _table()
    ions = generate_theoretical_ms2_ions([_fragment("ACGU")], config, base_masses)
    one_nt_ions = [ion for ion in ions if len(ion.ion_sequence) == 1]
    assert one_nt_ions
    base_loss_ions = generate_base_loss_ions(ions, config, table)
    assert all(len(row.ion_sequence) >= 2 for row in base_loss_ions)

    config_len3 = _config(base_loss_min_ion_length=3)
    base_loss_ions_len3 = generate_base_loss_ions(ions, config_len3, table)
    assert all(len(row.ion_sequence) >= 3 for row in base_loss_ions_len3)


def test_max_rows_truncates_generation():
    base_masses = _base_masses()
    table = _table()
    ions = generate_theoretical_ms2_ions([_fragment("ACGUACGUACGU")], _config(), base_masses)
    config = _config(base_loss_max_rows=5)
    base_loss_ions = generate_base_loss_ions(ions, config, table)
    assert len(base_loss_ions) == 5


def test_include_base_loss_false_returns_empty_list():
    base_masses = _base_masses()
    table = _table()
    ions = generate_theoretical_ms2_ions([_fragment("ACGU")], _config(), base_masses)
    config = _config(include_base_loss=False)
    assert generate_base_loss_ions(ions, config, table) == []
    assert generate_precursor_ions_for_base_loss([_fragment("ACGU")], config, base_masses) == []


def test_base_loss_include_precursor_false_omits_precursor_variants():
    base_masses = _base_masses()
    config = _config(base_loss_include_precursor=False)
    precursor_ions = generate_precursor_ions_for_base_loss([_fragment("ACGU")], config, base_masses)
    assert precursor_ions == []

    ions = generate_theoretical_ms2_ions([_fragment("ACGU")], config, base_masses)
    table = _table()
    base_loss_ions = generate_base_loss_ions(ions + precursor_ions, config, table)
    assert all(row.ion_type != "precursor" for row in base_loss_ions)


def test_precursor_mass_uses_fragment_terminal_form():
    base_masses = _base_masses()
    config = _config()
    default_ions = generate_precursor_ions_for_base_loss(
        [_fragment("ACGU", terminal_form="default")], config, base_masses,
    )
    residual_phosphate_ions = generate_precursor_ions_for_base_loss(
        [_fragment("ACGU", terminal_form="residual_phosphate")], config, base_masses,
    )
    default_mass = next(ion.theoretical_mass for ion in default_ions if ion.charge == 1)
    residual_mass = next(ion.theoretical_mass for ion in residual_phosphate_ions if ion.charge == 1)
    phosphate = base_masses["constants"]["phosphate"]
    assert abs((residual_mass - default_mass) - phosphate) < 1e-6


def test_match_base_loss_ions_finds_peak_within_tolerance():
    base_masses = _base_masses()
    config = _config()
    table = _table()
    ions = generate_theoretical_ms2_ions([_fragment("ACGU")], config, base_masses)
    base_loss_ions = generate_base_loss_ions(ions, config, table)
    target = base_loss_ions[0]
    spectrum = MS2SpectrumInfo(
        "S1", 1, 1.5, 1000.0, 1, None, 1, target.theoretical_mz, 1000.0, 1000.0,
        [(target.theoretical_mz, 1000.0)],
    )
    matches = match_base_loss_ions([spectrum], base_loss_ions, config)
    assert matches
    assert matches[0].best_ion_id == target.ion_id
    assert matches[0].base_loss_base == target.base_loss_base
    assert matches[0].base_loss_position == target.base_loss_position


if __name__ == "__main__":
    test_base_loss_variants_generated_for_each_backbone_and_precursor_ion()
    test_min_ion_length_excludes_short_ions()
    test_max_rows_truncates_generation()
    test_include_base_loss_false_returns_empty_list()
    test_base_loss_include_precursor_false_omits_precursor_variants()
    test_precursor_mass_uses_fragment_terminal_form()
    test_match_base_loss_ions_finds_peak_within_tolerance()
    print("base loss ion tests: OK")
