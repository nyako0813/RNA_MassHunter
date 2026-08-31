from pathlib import Path
from types import SimpleNamespace

from rna_masshunter.masses import calculate_unmodified_rna_mass, fragment_ion_series_offsets, load_base_masses
from rna_masshunter.models import Fragment
from rna_masshunter.ms2_annotation import generate_theoretical_ms2_ions


def _base_masses():
    return load_base_masses(Path(__file__).parent / "data" / "base_masses.yaml")


def _config(ion_series=None):
    ms2_annotation = {"use_theoretical_fragments": True, "min_ion_length": 1}
    if ion_series is not None:
        ms2_annotation["ion_series"] = ion_series
    return SimpleNamespace(instrument={"polarity": "negative"}, ms2_annotation=ms2_annotation)


def _fragment(sequence="ACGU"):
    return Fragment("F1", "target", sequence, 1, len(sequence), 1, len(sequence), "RNase_T1", 0, "default", 1000.0)


def test_default_ion_series_generates_d_w_a_z_with_correct_mass_offsets():
    base_masses = _base_masses()
    ions = generate_theoretical_ms2_ions([_fragment("ACGU")], _config(), base_masses)
    by_type_and_seq = {(ion.ion_type, ion.ion_sequence): ion for ion in ions}
    assert {ion.ion_type for ion in ions} == {"d", "w", "a", "z"}

    offsets = fragment_ion_series_offsets(base_masses)
    for cut in range(1, 4):
        prefix, suffix = "ACGU"[:cut], "ACGU"[cut:]
        prefix_mass = calculate_unmodified_rna_mass(prefix, base_masses, terminal_form="default")
        suffix_mass = calculate_unmodified_rna_mass(suffix, base_masses, terminal_form="default")

        d_ion = by_type_and_seq[("d", prefix)]
        a_ion = by_type_and_seq[("a", prefix)]
        w_ion = by_type_and_seq[("w", suffix)]
        z_ion = by_type_and_seq[("z", suffix)]

        assert d_ion.theoretical_mass == prefix_mass
        assert w_ion.theoretical_mass == suffix_mass
        assert abs(a_ion.theoretical_mass - (prefix_mass + offsets["a"])) < 1e-9
        assert abs(z_ion.theoretical_mass - (suffix_mass + offsets["z"])) < 1e-9

        # a/z are d/z minus a complete phosphate + water relative to d/w.
        assert abs((d_ion.theoretical_mass - a_ion.theoretical_mass) - (offsets["d"] - offsets["a"])) < 1e-9
        assert abs((w_ion.theoretical_mass - z_ion.theoretical_mass) - (offsets["w"] - offsets["z"])) < 1e-9


def test_complementary_pair_mass_conservation():
    # For a given cut, a(prefix) + w(suffix) should equal d(prefix) + z(suffix):
    # both offsets equal -(phosphate+water) applied to exactly one side.
    base_masses = _base_masses()
    ions = generate_theoretical_ms2_ions([_fragment("ACGU")], _config(), base_masses)
    by_key = {(ion.ion_type, ion.ion_start, ion.ion_end): ion.theoretical_mass for ion in ions}
    cut = 2
    a_mass = by_key[("a", 1, cut)]
    w_mass = by_key[("w", cut + 1, 4)]
    d_mass = by_key[("d", 1, cut)]
    z_mass = by_key[("z", cut + 1, 4)]
    assert abs((a_mass + w_mass) - (d_mass + z_mass)) < 1e-9


def test_ion_series_config_filters_generated_types():
    base_masses = _base_masses()
    ions = generate_theoretical_ms2_ions([_fragment("ACGU")], _config(ion_series=["d", "w"]), base_masses)
    assert ions
    assert {ion.ion_type for ion in ions} == {"d", "w"}


def test_unknown_ion_series_values_fall_back_to_default_four():
    base_masses = _base_masses()
    ions = generate_theoretical_ms2_ions([_fragment("ACGU")], _config(ion_series=["nonsense"]), base_masses)
    assert {ion.ion_type for ion in ions} == {"d", "w", "a", "z"}


if __name__ == "__main__":
    test_default_ion_series_generates_d_w_a_z_with_correct_mass_offsets()
    test_complementary_pair_mass_conservation()
    test_ion_series_config_filters_generated_types()
    test_unknown_ion_series_values_fall_back_to_default_four()
    print("ms2 d/w/a/z ion series tests: OK")
