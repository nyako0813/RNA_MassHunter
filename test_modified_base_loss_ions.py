from pathlib import Path
from types import SimpleNamespace

from rna_masshunter.base_loss_masses import build_base_loss_masses
from rna_masshunter.masses import load_base_masses
from rna_masshunter.models import Modification, MS2SpectrumInfo
from rna_masshunter.modified_fragment_ions import (
    generate_modified_base_loss_ions,
    generate_modified_precursor_ions_for_base_loss,
    generate_modified_theoretical_ions,
    match_modified_base_loss_ions,
)


def _base_masses():
    return load_base_masses(Path(__file__).parent / "data" / "base_masses.yaml")


def _config(**overrides):
    ms2_annotation = {
        "include_modified_fragment_ions": True,
        "modified_fragment_require_target_base": True,
        "modified_fragment_include_unmodified_counterparts": True,
        "modified_fragment_max_positions_per_candidate": 20,
        "modified_fragment_min_ion_length": 1,
        "modified_fragment_min_ion_length_for_localization": 2,
        "modified_fragment_max_rows": 100000,
        "mz_tolerance_ppm": 20,
        "include_base_loss": True,
        "base_loss_min_ion_length": 2,
        "base_loss_max_rows": 50000,
        "base_loss_include_precursor": True,
    }
    ms2_annotation.update(overrides)
    return SimpleNamespace(instrument={"polarity": "negative"}, ms2_annotation=ms2_annotation)


def _parent(modification_id="m1A", mass_shift=14.0156, modified_mass=None):
    # "ACGU" with an "A" modification: position 1 is the candidate site.
    return {
        "Spectrum_ID": "S1", "Candidate_Parent_Fragment_ID": "F1",
        "Candidate_Parent_Sequence": "ACGU", "Candidate_Parent_Start": 10,
        "Candidate_Parent_End": 13, "Candidate_Type": "modified",
        "Modification_ID": modification_id, "Modification_Name": "1-methyladenosine",
        "Modification_Target_Base": "A", "Modification_Mass_Shift": mass_shift,
        "Parent_Charge": 1,
        "Parent_Modified_Mass": modified_mass if modified_mass is not None else (1000.0 + mass_shift),
    }


def _table():
    base_masses = _base_masses()
    mods = [Modification(id="m1A", symbol="m1A", mass_shift_from_unmodified=14.0156, category="biological", target_bases=["A"])]
    return build_base_loss_masses(mods, base_masses)


def test_modified_position_uses_modification_base_loss_mass():
    config = _config()
    table = _table()
    parent = _parent()
    ions = generate_modified_theoretical_ions([parent], config, _base_masses())
    base_loss_ions = generate_modified_base_loss_ions(ions, config, table)
    assert base_loss_ions

    # An ion containing the modified position (position 1) should apply the
    # modification's own base-loss mass at that offset.
    modified_hits = [row for row in base_loss_ions if row["Base_Loss_Is_Modified_Base"]]
    assert modified_hits
    for row in modified_hits:
        assert row["Base_Loss_Base_Or_Modification_ID"] == "m1A"
        assert row["Base_Loss_Position_In_Parent"] == row["Candidate_Modification_Position_In_Parent"]


def test_unmodified_positions_use_canonical_base_loss_mass():
    config = _config()
    table = _table()
    parent = _parent()
    ions = generate_modified_theoretical_ions([parent], config, _base_masses())
    base_loss_ions = generate_modified_base_loss_ions(ions, config, table)
    unmodified_hits = [row for row in base_loss_ions if not row["Base_Loss_Is_Modified_Base"]]
    assert unmodified_hits
    for row in unmodified_hits:
        assert row["Base_Loss_Base_Or_Modification_ID"] in {"A", "U", "G", "C"}


def test_unresolvable_modification_id_is_skipped_without_error():
    config = _config()
    table = _table()  # only knows "m1A", "A", "U", "G", "C"
    parent = _parent(modification_id="mystery_mod_xyz")
    ions = generate_modified_theoretical_ions([parent], config, _base_masses())
    base_loss_ions = generate_modified_base_loss_ions(ions, config, table)
    # Every row touching the modified position (key = "mystery_mod_xyz") must
    # be skipped since it isn't in the lookup table; only unmodified-position
    # rows (canonical base keys) can survive.
    assert all(not row["Base_Loss_Is_Modified_Base"] for row in base_loss_ions)


def test_modified_precursor_ions_for_base_loss_generate_full_length_rows():
    config = _config()
    parent = _parent()
    precursor_rows = generate_modified_precursor_ions_for_base_loss([parent], config)
    assert precursor_rows
    for row in precursor_rows:
        assert row["Ion_Type"] == "precursor"
        assert row["Ion_Sequence"] == "ACGU"
        assert row["Ion_Start"] == 1 and row["Ion_End"] == 4
        assert row["Ion_Contains_Modification"] is True
        assert row["Theoretical_Mass"] == parent["Parent_Modified_Mass"]

    table = _table()
    base_loss_ions = generate_modified_base_loss_ions(precursor_rows, config, table)
    positions = {row["Base_Loss_Position_In_Parent"] for row in base_loss_ions}
    assert positions == {1, 2, 3, 4}
    modified_position_row = next(row for row in base_loss_ions if row["Base_Loss_Is_Modified_Base"])
    assert modified_position_row["Base_Loss_Position_In_Parent"] == modified_position_row["Candidate_Modification_Position_In_Parent"]


def test_base_loss_include_precursor_false_omits_precursor_rows():
    config = _config(base_loss_include_precursor=False)
    parent = _parent()
    assert generate_modified_precursor_ions_for_base_loss([parent], config) == []


def test_include_base_loss_false_returns_empty_list():
    config = _config(include_base_loss=False)
    table = _table()
    parent = _parent()
    ions = generate_modified_theoretical_ions([parent], config, _base_masses())
    assert generate_modified_base_loss_ions(ions, config, table) == []
    assert generate_modified_precursor_ions_for_base_loss([parent], config) == []


def test_match_modified_base_loss_ions_finds_peak_within_tolerance():
    config = _config()
    table = _table()
    parent = _parent()
    ions = generate_modified_theoretical_ions([parent], config, _base_masses())
    base_loss_ions = generate_modified_base_loss_ions(ions, config, table)
    target = next(row for row in base_loss_ions if row["Base_Loss_Is_Modified_Base"])
    spectrum = MS2SpectrumInfo(
        "S1", 1, 1.5, 1000.0, 1, None, 1, target["Theoretical_mz"], 1000.0, 1000.0,
        [(target["Theoretical_mz"], 1000.0)],
    )
    matches = match_modified_base_loss_ions([spectrum], base_loss_ions, config)
    assert matches
    assert any(row["Match_Status"] == "matched_modified_base_loss_ion" for row in matches)


if __name__ == "__main__":
    test_modified_position_uses_modification_base_loss_mass()
    test_unmodified_positions_use_canonical_base_loss_mass()
    test_unresolvable_modification_id_is_skipped_without_error()
    test_modified_precursor_ions_for_base_loss_generate_full_length_rows()
    test_base_loss_include_precursor_false_omits_precursor_rows()
    test_include_base_loss_false_returns_empty_list()
    test_match_modified_base_loss_ions_finds_peak_within_tolerance()
    print("modified base loss ion tests: OK")
