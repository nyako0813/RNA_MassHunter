from pathlib import Path

from rna_masshunter.base_loss_masses import build_base_loss_masses
from rna_masshunter.masses import load_base_masses
from rna_masshunter.modifications import load_modifications
from rna_masshunter.models import Modification

ROOT = Path(__file__).parent
BASE_MASSES = load_base_masses(ROOT / "data" / "base_masses.yaml")


def _mod(mod_id, target_bases, shift):
    return Modification(
        id=mod_id, symbol=mod_id, mass_shift_from_unmodified=shift,
        category="biological", target_bases=target_bases,
    )


def test_canonical_bases_use_exact_elemental_masses():
    table = build_base_loss_masses([], BASE_MASSES)
    assert table["A"] == 135.0545
    assert table["G"] == 151.0494
    assert table["C"] == 111.0433
    assert table["U"] == 112.0273


def test_pure_base_modification_adds_full_shift_to_canonical_base():
    mods = [_mod("m1A", ["A"], 14.0156)]
    table = build_base_loss_masses(mods, BASE_MASSES)
    assert table["m1A"] == 135.0545 + 14.0156


def test_pure_sugar_modification_keeps_canonical_base_mass():
    # Symbol ends in "m" and the stem is a bare canonical base letter:
    # the base itself is unmodified, only the sugar carries the methyl.
    mods = [_mod("Am", ["A"], 14.0156)]
    table = build_base_loss_masses(mods, BASE_MASSES)
    assert table["Am"] == 135.0545


def test_compound_base_and_sugar_modification_uses_base_component_shift_only():
    # m1Am's own curated shift (28.0313) mixes base (m1A: 14.0156) and
    # sugar (Am: 14.0156) contributions. Base-loss must use only the
    # base-side (m1A) contribution, matched via the naming convention,
    # not the compound's own mixed shift value.
    mods = [
        _mod("m1A", ["A"], 14.0156),
        _mod("Am", ["A"], 14.0156),
        _mod("m1Am", ["A"], 28.0313),
    ]
    table = build_base_loss_masses(mods, BASE_MASSES)
    assert table["m1Am"] == table["m1A"] == 135.0545 + 14.0156


def test_unresolvable_compound_modification_is_excluded_and_warns():
    # "m9Xm" implies a base component "m9X" that does not exist in the
    # catalog; it must be excluded rather than guessed.
    mods = [_mod("m9Xm", ["A"], 99.0)]
    warnings = []
    table = build_base_loss_masses(mods, BASE_MASSES, warnings)
    assert "m9Xm" not in table
    assert any(w["Message"].startswith("Compound base+sugar modification could not be decomposed") for w in warnings)


def test_modification_with_unresolvable_target_base_is_excluded_and_warns():
    mods = [_mod("weird_mod", [], 5.0)]
    warnings = []
    table = build_base_loss_masses(mods, BASE_MASSES, warnings)
    assert "weird_mod" not in table
    assert any("no resolvable canonical target base" in w["Message"] for w in warnings)


def test_pure_sugar_special_ids_keep_canonical_base_mass():
    mods = [_mod("Arp", ["A"], 176.0)]
    table = build_base_loss_masses(mods, BASE_MASSES)
    assert table["Arp"] == 135.0545


def test_real_catalog_resolves_every_entry_with_no_warnings():
    modifications = load_modifications(ROOT / "data" / "modifications.yaml")
    warnings = []
    table = build_base_loss_masses(modifications, BASE_MASSES, warnings)
    assert warnings == []
    assert len(table) == len(modifications) + 4  # + A, U, G, C
    # Spot-check a compound entry against its base-only counterpart.
    by_id = {mod.id: mod for mod in modifications}
    assert table["m1Am"] == table["m1A"]
    assert table["m4_4Cm"] == table["m4_4C"]
    assert table["Im"] == table["I"]
    assert table["Ym"] == table["Y"] == table["U"]
    assert "m1A" in by_id and "m4_4C" in by_id  # sanity: components exist
