from dataclasses import fields
from hashlib import sha256
from pathlib import Path

import pytest

from rna_masshunter.intact_rna_mass import (
    CcaPolicy,
    FivePrimeState,
    IntactRnaMassParameters,
    RnaTopology,
    ThreePrimeState,
    calculate_intact_rna_mass,
)
from rna_masshunter.masses import calculate_unmodified_rna_mass, load_base_masses

ROOT = Path(__file__).resolve().parent


def parameters(**overrides):
    values = {
        "five_prime_state": FivePrimeState.OH,
        "three_prime_state": ThreePrimeState.OH,
        "topology": RnaTopology.LINEAR,
        "cca_policy": CcaPolicy.AS_PROVIDED,
        "convert_t_to_u": False,
        "terminal_state_confirmed": False,
    }
    values.update(overrides)
    return IntactRnaMassParameters(**values)


def calculate(sequence, **overrides):
    return calculate_intact_rna_mass(sequence, parameters=parameters(**overrides))


def delta(left, right):
    left_counts = left.elemental_composition.to_dict()
    right_counts = right.elemental_composition.to_dict()
    return {
        element: left_counts.get(element, 0) - right_counts.get(element, 0)
        for element in sorted(set(left_counts) | set(right_counts))
        if left_counts.get(element, 0) != right_counts.get(element, 0)
    }


@pytest.mark.parametrize(
    "sequence,formula,expected_mass",
    [
        ("A", "C10H13N5O4", 267.09675391942),
        ("AC", "C19H25N8O11P1", 572.1380406548801),
        ("ACGU", "C38H48N15O26P3", 1223.21077771757),
    ],
)
def test_known_linear_oh_oh_compositions_and_masses(sequence, formula, expected_mass):
    result = calculate(sequence)
    assert result.formula == formula
    assert result.monoisotopic_neutral_mass == pytest.approx(expected_mass, abs=1e-9)


@pytest.mark.parametrize(
    "sequence,formula,expected_mass",
    [
        ("A", "C10H14N5O7P1", 347.06308480878),
        ("ACGU", "C38H49N15O29P4", 1303.17710860693),
    ],
)
def test_known_linear_five_prime_monophosphate_masses(sequence, formula, expected_mass):
    result = calculate(sequence, five_prime_state=FivePrimeState.MONOPHOSPHATE)
    assert result.formula == formula
    assert result.monoisotopic_neutral_mass == pytest.approx(expected_mass, abs=1e-9)


def test_lowercase_is_normalized_without_warning():
    result = calculate("acgu", terminal_state_confirmed=True)
    assert result.normalized_sequence == "ACGU"
    assert result.warnings == ()


def test_ascii_whitespace_is_removed():
    result = calculate("  A\tC\nG\rU\v\f")
    assert result.normalized_sequence == "ACGU"


@pytest.mark.parametrize("sequence", ["", " \n\t"])
def test_empty_sequence_is_rejected(sequence):
    with pytest.raises(ValueError, match="must not be empty"):
        calculate(sequence)


@pytest.mark.parametrize("sequence,bad", [("AT", "T"), ("AN", "N"), ("Am", "M"), ("A1", "1"), ("A*", r"\*")])
def test_noncanonical_symbols_are_rejected(sequence, bad):
    with pytest.raises(ValueError, match=bad):
        calculate(sequence)


def test_t_to_u_requires_explicit_option_and_records_provenance():
    result = calculate("atg", convert_t_to_u=True)
    assert result.normalized_sequence == "AUG"
    assert result.t_to_u_conversion_applied is True
    assert result.warnings == ("T_TO_U_CONVERSION_APPLIED", "TERMINAL_STATE_NOT_CONFIRMED")


def test_convert_option_without_t_records_no_conversion():
    result = calculate("AUG", convert_t_to_u=True, terminal_state_confirmed=True)
    assert result.t_to_u_conversion_applied is False
    assert result.warnings == ()


def test_non_ascii_whitespace_is_not_silently_removed():
    with pytest.raises(ValueError, match="invalid canonical RNA base"):
        calculate("A\u00a0C")


def test_input_string_is_unchanged_and_result_is_deterministic():
    source = " aC\nGU "
    snapshot = source
    first = calculate(source)
    second = calculate(source)
    assert source == snapshot
    assert first == second


def test_sequence_sha256_is_based_on_normalized_sequence():
    result = calculate(" a c g u ")
    assert result.sequence_sha256 == sha256(b"ACGU").hexdigest()


def test_cca_is_reported_as_provided():
    with_cca = calculate("ACGCCA")
    without_cca = calculate("ACG")
    assert with_cca.ends_with_cca is True
    assert without_cca.ends_with_cca is False
    assert with_cca.cca_policy == "AS_PROVIDED"


def test_cca_is_neither_added_nor_removed():
    assert calculate("ACG").normalized_sequence == "ACG"
    assert calculate("ACGCCA").normalized_sequence == "ACGCCA"


@pytest.mark.parametrize("sequence,bonds", [("A", 0), ("AC", 1), ("ACGU", 3)])
def test_linear_phosphodiester_bond_count_is_n_minus_one(sequence, bonds):
    assert calculate(sequence).phosphodiester_bond_count == bonds


def test_circular_phosphodiester_bond_count_is_n():
    assert calculate("AC", topology=RnaTopology.CIRCULAR).phosphodiester_bond_count == 2
    assert calculate("ACGU", topology=RnaTopology.CIRCULAR).phosphodiester_bond_count == 4


def test_formula_mass_matches_final_elemental_composition():
    result = calculate(
        "ACGU",
        five_prime_state=FivePrimeState.TRIPHOSPHATE,
        three_prime_state=ThreePrimeState.MONOPHOSPHATE,
        terminal_state_confirmed=True,
    )
    assert result.formula == result.elemental_composition.canonical_string()
    assert result.monoisotopic_neutral_mass == result.elemental_composition.exact_mass
    assert all(value >= 0 for value in result.elemental_composition.to_dict().values())


@pytest.mark.parametrize(
    "state,expected_delta",
    [
        (FivePrimeState.OH, {}),
        (FivePrimeState.MONOPHOSPHATE, {"H": 1, "O": 3, "P": 1}),
        (FivePrimeState.DIPHOSPHATE, {"H": 2, "O": 6, "P": 2}),
        (FivePrimeState.TRIPHOSPHATE, {"H": 3, "O": 9, "P": 3}),
    ],
)
def test_five_prime_states_are_direct_deltas_from_oh(state, expected_delta):
    reference = calculate("ACGU")
    candidate = calculate("ACGU", five_prime_state=state)
    assert delta(candidate, reference) == expected_delta


def test_five_prime_monophosphate_to_diphosphate_delta_is_hpo3():
    mono = calculate("ACGU", five_prime_state=FivePrimeState.MONOPHOSPHATE)
    diphosphate = calculate("ACGU", five_prime_state=FivePrimeState.DIPHOSPHATE)
    assert delta(diphosphate, mono) == {"H": 1, "O": 3, "P": 1}


def test_five_prime_monophosphate_to_triphosphate_delta_is_h2p2o6():
    mono = calculate("ACGU", five_prime_state=FivePrimeState.MONOPHOSPHATE)
    triphosphate = calculate("ACGU", five_prime_state=FivePrimeState.TRIPHOSPHATE)
    assert delta(triphosphate, mono) == {"H": 2, "O": 6, "P": 2}


@pytest.mark.parametrize(
    "state,expected_delta",
    [
        (ThreePrimeState.OH, {}),
        (ThreePrimeState.MONOPHOSPHATE, {"H": 1, "O": 3, "P": 1}),
        (ThreePrimeState.CYCLIC_PHOSPHATE, {"H": -1, "O": 2, "P": 1}),
    ],
)
def test_three_prime_states_are_direct_deltas_from_oh(state, expected_delta):
    reference = calculate("ACGU")
    candidate = calculate("ACGU", three_prime_state=state)
    assert delta(candidate, reference) == expected_delta


def test_three_prime_monophosphate_to_cyclic_delta_is_minus_water():
    mono = calculate("ACGU", three_prime_state=ThreePrimeState.MONOPHOSPHATE)
    cyclic = calculate("ACGU", three_prime_state=ThreePrimeState.CYCLIC_PHOSPHATE)
    assert delta(cyclic, mono) == {"H": -2, "O": -1}


def test_circular_minus_linear_delta_is_minus_water():
    linear = calculate("ACGU")
    circular = calculate("ACGU", topology=RnaTopology.CIRCULAR)
    assert circular.topology == "CIRCULAR"
    assert delta(circular, linear) == {"H": -2, "O": -1}


def test_single_residue_circular_rna_is_rejected():
    with pytest.raises(ValueError, match="at least two residues"):
        calculate("A", topology=RnaTopology.CIRCULAR)


@pytest.mark.parametrize(
    "five,three",
    [
        (FivePrimeState.MONOPHOSPHATE, ThreePrimeState.OH),
        (FivePrimeState.OH, ThreePrimeState.MONOPHOSPHATE),
        (FivePrimeState.OH, ThreePrimeState.CYCLIC_PHOSPHATE),
    ],
)
def test_circular_rna_rejects_independent_terminal_states(five, three):
    with pytest.raises(ValueError, match="does not have independent"):
        calculate("AC", topology=RnaTopology.CIRCULAR, five_prime_state=five, three_prime_state=three)


@pytest.mark.parametrize(
    "field,value,message",
    [
        ("five_prime_state", "UNKNOWN", "five_prime_state"),
        ("three_prime_state", "UNKNOWN", "three_prime_state"),
        ("topology", "BRANCHED", "topology"),
        ("cca_policy", "AUTO_ADD", "cca_policy"),
    ],
)
def test_unknown_chemistry_parameters_are_rejected(field, value, message):
    with pytest.raises(ValueError, match=message):
        calculate("AC", **{field: value})


def test_nonboolean_provenance_parameters_are_rejected():
    with pytest.raises(ValueError, match="convert_t_to_u"):
        calculate_intact_rna_mass("AC", parameters=parameters(convert_t_to_u="yes"))
    with pytest.raises(ValueError, match="terminal_state_confirmed"):
        calculate_intact_rna_mass("AC", parameters=parameters(terminal_state_confirmed="yes"))


def test_terminal_state_confirmation_is_provenance_not_calculation_gate():
    unconfirmed = calculate("AC")
    confirmed = calculate("AC", terminal_state_confirmed=True)
    assert unconfirmed.monoisotopic_neutral_mass == confirmed.monoisotopic_neutral_mass
    assert unconfirmed.terminal_state_confirmed is False
    assert unconfirmed.warnings == ("TERMINAL_STATE_NOT_CONFIRMED",)
    assert confirmed.terminal_state_confirmed is True
    assert confirmed.warnings == ()


def test_mass_type_and_formula_provenance_are_fixed_and_deterministic():
    result = calculate("ACGU")
    assert result.theoretical_mass_type == "MONOISOTOPIC_NEUTRAL"
    assert result.cca_policy == "AS_PROVIDED"
    assert result.formula == "C38H48N15O26P3"


def test_result_has_no_formal_score_fields():
    names = {field.name for field in fields(type(calculate("A")))}
    assert not any("formal" in name.lower() or "ranking" in name.lower() for name in names)


@pytest.mark.parametrize(
    "sequence,legacy_expected",
    [("A", 347.063064684), ("AC", 652.104364684), ("ACGU", 1303.177064684)],
)
def test_legacy_default_is_approximately_new_five_prime_monophosphate(sequence, legacy_expected):
    base_masses = load_base_masses(ROOT / "data" / "base_masses.yaml")
    legacy = calculate_unmodified_rna_mass(sequence, base_masses, terminal_form="default")
    new_oh = calculate(sequence)
    new_monophosphate = calculate(sequence, five_prime_state=FivePrimeState.MONOPHOSPHATE)
    assert legacy == pytest.approx(legacy_expected, abs=1e-12)
    assert legacy == pytest.approx(new_monophosphate.monoisotopic_neutral_mass, abs=0.001)
    assert legacy - new_oh.monoisotopic_neutral_mass == pytest.approx(79.96633088936, abs=0.001)


def test_valid_string_enum_values_are_normalized_at_api_boundary():
    result = calculate_intact_rna_mass(
        "AC",
        parameters=IntactRnaMassParameters(
            five_prime_state="monophosphate",
            three_prime_state="oh",
            topology="linear",
            cca_policy="as_provided",
        ),
    )
    assert result.five_prime_state == "MONOPHOSPHATE"
    assert result.three_prime_state == "OH"
    assert result.topology == "LINEAR"
