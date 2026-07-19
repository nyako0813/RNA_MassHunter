from __future__ import annotations

from dataclasses import asdict

import pytest

from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.intact_rna_average_mass import (
    AVERAGE_ATOMIC_MASSES,
    AVERAGE_ATOMIC_MASS_CONSTANT_SET,
    AVERAGE_ATOMIC_MASS_CONSTANT_SOURCE,
    AVERAGE_ATOMIC_MASS_CONSTANT_VERSION,
    PROTON_MASS_DA,
    ComparisonReferenceRole,
    MassDisplaySpecies,
    TheoreticalMassDefinition,
    calculate_average_neutral_mass_from_composition,
    calculate_intact_rna_average_mass,
)
from rna_masshunter.intact_rna_mass import (
    FivePrimeState,
    IntactRnaMassParameters,
    RnaTopology,
    ThreePrimeState,
    calculate_intact_rna_mass,
)

GLU_UUC = "GCUCCGGUAGUGUAGUCCGGCCAAUCAUUCCGGCCUUUCGAGCCGAAGACUCGGGUUCGAAUCCCGGCCGGAGCACCA"
LEU_UAA = "GCGAGGGUUGCCCAGCCAGGCCAAAGGCGCCAGACUUAAGAUCUGGUAUCGAAGGAUUUCGUGGGUUCGAAUCCCACCCCUCGCA"
LEU_UAG = "GCGAGGGUUGCCCAGCUAGGUCAAAGGCGAUGGGCUUAGGACCCAUUUUCGUAGGAAUUCGUGCGUUCGAAUCGCACCCCUCGCA"
PARAMETERS = IntactRnaMassParameters(five_prime_state=FivePrimeState.MONOPHOSPHATE)


def result(sequence=GLU_UUC, **kwargs):
    return calculate_intact_rna_average_mass(sequence, parameters=PARAMETERS, **kwargs)


def test_fixed_average_atomic_mass_constants_and_metadata():
    assert set(AVERAGE_ATOMIC_MASSES) == {"C", "H", "N", "O", "P", "S"}
    assert all(value > 0 for value in AVERAGE_ATOMIC_MASSES.values())
    assert AVERAGE_ATOMIC_MASS_CONSTANT_SET == "IUPAC_CONVENTIONAL_ATOMIC_WEIGHTS_FIXED"
    assert AVERAGE_ATOMIC_MASS_CONSTANT_VERSION == "RNA_MASSHUNTER_2026_01"
    assert "Fixed conventional" in AVERAGE_ATOMIC_MASS_CONSTANT_SOURCE


def test_composition_calculation_is_deterministic_and_does_not_mutate_input():
    composition = {"C": 2, "H": 6, "O": 1}
    before = composition.copy()
    first = calculate_average_neutral_mass_from_composition(composition)
    second = calculate_average_neutral_mass_from_composition(composition)
    assert first == second == 2 * 12.0107 + 6 * 1.00794 + 15.9994
    assert composition == before


def test_elemental_composition_is_supported_and_empty_is_zero():
    assert calculate_average_neutral_mass_from_composition(ElementalComposition({"P": 1})) == AVERAGE_ATOMIC_MASSES["P"]
    assert calculate_average_neutral_mass_from_composition({}) == 0.0


@pytest.mark.parametrize("composition", [{"Se": 1}, {"X": 1}])
def test_unknown_average_mass_element_is_rejected(composition):
    with pytest.raises(ValueError, match="unsupported average-mass element"):
        calculate_average_neutral_mass_from_composition(composition)


def test_negative_and_noninteger_element_counts_are_rejected():
    with pytest.raises(ValueError, match="negative element count"):
        calculate_average_neutral_mass_from_composition({"C": -1})
    with pytest.raises(ValueError, match="must be an integer"):
        calculate_average_neutral_mass_from_composition({"C": 1.5})


def test_average_calculation_is_independent_of_exact_mass_property():
    composition = ElementalComposition({"C": 10, "H": 20})
    average = calculate_average_neutral_mass_from_composition(composition)
    assert average != composition.exact_mass


def test_complete_sequence_terminal_chemistry_and_topology_are_reused_from_mono_layer():
    calculated = result()
    mono = calculated.monoisotopic_result
    assert mono.normalized_sequence == GLU_UUC
    assert mono.sequence_length == 78
    assert mono.ends_with_cca is True
    assert mono.five_prime_state == "MONOPHOSPHATE"
    assert mono.three_prime_state == "OH"
    assert mono.topology == "LINEAR"
    assert calculated.elemental_composition.to_dict() == {"C": 739, "H": 923, "N": 293, "O": 548, "P": 78}


def test_cca_state_changes_average_neutral_mass():
    full = result(GLU_UUC).average_neutral_molecular_mass_m
    cc = result(GLU_UUC[:-1]).average_neutral_molecular_mass_m
    assert full > cc
    assert full != cc


def test_five_prime_three_prime_and_topology_change_composition_based_average_mass():
    base = result().average_neutral_molecular_mass_m
    five_oh = calculate_intact_rna_average_mass(
        GLU_UUC, parameters=IntactRnaMassParameters(five_prime_state=FivePrimeState.OH)
    ).average_neutral_molecular_mass_m
    three_p = calculate_intact_rna_average_mass(
        GLU_UUC,
        parameters=IntactRnaMassParameters(
            five_prime_state=FivePrimeState.MONOPHOSPHATE,
            three_prime_state=ThreePrimeState.MONOPHOSPHATE,
        ),
    ).average_neutral_molecular_mass_m
    circular = calculate_intact_rna_average_mass(
        GLU_UUC,
        parameters=IntactRnaMassParameters(topology=RnaTopology.CIRCULAR),
    ).average_neutral_molecular_mass_m
    assert len({base, five_oh, three_p, circular}) == 4


def test_neutral_m_and_diagnostic_species_are_separate_and_metadata_is_explicit():
    calculated = result()
    assert calculated.average_m_plus_h == calculated.average_neutral_molecular_mass_m + PROTON_MASS_DA
    assert calculated.average_m_minus_h == calculated.average_neutral_molecular_mass_m - PROTON_MASS_DA
    assert calculated.average_m_plus_h != calculated.average_m_minus_h
    assert [item.theoretical_mass_definition for item in calculated.references] == [
        TheoreticalMassDefinition.AVERAGE_NEUTRAL_M,
        TheoreticalMassDefinition.AVERAGE_M_PLUS_H,
        TheoreticalMassDefinition.AVERAGE_M_MINUS_H,
        TheoreticalMassDefinition.MONOISOTOPIC_NEUTRAL_M,
    ]
    neutral, plus_h, minus_h, mono = calculated.references
    assert neutral.mass_display_species is MassDisplaySpecies.M
    assert neutral.protonation_adjustment_applied is False
    assert neutral.deprotonation_adjustment_applied is False
    assert plus_h.protonation_adjustment_applied is True
    assert minus_h.deprotonation_adjustment_applied is True
    assert mono.comparison_role is ComparisonReferenceRole.MONOISOTOPIC_DIAGNOSTIC_ONLY
    assert all(not item.charge_adjustment_applied for item in calculated.references)
    assert all(not item.electron_mass_adjustment_applied for item in calculated.references)
    assert all(not item.mz_conversion_applied for item in calculated.references)


def test_ion_mode_and_charge_descriptors_do_not_change_neutral_m():
    baseline = result().average_neutral_molecular_mass_m
    assert result(ion_mode="positive").average_neutral_molecular_mass_m == baseline
    assert result(ion_mode="negative").average_neutral_molecular_mass_m == baseline
    assert result(charge_state=-12).average_neutral_molecular_mass_m == baseline


def test_glu_uuc_reference_values_and_apex_diagnostics():
    calculated = result()
    assert calculated.monoisotopic_neutral_mass == pytest.approx(25082.289835, abs=1e-6)
    # Mongo Oligo is an external validation reference, never a production constant.
    mongo_average_m = 25094.067
    assert calculated.average_neutral_molecular_mass_m == pytest.approx(mongo_average_m, abs=0.3)
    assert calculated.average_neutral_molecular_mass_m > calculated.monoisotopic_neutral_mass
    assert calculated.average_neutral_molecular_mass_m - calculated.monoisotopic_neutral_mass == pytest.approx(11.533820397, abs=1e-6)
    apex = 25326.0
    neutral_delta = apex - calculated.average_neutral_molecular_mass_m
    plus_delta = apex - calculated.average_m_plus_h
    minus_delta = apex - calculated.average_m_minus_h
    assert neutral_delta == pytest.approx(232.176344156, abs=1e-6)
    assert plus_delta == pytest.approx(neutral_delta - PROTON_MASS_DA, abs=1e-9)
    assert minus_delta == pytest.approx(neutral_delta + PROTON_MASS_DA, abs=1e-9)
    assert len({neutral_delta, plus_delta, minus_delta}) == 3


def test_leu_average_references_and_deterministic_difference():
    uaa = result(LEU_UAA + "CCA")
    uag = result(LEU_UAG + "CCA")
    assert uaa.average_neutral_molecular_mass_m == pytest.approx(28370.6501, abs=0.2)
    assert uag.average_neutral_molecular_mass_m == pytest.approx(28365.5801, abs=0.2)
    assert uaa.average_neutral_molecular_mass_m - uag.average_neutral_molecular_mass_m == pytest.approx(5.07102, abs=1e-9)
    assert uaa.monoisotopic_neutral_mass == calculate_intact_rna_mass(LEU_UAA + "CCA", parameters=PARAMETERS).monoisotopic_neutral_mass
    assert uag.monoisotopic_neutral_mass == calculate_intact_rna_mass(LEU_UAG + "CCA", parameters=PARAMETERS).monoisotopic_neutral_mass


def test_result_is_repeatable_and_input_parameters_are_unchanged():
    before = asdict(PARAMETERS)
    assert result() == result()
    assert asdict(PARAMETERS) == before
