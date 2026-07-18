from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
from pathlib import Path

import pytest

from rna_masshunter.cca_tail_state import (
    CCA_TAIL_STATE_APPLIED_TO_CANDIDATE_FILTERING,
    CCA_TAIL_STATE_APPLIED_TO_FINAL_CONSENSUS,
    CCA_TAIL_STATE_APPLIED_TO_FORMAL_SCORE,
    CCA_TAIL_STATE_APPLIED_TO_RANKING,
    CCATailState,
    CCATailStatus,
    RegisteredSequenceCCAMode,
    build_cca_tail_variant,
    calculate_cca_tail_variant_mass,
    derive_core_sequence,
    generate_cca_tail_variants,
)
from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.intact_rna_mass import (
    FivePrimeState,
    IntactRnaMassParameters,
    ThreePrimeState,
    calculate_intact_rna_mass,
)
from rna_masshunter.sciex_sample_manifest import (
    SAMPLE_MANIFEST_APPLIED_TO_CANDIDATE_FILTERING,
    SAMPLE_MANIFEST_APPLIED_TO_FINAL_CONSENSUS,
    SAMPLE_MANIFEST_APPLIED_TO_FORMAL_SCORE,
    SAMPLE_MANIFEST_APPLIED_TO_RANKING,
    get_measurements_for_sample,
    get_rna_identity,
    load_sciex_sample_manifest,
    resolve_measurement_identity,
)
from rna_masshunter.structure_fragment import RNA_RESIDUE_COMPOSITIONS

ROOT = Path(__file__).parent
MANIFEST_PATH = ROOT / "data" / "sciex_sample_manifest.yaml"
LEU_UAA_SEQUENCE = "GCGAGGGUUGCCCAGCCAGGCCAAAGGCGCCAGACUUAAGAUCUGGUAUCGAAGGAUUUCGUGGGUUCGAAUCCCACCCCUCGCA"
LEU_UAG_SEQUENCE = "GCGAGGGUUGCCCAGCUAGGUCAAAGGCGAUGGGCUUAGGACCCAUUUUCGUAGGAAUUCGUGCGUUCGAAUCGCACCCCUCGCA"
GLU_UUC_SEQUENCE = "GCUCCGGUAGUGUAGUCCGGCCAAUCAUUCCGGCCUUUCGAGCCGAAGACUCGGGUUCGAAUCCCGGCCGGAGCACCA"
TAIL_EXPECTATIONS = {
    CCATailState.NONE: ("", "CCA", ("C", "C", "A"), 3),
    CCATailState.C: ("C", "CA", ("C", "A"), 2),
    CCATailState.CC: ("CC", "A", ("A",), 1),
    CCATailState.CCA: ("CCA", "", (), 0),
}


@pytest.mark.parametrize("state", list(TAIL_EXPECTATIONS))
def test_cca_completion_behavior(state):
    tail, missing, additions, count = TAIL_EXPECTATIONS[state]
    variant = build_cca_tail_variant("ACGU", RegisteredSequenceCCAMode.EXCLUDES_CCA, state)
    assert variant.candidate_cca_tail_sequence == tail
    assert variant.candidate_cca_tail_length == len(tail)
    assert variant.missing_cca_suffix == missing
    assert variant.required_added_nucleotides == additions
    assert variant.cca_completion_step_count == count
    assert variant.cca_completion_required is bool(missing)


@pytest.mark.parametrize("sequence", [LEU_UAA_SEQUENCE, LEU_UAG_SEQUENCE])
def test_excludes_cca_core_is_registered_sequence(sequence):
    before = sequence
    core = derive_core_sequence(sequence, RegisteredSequenceCCAMode.EXCLUDES_CCA)
    assert core == sequence
    assert sequence == before
    assert core.endswith("GCA")


@pytest.mark.parametrize("sequence", [LEU_UAA_SEQUENCE, LEU_UAG_SEQUENCE])
def test_excludes_cca_default_candidates_are_fixed_and_complete(sequence):
    variants = generate_cca_tail_variants(sequence, RegisteredSequenceCCAMode.EXCLUDES_CCA)
    assert [item.candidate_cca_tail_state for item in variants] == [
        CCATailState.NONE,
        CCATailState.C,
        CCATailState.CC,
        CCATailState.CCA,
    ]
    assert [item.complete_candidate_sequence_length for item in variants] == [85, 86, 87, 88]
    for variant, suffix in zip(variants, ("GCA", "GCAC", "GCACC", "GCACCA"), strict=True):
        assert variant.complete_candidate_sequence.endswith(suffix)


@pytest.mark.parametrize("sequence", [LEU_UAA_SEQUENCE, LEU_UAG_SEQUENCE])
def test_excludes_cca_gca_suffix_is_not_misread_as_tail(sequence):
    variants = generate_cca_tail_variants(sequence, RegisteredSequenceCCAMode.EXCLUDES_CCA)
    assert variants[0].core_sequence.endswith("GCA")
    assert variants[0].candidate_cca_tail_sequence == ""
    assert variants[-1].complete_candidate_sequence == sequence + "CCA"


def test_excludes_mode_does_not_split_incidental_cca_suffix():
    sequence = "ACGUCCA"
    assert derive_core_sequence(sequence, RegisteredSequenceCCAMode.EXCLUDES_CCA) == sequence


def test_includes_complete_cca_derives_core():
    assert GLU_UUC_SEQUENCE.endswith("CCA")
    assert derive_core_sequence(
        GLU_UUC_SEQUENCE, RegisteredSequenceCCAMode.INCLUDES_COMPLETE_CCA
    ) == GLU_UUC_SEQUENCE[:-3]


def test_includes_complete_cca_default_candidates_are_cca_then_cc():
    variants = generate_cca_tail_variants(
        GLU_UUC_SEQUENCE, RegisteredSequenceCCAMode.INCLUDES_COMPLETE_CCA
    )
    assert [item.candidate_cca_tail_state for item in variants] == [
        CCATailState.CCA,
        CCATailState.CC,
    ]
    assert variants[0].complete_candidate_sequence == GLU_UUC_SEQUENCE
    assert variants[1].complete_candidate_sequence == GLU_UUC_SEQUENCE[:-1]


def test_includes_complete_cca_explicit_c_and_none_candidates():
    variants = generate_cca_tail_variants(
        GLU_UUC_SEQUENCE,
        RegisteredSequenceCCAMode.INCLUDES_COMPLETE_CCA,
        candidate_states=(CCATailState.C, CCATailState.NONE),
    )
    assert [item.candidate_cca_tail_state for item in variants] == [
        CCATailState.NONE,
        CCATailState.C,
    ]
    assert variants[0].complete_candidate_sequence == GLU_UUC_SEQUENCE[:-3]
    assert variants[1].complete_candidate_sequence == GLU_UUC_SEQUENCE[:-2]


def test_includes_complete_cca_never_doubles_tail():
    variants = generate_cca_tail_variants(
        GLU_UUC_SEQUENCE,
        RegisteredSequenceCCAMode.INCLUDES_COMPLETE_CCA,
        candidate_states=tuple(CCATailState),
    )
    assert all(not item.complete_candidate_sequence.endswith("CCACCA") for item in variants)
    assert max(item.complete_candidate_sequence_length for item in variants) == len(GLU_UUC_SEQUENCE)


def test_includes_mode_rejects_non_cca_suffix():
    with pytest.raises(ValueError, match="must end with CCA"):
        derive_core_sequence("ACGU", RegisteredSequenceCCAMode.INCLUDES_COMPLETE_CCA)


@pytest.mark.parametrize("candidate_states", [None, (CCATailState.CCA,)])
def test_unknown_mode_never_infers_core_from_suffix(candidate_states):
    with pytest.raises(ValueError, match="UNKNOWN"):
        generate_cca_tail_variants("ACGUCCA", RegisteredSequenceCCAMode.UNKNOWN, candidate_states)


def test_invalid_mode_is_rejected():
    with pytest.raises(ValueError, match="unknown registered_sequence_cca_mode"):
        derive_core_sequence("ACGU", "INFER_FROM_SUFFIX")


def test_invalid_candidate_state_is_rejected():
    with pytest.raises(ValueError, match="unknown candidate_cca_tail_state"):
        build_cca_tail_variant("ACGU", RegisteredSequenceCCAMode.EXCLUDES_CCA, "CA")


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [
        (" ac gu \n", "ACGU"),
        ("acgu", "ACGU"),
    ],
)
def test_sequence_normalization_matches_canonical_rna_rules(sequence, expected):
    variant = build_cca_tail_variant(
        sequence, RegisteredSequenceCCAMode.EXCLUDES_CCA, CCATailState.NONE
    )
    assert variant.registered_sequence == expected
    assert variant.core_sequence == expected


@pytest.mark.parametrize("sequence", ["", " \n\t", "ACGT", "ACGX", "ACG-"])
def test_invalid_or_empty_sequence_is_rejected_without_t_to_u(sequence):
    with pytest.raises(ValueError, match="empty|invalid canonical RNA base"):
        derive_core_sequence(sequence, RegisteredSequenceCCAMode.EXCLUDES_CCA)


def test_non_string_sequence_is_rejected():
    with pytest.raises(TypeError, match="must be a string"):
        derive_core_sequence(["A", "C"], RegisteredSequenceCCAMode.EXCLUDES_CCA)


def test_hashes_and_repeated_generation_are_deterministic():
    first = generate_cca_tail_variants(LEU_UAA_SEQUENCE, RegisteredSequenceCCAMode.EXCLUDES_CCA)
    second = generate_cca_tail_variants(LEU_UAA_SEQUENCE, RegisteredSequenceCCAMode.EXCLUDES_CCA)
    assert first == second
    for variant in first:
        assert variant.registered_sequence_sha256 == sha256(LEU_UAA_SEQUENCE.encode("ascii")).hexdigest()
        assert variant.core_sequence_sha256 == sha256(variant.core_sequence.encode("ascii")).hexdigest()
        assert variant.complete_candidate_sequence_sha256 == sha256(
            variant.complete_candidate_sequence.encode("ascii")
        ).hexdigest()


def test_explicit_candidate_order_is_fixed_not_caller_or_enum_order():
    variants = generate_cca_tail_variants(
        LEU_UAA_SEQUENCE,
        RegisteredSequenceCCAMode.EXCLUDES_CCA,
        candidate_states=(CCATailState.CCA, CCATailState.C, CCATailState.NONE),
    )
    assert [item.candidate_cca_tail_state for item in variants] == [
        CCATailState.NONE,
        CCATailState.C,
        CCATailState.CCA,
    ]


def test_empty_or_string_candidate_state_collection_is_rejected():
    with pytest.raises(ValueError, match="must not be empty"):
        generate_cca_tail_variants(
            LEU_UAA_SEQUENCE, RegisteredSequenceCCAMode.EXCLUDES_CCA, ()
        )
    with pytest.raises(TypeError, match="iterable"):
        generate_cca_tail_variants(
            LEU_UAA_SEQUENCE, RegisteredSequenceCCAMode.EXCLUDES_CCA, "CCA"
        )


def test_candidate_mass_deltas_equal_polymer_residue_compositions():
    parameters = IntactRnaMassParameters(
        five_prime_state=FivePrimeState.OH,
        three_prime_state=ThreePrimeState.OH,
        terminal_state_confirmed=False,
    )
    variants = generate_cca_tail_variants(
        LEU_UAA_SEQUENCE, RegisteredSequenceCCAMode.EXCLUDES_CCA
    )
    masses = [calculate_cca_tail_variant_mass(item, parameters) for item in variants]
    values = [item.monoisotopic_neutral_mass for item in masses]
    deltas = [values[index + 1] - values[index] for index in range(3)]
    cytidine_residue_mass = ElementalComposition(RNA_RESIDUE_COMPOSITIONS["C"]).exact_mass
    adenosine_residue_mass = ElementalComposition(RNA_RESIDUE_COMPOSITIONS["A"]).exact_mass
    assert values == sorted(values)
    assert deltas[0] == pytest.approx(cytidine_residue_mass, abs=1e-9)
    assert deltas[1] == pytest.approx(cytidine_residue_mass, abs=1e-9)
    assert deltas[2] == pytest.approx(adenosine_residue_mass, abs=1e-9)
    assert deltas[0] == pytest.approx(deltas[1], abs=1e-9)


@pytest.mark.parametrize(
    ("state", "direct_sequence"),
    [
        (CCATailState.CCA, GLU_UUC_SEQUENCE),
        (CCATailState.CC, GLU_UUC_SEQUENCE[:-1]),
    ],
)
def test_glu_candidate_mass_matches_direct_complete_sequence(state, direct_sequence):
    parameters = IntactRnaMassParameters(terminal_state_confirmed=False)
    variant = build_cca_tail_variant(
        GLU_UUC_SEQUENCE, RegisteredSequenceCCAMode.INCLUDES_COMPLETE_CCA, state
    )
    candidate_mass = calculate_cca_tail_variant_mass(variant, parameters)
    direct_mass = calculate_intact_rna_mass(direct_sequence, parameters=parameters)
    assert candidate_mass.monoisotopic_neutral_mass == direct_mass.monoisotopic_neutral_mass
    assert candidate_mass.intact_mass_result.elemental_composition == direct_mass.elemental_composition


def test_mass_helper_does_not_mutate_or_conflate_terminal_parameters():
    parameters = IntactRnaMassParameters(
        five_prime_state=FivePrimeState.MONOPHOSPHATE,
        three_prime_state=ThreePrimeState.CYCLIC_PHOSPHATE,
        terminal_state_confirmed=False,
    )
    before = asdict(parameters)
    variant = build_cca_tail_variant(
        LEU_UAG_SEQUENCE, RegisteredSequenceCCAMode.EXCLUDES_CCA, CCATailState.CCA
    )
    result = calculate_cca_tail_variant_mass(variant, parameters)
    assert asdict(parameters) == before
    assert result.intact_mass_parameters is parameters
    assert result.intact_mass_result.five_prime_state == "MONOPHOSPHATE"
    assert result.intact_mass_result.three_prime_state == "CYCLIC_PHOSPHATE"
    assert result.charge_adjustment_applied is False
    assert result.mz_conversion_applied is False


@pytest.mark.parametrize("state", list(CCATailState))
def test_all_candidates_are_mass_only_without_biological_assignment(state):
    variant = build_cca_tail_variant(
        LEU_UAA_SEQUENCE, RegisteredSequenceCCAMode.EXCLUDES_CCA, state
    )
    assert variant.mass_match_only is True
    assert variant.cca_state_confirmed is False
    assert variant.repair_intermediate_assigned is False
    assert variant.biological_cause_assigned is False
    assert variant.rnase_t_assigned is False
    assert variant.structure_identity_assigned is False
    assert variant.formal_evidence is False


def test_manifest_registered_and_sample_cca_states_are_explicit():
    manifest = load_sciex_sample_manifest(MANIFEST_PATH)
    uaa = get_rna_identity(manifest, "TRNA_LEU_UAA")
    uag = get_rna_identity(manifest, "TRNA_LEU_UAG")
    glu = get_rna_identity(manifest, "TRNA_GLU_UUC")
    assert (uaa.registered_sequence_cca_mode, uaa.registered_cca_tail_state) == (
        RegisteredSequenceCCAMode.EXCLUDES_CCA,
        CCATailState.NONE,
    )
    assert (uag.registered_sequence_cca_mode, uag.registered_cca_tail_state) == (
        RegisteredSequenceCCAMode.EXCLUDES_CCA,
        CCATailState.NONE,
    )
    assert (glu.registered_sequence_cca_mode, glu.registered_cca_tail_state) == (
        RegisteredSequenceCCAMode.INCLUDES_COMPLETE_CCA,
        CCATailState.CCA,
    )
    assert uaa.sequence_length == uag.sequence_length == 85
    assert glu.sequence == GLU_UUC_SEQUENCE
    assert glu.sequence.endswith("CCA") and not glu.sequence.endswith("CCACCA")
    assert all(sample.sample_cca_tail_state is None for sample in manifest.samples)
    assert all(sample.sample_cca_tail_status is CCATailStatus.UNKNOWN for sample in manifest.samples)


def test_manifest_measurements_pairing_and_sequence_hash_linkage_are_unchanged():
    manifest = load_sciex_sample_manifest(MANIFEST_PATH)
    assert len(manifest.measurements) == 6
    expected = {
        "LEU_UAA_WT": {"LEU_UAA_WT_FULL", "LEU_UAA_WT_T1"},
        "GLU_UUC_WT": {"GLU_UUC_WT_FULL", "GLU_UUC_WT_P1_AP"},
        "LEU_UAG_WT": {"LEU_UAG_WT_FULL", "LEU_UAG_WT_T1"},
    }
    for sample_id, measurement_ids in expected.items():
        actual = get_measurements_for_sample(manifest, sample_id)
        assert {item.measurement_id for item in actual} == measurement_ids
        hashes = {
            resolve_measurement_identity(manifest, item.measurement_id).expected_sequence_sha256
            for item in actual
        }
        assert len(hashes) == 1 and None not in hashes


def test_cca_and_manifest_formal_flags_are_false():
    assert CCA_TAIL_STATE_APPLIED_TO_FORMAL_SCORE is False
    assert CCA_TAIL_STATE_APPLIED_TO_RANKING is False
    assert CCA_TAIL_STATE_APPLIED_TO_CANDIDATE_FILTERING is False
    assert CCA_TAIL_STATE_APPLIED_TO_FINAL_CONSENSUS is False
    assert SAMPLE_MANIFEST_APPLIED_TO_FORMAL_SCORE is False
    assert SAMPLE_MANIFEST_APPLIED_TO_RANKING is False
    assert SAMPLE_MANIFEST_APPLIED_TO_CANDIDATE_FILTERING is False
    assert SAMPLE_MANIFEST_APPLIED_TO_FINAL_CONSENSUS is False
