from __future__ import annotations

from dataclasses import FrozenInstanceError, asdict, replace
from pathlib import Path

import pytest

from rna_masshunter.cca_tail_state import (
    CCATailState,
    RegisteredSequenceCCAMode,
    generate_cca_tail_variants,
)
from rna_masshunter.intact_rna_candidate_generation import (
    CANDIDATE_GENERATION_APPLIED_TO_CANDIDATE_FILTERING,
    CANDIDATE_GENERATION_APPLIED_TO_FINAL_CONSENSUS,
    CANDIDATE_GENERATION_APPLIED_TO_FORMAL_SCORE,
    CANDIDATE_GENERATION_APPLIED_TO_RANKING,
    PRIMARY_DEFAULT,
    SECONDARY_TERMINAL_DIAGNOSTIC,
    AssumptionStatus,
    CandidateCategory,
    CandidateGenerationError,
    CandidateSetName,
    build_unmodified_intact_candidate,
    generate_candidates_for_measurement,
)
from rna_masshunter.intact_rna_mass import (
    CcaPolicy,
    FivePrimeState,
    IntactRnaMassParameters,
    RnaTopology,
    ThreePrimeState,
    calculate_intact_rna_mass,
)
from rna_masshunter.sciex_sample_manifest import (
    SequenceStatus,
    get_measurement,
    get_rna_identity,
    load_sciex_sample_manifest,
)

ROOT = Path(__file__).parent
MANIFEST_PATH = ROOT / "data" / "sciex_sample_manifest.yaml"


@pytest.fixture
def manifest():
    return load_sciex_sample_manifest(MANIFEST_PATH)


def candidate_states(candidates):
    return [item.cca_tail_state for item in candidates]


def test_glu_full_routes_correct_identity_and_default_order(manifest):
    candidates = generate_candidates_for_measurement(manifest, "GLU_UUC_WT_FULL")
    assert len(candidates) == 4
    assert all(item.rna_identity_id == "TRNA_GLU_UUC" for item in candidates)
    assert candidate_states(candidates) == [
        CCATailState.CCA,
        CCATailState.CC,
        CCATailState.CCA,
        CCATailState.CC,
    ]
    assert [item.candidate_set_name for item in candidates] == [
        CandidateSetName.PRIMARY_DEFAULT,
        CandidateSetName.PRIMARY_DEFAULT,
        CandidateSetName.SECONDARY_TERMINAL_DIAGNOSTIC,
        CandidateSetName.SECONDARY_TERMINAL_DIAGNOSTIC,
    ]
    assert [item.complete_sequence_length for item in candidates] == [78, 77, 78, 77]
    assert all(not item.complete_sequence.endswith("CCACCA") for item in candidates)


@pytest.mark.parametrize(
    ("measurement_id", "rna_identity_id", "expected_hash"),
    [
        (
            "LEU_UAA_WT_FULL",
            "TRNA_LEU_UAA",
            "71664e8092c4c48c9c31fbebd57ec9db5958f4ba29bfcf424ffc5a8fce26e72d",
        ),
        (
            "LEU_UAG_WT_FULL",
            "TRNA_LEU_UAG",
            "bb309562720cf5b7795103d75325791af4b9c27eaf03c4320cd48912e50d0dd6",
        ),
    ],
)
def test_leu_full_routes_identity_and_eight_candidates(
    manifest, measurement_id, rna_identity_id, expected_hash
):
    candidates = generate_candidates_for_measurement(manifest, measurement_id)
    assert len(candidates) == 8
    assert all(item.rna_identity_id == rna_identity_id for item in candidates)
    assert candidate_states(candidates[:4]) == list(CCATailState)
    assert candidate_states(candidates[4:]) == list(CCATailState)
    assert [item.complete_sequence_length for item in candidates[:4]] == [85, 86, 87, 88]
    assert [item.complete_sequence_length for item in candidates[4:]] == [85, 86, 87, 88]
    assert all(item.registered_sequence_sha256 == expected_hash for item in candidates)


@pytest.mark.parametrize(
    ("measurement_id", "primary_count"),
    [
        ("GLU_UUC_WT_FULL", 2),
        ("LEU_UAA_WT_FULL", 4),
        ("LEU_UAG_WT_FULL", 4),
    ],
)
def test_secondary_toggle_keeps_only_primary(manifest, measurement_id, primary_count):
    candidates = generate_candidates_for_measurement(
        manifest,
        measurement_id,
        include_secondary_terminal_state=False,
    )
    assert len(candidates) == primary_count
    assert all(item.candidate_set_name is CandidateSetName.PRIMARY_DEFAULT for item in candidates)
    assert all(item.candidate_priority == 1 for item in candidates)


def test_primary_terminal_default_is_assumed_5p_3oh_linear(manifest):
    candidates = generate_candidates_for_measurement(manifest, "LEU_UAA_WT_FULL")[:4]
    for candidate in candidates:
        assert candidate.five_prime_state is FivePrimeState.MONOPHOSPHATE
        assert candidate.five_prime_state_status is AssumptionStatus.ASSUMED
        assert candidate.three_prime_state is ThreePrimeState.OH
        assert candidate.three_prime_state_status is AssumptionStatus.ASSUMED
        assert candidate.topology is RnaTopology.LINEAR
        assert candidate.topology_status is AssumptionStatus.ASSUMED
        assert candidate.terminal_state_confirmed is False


def test_secondary_terminal_diagnostic_is_alternative_5oh_3oh_linear(manifest):
    candidates = generate_candidates_for_measurement(manifest, "LEU_UAA_WT_FULL")[4:]
    for candidate in candidates:
        assert candidate.five_prime_state is FivePrimeState.OH
        assert candidate.five_prime_state_status is AssumptionStatus.ALTERNATIVE
        assert candidate.three_prime_state is ThreePrimeState.OH
        assert candidate.three_prime_state_status is AssumptionStatus.ASSUMED
        assert candidate.topology is RnaTopology.LINEAR
        assert candidate.topology_status is AssumptionStatus.ASSUMED
        assert candidate.terminal_state_confirmed is False


@pytest.mark.parametrize("measurement_id", ["GLU_UUC_WT_FULL", "LEU_UAA_WT_FULL", "LEU_UAG_WT_FULL"])
def test_primary_precedes_secondary(manifest, measurement_id):
    candidates = generate_candidates_for_measurement(manifest, measurement_id)
    priorities = [item.candidate_priority for item in candidates]
    assert priorities == sorted(priorities)
    assert set(priorities) == {1, 2}


def test_candidate_ids_are_readable_exact_and_deterministic(manifest):
    candidates = generate_candidates_for_measurement(manifest, "GLU_UUC_WT_FULL")
    assert [item.candidate_id for item in candidates] == [
        "TRNA_GLU_UUC__CCA__5P__3OH__LINEAR",
        "TRNA_GLU_UUC__CC__5P__3OH__LINEAR",
        "TRNA_GLU_UUC__CCA__5OH__3OH__LINEAR",
        "TRNA_GLU_UUC__CC__5OH__3OH__LINEAR",
    ]
    assert all("/" not in item.candidate_id and "\\" not in item.candidate_id for item in candidates)
    assert candidates == generate_candidates_for_measurement(manifest, "GLU_UUC_WT_FULL")


@pytest.mark.parametrize("measurement_id", ["GLU_UUC_WT_FULL", "LEU_UAA_WT_FULL", "LEU_UAG_WT_FULL"])
def test_candidate_category_formula_and_mass_type(manifest, measurement_id):
    for candidate in generate_candidates_for_measurement(manifest, measurement_id):
        assert candidate.candidate_category is CandidateCategory.UNMODIFIED_INTACT_RNA_CCA_TERMINAL_STATE
        assert candidate.theoretical_formula.startswith("C")
        assert "P" in candidate.theoretical_formula
        assert candidate.theoretical_mass > 0
        assert candidate.theoretical_mass_type == "MONOISOTOPIC_NEUTRAL"


@pytest.mark.parametrize("measurement_id", ["GLU_UUC_WT_FULL", "LEU_UAA_WT_FULL", "LEU_UAG_WT_FULL"])
def test_mass_uses_complete_sequence_and_matches_direct_calculation(manifest, measurement_id):
    for candidate in generate_candidates_for_measurement(manifest, measurement_id):
        parameters = IntactRnaMassParameters(
            five_prime_state=candidate.five_prime_state,
            three_prime_state=candidate.three_prime_state,
            topology=candidate.topology,
            cca_policy=CcaPolicy.AS_PROVIDED,
            convert_t_to_u=False,
            terminal_state_confirmed=False,
        )
        direct = calculate_intact_rna_mass(candidate.complete_sequence, parameters=parameters)
        assert candidate.theoretical_mass == direct.monoisotopic_neutral_mass
        assert candidate.theoretical_formula == direct.formula
        assert candidate.complete_sequence_sha256 == direct.sequence_sha256


def test_5p_minus_5oh_mass_delta_is_explicit_phosphate_delta(manifest):
    candidates = generate_candidates_for_measurement(manifest, "GLU_UUC_WT_FULL")
    primary = candidates[0]
    secondary = candidates[2]
    assert primary.complete_sequence == secondary.complete_sequence
    assert primary.theoretical_mass - secondary.theoretical_mass == pytest.approx(
        79.966330889, abs=1e-9
    )


@pytest.mark.parametrize("measurement_id", ["GLU_UUC_WT_FULL", "LEU_UAA_WT_FULL", "LEU_UAG_WT_FULL"])
def test_no_charge_proton_mz_or_average_mass_adjustment(manifest, measurement_id):
    for candidate in generate_candidates_for_measurement(manifest, measurement_id):
        assert candidate.charge_correction_applied is False
        assert candidate.proton_correction_applied is False
        assert candidate.mz_conversion_applied is False
        assert candidate.average_mass_used is False


@pytest.mark.parametrize(
    "measurement_id",
    ["LEU_UAA_WT_T1", "LEU_UAG_WT_T1", "GLU_UUC_WT_P1_AP"],
)
def test_digest_measurements_are_rejected(manifest, measurement_id):
    with pytest.raises(CandidateGenerationError, match="not FULL_LENGTH"):
        generate_candidates_for_measurement(manifest, measurement_id)


def test_unknown_measurement_is_rejected(manifest):
    with pytest.raises(KeyError, match="unknown measurement_id"):
        generate_candidates_for_measurement(manifest, "UNKNOWN_MEASUREMENT")


def test_unknown_sequence_is_rejected_without_mutating_manifest(manifest):
    identity = get_rna_identity(manifest, "TRNA_LEU_UAA")
    changed_identity = replace(identity, sequence=None, sequence_status=SequenceStatus.UNKNOWN)
    changed_manifest = replace(
        manifest,
        rna_identities=tuple(
            changed_identity if item.rna_identity_id == identity.rna_identity_id else item
            for item in manifest.rna_identities
        ),
    )
    before = asdict(changed_manifest)
    with pytest.raises(CandidateGenerationError, match="no candidate-ready sequence"):
        generate_candidates_for_measurement(changed_manifest, "LEU_UAA_WT_FULL")
    assert asdict(changed_manifest) == before


def test_unknown_registered_cca_mode_is_rejected(manifest):
    identity = get_rna_identity(manifest, "TRNA_LEU_UAG")
    changed_identity = replace(
        identity,
        registered_sequence_cca_mode=RegisteredSequenceCCAMode.UNKNOWN,
        registered_cca_tail_state=None,
    )
    changed_manifest = replace(
        manifest,
        rna_identities=tuple(
            changed_identity if item.rna_identity_id == identity.rna_identity_id else item
            for item in manifest.rna_identities
        ),
    )
    with pytest.raises(CandidateGenerationError, match="UNKNOWN registered CCA mode"):
        generate_candidates_for_measurement(changed_manifest, "LEU_UAG_WT_FULL")


def test_filename_is_not_used_for_identity_routing(manifest):
    measurement = get_measurement(manifest, "LEU_UAA_WT_FULL")
    renamed = replace(measurement, source_file_name="looks-like-Glu-UUC.mzML")
    changed_manifest = replace(
        manifest,
        measurements=tuple(
            renamed if item.measurement_id == measurement.measurement_id else item
            for item in manifest.measurements
        ),
    )
    candidates = generate_candidates_for_measurement(changed_manifest, measurement.measurement_id)
    assert all(item.rna_identity_id == "TRNA_LEU_UAA" for item in candidates)
    assert all(item.registered_sequence_sha256.startswith("71664e80") for item in candidates)


@pytest.mark.parametrize(
    ("measurement_id", "expected_identity", "forbidden_identities"),
    [
        ("LEU_UAA_WT_FULL", "TRNA_LEU_UAA", {"TRNA_LEU_UAG", "TRNA_GLU_UUC"}),
        ("LEU_UAG_WT_FULL", "TRNA_LEU_UAG", {"TRNA_LEU_UAA", "TRNA_GLU_UUC"}),
        ("GLU_UUC_WT_FULL", "TRNA_GLU_UUC", {"TRNA_LEU_UAA", "TRNA_LEU_UAG"}),
    ],
)
def test_no_cross_identity_candidates(manifest, measurement_id, expected_identity, forbidden_identities):
    candidates = generate_candidates_for_measurement(manifest, measurement_id)
    assert {item.rna_identity_id for item in candidates} == {expected_identity}
    assert not ({item.rna_identity_id for item in candidates} & forbidden_identities)


@pytest.mark.parametrize("measurement_id", ["GLU_UUC_WT_FULL", "LEU_UAA_WT_FULL", "LEU_UAG_WT_FULL"])
def test_mass_equivalent_groups_are_formula_based_and_standalone(manifest, measurement_id):
    candidates = generate_candidates_for_measurement(manifest, measurement_id)
    for candidate in candidates:
        assert candidate.mass_equivalent_group_id.startswith("MEG__")
        assert candidate.candidate_ambiguity_count == 1
        assert candidate.mass_equivalent_candidate_ids == (candidate.candidate_id,)
    assert len({item.mass_equivalent_group_id for item in candidates}) == len(candidates)
    repeated = generate_candidates_for_measurement(manifest, measurement_id)
    assert [item.mass_equivalent_group_id for item in candidates] == [
        item.mass_equivalent_group_id for item in repeated
    ]


@pytest.mark.parametrize("measurement_id", ["GLU_UUC_WT_FULL", "LEU_UAA_WT_FULL", "LEU_UAG_WT_FULL"])
def test_assumption_count_is_four_by_documented_rule(manifest, measurement_id):
    candidates = generate_candidates_for_measurement(manifest, measurement_id)
    assert all(item.candidate_assumption_count == 4 for item in candidates)


@pytest.mark.parametrize("measurement_id", ["GLU_UUC_WT_FULL", "LEU_UAA_WT_FULL", "LEU_UAG_WT_FULL"])
def test_false_certainty_safeguards(manifest, measurement_id):
    for candidate in generate_candidates_for_measurement(manifest, measurement_id):
        assert candidate.mass_match_only is True
        assert candidate.unmodified_candidate is True
        assert candidate.cca_state_confirmed is False
        assert candidate.terminal_state_confirmed is False
        assert candidate.structure_identity_assigned is False
        assert candidate.position_assigned is False
        assert candidate.modification_assigned is False
        assert candidate.biological_cause_assigned is False
        assert candidate.rnase_t_assigned is False


def test_all_formal_flags_are_false(manifest):
    assert CANDIDATE_GENERATION_APPLIED_TO_FORMAL_SCORE is False
    assert CANDIDATE_GENERATION_APPLIED_TO_RANKING is False
    assert CANDIDATE_GENERATION_APPLIED_TO_CANDIDATE_FILTERING is False
    assert CANDIDATE_GENERATION_APPLIED_TO_FINAL_CONSENSUS is False
    for candidate in generate_candidates_for_measurement(manifest, "GLU_UUC_WT_FULL"):
        assert candidate.applied_to_formal_score is False
        assert candidate.applied_to_ranking is False
        assert candidate.applied_to_candidate_filtering is False
        assert candidate.applied_to_final_consensus is False


def test_manifest_identity_variant_and_terminal_set_are_not_mutated(manifest):
    identity = get_rna_identity(manifest, "TRNA_LEU_UAA")
    variant = generate_cca_tail_variants(
        identity.sequence,
        identity.registered_sequence_cca_mode,
    )[0]
    manifest_before = asdict(manifest)
    identity_before = asdict(identity)
    variant_before = asdict(variant)
    terminal_before = asdict(PRIMARY_DEFAULT)
    build_unmodified_intact_candidate(
        rna_identity_id=identity.rna_identity_id,
        sample_id="LEU_UAA_WT",
        measurement_id="LEU_UAA_WT_FULL",
        cca_variant=variant,
        terminal_candidate_set=PRIMARY_DEFAULT,
    )
    assert asdict(manifest) == manifest_before
    assert asdict(identity) == identity_before
    assert asdict(variant) == variant_before
    assert asdict(PRIMARY_DEFAULT) == terminal_before


def test_candidate_result_is_immutable(manifest):
    candidate = generate_candidates_for_measurement(manifest, "GLU_UUC_WT_FULL")[0]
    with pytest.raises(FrozenInstanceError):
        candidate.candidate_priority = 99


def test_invalid_secondary_flag_is_rejected(manifest):
    with pytest.raises(TypeError, match="must be boolean"):
        generate_candidates_for_measurement(
            manifest,
            "GLU_UUC_WT_FULL",
            include_secondary_terminal_state=1,
        )


def test_only_allowed_default_terminal_states_are_generated(manifest):
    candidates = generate_candidates_for_measurement(manifest, "LEU_UAA_WT_FULL")
    assert {item.five_prime_state for item in candidates} == {
        FivePrimeState.MONOPHOSPHATE,
        FivePrimeState.OH,
    }
    assert {item.three_prime_state for item in candidates} == {ThreePrimeState.OH}
    assert {item.topology for item in candidates} == {RnaTopology.LINEAR}
    assert FivePrimeState.DIPHOSPHATE not in {item.five_prime_state for item in candidates}
    assert FivePrimeState.TRIPHOSPHATE not in {item.five_prime_state for item in candidates}

@pytest.mark.parametrize(
    ("measurement_id", "expected_count"),
    [("GLU_UUC_WT_FULL", 4), ("LEU_UAA_WT_FULL", 8), ("LEU_UAG_WT_FULL", 8)],
)
def test_candidates_keep_four_mass_references_and_baseline_safeguards(
    manifest, measurement_id, expected_count
):
    candidates = generate_candidates_for_measurement(manifest, measurement_id)
    assert len(candidates) == expected_count
    for candidate in candidates:
        assert candidate.theoretical_monoisotopic_neutral_mass == candidate.theoretical_mass
        assert candidate.theoretical_average_neutral_molecular_mass_m > 0
        assert candidate.theoretical_average_m_plus_h > candidate.theoretical_average_neutral_molecular_mass_m
        assert candidate.theoretical_average_m_minus_h < candidate.theoretical_average_neutral_molecular_mass_m
        assert candidate.primary_reference_candidate.value == "AVERAGE_NEUTRAL_M"
        assert candidate.observed_output_species == "UNKNOWN"
        assert candidate.observed_output_species_confirmed is False
        assert candidate.candidate_role.value == "UNMODIFIED_REFERENCE_BASELINE"
        assert candidate.native_modifications_expected is True
        assert candidate.modification_mass_not_yet_applied is True
        assert candidate.biological_unmodified_state_assigned is False
        assert candidate.target_rna_identity_confirmed_by_mass is False
        assert candidate.co_captured_rna_excluded is False

def test_generate_candidates_for_measurement_respects_cca_tail_config(manifest):
    cca_tail_config = {
        "enabled": True,
        "excludes_cca_candidate_states": ["CCA"],
    }
    candidates = generate_candidates_for_measurement(
        manifest, "LEU_UAA_WT_FULL", cca_tail_config=cca_tail_config,
    )
    assert candidate_states(candidates) == [CCATailState.CCA, CCATailState.CCA]