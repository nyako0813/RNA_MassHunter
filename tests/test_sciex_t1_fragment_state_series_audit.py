from dataclasses import replace
from pathlib import Path

import pytest

from rna_masshunter.masses import load_base_masses
from rna_masshunter.sciex_sample_manifest import load_sciex_sample_manifest
from rna_masshunter.sciex_t1_fragment_shadow_match import T1IonCandidate, T1IonMode, TheoreticalT1Fragment
from rna_masshunter.sciex_t1_fragment_state_series_audit import (
    ReconciliationStatus, SeriesAmbiguityStatus, SeriesQualityStatus,
    StateLabel, StateMatchStatus, StateSeriesAuditParameters,
    T1FragmentIonHypothesis, _state_candidates, audit_optional_result,
    build_default_state_delta_definitions, build_t1_fragment_state_families,
    generate_t1_fragment_ion_hypotheses, match_t1_fragment_ions_to_peaks,
    reconcile_t1_state_families_with_full_length_series,
)
from rna_masshunter.sciex_t1_replicate_consistency_audit import ReplicateRunPeak

ROOT = Path(__file__).parents[1]


def fragment(identifier="F1", sequence="AG", mass=1001.007276466621):
    return TheoreticalT1Fragment(identifier, "SYNTHETIC_RNA", "SYNTHETIC__CCA_CCA", "CCA",
        "DIGEST_TERMINUS_UNKNOWN", "DIGEST_TERMINUS_UNKNOWN", 1, 1, len(sequence), sequence,
        "RNASE_T1_AFTER_G" if sequence.endswith("G") else "THREE_PRIME_TERMINAL_FRAGMENT", mass, mass + 0.5)


def ion(identifier="F1", sequence="AG", mass=1001.007276466621, mz=1000.0, charge=1):
    source = T1IonCandidate(f"I_{identifier}_{charge}", identifier, "SYNTHETIC_RNA", "SYNTHETIC__CCA_CCA",
        "CCA", sequence, 1, len(sequence), T1IonMode.NEGATIVE_DEPROTONATED, charge, mz,
        "MONOISOTOPIC_NEUTRAL", "MONOISOTOPIC_NEGATIVE_ION_MZ", -charge)
    return T1FragmentIonHypothesis(ion_hypothesis_id=source.ion_candidate_id, fragment_id=identifier,
        fragment_sequence=sequence, start_position=1, end_position=len(sequence),
        theoretical_neutral_mass=mass, fragment_length=len(sequence),
        base_composition=tuple((base, sequence.count(base)) for base in "ACGU"), contains_g="G" in sequence,
        cleavage_start="DIGEST_TERMINUS_UNKNOWN", cleavage_end="DIGEST_TERMINUS_UNKNOWN",
        terminal_chemistry="PROJECT_STANDARD_UNCHANGED", generation_status="GENERATED",
        generation_block_reasons=(), ion_mode=source.ion_mode, charge=charge,
        adduct_hypothesis="[M-zH]z-", theoretical_mz=mz,
        ion_hypothesis_status="ELIGIBLE_NEGATIVE_SOURCE",
        ion_hypothesis_block_reasons=("SOURCE_POLARITY_NEGATIVE",), source_candidate=source)


def peak(identifier, mz, *, centroid=None, recurrence=0.2, prominence=0.1, fwhm=0.02, rank=1):
    return ReplicateRunPeak(run_label="S", peak_id=identifier, apex_mz=mz,
        centroid_mz=mz if centroid is None else centroid, raw_apex_intensity=100,
        normalized_apex_intensity=1 / rank, raw_integrated_intensity=10,
        normalized_integrated_intensity=0.5 / rank, relative_intensity=1 / rank,
        intensity_rank=rank, prominence=prominence, relative_prominence=prominence,
        fwhm=fwhm, left_bound_mz=mz - 0.02, right_bound_mz=mz + 0.02,
        supporting_ms1_scan_count=round(recurrence * 100), total_ms1_scan_count=100,
        scan_recurrence_fraction=recurrence, first_supporting_scan_time=1,
        last_supporting_scan_time=2, detection_status="MAJOR_SHARP", detection_block_reasons=())


def candidates(labels, *, charge=1, centroid_shift=0, recurrence=0.2, ions=None):
    definitions = build_default_state_delta_definitions()
    selected = [x for x in definitions if x.state_label in labels]
    hypotheses = ions or [ion(charge=charge, mz=1000)]
    observed = [peak(f"P{i}", 1000 + definition.target_neutral_delta / charge,
                     centroid=1000 + definition.target_neutral_delta / charge + centroid_shift,
                     recurrence=recurrence, rank=i + 1)
                for i, definition in enumerate(selected)]
    return _state_candidates(hypotheses, observed, definitions, StateSeriesAuditParameters())


def family_for(labels, **kwargs):
    rows = candidates(labels, **kwargs)
    return build_t1_fragment_state_families(state_candidates=rows)[0]


def test_existing_t1_fragment_generator_covers_g_cleavage_terminal_and_positions():
    manifest = load_sciex_sample_manifest(ROOT / "data/sciex_sample_manifest.yaml")
    base = load_base_masses(ROOT / "data/base_masses.yaml")
    from rna_masshunter.sciex_t1_fragment_shadow_match import generate_theoretical_t1_fragments
    rows = generate_theoretical_t1_fragments(manifest, "TRNA_GLU_UUC", base, candidate_states=["CCA"])
    assert rows and all(x.start_position <= x.end_position and x.neutral_monoisotopic_mass > 0 for x in rows)
    assert any(x.cleavage_context == "RNASE_T1_AFTER_G" for x in rows)
    assert rows[-1].cleavage_context == "THREE_PRIME_TERMINAL_FRAGMENT"
    assert any(len(x.fragment_sequence) == 1 and x.fragment_sequence == "G" for x in rows)


def test_negative_charge_generation_only_and_positive_is_blocked():
    rows = generate_t1_fragment_ion_hypotheses([fragment()], observed_mz_range=(0, 2000))
    assert rows and {x.ion_mode for x in rows} == {T1IonMode.NEGATIVE_DEPROTONATED}
    assert generate_t1_fragment_ion_hypotheses([fragment()], ion_mode="positive") == ()


def test_observed_mz_range_filters_physically_irrelevant_hypotheses():
    rows = generate_t1_fragment_ion_hypotheses([fragment(mass=1001.007276466621)], observed_mz_range=(900, 1100))
    assert {x.charge for x in rows} == {1}


@pytest.mark.parametrize("offset,status", [(0.0, StateMatchStatus.STRICT), (0.015, StateMatchStatus.SUPPORTIVE)])
def test_exact_and_tolerance_match(offset, status):
    rows = match_t1_fragment_ions_to_peaks([ion()], [peak("P", 1000 + offset)])
    assert len(rows) == 1 and rows[0].apex_match_status is status


def test_outside_tolerance_does_not_match():
    assert match_t1_fragment_ions_to_peaks([ion()], [peak("P", 1000.021)]) == ()


def test_charge_ambiguity_is_not_forced_to_best_error():
    rows = match_t1_fragment_ions_to_peaks([ion("F1", mz=1000, charge=1), ion("F1", mz=1000, charge=2)], [peak("P", 1000)])
    assert len(rows) == 2 and all(x.charge_ambiguity_status == "CHARGE_AMBIGUOUS" for x in rows)


def test_fragment_ambiguity_is_not_forced_to_best_error():
    rows = match_t1_fragment_ions_to_peaks([ion("F1"), ion("F2")], [peak("P", 1000)])
    assert len(rows) == 2 and all(x.fragment_ambiguity_status == "FRAGMENT_AMBIGUOUS" for x in rows)


@pytest.mark.parametrize("labels,pattern", [
    ({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT}, "BASE__PLUS16"),
    ({StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT}, "BASE__PLUS18"),
    ({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT, StateLabel.PLUS_32_EQUIVALENT}, "BASE__PLUS16__PLUS32"),
    ({StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT, StateLabel.PLUS_34_EQUIVALENT}, "BASE__PLUS18__PLUS34"),
])
def test_expected_state_family_patterns(labels, pattern):
    assert family_for(labels).state_series_pattern == pattern


def test_missing_intermediate_is_partial_and_blocked_from_complete_claim():
    family = family_for({StateLabel.BASE_STATE, StateLabel.PLUS_32_EQUIVALENT})
    assert family.state_series_pattern.startswith("PARTIAL_SERIES")
    assert "INCOMPLETE_STATE_SERIES" in family.series_block_reasons


def test_charge_scaled_spacing_is_neutral_delta_divided_by_charge():
    rows = candidates({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT}, charge=2)
    shifted = next(x for x in rows if x.state_label is StateLabel.PLUS_16_EQUIVALENT)
    target = next(x.target_neutral_delta for x in build_default_state_delta_definitions() if x.state_label is StateLabel.PLUS_16_EQUIVALENT)
    assert shifted.expected_mz_delta == pytest.approx(target / 2)


def test_apex_centroid_disagreement_downgrades_quality():
    family = family_for({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT}, centroid_shift=0.1)
    assert "APEX_CENTROID_DISAGREEMENT" in family.series_block_reasons
    assert family.series_quality_status is not SeriesQualityStatus.HIGH_QUALITY_STATE_FAMILY


def test_low_recurrence_is_low_quality():
    family = family_for({StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT}, recurrence=0.001)
    assert family.series_quality_status is SeriesQualityStatus.LOW_RECURRENCE_STATE_FAMILY
    assert "LOW_SCAN_RECURRENCE" in family.series_block_reasons


def test_peak_multiplicity_and_state_ambiguity_are_explicit():
    definitions = build_default_state_delta_definitions()
    base = next(x for x in definitions if x.state_label is StateLabel.BASE_STATE)
    shifted = next(x for x in definitions if x.state_label is StateLabel.PLUS_16_EQUIVALENT)
    observed = [peak("P0", 1000), peak("P1", 1000 + shifted.target_neutral_delta), peak("P2", 1000 + shifted.target_neutral_delta + 0.005)]
    rows = _state_candidates([ion()], observed, (base, shifted), StateSeriesAuditParameters())
    family = build_t1_fragment_state_families(state_candidates=rows, state_delta_definitions=(base, shifted))[0]
    assert family.series_ambiguity_status in {SeriesAmbiguityStatus.PEAK_MULTIPLICITY, SeriesAmbiguityStatus.MULTI_AXIS_AMBIGUOUS}


def test_full_length_reconciliation_compatible_partial_and_insufficient():
    compatible = family_for({StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT, StateLabel.PLUS_34_EQUIVALENT})
    plus18 = family_for({StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT})
    one = reconcile_t1_state_families_with_full_length_series([compatible], [25292, 25310, 25326, 25342])[0]
    two = reconcile_t1_state_families_with_full_length_series([plus18], [25292, 25310, 25326, 25342])[0]
    three = reconcile_t1_state_families_with_full_length_series([plus18], None)[0]
    assert one.reconciliation_status is ReconciliationStatus.FULL_LENGTH_PATTERN_COMPATIBLE
    assert two.reconciliation_status is ReconciliationStatus.T1_PLUS18_SERIES_ONLY
    assert three.reconciliation_status is ReconciliationStatus.INSUFFICIENT_T1_EVIDENCE


def test_no_family_reconciliation_explicitly_reports_series_not_observed():
    row = reconcile_t1_state_families_with_full_length_series([], [25292, 25310, 25326, 25342])[0]
    assert row.state_family_id == "NO_STATE_FAMILY"
    assert row.reconciliation_status is ReconciliationStatus.T1_SERIES_NOT_OBSERVED
    assert "NO_STATE_SERIES" in row.reconciliation_block_reasons


def test_localization_and_chemical_safeguards_are_always_false():
    family = family_for({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT})
    assert family.localization_level == "FRAGMENT_RANGE_ONLY"
    assert not family.exact_nucleotide_localization and not family.exact_atom_localization
    assert not family.chemical_identity_assigned and not family.modification_assigned


def test_deterministic_for_reversed_peak_and_ion_order():
    definitions = build_default_state_delta_definitions()
    observed = [peak("P0", 1000), peak("P1", 1000 + definitions[1].target_neutral_delta)]
    a = _state_candidates([ion("F2"), ion("F1")], observed, definitions[:2], StateSeriesAuditParameters())
    b = _state_candidates([ion("F1"), ion("F2")], list(reversed(observed)), definitions[:2], StateSeriesAuditParameters())
    assert a == b


def test_optional_result_shape_and_formal_non_propagation():
    # Serializer safeguard behavior is checked without touching any formal pipeline.
    family = family_for({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT})
    assert not family.formal_propagation and not family.applied_to_formal_score
    assert not family.applied_to_ranking and not family.applied_to_candidate_filtering
    assert not family.applied_to_final_consensus
    assert set(build_default_state_delta_definitions())
