from pathlib import Path

import pytest

from rna_masshunter.masses import PROTON_MASS
from rna_masshunter.sciex_t1_fragment_state_series_audit import StateLabel, build_default_state_delta_definitions
from rna_masshunter.sciex_t1_replicate_consistency_audit import ReplicateRunPeak
from rna_masshunter.sciex_p1ap_nucleoside_state_audit import (
    MatchStatus, NucleosideCandidate, NucleosideCandidateClass,
    NucleosideIonHypothesis, P1APAuditParameters, RTCoelutionStatus,
    StateFamilyQuality, TargetRTSummary, _state_candidates, _target_id,
    build_nucleoside_state_families, generate_negative_nucleoside_ion_hypotheses,
    generate_nucleoside_candidates, generate_positive_nucleoside_ion_hypotheses,
    match_nucleoside_ions_to_peaks, reconcile_p1ap_with_t1_and_full_length,
)

ROOT = Path(__file__).parents[1]


def candidate(identifier="A", name="adenosine", base="A", mass=267.096753):
    return NucleosideCandidate(candidate_id=identifier, candidate_name=name, parent_base=base,
        candidate_class=NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE,
        modification_components=(), molecular_formula="C10H13N5O4",
        theoretical_neutral_mass=mass, mass_provenance="TEST", source_registry="TEST",
        structure_constraint_status="CANONICAL_FORMULA", candidate_block_reasons=(),
        eligible_for_mass_matching=True, source_candidate={})


def ion(identifier="A", name="adenosine", base="A", mz=268.104029, adduct="[M+H]+"):
    return NucleosideIonHypothesis(ion_hypothesis_id=f"I_{identifier}_{adduct}", candidate_id=identifier,
        candidate_name=name, parent_base=base, candidate_class=NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE,
        ion_mode="POSITIVE", charge=1, adduct_type=adduct, adduct_mass=PROTON_MASS,
        theoretical_mz=mz, ion_hypothesis_status="ELIGIBLE_POSITIVE_SOURCE",
        ion_hypothesis_block_reasons=("SOURCE_POLARITY_POSITIVE",))


def peak(identifier, mz, *, centroid=None, recurrence=0.2, prominence=0.1, fwhm=0.02,
         rank=1, quality="MAJOR_SHARP"):
    return ReplicateRunPeak(run_label="P1", peak_id=identifier, apex_mz=mz,
        centroid_mz=mz if centroid is None else centroid, raw_apex_intensity=100,
        normalized_apex_intensity=1 / rank, raw_integrated_intensity=10,
        normalized_integrated_intensity=0.5 / rank, relative_intensity=1 / rank,
        intensity_rank=rank, prominence=prominence, relative_prominence=prominence,
        fwhm=fwhm, left_bound_mz=mz - 0.02, right_bound_mz=mz + 0.02,
        supporting_ms1_scan_count=round(100 * recurrence), total_ms1_scan_count=100,
        scan_recurrence_fraction=recurrence, first_supporting_scan_time=1,
        last_supporting_scan_time=2, detection_status=quality, detection_block_reasons=())


def trace(target_id, apex=5.0, start=4.9, end=5.1, recurrence=0.2):
    return TargetRTSummary(target_id, start, end, apex, apex, round(recurrence * 100), 100,
        recurrence, end - start, "TARGETED_SCAN_LEVEL_RT_EVIDENCE")


def state_rows(labels, *, charge_mz=268.0, rt_apices=None, recurrence=0.2, centroid_shift=0,
               quality="MAJOR_SHARP", ions=None):
    definitions = build_default_state_delta_definitions()
    chosen = [x for x in definitions if x.state_label in labels]
    ion_rows = ions or [ion(mz=charge_mz)]
    peaks = [peak(f"P{i}", charge_mz + definition.target_neutral_delta,
                  centroid=charge_mz + definition.target_neutral_delta + centroid_shift,
                  recurrence=recurrence, rank=i + 1, quality=quality)
             for i, definition in enumerate(chosen)]
    rt = {}
    for ion_row in ion_rows:
        for i, definition in enumerate(chosen):
            apex = (rt_apices or {}).get(definition.state_label, 5.0)
            rt[_target_id(ion_row.ion_hypothesis_id, definition.state_label)] = trace(
                _target_id(ion_row.ion_hypothesis_id, definition.state_label), apex=apex,
                start=apex - 0.1, end=apex + 0.1, recurrence=recurrence)
    return _state_candidates(ion_rows, peaks, definitions, rt, P1APAuditParameters())


def family(labels, **kwargs):
    rows = state_rows(labels, **kwargs)
    return build_nucleoside_state_families(state_candidates=rows)[0]


def test_canonical_a_g_c_u_candidates_use_existing_formula_registry():
    rows = generate_nucleoside_candidates(project_root=ROOT, sequence="AGCU", modification_registry=[])
    canonical = {x.parent_base: x for x in rows if x.candidate_class is NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE}
    assert set(canonical) == set("AGCU")
    assert all(x.theoretical_neutral_mass and x.theoretical_neutral_mass > 240 for x in canonical.values())
    assert all(x.molecular_formula != "MODEL_NOT_DEFINED" for x in canonical.values())


def test_residual_phosphate_candidates_are_diagnostic_not_ion_targets():
    rows = generate_nucleoside_candidates(project_root=ROOT, sequence="AGCU", modification_registry=[])
    residual = [x for x in rows if x.candidate_class is NucleosideCandidateClass.MONOPHOSPHATE_RESIDUAL]
    assert residual and all(not x.eligible_for_mass_matching for x in residual)
    generated = generate_positive_nucleoside_ion_hypotheses(rows)
    assert all(x.candidate_class is not NucleosideCandidateClass.MONOPHOSPHATE_RESIDUAL for x in generated)


def test_positive_protonated_mz_and_negative_block():
    source = candidate(mass=267.096753)
    rows = generate_positive_nucleoside_ion_hypotheses([source])
    assert len(rows) == 1 and rows[0].adduct_type == "[M+H]+"
    assert rows[0].theoretical_mz == pytest.approx(267.096753 + PROTON_MASS)
    assert generate_negative_nucleoside_ion_hypotheses([source]) == ()


def test_optional_adducts_only_when_configured():
    source = candidate()
    default = generate_positive_nucleoside_ion_hypotheses([source])
    configured = generate_positive_nucleoside_ion_hypotheses([source], ion_config={"adducts": ("H", "Na", "K")})
    assert {x.adduct_type for x in default} == {"[M+H]+"}
    assert {x.adduct_type for x in configured} == {"[M+H]+", "[M+Na]+", "[M+K]+"}
    with pytest.raises(ValueError): generate_positive_nucleoside_ion_hypotheses([source], ion_config={"adducts": ("NH4",)})


@pytest.mark.parametrize("offset,status", [(0, MatchStatus.STRICT), (0.002, MatchStatus.SUPPORTIVE)])
def test_exact_and_tolerance_candidate_match(offset, status):
    # At m/z 268, 5 ppm=0.00134 and 10 ppm=0.00268.
    rows = match_nucleoside_ions_to_peaks([ion(mz=268)], [peak("P", 268 + offset)])
    assert len(rows) == 1 and rows[0].apex_match_status is status


def test_outside_tolerance_is_not_matched():
    assert match_nucleoside_ions_to_peaks([ion(mz=268)], [peak("P", 268.003)]) == ()


def test_broad_rt_trace_is_reported_as_possible_background_match():
    ion_row = ion(mz=268)
    target = _target_id(ion_row.ion_hypothesis_id, StateLabel.BASE_STATE)
    evidence = {target: TargetRTSummary(target, 1.0, 20.0, 5.0, 8.0, 50, 100, 0.5, 19.0, "TARGETED_SCAN_LEVEL_RT_EVIDENCE")}
    row = match_nucleoside_ions_to_peaks([ion_row], [peak("P", 268)], rt_evidence=evidence)[0]
    assert row.match_quality_status == "POSSIBLE_BACKGROUND_MATCH"
    assert "POSSIBLE_BACKGROUND_EXPLANATION" in row.match_block_reasons


def test_identity_ambiguity_is_preserved():
    ions = [ion("A", mz=268), ion("B", name="isobar", base="G", mz=268)]
    rows = match_nucleoside_ions_to_peaks(ions, [peak("P", 268)])
    assert len(rows) == 2 and all(x.identity_ambiguity_status == "IDENTITY_AMBIGUOUS" for x in rows)


def test_adduct_ambiguity_is_preserved():
    ions = [ion("A", mz=268, adduct="[M+H]+"), ion("A", mz=268, adduct="[M+Na]+")]
    rows = match_nucleoside_ions_to_peaks(ions, [peak("P", 268)])
    assert len(rows) == 2 and all(x.adduct_ambiguity_status == "ADDUCT_AMBIGUOUS" for x in rows)


@pytest.mark.parametrize("labels,expected", [
    ({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT}, {StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT}),
    ({StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT}, {StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT}),
    ({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT, StateLabel.PLUS_32_EQUIVALENT}, {StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT, StateLabel.PLUS_32_EQUIVALENT}),
    ({StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT, StateLabel.PLUS_34_EQUIVALENT}, {StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT, StateLabel.PLUS_34_EQUIVALENT}),
])
def test_state_family_patterns(labels, expected):
    assert set(family(labels).detected_state_labels) == expected


def test_missing_intermediate_is_incomplete():
    row = family({StateLabel.BASE_STATE, StateLabel.PLUS_32_EQUIVALENT})
    assert "INCOMPLETE_STATE_SERIES" in row.series_block_reasons


def test_isotope_alternative_downgrades_confidence():
    row = family({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT}, quality="ISOTOPE_OR_ENVELOPE_COMPONENT")
    assert row.isotope_compatibility == "POSSIBLE"
    assert row.series_quality_status is StateFamilyQuality.AMBIGUOUS_STATE_FAMILY


def test_rt_coelution_and_separation_are_evidence_not_identity():
    coeluting = family({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT}, rt_apices={StateLabel.BASE_STATE: 5.0, StateLabel.PLUS_16_EQUIVALENT: 5.03})
    distinct = family({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT}, rt_apices={StateLabel.BASE_STATE: 5.0, StateLabel.PLUS_16_EQUIVALENT: 8.0})
    assert coeluting.rt_coelution_status is RTCoelutionStatus.COELUTING
    assert distinct.rt_coelution_status is RTCoelutionStatus.DISTINCT_RETENTION
    assert not coeluting.chemical_identity_assigned and not distinct.chemical_identity_assigned


def test_low_recurrence_and_apex_centroid_disagreement_lower_quality():
    low = family({StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT}, recurrence=0.001)
    disagreement = family({StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT}, centroid_shift=0.02)
    assert low.series_quality_status is StateFamilyQuality.LOW_QUALITY_STATE_FAMILY
    assert "LOW_SCAN_RECURRENCE" in low.series_block_reasons
    assert "APEX_CENTROID_DISAGREEMENT" in disagreement.series_block_reasons
    assert disagreement.series_quality_status is not StateFamilyQuality.HIGH_QUALITY_STATE_FAMILY


def test_t1_reconciliation_reports_p1_evidence_without_t1_series():
    row = family({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT})
    recon = reconcile_p1ap_with_t1_and_full_length([row], t1_result="T1_SERIES_NOT_OBSERVED", full_length_series=[0, 18, 34, 50])[0]
    assert recon.p1ap_state_status == "P1AP_STATE_EVIDENCE_WITHOUT_T1_SERIES"
    assert recon.localization_status == "NUCLEOSIDE_CLASS_ONLY_T1_FRAGMENT_LOCALIZATION_UNSUPPORTED"


def test_full_length_reconciliation_compatible_partial_not_observed_and_insufficient():
    compatible = family({StateLabel.BASE_STATE, StateLabel.PLUS_18_EQUIVALENT, StateLabel.PLUS_34_EQUIVALENT})
    partial = family({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT, StateLabel.PLUS_18_EQUIVALENT})
    assert reconcile_p1ap_with_t1_and_full_length([compatible], full_length_series=[0, 18, 34, 50])[0].full_length_reconciliation_status == "FULL_LENGTH_DELTA_PATTERN_COMPATIBLE"
    assert reconcile_p1ap_with_t1_and_full_length([partial], full_length_series=[0, 18, 34, 50])[0].full_length_reconciliation_status == "PARTIALLY_COMPATIBLE_WITH_FULL_LENGTH_PATTERN"
    assert reconcile_p1ap_with_t1_and_full_length([], full_length_series=[0, 18, 34, 50])[0].full_length_reconciliation_status == "P1AP_STATE_PATTERN_NOT_OBSERVED"
    assert reconcile_p1ap_with_t1_and_full_length([], full_length_series=None)[0].full_length_reconciliation_status == "INSUFFICIENT_P1AP_EVIDENCE"


def test_structural_incompatibility_is_blocked_from_ion_generation():
    rows = generate_nucleoside_candidates(project_root=ROOT, sequence="A", modification_registry=[],
        structure_constraints=[{"candidate_id": "BAD", "candidate_name": "bad", "parent_base": "A",
                                "components": ["m1A", "m6A"], "status": "STRUCTURALLY_INCOMPATIBLE_COMBINATION"}])
    blocked = next(x for x in rows if x.candidate_id == "BAD")
    assert not blocked.eligible_for_mass_matching
    assert "STRUCTURALLY_INCOMPATIBLE_COMBINATION" in blocked.candidate_block_reasons
    assert all(x.candidate_id != "BAD" for x in generate_positive_nucleoside_ion_hypotheses(rows))


def test_localization_and_formal_safeguards_remain_false():
    row = family({StateLabel.BASE_STATE, StateLabel.PLUS_16_EQUIVALENT})
    assert not row.exact_nucleotide_localization and not row.exact_atom_localization
    assert not row.modification_assigned and not row.reaction_order_assigned
    assert not row.formal_propagation and not row.applied_to_formal_score
    assert not row.applied_to_ranking and not row.applied_to_candidate_filtering and not row.applied_to_final_consensus


def test_determinism_for_reversed_candidate_ion_and_peak_order():
    definitions = build_default_state_delta_definitions()[:2]
    ions = [ion("B", name="B", mz=268), ion("A", name="A", mz=268)]
    peaks = [peak("P0", 268), peak("P1", 268 + definitions[1].target_neutral_delta)]
    rt = {_target_id(x.ion_hypothesis_id, d.state_label): trace(_target_id(x.ion_hypothesis_id, d.state_label)) for x in ions for d in definitions}
    a = _state_candidates(ions, peaks, definitions, rt, P1APAuditParameters())
    b = _state_candidates(list(reversed(ions)), list(reversed(peaks)), definitions, rt, P1APAuditParameters())
    assert a == b
