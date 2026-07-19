from dataclasses import replace
from pathlib import Path

import pytest

from rna_masshunter.modifications import load_modifications
from rna_masshunter.sciex_intact_oxygen_water_state_audit import (
    PeakShapeSupportClass,
    StateEdgeType,
    StateRelationClass,
    StateSeriesPattern,
    audit_oxygen_water_state_series,
    oxygen_water_reference_provenance,
)
from rna_masshunter.sciex_intact_peak_family import (
    DeltaMassDefinition,
    PeakFamilyParameters,
    PeakFamilyPeak,
    PeakQualityClass,
    SciexIntactPeakFamilyResult,
    build_delta_reference_registry,
    generate_delta_pairs,
    match_delta_pairs,
)
from rna_masshunter.sciex_reconstructed_profile_registry import (
    DEFAULT_PROFILE_REGISTRY,
    analyze_loaded_profile_peak_families,
    load_profile_source,
    resolve_profile_source,
    route_profile_source_to_candidates,
)
from rna_masshunter.sciex_sample_manifest import load_sciex_sample_manifest

ROOT = Path(__file__).parent
GLU_RUNTIME = ROOT / ".cache/runtime_profiles/LC-MS_旧WT.xlsx"


def peak(index, mass, *, shoulder=False, duplicate=False, relative=0.5, fwhm=3.0, prominence=100.0):
    return PeakFamilyPeak(
        peak_id=f"P{index}", source_id="GLU", measurement_id="M", rna_identity="TRNA_GLU_UUC",
        apex_mass=mass, centroid_mass=mass + 0.1, apex_intensity=relative * 1000,
        integrated_intensity=relative * 2000, relative_apex_intensity=relative,
        relative_integrated_intensity=relative, left_boundary_mass=mass - 2,
        right_boundary_mass=mass + 2, peak_width_da=4.0, fwhm_da=fwhm,
        prominence=prominence, relative_prominence=0.1, sharpness_score=10.0,
        nearest_peak_separation_da=None, peak_overlap_fraction=0.0,
        peak_detection_status="DETECTED", peak_quality_class=PeakQualityClass.MAJOR_SHARP,
        selected_as_major_peak=True, possible_isotope_or_reconstruction_artifact=False,
        possible_shoulder=shoulder, possible_duplicate_peak=duplicate, possible_adduct=False,
        possible_output_convention_offset=False,
    )


def result(masses, *, shoulder_index=None, duplicate_index=None):
    peaks = tuple(peak(i, mass, shoulder=i == shoulder_index, duplicate=i == duplicate_index,
                       relative=(i + 1) / len(masses)) for i, mass in enumerate(masses))
    params = PeakFamilyParameters()
    pairs = generate_delta_pairs(peaks)
    refs = build_delta_reference_registry()
    matches = match_delta_pairs(pairs, refs, parameters=params)
    return SciexIntactPeakFamilyResult(
        "COMPLETED", "TEST", params, peaks, peaks, pairs, refs, matches, (), (),
    )


def test_reference_provenance_reuses_average_and_mono_references():
    refs = oxygen_water_reference_provenance(build_delta_reference_registry())
    by_id = {ref.reference_id: ref for ref in refs}
    assert by_id["O_AVERAGE_DELTA"].reference_delta_da == pytest.approx(15.9994)
    assert by_id["WATER_AVERAGE_DELTA_HYDRATION"].reference_delta_da == pytest.approx(18.01528)
    assert by_id["O_AVERAGE_DELTA"].delta_mass_definition is DeltaMassDefinition.AVERAGE_DELTA
    assert by_id["WATER_AVERAGE_DELTA_HYDRATION"].mass_definition_compatible is True
    assert by_id["O_MONOISOTOPIC_DELTA"].mass_definition_compatible is False
    assert by_id["WATER_MONOISOTOPIC_DELTA_HYDRATION"].mass_definition_compatible is False


@pytest.mark.parametrize("delta, expected", [
    (16.0, StateRelationClass.O_EQUIVALENT_STRICT),
    (16.8, StateRelationClass.O_EQUIVALENT_EXPLORATORY),
    (18.0, StateRelationClass.H2O_EQUIVALENT_STRICT),
])
def test_o_and_water_relation_classes(delta, expected):
    audit = audit_oxygen_water_state_series(result((100.0, 100.0 + delta)))
    relation = next(item for item in audit.relations if item.state_relation_class is not StateRelationClass.NO_MATCH)
    assert relation.state_relation_class is expected
    assert relation.observed_delta_mass_definition is DeltaMassDefinition.AVERAGE_DELTA
    assert relation.mass_definition_compatible is True


def test_no_match_and_isolated_peaks():
    audit = audit_oxygen_water_state_series(result((100.0, 1500.0)))
    assert audit.relations[0].state_relation_class is StateRelationClass.NO_MATCH
    assert not audit.edges
    assert [series.member_count for series in audit.series] == [1, 1]
    assert all(series.series_pattern is StateSeriesPattern.UNRESOLVED_STATE_SERIES for series in audit.series)


def test_edges_are_directed_lower_to_higher_without_reverse_duplicate():
    audit = audit_oxygen_water_state_series(result((118.0, 100.0)))
    assert len(audit.edges) == 1
    edge = audit.edges[0]
    assert edge.lower_apex_mass == 100.0 and edge.higher_apex_mass == 118.0
    assert edge.edge_type is StateEdgeType.PLUS_H2O_EQUIVALENT
    assert edge.reaction_direction_assigned is False
    assert edge.precursor_product_assigned is False
    assert edge.oxidation_direction_assigned is False


@pytest.mark.parametrize("masses,pattern,o_count,water_count", [
    ((100.0, 116.0), StateSeriesPattern.SINGLE_O_STEP, 1, 0),
    ((100.0, 116.0, 132.0), StateSeriesPattern.MULTIPLE_SEQUENTIAL_O_STEPS, 2, 0),
    ((100.0, 118.0), StateSeriesPattern.SINGLE_H2O_STEP, 0, 1),
    ((100.0, 118.0, 136.0), StateSeriesPattern.MULTIPLE_SEQUENTIAL_H2O_STEPS, 0, 2),
    ((100.0, 118.0, 134.0, 150.0), StateSeriesPattern.MIXED_H2O_AND_O_STEPS, 2, 1),
])
def test_connected_state_series_patterns(masses, pattern, o_count, water_count):
    series = max(audit_oxygen_water_state_series(result(masses)).series, key=lambda item: item.member_count)
    assert series.member_count == len(masses)
    assert series.series_pattern is pattern
    assert series.o_equivalent_edge_count == o_count
    assert series.h2o_equivalent_edge_count == water_count
    assert series.sequential_o_step_count == o_count
    assert series.sequential_h2o_step_count == water_count


@pytest.mark.parametrize("field,index", [("shoulder", 1), ("duplicate", 1)])
def test_shoulder_and_duplicate_pairs_are_excluded(field, index):
    kwargs = {f"{field}_index": index}
    audit = audit_oxygen_water_state_series(result((100.0, 116.0), **kwargs))
    assert not audit.edges
    assert len(audit.series) == 2


def test_peak_shape_relative_abundance_and_false_certainty_fields():
    series = max(audit_oxygen_water_state_series(result((100.0, 118.0, 134.0))).series,
                 key=lambda item: item.member_count)
    assert series.mass_ordered_apex_masses == (100.0, 118.0, 134.0)
    assert series.member_relative_apex_intensities == pytest.approx((1/3, 2/3, 1.0))
    assert series.all_members_independent_peaks is True
    assert series.peak_shape_support_class is PeakShapeSupportClass.STRONG_DISTINCT_PEAK_SUPPORT
    for name in ("oxidation_assigned", "hydration_assigned", "dehydration_assigned",
                 "thioamide_assigned", "structure_assigned", "reaction_direction_assigned",
                 "applied_to_formal_score", "applied_to_ranking",
                 "applied_to_candidate_filtering", "applied_to_final_consensus"):
        assert getattr(series, name) is False


def test_repeated_run_is_identical_and_input_is_not_mutated():
    source = result((100.0, 118.0, 134.0))
    before = repr(source)
    first = audit_oxygen_water_state_series(source)
    second = audit_oxygen_water_state_series(source)
    assert first == second
    assert repr(source) == before


@pytest.mark.skipif(not GLU_RUNTIME.is_file(), reason="runtime-only Glu workbook absent")
def test_glu_runtime_four_member_mixed_series_regression():
    manifest = load_sciex_sample_manifest(ROOT / "data/sciex_sample_manifest.yaml")
    source = resolve_profile_source(DEFAULT_PROFILE_REGISTRY, "GLU_UUC_WT_FULL_RECONSTRUCTED")
    loaded = load_profile_source(source, GLU_RUNTIME, sheet_name="旧WT_kenki_2")
    routing = route_profile_source_to_candidates(DEFAULT_PROFILE_REGISTRY, manifest, source.profile_source_id)
    family = analyze_loaded_profile_peak_families(
        loaded, routing, known_modifications=load_modifications(ROOT / "data/modifications.yaml"),
    )
    audit = audit_oxygen_water_state_series(family)
    target = next(series for series in audit.series
                  if series.member_count == 4 and 25291 <= series.lowest_mass <= 25293)
    assert target.member_apex_masses == pytest.approx((25292, 25310, 25326, 25342), abs=1.0)
    assert target.series_pattern is StateSeriesPattern.MIXED_H2O_AND_O_STEPS
    assert target.h2o_equivalent_edge_count == 1
    assert target.o_equivalent_edge_count >= 2
    assert target.highest_intensity_apex == pytest.approx(25326, abs=1.0)
    assert target.mass_span_da == pytest.approx(50, abs=1.0)
