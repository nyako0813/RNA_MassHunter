from dataclasses import replace
from pathlib import Path

import pytest

from rna_masshunter.modifications import load_modifications
from rna_masshunter.sciex_intact_cross_profile_audit import (
    ComparisonLayer,
    CrossProfileMassMatchClass,
    CrossProfileParameters,
    SelectedPeakClassification,
    ShapeSimilarityClass,
    audit_leu_cross_profiles,
    match_cross_profile_peaks,
)
from rna_masshunter.sciex_intact_peak_family import (
    PeakFamilyParameters,
    PeakFamilyPeak,
    PeakQualityClass,
    SciexIntactPeakFamilyResult,
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
UAA_RUNTIME = ROOT / ".cache/runtime_profiles/WT_LeuUAA(Full).txt"
UAG_RUNTIME = ROOT / ".cache/runtime_profiles/WT_LeuUAG(Full).txt"


def peak(prefix, index, mass, *, relative=0.5, integrated=None, fwhm=2.0, width=4.0,
         prominence=100.0, quality=PeakQualityClass.MAJOR_SHARP, centroid_offset=0.1):
    integrated = relative if integrated is None else integrated
    return PeakFamilyPeak(
        peak_id=f"{prefix}{index}", source_id=prefix, measurement_id=prefix,
        rna_identity=f"TRNA_LEU_{prefix}", apex_mass=mass,
        centroid_mass=None if centroid_offset is None else mass + centroid_offset,
        apex_intensity=relative * 1000, integrated_intensity=integrated * 2000,
        relative_apex_intensity=relative, relative_integrated_intensity=integrated,
        left_boundary_mass=mass - width/2 if width is not None else None,
        right_boundary_mass=mass + width/2 if width is not None else None,
        peak_width_da=width, fwhm_da=fwhm, prominence=prominence,
        relative_prominence=0.1, sharpness_score=10.0,
        nearest_peak_separation_da=None, peak_overlap_fraction=0.0,
        peak_detection_status="DETECTED", peak_quality_class=quality,
        selected_as_major_peak=True, possible_isotope_or_reconstruction_artifact=False,
        possible_shoulder=False, possible_duplicate_peak=False, possible_adduct=False,
        possible_output_convention_offset=False,
    )


def family_result(peaks):
    peaks = tuple(peaks)
    return SciexIntactPeakFamilyResult(
        "COMPLETED", "TEST", PeakFamilyParameters(), peaks, peaks, (), (), (), (), (),
    )


def match(uaa, uag, params=None):
    return match_cross_profile_peaks(
        uaa, uag, comparison_layer=ComparisonLayer.SELECTED_MAJOR_PEAK_COMPARISON,
        parameters=params,
    )


@pytest.mark.parametrize("difference,klass", [
    (0.5, CrossProfileMassMatchClass.STRICT),
    (0.8, CrossProfileMassMatchClass.EXPLORATORY),
])
def test_strict_and_exploratory_one_to_one_match(difference, klass):
    rows = match((peak("UAA", 1, 100),), (peak("UAG", 1, 100+difference),))
    assert len(rows) == 1 and rows[0].mass_match_class is klass


def test_unmatched_peaks_and_selected_classification():
    audit = audit_leu_cross_profiles(
        family_result((peak("UAA", 1, 100), peak("UAA", 2, 200))),
        family_result((peak("UAG", 1, 100), peak("UAG", 2, 300))),
    )
    classes = {status.classification for status in audit.selected_peak_statuses}
    assert SelectedPeakClassification.COMMON_SELECTED_PEAK in classes
    assert SelectedPeakClassification.UAA_ONLY_SELECTED_PEAK in classes
    assert SelectedPeakClassification.UAG_ONLY_SELECTED_PEAK in classes
    assert audit.uaa_summary.sample_specific_selected_peak_count == 1
    assert audit.uag_summary.sample_specific_selected_peak_count == 1


def test_assignment_is_unique_minimum_error_and_deterministic():
    uaa = (peak("UAA", 1, 100.0), peak("UAA", 2, 100.8))
    uag = (peak("UAG", 1, 100.2), peak("UAG", 2, 100.9))
    first = match(uaa, uag)
    second = match(tuple(reversed(uaa)), tuple(reversed(uag)))
    assert first == second
    assert len({row.uaa_peak_id for row in first}) == len(first)
    assert len({row.uag_peak_id for row in first}) == len(first)
    assert sum(row.apex_mass_difference_da for row in first) == pytest.approx(0.3)


def test_centroid_ratios_prominence_quality_and_relative_values():
    left = peak("UAA", 1, 100, relative=0.4, integrated=0.25, fwhm=2, width=4,
                prominence=50, quality=PeakQualityClass.MAJOR_SHARP, centroid_offset=0.1)
    right = peak("UAG", 1, 100.2, relative=0.8, integrated=0.5, fwhm=4, width=8,
                 prominence=100, quality=PeakQualityClass.MAJOR_BROAD, centroid_offset=0.3)
    row = match((left,), (right,))[0]
    assert row.centroid_mass_difference_da == pytest.approx(0.4)
    assert row.fwhm_ratio == pytest.approx(2)
    assert row.peak_width_ratio == pytest.approx(2)
    assert row.prominence_ratio == pytest.approx(2)
    assert row.relative_apex_intensity_ratio == pytest.approx(2)
    assert row.relative_integrated_intensity_ratio == pytest.approx(2)
    assert row.uaa_quality_class == "MAJOR_SHARP" and row.uag_quality_class == "MAJOR_BROAD"


@pytest.mark.parametrize("right_kwargs,expected", [
    ({}, ShapeSimilarityClass.HIGHLY_SIMILAR_PROFILE_PEAK),
    ({"relative": 0.2, "integrated": 0.2, "fwhm": 3, "width": 6},
     ShapeSimilarityClass.MODERATELY_SIMILAR_PROFILE_PEAK),
    ({"relative": 0.01, "integrated": 0.01, "fwhm": 20, "width": 40,
      "prominence": 1, "quality": PeakQualityClass.MAJOR_BROAD},
     ShapeSimilarityClass.MASS_MATCH_SHAPE_DIFFERENT),
    ({"fwhm": None, "width": None, "centroid_offset": None},
     ShapeSimilarityClass.MASS_ONLY_MATCH),
])
def test_shape_similarity_classes(right_kwargs, expected):
    left = peak("UAA", 1, 100)
    right = peak("UAG", 1, 100.1, **right_kwargs)
    assert match((left,), (right,))[0].shape_similarity_class is expected


def test_ambiguous_classification_when_two_candidates_are_equivalent():
    audit = audit_leu_cross_profiles(
        family_result((peak("UAA", 1, 100.0),)),
        family_result((peak("UAG", 1, 99.99), peak("UAG", 2, 100.01))),
    )
    assert audit.selected_major_matches[0].ambiguous_assignment is True
    assert any(status.classification is SelectedPeakClassification.AMBIGUOUS_CROSS_PROFILE_MATCH
               for status in audit.selected_peak_statuses)


def test_common_and_specific_intensity_fractions():
    audit = audit_leu_cross_profiles(
        family_result((peak("UAA", 1, 100, relative=0.75, integrated=0.6),
                       peak("UAA", 2, 200, relative=0.25, integrated=0.4))),
        family_result((peak("UAG", 1, 100, relative=0.5, integrated=0.7),
                       peak("UAG", 2, 300, relative=0.5, integrated=0.3))),
    )
    assert audit.uaa_summary.common_selected_apex_intensity_fraction == pytest.approx(0.75)
    assert audit.uaa_summary.sample_specific_apex_intensity_fraction == pytest.approx(0.25)
    assert audit.uaa_summary.common_selected_integrated_intensity_fraction == pytest.approx(0.6)
    assert audit.uag_summary.common_selected_integrated_intensity_fraction == pytest.approx(0.7)


def test_spearman_correlations_and_false_certainty_safeguards():
    uaa = tuple(peak("UAA", i, 100*i, relative=i/4, prominence=i*10, fwhm=i+1) for i in range(1, 4))
    uag = tuple(peak("UAG", i, 100*i, relative=i/4, prominence=i*10, fwhm=i+1) for i in range(1, 4))
    audit = audit_leu_cross_profiles(family_result(uaa), family_result(uag))
    assert audit.correlations.method == "SPEARMAN_RANK_CORRELATION"
    assert audit.correlations.apex_mass == pytest.approx(1)
    assert audit.correlations.relative_apex_intensity == pytest.approx(1)
    for record in (*audit.selected_major_matches, audit.uaa_summary, audit.uag_summary):
        assert record.shadow_analysis_only is True and record.mass_evidence_only is True
        for name in ("rna_identity_confirmed", "target_rna_identity_confirmed_by_mass",
                     "both_target_trnas_assigned", "common_peak_biological_identity_assigned",
                     "co_captured_rna_excluded", "reconstruction_artifact_excluded",
                     "background_component_excluded", "sequence_cocapture_ranking_performed",
                     "applied_to_formal_score", "applied_to_ranking",
                     "applied_to_candidate_filtering", "applied_to_final_consensus"):
            assert getattr(record, name) is False


def test_input_nonmutation():
    left = family_result((peak("UAA", 1, 100),))
    right = family_result((peak("UAG", 1, 100),))
    before = (repr(left), repr(right))
    audit_leu_cross_profiles(left, right)
    assert (repr(left), repr(right)) == before


@pytest.mark.skipif(not (UAA_RUNTIME.is_file() and UAG_RUNTIME.is_file()),
                    reason="runtime-only Leu profiles absent")
def test_leu_runtime_priority_common_peaks_and_specific_peaks_regression():
    manifest = load_sciex_sample_manifest(ROOT / "data/sciex_sample_manifest.yaml")
    modifications = load_modifications(ROOT / "data/modifications.yaml")
    results = []
    for source_id, path in (
        ("LEU_UAA_WT_FULL_RECONSTRUCTED", UAA_RUNTIME),
        ("LEU_UAG_WT_FULL_RECONSTRUCTED", UAG_RUNTIME),
    ):
        source = resolve_profile_source(DEFAULT_PROFILE_REGISTRY, source_id)
        loaded = load_profile_source(source, path)
        routing = route_profile_source_to_candidates(DEFAULT_PROFILE_REGISTRY, manifest, source_id)
        results.append(analyze_loaded_profile_peak_families(
            loaded, routing, known_modifications=modifications,
        ))
    audit = audit_leu_cross_profiles(*results)
    for target in (27946, 25054, 22797.5, 28766, 29693.5):
        assert any(abs(match.uaa_apex_mass-target) <= 1 and abs(match.uag_apex_mass-target) <= 1
                   for match in audit.selected_major_matches)
    assert audit.uaa_summary.common_selected_peak_count > 0
    assert audit.uaa_summary.sample_specific_selected_peak_count > 0
    assert audit.uag_summary.sample_specific_selected_peak_count > 0
