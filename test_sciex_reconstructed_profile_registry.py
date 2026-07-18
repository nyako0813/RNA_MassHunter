from __future__ import annotations

from dataclasses import FrozenInstanceError, asdict, replace
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest
from openpyxl import Workbook

from rna_masshunter.sciex_reconstructed_profile_registry import (
    DEFAULT_PROFILE_REGISTRY,
    PROFILE_REGISTRY_APPLIED_TO_CANDIDATE_FILTERING,
    PROFILE_REGISTRY_APPLIED_TO_FINAL_CONSENSUS,
    PROFILE_REGISTRY_APPLIED_TO_FORMAL_SCORE,
    PROFILE_REGISTRY_APPLIED_TO_RANKING,
    LoadedProfileSource,
    ProfileRegistryValidationError,
    ProfileRoutingStatus,
    ProfileShadowComparisonRow,
    ProfileSourceType,
    ReconstructedProfileType,
    ShadowComparisonStatus,
    UnknownMassMetadata,
    build_profile_registry,
    compare_loaded_profile_shadow,
    load_profile_source,
    resolve_profile_source,
    route_profile_source_to_candidates,
    validate_profile_source_against_manifest,
)
from rna_masshunter.sciex_sample_manifest import load_sciex_sample_manifest

ROOT = Path(__file__).parent
MANIFEST_PATH = ROOT / "data" / "sciex_sample_manifest.yaml"
EXPECTED = {
    "LEU_UAA_WT_FULL_RECONSTRUCTED": (
        "LEU_UAA_WT_FULL", "LEU_UAA_WT", "TRNA_LEU_UAA", "WT_LeuUAA(Full).txt"
    ),
    "LEU_UAA_WT_T1_MZ": (
        "LEU_UAA_WT_T1", "LEU_UAA_WT", "TRNA_LEU_UAA", "UAA-T1.txt"
    ),
    "LEU_UAG_WT_FULL_RECONSTRUCTED": (
        "LEU_UAG_WT_FULL", "LEU_UAG_WT", "TRNA_LEU_UAG", "WT_LeuUAG(Full).txt"
    ),
    "LEU_UAG_WT_T1_MZ": (
        "LEU_UAG_WT_T1", "LEU_UAG_WT", "TRNA_LEU_UAG", "UAG-T1.txt"
    ),
    "GLU_UUC_WT_FULL_RECONSTRUCTED": (
        "GLU_UUC_WT_FULL", "GLU_UUC_WT", "TRNA_GLU_UUC", None
    ),
    "GLU_UUC_WT_P1_AP_MZML": (
        "GLU_UUC_WT_P1_AP", "GLU_UUC_WT", "TRNA_GLU_UUC", "Nsd01.mzML"
    ),
}


@pytest.fixture
def registry():
    return DEFAULT_PROFILE_REGISTRY


@pytest.fixture
def manifest():
    return load_sciex_sample_manifest(MANIFEST_PATH)


def write_text(path: Path, header: str, values=((100.0, 1.0), (101.0, 2.0), (102.0, 1.0))):
    lines = [header] + [f"{axis}\t{intensity}" for axis, intensity in values]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_default_registry_has_six_unique_explicit_sources(registry):
    assert registry.schema_version == 1
    assert len(registry.sources) == 6
    assert len({item.profile_source_id for item in registry.sources}) == 6
    assert len({item.measurement_id for item in registry.sources}) == 6
    actual = {
        item.profile_source_id: (
            item.measurement_id, item.sample_id, item.rna_identity_id, item.source_file_name
        )
        for item in registry.sources
    }
    assert actual == EXPECTED


def test_default_registry_contains_no_absolute_paths(registry):
    for source in registry.sources:
        for value in (source.source_file_name, source.source_workbook_name):
            assert value is None or not Path(value).is_absolute()
            assert value is None or "/mnt/" not in value
            assert value is None or ":\\" not in value


def test_duplicate_source_and_measurement_are_rejected(registry):
    with pytest.raises(ProfileRegistryValidationError, match="duplicate profile_source_id"):
        build_profile_registry((*registry.sources, registry.sources[0]))
    duplicate_measurement = replace(
        registry.sources[1], profile_source_id="UNIQUE_SOURCE",
        measurement_id=registry.sources[0].measurement_id,
    )
    with pytest.raises(ProfileRegistryValidationError, match="duplicate measurement_id"):
        build_profile_registry((*registry.sources, duplicate_measurement))


@pytest.mark.parametrize(
    ("source_id", "expected_type", "profile_type"),
    [
        (
            "LEU_UAA_WT_FULL_RECONSTRUCTED",
            ProfileSourceType.TEXT_NEUTRAL_PROFILE,
            ReconstructedProfileType.NEUTRAL_MASS_PROFILE,
        ),
        (
            "LEU_UAA_WT_T1_MZ",
            ProfileSourceType.TEXT_MZ_PROFILE,
            ReconstructedProfileType.MZ_PROFILE,
        ),
        (
            "GLU_UUC_WT_FULL_RECONSTRUCTED",
            ProfileSourceType.EXCEL_SHEET_NEUTRAL_PROFILE,
            ReconstructedProfileType.NEUTRAL_MASS_PROFILE,
        ),
        (
            "GLU_UUC_WT_P1_AP_MZML",
            ProfileSourceType.MZML_DIGEST,
            ReconstructedProfileType.DIGEST_MZML,
        ),
    ],
)
def test_source_and_profile_types_are_explicit(registry, source_id, expected_type, profile_type):
    source = resolve_profile_source(registry, source_id)
    assert source.source_type is expected_type
    assert source.profile_type is profile_type


def test_invalid_source_type_profile_combination_is_rejected(registry):
    invalid = replace(
        registry.sources[0],
        source_type=ProfileSourceType.TEXT_MZ_PROFILE,
    )
    with pytest.raises(ProfileRegistryValidationError, match="TEXT_MZ_PROFILE requires"):
        build_profile_registry((invalid,))


def test_unknown_profile_source_is_rejected(registry):
    with pytest.raises(KeyError, match="unknown profile_source_id"):
        resolve_profile_source(registry, "UNKNOWN")


def test_all_sources_validate_against_manifest(registry, manifest):
    for source in registry.sources:
        result = validate_profile_source_against_manifest(
            registry, manifest, source.profile_source_id
        )
        assert result.valid is True
        assert result.measurement_id == source.measurement_id
        assert result.sample_id == source.sample_id
        assert result.rna_identity_id == source.rna_identity_id


def test_unknown_measurement_reference_is_rejected(registry, manifest):
    changed = replace(registry.sources[0], measurement_id="UNKNOWN_MEASUREMENT")
    changed_registry = replace(registry, sources=(changed, *registry.sources[1:]))
    with pytest.raises(KeyError, match="unknown measurement_id"):
        validate_profile_source_against_manifest(
            changed_registry, manifest, changed.profile_source_id
        )


@pytest.mark.parametrize("field", ["sample_id", "rna_identity_id"])
def test_manifest_identity_conflicts_are_rejected(registry, manifest, field):
    changed = replace(registry.sources[0], **{field: "WRONG"})
    changed_registry = replace(registry, sources=(changed, *registry.sources[1:]))
    with pytest.raises(ProfileRegistryValidationError, match=field):
        validate_profile_source_against_manifest(
            changed_registry, manifest, changed.profile_source_id
        )


def test_only_three_full_sources_are_eligible(registry):
    candidate_eligible = [
        item.profile_source_id for item in registry.sources
        if item.eligible_for_intact_candidate_generation
    ]
    comparison_eligible = [
        item.profile_source_id for item in registry.sources
        if item.eligible_for_intact_mass_comparison
    ]
    expected = [
        "LEU_UAA_WT_FULL_RECONSTRUCTED",
        "LEU_UAG_WT_FULL_RECONSTRUCTED",
        "GLU_UUC_WT_FULL_RECONSTRUCTED",
    ]
    assert candidate_eligible == expected
    assert comparison_eligible == expected


def test_digest_sources_are_ineligible_with_reasons(registry):
    expected = {
        "LEU_UAA_WT_T1_MZ": "RNASE_T1_DIGEST_MZ_PROFILE",
        "LEU_UAG_WT_T1_MZ": "RNASE_T1_DIGEST_MZ_PROFILE",
        "GLU_UUC_WT_P1_AP_MZML": "P1_AP_DIGEST_NOT_INTACT_RNA",
    }
    for source_id, reason in expected.items():
        source = resolve_profile_source(registry, source_id)
        assert source.eligible_for_intact_candidate_generation is False
        assert source.eligible_for_intact_mass_comparison is False
        assert source.eligibility_reason == reason


@pytest.mark.parametrize(
    ("source_id", "identity", "count"),
    [
        ("LEU_UAA_WT_FULL_RECONSTRUCTED", "TRNA_LEU_UAA", 8),
        ("LEU_UAG_WT_FULL_RECONSTRUCTED", "TRNA_LEU_UAG", 8),
        ("GLU_UUC_WT_FULL_RECONSTRUCTED", "TRNA_GLU_UUC", 4),
    ],
)
def test_full_sources_route_correct_identity_and_candidate_count(
    registry, manifest, source_id, identity, count
):
    result = route_profile_source_to_candidates(registry, manifest, source_id)
    assert result.status is ProfileRoutingStatus.ROUTED
    assert result.rna_identity_id == identity
    assert len(result.candidates) == count
    assert {item.rna_identity_id for item in result.candidates} == {identity}


@pytest.mark.parametrize(
    "source_id",
    ["LEU_UAA_WT_T1_MZ", "LEU_UAG_WT_T1_MZ", "GLU_UUC_WT_P1_AP_MZML"],
)
def test_digest_routing_is_skipped_with_zero_candidates(registry, manifest, source_id):
    result = route_profile_source_to_candidates(registry, manifest, source_id)
    assert result.status is ProfileRoutingStatus.SKIPPED
    assert result.candidates == ()


def test_filename_and_sheet_name_do_not_route_identity(registry, manifest):
    source = resolve_profile_source(registry, "LEU_UAA_WT_FULL_RECONSTRUCTED")
    misleading = replace(source, source_file_name="TRNA_GLU_UUC.xlsx")
    changed_registry = replace(
        registry,
        sources=tuple(misleading if item is source else item for item in registry.sources),
    )
    result = route_profile_source_to_candidates(
        changed_registry, manifest, source.profile_source_id
    )
    assert result.rna_identity_id == "TRNA_LEU_UAA"


def test_text_mass_header_loads_as_neutral(tmp_path, registry):
    source = resolve_profile_source(registry, "LEU_UAA_WT_FULL_RECONSTRUCTED")
    path = write_text(tmp_path / source.source_file_name, "Mass\tIntensity")
    loaded = load_profile_source(source, path)
    assert loaded.header == ("Mass", "Intensity")
    assert loaded.profile_type is ReconstructedProfileType.NEUTRAL_MASS_PROFILE
    assert loaded.row_count == 3
    assert loaded.coordinates == (100.0, 101.0, 102.0)
    assert loaded.median_spacing == 1.0


def test_text_mass_charge_header_remains_mz_without_charge_assumption(tmp_path, registry):
    source = resolve_profile_source(registry, "LEU_UAA_WT_T1_MZ")
    path = write_text(tmp_path / source.source_file_name, "Mass/Charge\tIntensity")
    loaded = load_profile_source(source, path)
    assert loaded.header == ("Mass/Charge", "Intensity")
    assert loaded.profile_type is ReconstructedProfileType.MZ_PROFILE
    assert loaded.coordinates == (100.0, 101.0, 102.0)
    assert source.eligible_for_intact_mass_comparison is False


def test_wrong_text_profile_type_is_rejected(tmp_path, registry):
    source = resolve_profile_source(registry, "LEU_UAA_WT_T1_MZ")
    path = write_text(tmp_path / source.source_file_name, "Mass\tIntensity")
    with pytest.raises(ProfileRegistryValidationError, match="profile type mismatch"):
        load_profile_source(source, path)


def test_excel_loader_uses_only_registered_sheet(tmp_path, registry):
    source = resolve_profile_source(registry, "GLU_UUC_WT_FULL_RECONSTRUCTED")
    path = tmp_path / source.source_workbook_name
    workbook = Workbook()
    selected = workbook.active
    selected.title = source.source_sheet_name
    selected.append(("Mass", "Intensity"))
    selected.append((20000.0, 1.0))
    selected.append((20000.5, 2.0))
    other = workbook.create_sheet("other")
    other.append(("Mass/Charge", "Intensity"))
    other.append((500.0, 10.0))
    workbook.save(path)
    loaded = load_profile_source(source, path)
    assert loaded.sheet_name == source.source_sheet_name
    assert loaded.profile_type is ReconstructedProfileType.NEUTRAL_MASS_PROFILE
    assert loaded.coordinates == (20000.0, 20000.5)
    with pytest.raises(ProfileRegistryValidationError, match="wrong sheet"):
        load_profile_source(source, path, sheet_name="other")


def test_mzml_loading_is_metadata_only(tmp_path, registry):
    source = resolve_profile_source(registry, "GLU_UUC_WT_P1_AP_MZML")
    path = tmp_path / source.source_file_name
    path.write_bytes(b"not parsed mzML content")
    loaded = load_profile_source(source, path)
    assert loaded.profile_type is ReconstructedProfileType.DIGEST_MZML
    assert loaded.input_status == "METADATA_ONLY"
    assert loaded.row_count == 0
    assert loaded.source_size_bytes == len(b"not parsed mzML content")
    assert loaded.source_sha256 == sha256(b"not parsed mzML content").hexdigest()


def test_runtime_filename_mismatch_is_warning_not_identity_change(tmp_path, registry):
    source = resolve_profile_source(registry, "LEU_UAA_WT_FULL_RECONSTRUCTED")
    path = write_text(tmp_path / "renamed.txt", "Mass\tIntensity")
    loaded = load_profile_source(source, path)
    assert loaded.source_name_matches is False
    assert loaded.warnings == ("SOURCE_NAME_MISMATCH",)
    assert loaded.profile_source_id == source.profile_source_id


def synthetic_loaded(source, candidates):
    centers = [item.theoretical_mass for item in candidates]
    axis = np.arange(min(centers) - 20.0, max(centers) + 20.0, 0.5)
    intensity = np.ones_like(axis)
    for center in centers:
        intensity += 100.0 * np.exp(-0.5 * ((axis - center) / 1.0) ** 2)
    return LoadedProfileSource(
        source.profile_source_id,
        Path(source.source_file_name or source.source_workbook_name or "profile"),
        True,
        "0" * 64,
        1,
        ("Mass", "Intensity"),
        ReconstructedProfileType.NEUTRAL_MASS_PROFILE,
        "SUPPORTED_INPUT",
        tuple(float(value) for value in axis),
        tuple(float(value) for value in intensity),
        len(axis),
        float(axis.min()),
        float(axis.max()),
        0.5,
        None,
        (),
    )


def test_neutral_profile_shadow_comparison_uses_existing_detector_and_tolerances(registry, manifest):
    routing = route_profile_source_to_candidates(
        registry, manifest, "GLU_UUC_WT_FULL_RECONSTRUCTED"
    )
    loaded = synthetic_loaded(routing.profile_source, routing.candidates)
    result = compare_loaded_profile_shadow(loaded, routing)
    assert result.status is ShadowComparisonStatus.COMPLETED
    assert result.detection_result is not None
    assert len(result.detection_result.peaks) >= 2
    assert result.strict_tolerance_da == 1.0
    assert result.exploratory_tolerance_da == 5.0
    assert result.detail_rows
    assert all(row.rna_identity_id == "TRNA_GLU_UUC" for row in result.detail_rows)


@pytest.mark.parametrize(
    "source_id",
    ["LEU_UAA_WT_T1_MZ", "LEU_UAG_WT_T1_MZ", "GLU_UUC_WT_P1_AP_MZML"],
)
def test_digest_shadow_comparison_is_skipped_with_no_details(registry, manifest, source_id):
    routing = route_profile_source_to_candidates(registry, manifest, source_id)
    source = routing.profile_source
    loaded = LoadedProfileSource(
        source.profile_source_id, Path(source.source_file_name or "digest"), True,
        "0" * 64, 0, (), source.profile_type, "METADATA_ONLY", (), (), 0,
        None, None, None, None, (),
    )
    result = compare_loaded_profile_shadow(loaded, routing)
    assert result.status is ShadowComparisonStatus.SKIPPED
    assert result.detail_rows == ()
    assert result.candidates == ()


def test_loaded_routing_identity_conflict_is_rejected(registry, manifest):
    routing = route_profile_source_to_candidates(
        registry, manifest, "LEU_UAA_WT_FULL_RECONSTRUCTED"
    )
    loaded = replace(
        synthetic_loaded(routing.profile_source, routing.candidates),
        profile_source_id="LEU_UAG_WT_FULL_RECONSTRUCTED",
    )
    with pytest.raises(ProfileRegistryValidationError, match="conflict"):
        compare_loaded_profile_shadow(loaded, routing)


def test_shadow_rows_keep_unknown_mass_metadata_and_false_certainty(registry, manifest):
    routing = route_profile_source_to_candidates(
        registry, manifest, "GLU_UUC_WT_FULL_RECONSTRUCTED"
    )
    result = compare_loaded_profile_shadow(
        synthetic_loaded(routing.profile_source, routing.candidates), routing
    )
    for row in result.detail_rows:
        assert row.observed_mass_type == "UNKNOWN"
        assert row.mass_definition_compatibility == "UNKNOWN"
        assert row.theoretical_mass_type == "MONOISOTOPIC_NEUTRAL"
        assert row.calibration_applied is False
        assert row.mass_match_only is True
        assert row.unmodified_candidate is True
        assert row.cca_state_confirmed is False
        assert row.terminal_state_confirmed is False
        assert row.structure_identity_assigned is False
        assert row.position_assigned is False
        assert row.modification_assigned is False
        assert row.biological_cause_assigned is False
        assert row.rnase_t_assigned is False
        assert row.applied_to_formal_score is False
        assert row.applied_to_ranking is False
        assert row.applied_to_candidate_filtering is False
        assert row.applied_to_final_consensus is False


def test_source_metadata_registry_manifest_and_candidates_are_immutable(registry, manifest):
    source = registry.sources[0]
    registry_before = asdict(registry)
    manifest_before = asdict(manifest)
    source_before = asdict(source)
    routing = route_profile_source_to_candidates(registry, manifest, source.profile_source_id)
    assert asdict(registry) == registry_before
    assert asdict(manifest) == manifest_before
    assert asdict(source) == source_before
    with pytest.raises(FrozenInstanceError):
        routing.candidates[0].candidate_priority = 9


def test_registry_and_rows_are_nonformal(registry):
    assert PROFILE_REGISTRY_APPLIED_TO_FORMAL_SCORE is False
    assert PROFILE_REGISTRY_APPLIED_TO_RANKING is False
    assert PROFILE_REGISTRY_APPLIED_TO_CANDIDATE_FILTERING is False
    assert PROFILE_REGISTRY_APPLIED_TO_FINAL_CONSENSUS is False
    assert all(item.observed_mass_type is UnknownMassMetadata.UNKNOWN for item in registry.sources)
    assert all(
        item.mass_definition_compatibility is UnknownMassMetadata.UNKNOWN
        for item in registry.sources
    )
