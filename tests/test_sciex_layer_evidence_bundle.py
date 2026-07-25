from copy import deepcopy
from dataclasses import replace
from enum import IntEnum
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pytest

import rna_masshunter.sciex_layer_evidence_bundle as bundle_serializer
from rna_masshunter.masses import load_base_masses
from rna_masshunter.sciex_intact_oxygen_water_state_audit import audit_oxygen_water_state_series
from rna_masshunter.sciex_intact_peak_family import (
    PeakFamilyParameters, PeakFamilyPeak, PeakQualityClass, SciexIntactPeakFamilyResult,
    build_delta_reference_registry, generate_delta_pairs, match_delta_pairs,
)
from rna_masshunter.sciex_layer_evidence_bundle import (
    BUNDLE_TYPE, LAYER_CONTRACTS, SAFEGUARDS, SCHEMA_VERSION,
    SERIALIZER_FORMAT_VERSION,
    LayerEvidenceBundleError, canonical_json_bytes,
    compare_layer_evidence_provenance, export_layer_evidence_bundle,
    load_layer_evidence_bundle, restore_layer_evidence_result,
    validate_layer_evidence_bundle,
)
from rna_masshunter.sciex_p1ap_nucleoside_ms2_identity_audit import (
    NucleosideCandidateMS2Summary, NucleosideMS2ProductMatch,
    P1APNucleosideMS2AuditResult, P1APNucleosideMS2Summary, ProcessedMS2Spectrum,
    audit_p1ap_nucleoside_ms2_identity, candidate_records_from_ms1_result,
)
from rna_masshunter.sciex_p1ap_nucleoside_state_audit import (
    NucleosideCandidateClass, P1APNucleosideStateAuditResult,
    audit_p1ap_nucleoside_state_series,
    match_nucleoside_ions_to_peaks,
)
from rna_masshunter.sciex_rna_cross_layer_evidence_reconciliation import (
    audit_rna_cross_layer_evidence_reconciliation,
)
from rna_masshunter.sciex_sample_manifest import load_sciex_sample_manifest
from rna_masshunter.sciex_t1_fragment_shadow_match import T1IonCandidate, T1IonMode
from rna_masshunter.sciex_t1_fragment_state_series_audit import (
    T1FragmentIonHypothesis, T1FragmentStateSeriesAuditResult,
    audit_t1_fragment_state_series, match_t1_fragment_ions_to_peaks,
)
from rna_masshunter.sciex_t1_replicate_consistency_audit import ReplicateRunPeak

ROOT = Path(__file__).parents[1]
COMMIT = "0e67175109a5c5bef9f3d57f31b8f66b18001c95"
CREATED = "2026-07-25T00:00:00Z"
RNA = {
    "name": "Mac_tRNA-Glu-UUC",
    "sequence": "GCUCCGGUAGUGUAGUCCGGCCAAUCAUUCCGGCCUUUCGAGCCGAAGACUCGGGUUCGAAUCCCGGCCGGAGCACCA",
    "anticodon": "UUC",
    "wobble_position": 37,
    "organism_group": "ARCHAEA",
    "species": "Methanosarcina acetivorans",
}


def _peak(identifier, mz):
    return ReplicateRunPeak(
        run_label="SYNTHETIC", peak_id=identifier, apex_mz=mz, centroid_mz=mz,
        raw_apex_intensity=100.0, normalized_apex_intensity=1.0,
        raw_integrated_intensity=10.0, normalized_integrated_intensity=1.0,
        relative_intensity=1.0, intensity_rank=1, prominence=0.1,
        relative_prominence=0.1, fwhm=0.02, left_bound_mz=mz - 0.02,
        right_bound_mz=mz + 0.02, supporting_ms1_scan_count=10,
        total_ms1_scan_count=100, scan_recurrence_fraction=0.1,
        first_supporting_scan_time=1.0, last_supporting_scan_time=2.0,
        detection_status="MAJOR_SHARP", detection_block_reasons=(),
    )


def _full_result():
    peaks = tuple(
        PeakFamilyPeak(
            peak_id=f"P{index}", source_id="FULL", measurement_id="SYNTHETIC",
            rna_identity=RNA["name"], apex_mass=mass, centroid_mass=mass + 0.1,
            apex_intensity=1000.0, integrated_intensity=2000.0,
            relative_apex_intensity=1.0, relative_integrated_intensity=1.0,
            left_boundary_mass=mass - 2.0, right_boundary_mass=mass + 2.0,
            peak_width_da=4.0, fwhm_da=3.0, prominence=100.0,
            relative_prominence=0.1, sharpness_score=10.0,
            nearest_peak_separation_da=None, peak_overlap_fraction=0.0,
            peak_detection_status="DETECTED", peak_quality_class=PeakQualityClass.MAJOR_SHARP,
            selected_as_major_peak=True, possible_isotope_or_reconstruction_artifact=False,
            possible_shoulder=False, possible_duplicate_peak=False, possible_adduct=False,
            possible_output_convention_offset=False,
        )
        for index, mass in enumerate((100.0, 118.0))
    )
    parameters = PeakFamilyParameters()
    pairs = generate_delta_pairs(peaks)
    references = build_delta_reference_registry()
    matches = match_delta_pairs(pairs, references, parameters=parameters)
    family = SciexIntactPeakFamilyResult(
        "COMPLETED", "FULL", parameters, peaks, peaks, pairs, references, matches, (), (),
    )
    return audit_oxygen_water_state_series(family)


def _t1_result(tmp_path):
    missing = tmp_path / "missing-t1.mzML"
    manifest = load_sciex_sample_manifest(ROOT / "data/sciex_sample_manifest.yaml")
    base_masses = load_base_masses(ROOT / "data/base_masses.yaml")
    empty = audit_t1_fragment_state_series(
        missing, RNA["sequence"], manifest=manifest, rna_identity_id="TRNA_GLU_UUC",
        base_masses=base_masses,
    )
    source = T1IonCandidate(
        "I_F1_1", "F1", "TRNA_GLU_UUC", "SYNTHETIC__CCA_CCA", "CCA", "ACG",
        1, 3, T1IonMode.NEGATIVE_DEPROTONATED, 1, 500.0,
        "MONOISOTOPIC_NEUTRAL", "MONOISOTOPIC_NEGATIVE_ION_MZ", -1,
    )
    ion = T1FragmentIonHypothesis(
        ion_hypothesis_id=source.ion_candidate_id, fragment_id="F1",
        fragment_sequence="ACG", start_position=1, end_position=3,
        theoretical_neutral_mass=501.007276466621, fragment_length=3,
        base_composition=(("A", 1), ("C", 1), ("G", 1), ("U", 0)),
        contains_g=True, cleavage_start="DIGEST_TERMINUS_UNKNOWN",
        cleavage_end="RNASE_T1_AFTER_G", terminal_chemistry="PROJECT_STANDARD_UNCHANGED",
        generation_status="GENERATED", generation_block_reasons=(),
        ion_mode=T1IonMode.NEGATIVE_DEPROTONATED, charge=1,
        adduct_hypothesis="[M-H]-", theoretical_mz=500.0,
        ion_hypothesis_status="ELIGIBLE_NEGATIVE_SOURCE",
        ion_hypothesis_block_reasons=(), source_candidate=source,
    )
    match = match_t1_fragment_ions_to_peaks((ion,), (_peak("T1P", 500.0),))[0]
    return replace(
        empty, ion_hypotheses=(ion,), fragment_matches=(match,),
        summary=replace(empty.summary, fragment_match_count=1, strict_match_count=1,
                        unambiguous_match_count=1),
    )


def _p1_ms1_result(tmp_path):
    empty = audit_p1ap_nucleoside_state_series(
        tmp_path / "missing-p1.mzML", project_root=ROOT, sequence="ACG",
    )
    ion = empty.ion_hypotheses[0]
    match = match_nucleoside_ions_to_peaks(
        (ion,), (_peak("P1P", ion.theoretical_mz),),
    )[0]
    return replace(
        empty, matches=(match,),
        summary=replace(empty.summary, candidate_match_count=1, strict_match_count=1,
                        unique_matched_peak_count=1, unambiguous_match_count=1),
    )


def _ms2_match():
    return NucleosideMS2ProductMatch(
        candidate_id="A", ms2_spectrum_id="MS2_1", product_ion_id="A__PROTONATED_BASE",
        product_ion_label="PROTONATED_BASE", product_ion_class="BASE_RELATED",
        theoretical_product_mz=136.06177, observed_product_mz=136.06177,
        delta_mz=0.0, absolute_delta_mz=0.0, ppm_error=0.0,
        observed_intensity=100.0, relative_intensity=1.0, intensity_rank=1,
        candidate_count_for_observed_peak=1, observed_peak_count_for_product_ion=1,
        match_ambiguity_status="UNAMBIGUOUS", match_quality_status="STRICT",
        match_block_reasons=(),
    )


def _ms2_candidate_summary(status):
    return NucleosideCandidateMS2Summary(
        candidate_id="A", candidate_name="adenosine", compatible_ms2_spectrum_count=1,
        usable_ms2_spectrum_count=1, collision_energy_count=1,
        theoretical_product_ion_count=3, recurrent_matched_product_ion_count=0,
        candidate_unique_recurrent_ion_count=0, shared_recurrent_ion_count=0,
        best_spectrum_id="MS2_1", best_spectrum_evidence_status="LOW",
        median_explained_intensity_fraction=0.1, median_top10_explained_fraction=0.1,
        median_mass_error=0.0, ms2_identity_evidence_status=status,
        ms2_identity_confidence="LOW",
        identity_ambiguity_before_ms2="IDENTITY_AMBIGUOUS",
        identity_ambiguity_after_ms2="CLASS_SUPPORTED_EXACT_IDENTITY_UNCONFIRMED",
        candidate_specific_ms2_rules_available=False, ms2_block_reasons=(),
    )


def _p1_ms2_result(tmp_path, *, with_spectrum=False):
    empty = audit_p1ap_nucleoside_ms2_identity(
        tmp_path / "missing-ms2.mzML", p1ap_ms1_result=None,
    )
    spectra = ()
    if with_spectrum:
        spectra = (ProcessedMS2Spectrum(
            ms2_spectrum_id="MS2_1", scan_time=1.0, raw_peak_count=1,
            positive_intensity_peak_count=1, zero_intensity_peak_count=0,
            negative_intensity_peak_count=0, filtered_peak_count=1,
            base_peak_mz=136.06177, base_peak_intensity=100.0, tic=100.0,
            mz_min=136.06177, mz_max=136.06177,
            profile_or_centroid_metadata="CENTROID", ms2_preprocessing_status="COMPLETED",
            ms2_preprocessing_block_reasons=(), peaks=((136.06177, 100.0),),
        ),)
    return replace(
        empty, spectrum_records=spectra, product_match_records=(_ms2_match(),),
        summary=replace(empty.summary, status="COMPLETED", product_match_count=1,
                        overall_evidence_status="PRODUCT_ION_EVIDENCE_PRESENT",
                        overall_confidence="LOW"),
    )


@pytest.fixture
def source_file(tmp_path):
    path = tmp_path / "sample.mzML"
    path.write_bytes(b"synthetic mzML source bytes")
    return path


@pytest.fixture
def results(tmp_path):
    return {
        "FULL": _full_result(),
        "T1": _t1_result(tmp_path),
        "P1AP_MS1": _p1_ms1_result(tmp_path),
        "P1AP_MS2": _p1_ms2_result(tmp_path),
    }


def _experiment(layer, *, independence="I1", shared="S1"):
    return {
        "condition_name": "SYNTHETIC", "digest_type": layer, "layer": layer,
        "independence_group": independence, "shared_source_group": shared,
    }


def _bundle(result, layer, source_file, **kwargs):
    return export_layer_evidence_bundle(
        result, layer=layer, source_path=source_file, rna=RNA,
        experiment=_experiment(layer), producer_commit=COMMIT,
        created_at_utc=kwargs.pop("created_at_utc", CREATED), **kwargs,
    )


def _rehash(bundle):
    payload = deepcopy(bundle)
    payload.pop("created_at_utc", None)
    payload.pop("canonical_payload_sha256", None)
    bundle["canonical_payload_sha256"] = sha256(canonical_json_bytes(payload)).hexdigest()


@pytest.mark.parametrize("layer", ["FULL", "T1", "P1AP_MS1", "P1AP_MS2"])
def test_layer_round_trip_returns_exact_production_class(layer, results, source_file):
    bundle = _bundle(results[layer], layer, source_file)
    restored = restore_layer_evidence_result(bundle)
    assert type(restored) is type(results[layer]) is LAYER_CONTRACTS[layer].result_class
    assert restored == results[layer]


def test_deterministic_serialization(results, source_file):
    left = _bundle(results["FULL"], "FULL", source_file)
    right = _bundle(results["FULL"], "FULL", source_file)
    assert canonical_json_bytes(left) == canonical_json_bytes(right)


def test_created_time_excluded_from_canonical_payload_sha(results, source_file):
    left = _bundle(results["FULL"], "FULL", source_file, created_at_utc="2026-01-01T00:00:00Z")
    right = _bundle(results["FULL"], "FULL", source_file, created_at_utc="2026-01-02T00:00:00Z")
    assert left["canonical_payload_sha256"] == right["canonical_payload_sha256"]


def test_canonical_payload_sha_detects_tampering(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    bundle["producer_commit"] = "tampered"
    with pytest.raises(LayerEvidenceBundleError, match="canonical payload SHA256"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_unknown_schema_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    bundle["schema_version"] = "999"
    with pytest.raises(LayerEvidenceBundleError, match="unknown schema_version"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_unknown_layer_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    bundle["layer"] = "UNKNOWN"
    with pytest.raises(LayerEvidenceBundleError, match="unknown layer"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_unknown_optional_result_key_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    bundle["optional_result_key"] = "identity_audit"
    with pytest.raises(LayerEvidenceBundleError, match="optional_result_key"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_malformed_json_rejected(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text("{not-json", encoding="utf-8")
    with pytest.raises(LayerEvidenceBundleError, match="malformed"):
        load_layer_evidence_bundle(path)


def test_missing_required_metadata_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    del bundle["source"]["sample_id"]
    with pytest.raises(LayerEvidenceBundleError, match="missing required fields"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


@pytest.mark.parametrize("layer,field", [
    ("FULL", "relations"), ("T1", "fragment_matches"),
    ("P1AP_MS1", "matches"), ("P1AP_MS2", "product_match_records"),
])
def test_empty_layer_result_rejected(layer, field, results, source_file):
    empty = replace(results[layer], **{field: ()})
    if layer == "P1AP_MS2":
        empty = replace(empty, candidate_summary_records=())
    with pytest.raises(LayerEvidenceBundleError, match="empty production result"):
        _bundle(empty, layer, source_file)


def test_p1ap_ms2_precursor_only_identity_summary_rejected(results, source_file):
    precursor_only = replace(
        results["P1AP_MS2"], product_match_records=(),
        candidate_summary_records=(_ms2_candidate_summary("MS2_PRECURSOR_COMPATIBLE_ONLY"),),
    )
    with pytest.raises(LayerEvidenceBundleError, match="empty production result"):
        _bundle(precursor_only, "P1AP_MS2", source_file)


def test_p1ap_ms2_substantive_identity_evidence_is_non_empty(results, source_file):
    identity = replace(
        results["P1AP_MS2"], product_match_records=(),
        candidate_summary_records=(_ms2_candidate_summary("MS2_SUPPORTS_CANONICAL_NUCLEOSIDE_CLASS"),),
    )
    bundle = _bundle(identity, "P1AP_MS2", source_file)
    assert bundle["validation"]["record_count"] == 1
    assert bundle["validation"]["non_empty"] is True


def test_t1_zero_state_family_is_valid_with_fragment_match(results, source_file):
    assert results["T1"].state_families == ()
    bundle = _bundle(results["T1"], "T1", source_file)
    assert bundle["validation"]["status"] == "PASSED"


def test_identity_audit_string_only_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    bundle["result"] = "identity_audit=PASSED"
    with pytest.raises(LayerEvidenceBundleError, match="encoded production object"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_wrong_source_sha_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    bundle["source"]["sha256"] = "0" * 64
    with pytest.raises(LayerEvidenceBundleError, match="SHA256 mismatch"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_wrong_source_size_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    bundle["source"]["size_bytes"] += 1
    with pytest.raises(LayerEvidenceBundleError, match="size mismatch"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_wrong_rna_sequence_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    expected = dict(RNA, sequence="DIFFERENT")
    with pytest.raises(LayerEvidenceBundleError, match="RNA sequence mismatch"):
        validate_layer_evidence_bundle(bundle, source_path=source_file, expected_rna=expected)


def test_wrong_anticodon_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    expected = dict(RNA, anticodon="AAA")
    with pytest.raises(LayerEvidenceBundleError, match="anticodon mismatch"):
        validate_layer_evidence_bundle(bundle, source_path=source_file, expected_rna=expected)


def test_sample_mismatch_is_reported(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file, sample_id="SAMPLE_A")
    with pytest.raises(LayerEvidenceBundleError, match="sample mismatch"):
        validate_layer_evidence_bundle(bundle, source_path=source_file, expected_sample_id="SAMPLE_B")


def test_04_new_t1_cannot_substitute_for_05_old_t1(results, tmp_path):
    old = tmp_path / "05 old T1.mzML"
    new = tmp_path / "04 new T1.mzML"
    old.write_bytes(b"same bytes")
    new.write_bytes(b"same bytes")
    bundle = _bundle(results["T1"], "T1", old, sample_id="OLD_T1")
    with pytest.raises(LayerEvidenceBundleError, match="cannot substitute"):
        validate_layer_evidence_bundle(bundle, source_path=new)


def test_safeguard_missing_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    del bundle["safeguards"]["formal_propagation"]
    with pytest.raises(LayerEvidenceBundleError, match="missing required fields|safeguard fields mismatch"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_safeguard_true_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    bundle["safeguards"]["formal_propagation"] = True
    with pytest.raises(LayerEvidenceBundleError, match="unsafe safeguard"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_raw_spectrum_array_rejected(results, source_file):
    bundle = _bundle(results["P1AP_MS2"], "P1AP_MS2", source_file)
    bundle["result"]["raw_spectrum_arrays"] = [[1.0, 2.0]]
    with pytest.raises(LayerEvidenceBundleError, match="raw spectrum/binary"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_binary_base64_payload_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    bundle["result"]["base64"] = "AAAA"
    with pytest.raises(LayerEvidenceBundleError, match="raw spectrum/binary"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_nonfinite_float_rejected_on_export(results, source_file):
    broken = replace(results["P1AP_MS2"], summary=replace(results["P1AP_MS2"].summary, runtime_seconds=float("nan")))
    with pytest.raises(LayerEvidenceBundleError, match="non-finite"):
        _bundle(broken, "P1AP_MS2", source_file)


def test_nonfinite_json_rejected_on_load(tmp_path):
    path = tmp_path / "nan.json"
    path.write_text('{"value":NaN}', encoding="utf-8")
    with pytest.raises(LayerEvidenceBundleError, match="non-finite"):
        load_layer_evidence_bundle(path)


def test_restored_objects_are_accepted_by_cross_layer_builder(results, source_file):
    restored = {
        layer: restore_layer_evidence_result(_bundle(result, layer, source_file))
        for layer, result in results.items()
    }
    output = audit_rna_cross_layer_evidence_reconciliation(
        full_length_result=restored["FULL"], t1_result=restored["T1"],
        p1ap_ms1_result=restored["P1AP_MS1"], p1ap_ms2_result=restored["P1AP_MS2"],
        runtime_context={"RNA_Identity": RNA["name"]},
    )
    assert output.nodes and output.consensus.formal_propagation is False


def test_p1ap_ms1_ms2_shared_source_grouping_preserved(results, source_file):
    ms1 = _bundle(results["P1AP_MS1"], "P1AP_MS1", source_file)
    ms2 = _bundle(results["P1AP_MS2"], "P1AP_MS2", source_file)
    comparison = compare_layer_evidence_provenance(ms1, ms2)
    assert comparison["same_raw_source"]
    assert comparison["shared_source_relationship"]
    assert comparison["same_independence_group"]
    assert comparison["p1ap_ms1_ms2_compatible"]
    assert not comparison["independent_support"]


def test_different_samples_are_not_counted_as_independent_support(results, tmp_path):
    left_source = tmp_path / "left.mzML"
    right_source = tmp_path / "right.mzML"
    left_source.write_bytes(b"left")
    right_source.write_bytes(b"right")
    left = export_layer_evidence_bundle(
        results["FULL"], layer="FULL", source_path=left_source, rna=RNA,
        experiment=_experiment("FULL", independence="I1", shared="S1"),
        producer_commit=COMMIT, created_at_utc=CREATED, sample_id="SAMPLE_A",
    )
    right = export_layer_evidence_bundle(
        results["T1"], layer="T1", source_path=right_source, rna=RNA,
        experiment=_experiment("T1", independence="I2", shared="S2"),
        producer_commit=COMMIT, created_at_utc=CREATED, sample_id="SAMPLE_B",
    )
    comparison = compare_layer_evidence_provenance(left, right)
    assert comparison["different_sample"]
    assert not comparison["independent_support"]


def test_formal_nonpropagation_preserved(results, source_file):
    for layer, result in results.items():
        bundle = _bundle(result, layer, source_file)
        assert bundle["safeguards"] == SAFEGUARDS
        restored = restore_layer_evidence_result(bundle)
        assert getattr(restored, "formal_propagation", False) is False


def test_ms2_peak_arrays_are_omitted_and_restore_empty(tmp_path, source_file):
    result = _p1_ms2_result(tmp_path, with_spectrum=True)
    bundle = _bundle(result, "P1AP_MS2", source_file)
    encoded = canonical_json_bytes(bundle)
    assert b'"peaks"' not in encoded
    restored = restore_layer_evidence_result(bundle)
    assert restored.spectrum_records[0].peaks == ()
    assert restored.summary == result.summary
    assert restored.product_match_records == result.product_match_records


def test_export_write_load_restore_without_raw_parser(results, source_file, tmp_path, monkeypatch):
    import rna_masshunter.sciex_profile_parser as profile_parser
    import rna_masshunter.sciex_intact_peak_detection as peak_detection

    def forbidden(*args, **kwargs):
        raise AssertionError("raw parser must not be invoked")

    monkeypatch.setattr(profile_parser, "parse_sciex_profile", forbidden)
    monkeypatch.setattr(peak_detection, "detect_sciex_intact_peaks", forbidden, raising=False)
    destination = tmp_path / "bundle.json"
    _bundle(results["FULL"], "FULL", source_file, output_path=destination)
    restored = load_layer_evidence_bundle(destination, source_path=source_file, restore=True)
    assert type(restored) is type(results["FULL"])


def test_schema_and_contract_fields_are_explicit(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    assert bundle["schema_version"] == SCHEMA_VERSION
    assert bundle["bundle_type"] == BUNDLE_TYPE
    assert bundle["optional_result_key"] == LAYER_CONTRACTS["FULL"].optional_result_key
    assert bundle["producer_name"] == LAYER_CONTRACTS["FULL"].producer_name
    assert bundle["validation"]["status"] == "PASSED"


def test_harness_schema_is_safely_rejected_as_unknown(tmp_path):
    path = tmp_path / "legacy.json"
    path.write_text(json.dumps({"schema_version": "rna-masshunter-cross-layer-bundle-v1"}), encoding="utf-8")
    with pytest.raises(LayerEvidenceBundleError, match="missing required fields|unknown schema"):
        load_layer_evidence_bundle(path)


_TRANSIENT_PROFILE_FIELDS = (
    "comparison_mz_grid", "comparison_raw_profile", "comparison_normalized_profile",
)


def _with_comparison_arrays(profile):
    return replace(
        profile,
        comparison_mz_grid=np.array([100.0, 100.5, 101.0], dtype=float),
        comparison_raw_profile=np.array([10.0, 20.0, 5.0], dtype=float),
        comparison_normalized_profile=np.array([0.5, 1.0, 0.25], dtype=float),
    )


def _find_encoded_dataclass(value, class_suffix):
    if isinstance(value, dict):
        if value.get("__bundle_type__") == "dataclass" and str(value.get("class", "")).endswith(class_suffix):
            return value
        for item in value.values():
            found = _find_encoded_dataclass(item, class_suffix)
            if found is not None:
                return found
    elif isinstance(value, list):
        for item in value:
            found = _find_encoded_dataclass(item, class_suffix)
            if found is not None:
                return found
    return None


def test_t1_nonempty_comparison_arrays_round_trip_as_transient(results, source_file):
    original = results["T1"]
    profile = _with_comparison_arrays(original.run_profile)
    result = replace(original, run_profile=profile)
    before = {name: getattr(profile, name).copy() for name in _TRANSIENT_PROFILE_FIELDS}

    bundle = _bundle(result, "T1", source_file)
    restored = restore_layer_evidence_result(bundle)

    assert type(restored) is type(result) is T1FragmentStateSeriesAuditResult
    assert restored == result
    assert all(getattr(restored.run_profile, name) is None for name in _TRANSIENT_PROFILE_FIELDS)
    assert all(np.array_equal(getattr(profile, name), before[name]) for name in _TRANSIENT_PROFILE_FIELDS)


def test_t1_transient_arrays_and_field_names_absent_from_canonical_payload(results, source_file):
    result = replace(results["T1"], run_profile=_with_comparison_arrays(results["T1"].run_profile))
    encoded = canonical_json_bytes(_bundle(result, "T1", source_file))
    assert all(name.encode() not in encoded for name in _TRANSIENT_PROFILE_FIELDS)
    assert b"100.5" not in encoded


def test_t1_with_transient_arrays_restores_into_cross_layer_builder(results, source_file):
    result = replace(results["T1"], run_profile=_with_comparison_arrays(results["T1"].run_profile))
    restored = restore_layer_evidence_result(_bundle(result, "T1", source_file))
    output = audit_rna_cross_layer_evidence_reconciliation(
        full_length_result=None, t1_result=restored, p1ap_ms1_result=None,
        p1ap_ms2_result=None, runtime_context={"RNA_Identity": RNA["name"]},
    )
    assert output.nodes
    assert any(node.layer == "T1" for node in output.nodes)


def test_unknown_numpy_array_remains_rejected(results, source_file):
    result = results["P1AP_MS1"]
    first = replace(result.candidates[0], source_candidate={"unexpected_profile": np.array([1.0, 2.0])})
    changed = replace(result, candidates=(first,) + result.candidates[1:])
    with pytest.raises(LayerEvidenceBundleError, match="numpy/raw arrays are forbidden"):
        _bundle(changed, "P1AP_MS1", source_file)


def test_manually_serialized_transient_field_is_rejected(results, source_file):
    result = replace(results["T1"], run_profile=_with_comparison_arrays(results["T1"].run_profile))
    bundle = _bundle(result, "T1", source_file)
    encoded_profile = _find_encoded_dataclass(bundle["result"], ".ReplicateRunPeakProfile")
    assert encoded_profile is not None
    encoded_profile["fields"]["comparison_mz_grid"] = {"__bundle_type__": "list", "items": [1.0]}
    with pytest.raises(LayerEvidenceBundleError, match="transient fields must not be serialized"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


@pytest.mark.parametrize("replacement", [
    {"base64": "AAAA"},
    {"binary_payload": "AAAA"},
    b"not-json-binary",
])
def test_malformed_transient_replacement_is_rejected(replacement, results, source_file):
    bundle = _bundle(results["T1"], "T1", source_file)
    encoded_profile = _find_encoded_dataclass(bundle["result"], ".ReplicateRunPeakProfile")
    assert encoded_profile is not None
    encoded_profile["fields"]["comparison_raw_profile"] = replacement
    with pytest.raises(LayerEvidenceBundleError, match="raw spectrum/binary|transient fields"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_p1ap_ms1_nonempty_comparison_arrays_round_trip_as_transient(results, source_file):
    original = results["P1AP_MS1"]
    nested = _with_comparison_arrays(original.run_profile.profile)
    result = replace(original, run_profile=replace(original.run_profile, profile=nested))
    bundle = _bundle(result, "P1AP_MS1", source_file)
    restored = restore_layer_evidence_result(bundle)

    assert type(restored) is type(result) is P1APNucleosideStateAuditResult
    assert restored == result
    assert all(getattr(restored.run_profile.profile, name) is None for name in _TRANSIENT_PROFILE_FIELDS)
    assert all(name.encode() not in canonical_json_bytes(bundle) for name in _TRANSIENT_PROFILE_FIELDS)


class _SyntheticIntEnum(IntEnum):
    ONE = 1


def _approved_test_enum(monkeypatch, enum_class):
    approved = dict(bundle_serializer._allowed_types())
    approved[bundle_serializer._type_name(enum_class)] = enum_class
    monkeypatch.setattr(bundle_serializer, "_allowed_types", lambda: approved)


def test_str_enum_uses_type_tag_and_restores_exact_class():
    encoded = bundle_serializer._encode(NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE)
    assert encoded == {
        "__bundle_type__": "enum",
        "class": bundle_serializer._type_name(NucleosideCandidateClass),
        "value": NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE.value,
    }
    restored = bundle_serializer._decode(encoded)
    assert type(restored) is NucleosideCandidateClass
    assert restored.value == NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE.value


def test_primitive_string_int_and_bool_remain_exact_primitive_types():
    for value in ("plain", 7, True, False):
        encoded = bundle_serializer._encode(value)
        restored = bundle_serializer._decode(encoded)
        assert encoded is value or encoded == value
        assert type(restored) is type(value)


def test_int_enum_uses_type_tag_and_restores_exact_class(monkeypatch):
    _approved_test_enum(monkeypatch, _SyntheticIntEnum)
    encoded = bundle_serializer._encode(_SyntheticIntEnum.ONE)
    assert encoded["__bundle_type__"] == "enum"
    assert encoded["value"] == 1
    restored = bundle_serializer._decode(encoded)
    assert type(restored) is _SyntheticIntEnum
    assert restored.value == 1


def test_unapproved_enum_type_is_rejected():
    with pytest.raises(LayerEvidenceBundleError, match="unapproved enum type"):
        bundle_serializer._encode(_SyntheticIntEnum.ONE)


def test_invalid_enum_type_tag_is_rejected():
    encoded = bundle_serializer._encode(NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE)
    encoded["class"] = bundle_serializer._type_name(P1APNucleosideStateAuditResult)
    with pytest.raises(LayerEvidenceBundleError, match="malformed enum encoding"):
        bundle_serializer._decode(encoded)


def test_invalid_enum_value_is_rejected():
    encoded = bundle_serializer._encode(NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE)
    encoded["value"] = "NOT_A_CANDIDATE_CLASS"
    with pytest.raises(LayerEvidenceBundleError, match="invalid enum value"):
        bundle_serializer._decode(encoded)


def test_unknown_enum_class_is_rejected():
    encoded = bundle_serializer._encode(NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE)
    encoded["class"] = "unknown.module.MissingEnum"
    with pytest.raises(LayerEvidenceBundleError, match="unapproved encoded class"):
        bundle_serializer._decode(encoded)


def test_legacy_schema_1_bundle_without_serializer_format_is_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    del bundle["serializer_format_version"]
    _rehash(bundle)
    with pytest.raises(LayerEvidenceBundleError, match="missing required fields"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_unknown_serializer_format_is_rejected(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    bundle["serializer_format_version"] = "legacy-primitive-first"
    _rehash(bundle)
    with pytest.raises(LayerEvidenceBundleError, match="unknown serializer_format_version"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_tampered_nested_enum_replaced_by_plain_string_is_rejected(results, source_file):
    bundle = _bundle(results["P1AP_MS1"], "P1AP_MS1", source_file)
    encoded_match = _find_encoded_dataclass(bundle["result"], ".NucleosidePeakMatch")
    assert encoded_match is not None
    encoded_match["fields"]["candidate_class"] = (
        NucleosideCandidateClass.NEUTRAL_NUCLEOSIDE.value
    )
    _rehash(bundle)
    with pytest.raises(LayerEvidenceBundleError, match="enum field type mismatch"):
        validate_layer_evidence_bundle(bundle, source_path=source_file)


def test_p1ap_ms1_nested_candidate_class_exact_type_and_ms2_setup(results, source_file):
    original = results["P1AP_MS1"]
    restored = restore_layer_evidence_result(_bundle(original, "P1AP_MS1", source_file))
    assert type(restored) is P1APNucleosideStateAuditResult
    assert type(restored.matches[0].candidate_class) is NucleosideCandidateClass
    assert restored.matches[0].candidate_class.value
    records = candidate_records_from_ms1_result(restored)
    assert records
    assert type(records[0].candidate_class) is str


def test_full_nested_enums_restore_exact_types(results, source_file):
    original = results["FULL"]
    restored = restore_layer_evidence_result(_bundle(original, "FULL", source_file))
    assert type(restored) is type(original)
    assert type(restored.references[0].reference_category) is type(
        original.references[0].reference_category
    )
    assert type(restored.relations[0].state_relation_class) is type(
        original.relations[0].state_relation_class
    )


def test_t1_nested_enums_restore_exact_types(results, source_file):
    original = results["T1"]
    restored = restore_layer_evidence_result(_bundle(original, "T1", source_file))
    assert type(restored) is type(original)
    assert type(restored.ion_hypotheses[0].ion_mode) is type(
        original.ion_hypotheses[0].ion_mode
    )
    assert type(restored.fragment_matches[0].ion_mode) is type(
        original.fragment_matches[0].ion_mode
    )


def test_p1ap_ms2_synthetic_nested_round_trip_after_enum_safe_ms1_setup(
    results, source_file,
):
    restored_ms1 = restore_layer_evidence_result(
        _bundle(results["P1AP_MS1"], "P1AP_MS1", source_file)
    )
    assert type(restored_ms1.matches[0].candidate_class) is NucleosideCandidateClass
    records = candidate_records_from_ms1_result(restored_ms1)
    ms2 = replace(results["P1AP_MS2"], candidates=records)
    restored_ms2 = restore_layer_evidence_result(_bundle(ms2, "P1AP_MS2", source_file))
    assert type(restored_ms2) is P1APNucleosideMS2AuditResult
    assert restored_ms2.candidates == records
    assert all(type(record.candidate_class) is str for record in restored_ms2.candidates)


def test_serializer_format_is_explicit(results, source_file):
    bundle = _bundle(results["FULL"], "FULL", source_file)
    assert bundle["schema_version"] == "1.0"
    assert bundle["serializer_format_version"] == SERIALIZER_FORMAT_VERSION
