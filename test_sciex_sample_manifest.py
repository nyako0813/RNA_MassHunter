from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path
import re

import pytest
import yaml

from rna_masshunter.sciex_sample_manifest import (
    AnalyteLevel,
    AnticodonOrientation,
    CCAStatus,
    DigestionEnzyme,
    ExperimentType,
    IdentityStatus,
    ManifestValidationError,
    PhosphataseTreatment,
    SAMPLE_MANIFEST_APPLIED_TO_CANDIDATE_FILTERING,
    SAMPLE_MANIFEST_APPLIED_TO_FINAL_CONSENSUS,
    SAMPLE_MANIFEST_APPLIED_TO_FORMAL_SCORE,
    SAMPLE_MANIFEST_APPLIED_TO_RANKING,
    SequenceStatus,
    TerminalStateStatus,
    filter_measurements_by_experiment_type,
    get_measurement,
    get_measurements_for_sample,
    get_rna_identity,
    get_paired_measurements,
    get_samples_for_rna_identity,
    load_sciex_sample_manifest,
    resolve_measurement_identity,
    resolve_measurement_path,
    validate_measurement_against_observed_metadata,
    validate_measurement_files,
)

ROOT = Path(__file__).parent
MANIFEST_PATH = ROOT / "data" / "sciex_sample_manifest.yaml"
LEU_UAA_SEQUENCE = "GCGAGGGUUGCCCAGCCAGGCCAAAGGCGCCAGACUUAAGAUCUGGUAUCGAAGGAUUUCGUGGGUUCGAAUCCCACCCCUCGCA"
LEU_UAG_SEQUENCE = "GCGAGGGUUGCCCAGCUAGGUCAAAGGCGAUGGGCUUAGGACCCAUUUUCGUAGGAAUUCGUGCGUUCGAAUCGCACCCCUCGCA"

EXPECTED_MEASUREMENTS = {
    "LEU_UAA_WT_FULL": ("LEU_UAA_WT", "leu_uaa_full", "12Old UAA.mzML"),
    "LEU_UAA_WT_T1": ("LEU_UAA_WT", "leu_uaa_t1", "03_LeuUAA_T1.mzML"),
    "GLU_UUC_WT_FULL": ("GLU_UUC_WT", "glu_uuc_full", "01kenki.mzML"),
    "GLU_UUC_WT_P1_AP": ("GLU_UUC_WT", "glu_uuc_p1_ap", "Nsd01.mzML"),
    "LEU_UAG_WT_FULL": ("LEU_UAG_WT", "leu_uag_full", "11Old UAG.mzML"),
    "LEU_UAG_WT_T1": ("LEU_UAG_WT", "leu_uag_t1", "04_LeuUAG_T1.mzML"),
}


@pytest.fixture
def manifest():
    return load_sciex_sample_manifest(MANIFEST_PATH)


@pytest.fixture
def payload():
    return yaml.safe_load(MANIFEST_PATH.read_text(encoding="utf-8"))


def write_payload(tmp_path: Path, payload: dict) -> Path:
    path = tmp_path / "manifest.yaml"
    path.write_text(yaml.safe_dump(payload, sort_keys=False, allow_unicode=True), encoding="utf-8")
    return path


def confirmed_identity(payload: dict, sequence: str = "ac gu cca") -> dict:
    raw = payload["rna_identities"][0]
    raw.update({
        "anticodon": "ACG",
        "anticodon_orientation": "FIVE_TO_THREE",
        "wobble_sequence_index_1based": 1,
        "anticodon_start_index_1based": 1,
        "anticodon_end_index_1based": 3,
        "sequence": sequence,
        "sequence_source": "curated reference accession TEST:1",
        "sequence_status": "CONFIRMED",
        "sequence_orientation": "FIVE_TO_THREE",
        "sequence_alphabet": "RNA",
        "sequence_sha256": None,
        "sequence_length": None,
        "ends_with_cca": None,
        "cca_status": "CONFIRMED_PRESENT",
        "mature_rna_status": "MATURE_RNA",
    })
    return raw


def test_valid_manifest_load_and_schema_version(manifest):
    assert manifest.schema_version == 1
    assert len(manifest.rna_identities) == 3
    assert len(manifest.samples) == 3
    assert len(manifest.measurements) == 6


def test_load_is_deterministic():
    assert load_sciex_sample_manifest(MANIFEST_PATH) == load_sciex_sample_manifest(MANIFEST_PATH)


def test_missing_file_has_useful_error(tmp_path):
    with pytest.raises(ManifestValidationError, match="manifest not found"):
        load_sciex_sample_manifest(tmp_path / "missing.yaml")


def test_malformed_yaml_has_useful_error(tmp_path):
    path = tmp_path / "bad.yaml"
    path.write_text("measurements: [\n", encoding="utf-8")
    with pytest.raises(ManifestValidationError, match="cannot load manifest"):
        load_sciex_sample_manifest(path)


def test_unknown_root_field_is_rejected(tmp_path, payload):
    payload["unexpected"] = True
    with pytest.raises(ManifestValidationError, match="unknown field"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


@pytest.mark.parametrize(
    ("section", "field", "value"),
    [
        ("rna_identities", "sequence_status", "GUESSED"),
        ("samples", "identity_status", "CERTAIN"),
        ("measurements", "experiment_type", "P1_ONLY"),
        ("measurements", "expected_analyte_level", "PEPTIDE"),
    ],
)
def test_invalid_enums_are_rejected(tmp_path, payload, section, field, value):
    payload[section][0][field] = value
    with pytest.raises(ManifestValidationError, match="unknown value"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


@pytest.mark.parametrize(
    ("section", "reference_field", "value", "message"),
    [
        ("samples", "rna_identity_id", "MISSING_RNA", "unknown rna_identity_id"),
        ("measurements", "sample_id", "MISSING_SAMPLE", "unknown sample_id"),
    ],
)
def test_unknown_references_are_rejected(tmp_path, payload, section, reference_field, value, message):
    payload[section][0][reference_field] = value
    with pytest.raises(ManifestValidationError, match=message):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


@pytest.mark.parametrize(
    ("section", "id_field", "message"),
    [
        ("rna_identities", "rna_identity_id", "duplicate rna_identity_id"),
        ("samples", "sample_id", "duplicate sample_id"),
        ("measurements", "measurement_id", "duplicate measurement_id"),
        ("measurements", "input_alias", "duplicate input_alias"),
    ],
)
def test_duplicate_ids_and_aliases_are_rejected(tmp_path, payload, section, id_field, message):
    duplicate = deepcopy(payload[section][0])
    if id_field != "input_alias":
        duplicate[id_field] = payload[section][0][id_field]
    else:
        duplicate["measurement_id"] = "UNIQUE_MEASUREMENT"
    payload[section].append(duplicate)
    with pytest.raises(ManifestValidationError, match=message):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


def test_current_measurement_ids_aliases_and_filenames_are_exact(manifest):
    actual = {
        item.measurement_id: (item.sample_id, item.input_alias, item.source_file_name)
        for item in manifest.measurements
    }
    assert actual == EXPECTED_MEASUREMENTS
    assert len({item.input_alias for item in manifest.measurements}) == 6


@pytest.mark.parametrize("measurement_id", EXPECTED_MEASUREMENTS)
def test_each_current_measurement_exists(manifest, measurement_id):
    assert get_measurement(manifest, measurement_id).measurement_id == measurement_id


def test_each_full_and_digest_belongs_to_same_sample(manifest):
    assert {item.measurement_id for item in get_measurements_for_sample(manifest, "LEU_UAA_WT")} == {
        "LEU_UAA_WT_FULL", "LEU_UAA_WT_T1"
    }
    assert {item.measurement_id for item in get_measurements_for_sample(manifest, "GLU_UUC_WT")} == {
        "GLU_UUC_WT_FULL", "GLU_UUC_WT_P1_AP"
    }
    assert {item.measurement_id for item in get_measurements_for_sample(manifest, "LEU_UAG_WT")} == {
        "LEU_UAG_WT_FULL", "LEU_UAG_WT_T1"
    }


def test_leu_isoacceptors_are_distinct_and_never_cross_paired(manifest):
    uaa = resolve_measurement_identity(manifest, "LEU_UAA_WT_T1")
    uag = resolve_measurement_identity(manifest, "LEU_UAG_WT_FULL")
    assert (uaa.rna_identity_id, uaa.expected_anticodon) == ("TRNA_LEU_UAA", "UAA")
    assert (uag.rna_identity_id, uag.expected_anticodon) == ("TRNA_LEU_UAG", "UAG")
    assert [item.measurement_id for item in get_paired_measurements(manifest, "LEU_UAA_WT_T1")] == [
        "LEU_UAA_WT_FULL"
    ]


def test_full_length_treatment_and_analyte(manifest):
    measurements = filter_measurements_by_experiment_type(manifest, ExperimentType.FULL_LENGTH)
    assert len(measurements) == 3
    assert all(item.digestion_enzyme is DigestionEnzyme.NONE for item in measurements)
    assert all(item.phosphatase_treatment is PhosphataseTreatment.NONE for item in measurements)
    assert all(item.expected_analyte_level is AnalyteLevel.INTACT_RNA for item in measurements)


def test_t1_treatment_and_analyte(manifest):
    measurements = filter_measurements_by_experiment_type(manifest, ExperimentType.RNASE_T1_DIGEST)
    assert len(measurements) == 2
    assert all(item.digestion_enzyme is DigestionEnzyme.RNASE_T1 for item in measurements)
    assert all(item.phosphatase_treatment is PhosphataseTreatment.NONE for item in measurements)
    assert all(item.expected_analyte_level is AnalyteLevel.OLIGONUCLEOTIDE for item in measurements)


def test_p1_ap_treatment_and_analyte(manifest):
    measurements = filter_measurements_by_experiment_type(manifest, ExperimentType.NUCLEASE_P1_AP_DIGEST)
    assert [item.measurement_id for item in measurements] == ["GLU_UUC_WT_P1_AP"]
    item = measurements[0]
    assert item.digestion_enzyme is DigestionEnzyme.NUCLEASE_P1
    assert item.phosphatase_treatment is PhosphataseTreatment.ALKALINE_PHOSPHATASE
    assert item.expected_analyte_level is AnalyteLevel.NUCLEOSIDE


@pytest.mark.parametrize(
    ("experiment_index", "field", "value"),
    [
        (0, "digestion_enzyme", "RNASE_T1"),
        (1, "digestion_enzyme", "NUCLEASE_P1"),
        (3, "phosphatase_treatment", "NONE"),
        (0, "expected_analyte_level", "NUCLEOSIDE"),
        (1, "expected_analyte_level", "INTACT_RNA"),
    ],
)
def test_inconsistent_experiment_treatment_is_rejected(tmp_path, payload, experiment_index, field, value):
    payload["measurements"][experiment_index][field] = value
    with pytest.raises(ManifestValidationError, match="requires digestion"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


def test_confirmed_canonical_sequence_is_normalized_and_derived(tmp_path, payload):
    confirmed_identity(payload)
    identity = load_sciex_sample_manifest(write_payload(tmp_path, payload)).rna_identities[0]
    assert identity.sequence == "ACGUCCA"
    assert identity.sequence_length == 7
    assert identity.sequence_sha256 == sha256(b"ACGUCCA").hexdigest()
    assert identity.ends_with_cca is True
    assert identity.sequence_status is SequenceStatus.CONFIRMED


def test_unknown_sequence_is_allowed_without_derived_values(tmp_path, payload):
    raw = payload["rna_identities"][1]
    raw.update({
        "wobble_sequence_index_1based": None,
        "anticodon_start_index_1based": None,
        "anticodon_end_index_1based": None,
        "sequence": None,
        "sequence_source": None,
        "sequence_status": "UNKNOWN",
        "sequence_orientation": "UNKNOWN",
        "sequence_alphabet": "UNKNOWN",
        "registered_sequence_cca_mode": "UNKNOWN",
        "registered_cca_tail_state": None,
        "cca_status": "UNKNOWN",
        "sequence_notes": None,
    })
    identity = load_sciex_sample_manifest(write_payload(tmp_path, payload)).rna_identities[1]
    assert identity.sequence is None
    assert identity.sequence_sha256 is None
    assert identity.sequence_length is None and identity.ends_with_cca is None


@pytest.mark.parametrize("sequence", ["ACGX", "ACGT", "ACG-", ""])
def test_invalid_sequence_is_rejected_without_implicit_t_to_u(tmp_path, payload, sequence):
    confirmed_identity(payload, sequence)
    with pytest.raises(ManifestValidationError, match="sequence"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


def test_supplied_sequence_hash_and_length_must_match(tmp_path, payload):
    raw = confirmed_identity(payload)
    raw["sequence_sha256"] = "0" * 64
    with pytest.raises(ManifestValidationError, match="sequence_sha256"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))
    raw["sequence_sha256"] = None
    raw["sequence_length"] = 99
    with pytest.raises(ManifestValidationError, match="sequence_length"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


def test_cca_consistency_is_deterministic(tmp_path, payload):
    raw = confirmed_identity(payload)
    raw["ends_with_cca"] = True
    identity = load_sciex_sample_manifest(write_payload(tmp_path, payload)).rna_identities[0]
    assert identity.ends_with_cca is True
    assert identity.cca_status is CCAStatus.CONFIRMED_PRESENT


@pytest.mark.parametrize(
    ("sequence", "ends_with_cca", "cca_status"),
    [
        ("ACGUCCA", False, "CONFIRMED_PRESENT"),
        ("ACGUCCA", True, "CONFIRMED_ABSENT"),
        ("ACGU", False, "CONFIRMED_PRESENT"),
    ],
)
def test_cca_contradictions_are_rejected(tmp_path, payload, sequence, ends_with_cca, cca_status):
    raw = confirmed_identity(payload, sequence)
    raw["ends_with_cca"] = ends_with_cca
    raw["cca_status"] = cca_status
    with pytest.raises(ManifestValidationError, match="ends_with_cca|cca_status"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("sequence_source", None),
        ("sequence_orientation", "UNKNOWN"),
        ("sequence_alphabet", "UNKNOWN"),
        ("cca_status", "LIKELY_PRESENT"),
    ],
)
def test_confirmed_sequence_requires_provenance_and_explicit_state(tmp_path, payload, field, value):
    raw = confirmed_identity(payload)
    raw[field] = value
    with pytest.raises(ManifestValidationError, match="CONFIRMED sequence"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


def test_current_anticodon_orientations_reflect_only_supplied_information(manifest):
    assert get_rna_identity(manifest, "TRNA_LEU_UAA").anticodon_orientation is AnticodonOrientation.FIVE_TO_THREE
    assert get_rna_identity(manifest, "TRNA_LEU_UAG").anticodon_orientation is AnticodonOrientation.FIVE_TO_THREE
    assert get_rna_identity(manifest, "TRNA_GLU_UUC").anticodon_orientation is AnticodonOrientation.UNCONFIRMED


def test_terminal_states_remain_unknown_and_digest_does_not_reuse_them(manifest):
    for item in manifest.measurements:
        assert item.five_prime_state is None and item.three_prime_state is None
        assert item.terminal_state_status is TerminalStateStatus.UNKNOWN
        assert item.terminal_state_required is (item.experiment_type is ExperimentType.FULL_LENGTH)


def test_digest_terminal_state_metadata_is_rejected(tmp_path, payload):
    payload["measurements"][1]["five_prime_state"] = "MONOPHOSPHATE"
    with pytest.raises(ManifestValidationError, match="cannot reuse intact terminal-state"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


def test_lookup_sample_to_measurements_is_ordered(manifest):
    assert [item.measurement_id for item in get_measurements_for_sample(manifest, "LEU_UAA_WT")] == [
        "LEU_UAA_WT_FULL", "LEU_UAA_WT_T1"
    ]


def test_lookup_identity_to_samples(manifest):
    assert [item.sample_id for item in get_samples_for_rna_identity(manifest, "TRNA_GLU_UUC")] == [
        "GLU_UUC_WT"
    ]


@pytest.mark.parametrize(
    ("measurement_id", "paired_id"),
    [
        ("LEU_UAA_WT_FULL", "LEU_UAA_WT_T1"),
        ("LEU_UAA_WT_T1", "LEU_UAA_WT_FULL"),
        ("GLU_UUC_WT_FULL", "GLU_UUC_WT_P1_AP"),
        ("GLU_UUC_WT_P1_AP", "GLU_UUC_WT_FULL"),
        ("LEU_UAG_WT_FULL", "LEU_UAG_WT_T1"),
        ("LEU_UAG_WT_T1", "LEU_UAG_WT_FULL"),
    ],
)
def test_bidirectional_full_digest_pairing(manifest, measurement_id, paired_id):
    assert [item.measurement_id for item in get_paired_measurements(manifest, measurement_id)] == [paired_id]


def test_experiment_filter_accepts_enum_value_and_preserves_order(manifest):
    assert [item.measurement_id for item in filter_measurements_by_experiment_type(manifest, "FULL_LENGTH")] == [
        "LEU_UAA_WT_FULL", "GLU_UUC_WT_FULL", "LEU_UAG_WT_FULL"
    ]


@pytest.mark.parametrize(
    ("function", "unknown"),
    [
        (get_measurements_for_sample, "UNKNOWN_SAMPLE"),
        (get_samples_for_rna_identity, "UNKNOWN_IDENTITY"),
        (get_paired_measurements, "UNKNOWN_MEASUREMENT"),
    ],
)
def test_lookup_unknown_ids_raise_key_error(manifest, function, unknown):
    with pytest.raises(KeyError, match="unknown"):
        function(manifest, unknown)


def test_source_filename_and_optional_hash_are_preserved(manifest):
    measurement = get_measurement(manifest, "GLU_UUC_WT_P1_AP")
    assert measurement.source_file_name == "Nsd01.mzML"
    assert measurement.source_file_sha256 is None


def test_filename_does_not_confirm_sample_identity(manifest):
    metadata = resolve_measurement_identity(manifest, "LEU_UAG_WT_FULL")
    assert metadata.expected_source_file_name == "11Old UAG.mzML"
    assert metadata.identity_status is IdentityStatus.UNCONFIRMED
    assert metadata.expected_sequence_sha256 == sha256(LEU_UAG_SEQUENCE.encode("ascii")).hexdigest()


def test_tracked_yaml_contains_no_absolute_windows_or_wsl_paths():
    text = MANIFEST_PATH.read_text(encoding="utf-8")
    assert "/mnt/" not in text and "/home/" not in text
    assert re.search(r"[A-Za-z]:[\\/]", text) is None


@pytest.mark.parametrize("absolute", ["/home/user/input.mzML", "/mnt/c/input.mzML", r"C:\\Data\\input.mzML"])
def test_absolute_path_hints_are_rejected(tmp_path, payload, absolute):
    payload["measurements"][0]["source_file_path_hint"] = absolute
    with pytest.raises(ManifestValidationError, match="absolute paths are forbidden"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


def test_resolve_measurement_path_returns_path_without_mutation(manifest):
    mapping = {"LEU_UAA_WT_FULL": "測定 データ/12Old UAA.mzML"}
    before = dict(mapping)
    manifest_before = asdict(manifest)
    result = resolve_measurement_path(manifest, "LEU_UAA_WT_FULL", mapping)
    assert result == Path("測定 データ/12Old UAA.mzML")
    assert isinstance(result, Path)
    assert mapping == before
    assert asdict(manifest) == manifest_before


def test_resolve_measurement_path_does_not_open_file(manifest, monkeypatch):
    def fail_open(*args, **kwargs):
        raise AssertionError("resolve_measurement_path opened a file")

    monkeypatch.setattr(Path, "open", fail_open)
    assert resolve_measurement_path(manifest, "LEU_UAA_WT_FULL", {"LEU_UAA_WT_FULL": "missing.mzML"}) == Path(
        "missing.mzML"
    )


def test_resolve_missing_mapping_and_unknown_measurement(manifest):
    with pytest.raises(KeyError, match="missing path mapping"):
        resolve_measurement_path(manifest, "LEU_UAA_WT_FULL", {})
    with pytest.raises(KeyError, match="unknown measurement_id"):
        resolve_measurement_path(manifest, "DOES_NOT_EXIST", {"DOES_NOT_EXIST": "x.mzML"})


def test_validate_measurement_files_checks_existence_filename_and_unicode_space_path(manifest, tmp_path):
    directory = tmp_path / "日本語 data"
    directory.mkdir()
    correct = directory / "12Old UAA.mzML"
    correct.write_bytes(b"metadata only")
    wrong = directory / "wrong.mzML"
    wrong.write_bytes(b"metadata only")
    mapping = {
        "LEU_UAA_WT_FULL": correct,
        "LEU_UAA_WT_T1": wrong,
    }
    before = dict(mapping)
    results = {item.measurement_id: item for item in validate_measurement_files(manifest, mapping)}
    assert results["LEU_UAA_WT_FULL"].exists is True
    assert results["LEU_UAA_WT_FULL"].is_file is True
    assert results["LEU_UAA_WT_FULL"].source_file_name_matches is True
    assert results["LEU_UAA_WT_FULL"].warnings == ()
    assert results["LEU_UAA_WT_T1"].source_file_name_matches is False
    assert "SOURCE_FILE_NAME_MISMATCH" in results["LEU_UAA_WT_T1"].warnings
    assert results["GLU_UUC_WT_FULL"].warnings == ("MISSING_PATH_MAPPING",)
    assert mapping == before


def test_validate_measurement_files_reports_nonexistent_path(manifest, tmp_path):
    result = validate_measurement_files(
        manifest, {"LEU_UAA_WT_FULL": tmp_path / "12Old UAA.mzML"}
    )[0]
    assert result.exists is False and result.is_file is False
    assert result.source_file_name_matches is True
    assert result.warnings == ("PATH_DOES_NOT_EXIST",)


def test_identity_metadata_and_observed_validation_are_pure(manifest):
    observed = {
        "sample_id": "LEU_UAA_WT",
        "rna_identity_id": "TRNA_LEU_UAA",
        "anticodon": "UAA",
        "experiment_type": "FULL_LENGTH",
        "source_file_name": "12Old UAA.mzML",
    }
    before = dict(observed)
    assert validate_measurement_against_observed_metadata(manifest, "LEU_UAA_WT_FULL", observed).valid is True
    observed_bad = {**observed, "anticodon": "UAG"}
    result = validate_measurement_against_observed_metadata(manifest, "LEU_UAA_WT_FULL", observed_bad)
    assert result.valid is False and result.errors == ("anticodon_mismatch:UAG!=UAA",)
    assert observed == before


def test_manifest_is_explicitly_non_formal():
    assert SAMPLE_MANIFEST_APPLIED_TO_FORMAL_SCORE is False
    assert SAMPLE_MANIFEST_APPLIED_TO_RANKING is False
    assert SAMPLE_MANIFEST_APPLIED_TO_CANDIDATE_FILTERING is False
    assert SAMPLE_MANIFEST_APPLIED_TO_FINAL_CONSENSUS is False


def test_manifest_module_has_no_config_or_import_time_load_reference():
    source = (ROOT / "rna_masshunter" / "sciex_sample_manifest.py").read_text(encoding="utf-8")
    assert "config.yaml" not in source
    assert "rna_masshunter.config" not in source


@pytest.mark.parametrize(
    ("rna_identity_id", "sequence", "anticodon", "expected_hash"),
    [
        (
            "TRNA_LEU_UAA",
            LEU_UAA_SEQUENCE,
            "UAA",
            "71664e8092c4c48c9c31fbebd57ec9db5958f4ba29bfcf424ffc5a8fce26e72d",
        ),
        (
            "TRNA_LEU_UAG",
            LEU_UAG_SEQUENCE,
            "UAG",
            "bb309562720cf5b7795103d75325791af4b9c27eaf03c4320cd48912e50d0dd6",
        ),
    ],
)
def test_user_provided_leu_sequences_are_registered_exactly(
    manifest, rna_identity_id, sequence, anticodon, expected_hash
):
    identity = get_rna_identity(manifest, rna_identity_id)
    assert identity.sequence == sequence
    assert identity.sequence_length == 85
    assert set(identity.sequence) <= set("ACGU")
    assert identity.sequence[36] == "U"
    assert identity.sequence[36:39] == anticodon
    assert identity.anticodon == anticodon
    assert identity.anticodon_orientation is AnticodonOrientation.FIVE_TO_THREE
    assert identity.wobble_sequence_index_1based == 37
    assert identity.anticodon_start_index_1based == 37
    assert identity.anticodon_end_index_1based == 39
    assert identity.sequence_sha256 == expected_hash
    assert identity.sequence_source == "USER_PROVIDED"
    assert identity.sequence_status is SequenceStatus.CONFIRMED
    assert identity.ends_with_cca is False
    assert identity.sequence[-3:] == "GCA"
    assert identity.cca_status is CCAStatus.CONFIRMED_ABSENT
    assert identity.mature_rna_status.value == "UNKNOWN"
    assert identity.intron_status.value == "UNKNOWN"
    assert "sample-level mature CCA state remains unconfirmed" in identity.sequence_notes


def test_leu_sequences_and_hashes_are_distinct(manifest):
    uaa = get_rna_identity(manifest, "TRNA_LEU_UAA")
    uag = get_rna_identity(manifest, "TRNA_LEU_UAG")
    assert uaa.sequence != uag.sequence
    assert uaa.anticodon != uag.anticodon
    assert uaa.sequence_sha256 != uag.sequence_sha256


@pytest.mark.parametrize(
    ("rna_identity_id", "full_id", "digest_id"),
    [
        ("TRNA_LEU_UAA", "LEU_UAA_WT_FULL", "LEU_UAA_WT_T1"),
        ("TRNA_LEU_UAG", "LEU_UAG_WT_FULL", "LEU_UAG_WT_T1"),
    ],
)
def test_full_and_digest_resolve_to_same_registered_sequence(
    manifest, rna_identity_id, full_id, digest_id
):
    identity = get_rna_identity(manifest, rna_identity_id)
    full = resolve_measurement_identity(manifest, full_id)
    digest = resolve_measurement_identity(manifest, digest_id)
    assert full.rna_identity_id == digest.rna_identity_id == rna_identity_id
    assert full.expected_sequence_sha256 == digest.expected_sequence_sha256 == identity.sequence_sha256


def test_measurement_filename_change_does_not_select_sequence(tmp_path, payload):
    expected = load_sciex_sample_manifest(write_payload(tmp_path, deepcopy(payload)))
    expected_hash = get_rna_identity(expected, "TRNA_LEU_UAA").sequence_sha256
    payload["measurements"][0]["source_file_name"] = "unrelated-name.mzML"
    changed = load_sciex_sample_manifest(write_payload(tmp_path, payload))
    assert get_rna_identity(changed, "TRNA_LEU_UAA").sequence_sha256 == expected_hash


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"wobble_sequence_index_1based": 0}, "positive 1-based"),
        (
            {
                "wobble_sequence_index_1based": 84,
                "anticodon_start_index_1based": 84,
                "anticodon_end_index_1based": 86,
            },
            "exceeds sequence length",
        ),
        ({"anticodon": "UA"}, "three canonical RNA bases"),
        ({"anticodon": "UAG"}, "does not match"),
        ({"anticodon_end_index_1based": 40}, "span must contain exactly three"),
        ({"wobble_sequence_index_1based": 36}, "must equal anticodon_start"),
    ],
)
def test_anticodon_sequence_index_validation_rejects_inconsistency(
    tmp_path, payload, updates, message
):
    before = deepcopy(payload)
    payload["rna_identities"][0].update(updates)
    with pytest.raises(ManifestValidationError, match=message):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))
    assert before["rna_identities"][2] == payload["rna_identities"][2]


def test_registered_sequence_requires_all_explicit_1_based_indexes(tmp_path, payload):
    payload["rna_identities"][0]["anticodon_start_index_1based"] = None
    with pytest.raises(ManifestValidationError, match="requires wobble and anticodon 1-based indexes"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


def test_unknown_sequence_rejects_orphan_sequence_indexes(tmp_path, payload):
    raw = payload["rna_identities"][1]
    raw.update({
        "sequence": None,
        "sequence_source": None,
        "sequence_status": "UNKNOWN",
        "sequence_orientation": "UNKNOWN",
        "sequence_alphabet": "UNKNOWN",
        "registered_sequence_cca_mode": "UNKNOWN",
        "registered_cca_tail_state": None,
        "cca_status": "UNKNOWN",
        "anticodon_start_index_1based": None,
        "anticodon_end_index_1based": None,
    })
    with pytest.raises(ManifestValidationError, match="sequence indexes require a registered sequence"):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


def test_confirmed_user_sequence_does_not_confirm_mature_processing_state(manifest):
    for rna_identity_id in ("TRNA_LEU_UAA", "TRNA_LEU_UAG"):
        identity = get_rna_identity(manifest, rna_identity_id)
        assert identity.sequence_status is SequenceStatus.CONFIRMED
        assert identity.mature_rna_status.value == "UNKNOWN"
        assert identity.intron_status.value == "UNKNOWN"


def test_provided_sequences_are_serialized_without_cca_or_other_editing(payload, manifest):
    for index, expected in ((0, LEU_UAA_SEQUENCE), (2, LEU_UAG_SEQUENCE)):
        assert payload["rna_identities"][index]["sequence"] == expected
        loaded = manifest.rna_identities[index].sequence
        assert loaded == expected
        assert len(loaded) == 85
        assert loaded[-3:] == "GCA"
        assert not loaded.endswith("CCA")


def test_glu_uuc_identity_uses_confirmed_existing_config_sequence(manifest):
    identity = get_rna_identity(manifest, "TRNA_GLU_UUC")
    expected = "GCUCCGGUAGUGUAGUCCGGCCAAUCAUUCCGGCCUUUCGAGCCGAAGACUCGGGUUCGAAUCCCGGCCGGAGCACCA"
    assert identity.display_name == "tRNA Glu-UUC"
    assert identity.anticodon == "UUC"
    assert identity.anticodon_orientation is AnticodonOrientation.UNCONFIRMED
    assert identity.sequence == expected
    assert identity.sequence_source == "USER_PROVIDED_AND_EXISTING_CONFIG_CONFIRMED"
    assert identity.sequence_status is SequenceStatus.CONFIRMED
    assert identity.wobble_sequence_index_1based == 37
    assert identity.anticodon_start_index_1based == 37
    assert identity.anticodon_end_index_1based == 39
    assert identity.sequence[36:39] == "UUC"
    assert identity.ends_with_cca is True


def test_sequence_registration_does_not_change_measurements_or_formal_flags(manifest):
    assert len(manifest.measurements) == 6
    assert {item.measurement_id for item in manifest.measurements} == set(EXPECTED_MEASUREMENTS)
    assert SAMPLE_MANIFEST_APPLIED_TO_FORMAL_SCORE is False
    assert SAMPLE_MANIFEST_APPLIED_TO_RANKING is False
    assert SAMPLE_MANIFEST_APPLIED_TO_CANDIDATE_FILTERING is False
    assert SAMPLE_MANIFEST_APPLIED_TO_FINAL_CONSENSUS is False


@pytest.mark.parametrize(
    ("identity_index", "updates", "message"),
    [
        (0, {"registered_cca_tail_state": "CCA"}, "EXCLUDES_CCA requires"),
        (
            1,
            {"sequence": "GCUCCGGUAGUGUAGUCCGGCCAAUCAUUCCGGCCUUUCGAGCCGAAGACUCGGGUUCGAAUCCCGGCCGGAGCACC"},
            "requires sequence ending in CCA",
        ),
        (
            2,
            {"registered_sequence_cca_mode": "UNKNOWN", "registered_cca_tail_state": "NONE"},
            "UNKNOWN registered CCA mode requires null",
        ),
    ],
)
def test_registered_cca_mode_and_state_inconsistency_is_rejected(
    tmp_path, payload, identity_index, updates, message
):
    payload["rna_identities"][identity_index].update(updates)
    with pytest.raises(ManifestValidationError, match=message):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        (
            {"sample_cca_tail_state": "CCA", "sample_cca_tail_status": "UNKNOWN"},
            "UNKNOWN sample CCA status requires null",
        ),
        (
            {"sample_cca_tail_state": None, "sample_cca_tail_status": "CONFIRMED"},
            "CONFIRMED sample CCA status requires",
        ),
    ],
)
def test_sample_cca_state_remains_separate_and_explicit(tmp_path, payload, updates, message):
    payload["samples"][0].update(updates)
    with pytest.raises(ManifestValidationError, match=message):
        load_sciex_sample_manifest(write_payload(tmp_path, payload))
