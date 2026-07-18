"""Validated SCIEX sample manifest and RNA sequence-registry models.

This module is deliberately independent from production routing and formal scoring.
Loading is explicit; importing the module never reads configuration or data files.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from hashlib import sha256
from pathlib import Path, PureWindowsPath
from typing import Any, Mapping, TypeVar

import yaml

SCHEMA_VERSION = 1
SAMPLE_MANIFEST_APPLIED_TO_FORMAL_SCORE = False
SAMPLE_MANIFEST_APPLIED_TO_RANKING = False
SAMPLE_MANIFEST_APPLIED_TO_CANDIDATE_FILTERING = False
SAMPLE_MANIFEST_APPLIED_TO_FINAL_CONSENSUS = False
_ASCII_WHITESPACE = frozenset(" \t\n\r\v\f")
_HEX_DIGITS = frozenset("0123456789abcdef")


class ManifestValidationError(ValueError):
    """Raised when a SCIEX sample manifest violates its schema."""


class _TextEnum(str, Enum):
    pass


class RNAType(_TextEnum):
    TRNA = "TRNA"


class AnticodonOrientation(_TextEnum):
    FIVE_TO_THREE = "FIVE_TO_THREE"
    THREE_TO_FIVE = "THREE_TO_FIVE"
    UNCONFIRMED = "UNCONFIRMED"


class SequenceStatus(_TextEnum):
    CONFIRMED = "CONFIRMED"
    PROVISIONAL = "PROVISIONAL"
    UNKNOWN = "UNKNOWN"


class SequenceOrientation(_TextEnum):
    FIVE_TO_THREE = "FIVE_TO_THREE"
    UNKNOWN = "UNKNOWN"


class SequenceAlphabet(_TextEnum):
    RNA = "RNA"
    UNKNOWN = "UNKNOWN"


class CCAStatus(_TextEnum):
    CONFIRMED_PRESENT = "CONFIRMED_PRESENT"
    CONFIRMED_ABSENT = "CONFIRMED_ABSENT"
    LIKELY_PRESENT = "LIKELY_PRESENT"
    LIKELY_ABSENT = "LIKELY_ABSENT"
    UNKNOWN = "UNKNOWN"


class IntronStatus(_TextEnum):
    CONFIRMED_PRESENT = "CONFIRMED_PRESENT"
    CONFIRMED_ABSENT = "CONFIRMED_ABSENT"
    UNKNOWN = "UNKNOWN"


class MatureRNAStatus(_TextEnum):
    MATURE_RNA = "MATURE_RNA"
    GENOMIC_SEQUENCE = "GENOMIC_SEQUENCE"
    UNKNOWN = "UNKNOWN"


class IdentityStatus(_TextEnum):
    CONFIRMED = "CONFIRMED"
    PROVISIONAL = "PROVISIONAL"
    UNCONFIRMED = "UNCONFIRMED"


class NativeOrTranscript(_TextEnum):
    NATIVE = "NATIVE"
    TRANSCRIPT = "TRANSCRIPT"
    UNKNOWN = "UNKNOWN"


class ExperimentType(_TextEnum):
    FULL_LENGTH = "FULL_LENGTH"
    RNASE_T1_DIGEST = "RNASE_T1_DIGEST"
    NUCLEASE_P1_AP_DIGEST = "NUCLEASE_P1_AP_DIGEST"


class InputRole(_TextEnum):
    PRIMARY = "PRIMARY"
    SUPPORTING = "SUPPORTING"


class ProfileType(_TextEnum):
    MZML = "MZML"
    UNKNOWN = "UNKNOWN"


class DigestionEnzyme(_TextEnum):
    NONE = "NONE"
    RNASE_T1 = "RNASE_T1"
    NUCLEASE_P1 = "NUCLEASE_P1"


class PhosphataseTreatment(_TextEnum):
    NONE = "NONE"
    ALKALINE_PHOSPHATASE = "ALKALINE_PHOSPHATASE"


class AnalyteLevel(_TextEnum):
    INTACT_RNA = "INTACT_RNA"
    OLIGONUCLEOTIDE = "OLIGONUCLEOTIDE"
    NUCLEOSIDE = "NUCLEOSIDE"


class TerminalStateStatus(_TextEnum):
    CONFIRMED = "CONFIRMED"
    ASSUMED = "ASSUMED"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class RNAIdentity:
    rna_identity_id: str
    display_name: str
    rna_type: RNAType
    amino_acid: str
    anticodon: str
    anticodon_orientation: AnticodonOrientation
    organism: str | None
    strain: str | None
    gene_or_locus: str | None
    sequence: str | None
    sequence_source: str | None
    sequence_status: SequenceStatus
    sequence_orientation: SequenceOrientation
    sequence_alphabet: SequenceAlphabet
    sequence_sha256: str | None
    sequence_length: int | None
    ends_with_cca: bool | None
    cca_status: CCAStatus
    intron_status: IntronStatus
    mature_rna_status: MatureRNAStatus


@dataclass(frozen=True)
class BiologicalSample:
    sample_id: str
    display_name: str
    rna_identity_id: str
    condition: str | None
    biological_source: str | None
    purification_method: str | None
    native_or_transcript: NativeOrTranscript
    sample_notes: str | None
    identity_status: IdentityStatus


@dataclass(frozen=True)
class SCIEXMeasurement:
    measurement_id: str
    sample_id: str
    experiment_type: ExperimentType
    input_role: InputRole
    input_alias: str
    source_file_name: str
    source_file_path_hint: str | None
    source_file_sha256: str | None
    profile_type: ProfileType
    digestion_enzyme: DigestionEnzyme
    phosphatase_treatment: PhosphataseTreatment
    expected_analyte_level: AnalyteLevel
    sequence_required: bool
    terminal_state_required: bool
    five_prime_state: str | None
    three_prime_state: str | None
    terminal_state_status: TerminalStateStatus
    enabled: bool
    notes: str | None


@dataclass(frozen=True)
class SCIEXSampleManifest:
    schema_version: int
    rna_identities: tuple[RNAIdentity, ...]
    samples: tuple[BiologicalSample, ...]
    measurements: tuple[SCIEXMeasurement, ...]


@dataclass(frozen=True)
class MeasurementFileValidation:
    measurement_id: str
    resolved_path: Path | None
    exists: bool
    is_file: bool
    source_file_name_matches: bool
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class MeasurementIdentityMetadata:
    measurement_id: str
    sample_id: str
    sample_family: str
    rna_identity_id: str
    expected_anticodon: str
    expected_sequence_sha256: str | None
    expected_experiment_type: ExperimentType
    expected_source_file_name: str
    identity_status: IdentityStatus


@dataclass(frozen=True)
class ObservedMetadataValidation:
    valid: bool
    errors: tuple[str, ...]


_EnumT = TypeVar("_EnumT", bound=Enum)
_IDENTITY_FIELDS = frozenset({
    "rna_identity_id", "display_name", "rna_type", "amino_acid", "anticodon",
    "anticodon_orientation", "organism", "strain", "gene_or_locus", "sequence",
    "sequence_source", "sequence_status", "sequence_orientation", "sequence_alphabet",
    "sequence_sha256", "sequence_length", "ends_with_cca", "cca_status", "intron_status",
    "mature_rna_status",
})
_SAMPLE_FIELDS = frozenset({
    "sample_id", "display_name", "rna_identity_id", "condition", "biological_source",
    "purification_method", "native_or_transcript", "sample_notes", "identity_status",
})
_MEASUREMENT_FIELDS = frozenset({
    "measurement_id", "sample_id", "experiment_type", "input_role", "input_alias",
    "source_file_name", "source_file_path_hint", "source_file_sha256", "profile_type",
    "digestion_enzyme", "phosphatase_treatment", "expected_analyte_level",
    "sequence_required", "terminal_state_required", "five_prime_state", "three_prime_state",
    "terminal_state_status", "enabled", "notes",
})
_EXPECTED_TREATMENT = {
    ExperimentType.FULL_LENGTH: (DigestionEnzyme.NONE, PhosphataseTreatment.NONE, AnalyteLevel.INTACT_RNA),
    ExperimentType.RNASE_T1_DIGEST: (
        DigestionEnzyme.RNASE_T1,
        PhosphataseTreatment.NONE,
        AnalyteLevel.OLIGONUCLEOTIDE,
    ),
    ExperimentType.NUCLEASE_P1_AP_DIGEST: (
        DigestionEnzyme.NUCLEASE_P1,
        PhosphataseTreatment.ALKALINE_PHOSPHATASE,
        AnalyteLevel.NUCLEOSIDE,
    ),
}


def _mapping(value: Any, location: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ManifestValidationError(f"{location}: expected a mapping")
    return value


def _strict_fields(raw: Mapping[str, Any], allowed: frozenset[str], location: str) -> None:
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ManifestValidationError(f"{location}: unknown field(s): {', '.join(unknown)}")


def _required_text(raw: Mapping[str, Any], field: str, location: str) -> str:
    value = raw.get(field)
    if not isinstance(value, str) or not value.strip():
        raise ManifestValidationError(f"{location}.{field}: expected non-empty text")
    return value.strip()


def _optional_text(raw: Mapping[str, Any], field: str, location: str) -> str | None:
    value = raw.get(field)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ManifestValidationError(f"{location}.{field}: expected non-empty text or null")
    return value.strip()


def _enum(raw: Mapping[str, Any], field: str, enum_type: type[_EnumT], location: str) -> _EnumT:
    value = raw.get(field)
    try:
        return enum_type(value)
    except (TypeError, ValueError) as exc:
        allowed = ", ".join(item.value for item in enum_type)
        raise ManifestValidationError(
            f"{location}.{field}: unknown value {value!r}; expected one of {allowed}"
        ) from exc


def _boolean(raw: Mapping[str, Any], field: str, location: str) -> bool:
    value = raw.get(field)
    if not isinstance(value, bool):
        raise ManifestValidationError(f"{location}.{field}: expected boolean")
    return value


def _optional_sha256(raw: Mapping[str, Any], field: str, location: str) -> str | None:
    value = _optional_text(raw, field, location)
    if value is None:
        return None
    normalized = value.lower()
    if len(normalized) != 64 or any(character not in _HEX_DIGITS for character in normalized):
        raise ManifestValidationError(f"{location}.{field}: expected a 64-character SHA-256")
    return normalized


def _normalize_sequence(value: Any, location: str) -> str:
    if not isinstance(value, str):
        raise ManifestValidationError(f"{location}.sequence: expected text or null")
    normalized = "".join(base for base in value if base not in _ASCII_WHITESPACE).upper()
    if not normalized:
        raise ManifestValidationError(f"{location}.sequence: must not be empty")
    for position, base in enumerate(normalized, 1):
        if base not in "ACGU":
            raise ManifestValidationError(
                f"{location}.sequence: invalid canonical RNA base at position {position}: {base!r}"
            )
    return normalized


def _parse_identity(raw_value: Any, index: int) -> RNAIdentity:
    location = f"rna_identities[{index}]"
    raw = _mapping(raw_value, location)
    _strict_fields(raw, _IDENTITY_FIELDS, location)
    status = _enum(raw, "sequence_status", SequenceStatus, location)
    sequence_value = raw.get("sequence")
    sequence = None if sequence_value is None else _normalize_sequence(sequence_value, location)
    if status is SequenceStatus.UNKNOWN and sequence is not None:
        raise ManifestValidationError(f"{location}: UNKNOWN sequence_status requires null sequence")
    if status is not SequenceStatus.UNKNOWN and sequence is None:
        raise ManifestValidationError(f"{location}: {status.value} sequence_status requires a sequence")

    orientation = _enum(raw, "sequence_orientation", SequenceOrientation, location)
    alphabet = _enum(raw, "sequence_alphabet", SequenceAlphabet, location)
    sequence_source = _optional_text(raw, "sequence_source", location)
    cca_status = _enum(raw, "cca_status", CCAStatus, location)
    mature_status = _enum(raw, "mature_rna_status", MatureRNAStatus, location)
    if status is SequenceStatus.CONFIRMED:
        if sequence_source is None:
            raise ManifestValidationError(f"{location}: CONFIRMED sequence requires sequence_source")
        if orientation is not SequenceOrientation.FIVE_TO_THREE:
            raise ManifestValidationError(f"{location}: CONFIRMED sequence requires FIVE_TO_THREE orientation")
        if alphabet is not SequenceAlphabet.RNA:
            raise ManifestValidationError(f"{location}: CONFIRMED sequence requires RNA alphabet")
        if mature_status is MatureRNAStatus.UNKNOWN:
            raise ManifestValidationError(f"{location}: CONFIRMED sequence requires mature/genomic status")
        if cca_status not in {CCAStatus.CONFIRMED_PRESENT, CCAStatus.CONFIRMED_ABSENT}:
            raise ManifestValidationError(f"{location}: CONFIRMED sequence requires confirmed CCA status")

    computed_hash = sha256(sequence.encode("ascii")).hexdigest() if sequence is not None else None
    computed_length = len(sequence) if sequence is not None else None
    computed_cca = sequence.endswith("CCA") if sequence is not None else None
    supplied_hash = _optional_sha256(raw, "sequence_sha256", location)
    supplied_length = raw.get("sequence_length")
    supplied_cca = raw.get("ends_with_cca")
    if supplied_hash is not None and supplied_hash != computed_hash:
        raise ManifestValidationError(f"{location}.sequence_sha256: does not match normalized sequence")
    if supplied_length is not None and (not isinstance(supplied_length, int) or isinstance(supplied_length, bool)):
        raise ManifestValidationError(f"{location}.sequence_length: expected integer or null")
    if supplied_length is not None and supplied_length != computed_length:
        raise ManifestValidationError(f"{location}.sequence_length: does not match normalized sequence")
    if supplied_cca is not None and not isinstance(supplied_cca, bool):
        raise ManifestValidationError(f"{location}.ends_with_cca: expected boolean or null")
    if supplied_cca is not None and supplied_cca != computed_cca:
        raise ManifestValidationError(f"{location}.ends_with_cca: does not match normalized sequence")
    if computed_cca is True and cca_status is CCAStatus.CONFIRMED_ABSENT:
        raise ManifestValidationError(f"{location}.cca_status: sequence ends with CCA")
    if computed_cca is False and cca_status is CCAStatus.CONFIRMED_PRESENT:
        raise ManifestValidationError(f"{location}.cca_status: sequence does not end with CCA")

    anticodon = _required_text(raw, "anticodon", location).upper()
    if len(anticodon) != 3 or any(base not in "ACGU" for base in anticodon):
        raise ManifestValidationError(f"{location}.anticodon: expected three canonical RNA bases")
    return RNAIdentity(
        rna_identity_id=_required_text(raw, "rna_identity_id", location),
        display_name=_required_text(raw, "display_name", location),
        rna_type=_enum(raw, "rna_type", RNAType, location),
        amino_acid=_required_text(raw, "amino_acid", location),
        anticodon=anticodon,
        anticodon_orientation=_enum(raw, "anticodon_orientation", AnticodonOrientation, location),
        organism=_optional_text(raw, "organism", location),
        strain=_optional_text(raw, "strain", location),
        gene_or_locus=_optional_text(raw, "gene_or_locus", location),
        sequence=sequence,
        sequence_source=sequence_source,
        sequence_status=status,
        sequence_orientation=orientation,
        sequence_alphabet=alphabet,
        sequence_sha256=computed_hash,
        sequence_length=computed_length,
        ends_with_cca=computed_cca,
        cca_status=cca_status,
        intron_status=_enum(raw, "intron_status", IntronStatus, location),
        mature_rna_status=mature_status,
    )


def _parse_sample(raw_value: Any, index: int) -> BiologicalSample:
    location = f"samples[{index}]"
    raw = _mapping(raw_value, location)
    _strict_fields(raw, _SAMPLE_FIELDS, location)
    return BiologicalSample(
        sample_id=_required_text(raw, "sample_id", location),
        display_name=_required_text(raw, "display_name", location),
        rna_identity_id=_required_text(raw, "rna_identity_id", location),
        condition=_optional_text(raw, "condition", location),
        biological_source=_optional_text(raw, "biological_source", location),
        purification_method=_optional_text(raw, "purification_method", location),
        native_or_transcript=_enum(raw, "native_or_transcript", NativeOrTranscript, location),
        sample_notes=_optional_text(raw, "sample_notes", location),
        identity_status=_enum(raw, "identity_status", IdentityStatus, location),
    )


def _parse_measurement(raw_value: Any, index: int) -> SCIEXMeasurement:
    location = f"measurements[{index}]"
    raw = _mapping(raw_value, location)
    _strict_fields(raw, _MEASUREMENT_FIELDS, location)
    experiment = _enum(raw, "experiment_type", ExperimentType, location)
    digestion = _enum(raw, "digestion_enzyme", DigestionEnzyme, location)
    phosphatase = _enum(raw, "phosphatase_treatment", PhosphataseTreatment, location)
    analyte = _enum(raw, "expected_analyte_level", AnalyteLevel, location)
    expected = _EXPECTED_TREATMENT[experiment]
    if (digestion, phosphatase, analyte) != expected:
        raise ManifestValidationError(
            f"{location}: {experiment.value} requires digestion={expected[0].value}, "
            f"phosphatase={expected[1].value}, analyte={expected[2].value}"
        )
    file_name = _required_text(raw, "source_file_name", location)
    if (
        Path(file_name).name != file_name
        or PureWindowsPath(file_name).name != file_name
        or Path(file_name).is_absolute()
        or PureWindowsPath(file_name).is_absolute()
    ):
        raise ManifestValidationError(
            f"{location}.source_file_name: filename only; absolute/path values are forbidden"
        )
    path_hint = _optional_text(raw, "source_file_path_hint", location)
    if path_hint is not None and (
        Path(path_hint).is_absolute() or PureWindowsPath(path_hint).is_absolute()
    ):
        raise ManifestValidationError(f"{location}.source_file_path_hint: absolute paths are forbidden")
    terminal_status = _enum(raw, "terminal_state_status", TerminalStateStatus, location)
    five_prime = _optional_text(raw, "five_prime_state", location)
    three_prime = _optional_text(raw, "three_prime_state", location)
    terminal_required = _boolean(raw, "terminal_state_required", location)
    if experiment is not ExperimentType.FULL_LENGTH:
        has_terminal_metadata = (
            terminal_required
            or five_prime is not None
            or three_prime is not None
            or terminal_status is not TerminalStateStatus.UNKNOWN
        )
        if has_terminal_metadata:
            raise ManifestValidationError(
                f"{location}: digest measurements cannot reuse intact terminal-state metadata"
            )
    if terminal_status is TerminalStateStatus.CONFIRMED and (five_prime is None or three_prime is None):
        raise ManifestValidationError(f"{location}: CONFIRMED terminal state requires both terminal values")
    return SCIEXMeasurement(
        measurement_id=_required_text(raw, "measurement_id", location),
        sample_id=_required_text(raw, "sample_id", location),
        experiment_type=experiment,
        input_role=_enum(raw, "input_role", InputRole, location),
        input_alias=_required_text(raw, "input_alias", location),
        source_file_name=file_name,
        source_file_path_hint=path_hint,
        source_file_sha256=_optional_sha256(raw, "source_file_sha256", location),
        profile_type=_enum(raw, "profile_type", ProfileType, location),
        digestion_enzyme=digestion,
        phosphatase_treatment=phosphatase,
        expected_analyte_level=analyte,
        sequence_required=_boolean(raw, "sequence_required", location),
        terminal_state_required=terminal_required,
        five_prime_state=five_prime,
        three_prime_state=three_prime,
        terminal_state_status=terminal_status,
        enabled=_boolean(raw, "enabled", location),
        notes=_optional_text(raw, "notes", location),
    )


def _unique(items: tuple[Any, ...], field: str, label: str) -> None:
    seen: set[str] = set()
    for item in items:
        value = getattr(item, field)
        if value in seen:
            raise ManifestValidationError(f"duplicate {label}: {value}")
        seen.add(value)


def load_sciex_sample_manifest(path: str | Path) -> SCIEXSampleManifest:
    """Load and fully validate a manifest using PyYAML's safe loader."""
    source = Path(path)
    if not source.is_file():
        raise ManifestValidationError(f"manifest not found: {source}")
    try:
        payload = yaml.safe_load(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise ManifestValidationError(f"cannot load manifest {source}: {exc}") from exc
    root = _mapping(payload, "manifest")
    _strict_fields(root, frozenset({"schema_version", "rna_identities", "samples", "measurements"}), "manifest")
    if root.get("schema_version") != SCHEMA_VERSION:
        raise ManifestValidationError(f"unsupported schema_version: expected {SCHEMA_VERSION}")
    for field in ("rna_identities", "samples", "measurements"):
        if not isinstance(root.get(field), list) or not root[field]:
            raise ManifestValidationError(f"manifest.{field}: expected a non-empty list")
    identities = tuple(_parse_identity(value, index) for index, value in enumerate(root["rna_identities"], 1))
    samples = tuple(_parse_sample(value, index) for index, value in enumerate(root["samples"], 1))
    measurements = tuple(_parse_measurement(value, index) for index, value in enumerate(root["measurements"], 1))
    _unique(identities, "rna_identity_id", "rna_identity_id")
    _unique(samples, "sample_id", "sample_id")
    _unique(measurements, "measurement_id", "measurement_id")
    _unique(measurements, "input_alias", "input_alias")
    identity_ids = {item.rna_identity_id for item in identities}
    for sample in samples:
        if sample.rna_identity_id not in identity_ids:
            raise ManifestValidationError(
                f"sample {sample.sample_id}: unknown rna_identity_id {sample.rna_identity_id}"
            )
    sample_ids = {item.sample_id for item in samples}
    for measurement in measurements:
        if measurement.sample_id not in sample_ids:
            raise ManifestValidationError(
                f"measurement {measurement.measurement_id}: unknown sample_id {measurement.sample_id}"
            )
    return SCIEXSampleManifest(SCHEMA_VERSION, identities, samples, measurements)


def get_measurement(manifest: SCIEXSampleManifest, measurement_id: str) -> SCIEXMeasurement:
    for measurement in manifest.measurements:
        if measurement.measurement_id == measurement_id:
            return measurement
    raise KeyError(f"unknown measurement_id: {measurement_id}")


def get_sample(manifest: SCIEXSampleManifest, sample_id: str) -> BiologicalSample:
    for sample in manifest.samples:
        if sample.sample_id == sample_id:
            return sample
    raise KeyError(f"unknown sample_id: {sample_id}")


def get_rna_identity(manifest: SCIEXSampleManifest, rna_identity_id: str) -> RNAIdentity:
    for identity in manifest.rna_identities:
        if identity.rna_identity_id == rna_identity_id:
            return identity
    raise KeyError(f"unknown rna_identity_id: {rna_identity_id}")


def get_measurements_for_sample(manifest: SCIEXSampleManifest, sample_id: str) -> tuple[SCIEXMeasurement, ...]:
    get_sample(manifest, sample_id)
    return tuple(item for item in manifest.measurements if item.sample_id == sample_id)


def get_samples_for_rna_identity(manifest: SCIEXSampleManifest, rna_identity_id: str) -> tuple[BiologicalSample, ...]:
    get_rna_identity(manifest, rna_identity_id)
    return tuple(item for item in manifest.samples if item.rna_identity_id == rna_identity_id)


def filter_measurements_by_experiment_type(
    manifest: SCIEXSampleManifest, experiment_type: ExperimentType | str
) -> tuple[SCIEXMeasurement, ...]:
    try:
        selected = experiment_type if isinstance(experiment_type, ExperimentType) else ExperimentType(experiment_type)
    except ValueError as exc:
        raise KeyError(f"unknown experiment_type: {experiment_type}") from exc
    return tuple(item for item in manifest.measurements if item.experiment_type is selected)


def get_paired_measurements(
    manifest: SCIEXSampleManifest, measurement_id: str
) -> tuple[SCIEXMeasurement, ...]:
    """Return opposite-role measurements from the same biological sample only."""
    selected = get_measurement(manifest, measurement_id)
    same_sample = get_measurements_for_sample(manifest, selected.sample_id)
    if selected.experiment_type is ExperimentType.FULL_LENGTH:
        return tuple(item for item in same_sample if item.experiment_type is not ExperimentType.FULL_LENGTH)
    return tuple(item for item in same_sample if item.experiment_type is ExperimentType.FULL_LENGTH)


def resolve_measurement_path(
    manifest: SCIEXSampleManifest,
    measurement_id: str,
    path_mapping: Mapping[str, str | Path],
) -> Path:
    """Resolve a runtime mapping without opening a file or mutating either input."""
    get_measurement(manifest, measurement_id)
    if measurement_id not in path_mapping:
        raise KeyError(f"missing path mapping for measurement_id: {measurement_id}")
    value = path_mapping[measurement_id]
    if not isinstance(value, (str, Path)):
        raise TypeError(f"path mapping for {measurement_id} must be str or Path")
    return Path(value)


def validate_measurement_files(
    manifest: SCIEXSampleManifest,
    path_mapping: Mapping[str, str | Path],
) -> tuple[MeasurementFileValidation, ...]:
    """Perform read-only path/filename checks; file contents are never opened."""
    results: list[MeasurementFileValidation] = []
    for measurement in manifest.measurements:
        warnings: list[str] = []
        if measurement.measurement_id not in path_mapping:
            results.append(MeasurementFileValidation(
                measurement.measurement_id, None, False, False, False, ("MISSING_PATH_MAPPING",)
            ))
            continue
        path = resolve_measurement_path(manifest, measurement.measurement_id, path_mapping)
        exists = path.exists()
        is_file = path.is_file()
        filename_matches = path.name == measurement.source_file_name
        if not exists:
            warnings.append("PATH_DOES_NOT_EXIST")
        elif not is_file:
            warnings.append("PATH_IS_NOT_FILE")
        if not filename_matches:
            warnings.append("SOURCE_FILE_NAME_MISMATCH")
        results.append(MeasurementFileValidation(
            measurement.measurement_id, path, exists, is_file, filename_matches, tuple(warnings)
        ))
    return tuple(results)


def resolve_measurement_identity(
    manifest: SCIEXSampleManifest, measurement_id: str
) -> MeasurementIdentityMetadata:
    measurement = get_measurement(manifest, measurement_id)
    sample = get_sample(manifest, measurement.sample_id)
    identity = get_rna_identity(manifest, sample.rna_identity_id)
    return MeasurementIdentityMetadata(
        measurement_id=measurement.measurement_id,
        sample_id=sample.sample_id,
        sample_family=sample.sample_id,
        rna_identity_id=identity.rna_identity_id,
        expected_anticodon=identity.anticodon,
        expected_sequence_sha256=identity.sequence_sha256,
        expected_experiment_type=measurement.experiment_type,
        expected_source_file_name=measurement.source_file_name,
        identity_status=sample.identity_status,
    )


def validate_measurement_against_observed_metadata(
    manifest: SCIEXSampleManifest,
    measurement_id: str,
    observed_metadata: Mapping[str, Any],
) -> ObservedMetadataValidation:
    """Compare explicit observed labels without promoting identity confidence."""
    expected = resolve_measurement_identity(manifest, measurement_id)
    comparisons = {
        "rna_identity_id": expected.rna_identity_id,
        "anticodon": expected.expected_anticodon,
        "experiment_type": expected.expected_experiment_type.value,
        "source_file_name": expected.expected_source_file_name,
        "sample_id": expected.sample_id,
    }
    errors = tuple(
        f"{field}_mismatch:{observed_metadata[field]}!={value}"
        for field, value in comparisons.items()
        if field in observed_metadata and str(observed_metadata[field]) != str(value)
    )
    return ObservedMetadataValidation(not errors, errors)
