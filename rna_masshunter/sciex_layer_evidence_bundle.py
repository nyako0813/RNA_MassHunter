"""Safe JSON bundles for production SCIEX layer-evidence audit results.

A layer-evidence bundle moves an already-produced shadow audit result between
processes without re-reading raw mzML.  Raw peak arrays, complete spectrum arrays,
and binary payloads are deliberately excluded: cross-layer reconciliation needs
summaries and evidence records, not a second copy of raw observations.

Layer non-empty contracts are intentionally strict.  FULL needs relations plus a
state series or an eligible state relation; T1 needs a fragment match or valid
state evidence (zero state families alone is not absence); P1AP_MS1 needs an MS1
match; P1AP_MS2 needs a product match or substantive identity evidence.  Candidate
or precursor lists and identity-audit strings alone never satisfy the contract.

RNA identity (name, sequence, anticodon, taxonomy) and sample identity are separate
provenance axes.  P1/AP MS1 and MS2 from one raw run share both source and
independence groups and therefore are not independent support.  Bundles preserve
explicit shadow-only safeguards and cannot propagate into formal score, ranking,
candidate filtering, localization, identity assignment, or final consensus.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import MISSING, dataclass, fields, is_dataclass
from datetime import datetime, timezone
from enum import Enum
from hashlib import sha256
from importlib import import_module
from functools import lru_cache
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping, get_type_hints

SCHEMA_VERSION = "1.0"
SERIALIZER_FORMAT_VERSION = "enum-tagged-v2"
BUNDLE_TYPE = "RNA_MASSHUNTER_LAYER_EVIDENCE"

SAFEGUARDS = {
    "formal_propagation": False,
    "chemical_identity_assigned": False,
    "modification_assigned": False,
    "exact_candidate_identity_confirmed": False,
    "exact_isomer_identity_confirmed": False,
    "exact_nucleotide_localization": False,
    "exact_atom_localization": False,
    "reaction_order_assigned": False,
    "applied_to_formal_score": False,
    "applied_to_ranking": False,
    "applied_to_candidate_filtering": False,
    "applied_to_final_consensus": False,
    "shadow_analysis_only": True,
    "cross_layer_reconciliation_only": True,
}

_REQUIRED_TOP_LEVEL = {
    "schema_version", "serializer_format_version", "bundle_type", "layer", "optional_result_key",
    "producer_name", "producer_module", "producer_commit", "created_at_utc",
    "source", "rna", "experiment", "result", "safeguards", "validation",
    "canonical_payload_sha256",
}
_REQUIRED_SOURCE = {
    "path", "basename", "size_bytes", "mtime_ns", "sha256", "run_id",
    "sample_id", "biological_sample_id",
}
_REQUIRED_RNA = {
    "name", "sequence", "anticodon", "wobble_position", "organism_group", "species",
}
_REQUIRED_EXPERIMENT = {
    "condition_name", "digest_type", "layer", "independence_group", "shared_source_group",
}
_REQUIRED_VALIDATION = {
    "status", "non_empty", "record_count", "required_fields_present",
    "provenance_verified", "safeguards_verified", "validation_messages",
}

_RAW_FIELD_NAMES = {
    "raw_peaks", "raw_peak_arrays", "spectrum_arrays", "raw_spectrum_arrays",
    "mz_array", "m/z array", "intensity_array", "intensity array", "binary", "base64",
    "binary_payload", "binary_data_array",
}

# These fields are derived in-memory profiles, not cross-layer evidence.  They
# are omitted by exact production class and field name; every other NumPy array
# remains forbidden.  Missing fields use the dataclass default, except for the
# existing ProcessedMS2Spectrum.peaks contract, whose constructor requires ().
TRANSIENT_DATACLASS_FIELDS = {
    "rna_masshunter.sciex_p1ap_nucleoside_ms2_identity_audit.ProcessedMS2Spectrum": frozenset({"peaks"}),
    "rna_masshunter.sciex_t1_replicate_consistency_audit.ReplicateRunPeakProfile": frozenset({
        "comparison_mz_grid", "comparison_raw_profile", "comparison_normalized_profile",
    }),
}
_TRANSIENT_RESTORE_VALUES = {
    "rna_masshunter.sciex_p1ap_nucleoside_ms2_identity_audit.ProcessedMS2Spectrum": {"peaks": ()},
}

_ALLOWED_TYPE_MODULES = {
    "rna_masshunter.intact_rna_average_mass",
    "rna_masshunter.sciex_intact_peak_family",
    "rna_masshunter.sciex_intact_oxygen_water_state_audit",
    "rna_masshunter.sciex_mzml_source_metadata_audit",
    "rna_masshunter.sciex_t1_fragment_delta_audit",
    "rna_masshunter.sciex_t1_fragment_shadow_match",
    "rna_masshunter.sciex_t1_profile_peak_audit",
    "rna_masshunter.sciex_t1_replicate_consistency_audit",
    "rna_masshunter.sciex_t1_fragment_state_series_audit",
    "rna_masshunter.sciex_p1ap_nucleoside_state_audit",
    "rna_masshunter.sciex_p1ap_nucleoside_ms2_identity_audit",
}


class LayerEvidenceBundleError(ValueError):
    """Raised when a layer-evidence bundle violates its production contract."""


@dataclass(frozen=True)
class LayerContract:
    layer: str
    optional_result_key: str
    result_module: str
    result_class_name: str
    producer_module: str
    producer_name: str
    required_result_fields: tuple[str, ...]

    @property
    def result_class(self) -> type:
        return getattr(import_module(self.result_module), self.result_class_name)


LAYER_CONTRACTS: dict[str, LayerContract] = {
    "FULL": LayerContract(
        "FULL", "sciex_intact_oxygen_water_state_audit",
        "rna_masshunter.sciex_intact_oxygen_water_state_audit", "OxygenWaterStateAuditResult",
        "rna_masshunter.sciex_intact_oxygen_water_state_audit", "audit_oxygen_water_state_series",
        ("source_id", "status", "references", "relations", "edges", "series", "algorithm_version"),
    ),
    "T1": LayerContract(
        "T1", "sciex_t1_fragment_state_series_audit",
        "rna_masshunter.sciex_t1_fragment_state_series_audit", "T1FragmentStateSeriesAuditResult",
        "rna_masshunter.sciex_t1_fragment_state_series_audit", "audit_t1_fragment_state_series",
        ("parameters", "run_profile", "run_summary", "fragment_matches", "state_families", "summary"),
    ),
    "P1AP_MS1": LayerContract(
        "P1AP_MS1", "sciex_p1ap_nucleoside_state_audit",
        "rna_masshunter.sciex_p1ap_nucleoside_state_audit", "P1APNucleosideStateAuditResult",
        "rna_masshunter.sciex_p1ap_nucleoside_state_audit", "audit_p1ap_nucleoside_state_series",
        ("parameters", "run_profile", "run_summary", "candidates", "matches", "state_families", "summary"),
    ),
    "P1AP_MS2": LayerContract(
        "P1AP_MS2", "sciex_p1ap_nucleoside_ms2_identity_audit",
        "rna_masshunter.sciex_p1ap_nucleoside_ms2_identity_audit", "P1APNucleosideMS2AuditResult",
        "rna_masshunter.sciex_p1ap_nucleoside_ms2_identity_audit", "audit_p1ap_nucleoside_ms2_identity",
        ("candidates", "precursor_records", "product_match_records", "candidate_summary_records", "summary"),
    ),
}
_CONTRACTS_BY_KEY = {contract.optional_result_key: contract for contract in LAYER_CONTRACTS.values()}


def compute_file_sha256(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA256 of a file without loading it all into memory."""
    digest = sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    """Serialize JSON deterministically, rejecting NaN and Infinity."""
    try:
        text = json.dumps(
            value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise LayerEvidenceBundleError(f"value is not canonical JSON: {exc}") from exc
    return text.encode("utf-8")


def _type_name(value: type) -> str:
    return f"{value.__module__}.{value.__qualname__}"


@lru_cache(maxsize=1)
def _allowed_types() -> dict[str, type]:
    output: dict[str, type] = {}
    for module_name in sorted(_ALLOWED_TYPE_MODULES):
        module = import_module(module_name)
        for value in vars(module).values():
            if not isinstance(value, type) or value.__module__ != module_name:
                continue
            if is_dataclass(value) or issubclass(value, Enum):
                output[_type_name(value)] = value
    return output


def _dangerous_field(name: Any) -> bool:
    text = str(name).strip().lower()
    return text in _RAW_FIELD_NAMES or text.endswith("spectrum_arrays") or text.endswith("binary_payload")


def _sort_encoded(values: list[Any]) -> list[Any]:
    return sorted(values, key=canonical_json_bytes)


def _reject_forbidden_payload(value: Any) -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            if _dangerous_field(key):
                raise LayerEvidenceBundleError(f"raw spectrum/binary field is forbidden: {key}")
            _reject_forbidden_payload(item)
    elif isinstance(value, list):
        for item in value:
            _reject_forbidden_payload(item)
    elif isinstance(value, float) and not math.isfinite(value):
        raise LayerEvidenceBundleError("non-finite float is forbidden")


def _encode(value: Any, *, field_name: str | None = None) -> Any:
    if field_name and _dangerous_field(field_name):
        raise LayerEvidenceBundleError(f"raw spectrum/binary field is forbidden: {field_name}")
    if value is None:
        return None
    if isinstance(value, Enum):
        class_name = _type_name(type(value))
        if class_name not in _allowed_types():
            raise LayerEvidenceBundleError(f"unapproved enum type: {class_name}")
        return {"__bundle_type__": "enum", "class": class_name, "value": _encode(value.value)}
    if isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise LayerEvidenceBundleError("non-finite float is forbidden")
        return value
    if is_dataclass(value) and not isinstance(value, type):
        class_name = _type_name(type(value))
        if class_name not in _allowed_types():
            raise LayerEvidenceBundleError(f"unapproved dataclass type: {class_name}")
        encoded_fields: dict[str, Any] = {}
        transient_fields = TRANSIENT_DATACLASS_FIELDS.get(class_name, frozenset())
        for item in sorted(fields(value), key=lambda row: row.name):
            if item.name in transient_fields:
                continue
            encoded_fields[item.name] = _encode(getattr(value, item.name), field_name=item.name)
        return {"__bundle_type__": "dataclass", "class": class_name, "fields": encoded_fields}
    if isinstance(value, Path):
        return {"__bundle_type__": "path", "value": str(value)}
    if isinstance(value, datetime):
        if value.tzinfo is None:
            raise LayerEvidenceBundleError("naive datetime is forbidden")
        return {"__bundle_type__": "datetime", "value": value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")}
    if isinstance(value, tuple):
        return {"__bundle_type__": "tuple", "items": [_encode(item) for item in value]}
    if isinstance(value, list):
        return {"__bundle_type__": "list", "items": [_encode(item) for item in value]}
    if isinstance(value, (set, frozenset)):
        items = _sort_encoded([_encode(item) for item in value])
        return {"__bundle_type__": "frozenset" if isinstance(value, frozenset) else "set", "items": items}
    if isinstance(value, Mapping):
        items = []
        for key, item in value.items():
            if _dangerous_field(key):
                raise LayerEvidenceBundleError(f"raw spectrum/binary field is forbidden: {key}")
            items.append([_encode(key), _encode(item, field_name=str(key))])
        items.sort(key=lambda pair: canonical_json_bytes(pair[0]))
        return {"__bundle_type__": "dict", "items": items}
    if type(value).__module__ == "numpy" and hasattr(value, "item"):
        if type(value).__name__ == "ndarray":
            raise LayerEvidenceBundleError("numpy/raw arrays are forbidden")
        return _encode(value.item())
    if isinstance(value, (bytes, bytearray, memoryview)):
        raise LayerEvidenceBundleError("binary payload is forbidden")
    raise LayerEvidenceBundleError(f"unsupported value type: {_type_name(type(value))}")


def _decode(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise LayerEvidenceBundleError("non-finite float is forbidden")
        return value
    if isinstance(value, list):
        return [_decode(item) for item in value]
    if not isinstance(value, dict):
        raise LayerEvidenceBundleError(f"unsupported JSON node: {type(value).__name__}")
    kind = value.get("__bundle_type__")
    if kind is None:
        return {key: _decode(item) for key, item in value.items()}
    if kind in {"tuple", "list", "set", "frozenset"}:
        if set(value) != {"__bundle_type__", "items"} or not isinstance(value["items"], list):
            raise LayerEvidenceBundleError(f"malformed {kind} encoding")
        items = [_decode(item) for item in value["items"]]
        return tuple(items) if kind == "tuple" else items if kind == "list" else set(items) if kind == "set" else frozenset(items)
    if kind == "dict":
        if set(value) != {"__bundle_type__", "items"} or not isinstance(value["items"], list):
            raise LayerEvidenceBundleError("malformed dict encoding")
        output = {}
        for pair in value["items"]:
            if not isinstance(pair, list) or len(pair) != 2:
                raise LayerEvidenceBundleError("malformed dict item")
            key = _decode(pair[0])
            if key in output:
                raise LayerEvidenceBundleError("duplicate decoded dict key")
            output[key] = _decode(pair[1])
        return output
    if kind == "path":
        if set(value) != {"__bundle_type__", "value"} or not isinstance(value["value"], str):
            raise LayerEvidenceBundleError("malformed path encoding")
        return Path(value["value"])
    if kind == "datetime":
        if set(value) != {"__bundle_type__", "value"} or not isinstance(value["value"], str):
            raise LayerEvidenceBundleError("malformed datetime encoding")
        try:
            parsed = datetime.fromisoformat(value["value"].replace("Z", "+00:00"))
        except ValueError as exc:
            raise LayerEvidenceBundleError("malformed datetime") from exc
        if parsed.tzinfo is None:
            raise LayerEvidenceBundleError("datetime lacks timezone")
        return parsed.astimezone(timezone.utc)
    if kind in {"enum", "dataclass"}:
        class_name = value.get("class")
        approved = _allowed_types()
        cls = approved.get(class_name)
        if cls is None:
            raise LayerEvidenceBundleError(f"unapproved encoded class: {class_name}")
        if kind == "enum":
            if set(value) != {"__bundle_type__", "class", "value"} or not issubclass(cls, Enum):
                raise LayerEvidenceBundleError("malformed enum encoding")
            try:
                return cls(_decode(value["value"]))
            except (TypeError, ValueError) as exc:
                raise LayerEvidenceBundleError(f"invalid enum value for {class_name}") from exc
        if set(value) != {"__bundle_type__", "class", "fields"} or not is_dataclass(cls) or not isinstance(value["fields"], dict):
            raise LayerEvidenceBundleError("malformed dataclass encoding")
        definitions = {item.name: item for item in fields(cls)}
        unknown = set(value["fields"]) - set(definitions)
        if unknown:
            raise LayerEvidenceBundleError(f"unknown fields for {class_name}: {sorted(unknown)}")
        serialized_transient = set(value["fields"]) & set(TRANSIENT_DATACLASS_FIELDS.get(class_name, ()))
        if serialized_transient:
            raise LayerEvidenceBundleError(
                f"transient fields must not be serialized for {class_name}: {sorted(serialized_transient)}"
            )
        kwargs = {}
        type_hints = get_type_hints(cls)
        for name, definition in definitions.items():
            if not definition.init:
                continue
            if name in value["fields"]:
                decoded = _decode(value["fields"][name])
                expected_type = type_hints.get(name)
                if (
                    isinstance(expected_type, type)
                    and issubclass(expected_type, Enum)
                    and type(decoded) is not expected_type
                ):
                    raise LayerEvidenceBundleError(
                        f"enum field type mismatch for {class_name}.{name}: "
                        f"expected {_type_name(expected_type)}, got {_type_name(type(decoded))}"
                    )
                kwargs[name] = decoded
            elif name in _TRANSIENT_RESTORE_VALUES.get(class_name, {}):
                kwargs[name] = _TRANSIENT_RESTORE_VALUES[class_name][name]
            elif definition.default is MISSING and definition.default_factory is MISSING:
                raise LayerEvidenceBundleError(f"missing field for {class_name}: {name}")
        try:
            return cls(**kwargs)
        except (TypeError, ValueError) as exc:
            raise LayerEvidenceBundleError(f"cannot restore {class_name}: {exc}") from exc
    raise LayerEvidenceBundleError(f"unknown encoded type: {kind}")


def _require_mapping(value: Any, name: str, keys: set[str]) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise LayerEvidenceBundleError(f"{name} must be an object")
    missing = keys - set(value)
    if missing:
        raise LayerEvidenceBundleError(f"{name} missing required fields: {sorted(missing)}")
    return value


def _contract_for_bundle(bundle: Mapping[str, Any]) -> LayerContract:
    layer = bundle.get("layer")
    contract = LAYER_CONTRACTS.get(layer)
    if contract is None:
        raise LayerEvidenceBundleError(f"unknown layer: {layer}")
    key = bundle.get("optional_result_key")
    if key not in _CONTRACTS_BY_KEY or key != contract.optional_result_key:
        raise LayerEvidenceBundleError(f"unknown or incompatible optional_result_key: {key}")
    if bundle.get("producer_module") != contract.producer_module or bundle.get("producer_name") != contract.producer_name:
        raise LayerEvidenceBundleError("producer does not match layer contract")
    return contract


def _required_fields_present(result: Any, contract: LayerContract) -> bool:
    return all(hasattr(result, name) for name in contract.required_result_fields)


def _substantive_ms2_identity(record: Any) -> bool:
    status = str(getattr(record, "ms2_identity_evidence_status", "") or "").upper()
    excluded = {"", "MS2_PRECURSOR_COMPATIBLE_ONLY", "MS2_INSUFFICIENT", "NO_MS2_EVIDENCE"}
    return status not in excluded


def _record_metrics(result: Any, contract: LayerContract) -> tuple[int, bool]:
    if contract.layer == "FULL":
        references = tuple(result.references or ())
        relations = tuple(result.relations or ())
        edges = tuple(result.edges or ())
        series = tuple(result.series or ())
        eligible_relation = any(bool(getattr(row, "eligible_for_state_series", False)) for row in relations)
        meaningful_series = any(
            int(getattr(row, "member_count", 0) or 0) > 1
            and (int(getattr(row, "o_equivalent_edge_count", 0) or 0)
                 + int(getattr(row, "h2o_equivalent_edge_count", 0) or 0)) > 0
            for row in series
        )
        return len(references) + len(relations) + len(edges) + len(series), bool(relations and (meaningful_series or eligible_relation))
    if contract.layer == "T1":
        matches = tuple(result.fragment_matches or ())
        families = tuple(result.state_families or ())
        reconciliations = tuple(getattr(result, "reconciliations", ()) or ())
        valid_state = bool(families) or any(
            str(getattr(getattr(row, "reconciliation_status", ""), "value", getattr(row, "reconciliation_status", "")))
            not in {"", "T1_SERIES_NOT_OBSERVED", "INSUFFICIENT_T1_EVIDENCE"}
            for row in reconciliations
        )
        return len(matches) + len(families), bool(matches or valid_state)
    if contract.layer == "P1AP_MS1":
        matches = tuple(result.matches or ())
        families = tuple(result.state_families or ())
        return len(matches) + len(families), bool(matches)
    matches = tuple(result.product_match_records or ())
    summaries = tuple(result.candidate_summary_records or ())
    substantive = sum(_substantive_ms2_identity(row) for row in summaries)
    return len(matches) + len(summaries), bool(matches or substantive)


def _validate_safeguards(result: Any, encoded_result: Any, top: Mapping[str, Any]) -> None:
    if set(top) != set(SAFEGUARDS):
        missing = set(SAFEGUARDS) - set(top)
        extra = set(top) - set(SAFEGUARDS)
        raise LayerEvidenceBundleError(f"safeguard fields mismatch; missing={sorted(missing)} extra={sorted(extra)}")
    for name, expected in SAFEGUARDS.items():
        if top[name] is not expected:
            raise LayerEvidenceBundleError(f"unsafe safeguard value: {name}={top[name]!r}")

    def visit(value: Any) -> None:
        if isinstance(value, dict):
            if value.get("__bundle_type__") == "dataclass" and isinstance(value.get("fields"), dict):
                for name, expected in SAFEGUARDS.items():
                    if name in value["fields"] and value["fields"][name] is not expected:
                        raise LayerEvidenceBundleError(f"result safeguard conflicts: {name}")
            for item in value.values():
                visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)
    visit(encoded_result)
    def visit_restored(value: Any) -> None:
        if is_dataclass(value) and not isinstance(value, type):
            for definition in fields(value):
                item = getattr(value, definition.name)
                if definition.name in SAFEGUARDS and item is not SAFEGUARDS[definition.name]:
                    raise LayerEvidenceBundleError(f"restored result safeguard conflicts: {definition.name}")
                visit_restored(item)
        elif isinstance(value, Mapping):
            for name, item in value.items():
                if name in SAFEGUARDS and item is not SAFEGUARDS[name]:
                    raise LayerEvidenceBundleError(f"restored result safeguard conflicts: {name}")
                visit_restored(item)
        elif isinstance(value, (tuple, list, set, frozenset)):
            for item in value:
                visit_restored(item)
    visit_restored(result)


def _canonical_payload(bundle: Mapping[str, Any]) -> dict[str, Any]:
    payload = deepcopy(dict(bundle))
    payload.pop("created_at_utc", None)
    payload.pop("canonical_payload_sha256", None)
    return payload


def _canonical_payload_sha(bundle: Mapping[str, Any]) -> str:
    return sha256(canonical_json_bytes(_canonical_payload(bundle))).hexdigest()


def _parse_created_at(value: Any) -> None:
    if not isinstance(value, str) or not value.endswith("Z"):
        raise LayerEvidenceBundleError("created_at_utc must be a UTC ISO 8601 string")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise LayerEvidenceBundleError("created_at_utc is malformed") from exc
    if parsed.utcoffset() != timezone.utc.utcoffset(parsed):
        raise LayerEvidenceBundleError("created_at_utc is not UTC")


def validate_layer_evidence_bundle(
    bundle: Mapping[str, Any], *, source_path: str | Path | None = None,
    expected_rna: Mapping[str, Any] | None = None, expected_sample_id: str | None = None,
) -> dict[str, Any]:
    """Validate structure, identity, source bytes, non-empty evidence, and safeguards.

    ``source_path`` defaults to the path stored in the bundle and is always hashed;
    callers may supply a relocated source path with the same basename and bytes.
    """
    if not isinstance(bundle, Mapping):
        raise LayerEvidenceBundleError("bundle must be an object")
    if bundle.get("schema_version") != SCHEMA_VERSION:
        raise LayerEvidenceBundleError(f"unknown schema_version: {bundle.get('schema_version')}")
    top = _require_mapping(bundle, "bundle", _REQUIRED_TOP_LEVEL)
    if top["serializer_format_version"] != SERIALIZER_FORMAT_VERSION:
        raise LayerEvidenceBundleError(
            f"unknown serializer_format_version: {top['serializer_format_version']}"
        )
    if top["bundle_type"] != BUNDLE_TYPE:
        raise LayerEvidenceBundleError(f"unknown bundle_type: {top['bundle_type']}")
    contract = _contract_for_bundle(top)
    if not isinstance(top.get("producer_commit"), str) or not top["producer_commit"].strip():
        raise LayerEvidenceBundleError("producer_commit is required")
    _parse_created_at(top["created_at_utc"])
    source = _require_mapping(top["source"], "source", _REQUIRED_SOURCE)
    rna = _require_mapping(top["rna"], "rna", _REQUIRED_RNA)
    experiment = _require_mapping(top["experiment"], "experiment", _REQUIRED_EXPERIMENT)
    for name in ("path", "basename", "sha256", "run_id", "sample_id", "biological_sample_id"):
        if not isinstance(source[name], str) or not source[name].strip():
            raise LayerEvidenceBundleError(f"source {name} must be non-empty")
    for name in ("name", "sequence", "anticodon", "organism_group", "species"):
        if not isinstance(rna[name], str) or not rna[name].strip():
            raise LayerEvidenceBundleError(f"rna {name} must be non-empty")
    for name in ("condition_name", "digest_type", "independence_group", "shared_source_group"):
        if not isinstance(experiment[name], str) or not experiment[name].strip():
            raise LayerEvidenceBundleError(f"experiment {name} must be non-empty")
    stored_validation = _require_mapping(top["validation"], "validation", _REQUIRED_VALIDATION)
    if experiment["layer"] != contract.layer:
        raise LayerEvidenceBundleError("experiment layer mismatch")
    if Path(str(source["path"])).name != source["basename"]:
        raise LayerEvidenceBundleError("source path basename mismatch")
    actual_path = Path(source_path) if source_path is not None else Path(str(source["path"]))
    if actual_path.name != source["basename"]:
        expected = source["basename"]
        actual = actual_path.name
        if {expected, actual} == {"04 new T1.mzML", "05 old T1.mzML"}:
            raise LayerEvidenceBundleError("04 new T1 cannot substitute for 05 old T1")
        raise LayerEvidenceBundleError("source path basename mismatch")
    try:
        stat = actual_path.stat()
    except OSError as exc:
        raise LayerEvidenceBundleError(f"source unavailable: {actual_path}") from exc
    if stat.st_size != source["size_bytes"]:
        raise LayerEvidenceBundleError("source size mismatch")
    if compute_file_sha256(actual_path) != source["sha256"]:
        raise LayerEvidenceBundleError("source SHA256 mismatch")
    if expected_rna is not None:
        for name in _REQUIRED_RNA:
            if name in expected_rna and rna[name] != expected_rna[name]:
                label = "RNA sequence" if name == "sequence" else "anticodon" if name == "anticodon" else f"RNA {name}"
                raise LayerEvidenceBundleError(f"{label} mismatch")
    if expected_sample_id is not None and source["sample_id"] != expected_sample_id:
        raise LayerEvidenceBundleError(
            f"sample mismatch: expected {expected_sample_id!r}, observed {source['sample_id']!r}"
        )
    if not isinstance(top["result"], Mapping):
        raise LayerEvidenceBundleError("result must be an encoded production object")
    _reject_forbidden_payload(top["result"])
    result = _decode(top["result"])
    if type(result) is not contract.result_class:
        raise LayerEvidenceBundleError(
            f"restored result class mismatch: expected {_type_name(contract.result_class)}, got {_type_name(type(result))}"
        )
    required_present = _required_fields_present(result, contract)
    if not required_present:
        raise LayerEvidenceBundleError("production result missing required fields")
    record_count, non_empty = _record_metrics(result, contract)
    if not non_empty:
        raise LayerEvidenceBundleError(f"empty production result for layer {contract.layer}")
    _validate_safeguards(result, top["result"], _require_mapping(top["safeguards"], "safeguards", set(SAFEGUARDS)))
    if top["canonical_payload_sha256"] != _canonical_payload_sha(top):
        raise LayerEvidenceBundleError("canonical payload SHA256 mismatch")
    expected_validation = {
        "status": "PASSED", "non_empty": True, "record_count": record_count,
        "required_fields_present": True, "provenance_verified": True,
        "safeguards_verified": True,
    }
    for name, expected in expected_validation.items():
        if stored_validation[name] != expected:
            raise LayerEvidenceBundleError(f"validation field mismatch: {name}")
    if not isinstance(stored_validation["validation_messages"], list):
        raise LayerEvidenceBundleError("validation_messages must be a list")
    return dict(stored_validation)


def restore_layer_evidence_result(
    bundle: Mapping[str, Any], *, source_path: str | Path | None = None,
    expected_rna: Mapping[str, Any] | None = None, expected_sample_id: str | None = None,
) -> Any:
    """Validate then restore the exact production root class; no raw parser is invoked."""
    validate_layer_evidence_bundle(
        bundle, source_path=source_path, expected_rna=expected_rna,
        expected_sample_id=expected_sample_id,
    )
    contract = _contract_for_bundle(bundle)
    result = _decode(bundle.get("result"))
    if type(result) is not contract.result_class:
        raise LayerEvidenceBundleError("restored result class does not match layer contract")
    _validate_safeguards(result, bundle.get("result"), bundle.get("safeguards", {}))
    return result


def export_layer_evidence_bundle(
    result: Any, *, layer: str, source_path: str | Path, rna: Mapping[str, Any],
    experiment: Mapping[str, Any], producer_commit: str, output_path: str | Path | None = None,
    created_at_utc: str | None = None, run_id: str | None = None,
    sample_id: str | None = None, biological_sample_id: str | None = None,
) -> dict[str, Any]:
    """Export a production audit result as canonical UTF-8 JSON.

    The write is atomic when ``output_path`` is provided.  The returned mapping is
    already fully validated against the source file and layer non-empty contract.
    """
    contract = LAYER_CONTRACTS.get(layer)
    if contract is None:
        raise LayerEvidenceBundleError(f"unknown layer: {layer}")
    if type(result) is not contract.result_class:
        raise LayerEvidenceBundleError(f"result is not {contract.result_class_name}")
    if not producer_commit or not str(producer_commit).strip():
        raise LayerEvidenceBundleError("producer_commit is required")
    source_file = Path(source_path)
    stat = source_file.stat()
    source = {
        "path": str(source_file), "basename": source_file.name, "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns, "sha256": compute_file_sha256(source_file),
        "run_id": str(run_id or source_file.stem),
        "sample_id": str(sample_id or source_file.stem),
        "biological_sample_id": str(biological_sample_id or sample_id or source_file.stem),
    }
    rna_payload = {name: rna.get(name) for name in sorted(_REQUIRED_RNA)}
    missing_rna = [name for name, value in rna_payload.items() if value is None]
    if missing_rna:
        raise LayerEvidenceBundleError(f"rna missing required fields: {missing_rna}")
    experiment_payload = {name: experiment.get(name) for name in sorted(_REQUIRED_EXPERIMENT)}
    experiment_payload["layer"] = layer
    missing_experiment = [name for name, value in experiment_payload.items() if value is None]
    if missing_experiment:
        raise LayerEvidenceBundleError(f"experiment missing required fields: {missing_experiment}")
    encoded_result = _encode(result)
    record_count, non_empty = _record_metrics(result, contract)
    if not non_empty:
        raise LayerEvidenceBundleError(f"empty production result for layer {layer}")
    _validate_safeguards(result, encoded_result, SAFEGUARDS)
    created = created_at_utc or datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    bundle: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION, "serializer_format_version": SERIALIZER_FORMAT_VERSION,
        "bundle_type": BUNDLE_TYPE, "layer": layer,
        "optional_result_key": contract.optional_result_key,
        "producer_name": contract.producer_name, "producer_module": contract.producer_module,
        "producer_commit": str(producer_commit), "created_at_utc": created,
        "source": source, "rna": rna_payload, "experiment": experiment_payload,
        "result": encoded_result, "safeguards": dict(SAFEGUARDS),
        "validation": {
            "status": "PASSED", "non_empty": True, "record_count": record_count,
            "required_fields_present": _required_fields_present(result, contract),
            "provenance_verified": True, "safeguards_verified": True,
            "validation_messages": [f"{layer} production result validated"],
        },
    }
    bundle["canonical_payload_sha256"] = _canonical_payload_sha(bundle)
    validate_layer_evidence_bundle(bundle, source_path=source_file, expected_rna=rna_payload,
                                   expected_sample_id=source["sample_id"])
    if output_path is not None:
        destination = Path(output_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(destination.name + f".tmp-{os.getpid()}")
        try:
            temporary.write_bytes(canonical_json_bytes(bundle))
            os.replace(temporary, destination)
        finally:
            if temporary.exists():
                temporary.unlink()
    return bundle


def _reject_constant(value: str) -> None:
    raise LayerEvidenceBundleError(f"non-finite JSON constant is forbidden: {value}")


def load_layer_evidence_bundle(
    path: str | Path, *, source_path: str | Path | None = None,
    expected_rna: Mapping[str, Any] | None = None, expected_sample_id: str | None = None,
    restore: bool = False,
) -> dict[str, Any] | Any:
    """Load and validate a schema-1.0 bundle, optionally restoring its result."""
    try:
        raw = Path(path).read_text(encoding="utf-8")
        bundle = json.loads(raw, parse_constant=_reject_constant)
    except LayerEvidenceBundleError:
        raise
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise LayerEvidenceBundleError(f"malformed or unreadable JSON bundle: {exc}") from exc
    if not isinstance(bundle, dict):
        raise LayerEvidenceBundleError("bundle JSON root must be an object")
    validate_layer_evidence_bundle(
        bundle, source_path=source_path, expected_rna=expected_rna,
        expected_sample_id=expected_sample_id,
    )
    return restore_layer_evidence_result(
        bundle, source_path=source_path, expected_rna=expected_rna,
        expected_sample_id=expected_sample_id,
    ) if restore else bundle


def compare_layer_evidence_provenance(
    left: Mapping[str, Any], right: Mapping[str, Any],
) -> dict[str, Any]:
    """Compare RNA, sample, source, condition, sharing, and independence metadata."""
    left_source = _require_mapping(left.get("source"), "left source", _REQUIRED_SOURCE)
    right_source = _require_mapping(right.get("source"), "right source", _REQUIRED_SOURCE)
    left_rna = _require_mapping(left.get("rna"), "left rna", _REQUIRED_RNA)
    right_rna = _require_mapping(right.get("rna"), "right rna", _REQUIRED_RNA)
    left_exp = _require_mapping(left.get("experiment"), "left experiment", _REQUIRED_EXPERIMENT)
    right_exp = _require_mapping(right.get("experiment"), "right experiment", _REQUIRED_EXPERIMENT)
    same_rna_identity = left_rna["name"] == right_rna["name"]
    same_sequence = left_rna["sequence"] == right_rna["sequence"]
    same_anticodon = left_rna["anticodon"] == right_rna["anticodon"]
    same_taxonomy = (
        left_rna["organism_group"] == right_rna["organism_group"]
        and left_rna["species"] == right_rna["species"]
    )
    same_sample = left_source["sample_id"] == right_source["sample_id"]
    same_raw_source = (
        left_source["sha256"] == right_source["sha256"]
        and left_source["run_id"] == right_source["run_id"]
    )
    shared_source = (
        left_exp["shared_source_group"] == right_exp["shared_source_group"]
        and same_raw_source
    )
    same_independence_group = left_exp["independence_group"] == right_exp["independence_group"]
    p1_pair = {left.get("layer"), right.get("layer")} == {"P1AP_MS1", "P1AP_MS2"}
    p1_pair_compatible = not p1_pair or (same_raw_source and shared_source and same_independence_group)
    return {
        "same_rna_identity": same_rna_identity,
        "same_sequence": same_sequence,
        "same_anticodon": same_anticodon,
        "same_organism_species": same_taxonomy,
        "compatible_condition": left_exp["condition_name"] == right_exp["condition_name"],
        "same_sample": same_sample,
        "different_sample": not same_sample,
        "same_raw_source": same_raw_source,
        "shared_source_relationship": shared_source,
        "same_independence_group": same_independence_group,
        "independent_support": not same_independence_group and same_sample,
        "p1ap_ms1_ms2_compatible": p1_pair_compatible,
        "rna_compatible": same_rna_identity and same_sequence and same_anticodon and same_taxonomy,
    }
