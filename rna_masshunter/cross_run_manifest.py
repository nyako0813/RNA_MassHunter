"""Explicit cross-run manifest validation and replicate-independence classification."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import yaml

REQUIRED_FIELDS = (
    "run_id", "mzml_path", "sample_id", "biological_replicate_id",
    "sample_preparation_id", "digestion_id", "technical_replicate_id",
    "acquisition_batch_id", "instrument_method_id", "condition", "enzyme",
    "sequence_id", "organism",
)
OPTIONAL_FIELDS = ("notes",)
INDEPENDENCE_ORDER = {
    "UNKNOWN_INDEPENDENCE": 0, "SAME_INJECTION": 1, "TECHNICAL_REPLICATE": 2,
    "INDEPENDENT_INJECTION": 3, "INDEPENDENT_DIGESTION": 4,
    "INDEPENDENT_SAMPLE_PREPARATION": 5, "BIOLOGICAL_REPLICATE": 6,
}

@dataclass(frozen=True)
class CrossRunManifest:
    schema_version: int
    runs: tuple[dict[str, Any], ...]
    source_path: Path

class ManifestValidationError(ValueError):
    pass

def _text(value: Any) -> str:
    return str(value).strip() if value is not None else ""

def load_cross_run_manifest(path: str | Path, *, require_files: bool = True) -> CrossRunManifest:
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise ManifestValidationError(f"manifest_not_found:{source}")
    payload = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    if payload.get("schema_version") != 1:
        raise ManifestValidationError("unsupported_schema_version: expected 1")
    raw_runs = payload.get("runs")
    if not isinstance(raw_runs, list) or not raw_runs:
        raise ManifestValidationError("runs must be a non-empty list")
    ids: set[str] = set(); paths: set[str] = set(); runs = []
    for index, raw in enumerate(raw_runs, 1):
        if not isinstance(raw, dict):
            raise ManifestValidationError(f"run[{index}]:not_a_mapping")
        missing = [field for field in REQUIRED_FIELDS if not _text(raw.get(field))]
        if missing:
            raise ManifestValidationError(f"run[{index}]:missing_metadata:{','.join(missing)}")
        run_id = _text(raw["run_id"])
        if run_id in ids:
            raise ManifestValidationError(f"duplicate_run_id:{run_id}")
        mzml = Path(_text(raw["mzml_path"])).expanduser()
        if not mzml.is_absolute():
            mzml = (source.parent / mzml).resolve()
        else:
            mzml = mzml.resolve()
        normalized = str(mzml)
        if normalized in paths:
            raise ManifestValidationError(f"duplicate_mzml_path:{normalized}")
        if require_files and not mzml.is_file():
            raise ManifestValidationError(f"missing_mzml_file:{run_id}:{normalized}")
        item = {field: _text(raw.get(field)) for field in REQUIRED_FIELDS + OPTIONAL_FIELDS}
        item["mzml_path"] = normalized
        ids.add(run_id); paths.add(normalized); runs.append(item)
    return CrossRunManifest(1, tuple(runs), source)

def classify_run_independence(left: dict[str, Any], right: dict[str, Any]) -> str:
    """Classify a run pair from strongest biological distinction downwards."""
    if _text(left.get("run_id")) == _text(right.get("run_id")) or (
        _text(left.get("mzml_path")) and _text(left.get("mzml_path")) == _text(right.get("mzml_path"))
    ):
        return "SAME_INJECTION"
    required = ("sample_id", "biological_replicate_id", "sample_preparation_id", "digestion_id", "technical_replicate_id")
    if any(not _text(left.get(x)) or not _text(right.get(x)) for x in required):
        return "UNKNOWN_INDEPENDENCE"
    if left["biological_replicate_id"] != right["biological_replicate_id"] or left["sample_id"] != right["sample_id"]:
        return "BIOLOGICAL_REPLICATE"
    if left["sample_preparation_id"] != right["sample_preparation_id"]:
        return "INDEPENDENT_SAMPLE_PREPARATION"
    if left["digestion_id"] != right["digestion_id"]:
        return "INDEPENDENT_DIGESTION"
    if left["technical_replicate_id"] != right["technical_replicate_id"]:
        return "TECHNICAL_REPLICATE"
    return "INDEPENDENT_INJECTION"

def strongest_independence(runs: list[dict[str, Any]]) -> str:
    if len(runs) < 2:
        return "SAME_INJECTION" if runs else "UNKNOWN_INDEPENDENCE"
    levels = [classify_run_independence(runs[i], runs[j]) for i in range(len(runs)) for j in range(i + 1, len(runs))]
    return max(levels, key=lambda x: INDEPENDENCE_ORDER[x])
