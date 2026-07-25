"""Command-line entry point for the production SCIEX bundle-only runner.

Example::

    python -m rna_masshunter.sciex_cross_layer_bundle_cli \
      --full-bundle bundle_FULL.json --t1-bundle bundle_T1.json \
      --p1ap-ms1-bundle bundle_P1AP_MS1.json \
      --p1ap-ms2-bundle bundle_P1AP_MS2.json \
      --aggregate-json cross_layer_aggregate.json --excel-output cross_layer.xlsx

All four production bundles are explicit and required for an active run.  The
CLI never discovers bundles or reparses raw mzML, refuses output overwrites, and
supports ``--dry-run`` validation without writes.  FULL node-level provenance
remains limited and is reported as a warning while bundle-level provenance stays
verified.  Results are shadow-only and never propagate into formal analysis.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import os
from pathlib import Path
import sys
import traceback
from typing import Any, Mapping, Sequence

import yaml

from rna_masshunter.sciex_cross_layer_bundle_runner import (
    CrossLayerBundleRunResult,
    CrossLayerBundleRunnerError,
    REQUIRED_LAYERS,
    run_cross_layer_from_bundles,
)
from rna_masshunter.sciex_layer_evidence_bundle import (
    LayerEvidenceBundleError,
    canonical_json_bytes,
)

EXIT_SUCCESS = 0
EXIT_ARGUMENT_CONFIG = 2
EXIT_BUNDLE_VALIDATION = 3
EXIT_COMPATIBILITY = 4
EXIT_OUTPUT = 5
EXIT_AGGREGATION = 6

_LAYER_ARGUMENTS = {
    "FULL": "full_bundle",
    "T1": "t1_bundle",
    "P1AP_MS1": "p1ap_ms1_bundle",
    "P1AP_MS2": "p1ap_ms2_bundle",
}
_CONFIG_KEYS = {"enabled", "bundles", "output", "overwrite"}
_CONFIG_BUNDLE_KEYS = {"full", "t1", "p1ap_ms1", "p1ap_ms2"}
_CONFIG_OUTPUT_KEYS = {"aggregate_json", "excel", "summary_json"}


class CrossLayerBundleCLIError(ValueError):
    """Categorized CLI error with a stable process exit code."""

    def __init__(self, category: str, message: str, exit_code: int):
        super().__init__(message)
        self.category = category
        self.exit_code = exit_code


@dataclass(frozen=True)
class CrossLayerBundleCLIConfig:
    enabled: bool
    bundles: dict[str, Path | None]
    aggregate_json: Path | None
    excel_output: Path | None
    summary_json: Path | None
    overwrite: bool = False
    source_path: Path | None = None


def _mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise CrossLayerBundleCLIError("CONFIG_ERROR", f"{name} must be a mapping", EXIT_ARGUMENT_CONFIG)
    return value


def _unknown_keys(value: Mapping[str, Any], allowed: set[str], name: str) -> None:
    unknown = sorted(set(map(str, value)) - allowed)
    if unknown:
        raise CrossLayerBundleCLIError(
            "CONFIG_ERROR", f"unknown {name} keys: {unknown}", EXIT_ARGUMENT_CONFIG,
        )


def _config_path(value: Any, base: Path, label: str) -> Path | None:
    if value is None or str(value).strip() == "":
        return None
    if not isinstance(value, (str, Path)):
        raise CrossLayerBundleCLIError("CONFIG_ERROR", f"{label} must be a path string or null", EXIT_ARGUMENT_CONFIG)
    path = Path(value).expanduser()
    return path if path.is_absolute() else base / path


def load_cross_layer_bundle_cli_config(path: str | Path) -> CrossLayerBundleCLIConfig:
    """Load only the dedicated runner section; relative paths use the config directory."""
    source = Path(path).expanduser().resolve()
    try:
        payload = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise CrossLayerBundleCLIError(
            "CONFIG_ERROR", f"cannot read config {source}: {exc}", EXIT_ARGUMENT_CONFIG,
        ) from exc
    root = _mapping(payload, "config")
    section = root.get("cross_layer_bundle_runner")
    if section is None:
        return CrossLayerBundleCLIConfig(False, {layer: None for layer in REQUIRED_LAYERS}, None, None, None, source_path=source)
    section = _mapping(section, "cross_layer_bundle_runner")
    _unknown_keys(section, _CONFIG_KEYS, "cross_layer_bundle_runner")
    bundles = _mapping(section.get("bundles") or {}, "cross_layer_bundle_runner.bundles")
    output = _mapping(section.get("output") or {}, "cross_layer_bundle_runner.output")
    _unknown_keys(bundles, _CONFIG_BUNDLE_KEYS, "bundle")
    _unknown_keys(output, _CONFIG_OUTPUT_KEYS, "output")
    enabled = section.get("enabled", False)
    overwrite = section.get("overwrite", False)
    if type(enabled) is not bool or type(overwrite) is not bool:
        raise CrossLayerBundleCLIError("CONFIG_ERROR", "enabled and overwrite must be boolean", EXIT_ARGUMENT_CONFIG)
    if overwrite:
        raise CrossLayerBundleCLIError("CONFIG_ERROR", "overwrite=true is unsupported; outputs are immutable", EXIT_ARGUMENT_CONFIG)
    base = source.parent
    return CrossLayerBundleCLIConfig(
        enabled=enabled,
        bundles={
            "FULL": _config_path(bundles.get("full"), base, "bundles.full"),
            "T1": _config_path(bundles.get("t1"), base, "bundles.t1"),
            "P1AP_MS1": _config_path(bundles.get("p1ap_ms1"), base, "bundles.p1ap_ms1"),
            "P1AP_MS2": _config_path(bundles.get("p1ap_ms2"), base, "bundles.p1ap_ms2"),
        },
        aggregate_json=_config_path(output.get("aggregate_json"), base, "output.aggregate_json"),
        excel_output=_config_path(output.get("excel"), base, "output.excel"),
        summary_json=_config_path(output.get("summary_json"), base, "output.summary_json"),
        overwrite=False,
        source_path=source,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m rna_masshunter.sciex_cross_layer_bundle_cli",
        description=(
            "Validate four explicit production SCIEX bundles and run bundle-only "
            "cross-layer reconciliation without reading raw mzML."
        ),
        epilog=(
            "Outputs never overwrite existing files. --dry-run performs validation, "
            "restore, compatibility, and aggregation without writing files."
        ),
    )
    parser.add_argument("--config", help="YAML file containing only an optional cross_layer_bundle_runner section.")
    parser.add_argument("--full-bundle", help="Explicit FULL production bundle JSON.")
    parser.add_argument("--t1-bundle", help="Explicit T1 production bundle JSON.")
    parser.add_argument("--p1ap-ms1-bundle", help="Explicit P1/AP MS1 production bundle JSON.")
    parser.add_argument("--p1ap-ms2-bundle", help="Explicit P1/AP MS2 production bundle JSON.")
    parser.add_argument("--aggregate-json", help="Aggregate JSON output path (must not exist).")
    parser.add_argument("--excel-output", help="XL-only Excel output path (must not exist).")
    parser.add_argument("--summary-json", help="Optional CLI summary JSON path (must not exist).")
    parser.add_argument("--dry-run", action="store_true", help="Validate and aggregate without writing any output.")
    parser.add_argument("--debug", action="store_true", help="Show a traceback for unexpected failures.")
    return parser


def _cli_path(value: str | None) -> Path | None:
    return Path(value).expanduser() if value else None


def _explicit_mode(args: argparse.Namespace) -> bool:
    return any(
        getattr(args, name) is not None
        for name in (*_LAYER_ARGUMENTS.values(), "aggregate_json", "excel_output", "summary_json")
    ) or bool(args.dry_run)


def _effective_settings(args: argparse.Namespace) -> tuple[bool, dict[str, Path | None], Path | None, Path | None, Path | None]:
    config = (
        load_cross_layer_bundle_cli_config(args.config)
        if args.config else
        CrossLayerBundleCLIConfig(False, {layer: None for layer in REQUIRED_LAYERS}, None, None, None)
    )
    explicit = _explicit_mode(args)
    active = explicit or config.enabled or not bool(args.config)
    bundles = dict(config.bundles)
    for layer, argument in _LAYER_ARGUMENTS.items():
        override = _cli_path(getattr(args, argument))
        if override is not None:
            bundles[layer] = override
    aggregate = _cli_path(args.aggregate_json) or config.aggregate_json
    excel = _cli_path(args.excel_output) or config.excel_output
    summary = _cli_path(args.summary_json) or config.summary_json
    return active, bundles, aggregate, excel, summary


def _absolute_identity(path: Path) -> Path:
    try:
        return path.resolve(strict=False)
    except OSError:
        return Path(os.path.abspath(path))


def _preflight(
    bundles: Mapping[str, Path | None], outputs: Sequence[Path | None], *, dry_run: bool,
) -> dict[str, Path]:
    missing = [layer for layer in REQUIRED_LAYERS if bundles.get(layer) is None]
    if missing:
        raise CrossLayerBundleCLIError(
            "ARGUMENT_ERROR", f"missing required bundle arguments: {missing}", EXIT_ARGUMENT_CONFIG,
        )
    inputs = {layer: Path(bundles[layer]) for layer in REQUIRED_LAYERS}  # type: ignore[arg-type]
    for layer, path in inputs.items():
        if path.suffix.lower() == ".mzml":
            raise CrossLayerBundleCLIError(
                "ARGUMENT_ERROR", f"{layer} input is raw mzML, not a production bundle: {path}", EXIT_ARGUMENT_CONFIG,
            )
        if not path.is_file():
            raise CrossLayerBundleCLIError(
                "BUNDLE_VALIDATION_ERROR", f"{layer} bundle file is unavailable: {path}", EXIT_BUNDLE_VALIDATION,
            )
    if not dry_run and not any(path is not None for path in outputs):
        raise CrossLayerBundleCLIError(
            "ARGUMENT_ERROR", "at least one output is required unless --dry-run is used", EXIT_ARGUMENT_CONFIG,
        )
    actual_outputs = [path for path in outputs if path is not None]
    identities: dict[Path, str] = {}
    for layer, path in inputs.items():
        identities[_absolute_identity(path)] = f"input:{layer}"
    for index, path in enumerate(actual_outputs):
        identity = _absolute_identity(path)
        if identity in identities:
            raise CrossLayerBundleCLIError(
                "OUTPUT_ERROR", f"input/output path collision: {path} conflicts with {identities[identity]}", EXIT_OUTPUT,
            )
        if identity in {_absolute_identity(other) for other in actual_outputs[:index]}:
            raise CrossLayerBundleCLIError("OUTPUT_ERROR", f"output path collision: {path}", EXIT_OUTPUT)
        if path.exists():
            raise CrossLayerBundleCLIError("OUTPUT_ERROR", f"output file already exists: {path}", EXIT_OUTPUT)
        if path.parent.exists() and not path.parent.is_dir():
            raise CrossLayerBundleCLIError("OUTPUT_ERROR", f"malformed output directory: {path.parent}", EXIT_OUTPUT)
    return inputs


def _summary(run: CrossLayerBundleRunResult, *, outputs: Mapping[str, Path | None], dry_run: bool) -> dict[str, Any]:
    report = run.compatibility_report
    return {
        "status": "DRY_RUN_PASSED" if dry_run else "PASSED",
        "dry_run": dry_run,
        "bundle_paths": run.input_bundle_paths,
        "bundle_sha256": run.input_bundle_sha256,
        "compatibility": {
            "compatible": report.compatible,
            "blocking_errors": list(report.blocking_errors),
            "bundle_level_provenance_verified": report.bundle_level_provenance_verified,
            "node_level_provenance": report.node_level_provenance,
            "shared_source_relationships": list(report.shared_source_relationships),
            "independence_groups": report.independence_groups,
        },
        "warnings": list(run.warnings),
        "aggregate_counts": run.aggregate_counts,
        "consensus": run.consensus_status,
        "confidence": run.confidence,
        "safeguards": run.safeguard_summary,
        "outputs": {
            "aggregate_json": str(outputs.get("aggregate_json")) if outputs.get("aggregate_json") else None,
            "excel": str(outputs.get("excel")) if outputs.get("excel") else None,
            "summary_json": str(outputs.get("summary_json")) if outputs.get("summary_json") else None,
        },
    }


def _atomic_summary(path: Path, payload: Mapping[str, Any]) -> Path:
    if path.exists():
        raise CrossLayerBundleCLIError("OUTPUT_ERROR", f"output file already exists: {path}", EXIT_OUTPUT)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    if temporary.exists():
        raise CrossLayerBundleCLIError("OUTPUT_ERROR", f"summary staging file already exists: {temporary}", EXIT_OUTPUT)
    try:
        with temporary.open("xb") as handle:
            handle.write(canonical_json_bytes(payload) + b"\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
    return path


def _staged_path(path: Path) -> Path:
    suffix = ".xlsx" if path.suffix.lower() == ".xlsx" else ""
    return path.with_name(f".{path.name}.cli-stage-{os.getpid()}{suffix}")


def _classify_exception(exc: Exception) -> CrossLayerBundleCLIError:
    if isinstance(exc, CrossLayerBundleCLIError):
        return exc
    if isinstance(exc, LayerEvidenceBundleError):
        return CrossLayerBundleCLIError("BUNDLE_VALIDATION_ERROR", str(exc), EXIT_BUNDLE_VALIDATION)
    text = str(exc)
    if isinstance(exc, CrossLayerBundleRunnerError):
        if "incompatible bundle set" in text:
            return CrossLayerBundleCLIError("COMPATIBILITY_ERROR", text, EXIT_COMPATIBILITY)
        if any(token in text for token in ("output", "directory", "workbook", "Excel", "staging")):
            return CrossLayerBundleCLIError("OUTPUT_ERROR", text, EXIT_OUTPUT)
        if any(token in text for token in ("bundle file", "missing required layer", "duplicate layer", "unsupported layer", "declared layer")):
            return CrossLayerBundleCLIError("BUNDLE_VALIDATION_ERROR", text, EXIT_BUNDLE_VALIDATION)
        return CrossLayerBundleCLIError("AGGREGATION_ERROR", text, EXIT_AGGREGATION)
    return CrossLayerBundleCLIError("AGGREGATION_ERROR", f"{type(exc).__name__}: {exc}", EXIT_AGGREGATION)


def execute(args: argparse.Namespace) -> dict[str, Any]:
    active, bundles, aggregate, excel, summary_path = _effective_settings(args)
    if not active:
        return {
            "status": "DISABLED", "dry_run": False, "bundle_paths": {},
            "bundle_sha256": {}, "compatibility": None, "warnings": [],
            "aggregate_counts": {}, "consensus": None, "confidence": None,
            "safeguards": None,
            "outputs": {"aggregate_json": None, "excel": None, "summary_json": None},
        }
    if not args.dry_run and aggregate is None and excel is None:
        raise CrossLayerBundleCLIError(
            "ARGUMENT_ERROR",
            "at least --aggregate-json or --excel-output is required unless --dry-run is used",
            EXIT_ARGUMENT_CONFIG,
        )
    if args.dry_run and any(path is not None for path in (aggregate, excel, summary_path)):
        # Paths may be supplied for parity with a real run, but dry-run never writes them.
        pass
    inputs = _preflight(bundles, (aggregate, excel, summary_path), dry_run=args.dry_run)
    if args.dry_run:
        run = run_cross_layer_from_bundles(inputs)
        return _summary(
            run,
            outputs={"aggregate_json": None, "excel": None, "summary_json": None},
            dry_run=True,
        )

    desired = {"aggregate_json": aggregate, "excel": excel, "summary_json": summary_path}
    stages = {name: _staged_path(path) for name, path in desired.items() if path is not None}
    for stage in stages.values():
        if stage.exists():
            raise CrossLayerBundleCLIError("OUTPUT_ERROR", f"CLI staging file already exists: {stage}", EXIT_OUTPUT)
    try:
        run = run_cross_layer_from_bundles(
            inputs,
            output_json_path=stages.get("aggregate_json"),
            output_excel_path=stages.get("excel"),
        )
        summary_payload = _summary(run, outputs=desired, dry_run=False)
        if "summary_json" in stages:
            _atomic_summary(stages["summary_json"], summary_payload)
        committed: list[Path] = []
        try:
            for name, stage in stages.items():
                destination = desired[name]
                assert destination is not None
                destination.parent.mkdir(parents=True, exist_ok=True)
                os.replace(stage, destination)
                committed.append(destination)
        except Exception:
            for destination in committed:
                if destination.exists():
                    destination.unlink()
            raise
        return summary_payload
    except Exception:
        for stage in stages.values():
            if stage.exists():
                stage.unlink()
        raise


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        summary = execute(args)
    except Exception as exc:
        error = _classify_exception(exc)
        print(json.dumps({
            "status": "FAILED", "error_category": error.category,
            "exit_code": error.exit_code, "blocking_reason": str(error),
        }, ensure_ascii=False, sort_keys=True), file=sys.stderr)
        if args.debug:
            traceback.print_exc()
        return error.exit_code
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return EXIT_SUCCESS


if __name__ == "__main__":
    raise SystemExit(main())
