from pathlib import Path
from typing import Any

import yaml

from rna_masshunter.models import RunConfig
from rna_masshunter.warnings_manager import add_warning


DEFAULT_CONFIG: dict[str, dict[str, Any]] = {
    "project": {"name": "RNA_MassHunter_v2", "output_dir": "output", "log_dir": "logs", "cache_dir": ".cache"},
    "input": {"raw_path": "", "mzml_path": "", "msconvert_path": ""},
    "organism": {"group": "archaea", "species": "", "rule_set": "archaea_general"},
    "sequence": {"name": "target_tRNA", "type": "RNA", "sequence": "", "anticodon": "", "wobble_position": 34},
    "experiment": {"condition_name": "wild_type"},
    "instrument": {"polarity": "negative", "ms1_tolerance_ppm": 10, "ms2_tolerance_da": 0.02},
    "reconstruction": {
        "enabled": True,
        "rt_min": None,
        "rt_max": None,
        "mz_min": 500,
        "mz_max": 3000,
        "intensity_threshold": 1000,
        "min_charge": 5,
        "max_charge": 40,
        "min_charge_states": 3,
        "mass_cluster_tolerance_da": 1.0,
    },
    "digestion": {
        "enabled": True,
        "enzyme": "RNase_T1",
        "missed_cleavages": 1,
        "min_length": 2,
        "include_terminal_forms": True,
        "allow_partial_digestion": True,
        "allow_nonspecific_cleavage": False,
    },
    "alkaline_phosphatase": {
        "enabled": True,
        "assume_complete": False,
        "allow_residual_phosphate": True,
        "allow_cyclic_phosphate": True,
    },
    "peak_filtering": {
        "major_intensity_threshold": 25000,
        "minor_intensity_threshold": 5000,
        "trace_intensity_threshold": 1000,
        "report_trace_peaks": True,
        "use_trace_peaks_for_final_call": True,
    },
    "performance": {"cache_enabled": True, "checkpoint_enabled": True},
    "reporting": {"excel_output": True},
}


def _merge_defaults(data: dict[str, Any], warnings: list[dict[str, Any]] | None) -> dict[str, Any]:
    merged = {}
    for section, defaults in DEFAULT_CONFIG.items():
        current = data.get(section, {})
        if not isinstance(current, dict):
            current = {}
            if warnings is not None:
                add_warning(warnings, "WARNING", "config", f"Config section '{section}' was not a mapping; defaults were used.")
        missing = sorted(set(defaults) - set(current))
        if missing and warnings is not None:
            add_warning(warnings, "WARNING", "config", f"Config section '{section}' missing keys filled by defaults.", missing)
        merged[section] = {**defaults, **current}
    merged["raw"] = data
    return merged


def load_config(config_path: str | Path, warnings: list[dict[str, Any]] | None = None) -> RunConfig:
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"config.yaml not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    merged = _merge_defaults(data, warnings)
    return RunConfig(**merged)


def validate_config(config: RunConfig, warnings: list[dict[str, Any]] | None = None) -> None:
    polarity = str(config.instrument.get("polarity", "")).lower()
    if polarity not in {"negative", "positive"} and warnings is not None:
        add_warning(warnings, "ERROR", "config", "instrument.polarity must be 'negative' or 'positive'.", polarity)
    if not config.input.get("mzml_path") and not config.input.get("raw_path") and warnings is not None:
        add_warning(warnings, "WARNING", "config", "input.mzml_path and input.raw_path are empty; sample config mode.")
    if not config.sequence.get("sequence") and warnings is not None:
        add_warning(warnings, "WARNING", "config", "sequence.sequence is empty; theoretical mass will be skipped.")


def resolve_paths(config: RunConfig, project_root: str | Path) -> RunConfig:
    root = Path(project_root)
    for section_name in ("project", "input"):
        section = getattr(config, section_name)
        for key, value in list(section.items()):
            if not value or not isinstance(value, str):
                continue
            if key.endswith("_dir") or key.endswith("_path"):
                path = Path(value)
                if not path.is_absolute():
                    section[key] = str(root / path)
    return config
