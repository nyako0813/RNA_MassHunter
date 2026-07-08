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
        "max_charge_state_peak_rows": 100000,
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
        "enabled": False,
        "assume_complete": False,
        "allow_residual_phosphate": True,
        "allow_cyclic_phosphate": True,
    },
    "fragment_mapping": {
        "enabled": True,
        "mz_tolerance_ppm": 10,
        "min_charge": 1,
        "max_charge": 8,
        "polarity": "auto",
        "use_peak_tiers": True,
        "include_trace_peaks": True,
        "max_matches_per_fragment": 20,
        "min_fragment_length_for_filtered": 3,
        "filtered_peak_tiers": ["Major", "Minor"],
        "filtered_confidence": ["High", "Medium"],
        "summary_best_match_by": "fragment_id",
    },
    "peak_filtering": {
        "major_intensity_threshold": 25000,
        "minor_intensity_threshold": 5000,
        "trace_intensity_threshold": 1000,
        "report_trace_peaks": True,
        "use_trace_peaks_for_final_call": True,
    },
    "performance": {"cache_enabled": True, "checkpoint_enabled": True},
    "reporting": {
        "excel_output": True,
        "max_excel_rows_per_sheet": 100000,
        "truncate_large_sheets": True,
        "max_charge_state_peak_rows": 100000,
    },
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


def _as_bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


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

    digestion_enabled = _as_bool(config.digestion.get("enabled"), True)
    reconstruction_enabled = _as_bool(config.reconstruction.get("enabled"), True)
    fragment_mapping_enabled = _as_bool(config.fragment_mapping.get("enabled"), True)
    if not digestion_enabled and not reconstruction_enabled and warnings is not None:
        add_warning(
            warnings,
            "WARNING",
            "config",
            "Both digestion.enabled and reconstruction.enabled are false; no fragment or intact-mass analysis will be performed.",
        )
    if not digestion_enabled and fragment_mapping_enabled and warnings is not None:
        add_warning(
            warnings,
            "INFO",
            "config",
            "digestion.enabled is false, so fragment_mapping.enabled will be ignored for this run.",
        )

    ap_enabled = _as_bool(config.alkaline_phosphatase.get("enabled"), False)
    ap_complete = _as_bool(config.alkaline_phosphatase.get("assume_complete"), False)
    if not ap_enabled and ap_complete and warnings is not None:
        add_warning(
            warnings,
            "WARNING",
            "config",
            "alkaline_phosphatase.assume_complete is true but alkaline_phosphatase.enabled is false; assume_complete will be ignored.",
        )


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