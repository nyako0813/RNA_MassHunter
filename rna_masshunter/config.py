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
        "digestion_mode": None,
        "missed_cleavages": 1,
        "min_length": 2,
        "max_length": None,
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
    "modification_search": {
        "enabled": True,
        "source": {"use_fragments": True, "use_intact": False},
        "mz_tolerance_ppm": 10,
        "max_candidates_per_match": 10,
        "include_isobaric_modifications": False,
        "isobaric_mass_shift_tolerance_da": 1e-6,
        "require_base_compatibility": True,
        "allow_unknown_position_within_fragment": True,
        "report_unmodified_explained": False,
        "min_peak_tier": ["Major", "Minor"],
        "min_confidence": ["High", "Medium"],
    },
    "p1_annotation": {
        "enabled": True,
        "include_unmatched_peaks": True,
        "include_modified_monomers": True,
        "include_phosphate_forms": True,
        "mz_tolerance_ppm": 10,
        "min_intensity": 0,
        "charge_states": [1],
    },
    "ms2_annotation": {
        "enabled": True,
        "mz_tolerance_ppm": 20,
        "min_peak_intensity": 10,
        "min_relative_intensity_percent": 1.0,
        "max_peaks_per_spectrum": 500,
        "precursor_match_tolerance_ppm": 20,
        "constrain_by_precursor": True,
        "fallback_to_all_ions_if_no_precursor_match": False,
        "use_theoretical_fragments": True,
        "include_neutral_loss": False,
        "include_base_loss": False,
        "min_ion_length": 1,
        "max_ion_length": None,
        "min_ion_length_for_evidence": 2,
        "max_ms2_match_rows": 100000,
        "max_unmatched_peaks": 50000,
        "output_unmatched_peaks": True,
        "output_low_intensity_peaks": False,
        "output_all_peak_annotations": False,
        "include_modified_precursor_candidates": True,
        "modified_precursor_source": "known_modifications",
        "modified_precursor_max_mods_per_fragment": 1,
        "modified_precursor_require_base_compatibility": True,
        "modified_precursor_include_isobaric": False,
        "modified_precursor_mass_shift_tolerance_da": 1e-6,
        "modified_precursor_max_candidates_per_spectrum": 20,
        "include_modified_fragment_ions": True,
        "modified_fragment_ion_source": "modified_precursor_candidates",
        "modified_fragment_max_positions_per_candidate": 20,
        "modified_fragment_require_target_base": True,
        "modified_fragment_include_unmodified_counterparts": True,
        "modified_fragment_min_ion_length": 1,
        "modified_fragment_min_ion_length_for_localization": 2,
        "modified_fragment_output_all_position_candidates": True,
        "modified_fragment_max_rows": 100000,
        "annotate_position_discriminating_ions": True,
    },
    "modification_evidence_ranking": {
        "enabled": True,
        "use_ms1_fragment_evidence": True,
        "use_known_modification_candidates": True,
        "use_ms2_precursor_evidence": True,
        "use_ms2_modified_ion_evidence": True,
        "use_localization_evidence": True,
        "use_organism_rules": True,
        "use_trna_context": True,
        "require_ms2_evidence_for_high_confidence": True,
        "use_biological_context": True,
        "cap_context_only_confidence": "Medium",
        "require_ms_evidence_for_context_boosted_high": True,
        "enable_ambiguity_grouping": True,
        "ambiguity_group_by": ["Spectrum_ID", "Parent_Fragment_ID", "Modification_ID"],
        "require_position_discriminating_ions_for_localization_confidence": True,
        "min_discriminating_ions_for_position_support": 1,
        "min_informative_discriminating_ions_for_high": 2,
        "collapse_ambiguous_positions_in_summary": True,
        "min_final_score_to_report": 0,
        "max_ranked_candidates": 10000,
        "weights": {
            "ms1_fragment_match": 1.0, "known_modification_candidate": 1.0,
            "ms2_precursor_rescue": 2.0, "ms2_modified_ion_match": 2.0,
            "localization_weak": 1.0, "localization_moderate": 3.0, "localization_strong": 5.0,
            "organism_rule_supported": 1.5, "trna_context_supported": 1.0,
            "low_information_penalty": -1.0, "ambiguous_position_penalty": -1.0,
            "isobaric_precursor_penalty": -2.0,
            "ambiguity_penalty": -1.5, "non_discriminating_ion_penalty": -0.5,
            "curation_manually_checked": 0.5, "source_user_pdf": 0.5,
            "detectability_ms2_supported": 0.5,
        },
    },
    "biological_context": {
        "enabled": True,
        "organism_group": "", "organism_species": "",
        "trna_name": "", "trna_type": "", "anticodon": "",
        "focus_positions": [], "focus_position_window": 2,
        "priority_modifications": [], "priority_keywords": [],
        "boost": {
            "organism_rule_supported": 1.0, "pathway_supported": 1.0,
            "priority_modification": 1.5, "priority_keyword_match": 0.75,
            "focus_position_match": 1.0, "focus_position_nearby": 0.5,
            "trna_context_supported": 1.0,
        },
        "penalties": {"organism_context_conflict": -2.0, "unrelated_modification_family": -0.5},
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
