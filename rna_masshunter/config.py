from pathlib import Path
from typing import Any

import yaml

from rna_masshunter.models import RunConfig
from rna_masshunter.sciex_delta_mass_cluster_audit import DeltaMassClusterParameters
from rna_masshunter.sciex_spacing_resolution_audit import SpacingResolutionParameters
from rna_masshunter.sciex_relation_evidence_quality_audit import RelationEvidenceQualityParameters
from rna_masshunter.warnings_manager import add_warning


DEFAULT_CONFIG: dict[str, dict[str, Any]] = {
    "analysis": {"mode": "full"},
    "project": {"name": "RNA_MassHunter_v2", "output_dir": "output", "log_dir": "logs", "cache_dir": ".cache"},
    "input": {"raw_path": "", "mzml_path": "", "msconvert_path": ""},
    "organism": {"group": "archaea", "species": "", "rule_set": "archaea_general"},
    "sequence": {"name": "target_tRNA", "type": "RNA", "sequence": "", "anticodon": "", "wobble_position": 34},
    "experiment": {"condition_name": "wild_type"},
    "instrument": {"polarity": "negative", "ms1_tolerance_ppm": 10, "ms2_tolerance_da": 0.02},
    "sciex_profile": {
        "enabled": False,
        "path": None,
        "intact_peak_detection": {"enabled": True},
        "intact_mass_comparison": {
            "enabled": True,
            "strict_tolerance_da": 1.0,
            "broad_tolerance_da": 5.0,
        },
        "delta_mass_cluster_audit": {
            "enabled": True,
            "cluster_tolerance_da": 0.5,
            "duplicate_apex_tolerance_da": 0.25,
            "isotope_spacing_da": 1.003355,
            "isotope_spacing_tolerance_da": 0.15,
            "integer_spacing_tolerance_da": 0.15,
            "minimum_cluster_size": 2,
            "max_pair_spacing_da": 200.0,
            "max_pair_rows": 20000,
        },
        "spacing_resolution_audit": {
            "enabled": True,
            "minimum_spacing_sample_count": 20,
            "quantization_tolerance_da": 0.02,
            "distinguishability_margin_factor": 2.0,
            "maximum_spacing_multiple": 10,
        },
        "relation_evidence_quality_audit": {
            "enabled": True,
            "high_error_fraction_threshold": 0.25,
            "low_error_fraction_threshold": 0.75,
            "minimum_recurrent_support_pairs": 2,
            "minimum_interpretable_resolution_margin": 2.0,
        },
        "cross_layer_evidence_reconciliation": {
            "enabled": True,
        },
    },
    "cca_tail": {
        "enabled": True,
        "excludes_cca_candidate_states": ["NONE", "C", "CC", "CCA"],
        "includes_cca_candidate_states": ["CCA", "CC"],
        "default_status": "ASSUMED",
    },
    "cca_processing": {
    "enabled": True,
    },
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
        "intact_reconstruction": {
            "engine": "legacy_cluster",
            "compare_with_legacy": False,
            "min_charge_states_for_reliable": 3,
            "min_charge_states_for_review": 2,
            "require_contiguous_charge_states": True,
            "max_neutral_mass_sd_da": 0.5,
            "max_neutral_mass_range_da": 1.5,
            "max_mass_error_ppm": 20,
            "max_envelope_internal_error_ppm": 20,
            "min_relative_intensity_percent": 0.5,
            "min_relative_envelope_intensity_percent_for_reliable": 1.0,
            "min_relative_envelope_intensity_percent_for_review": 0.1,
            "max_competing_envelopes": 3,
            "comparison_ready_statuses": ["Reliable", "Review"],
            "comparison_ready_tiers": {
                "strict": ["Tier_1_high_quality"],
                "review": ["Tier_1_high_quality", "Tier_2_supported"],
            },
            "quality_tiers": {
                "tier1_min_charge_states": 3,
                "tier1_min_consecutive_charge_states": 3,
                "tier1_min_local_envelope_relative_intensity_percent": 1.0,
                "tier2_min_charge_states": 2,
                "tier2_min_consecutive_charge_states": 2,
                "tier2_min_local_envelope_relative_intensity_percent": 0.1,
            },
            "max_rt_range_min_for_reliable": 0.15,
            "max_rt_range_min_for_review": 0.30,
            "allow_trace_only_reliable": False,
            "search_mode": "untargeted",
            "reference_masses": [],
            "reference_mass_tolerance_ppm": 20,
            "neutral_mass_range": {"enabled": True, "min_da": 20000, "max_da": 30000},
            "target_review_mass_range": {"enabled": False, "min_da": None, "max_da": None},
            "envelope_grouping": {
                "enabled": True,
                "mass_tolerance_da": 1.0,
                "rt_tolerance_min": 0.15,
                "min_shared_peak_fraction": 0.5,
                "min_shared_charge_fraction": 0.5,
                "require_peak_overlap": True,
            },
            "engine_comparison": {
                "mass_tolerance_ppm": 20,
                "rt_tolerance_min": 0.15,
                "min_shared_charge_fraction": 0.5,
                "require_mass_match": True,
            },
            "competitive_assignment": {
                "enabled": True,
                "rt_tolerance_min": 0.15,
                "close_score_margin": 5.0,
                "dry_run": True,
                "min_independent_peak_fraction": 0.5,
                "min_independent_charge_states": 2,
                "allow_shared_peaks_between_selected": False,
                "minimum_score_margin_for_exclusive_selection": 1.0,
                "apply_to_comparison_ready": False,
                "sensitivity_analysis": {
                    "enabled": True,
                    "scenarios": ["strict", "balanced", "sensitive", "permissive"],
                },
                "audit_masses": {"enabled": False, "tolerance_da": 2.0, "masses": []},
                "evidence_score_config_version": "MVP-5.9.8a-v1",
                "score_weights": {
                    "charge_count": 12.0,
                    "consecutive_charge_run": 10.0,
                    "charge_coverage": 8.0,
                    "local_relative_intensity": 5.0,
                    "supporting_scan_count": 3.0,
                    "rt_consistency": 5.0,
                    "extension_support": 2.0,
                    "split_merge_support": 2.0,
                    "internal_error_penalty": 10.0,
                    "neutral_mass_sd_penalty": 5.0,
                    "neutral_mass_range_penalty": 5.0,
                    "rt_range_penalty": 5.0,
                    "peak_sharing_penalty": 12.0,
                    "peak_usage_penalty": 5.0,
                    "charge_gap_penalty": 8.0,
                    "severe_limiting_factor_penalty": 15.0,
                },
            },
            "rt_localized": {
                "enabled": True,
                "rt_window_min": 0.10,
                "rt_step_min": 0.05,
                "min_scans_per_window": 1,
                "peak_aggregation": "max",
                "mz_merge_tolerance_ppm": 10,
                "adjacent_charge_mz_tolerance_ppm": 20,
                "max_charge_gap": 1,
                "min_charge_states": 2,
                "min_consecutive_charge_states": 2,
                "require_consecutive_for_candidate": True,
                "min_local_relative_peak_intensity_percent": 0.1,
                "neutral_mass_estimator": "intensity_weighted_mean",
                "merge_across_windows": {
                    "enabled": True,
                    "mass_tolerance_ppm": 10,
                    "rt_overlap_required": True,
                    "min_shared_charge_fraction": 0.5,
                },
                "charge_extension": {
                    "enabled": True,
                    "max_extension_charges": 2,
                    "weak_peak_tolerance_ppm": 30,
                    "weak_peak_min_local_relative_percent": 0.01,
                    "add_weak_peaks_to_envelope": False,
                },
                "split_envelope_merge": {
                    "enabled": True,
                    "mass_tolerance_ppm": 20,
                    "rt_tolerance_min": 0.10,
                    "max_charge_gap": 1,
                },
                "peak_sharing": {
                    "high_usage_threshold": 5,
                    "max_highly_shared_fraction_for_tier1": 0.25,
                    "max_highly_shared_fraction_for_tier2": 0.50,
                },
            },
            "mass_spectrum_output": {
                "enabled": True,
                "representatives_only": True,
                "comparison_ready_only": False,
                "include_qc_ineligible": True,
                "intensity_method": "total_supporting_intensity",
                "normalize_to_percent": True,
                "bin_width_da": None,
                "minimum_quality_tier": "Tier_3_weak",
                "assignment_filter": "none",
            },
        },
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
    "unknown_modification_search": {
        "enabled": True,
        "mz_tolerance_ppm": 10,
        "max_candidates_per_match": 10,
        "min_peak_tier": ["Major", "Minor"],
        "min_confidence": ["High", "Medium"],
        "include_known_modification_composites": True,
        "candidate_deltas": [
            {"label": "+O (oxidation)", "elements": {"O": 1}},
            {"label": "+S (thiolation)", "elements": {"S": 1}},
            {"label": "O→S substitution (e.g. 2-/4-thiouridine)", "elements": {"O": -1, "S": 1}},
            {"label": "+CH2 (methylation)", "elements": {"C": 1, "H": 2}},
            {"label": "-H2O (dehydration)", "elements": {"H": -2, "O": -1}},
            {"label": "+CH2O2 (formylation-like, +CHO2H)", "elements": {"C": 1, "H": 2, "O": 2}},
        ],
    },
    "p1_sap_dinucleotide": {
        "enabled": True,
        "candidate_generation": {
            "max_modifications_per_side": 3,
            "max_composite_states_per_position": 64,
            "max_candidate_count": 100000,
            "include_normal_phosphate": True,
            "include_phosphorothioate": True,
            "charges": [1],
            "polarity": "auto",
        },
        "search": {"mz_min": 100, "mz_max": 1000, "tolerance_ppm": 10},
        "mass_accuracy": {"strong_ppm": 2, "moderate_ppm": 5, "search_ppm": 10},
        "feature_quality": {
            "min_spectrum_count": 2, "min_profile_point_count": 2,
            "max_rt_gap_min": 0.08, "background_window_rt_min": 0.5,
            "background_mz_tolerance_ppm": 10,
        },
        "isotope": {"enabled": True, "tolerance_ppm": 20, "require_same_scan": True},
        "ms2_provenance": {"enabled": True},
        "targets": [],
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
        "ion_series": ["d", "w", "a", "z"],
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
        optional_section_absent = section == "sciex_profile" and section not in data
        reported_missing = [
            key for key in missing
            if not (
                section == "sciex_profile"
                and key in {"intact_mass_comparison", "delta_mass_cluster_audit", "spacing_resolution_audit", "relation_evidence_quality_audit", "cross_layer_evidence_reconciliation"}
            )
        ]
        if reported_missing and warnings is not None and not optional_section_absent:
            add_warning(warnings, "WARNING", "config", f"Config section '{section}' missing keys filled by defaults.", reported_missing)
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
    analysis_mode = str((config.analysis or {}).get("mode", "full") or "full").lower()
    if analysis_mode not in {"full", "intact_only"}:
        raise ValueError("analysis.mode must be one of: full, intact_only")
    config.analysis["mode"] = analysis_mode
    polarity = str(config.instrument.get("polarity", "")).lower()
    if polarity not in {"negative", "positive"} and warnings is not None:
        add_warning(warnings, "ERROR", "config", "instrument.polarity must be 'negative' or 'positive'.", polarity)

    sciex_profile = config.sciex_profile or {}
    if _as_bool(sciex_profile.get("enabled"), False):
        profile_path = sciex_profile.get("path")
        if profile_path is None:
            raise ValueError("sciex_profile.path is required when sciex_profile.enabled=true")
        if isinstance(profile_path, str) and not profile_path.strip():
            raise ValueError("sciex_profile.path must not be empty when sciex_profile.enabled=true")
        detection = sciex_profile.get("intact_peak_detection")
        if not isinstance(detection, dict):
            raise ValueError("sciex_profile.intact_peak_detection must be a mapping")
        comparison = sciex_profile.get(
            "intact_mass_comparison", DEFAULT_CONFIG["sciex_profile"]["intact_mass_comparison"]
        )
        if not isinstance(comparison, dict):
            raise ValueError("sciex_profile.intact_mass_comparison must be a mapping")
        if _as_bool(comparison.get("enabled"), True):
            strict = comparison.get("strict_tolerance_da", 1.0)
            broad = comparison.get("broad_tolerance_da", 5.0)
            if isinstance(strict, bool) or not isinstance(strict, (int, float)) or not float("-inf") < strict < float("inf") or strict <= 0:
                raise ValueError("sciex_profile.intact_mass_comparison.strict_tolerance_da must be positive")
            if isinstance(broad, bool) or not isinstance(broad, (int, float)) or not float("-inf") < broad < float("inf") or broad < strict:
                raise ValueError("sciex_profile.intact_mass_comparison.broad_tolerance_da must be >= strict_tolerance_da")

        cluster_audit = sciex_profile.get(
            "delta_mass_cluster_audit",
            DEFAULT_CONFIG["sciex_profile"]["delta_mass_cluster_audit"],
        )
        if not isinstance(cluster_audit, dict):
            raise ValueError("sciex_profile.delta_mass_cluster_audit must be a mapping")
        if _as_bool(cluster_audit.get("enabled"), True):
            try:
                DeltaMassClusterParameters.from_mapping(cluster_audit)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid sciex_profile.delta_mass_cluster_audit: {exc}"
                ) from exc

        resolution_audit = sciex_profile.get(
            "spacing_resolution_audit",
            DEFAULT_CONFIG["sciex_profile"]["spacing_resolution_audit"],
        )
        if not isinstance(resolution_audit, dict):
            raise ValueError("sciex_profile.spacing_resolution_audit must be a mapping")
        if _as_bool(resolution_audit.get("enabled"), True):
            try:
                SpacingResolutionParameters.from_mapping(resolution_audit)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid sciex_profile.spacing_resolution_audit: {exc}"
                ) from exc

        evidence_audit = sciex_profile.get(
            "relation_evidence_quality_audit",
            DEFAULT_CONFIG["sciex_profile"]["relation_evidence_quality_audit"],
        )
        if not isinstance(evidence_audit, dict):
            raise ValueError("sciex_profile.relation_evidence_quality_audit must be a mapping")
        if _as_bool(evidence_audit.get("enabled"), True):
            try:
                RelationEvidenceQualityParameters.from_mapping(evidence_audit)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid sciex_profile.relation_evidence_quality_audit: {exc}"
                ) from exc

        cross_layer_audit = sciex_profile.get(
            "cross_layer_evidence_reconciliation",
            DEFAULT_CONFIG["sciex_profile"]["cross_layer_evidence_reconciliation"],
        )
        if not isinstance(cross_layer_audit, dict):
            raise ValueError("sciex_profile.cross_layer_evidence_reconciliation must be a mapping")

    if not config.input.get("mzml_path") and not config.input.get("raw_path") and warnings is not None:
        add_warning(warnings, "WARNING", "config", "input.mzml_path and input.raw_path are empty; sample config mode.")
    if not config.sequence.get("sequence") and warnings is not None:
        add_warning(warnings, "WARNING", "config", "sequence.sequence is empty; theoretical mass will be skipped.")

    digestion_enabled = _as_bool(config.digestion.get("enabled"), True)
    reconstruction_enabled = _as_bool(config.reconstruction.get("enabled"), True)
    if analysis_mode == "intact_only" and not reconstruction_enabled:
        raise ValueError("analysis.mode=intact_only requires intact_reconstruction.enabled=true")
    intact_config = config.reconstruction.get("intact_reconstruction") or {}
    engine = str(intact_config.get("engine") or "legacy_cluster")
    if engine not in {"legacy_cluster", "rt_localized"}:
        raise ValueError("intact_reconstruction.engine must be one of: legacy_cluster, rt_localized")
    rt_localized = intact_config.get("rt_localized") or {}
    aggregation = str(rt_localized.get("peak_aggregation") or "max")
    if aggregation not in {"max", "sum", "mean"}:
        raise ValueError("intact_reconstruction.rt_localized.peak_aggregation must be one of: max, sum, mean")
    estimator = str(rt_localized.get("neutral_mass_estimator") or "intensity_weighted_mean")
    if estimator not in {"unweighted_mean", "intensity_weighted_mean", "median"}:
        raise ValueError("intact_reconstruction.rt_localized.neutral_mass_estimator must be one of: unweighted_mean, intensity_weighted_mean, median")
    spectrum_config = intact_config.get("mass_spectrum_output") or {}
    intensity_method = str(spectrum_config.get("intensity_method") or "total_supporting_intensity")
    if intensity_method not in {"total_supporting_intensity", "mean_supporting_intensity", "max_supporting_intensity"}:
        raise ValueError("intact_reconstruction.mass_spectrum_output.intensity_method must be one of: total_supporting_intensity, mean_supporting_intensity, max_supporting_intensity")
    assignment_filter = str(spectrum_config.get("assignment_filter") or "none").lower()
    if assignment_filter not in {"none", "strict", "review", "balanced_selected", "all"}:
        raise ValueError("intact_reconstruction.mass_spectrum_output.assignment_filter must be one of: none, strict, review, balanced_selected, all")
    minimum_quality_tier = str(spectrum_config.get("minimum_quality_tier") or "Tier_3_weak")
    if minimum_quality_tier not in {"Tier_1_high_quality", "Tier_2_supported", "Tier_3_weak", "all"}:
        raise ValueError("intact_reconstruction.mass_spectrum_output.minimum_quality_tier must be one of: Tier_1_high_quality, Tier_2_supported, Tier_3_weak, all")
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
    for section_name in ("project", "input", "sciex_profile"):
        section = getattr(config, section_name)
        for key, value in list(section.items()):
            if not value or not isinstance(value, str):
                continue
            if (
                key.endswith("_dir")
                or key.endswith("_path")
                or (section_name == "sciex_profile" and key == "path")
            ):
                path = Path(value)
                if not path.is_absolute():
                    section[key] = str(root / path)
    return config
