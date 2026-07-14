import argparse
from pathlib import Path
import time
import tracemalloc

from rna_masshunter.config import load_config, validate_config, resolve_paths
from rna_masshunter.composite_modification_audit import (
    append_composite_diagnostics, build_composite_modification_audit,
)
from rna_masshunter.audit_policy import (
    AUDIT_LEVELS, AUDIT_STATUS_COLUMNS, AuditPolicy, append_audit_level_diagnostics,
    audit_status_row, sheet_category,
)
from rna_masshunter.biological_context import biological_context_priority_rows
from rna_masshunter.biological_position_prior import evaluate_biological_position_priors, load_position_prior_rules
from rna_masshunter.conversion import prepare_input_file
from rna_masshunter.digestion import digest_sequence
from rna_masshunter.excel_report import write_excel_report
from rna_masshunter.evidence_ranking import build_ambiguity_groups, build_modification_evidence_ranking
from rna_masshunter.intact_reconstruction import reconstruct_intact_masses
from rna_masshunter.logging_utils import setup_logger
from rna_masshunter.masses import calculate_unmodified_rna_mass, load_base_masses
from rna_masshunter.modification_search import known_modification_candidate_rows, search_known_modifications, summarize_known_modification_candidates
from rna_masshunter.modifications import load_modifications, validate_modifications
from rna_masshunter.ms1_mapping import map_fragments_to_ms1_peaks
from rna_masshunter.ms1_match_truncation_audit import (
    append_diagnostic_shadow_columns,
    append_top_shadow_columns,
    build_ms1_truncation_audit,
)
from rna_masshunter.ms1_selection_strategy_audit import (
    append_selection_diagnostic_columns,
    append_top_selection_columns,
    build_ms1_selection_strategy_audit,
)
from rna_masshunter.ms1_top50_dedup_audit import (
    append_top50_diagnostic_columns,
    append_top50_shadow_columns,
    build_ms1_top50_dedup_audit,
)
from rna_masshunter.ms1_cross_fragment_ambiguity import (
    append_crossfrag_diagnostic_columns,
    append_crossfrag_top_columns,
    build_ms1_cross_fragment_ambiguity_audit,
)
from rna_masshunter.ms2_annotation import annotate_ms2
from rna_masshunter.ms2_ambiguous_peak_audit import (
    build_ambiguous_peak_audit, build_ambiguity_summary, build_ambiguity_diagnostics,
)
from rna_masshunter.ms2_identity_evidence import build_ms2_modification_identity
from rna_masshunter.ms2_unmatched_audit import build_unmatched_ion_summary
from rna_masshunter.ms2_zero_intensity_audit import (
    build_zero_intensity_audit, update_top50_affected as update_zero_top50_affected,
)
from rna_masshunter.ms2_effective_ambiguity import (
    build_effective_ambiguity, update_top50_affected as update_effective_top50_affected,
)
from rna_masshunter.position_mapper import build_position_map
from rna_masshunter.review_dashboard import build_review_dashboard_results
from rna_masshunter.mzml_diagnostics import run_mzml_diagnostics
from rna_masshunter.pathway_loader import load_pathways, validate_pathways
from rna_masshunter.p1_annotation import build_p1_optional_results, is_p1_enabled
from rna_masshunter.peak_filtering import classify_peak_tiers
from rna_masshunter.peak_picking import extract_ms1_peaks
from rna_masshunter.rule_loader import load_rule_set, validate_rule_set
from rna_masshunter.startup_check import run_startup_check
from rna_masshunter.warnings_manager import add_warning


def _as_bool(value, default: bool = True) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _record_workflow_step(
    rows,
    analysis_mode: str,
    step_name: str,
    status: str,
    enabled_by_config,
    executed: bool,
    skip_reason: str = "",
    output_sheets: str = "",
    notes: str = "",
) -> None:
    rows.append(
        {
            "Analysis_Mode": analysis_mode,
            "Step_Name": step_name,
            "Step_Status": status,
            "Enabled_By_Config": enabled_by_config,
            "Executed": executed,
            "Skip_Reason": skip_reason,
            "Output_Sheets": output_sheets,
            "Notes": notes,
        }
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run RNA_MassHunter analysis.")
    parser.add_argument(
        "--config",
        default=None,
        help="Path to a YAML config file. Defaults to config.yaml beside main.py.",
    )
    parser.add_argument(
        "--audit-level",
        choices=AUDIT_LEVELS,
        default="full",
        help="Shadow audit output level: standard, audit summaries, or full detail (default: full).",
    )
    return parser.parse_args(argv)


def resolve_config_path(project_root: Path, configured_path: str | None) -> Path:
    if configured_path is None:
        return project_root / "config.yaml"
    path = Path(configured_path).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    return path.resolve()


def main(argv: list[str] | None = None) -> None:
    project_root = Path(__file__).resolve().parent
    args = parse_args(argv)
    audit_policy = AuditPolicy.from_level(args.audit_level)
    audit_status_rows = []
    config_path = resolve_config_path(project_root, args.config)
    logger = setup_logger(project_root / "logs")
    warnings = []
    logger.info("RNA_MassHunter_v2 MVP-5 started")

    config = load_config(config_path, warnings=warnings)
    validate_config(config, warnings=warnings)
    config = resolve_paths(config, project_root)
    analysis_mode = str((config.analysis or {}).get("mode") or "full")
    intact_only = analysis_mode == "intact_only"
    workflow_rows = []

    _record_workflow_step(workflow_rows, analysis_mode, "config_load_and_validation", "executed", True, True, output_sheets="Input_parameters")
    run_startup_check(project_root, config, logger, warnings, config_path=config_path)
    _record_workflow_step(workflow_rows, analysis_mode, "startup_check", "executed", True, True, output_sheets="Warnings")

    modifications = []
    rule_set = {}
    pathways = []
    if intact_only:
        _record_workflow_step(workflow_rows, analysis_mode, "modification_dictionary", "skipped_by_analysis_mode", True, False, "intact_only skips modification resources")
        _record_workflow_step(workflow_rows, analysis_mode, "organism_rule_set", "skipped_by_analysis_mode", True, False, "intact_only skips biological rule resources")
        _record_workflow_step(workflow_rows, analysis_mode, "pathway_resources", "skipped_by_analysis_mode", True, False, "intact_only skips pathway resources")
    else:
        modifications = load_modifications(project_root / "data" / "modifications.yaml", warnings=warnings)
        validate_modifications(modifications, warnings=warnings)
        _record_workflow_step(workflow_rows, analysis_mode, "modification_dictionary", "executed", True, True, notes=f"entries={len(modifications)}")

        rule_set = load_rule_set(
            project_root / "data" / "rule_sets",
            config.organism.get("rule_set", "methanosarcina_acetivorans"),
            warnings=warnings,
        )
        validate_rule_set(rule_set, warnings=warnings)
        _record_workflow_step(workflow_rows, analysis_mode, "organism_rule_set", "executed", True, True, notes=str(rule_set.get("id") or rule_set.get("name") or ""))

        pathways = load_pathways(project_root / "data" / "pathways", warnings=warnings)
        validate_pathways(pathways, warnings=warnings)
        _record_workflow_step(workflow_rows, analysis_mode, "pathway_resources", "executed", True, True, notes=f"files={len(pathways)}")

    mzml_path = prepare_input_file(config, logger, warnings)
    diagnostics = {}
    peaks = []
    tier_result = classify_peak_tiers([], config.peak_filtering, warnings=warnings)
    intact_results = []
    charge_state_peaks = []
    intact_engine_artifacts = {}

    if mzml_path:
        diagnostics = run_mzml_diagnostics(mzml_path, logger, warnings)
        _record_workflow_step(workflow_rows, analysis_mode, "mzML_diagnostics", "executed", True, True, output_sheets="mzML_diagnostics")
        peaks = extract_ms1_peaks(mzml_path, config.reconstruction, warnings=warnings)
        tier_result = classify_peak_tiers(peaks, config.peak_filtering, warnings=warnings)
        _record_workflow_step(workflow_rows, analysis_mode, "MS1_peak_extraction", "executed", True, True, output_sheets="Charge_state_peaks", notes=f"peaks={len(peaks)}")
    else:
        add_warning(warnings, "WARNING", "main", "No mzML input was provided; report will be written without mzML-derived results.")
        diagnostics = {"Warnings": "No mzML input was provided."}
        _record_workflow_step(workflow_rows, analysis_mode, "mzML_diagnostics", "unavailable", True, False, "no_mzML_input", output_sheets="mzML_diagnostics")
        _record_workflow_step(workflow_rows, analysis_mode, "MS1_peak_extraction", "unavailable", True, False, "no_mzML_input")

    base_masses = load_base_masses(project_root / "data" / "base_masses.yaml", warnings=warnings)
    theoretical_mass = None
    theoretical_fragments = []
    fragment_ms1_matches = []
    known_modification_candidates = []
    known_modification_summary = []
    ranking_rows = []
    ms1_audit_context = {}
    sequence = (config.sequence.get("sequence", "") or "").upper().replace("T", "U")
    digestion_enabled = _as_bool(config.digestion.get("enabled"), True)
    fragment_mapping_enabled = _as_bool(config.fragment_mapping.get("enabled"), True)

    if sequence:
        theoretical_mass = calculate_unmodified_rna_mass(sequence, base_masses, warnings=warnings)
        _record_workflow_step(workflow_rows, analysis_mode, "unmodified_theoretical_mass_annotation", "executed", True, True, notes=f"mass={theoretical_mass}")
        position_map = build_position_map(sequence, config.sequence.get("wobble_position", 34))
        if intact_only:
            _record_workflow_step(workflow_rows, analysis_mode, "digestion", "skipped_by_analysis_mode", digestion_enabled, False, "intact_only skips digestion", output_sheets="")
        elif digestion_enabled:
            theoretical_fragments = digest_sequence(
                target_id=config.sequence.get("name", "target_tRNA"),
                sequence=sequence,
                position_map=position_map,
                config=config,
                base_masses=base_masses,
                warnings=warnings,
            )
            _record_workflow_step(workflow_rows, analysis_mode, "digestion", "executed", True, True, output_sheets="Theoretical_fragments", notes=f"fragments={len(theoretical_fragments)}")
        else:
            add_warning(
                warnings,
                "INFO",
                "main",
                "digestion.enabled is false; theoretical fragments and Fragment_MS1 mapping were skipped. Intact mass reconstruction can still run.",
                {"target_id": config.sequence.get("name", "target_tRNA"), "sequence_length": len(sequence)},
            )
            _record_workflow_step(workflow_rows, analysis_mode, "digestion", "disabled_by_config", False, False, "digestion.enabled=false")
    else:
        add_warning(warnings, "WARNING", "masses", "config.sequence.sequence is empty; theoretical mass and theoretical fragments were not calculated.")
        _record_workflow_step(workflow_rows, analysis_mode, "unmodified_theoretical_mass_annotation", "unavailable", True, False, "empty_sequence")
        _record_workflow_step(workflow_rows, analysis_mode, "digestion", "unavailable", digestion_enabled, False, "empty_sequence")

    if intact_only:
        _record_workflow_step(workflow_rows, analysis_mode, "fragment_MS1_mapping", "skipped_by_analysis_mode", fragment_mapping_enabled, False, "intact_only skips fragment mapping")
    elif digestion_enabled and theoretical_fragments and peaks and fragment_mapping_enabled:
        fragment_ms1_matches = map_fragments_to_ms1_peaks(
            theoretical_fragments,
            peaks,
            config,
            warnings=warnings,
            audit_context=ms1_audit_context,
        )
        _record_workflow_step(workflow_rows, analysis_mode, "fragment_MS1_mapping", "executed", True, True, output_sheets="Fragment_MS1_matches; Fragment_MS1_filtered; Fragment_MS1_summary", notes=f"matches={len(fragment_ms1_matches)}")
    elif digestion_enabled and fragment_mapping_enabled:
        fragment_ms1_matches = map_fragments_to_ms1_peaks(
            theoretical_fragments,
            peaks,
            config,
            warnings=warnings,
            audit_context=ms1_audit_context,
        )
        _record_workflow_step(workflow_rows, analysis_mode, "fragment_MS1_mapping", "executed", True, True, output_sheets="Fragment_MS1_matches; Fragment_MS1_filtered; Fragment_MS1_summary", notes=f"matches={len(fragment_ms1_matches)}")
    else:
        _record_workflow_step(workflow_rows, analysis_mode, "fragment_MS1_mapping", "disabled_by_config", fragment_mapping_enabled, False, "digestion or fragment_mapping disabled")

    reconstruction_enabled = _as_bool(config.reconstruction.get("enabled"), True)
    if reconstruction_enabled:
        intact_results, charge_state_peaks, intact_engine_artifacts = reconstruct_intact_masses(
            tier_result,
            config.reconstruction,
            config.instrument,
            theoretical_mass,
            warnings=warnings,
        )
        _record_workflow_step(
            workflow_rows,
            analysis_mode,
            "intact_reconstruction_QC_grouping_comparison",
            "executed",
            True,
            True,
            output_sheets="Intact_mass_reconstruction; Charge_state_peaks; Intact_Reconstruction_QC; Intact_Reconstruction_Diag; Intact_Envelope_Groups; Intact_Comparison_Candidates; Target_Review_Candidates; Reconstructed_Mass_Spectrum; RT_Engine_QC_Summary; Intact_Competition_Groups; Intact_Competition_Scores",
            notes=f"candidates={len(intact_results)}",
        )
    else:
        _record_workflow_step(workflow_rows, analysis_mode, "intact_reconstruction_QC_grouping_comparison", "disabled_by_config", False, False, "reconstruction.enabled=false")

    if intact_only:
        _record_workflow_step(workflow_rows, analysis_mode, "known_modification_search", "skipped_by_analysis_mode", _as_bool(config.modification_search.get("enabled"), True), False, "intact_only skips modification search")
    elif _as_bool(config.modification_search.get("enabled"), True):
        known_modification_candidates = search_known_modifications(
            fragment_ms1_matches=fragment_ms1_matches,
            intact_results=intact_results,
            modifications=modifications,
            config=config,
            warnings=warnings,
        )
        known_modification_summary = summarize_known_modification_candidates(known_modification_candidates)
        _record_workflow_step(workflow_rows, analysis_mode, "known_modification_search", "executed", True, True, output_sheets="Known_Modification_Candidates; Known_Modification_Summary", notes=f"candidates={len(known_modification_candidates)}")
    else:
        _record_workflow_step(workflow_rows, analysis_mode, "known_modification_search", "disabled_by_config", False, False, "modification_search.enabled=false")

    optional_results = {"Workflow_Summary": workflow_rows}
    if intact_engine_artifacts.get("rt_envelope_diagnostics") is not None:
        optional_results["RT_Envelope_Diagnostics"] = intact_engine_artifacts.get("rt_envelope_diagnostics", [])
    if intact_engine_artifacts.get("missing_charge_diagnostics") is not None:
        optional_results["Missing_Charge_Diagnostics"] = intact_engine_artifacts.get("missing_charge_diagnostics", [])
    if intact_engine_artifacts.get("engine_comparison") is not None:
        optional_results["Intact_Engine_Comparison"] = intact_engine_artifacts.get("engine_comparison", [])
    if intact_only:
        skipped_optional_steps = [
            ("MS2_annotation", _as_bool(config.ms2_annotation.get("enabled"), True)),
            ("P1_annotation", is_p1_enabled(config)),
            ("modification_evidence_ranking", _as_bool(config.modification_evidence_ranking.get("enabled"), True)),
            ("biological_context", _as_bool(config.biological_context.get("enabled"), True)),
            ("review_dashboard", True),
        ]
        for step_name, enabled in skipped_optional_steps:
            _record_workflow_step(workflow_rows, analysis_mode, step_name, "skipped_by_analysis_mode", enabled, False, "intact_only skips downstream modification review")
    else:
        optional_results.update(annotate_ms2(
            mzml_path=str(mzml_path) if mzml_path else None,
            theoretical_fragments=theoretical_fragments,
            config=config,
            base_masses=base_masses,
            modifications=modifications,
            warnings=warnings,
        ))
        _record_workflow_step(workflow_rows, analysis_mode, "MS2_annotation", "executed" if _as_bool(config.ms2_annotation.get("enabled"), True) else "disabled_by_config", _as_bool(config.ms2_annotation.get("enabled"), True), _as_bool(config.ms2_annotation.get("enabled"), True), output_sheets="MS2_*")
        if is_p1_enabled(config):
            optional_results.update(build_p1_optional_results(config, tier_result, base_masses, modifications))
            _record_workflow_step(workflow_rows, analysis_mode, "P1_annotation", "executed", True, True, output_sheets="P1_*")
        else:
            _record_workflow_step(workflow_rows, analysis_mode, "P1_annotation", "disabled_by_config", False, False, "p1_annotation.enabled=false")
        if _as_bool(config.modification_evidence_ranking.get("enable_ambiguity_grouping"), True):
            optional_results["Modification_Ambiguity_Groups"] = build_ambiguity_groups(
                optional_results.get("MS2_Modification_Localization_Evidence", []),
                optional_results.get("MS2_Modified_Ion_Matches", []),
            )
        else:
            optional_results["Modification_Ambiguity_Groups"] = []
        ranking_rows, ranking_summary = build_modification_evidence_ranking(
            config=config,
            modifications=modifications,
            theoretical_fragments=theoretical_fragments,
            fragment_ms1_matches=fragment_ms1_matches,
            known_candidates=known_modification_candidates,
            ms2_results=optional_results,
            rule_set=rule_set,
            pathways=pathways,
        )
        position_prior_rules = load_position_prior_rules(project_root / "data" / "modification_position_priors.yaml")
        ranking_rows, position_prior_rows, plausibility_rows, position_diagnostics = evaluate_biological_position_priors(
            config, ranking_rows, modifications, position_prior_rules
        )
        ranking_rows, identity_rows, identity_assignment_rows = build_ms2_modification_identity(
            ranking_rows,
            optional_results.get("MS2_Modified_Ion_Matches", []),
            optional_results.get("MS2_Modification_Localization_Evidence", []),
            optional_results.get("Modification_Ambiguity_Groups", []),
            enabled=_as_bool(config.ms2_annotation.get("enabled"), True),
            return_assignments=True,
        )
        optional_results["Modification_Evidence_Summary"] = ranking_summary
        optional_results["Modification_Evidence_Ranking"] = ranking_rows
        optional_results["Modification_Position_Priors"] = position_prior_rows
        optional_results["MS2_Biological_Plausibility"] = plausibility_rows
        optional_results["Biological_Prior_Diagnostics"] = position_diagnostics
        optional_results["MS2_Modification_Identity"] = identity_rows
        optional_results["MS2_Identity_Peak_Assignments"] = identity_assignment_rows
        optional_results["MS2_Unmatched_Ion_Summary"] = build_unmatched_ion_summary(
            ranking_rows, optional_results.get("MS2_Unmatched_Ion_Audit", []),
        )
        ambiguity_context = optional_results.pop("_MS2_Ambiguous_Audit_Context", {})
        zero_context = optional_results.pop("_MS2_Zero_Intensity_Audit_Context", {})
        zero_summary = []
        effective_summary = []
        if audit_policy.run_shadow_audits:
            tracemalloc.start()
            audit_started = time.perf_counter()
            ambiguous_clusters, ambiguous_peak_details = build_ambiguous_peak_audit(
                ambiguity_context.get("spectra", []),
                optional_results.get("MS2_Unmatched_Ion_Audit", []),
                optional_results.get("MS2_Modified_Theoretical_Ions", []),
                ranking_rows, identity_assignment_rows,
                enabled=_as_bool(config.ms2_annotation.get("enabled"), True),
            )
            ambiguity_summary = build_ambiguity_summary(ranking_rows, ambiguous_clusters, ambiguous_peak_details)
            ambiguity_diagnostics = build_ambiguity_diagnostics(
                ambiguous_clusters, ambiguous_peak_details, ambiguity_summary,
                enabled=_as_bool(config.ms2_annotation.get("enabled"), True),
            )[0]
            ambiguity_runtime = time.perf_counter() - audit_started
            _, ambiguity_peak_bytes = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            audit_status_rows.append(audit_status_row(
                "MS2_ambiguous_peak", "MS2", audit_policy, True, True, True,
                ambiguity_runtime, ambiguity_peak_bytes / (1024 * 1024),
            ))
            optional_results["MS2_Ambiguous_Peak_Clusters"] = ambiguous_clusters
            optional_results["MS2_Ambiguous_Peak_Detail"] = ambiguous_peak_details
            optional_results["MS2_Ambiguity_Summary"] = ambiguity_summary
            unmatched_diagnostics = optional_results.get("MS2_Unmatched_Ion_Diagnostics") or [{}]
            unmatched_diagnostics[0].update(ambiguity_diagnostics)

            zero_enabled = _as_bool(config.ms2_annotation.get("Enable_MS2_Zero_Intensity_Audit"), True)
            nonzero_simulation = _as_bool(config.ms2_annotation.get("Enable_Nonzero_Shadow_Simulation"), True)
            report_row_limit = int(getattr(config, "reporting", {}).get("max_excel_rows_per_sheet", 100000) or 100000)
            configured_zero_limit = int(config.ms2_annotation.get("max_zero_intensity_detail_rows", report_row_limit) or report_row_limit)
            max_zero_detail_rows = min(configured_zero_limit, report_row_limit)
            tracemalloc.start()
            audit_started = time.perf_counter()
            zero_spectra, zero_detail, zero_summary, zero_candidates, zero_diagnostics = build_zero_intensity_audit(
                zero_context, ranking_rows, optional_results.get("MS2_Ion_Matches", []),
                optional_results.get("MS2_Modified_Ion_Matches", []), identity_assignment_rows,
                optional_results.get("MS2_Modification_Localization_Evidence", []),
                ambiguous_clusters, ambiguous_peak_details, enabled=zero_enabled,
                nonzero_simulation=nonzero_simulation, max_detail_rows=max_zero_detail_rows,
            )
            zero_runtime = time.perf_counter() - audit_started
            _, zero_peak_bytes = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            audit_status_rows.append(audit_status_row(
                "MS2_zero_intensity", "MS2", audit_policy, True, True, True,
                zero_runtime, zero_peak_bytes / (1024 * 1024),
            ))
            optional_results["MS2_Zero_Intensity_Spectra"] = zero_spectra
            optional_results["MS2_Zero_Intensity_Detail"] = zero_detail
            optional_results["MS2_Zero_Intensity_Summary"] = zero_summary
            optional_results["_MS2_Zero_Intensity_Candidate_Summary"] = zero_candidates

            effective_enabled = _as_bool(config.ms2_annotation.get("Enable_MS2_Effective_Ambiguity_Audit"), True)
            tracemalloc.start()
            audit_started = time.perf_counter()
            effective_clusters, effective_detail, effective_summary, effective_candidates, effective_diagnostics = build_effective_ambiguity(
                ambiguous_clusters, ambiguous_peak_details, ambiguity_summary,
                optional_results.get("MS2_Ion_Matches", []), optional_results.get("MS2_Modified_Ion_Matches", []),
                identity_assignment_rows, optional_results.get("MS2_Modification_Localization_Evidence", []),
                ranking_rows, zero_context, enabled=effective_enabled, max_detail_rows=report_row_limit,
            )
            effective_runtime = time.perf_counter() - audit_started
            _, effective_peak_bytes = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            audit_status_rows.append(audit_status_row(
                "MS2_effective_ambiguity", "MS2", audit_policy, True, True, True,
                effective_runtime, effective_peak_bytes / (1024 * 1024),
            ))
            optional_results["MS2_Effective_Ambiguity"] = effective_clusters
            optional_results["MS2_Effective_Ambig_Detail"] = effective_detail
            optional_results["MS2_Effective_Ambig_Summary"] = effective_summary
            optional_results["_MS2_Effective_Ambiguity_Candidate_Summary"] = effective_candidates
            unmatched_diagnostics[0].update(zero_diagnostics[0])
            unmatched_diagnostics[0].update(effective_diagnostics[0])
            optional_results["MS2_Unmatched_Ion_Diagnostics"] = unmatched_diagnostics
        _record_workflow_step(workflow_rows, analysis_mode, "modification_evidence_ranking", "executed", _as_bool(config.modification_evidence_ranking.get("enabled"), True), True, output_sheets="Modification_Evidence_Summary; Modification_Evidence_Ranking; Modification_Ambiguity_Groups", notes=f"ranked={len(ranking_rows)}")
        optional_results["Biological_Context_Priorities"] = biological_context_priority_rows(config)
        optional_results["Context_Supported_Candidates"] = [
            row for row in ranking_rows if float(row.get("Biological_Context_Score") or 0.0) > 0
        ]
        _record_workflow_step(workflow_rows, analysis_mode, "biological_context", "executed" if _as_bool(config.biological_context.get("enabled"), True) else "disabled_by_config", _as_bool(config.biological_context.get("enabled"), True), True, output_sheets="Biological_Context_Priorities; Context_Supported_Candidates")
        review_results = build_review_dashboard_results(optional_results, config)
        optional_results.update(review_results)
        if audit_policy.run_shadow_audits:
            update_zero_top50_affected(zero_summary, review_results.get("Top_Modification_Candidates", []))
            update_effective_top50_affected(effective_summary, review_results.get("Top_Modification_Candidates", []))
        optional_results.pop("_MS2_Zero_Intensity_Candidate_Summary", None)
        optional_results.pop("_MS2_Effective_Ambiguity_Candidate_Summary", None)
        _record_workflow_step(workflow_rows, analysis_mode, "review_dashboard", "executed", True, True, output_sheets="Review_*")

    ms1_audit_enabled = _as_bool(config.fragment_mapping.get("Enable_MS1_Truncation_Audit"), True)
    ms1_shadow_enabled = _as_bool(config.fragment_mapping.get("Enable_MS1_Expanded_Shadow_Simulation"), True)
    if audit_policy.run_shadow_audits and not intact_only and ms1_audit_enabled and ms1_shadow_enabled and ms1_audit_context:
        tracemalloc.start()
        ms1_shadow_started = time.perf_counter()
        ms1_audit = build_ms1_truncation_audit(
            context=ms1_audit_context,
            config=config,
            modifications=modifications,
            intact_results=intact_results,
            baseline_matches=fragment_ms1_matches,
            baseline_candidates=known_modification_candidates,
            baseline_ranking=ranking_rows,
            ms2_results=optional_results,
            rule_set=rule_set,
            pathways=pathways,
        )
        ms1_shadow_seconds = time.perf_counter() - ms1_shadow_started
        _, ms1_shadow_peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        ms1_audit["summary"]["Shadow_Additional_Time_Seconds"] = ms1_shadow_seconds
        ms1_audit["summary"]["Shadow_Peak_Tracked_Memory_MiB"] = ms1_shadow_peak_bytes / (1024 * 1024)
        ms1_audit["summary_rows"][0]["Shadow_Additional_Time_Seconds"] = ms1_shadow_seconds
        ms1_audit["summary_rows"][0]["Shadow_Peak_Tracked_Memory_MiB"] = ms1_shadow_peak_bytes / (1024 * 1024)
        optional_results["MS1_Truncation_Audit"] = ms1_audit["audit_rows"]
        optional_results["MS1_Truncation_Detail"] = ms1_audit["detail_rows"]
        optional_results["MS1_Truncation_Summary"] = ms1_audit["summary_rows"]
        optional_results["Top_Modification_Candidates"] = append_top_shadow_columns(
            optional_results.get("Top_Modification_Candidates"), ms1_audit
        )
        optional_results["MS2_Unmatched_Ion_Diagnostics"] = append_diagnostic_shadow_columns(
            optional_results.get("MS2_Unmatched_Ion_Diagnostics"), ms1_audit
        )
        audit_status_rows.append(audit_status_row(
            "MS1_match_truncation", "MS1", audit_policy, True, True, True,
            ms1_shadow_seconds, ms1_shadow_peak_bytes / (1024 * 1024),
        ))

    selection_audit_enabled = _as_bool(config.fragment_mapping.get("Enable_MS1_Selection_Strategy_Audit"), True)
    selection_apply_formal = _as_bool(config.fragment_mapping.get("Apply_MS1_Selection_Strategy_To_Formal_Result"), False)
    if audit_policy.run_shadow_audits and not intact_only and selection_audit_enabled and not selection_apply_formal and ms1_audit_context:
        tracemalloc.start()
        selection_started = time.perf_counter()
        selection_audit = build_ms1_selection_strategy_audit(
            context=ms1_audit_context, config=config, modifications=modifications,
            intact_results=intact_results, baseline_matches=fragment_ms1_matches,
            baseline_candidates=known_modification_candidates, baseline_ranking=ranking_rows,
            ms2_results=optional_results, rule_set=rule_set, pathways=pathways,
        )
        selection_seconds = time.perf_counter() - selection_started
        _, selection_peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        selection_audit["summary"]["Strategy_Audit_Additional_Time_Seconds"] = selection_seconds
        selection_audit["summary"]["Strategy_Audit_Peak_Tracked_Memory_MiB"] = selection_peak_bytes / (1024 * 1024)
        selection_audit["summary_rows"][0]["Strategy_Audit_Additional_Time_Seconds"] = selection_seconds
        selection_audit["summary_rows"][0]["Strategy_Audit_Peak_Tracked_Memory_MiB"] = selection_peak_bytes / (1024 * 1024)
        optional_results["MS1_Selection_Strategy"] = selection_audit["strategy_rows"]
        optional_results["MS1_Selection_Detail"] = selection_audit["detail_rows"]
        optional_results["MS1_Selection_Summary"] = selection_audit["summary_rows"]
        optional_results["Top_Modification_Candidates"] = append_top_selection_columns(
            optional_results.get("Top_Modification_Candidates"), selection_audit
        )
        optional_results["MS2_Unmatched_Ion_Diagnostics"] = append_selection_diagnostic_columns(
            optional_results.get("MS2_Unmatched_Ion_Diagnostics"), selection_audit
        )
        audit_status_rows.append(audit_status_row(
            "MS1_selection_strategy", "MS1", audit_policy, True, True, True,
            selection_seconds, selection_peak_bytes / (1024 * 1024),
        ))

    top50_audit = None
    top50_audit_enabled = _as_bool(config.fragment_mapping.get("Enable_MS1_Tier_Top50_Full_Shadow"), True)
    dedup_audit_enabled = _as_bool(config.fragment_mapping.get("Enable_MS1_Physical_Peak_Dedup_Audit"), True)
    top50_apply_formal = _as_bool(config.fragment_mapping.get("Apply_MS1_Tier_Top50_To_Formal_Result"), False)
    dedup_apply_formal = _as_bool(config.fragment_mapping.get("Apply_MS1_Dedup_To_Formal_Result"), False)
    if audit_policy.run_shadow_audits and not intact_only and top50_audit_enabled and dedup_audit_enabled and not top50_apply_formal and not dedup_apply_formal and ms1_audit_context:
        tracemalloc.start()
        top50_audit = build_ms1_top50_dedup_audit(
            context=ms1_audit_context, config=config, modifications=modifications,
            intact_results=intact_results, baseline_matches=fragment_ms1_matches,
            baseline_candidates=known_modification_candidates, baseline_ranking=ranking_rows,
            ms2_results=optional_results, rule_set=rule_set, pathways=pathways,
        )
        _, top50_peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        top50_audit["summary"]["Audit_Peak_Tracked_Memory_MiB"] = top50_peak_bytes / (1024 * 1024)
        top50_audit["summary_rows"][0]["Audit_Peak_Tracked_Memory_MiB"] = top50_peak_bytes / (1024 * 1024)
        optional_results["MS1_Top50_Shadow"] = top50_audit["top50_rows"]
        optional_results["MS1_Peak_Dedup_Detail"] = top50_audit["detail_rows"]
        optional_results["MS1_Top50_Dedup_Summary"] = top50_audit["summary_rows"]
        optional_results["Top_Modification_Candidates"] = append_top50_shadow_columns(
            optional_results.get("Top_Modification_Candidates"), top50_audit
        )
        optional_results["MS2_Unmatched_Ion_Diagnostics"] = append_top50_diagnostic_columns(
            optional_results.get("MS2_Unmatched_Ion_Diagnostics"), top50_audit
        )
        audit_status_rows.append(audit_status_row(
            "MS1_top50_physical_peak", "MS1", audit_policy, True, True, True,
            float(top50_audit["summary"].get("Top50_Shadow_Additional_Time_Seconds") or 0),
            top50_peak_bytes / (1024 * 1024),
        ))

    crossfrag_enabled = _as_bool(config.fragment_mapping.get("Enable_MS1_Cross_Fragment_Ambiguity_Audit"), True)
    crossfrag_apply_formal = _as_bool(config.fragment_mapping.get("Apply_MS1_Cross_Fragment_Ambiguity_To_Formal_Result"), False)
    if audit_policy.run_shadow_audits and not intact_only and crossfrag_enabled and not crossfrag_apply_formal and ms1_audit_context:
        tracemalloc.start()
        crossfrag_audit = build_ms1_cross_fragment_ambiguity_audit(
            context=ms1_audit_context, config=config, modifications=modifications,
            intact_results=intact_results, baseline_candidates=known_modification_candidates,
            baseline_ranking=ranking_rows, ms2_results=optional_results,
            rule_set=rule_set, pathways=pathways, top50_audit=top50_audit,
        )
        _, crossfrag_peak_bytes = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        crossfrag_audit["summary"]["Audit_Peak_Tracked_Memory_MiB"] = crossfrag_peak_bytes / (1024 * 1024)
        crossfrag_audit["summary_rows"][0]["Audit_Peak_Tracked_Memory_MiB"] = crossfrag_peak_bytes / (1024 * 1024)
        optional_results["MS1_CrossFrag_Ambiguity"] = crossfrag_audit["ambiguity_rows"]
        optional_results["MS1_CrossFrag_Detail"] = crossfrag_audit["detail_rows"]
        optional_results["MS1_CrossFrag_Summary"] = crossfrag_audit["summary_rows"]
        optional_results["Top_Modification_Candidates"] = append_crossfrag_top_columns(
            optional_results.get("Top_Modification_Candidates"), crossfrag_audit
        )
        optional_results["MS2_Unmatched_Ion_Diagnostics"] = append_crossfrag_diagnostic_columns(
            optional_results.get("MS2_Unmatched_Ion_Diagnostics"), crossfrag_audit
        )
        audit_status_rows.append(audit_status_row(
            "MS1_cross_fragment", "MS1", audit_policy, True, True, True,
            float(crossfrag_audit["summary"].get("Ambiguity_Grouping_Time_Seconds") or 0)
            + float(crossfrag_audit["summary"].get("Shadow_Weighting_Time_Seconds") or 0),
            crossfrag_peak_bytes / (1024 * 1024),
        ))

    composite_audit = None
    if audit_policy.run_shadow_audits and sequence and not intact_only:
        composite_audit = build_composite_modification_audit(
            project_root, sequence, modifications, base_masses, audit_mode=audit_policy.level,
        )
        optional_results.update(composite_audit.sheets)
        audit_status_rows.extend([
            audit_status_row("Composite modification constraints", "Modification", audit_policy, True, True, audit_policy.include_detail, composite_audit.runtime_seconds, composite_audit.peak_memory_mb),
            audit_status_row("Backbone modification candidates", "Backbone", audit_policy, True, True, audit_policy.include_detail, 0.0, composite_audit.peak_memory_mb),
            audit_status_row("Cleavage blocking constraints", "Digestion", audit_policy, True, True, audit_policy.include_detail, 0.0, composite_audit.peak_memory_mb),
        ])

    audit_specs = (
        ("MS2_ambiguous_peak", "MS2"), ("MS2_zero_intensity", "MS2"),
        ("MS2_effective_ambiguity", "MS2"), ("MS1_match_truncation", "MS1"),
        ("MS1_selection_strategy", "MS1"), ("MS1_top50_physical_peak", "MS1"),
        ("MS1_cross_fragment", "MS1"),
        ("Composite modification constraints", "Modification"),
        ("Backbone modification candidates", "Backbone"),
        ("Cleavage blocking constraints", "Digestion"),
    )
    recorded = {row["Audit_Name"] for row in audit_status_rows}
    for audit_name, category in audit_specs:
        if audit_name not in recorded:
            audit_status_rows.append(audit_status_row(
                audit_name, category, audit_policy, False, True, True,
                reason=f"not executed for audit_level={audit_policy.level} or unavailable/disabled input",
            ))
    if audit_policy.include_summary:
        optional_results["Audit_Status"] = [
            {column: row.get(column, "") for column in AUDIT_STATUS_COLUMNS}
            for row in audit_status_rows
        ]
    shadow_sheet_count = sum(
        1 for name in optional_results
        if (sheet_category(name) or "").startswith("AUDIT_")
        and audit_policy.includes_category(sheet_category(name))
    )
    optional_results["MS2_Unmatched_Ion_Diagnostics"] = append_audit_level_diagnostics(
        optional_results.get("MS2_Unmatched_Ion_Diagnostics"), audit_policy,
        audit_status_rows, shadow_sheet_count,
    )
    optional_results["MS2_Unmatched_Ion_Diagnostics"] = append_composite_diagnostics(
        optional_results.get("MS2_Unmatched_Ion_Diagnostics"), composite_audit,
    )
    report_path = write_excel_report(
        output_dir=Path(config.project["output_dir"]),
        config=config,
        diagnostics=diagnostics,
        intact_results=intact_results,
        charge_state_peaks=charge_state_peaks,
        warnings=warnings,
        modifications=modifications,
        rule_set=rule_set,
        pathways=pathways,
        theoretical_fragments=theoretical_fragments,
        fragment_ms1_matches=fragment_ms1_matches,
        known_modification_candidates=known_modification_candidate_rows(known_modification_candidates),
        known_modification_summary=known_modification_summary,
        optional_results=optional_results,
        audit_policy=audit_policy,
    )
    logger.info("Excel report written: %s", report_path)
    logger.info("RNA_MassHunter_v2 MVP-5 finished")

if __name__ == "__main__":
    main()
