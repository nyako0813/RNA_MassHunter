from pathlib import Path

from rna_masshunter.config import load_config, validate_config, resolve_paths
from rna_masshunter.biological_context import biological_context_priority_rows
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
from rna_masshunter.ms2_annotation import annotate_ms2
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


def main() -> None:
    project_root = Path(__file__).resolve().parent
    logger = setup_logger(project_root / "logs")
    warnings = []
    logger.info("RNA_MassHunter_v2 MVP-5 started")

    config = load_config(project_root / "config.yaml", warnings=warnings)
    validate_config(config, warnings=warnings)
    config = resolve_paths(config, project_root)
    analysis_mode = str((config.analysis or {}).get("mode") or "full")
    intact_only = analysis_mode == "intact_only"
    workflow_rows = []

    _record_workflow_step(workflow_rows, analysis_mode, "config_load_and_validation", "executed", True, True, output_sheets="Input_parameters")
    run_startup_check(project_root, config, logger, warnings)
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
        )
        _record_workflow_step(workflow_rows, analysis_mode, "fragment_MS1_mapping", "executed", True, True, output_sheets="Fragment_MS1_matches; Fragment_MS1_filtered; Fragment_MS1_summary", notes=f"matches={len(fragment_ms1_matches)}")
    elif digestion_enabled and fragment_mapping_enabled:
        fragment_ms1_matches = map_fragments_to_ms1_peaks(
            theoretical_fragments,
            peaks,
            config,
            warnings=warnings,
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
            output_sheets="Intact_mass_reconstruction; Charge_state_peaks; Intact_Reconstruction_QC; Intact_Reconstruction_Diag; Intact_Envelope_Groups; Intact_Comparison_Candidates; Target_Review_Candidates; Reconstructed_Mass_Spectrum",
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
        optional_results["Modification_Evidence_Summary"] = ranking_summary
        optional_results["Modification_Evidence_Ranking"] = ranking_rows
        _record_workflow_step(workflow_rows, analysis_mode, "modification_evidence_ranking", "executed", _as_bool(config.modification_evidence_ranking.get("enabled"), True), True, output_sheets="Modification_Evidence_Summary; Modification_Evidence_Ranking; Modification_Ambiguity_Groups", notes=f"ranked={len(ranking_rows)}")
        optional_results["Biological_Context_Priorities"] = biological_context_priority_rows(config)
        optional_results["Context_Supported_Candidates"] = [
            row for row in ranking_rows if float(row.get("Biological_Context_Score") or 0.0) > 0
        ]
        _record_workflow_step(workflow_rows, analysis_mode, "biological_context", "executed" if _as_bool(config.biological_context.get("enabled"), True) else "disabled_by_config", _as_bool(config.biological_context.get("enabled"), True), True, output_sheets="Biological_Context_Priorities; Context_Supported_Candidates")
        optional_results.update(build_review_dashboard_results(optional_results, config))
        _record_workflow_step(workflow_rows, analysis_mode, "review_dashboard", "executed", True, True, output_sheets="Review_*")

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
    )
    logger.info("Excel report written: %s", report_path)
    logger.info("RNA_MassHunter_v2 MVP-5 finished")

if __name__ == "__main__":
    main()
