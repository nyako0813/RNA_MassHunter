from pathlib import Path

from rna_masshunter.config import load_config, validate_config, resolve_paths
from rna_masshunter.conversion import prepare_input_file
from rna_masshunter.digestion import digest_sequence
from rna_masshunter.excel_report import write_excel_report
from rna_masshunter.intact_reconstruction import reconstruct_intact_masses
from rna_masshunter.logging_utils import setup_logger
from rna_masshunter.masses import calculate_unmodified_rna_mass, load_base_masses
from rna_masshunter.modification_search import known_modification_candidate_rows, search_known_modifications, summarize_known_modification_candidates
from rna_masshunter.modifications import load_modifications, validate_modifications
from rna_masshunter.ms1_mapping import map_fragments_to_ms1_peaks
from rna_masshunter.ms2_annotation import annotate_ms2
from rna_masshunter.position_mapper import build_position_map
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


def main() -> None:
    project_root = Path(__file__).resolve().parent
    logger = setup_logger(project_root / "logs")
    warnings = []
    logger.info("RNA_MassHunter_v2 MVP-5 started")

    config = load_config(project_root / "config.yaml", warnings=warnings)
    validate_config(config, warnings=warnings)
    config = resolve_paths(config, project_root)

    run_startup_check(project_root, config, logger, warnings)

    modifications = load_modifications(project_root / "data" / "modifications.yaml", warnings=warnings)
    validate_modifications(modifications, warnings=warnings)

    rule_set = load_rule_set(
        project_root / "data" / "rule_sets",
        config.organism.get("rule_set", "methanosarcina_acetivorans"),
        warnings=warnings,
    )
    validate_rule_set(rule_set, warnings=warnings)

    pathways = load_pathways(project_root / "data" / "pathways", warnings=warnings)
    validate_pathways(pathways, warnings=warnings)

    mzml_path = prepare_input_file(config, logger, warnings)
    diagnostics = {}
    peaks = []
    tier_result = classify_peak_tiers([], config.peak_filtering, warnings=warnings)
    intact_results = []
    charge_state_peaks = []

    if mzml_path:
        diagnostics = run_mzml_diagnostics(mzml_path, logger, warnings)
        peaks = extract_ms1_peaks(mzml_path, config.reconstruction, warnings=warnings)
        tier_result = classify_peak_tiers(peaks, config.peak_filtering, warnings=warnings)
    else:
        add_warning(warnings, "WARNING", "main", "No mzML input was provided; MVP-3 report will be written without mzML-derived results.")
        diagnostics = {"Warnings": "No mzML input was provided."}

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
        position_map = build_position_map(sequence, config.sequence.get("wobble_position", 34))
        if digestion_enabled:
            theoretical_fragments = digest_sequence(
                target_id=config.sequence.get("name", "target_tRNA"),
                sequence=sequence,
                position_map=position_map,
                config=config,
                base_masses=base_masses,
                warnings=warnings,
            )
        else:
            add_warning(
                warnings,
                "INFO",
                "main",
                "digestion.enabled is false; theoretical fragments and Fragment_MS1 mapping were skipped. Intact mass reconstruction can still run.",
                {"target_id": config.sequence.get("name", "target_tRNA"), "sequence_length": len(sequence)},
            )
    else:
        add_warning(warnings, "WARNING", "masses", "config.sequence.sequence is empty; theoretical mass and theoretical fragments were not calculated.")

    if digestion_enabled and theoretical_fragments and peaks and fragment_mapping_enabled:
        fragment_ms1_matches = map_fragments_to_ms1_peaks(
            theoretical_fragments,
            peaks,
            config,
            warnings=warnings,
        )
    elif digestion_enabled and fragment_mapping_enabled:
        fragment_ms1_matches = map_fragments_to_ms1_peaks(
            theoretical_fragments,
            peaks,
            config,
            warnings=warnings,
        )

    if config.reconstruction.get("enabled", True):
        intact_results, charge_state_peaks = reconstruct_intact_masses(
            tier_result,
            config.reconstruction,
            config.instrument,
            theoretical_mass,
            warnings=warnings,
        )

    if _as_bool(config.modification_search.get("enabled"), True):
        known_modification_candidates = search_known_modifications(
            fragment_ms1_matches=fragment_ms1_matches,
            intact_results=intact_results,
            modifications=modifications,
            config=config,
            warnings=warnings,
        )
        known_modification_summary = summarize_known_modification_candidates(known_modification_candidates)

    optional_results = {}
    optional_results.update(annotate_ms2(
        mzml_path=str(mzml_path) if mzml_path else None,
        theoretical_fragments=theoretical_fragments,
        config=config,
        base_masses=base_masses,
        warnings=warnings,
    ))
    if is_p1_enabled(config):
        optional_results.update(build_p1_optional_results(config, tier_result, base_masses, modifications))

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