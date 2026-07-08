from pathlib import Path

from rna_masshunter.config import load_config, validate_config, resolve_paths
from rna_masshunter.conversion import prepare_input_file
from rna_masshunter.digestion import digest_sequence
from rna_masshunter.excel_report import write_excel_report
from rna_masshunter.intact_reconstruction import reconstruct_intact_masses
from rna_masshunter.logging_utils import setup_logger
from rna_masshunter.masses import calculate_unmodified_rna_mass, load_base_masses
from rna_masshunter.modifications import load_modifications, validate_modifications
from rna_masshunter.ms1_mapping import map_fragments_to_ms1_peaks
from rna_masshunter.position_mapper import build_position_map
from rna_masshunter.mzml_diagnostics import run_mzml_diagnostics
from rna_masshunter.pathway_loader import load_pathways, validate_pathways
from rna_masshunter.peak_filtering import classify_peak_tiers
from rna_masshunter.peak_picking import extract_ms1_peaks
from rna_masshunter.rule_loader import load_rule_set, validate_rule_set
from rna_masshunter.startup_check import run_startup_check
from rna_masshunter.warnings_manager import add_warning


def main() -> None:
    project_root = Path(__file__).resolve().parent
    logger = setup_logger(project_root / "logs")
    warnings = []
    logger.info("RNA_MassHunter_v2 MVP-3 started")

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
    sequence = (config.sequence.get("sequence", "") or "").upper().replace("T", "U")
    if sequence:
        theoretical_mass = calculate_unmodified_rna_mass(sequence, base_masses, warnings=warnings)
        position_map = build_position_map(sequence, config.sequence.get("wobble_position", 34))
        if config.digestion.get("enabled", True):
            theoretical_fragments = digest_sequence(
                target_id=config.sequence.get("name", "target_tRNA"),
                sequence=sequence,
                position_map=position_map,
                config=config,
                base_masses=base_masses,
                warnings=warnings,
            )
    else:
        add_warning(warnings, "WARNING", "masses", "config.sequence.sequence is empty; theoretical mass and theoretical fragments were not calculated.")

    if theoretical_fragments and peaks and config.fragment_mapping.get("enabled", True):
        fragment_ms1_matches = map_fragments_to_ms1_peaks(
            theoretical_fragments,
            peaks,
            config,
            warnings=warnings,
        )
    elif config.fragment_mapping.get("enabled", True):
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
    )
    logger.info("Excel report written: %s", report_path)
    logger.info("RNA_MassHunter_v2 MVP-3 finished")


if __name__ == "__main__":
    main()
