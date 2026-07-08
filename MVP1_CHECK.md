# RNA_MassHunter MVP-1 Check

## Date

2026-07-05

## Environment

- OS: WSL Ubuntu 24.04
- Python: 3.12.3
- Project path: /home/nyako/projects/RNA_MassHunter
- GitHub repository: nyako0813/RNA_MassHunter

## Confirmed

- MVP-1 code was integrated into the WSL-side repository.
- `python main.py` finished successfully.
- `config.yaml` was loaded.
- `data/modifications.yaml` was found and loaded.
- `data/rule_sets/` was found.
- `data/pathways/` was found.
- `data/base_masses.yaml` was found.
- `data/genotypes/knockout_template.yaml` was found.
- Startup check passed.
- Excel report was created in `output/`.
- Log file was created in `logs/`.
- GitHub push was completed.

## Confirmed Excel sheets

- Index
- Run_summary
- Input_parameters
- mzML_diagnostics
- Intact_mass_reconstruction
- Charge_state_peaks
- Warnings

## Fixed during MVP-1 check

- Added `lxml` to requirements because pyteomics depends on it.
- Fixed rule_set `inherits` handling for list-style inheritance.
- Fixed pathway loader to accept `pathway_id`.
- Added Excel Index sheet.
- Added links from Index to each sheet.
- Added Back to Index links from each sheet.
- Changed Index sheet so the sheet name itself is the link.
- Fixed Run_summary so `organism.rule_set` is displayed.

## Current expected warnings

The current `config.yaml` is still a sample config, so the following warnings are expected:

- `input.mzml_path` and `input.raw_path` are empty.
- `sequence.sequence` is empty.
- No mzML or raw input is configured.
- Theoretical mass is not calculated.

## Important notes

- `data/modifications.yaml` is the PDF-derived curated modification dictionary.
- `base_masses.yaml` currently contains Mongo Oligo-compatible placeholder values and must be checked before final quantitative use.
- MVP-1 does not perform RNase digestion.
- MVP-1 does not perform modification candidate search.
- MVP-1 does not perform unknown modification search.
- MVP-1 does not perform co-purified tRNA search.
- MVP-1 does not perform MS2 annotation.
- MVP-1 does not perform mass balance check.

## Next steps

1. Prepare an actual mzML file.
2. Edit `config.yaml` with:
   - mzML path
   - target tRNA sequence
   - polarity
   - charge range
   - optional RT range
3. Run `python main.py`.
4. Check `mzML_diagnostics`.
5. Check `Intact_mass_reconstruction`.
6. Use the results to decide MVP-2 priorities.
