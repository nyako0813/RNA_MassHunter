# RNA_MassHunter_v2

RNA_MassHunter_v2 is an MVP workflow for RNA/tRNA LC-MS analysis. Current MVP-5 functionality includes YAML configuration loading, mzML diagnostics, MS1 peak extraction, intact mass reconstruction, RNase digestion fragment generation, MS1 fragment matching, filtered fragment summaries, known modification candidate search, P1 peak annotation, MS2 theoretical ion annotation, and Excel output.

## Setup

Create and activate a virtual environment, then install dependencies.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## Configure

Edit `config.yaml` before real analysis.

- Set `input.mzml_path` to an existing mzML file, or set `input.raw_path`.
- WIFF/WIFF2 conversion requires `input.msconvert_path`.
- Fill `sequence.sequence` when theoretical unmodified RNA mass or fragment masses are needed.
- Set `organism.group`, `organism.species`, and `organism.rule_set` for the target RNA.
- The default sample config intentionally leaves input paths and sequence empty, so `python main.py` can start and produce an empty report with warnings.

## Analysis modes

### Intact/full-length reconstruction only

Use this when the experiment is not digested and you only want full-length charge-state reconstruction.

```yaml
digestion:
  enabled: false

fragment_mapping:
  enabled: false

reconstruction:
  enabled: true
```

In this mode, RNA_MassHunter still calculates the theoretical full-length unmodified mass when `sequence.sequence` is present, then performs intact mass reconstruction from MS1 charge states. Theoretical fragment generation and `Fragment_MS1_*` matching are skipped.

### RNase digestion + MS1 fragment matching

Use this when the sample was digested, for example with RNase A.

```yaml
digestion:
  enabled: true
  enzyme: RNase_A
  missed_cleavages: 1

fragment_mapping:
  enabled: true
```

Supported built-in enzyme names currently include `RNase_T1`, `RNase_A`, `RNase_T2`, `Nuclease_P1`, `Benzonase`, and `U_specific_RNase`.

### T1 digestion fragment analysis

Use this when RNase_T1 digestion fragments are the primary analysis target and full-length intact reconstruction is not needed.

```yaml
reconstruction:
  enabled: false

digestion:
  enabled: true
  enzyme: RNase_T1

fragment_mapping:
  enabled: true

modification_search:
  enabled: true
```

In this mode, the Excel report focuses on theoretical fragments, MS1 fragment matches, filtered fragment summaries, and known modification candidates. Full-length intact reconstruction sheets are omitted when `reconstruction.enabled` is false.

### MVP-5 MS2 annotation

MS2 annotation is an evidence table, not de novo sequencing. RNA_MassHunter reads MS2 spectra from mzML and matches observed MS2 peaks against theoretical RNA fragment ions generated from the known sequence and digestion fragments. The initial MVP-5 ion series is unmodified `c`/`y` ions. Neutral losses, base losses, and strict modification-site localization are deferred to later MVPs.

```yaml
ms2_annotation:
  enabled: true
  mz_tolerance_ppm: 20
  min_peak_intensity: 0
  max_peaks_per_spectrum: 500
  use_theoretical_fragments: true
```

When an mzML file has no MS2 spectra, the run still succeeds and the Excel report records `Total_MS2_Spectra = 0` with annotation skipped/no MS2 status.

## Alkaline phosphatase setting

Record whether the sample was treated with alkaline phosphatase.

No AP treatment:

```yaml
alkaline_phosphatase:
  enabled: false
```

AP treatment, assumed complete:

```yaml
alkaline_phosphatase:
  enabled: true
  assume_complete: true
```

AP treatment, not assumed complete:

```yaml
alkaline_phosphatase:
  enabled: true
  assume_complete: false
  allow_residual_phosphate: true
  allow_cyclic_phosphate: true
```

When AP is enabled and complete, only dephosphorylated terminal forms are generated for fragments. When AP is enabled but incomplete, residual phosphate and cyclic phosphate terminal forms can also be retained.

## Run

```bash
python main.py
```

Reports are written to `output/`, logs to `logs/`, and cache/checkpoint files to `.cache/`.

## Current Features

- Reads `config.yaml`.
- Reads and validates `data/modifications.yaml`.
- Reads rule sets from `data/rule_sets/*.yaml`, including `inherits` with child rule override by `position_rule.id`.
- Reads pathways from `data/pathways/*.yaml`.
- Runs startup checks and logs warnings/errors.
- Confirms mzML input and runs mzML diagnostics using pyteomics.
- Extracts centroid MS1 peaks from mzML.
- Classifies peaks as Major, Minor, Trace, or below reporting threshold.
- Reconstructs intact masses from charge states.
- Generates RNase theoretical fragments when digestion is enabled.
- Matches theoretical fragments to MS1 peaks when fragment mapping is enabled.
- Searches known modification candidates from `data/modifications.yaml` by comparing observed fragment neutral-mass shifts against curated modification mass shifts.
- Annotates MS2 spectra against theoretical unmodified c/y fragment ions when `ms2_annotation.enabled` is true.
- Writes Excel output with applicable sheets such as `Run_summary`, `Input_parameters`, `mzML_diagnostics`, `Intact_mass_reconstruction`, `Charge_state_peaks`, `Theoretical_fragments`, `Fragment_MS1_matches`, `Fragment_MS1_filtered`, `Fragment_MS1_summary`, `Known_Modification_Candidates`, `Known_Modification_Summary`, `P1_*`, `MS2_Summary`, `MS2_Spectra`, `MS2_Theoretical_Ions`, `MS2_Ion_Matches`, `MS2_Unmatched_Peaks`, and `Warnings`.

## Planned Later Features

- modified-fragment combination ranking beyond single known shifts
- unknown modification candidate generation
- mass balance check
- manual comparison

## Data Notes

`data/modifications.yaml` is treated as the user-confirmed PDF-derived modification dictionary when present. MVP-4 loads and validates it, then reports known single-modification candidates when observed fragment masses can be explained by dictionary mass shifts. Unknown modification search is still deferred.

MS1 mass differences alone cannot distinguish isobaric modifications such as pseudouridine. By default, modifications with `mass_shift_from_unmodified = 0` are excluded from `Known_Modification_Candidates`. Set `modification_search.include_isobaric_modifications: true` only when you explicitly want to report those mass-neutral candidates. RNA_MassHunter does not use special priority logic for standard position 34; Excel reports use sequence start/end positions for fragment review.

`data/base_masses.yaml` contains Mongo Oligo-compatible placeholder masses. Confirm these values against the final Mongo Oligo settings before quantitative final analysis.
