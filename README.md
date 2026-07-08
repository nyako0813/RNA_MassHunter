# RNA_MassHunter_v2

RNA_MassHunter_v2 is an MVP workflow for RNA/tRNA LC-MS analysis. Current MVP-3.1 functionality includes YAML configuration loading, mzML diagnostics, MS1 peak extraction, intact mass reconstruction, RNase digestion fragment generation, MS1 fragment matching, filtered fragment summaries, and Excel output.

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
- Writes Excel output with `Run_summary`, `Input_parameters`, `mzML_diagnostics`, `Intact_mass_reconstruction`, `Charge_state_peaks`, `Theoretical_fragments`, `Fragment_MS1_matches`, `Fragment_MS1_filtered`, `Fragment_MS1_summary`, and `Warnings`.

## Planned Later Features

- modification candidate ranking
- unknown modification candidate generation
- MS2 annotation
- mass balance check
- manual comparison

## Data Notes

`data/modifications.yaml` is treated as the user-confirmed PDF-derived modification dictionary when present. MVP-3.1 loads and validates it; full modification candidate search is planned for MVP-4.

`data/base_masses.yaml` contains Mongo Oligo-compatible placeholder masses. Confirm these values against the final Mongo Oligo settings before quantitative final analysis.
