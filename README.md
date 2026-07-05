# RNA_MassHunter_v2

RNA_MassHunter_v2 is an MVP-1 workflow for RNA/tRNA LC-MS/MS analysis. This first version focuses on stable startup, YAML loading, mzML diagnostics, MS1 peak extraction, tiered peak reporting, intact mass reconstruction, and Excel output.

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
- Fill `sequence.sequence` when theoretical unmodified RNA mass is needed.
- The default sample config intentionally leaves input paths and sequence empty, so `python main.py` can start and produce an empty report with warnings.

## Run

```bash
python main.py
```

Reports are written to `output/`, logs to `logs/`, and cache/checkpoint files to `.cache/`.

## MVP-1 Features

- Reads `config.yaml`.
- Reads and validates `data/modifications.yaml`.
- Reads rule sets from `data/rule_sets/*.yaml`, including `inherits` with child rule override by `position_rule.id`.
- Reads pathways from `data/pathways/*.yaml`.
- Runs startup checks and logs warnings/errors.
- Confirms mzML input and runs mzML diagnostics using pyteomics.
- Extracts centroid MS1 peaks from mzML.
- Classifies peaks as Major, Minor, Trace, or below reporting threshold.
- Reconstructs intact masses from charge states.
- Writes Excel output with `Run_summary`, `Input_parameters`, `mzML_diagnostics`, `Intact_mass_reconstruction`, `Charge_state_peaks`, and `Warnings`.

## Not Implemented In MVP-1

The following modules are intentionally stubs for MVP-2 or later:

- RNase digestion
- AP-aware theoretical fragment generation
- co-purified tRNA candidate search
- MS1 fragment mapping
- modification candidate ranking
- unknown modification candidate generation
- MS2 annotation
- mass balance check
- manual comparison

## Data Notes

`data/modifications.yaml` is treated as the user-confirmed PDF-derived modification dictionary when present. MVP-1 loads and validates it but does not perform modification candidate search.

`data/base_masses.yaml` contains Mongo Oligo-compatible placeholder masses. Confirm these values against the final Mongo Oligo settings before quantitative final analysis.
