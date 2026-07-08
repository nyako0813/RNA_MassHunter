# RNA_MassHunter MVP-3 Check

## Date

2026-07-08

## Environment

- OS: WSL Ubuntu
- Python: 3.12.3
- Project path: /home/nyako/projects/RNA_MassHunter
- GitHub repository: nyako0813/RNA_MassHunter

## Confirmed before MVP-3

- MVP-1 was implemented and pushed.
- MVP-2 theoretical fragment generation was implemented and pushed.
- `MVP1_CHECK.md` and `MVP2_CHECK.md` were present.
- Working tree was clean before MVP-3 implementation.

## MVP-3 implemented features

- Added `fragment_mapping` settings to `config.yaml`.
- Added default `fragment_mapping` settings to `rna_masshunter/config.py`.
- Added `FragmentMS1Match` dataclass.
- Implemented MS1 fragment mapping in `rna_masshunter/ms1_mapping.py`.
- Added theoretical m/z calculation from neutral fragment mass and charge.
- Supported negative mode and positive mode m/z calculation.
- Added ppm error calculation.
- Added charge-state search from configured `min_charge` to `max_charge`.
- Matched unmodified theoretical fragments against MS1 peaks using ppm tolerance.
- Added simple confidence labels:
  - High
  - Medium
  - Low
- Added `Fragment_MS1_matches` Excel sheet.
- Added Index link to `Fragment_MS1_matches`.
- Added Back-to-Index link from `Fragment_MS1_matches`.
- Added Run_summary `Fragment MS1 matches` count.

## Test 1: sample config with empty sequence and no mzML

Confirmed:

- `python main.py` finished successfully.
- Excel report was created.
- `Theoretical_fragments` sheet was created.
- `Fragment_MS1_matches` sheet was created.
- Both sheets were empty when no sequence and no mzML were provided.
- Existing MVP-1 and MVP-2 sheets were preserved.
- Index links and Back-to-Index links were preserved.

## Expected behavior without mzML

When mzML is not configured:

- MS1-derived results are empty.
- `Fragment_MS1_matches` is empty.
- The program still creates a complete Excel report.
- Warnings are recorded, but this is expected for sample config mode.

## Files changed in MVP-3

- `config.yaml`
- `main.py`
- `rna_masshunter/config.py`
- `rna_masshunter/excel_report.py`
- `rna_masshunter/models.py`
- `rna_masshunter/ms1_mapping.py`

## Notes

- `data/modifications.yaml` was not modified.
- `data/rule_sets/` was not modified.
- `data/pathways/` was not modified.
- MVP-3 maps unmodified theoretical fragments to MS1 peaks.
- MVP-3 does not yet perform modified fragment candidate generation.
- MVP-3 does not yet perform unknown modification search.
- MVP-3 does not yet perform MS2 annotation.
- MVP-3 does not yet perform mass balance check.
- MVP-3 does not yet perform detailed co-purified tRNA analysis.

## Next step

The next recommended step is to test MVP-3 with a real mzML file and a real tRNA sequence.

After real-data testing, MVP-4 should focus on modification candidate generation by comparing observed fragment masses against the modification dictionary and position/rule-set constraints.
