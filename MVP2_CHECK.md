# RNA_MassHunter MVP-2 Check

## Date

2026-07-08

## Environment

- OS: WSL Ubuntu
- Python: 3.12.3
- Project path: /home/nyako/projects/RNA_MassHunter
- GitHub repository: nyako0813/RNA_MassHunter

## Confirmed before MVP-2

- MVP-1 was already implemented and pushed.
- `python main.py` finished successfully on the new PC/WSL environment.
- `MVP1_CHECK.md` was present.
- Working tree was clean before implementation.

## MVP-2 implemented features

- RNase digestion basic logic.
- Enzyme rules for:
  - RNase_T1
  - RNase_A
  - RNase_T2
  - Nuclease_P1
  - Benzonase
  - U_specific_RNase
- Theoretical RNA fragment generation from `config.sequence.sequence`.
- Missed cleavage handling.
- Minimum fragment length filtering.
- Alkaline phosphatase-aware terminal forms:
  - dephosphorylated
  - residual_phosphate
  - cyclic_phosphate
- Fragment unmodified mass calculation.
- Sequence start/end positions.
- Standard tRNA position mapping using `wobble_position`.
- Excel `Theoretical_fragments` sheet.
- Index link to `Theoretical_fragments`.
- Back-to-Index link from `Theoretical_fragments`.
- Run_summary `Theoretical fragments` count.

## Test 1: sample config with empty sequence

Confirmed:

- `python main.py` finished successfully.
- Excel report was created.
- `Theoretical_fragments` sheet was created.
- `Theoretical_fragments` sheet had headers and 0 data rows.
- Existing MVP-1 sheets were preserved.
- Index links and Back-to-Index links worked.

## Test 2: short test sequence

Temporary test sequence:

```text
GGGAAAUUUGGCUAGC
```

Confirmed:

- `python main.py` finished successfully.
- Excel report was created.
- `Theoretical_fragments` sheet contained 24 data rows.
- RNase_T1 cleavage after G was reflected.
- Missed cleavages 0 and 1 were present.
- Terminal forms were present:
  - dephosphorylated
  - residual_phosphate
  - cyclic_phosphate
- Index link to `Theoretical_fragments` worked.
- Back-to-Index link from `Theoretical_fragments` worked.
- Run_summary `Theoretical fragments` count was 24.

## Post-test cleanup

- Temporary test sequence was removed from `config.yaml`.
- `config.sequence.sequence` was restored to an empty string.
- MVP-2 digestion and alkaline phosphatase configuration remained in `config.yaml`.

## Scope intentionally deferred

- MS1 fragment mapping was not implemented in MVP-2.
- Modification candidate generation was not implemented in MVP-2.
- Unknown modification search was not implemented in MVP-2.
- MS2 annotation was not implemented in MVP-2.
- Co-purified tRNA detailed analysis was not implemented in MVP-2.
- Mass balance check was not implemented in MVP-2.

## Notes

- `data/modifications.yaml` was preserved.
- `data/rule_sets` was preserved.
- `data/pathways` was preserved.
- Commit and push were performed after MVP-2 implementation before this check note was added.
