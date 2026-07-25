# SCIEX Cross-Layer Bundle Operational Workflow

This runbook covers the production four-layer handoff to the bundle-only cross-layer runner. The runner reconciles already-produced evidence; it never opens raw mzML and never replaces a layer producer.

## 1. Fixed layer set and source plan

Prepare one production result for each layer in this order: `FULL`, `T1`, `P1AP_MS1`, and `P1AP_MS2`. Record the raw source path, source SHA256, RNA identity, sample identity, digest/condition, producer commit, and intended independence group before execution. `04 new T1.mzML` and `05 old T1.mzML` are distinct samples; a bundle from one is never a substitute for the other. When P1/AP MS1 and MS2 originate from the same run, give them the same shared-source and independence groups.

## 2. Generate through the production producer chain

Run every layer producer in a separate process with its production settings. Pass the exact returned production result object directly to `export_layer_evidence_bundle(...)` with the matching layer, source, RNA, experiment, and producer commit. Do not construct result objects from report rows, edit JSON by hand, reuse Antigravity output, or feed an older serializer bundle to the runner. The export API rejects the wrong result class and empty evidence. The general cross-layer CLI does not invoke producers.

The four bundle JSON objects must retain their schema/serializer versions, layer, producer function/module, optional result key, source SHA256, RNA and experiment identities, non-empty record validation, safeguards, and `validation.status: PASSED`. Raw spectrum arrays and binary payloads are prohibited.

## 3. Validate and persist immediately

After each producer finishes:

1. Export its JSON to a new temporary path.
2. Load it with `load_layer_evidence_bundle(...)`, supplying the source and expected identity where available.
3. Confirm `validation.status` is `PASSED`, `non_empty` is true, and `record_count` is positive.
4. Compute the file SHA256.
5. Copy it to a repository-external persistent directory using a new filename.
6. Recompute SHA256 on the persistent copy and require an exact match.
7. Record path, size, SHA256, record count, source provenance, peak RSS, and elapsed time before starting the next layer.

Do not advance while a bundle exists only in `/tmp`. Never report an empty bundle as successful.

## 4. Persistent directory layout

Use immutable, run-specific paths, for example:

```text
valid_resume/2026-07-25_run-001/
  bundle_FULL.json
  bundle_T1.json
  bundle_P1AP_MS1.json
  bundle_P1AP_MS2.json
  checksums.sha256
```

Do not place production bundles, raw data, or generated workbooks under version control. Retain the source manifest beside the persistent copies when local data-governance policy permits it.

## 5. Validate the four-bundle set

Run the CLI dry run before requesting output:

```bash
PYTHONPATH=. python -m rna_masshunter.sciex_cross_layer_bundle_cli \
  --full-bundle /persistent/bundle_FULL.json \
  --t1-bundle /persistent/bundle_T1.json \
  --p1ap-ms1-bundle /persistent/bundle_P1AP_MS1.json \
  --p1ap-ms2-bundle /persistent/bundle_P1AP_MS2.json \
  --dry-run
```

A passing dry run validates, restores, checks compatibility, and aggregates without writing. Confirm bundle-level provenance is verified, inspect all compatibility warnings, and preserve the FULL node-provenance warning as a known limitation rather than suppressing it.

## 6. Create outputs once

Choose paths that do not exist, then run the command documented in [README](../README.md#sciex-cross-layer-evidence-bundles). At least aggregate JSON or Excel is required outside dry-run; summary JSON is optional. The CLI stages outputs and refuses overwrite, input/output collision, duplicate output paths, and raw mzML inputs.

The aggregate JSON is the machine-readable record. The six Excel sheets are `XL_Nodes`, `XL_Edges`, `XL_Hypotheses`, `XL_Layer_Summary`, `XL_Consensus`, and `XL_Next_Evidence`. All remain shadow-only.

## 7. Interpret safeguards and uncertainty

P1/AP MS1 and MS2 from a shared raw source are non-independent support. RNA identity and sample identity are separate checks. An ambiguous consensus is a successful reconciliation outcome. `LOW` confidence means the available evidence is insufficient. Neither state licenses exact chemical identity, structural isomer identity, nucleotide position, atom localization, reaction order, formal score changes, ranking changes, filtering, or final-consensus propagation.

## 8. Failure recovery

Use the CLI exit category to locate the failure: argument/configuration (2), bundle validation (3), compatibility (4), output (5), or aggregation (6). Keep all valid persistent inputs unchanged. Correct metadata only by regenerating through the production producer and serializer; never patch bundle JSON or migrate an old Antigravity bundle. For an output collision, select a new output path. Re-run dry-run after any regenerated input, compare its new SHA256 with the run manifest, and create final outputs only after all four inputs pass together.

If a process is interrupted, remove only known incomplete staging files according to local operator policy. Do not delete untracked raw measurements or valid bundles and do not weaken thresholds to obtain a pass.

## 9. Windows and WSL path handling

Quote every path containing spaces. From WSL, a persistent Windows directory such as `C:\Users\nyako\Documents\RNA_MassHunter_E2E_WIP\valid_resume` is addressed as `/mnt/c/Users/nyako/Documents/RNA_MassHunter_E2E_WIP/valid_resume`. Run Python and the repository from one environment so path identity checks are consistent. Do not mix a Windows path string into a WSL Python invocation.

Example:

```bash
cd /home/nyako/projects/RNA_MassHunter
PYTHONPATH=. .venv/bin/python -m rna_masshunter.sciex_cross_layer_bundle_cli \
  --config "/mnt/c/Users/nyako/Documents/RNA_MassHunter_E2E_WIP/run-config.yaml" \
  --dry-run
```

## 10. Audit checklist

Before sign-off, preserve the four persistent bundle paths and SHA256 values, source provenance, validation counts, producer commit, RNA/sample identities, shared-source groups, compatibility report, warnings, aggregate counts, consensus/confidence, output hashes, elapsed time, and memory observations. Confirm no raw arrays appear in JSON, every safeguard remains shadow-only, all Excel sheet names are at most 31 characters, repository `config.yaml` is unchanged, and no bundle, raw file, output, log, cache, backup, or virtual environment is staged.
