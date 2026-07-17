# RNA_MassHunter_v2

RNA_MassHunter_v2 is an MVP workflow for RNA/tRNA LC-MS analysis. Current MVP-5.1 functionality includes YAML configuration loading, mzML diagnostics, MS1 peak extraction, intact mass reconstruction, RNase digestion fragment generation, MS1 fragment matching, filtered fragment summaries, known modification candidate search, P1 peak annotation, precursor-constrained MS2 theoretical ion annotation, and Excel output.

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

### MVP-5.1 MS2 annotation

MS2 annotation is an evidence table, not de novo sequencing. RNA_MassHunter reads MS2 spectra from mzML, first matches each precursor against theoretical digestion fragments, then matches observed MS2 peaks against unmodified `c`/`y` ions from the precursor-compatible parent fragment. Neutral losses, base losses, and strict modification-site localization are deferred to later MVPs.

```yaml
ms2_annotation:
  enabled: true
  mz_tolerance_ppm: 20
  min_peak_intensity: 10
  min_relative_intensity_percent: 1.0
  max_peaks_per_spectrum: 500
  constrain_by_precursor: true
  fallback_to_all_ions_if_no_precursor_match: false
  output_all_peak_annotations: false
```

`MS2_Ion_Matches` reports only matched or multiple-candidate peaks. Unmatched peaks are kept in `MS2_Unmatched_Peaks` and capped by `ms2_annotation.max_unmatched_peaks`. One-nucleotide ions are marked as low-information evidence with `Informative_Ion = false`; parent-fragment evidence is summarized in `MS2_Fragment_Evidence`. When an mzML file has no MS2 spectra, the run still succeeds and the Excel report records `Total_MS2_Spectra = 0` with annotation skipped/no MS2 status.

MVP-5.2 can also match an MS2 precursor to a theoretical digestion fragment plus one compatible known-modification mass shift. Enable this with `ms2_annotation.include_modified_precursor_candidates`. Base compatibility is required and mass-neutral/isobaric modifications are excluded by default; at most one modification is considered per fragment. These candidates indicate that the precursor mass may be explained by a modification, but do not localize that modification. Fragment-ion matching remains an unmodified c/y-ion approximation in MVP-5.2, and this limitation is recorded in the candidate and evidence sheets.

MVP-5.3 extends each modified precursor candidate into position-specific modified c/y ions. Candidate positions are the matching target bases in the parent fragment; ions containing a candidate position receive the known modification mass shift, while unmodified counterparts are retained for comparison. `MS2_Modification_Localization_Evidence` summarizes support for each possible position. This is localization evidence, not definitive site assignment: one-nucleotide ions are low-information, strong evidence requires at least three informative matches with both c- and y-ion support, and similarly supported positions are explicitly marked `ambiguous-multiple-positions`.

MVP-5.4 adds evidence ranking for review prioritization, not definitive modification assignment. `Modification_Evidence_Ranking` integrates MS1 fragment matches, known-modification candidates, modified precursor evidence, modified-ion matches, localization evidence, and available organism/tRNA context into `Final_Score` and `Final_Confidence`. Isobaric candidates cannot become highly confident from precursor evidence alone; ambiguous localization and evidence dominated by one-nucleotide ions receive explicit penalties and warnings. `Modification_Evidence_Summary` provides run-level ranking counts.

MVP-5.4.1 calibrates confidence conservatively. High confidence requires Moderate/Strong localization or at least two informative modified-ion matches with both c- and y-series support. A precursor plus a single informative modified ion with Weak localization is retained as a Medium review candidate. `Confidence_Limiting_Factor` explains restrictions such as weak localization, single-ion support, one-sided ion series, or precursor-only evidence. Very High additionally requires Strong localization, at least three informative ions, both ion series, and good precursor/fragment-ion mass errors.

MVP-5.4.2 groups candidate positions that share a spectrum, parent fragment, and modification. `Modification_Ambiguity_Groups` distinguishes resolved, partially resolved, ambiguous, and unsupported groups without deleting individual candidates. Modified-ion matches are marked as position-discriminating only when their ion range supports one candidate position and excludes the alternatives; ions covering multiple candidate positions do not localize the modification. Ranking confidence therefore requires position-discriminating evidence, and unresolved groups receive explicit score penalties and limiting factors. Evidence ranking remains candidate prioritization rather than modification or site confirmation.

MVP-5.5 treats the user-curated PDF modification workbook as the authoritative import source for PDF-confirmed symbols, monoisotopic nucleoside masses, and mass shifts. Run `tools/import_curated_modifications.py` with a local workbook path to generate `data/modifications.yaml` and an Excel/TSV diff report; the source workbook is not committed unless redistribution rights are confirmed. Generated records include detectability, candidate policy, source priority, curation status, chemical group, and near-isobaric metadata. Pseudouridine is excluded from blind mass search but remains available to position-rule, literature, or user-specified workflows. Trimethylation and acetylation remain chemically distinct while sharing a near-isobaric review flag. Ranking uses curated/source metadata only as modest supporting evidence and never forces a modification call.

MVP-5.6 adds generic, user-configured biological context prioritization. Optional `priority_modifications`, `priority_keywords`, and `focus_positions` can move relevant candidates higher for review; all are empty in the sample config, and focus positions are never fixed to a built-in standard position. Organism, pathway, tRNA type, and anticodon boosts are applied only when explicit loaded rules support them. `Biological_Context_Priorities` records the settings used, while `Context_Supported_Candidates` lists positively boosted candidates. For example, a user may configure a modification family or research keyword such as `cnm5U` or `thioamide`, but no specific theme is hard-coded. Biological context prioritizes review and cannot establish modification identity or raise context-only candidates above the configured confidence cap without mass/MS evidence.

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

Use an alternate config without replacing the repository `config.yaml`:

```bash
python main.py --config path/to/alternate.yaml
```

Relative config paths are resolved from the current working directory; absolute paths are also accepted. Reports are written to `output/`, logs to `logs/`, and cache/checkpoint files to `.cache/`.

An optional SCIEX text profile can be routed to the shadow audit with:

```yaml
sciex_profile:
  enabled: false
  path: null
  intact_peak_detection:
    enabled: true
  intact_mass_comparison:
    enabled: true
    strict_tolerance_da: 1.0
    broad_tolerance_da: 5.0
  delta_mass_cluster_audit:
    enabled: true
    cluster_tolerance_da: 0.5
    duplicate_apex_tolerance_da: 0.25
    isotope_spacing_da: 1.003355
    isotope_spacing_tolerance_da: 0.15
    integer_spacing_tolerance_da: 0.15
    minimum_cluster_size: 2
    max_pair_spacing_da: 200.0
    max_pair_rows: 20000
```

`Mass`/`Intensity` profiles are eligible for intact neutral-mass peak detection. `Mass/Charge`/`Intensity` profiles remain m/z data and produce parser diagnostics only; no charge-1 conversion is assumed. When detection completes with at least one peak, the optional intact-mass comparison records proximity to the unmodified theoretical RNA mass and, when available, the nearest existing reconstructed intact mass. It does not perform modification lookup, chemical assignment, or molecular identity assignment. SCIEX parser, peak, and comparison sheets are hidden at the `standard` audit level and included at `audit`/`full`. The summary sheet uses the Excel-safe alias `SCIEX_Intact_Mass_Comp_Summary` because the descriptive name exceeds Excel's 31-character sheet-name limit. These shadow results do not affect formal scores, ranking, or candidate filtering.

When SCIEX routing is enabled, a lightweight input-identity shadow audit also compares conservative tRNA identity tokens from the source filename with the configured sequence name and anticodon. NFKC normalization, case folding, separator handling, known amino-acid three-letter codes, and explicit three-base RNA anticodons are used; generic filename words and replicate/run labels are ignored. A clear conflict adds one warning and marks the mass-comparison summary as biologically uninterpretable without stopping or changing its numerical calculation. `SCIEX_Input_Identity_Audit` is available at `audit` and `full`, hidden at `standard`, and never affects formal scoring, ranking, filtering, or molecular identity assignment.

The optional delta-mass cluster shadow audit groups comparison rows with a deterministic complete-link-style span bound, diagnoses duplicate-like apex proximity, and reports integer, isotope-like, and recurrent numerical spacing candidates. Pair analysis is limited by mass spacing and a deterministic row cap. The labels are numerical diagnostics only: no modification, isotope number, adduct, formula, charge, or molecular identity is assigned. `SCIEX_Delta_Mass_Clusters`, `SCIEX_Delta_Mass_Clust_Summary`, and `SCIEX_Delta_Mass_Relations` are available at `audit` and `full` and hidden at `standard`; identity conflicts disable biological interpretation but do not disable clustering or alter comparison values.

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
- Annotates MS2 spectra against precursor-compatible theoretical unmodified c/y fragment ions when `ms2_annotation.enabled` is true.
- Writes Excel output with applicable sheets such as `Run_summary`, `Input_parameters`, `mzML_diagnostics`, `Intact_mass_reconstruction`, `Charge_state_peaks`, `Theoretical_fragments`, `Fragment_MS1_matches`, `Fragment_MS1_filtered`, `Fragment_MS1_summary`, `Known_Modification_Candidates`, `Known_Modification_Summary`, `P1_*`, `MS2_Summary`, `MS2_Spectra`, `MS2_Parent_Candidates`, `MS2_Theoretical_Ions`, `MS2_Ion_Matches`, `MS2_Unmatched_Peaks`, `MS2_Fragment_Evidence`, and `Warnings`.

## Planned Later Features

- modified-fragment combination ranking beyond single known shifts
- unknown modification candidate generation
- mass balance check
- manual comparison


## MVP-5.9 Intact Reconstruction QC

MVP-5.9 adds quality-control diagnostics for full-length intact mass reconstruction. Existing intact reconstruction results are preserved, and additional QC columns plus `Intact_Reconstruction_QC` and `Intact_Reconstruction_Diag` sheets report charge-state support, charge continuity, neutral-mass spread, mass-error limits, intensity support, competing envelopes, limiting factors, reconstruction status, and a `Comparison_Ready` flag.

`Comparison_Ready` means the intact mass candidate has sufficient reconstruction quality for cross-condition full-length mass comparison. It does not confirm modification identity, modification position, or biological causality.


## MVP-5.9.1 Intact Envelope Quality Separation

MVP-5.9.1 separates charge-envelope reconstruction quality from proximity to the unmodified theoretical RNA mass. A modified full-length RNA can be far from the unmodified theoretical mass, so `Reliable` is driven by internal envelope consistency, charge-state support, RT consistency, and signal abundance rather than by `Unmodified_Theory_Delta` alone.

`Unmodified_Theory_Delta_Da` and `Unmodified_Theory_Delta_ppm` are annotation fields that can include total modification mass shifts. Optional reference masses can be configured to annotate agreement with external deconvolution results, but reference agreement does not confirm modification identity. `Comparison_Ready` indicates reconstruction quality suitable for cross-condition intact-mass comparison, not modification assignment.

## MVP-5.9.2 Neutral Mass Search Range

MVP-5.9.2 adds an SCiex-style absolute neutral mass search range for intact reconstruction. By default, full-length reconstruction review uses `intact_reconstruction.neutral_mass_range` with `enabled: true`, `min_da: 20000`, and `max_da: 30000`; candidates outside that range are flagged with `In_Neutral_Mass_Search_Range = false`, are not `Comparison_Ready`, and are excluded from the in-range dominant envelope used for full-length RNA review.

## MVP-5.9.3 Dominant Intact Envelope Selection

MVP-5.9.3 separates raw in-range intensity dominance from the Dominant Intact Envelope used for review. The neutral mass search range remains the absolute reconstruction range and can be changed for the analysis target; the default `20000-30000 Da` is only an initial setting. The strongest candidate inside that range is not automatically treated as the best intact envelope; `Dominant_Intact_Eligible` is selected from QC-eligible strict or review candidates using charge support, continuity, RT consistency, mass spread, internal error, and intensity. Optional `target_review_mass_range` settings can prioritize a specific RNA mass window without filtering reconstruction candidates. SCiex or other reference masses are optional external validation annotations: reconstruction, QC eligibility, `Comparison_Ready`, and dominant intact envelope selection work without reference data, and reference agreement is not required for successful reconstruction or modification identity.

MVP-5.9.5 adds `analysis.mode` with `full` as the default and `intact_only` for runs that execute mzML diagnostics, MS1 peak extraction, intact mass reconstruction, QC, grouping, comparison candidates, and common Excel output while skipping digestion, fragment mapping, modification search, MS2 annotation, evidence ranking, biological context, and review dashboard steps. The report includes `Workflow_Summary` to show executed and skipped steps. It also adds `Reconstructed_Mass_Spectrum`, which exports neutral-mass spectrum points with `Reconstructed_Envelope_Intensity` using configurable `intact_reconstruction.mass_spectrum_output` settings; the default intensity is `Total_Supporting_Intensity` normalized to 100 percent within the exported spectrum.

MVP-5.9.6 adds the optional `intact_reconstruction.engine: rt_localized` reconstruction engine while keeping `legacy_cluster` available. The RT-localized engine builds charge envelopes only from peaks in the same configurable RT window, uses anchor peaks to directly search predicted neighboring-charge m/z values, diagnoses charge gaps and missing charges, and reports local intensity and internal envelope mass consistency. RT windows, charge requirements, m/z tolerances, and neutral-mass estimator are configurable. Candidate generation does not use unmodified theoretical mass or reference masses; outputs remain deconvolution candidates, not absolute molecular identifications.

MVP-5.9.7 adds RT-localized engine precision diagnostics: per-candidate strict/review pass-fail matrices, generic quality tiers (`Tier_1_high_quality` through `Tier_4_rejected`), charge extension checks, split-envelope merge diagnostics, peak-sharing metrics, and criteria-based legacy-vs-RT engine matching. `Comparison_Ready` is now tier-driven by configurable `comparison_ready_tiers`; reference masses and target review ranges remain annotation/review aids and do not improve QC tier, strict/review eligibility, or dominant intact envelope selection. The Excel report includes `RT_Engine_QC_Summary` for tier counts, failure reasons, missing/extension charge status, split envelopes, peak usage, and engine comparison status.

MVP-5.9.8a adds diagnostic-only competition grouping for intact candidates that share supporting local peaks. Candidates are connected into competition groups by shared local-peak evidence and RT proximity, then scored with envelope-internal evidence only, including charge support, RT consistency, local intensity, extension/split support, internal error, mass spread, peak sharing, peak usage, charge gaps, and severe limiting factors. This stage does not exclude candidates, change `Comparison_Ready`, alter quality tiers, or change reconstructed mass spectrum selection; the score distribution is intended to guide a later peak-assignment design. Reference masses and target review ranges remain annotation-only and do not affect competition grouping, evidence score, or rank.

## MVP-5.9.4 Intact Envelope Grouping

MVP-5.9.4 records supporting peak identities for intact reconstruction candidates, merges exact duplicate peak-set candidates, and groups overlapping envelopes using configurable mass, RT, shared-peak, and shared-charge thresholds. Condition comparison uses one QC-selected group representative per envelope group, exported in `Intact_Comparison_Candidates`, while `Target_Review_Candidates` is an optional review aid when a target review range is configured. Reference masses remain annotation only and do not affect grouping, QC, or representative selection.


## Data Notes

`data/modifications.yaml` is treated as the user-confirmed PDF-derived modification dictionary when present. MVP-4 loads and validates it, then reports known single-modification candidates when observed fragment masses can be explained by dictionary mass shifts. Unknown modification search is still deferred.

MS1 mass differences alone cannot distinguish isobaric modifications such as pseudouridine. By default, modifications with `mass_shift_from_unmodified = 0` are excluded from `Known_Modification_Candidates`. Set `modification_search.include_isobaric_modifications: true` only when you explicitly want to report those mass-neutral candidates. RNA_MassHunter does not use special priority logic for standard position 34; Excel reports use sequence start/end positions for fragment review.

`data/base_masses.yaml` contains Mongo Oligo-compatible placeholder masses. Confirm these values against the final Mongo Oligo settings before quantitative final analysis.

## MVP-5.8: Candidate comparison report

Multiple RNA_MassHunter Excel reports can be compared across experimental conditions using `tools/compare_reports.py`.

```bash
python tools/compare_reports.py --input Run1=output/RNA_MassHunter_MVP5_run1.xlsx --input Run2=output/RNA_MassHunter_MVP5_run2.xlsx --output output/RNA_MassHunter_Comparison.xlsx
```

Each input may be written as `Condition=path`. When a condition name is omitted, the input workbook filename is used.

The comparison workbook contains:

- `Comparison_Summary`
- `Candidate_Comparison`
- `Condition_Presence_Matrix`
- `Review_Priority_Changes`
- `Ambiguity_Comparison`
- `Candidate_Delta_Summary`

The report compares candidate presence, review priority, final confidence, final score, candidate positions, and ambiguity status between conditions. It is intended for review prioritization and does not establish modification identity or experimental causality.

## MVP-5.8 Candidate Comparison Report

MVP-5.8 adds

MVP-5.9.8c compares diagnostic dry-run competitive assignment across strict, balanced, sensitive, and permissive threshold scenarios. It reports candidates selected consistently across scenarios and supports optional audit masses for result review only. Audit masses do not affect grouping, evidence scoring, assignment, quality tiers, Comparison Ready, representatives, dominant selection, or reconstructed spectrum output. This sensitivity analysis is intended for threshold review before any production connection of assignment results.

MVP-5.9.8d adds assignment strict/review eligibility and can optionally apply it to formal Comparison Ready and representative selection with `competitive_assignment.apply_to_comparison_ready`. Pre-assignment values are retained, ambiguous and threshold-sensitive candidates remain available in `Assignment_Ambiguous_Candidates`, and reconstructed spectrum assignment filtering supports `none`, `strict`, `review`, `balanced_selected`, and `all` (`review` is recommended for assignment-aware review). The connection is disabled by default for backward compatibility. Reference, audit, target-review, and theoretical masses do not affect assignment eligibility or selection.
## Audit output levels

`--audit-level` separates routine reports from research/development shadow-audit output without changing formal matching, candidates, scores, confidence, ranks, identity, localization, or warnings. The default remains `full` for backward compatibility.

- `standard` is intended for routine analysis. It runs the formal workflow, omits detailed shadow sheets and Top shadow columns, and records `not_run` audit status in Diagnostics without executing the heavy MS1/MS2 shadow builders.
- `audit` is intended for quality review. It runs the shadow audits, writes their Summary sheets and `Audit_Status`, and omits group/detail sheets. Some existing audits build Summary and Detail together, so this level currently saves Excel size more reliably than computation time.
- `full` is intended for research and development. It runs and writes every existing Summary, group, and Detail audit and retains all Top and Diagnostics shadow columns.

Examples:

```bash
python main.py --audit-level standard
python main.py --audit-level audit
python main.py --audit-level full
python main.py --config config.sample.yaml --audit-level standard
```

Running `python main.py` is equivalent to `--audit-level full`. `standard` normally provides the largest runtime and memory saving because shadow calculations are skipped; `audit` and `full` may have similar compute cost where an audit currently creates its Summary and Detail in one deterministic builder. The selected level is reported in `Run_summary` and Diagnostics. Audit results remain disconnected from formal results at every level.

### MS/MS biological position prior (shadow)

`ms2_annotation.biological_position_prior` evaluates candidate positions using input-sequence 1-based numbering, parent-base compatibility, and structural/isobaric alternatives. Canonical landmarks must be configured explicitly (the configured `sequence.wobble_position` is also accepted as the `wobble` landmark); no Sprinzl numbering is assumed. The results are diagnostic-only by default (`apply_to_final_score: false`) and do not change `Final_Score`, `Final_Confidence`, rank, candidate inclusion, or review priority. Review `Modification_Position_Priors` and `MS2_Biological_Plausibility` in the Excel report.

### MS/MS identity physical peak sharing (shadow)

A candidate match ID includes theoretical-ion metadata and is not the same as a physical observed peak. Identity shadow evaluation audits physical peak sharing across candidates. Peaks shared by structural isomers are group-level evidence and are not counted repeatedly as individual structure evidence. Group-level localization status is a ceiling for candidate-level position resolution. These fields remain shadow-only and do not change the formal score, confidence, rank, review priority, matching, or localization results.

### MS/MS unmatched modified-ion audit (shadow)

`MS2_Unmatched_Ion_Audit` and `MS2_Unmatched_Ion_Summary` classify why modified theoretical ions lack an existing observed match. An unmatched ion does not establish that a modification is absent: scan range, formal tolerance, configured intensity thresholds, post-threshold filtering, and fragmentation coverage are reported separately. The audit is not a hard filter, does not alter formal matching, and is disconnected from `Final_Score`, `Final_Confidence`, rank, candidate inclusion, localization, and review priority.

### MS/MS ambiguous nearby-peak clusters (shadow)

`MS2_Ambiguous_Peak_Clusters`, `MS2_Ambiguous_Peak_Detail`, and `MS2_Ambiguity_Summary` audit ambiguous nearby peaks as deterministic physical-peak clusters. Physical peak sharing across candidates is reported separately from competition among theoretical ions. An ambiguous peak neither confirms nor excludes a modification; this shadow audit is not a hard filter, does not alter formal matching, and is disconnected from `Final_Score`, confidence, rank, candidate inclusion, localization, and review priority.

In `MS2_Ambiguous_Peak_Clusters`, `Best_Peak_*` means the diagnostic raw peak closest to the theoretical m/z within the audit window (absolute m/z error, then higher intensity, then lower m/z). It may be a zero-intensity raw peak and is not the formal annotation best match; formal matches are selected only from intensity-filtered annotation-input peaks.

### MS/MS effective ambiguity stages (shadow)

`MS2_Effective_Ambiguity`, `MS2_Effective_Ambig_Detail`, and `MS2_Effective_Ambig_Summary` separate raw-window multiplicity from ambiguity that remains among positive-intensity peaks, positive peaks inside formal tolerance, and assignments already present in formal match tables. The strongest applicable stage is reported as `formal_match > formal_tolerance > positive_intensity > raw_only > none`. This classification reuses existing results without rematching and is not applied to `Final_Score`, confidence, rank, matching, identity, localization, or review priority.


### Fragment MS1 truncation audit (shadow)

`MS1_Truncation_Audit`, `MS1_Truncation_Detail`, and `MS1_Truncation_Summary` capture every pre-truncation Fragment MS1 match and compare the formal `fragment_mapping.max_matches_per_fragment` selection with an unlimited shadow. The formal selection remains `abs(mass_error_ppm)` ascending, then intensity descending; stable ties retain charge-ascending and input-peak order. The audit also compares filter-first, unique-physical-peak, charge-balanced, and tier-first diagnostic selections.

The unlimited shadow reuses the existing fragment filter, known-modification search, and evidence-ranking functions, then reports potential candidate-key, score, confidence, rank, Top-50, and cnm5U 36/37/38 changes. It appends diagnostic columns to `Top_Modification_Candidates` and `MS2_Unmatched_Ion_Diagnostics`. All `Applied_To_Final_Score` fields are `false`: the audit does not replace formal matches, change the configured 20-match cap, or propagate into formal candidates, scores, confidence, rank, localization, or review priority. Missing configuration keys use code defaults: audit enabled, unlimited expanded shadow enabled, and final-score application disabled.

`MS1_Selection_Strategy`, `MS1_Selection_Detail`, and `MS1_Selection_Summary` add a second non-mutating A/B audit across `current`, `filter_first`, `tier_then_error`, and `unlimited`. Filter-first applies the exact `Fragment_MS1_filtered` conditions before the existing error/intensity top-20 ordering. Tier-then-error prioritizes filter pass, peak tier, confidence, absolute ppm error, intensity, charge, and a deterministic physical-peak/input tie-break. The default recommendation and readiness diagnosis are shadow metadata only; `Apply_MS1_Selection_Strategy_To_Formal_Result` defaults to `false`.

`MS1_Top50_Shadow`, `MS1_Peak_Dedup_Detail`, and `MS1_Top50_Dedup_Summary` perform a full downstream tier-then-error top-50 simulation and an exact physical-peak ID audit. Exact-ID deduplication keeps scan/peak-index identity separate from near-m/z grouping, compares fragment-level charge deduplication with a diagnostic global constraint, and rebuilds known-modification candidates, evidence ranking, and review-dashboard shadows. Code defaults enable both audits with a 50-match shadow limit while all formal-application flags remain `false`.

`MS1_CrossFrag_Ambiguity`, `MS1_CrossFrag_Detail`, and `MS1_CrossFrag_Summary` treat one physical MS1 peak assigned to multiple theoretical fragments as assignment ambiguity rather than automatic duplication. The shadow compares full-count, winner-take-all, equal-fraction, quality-weighted, ambiguity-flag-only, and overlapping-fragment-family weighting. Deterministic best-assignment ordering uses formal-filter pass, tier, confidence, absolute ppm error, intensity, fragment length, missed cleavage, and fragment ID/order. Code defaults enable `Enable_MS1_Cross_Fragment_Ambiguity_Audit`, keep `Apply_MS1_Cross_Fragment_Ambiguity_To_Formal_Result=false`, and append only shadow columns to Top candidates and Diagnostics; no formal candidate, score, confidence, rank, localization, or review result is changed.

## Phase 1 constraint-aware composite modifications (shadow)

The legacy modification search remains the formal model and represents known entries primarily as single mass shifts. That model cannot by itself prove that multiple shifts occupy compatible atoms, prevent parent/finished-derivative double counting, or represent a modification on an inter-nucleotide bond. Phase 1 therefore adds a separate, non-propagating composition audit. It does not change formal candidate generation, digestion, fragment masses, scoring, confidence, rank, localization, or review priority.

`data/nucleoside_slots.yaml` defines immutable initial slot/state values for A, G, C, and U. `data/modification_transforms_v2.yaml` defines nucleoside state transitions, occupied slots, prerequisites, forbidden states, included/superseded components, evidence, and signed elemental-composition deltas. The bounded deterministic composer tries at most three transformation definitions per candidate. A finished derived transformation counts as one search component and lists the parent components that it already includes; this prevents their masses from being added again. Valid states and rejected attempts with structured reasons are reported separately, and composition-identical distinct states are grouped as isomers. Exact mass deltas are calculated from C/H/N/O/P/S/Se composition, not copied from legacy mass shifts.

The side-chain thioamide and oxidation definitions are explicitly hypothetical Phase 1 states. A generated candidate means that the declared slot/state constraints are internally consistent; it does not establish that the structure exists in a sample. New transformations should normally be added by extending the two YAML schemas with `from_state`, `to_state`, `occupies_slots`, `requires`/`forbids`, composition, evidence, and legacy mapping. Python changes should only be needed for a genuinely new constraint concept.

`data/backbone_modifications.yaml` represents phosphorothioate on bond `left_position_right_position` as an O-to-S delta. Its stereochemistry is `unknown` when unspecified. The cleavage shadow evaluates RNase T1, RNase A, Nuclease P1, and RNase T2 without changing `find_cleavage_sites` or formal digestion. A phosphorothioate-blocked site joins adjacent fragments and contributes its bond mass once. `Blocked_Cleavage_Count` and `Cleavage_Origin=phosphorothioate_blocked` remain distinct from `Stochastic_Missed_Cleavage_Count`; mixed cases are labelled separately. T1/P1 blocking is recorded as the research hypothesis under test, while A/T2 remain `potentially_blocked` with unknown evidence rather than being asserted as universal rules.

Output by audit level is:

- `standard`: the composite builder is not run, no five composite sheets are written, and Diagnostics reports `not_run` rather than zero.
- `audit`: `Composite_Mod_Summary`, the summary-equivalent `Cleavage_Block_Audit`, and three `Audit_Status` entries are written; detail sheets are suppressed.
- `full`: `Composite_Mod_Candidates`, `Composite_Mod_Invalid`, `Composite_Mod_Summary`, `Backbone_Mod_Candidates`, and `Cleavage_Block_Audit` are written.

All new rows use `Applied_To_Formal_Result=false`; `Formal_Change_Ready` is also false in Phase 1. `data/examples/composite_structure_example.yaml` is an `example_only` U37/two-bond hypothesis fixture for tests and shadow exploration. It is not read from `config.yaml`, is not applied automatically to Mac data, and is not evidence that the example structure occurs in any biological sample.

## Generic P1 + AP dinucleotide chemical-state audit (shadow)

The generic v2 P1 + AP system separates sequence-aware candidate construction, measured-feature auditing, and evidence interpretation. Candidate construction uses every adjacent 5′→3′ bond in an arbitrary input RNA, the schema-defined nucleoside transformations and position hypotheses that match that RNA, full elemental compositions, explicitly enabled linkage states, configured charges, and configured polarity. Composition-identical assignments are grouped without treating the first assignment as confirmed; every structural assignment remains available in a separate table.

Low-m/z dinucleotide ranges are configured independently from intact reconstruction. Acquisition, MS1 extraction, dinucleotide search, and MS2 product-ion ranges are reported separately. Observable groups are searched whether or not a target list is supplied. Raw matching uses sorted m/z indexes, same-scan profile points are reconstructed as spectrum-level local peaks, all disconnected RT features are retained, and the existing v1.1 background and same-scan CHNOSP isotope algorithms are reused. Sulfur evaluation is driven by the sulfur count in the complete composition and cannot by itself establish a phosphorothioate linkage.

Mass accuracy, background, isotope compatibility, physical-feature competition, precursor-compatible MS2 provenance, feature classification, and group interpretation use the same generic logic for every group. MS2 product ions below the MS1 extraction range remain in provenance counts. No dinucleotide fragmentation score or diagnostic ion is generated until a validated model exists. Position and original-bond localization always remain disabled after P1 + AP preparation.

Optional `p1_sap_dinucleotide.targets` entries contain an arbitrary label, theoretical m/z, and tolerance. They filter already-computed groups and features for display only; zero, one, or many targets produce identical generic analysis results. No target-specific JSON keys or processing functions are used.

At `audit` level, `P1_SAP_Dinuc_Summary`, `P1_SAP_Dinuc_Groups`, and `P1_SAP_Dinuc_Targets` are written. At `full` level, assignment, spectrum-peak, feature, isotope, competition, and MS2 detail sheets are added. JSON is stored under the generic `dinucleotide_audit` object. All rows use `Applied_To_Formal_Result=false`, `Formal_Change_Ready=false`, and `Formal_Result_Changed=false`.
