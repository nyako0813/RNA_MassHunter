# Curated modification import

The curated workbook is a local source input and is not committed to this repository unless its redistribution rights are explicitly confirmed. Place the workbook anywhere accessible locally, then run:

```bash
python tools/import_curated_modifications.py \
  --input /path/to/RNA_MassHunter_v2_PDF_modifications_v0_3_curated.xlsx \
  --output data/modifications.yaml \
  --compare data/modifications.yaml \
  --report output/curated_modification_import_report.xlsx
```

`PDF_modifications_v0_1` is the authoritative biological-modification input. PDF-confirmed symbol, monoisotopic nucleoside mass, and unmodified mass shift take precedence over older code values. `Non_PDF_supporting_entries` is excluded by default because it contains constants, adducts, terminal forms, and artifacts; add `--include-supporting-entries` only when those records are intentionally required in the generated YAML.

The report compares the generated records with the existing YAML before overwrite. Review every `conflict` and `missing_in_excel` row before committing a generated dictionary.
