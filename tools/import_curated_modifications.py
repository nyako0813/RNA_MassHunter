#!/usr/bin/env python3
"""Import the curated PDF modification workbook into RNA_MassHunter YAML."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from openpyxl.styles import Font, PatternFill


PRIMARY_SHEET = "PDF_modifications_v0_1"
SUPPORTING_SHEET = "Non_PDF_supporting_entries"
DECISIONS_SHEET = "Curation_decisions_v0_3"
BOOL_COLUMNS = (
    "candidate_policy_include_by_mass_search",
    "candidate_policy_include_if_position_rule_exists",
    "candidate_policy_include_if_literature_supported",
    "candidate_policy_include_if_user_specified",
)


def _clean(value: Any, default: Any = "") -> Any:
    if pd.isna(value):
        return default
    return value.item() if hasattr(value, "item") else value


def _bool(value: Any, default: bool = False) -> bool:
    if pd.isna(value):
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _bases(value: Any) -> list[str]:
    text = str(_clean(value, "")).replace(";", ",")
    return [item.strip().upper().replace("T", "U") for item in text.split(",") if item.strip()]


def _detectability(value: Any) -> bool | str:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"limited", "conditional", "position_only"}:
        return "limited"
    return _bool(value)


def _source(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "source_file": str(_clean(record.get("source_file"), "")),
        "source_page": _clean(record.get("source_page"), ""),
        "source_field": str(_clean(record.get("source_field"), "")),
        "source_priority": str(_clean(record.get("source_priority"), "")),
        "curation_status": str(_clean(record.get("curation_status"), "")),
        "notes": str(_clean(record.get("notes"), "")),
    }


def curated_record(record: dict[str, Any]) -> dict[str, Any]:
    mod_id = str(_clean(record.get("id"), "")).strip()
    name = str(_clean(record.get("common_name"), mod_id)).strip()
    shift = float(_clean(record.get("mass_shift_from_unmodified"), 0.0))
    source = _source(record)
    chemical_group = str(_clean(record.get("chemical_group"), ""))
    isobaric_group = str(_clean(record.get("isobaric_group"), ""))
    near_group = str(_clean(record.get("near_isobaric_group"), ""))
    lower_identity = f"{mod_id} {name} {chemical_group}".lower()
    if "trimethyl" in lower_identity:
        isobaric_group = "trimethylation_group"
        near_group = near_group or "near_isobaric_42Da_group"
    elif "acetyl" in lower_identity:
        isobaric_group = "acetylation_group"
        near_group = near_group or "near_isobaric_42Da_group"
    policy = {
        "include_by_mass_search": _bool(record.get(BOOL_COLUMNS[0]), True),
        "include_if_position_rule_exists": _bool(record.get(BOOL_COLUMNS[1]), True),
        "include_if_literature_supported": _bool(record.get(BOOL_COLUMNS[2]), True),
        "include_if_user_specified": _bool(record.get(BOOL_COLUMNS[3]), True),
    }
    detectability = {
        "ms1": _detectability(record.get("detectability_ms1")),
        "ms2": _detectability(record.get("detectability_ms2")),
    }
    if mod_id.upper() in {"Y", "PSI", "PSEUDOURIDINE"} or name.lower() == "pseudouridine":
        shift = 0.0
        detectability["ms1"] = False
        if detectability["ms2"] is False:
            detectability["ms2"] = "limited"
        policy = {
            "include_by_mass_search": False,
            "include_if_position_rule_exists": True,
            "include_if_literature_supported": True,
            "include_if_user_specified": True,
        }
    result = {
        "id": mod_id,
        "symbol": str(_clean(record.get("symbol"), mod_id)),
        "name": name,
        "short_name": str(_clean(record.get("symbol"), mod_id)),
        "target_bases": _bases(record.get("target_bases") or record.get("base")),
        "base": str(_clean(record.get("base"), "")),
        "modified_nucleoside_mass_mono": float(_clean(record.get("modified_nucleoside_mass_mono"), 0.0)),
        "mass_shift_from_unmodified": shift,
        "mass_basis": str(_clean(record.get("mass_basis"), "")),
        "category": str(_clean(record.get("category"), "biological")),
        "chemical_group": chemical_group,
        "isobaric_group": isobaric_group,
        "near_isobaric_group": near_group,
        "detectability": detectability,
        "candidate_policy": policy,
        "source": source,
        "source_priority": source["source_priority"],
        "curation_status": source["curation_status"],
        "sources": [{
            "type": "user_pdf", "file": source["source_file"], "page": source["source_page"],
            "source_field": source["source_field"], "source_priority": source["source_priority"],
            "curation_status": source["curation_status"],
        }],
        "curation": {"status": source["curation_status"], "notes": source["notes"]},
    }
    return result


def supporting_record(record: dict[str, Any]) -> dict[str, Any]:
    bases = _bases(record.get("target_bases"))
    return {
        "id": str(_clean(record.get("id"), "")), "symbol": str(_clean(record.get("short_name"), "")),
        "name": str(_clean(record.get("name"), "")), "short_name": str(_clean(record.get("short_name"), "")),
        "target_bases": bases, "base": "", "modified_nucleoside_mass_mono": None,
        "mass_shift_from_unmodified": float(_clean(record.get("mass_shift"), 0.0)), "mass_basis": "mass_delta",
        "category": str(_clean(record.get("category"), "supporting")), "chemical_group": "",
        "isobaric_group": "", "near_isobaric_group": "", "detectability": {"ms1": True, "ms2": "limited"},
        "candidate_policy": {"include_by_mass_search": False, "include_if_position_rule_exists": False,
                             "include_if_literature_supported": False, "include_if_user_specified": True},
        "source": {"source_file": "", "source_page": "", "source_field": SUPPORTING_SHEET,
                   "source_priority": str(_clean(record.get("source_priority"), "")), "curation_status": "supporting",
                   "notes": str(_clean(record.get("notes"), ""))},
        "source_priority": str(_clean(record.get("source_priority"), "")), "curation_status": "supporting",
    }


def _existing_records(path: Path) -> dict[str, dict[str, Any]]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    records = data.get("modifications", data if isinstance(data, list) else [])
    return {str(row.get("id") or row.get("symbol") or ""): row for row in records if isinstance(row, dict)}


def _json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def diff_rows(old: dict[str, dict[str, Any]], new: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    new_lookup = {row["id"]: row for row in new}
    for mod_id in sorted(set(old) | set(new_lookup)):
        before, after = old.get(mod_id), new_lookup.get(mod_id)
        if before is None:
            status, notes = "added", "Present only in curated workbook."
        elif after is None:
            status, notes = "missing_in_excel", "Existing YAML entry is not present in the primary curated sheet."
        else:
            mass_changed = abs(float(before.get("mass_shift_from_unmodified", 0.0)) - float(after.get("mass_shift_from_unmodified", 0.0))) > 1e-6
            bases_changed = sorted(before.get("target_bases", [])) != sorted(after.get("target_bases", []))
            policy_changed = before.get("candidate_policy", {}) != after.get("candidate_policy", {})
            both_checked = (before.get("curation", {}) or {}).get("status") == "manually_checked" and after.get("curation_status") == "manually_checked"
            status = "conflict" if mass_changed and both_checked else "updated" if mass_changed or bases_changed or policy_changed else "unchanged"
            notes = "Authoritative curated workbook differs from manually checked YAML mass." if status == "conflict" else "Metadata differs." if status == "updated" else ""
        rows.append({
            "id": mod_id, "status": status,
            "old_mass_shift": before.get("mass_shift_from_unmodified", "") if before else "",
            "new_mass_shift": after.get("mass_shift_from_unmodified", "") if after else "",
            "old_target_bases": ",".join(before.get("target_bases", [])) if before else "",
            "new_target_bases": ",".join(after.get("target_bases", [])) if after else "",
            "old_candidate_policy": _json(before.get("candidate_policy", {})) if before else "",
            "new_candidate_policy": _json(after.get("candidate_policy", {})) if after else "",
            "notes": notes,
        })
    return rows


def write_report(path: Path, rows: list[dict[str, Any]], decisions: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".tsv":
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        return
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        pd.DataFrame(rows).to_excel(writer, sheet_name="Import_Diff", index=False)
        decisions.to_excel(writer, sheet_name="Curation_Decisions", index=False)
        for sheet in writer.book.worksheets:
            sheet.freeze_panes = "A2"
            sheet.auto_filter.ref = sheet.dimensions
            for cell in sheet[1]:
                cell.font = Font(name=cell.font.name or "Calibri", size=cell.font.sz or 11, bold=True, color="FFFFFF")
                cell.fill = PatternFill(fill_type="solid", fgColor="0F766E")
            for column in sheet.columns:
                width = min(60, max(12, max(len(str(cell.value or "")) for cell in column) + 2))
                sheet.column_dimensions[column[0].column_letter].width = width


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Curated .xlsx workbook")
    parser.add_argument("--output", type=Path, required=True, help="Generated modifications.yaml")
    parser.add_argument("--report", type=Path, required=True, help="Diff report (.xlsx or .tsv)")
    parser.add_argument("--compare", type=Path, help="Existing YAML to compare; defaults to --output before overwrite")
    parser.add_argument("--include-supporting-entries", action="store_true", help="Include terminal/adduct/artifact/constants from the supporting sheet")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.input.exists():
        raise FileNotFoundError(args.input)
    primary = pd.read_excel(args.input, sheet_name=PRIMARY_SHEET)
    required = {"id", "symbol", "common_name", "target_bases", "mass_shift_from_unmodified"}
    missing = sorted(required - set(primary.columns))
    if missing:
        raise ValueError(f"{PRIMARY_SHEET} is missing required columns: {missing}")
    records = [curated_record(row) for row in primary.to_dict(orient="records") if str(_clean(row.get("id"), "")).strip()]
    if args.include_supporting_entries:
        supporting = pd.read_excel(args.input, sheet_name=SUPPORTING_SHEET)
        records.extend(supporting_record(row) for row in supporting.to_dict(orient="records") if str(_clean(row.get("id"), "")).strip())
    compare_path = args.compare or args.output
    old = _existing_records(compare_path)
    decisions = pd.read_excel(args.input, sheet_name=DECISIONS_SHEET)
    payload = {
        "schema_version": "RNA_MassHunter_modifications_v0.3_curated",
        "description": "Generated from the user-curated PDF modification workbook. PDF-confirmed mass values are authoritative.",
        "source_workbook": args.input.name,
        "modifications": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, allow_unicode=True, sort_keys=False, width=120)
    report = diff_rows(old, records)
    write_report(args.report, report, decisions)
    counts = pd.Series([row["status"] for row in report]).value_counts().to_dict()
    print(f"Imported {len(records)} modifications -> {args.output}")
    print(f"Diff report -> {args.report}; status counts: {counts}")


if __name__ == "__main__":
    main()
