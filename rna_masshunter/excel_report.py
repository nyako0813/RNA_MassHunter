from dataclasses import asdict, is_dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
from openpyxl.utils import get_column_letter


EXCEL_MAX_ROWS = 1_048_576
DATA_START_ROW = 3
EXCEL_DATA_ROW_LIMIT = EXCEL_MAX_ROWS - DATA_START_ROW


INTACT_COLUMNS = [
    "Observed_Mass",
    "Charge_State_Count",
    "Charge_States",
    "Supporting_Peak_Count",
    "Total_Intensity",
    "Theoretical_Mass",
    "Mass_Error_Da",
    "Mass_Error_ppm",
    "Assignment",
    "Confidence",
    "Warnings",
]

CHARGE_COLUMNS = ["Cluster_ID", "mz", "Intensity", "RT", "Scan_ID", "Charge", "Neutral_Mass", "Peak_Tier"]

THEORETICAL_FRAGMENT_COLUMNS = [
    "Fragment_ID",
    "Target_ID",
    "Sequence",
    "Length",
    "Start",
    "End",
    "Standard_Start",
    "Standard_End",
    "Enzyme",
    "Missed_Cleavages",
    "Terminal_Form",
    "Unmodified_Mass",
    "Warnings",
]

FRAGMENT_MS1_MATCH_COLUMNS = [
    "Match_ID",
    "Fragment_ID",
    "Target_ID",
    "Sequence",
    "Start",
    "End",
    "Standard_Start",
    "Standard_End",
    "Enzyme",
    "Missed_Cleavages",
    "Terminal_Form",
    "Fragment_Mass",
    "Charge",
    "Theoretical_mz",
    "Observed_mz",
    "Mass_Error_Da",
    "Mass_Error_ppm",
    "Intensity",
    "RT",
    "Scan_ID",
    "Peak_Tier",
    "Confidence",
    "Warnings",
]

SHEET_DESCRIPTIONS = {
    "Run_summary": "Run-level summary for this RNA_MassHunter MVP-3 report.",
    "Input_parameters": "Flattened parameters loaded from config.yaml.",
    "mzML_diagnostics": "mzML scan counts, ranges, precursor metadata, and warnings.",
    "Intact_mass_reconstruction": "Reconstructed intact mass clusters and mass errors.",
    "Charge_state_peaks": "Peak and charge-state evidence supporting reconstructed masses.",
    "Theoretical_fragments": "Theoretical RNase digestion fragments and terminal forms.",
    "Fragment_MS1_matches": "MS1 peak matches for unmodified theoretical fragments.",
    "Warnings": "Warnings and errors recorded during startup, loading, and analysis.",
}


def _flatten_dict(data: dict[str, Any], prefix: str = "") -> list[dict[str, Any]]:
    rows = []
    for key, value in data.items():
        full_key = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            rows.extend(_flatten_dict(value, full_key))
        else:
            rows.append({"Parameter": full_key, "Value": value})
    return rows


def _autosize_and_freeze(writer: pd.ExcelWriter) -> None:
    for worksheet in writer.book.worksheets:
        worksheet.freeze_panes = "A2" if worksheet.title == "Index" else "A4"
        for column_cells in worksheet.columns:
            max_length = 0
            column = get_column_letter(column_cells[0].column)
            for cell in column_cells:
                value = "" if cell.value is None else str(cell.value)
                max_length = max(max_length, min(len(value), 60))
            worksheet.column_dimensions[column].width = max(10, max_length + 2)


def _sheet_link(sheet_name: str, cell: str = "A1") -> str:
    return f"#'{sheet_name}'!{cell}"


def _coerce_to_frame(value: Any) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, list):
        return pd.DataFrame(value)
    if isinstance(value, dict):
        return pd.DataFrame([value])
    return pd.DataFrame([{"Value": value}])


def _add_index_and_backlinks(writer: pd.ExcelWriter, sheet_names: list[str]) -> None:
    workbook = writer.book
    index_sheet = workbook["Index"]
    for row_index, sheet_name in enumerate(sheet_names, start=2):
        link_cell = index_sheet.cell(row=row_index, column=1)
        link_cell.value = sheet_name
        link_cell.hyperlink = _sheet_link(sheet_name, "A1")
        link_cell.style = "Hyperlink"

    for sheet_name in sheet_names:
        worksheet = workbook[sheet_name]
        worksheet["A1"] = "← Back to Index"
        worksheet["A1"].hyperlink = _sheet_link("Index", "A1")
        worksheet["A1"].style = "Hyperlink"


def _fragment_rows(theoretical_fragments: list[Any]) -> list[dict[str, Any]]:
    rows = []
    for item in theoretical_fragments:
        raw = asdict(item) if is_dataclass(item) else dict(item)
        fragment_warnings = raw.get("warnings", [])
        if isinstance(fragment_warnings, list):
            fragment_warnings = "; ".join(map(str, fragment_warnings))
        rows.append(
            {
                "Fragment_ID": raw.get("fragment_id"),
                "Target_ID": raw.get("target_id"),
                "Sequence": raw.get("sequence"),
                "Length": len(raw.get("sequence") or ""),
                "Start": raw.get("start"),
                "End": raw.get("end"),
                "Standard_Start": raw.get("standard_start"),
                "Standard_End": raw.get("standard_end"),
                "Enzyme": raw.get("enzyme"),
                "Missed_Cleavages": raw.get("missed_cleavages"),
                "Terminal_Form": raw.get("terminal_form"),
                "Unmodified_Mass": raw.get("unmodified_mass"),
                "Warnings": fragment_warnings,
            }
        )
    return rows


def _fragment_ms1_match_rows(fragment_ms1_matches: list[Any]) -> list[dict[str, Any]]:
    rows = []
    for item in fragment_ms1_matches:
        raw = asdict(item) if is_dataclass(item) else dict(item)
        match_warnings = raw.get("warnings", [])
        if isinstance(match_warnings, list):
            match_warnings = "; ".join(map(str, match_warnings))
        rows.append(
            {
                "Match_ID": raw.get("match_id"),
                "Fragment_ID": raw.get("fragment_id"),
                "Target_ID": raw.get("target_id"),
                "Sequence": raw.get("sequence"),
                "Start": raw.get("start"),
                "End": raw.get("end"),
                "Standard_Start": raw.get("standard_start"),
                "Standard_End": raw.get("standard_end"),
                "Enzyme": raw.get("enzyme"),
                "Missed_Cleavages": raw.get("missed_cleavages"),
                "Terminal_Form": raw.get("terminal_form"),
                "Fragment_Mass": raw.get("fragment_mass"),
                "Charge": raw.get("charge"),
                "Theoretical_mz": raw.get("theoretical_mz"),
                "Observed_mz": raw.get("observed_mz"),
                "Mass_Error_Da": raw.get("mass_error_da"),
                "Mass_Error_ppm": raw.get("mass_error_ppm"),
                "Intensity": raw.get("intensity"),
                "RT": raw.get("rt"),
                "Scan_ID": raw.get("scan_id"),
                "Peak_Tier": raw.get("peak_tier"),
                "Confidence": raw.get("confidence"),
                "Warnings": match_warnings,
            }
        )
    return rows


def _as_bool(value: Any, default: bool = True) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _as_positive_int(value: Any, default: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed > 0 else default


def _append_excel_warning(
    warnings: list[dict[str, Any]],
    sheet_name: str,
    original_rows: int,
    written_rows: int,
) -> None:
    warnings.append(
        {
            "Timestamp": datetime.now().isoformat(timespec="seconds"),
            "Level": "WARNING",
            "Source": "excel_report",
            "Message": "Excel sheet was truncated because it exceeded max_excel_rows_per_sheet.",
            "Context": {"sheet": sheet_name, "original_rows": original_rows, "written_rows": written_rows},
        }
    )


def _truncate_frame_if_needed(
    sheet_name: str,
    frame: pd.DataFrame,
    max_rows: int,
    truncate_large_sheets: bool,
    warnings: list[dict[str, Any]],
    truncations: list[dict[str, Any]],
) -> pd.DataFrame:
    original_rows = len(frame)
    safe_limit = min(max_rows, EXCEL_DATA_ROW_LIMIT)
    if original_rows <= safe_limit:
        return frame

    if truncate_large_sheets:
        written_rows = safe_limit
    else:
        written_rows = EXCEL_DATA_ROW_LIMIT
    written_rows = min(written_rows, original_rows)
    _append_excel_warning(warnings, sheet_name, original_rows, written_rows)
    truncations.append({"sheet": sheet_name, "original_rows": original_rows, "written_rows": written_rows})
    return frame.head(written_rows).copy()


def _truncation_summary(truncations: list[dict[str, Any]]) -> str:
    if not truncations:
        return "None"
    return "; ".join(
        f"{item['sheet']}: {item['original_rows']} -> {item['written_rows']}" for item in truncations
    )


def write_excel_report(
    output_dir: str | Path,
    config,
    diagnostics: dict[str, Any],
    intact_results: list[Any],
    charge_state_peaks: list[dict[str, Any]],
    warnings: list[dict[str, Any]],
    modifications: list[Any] | None = None,
    rule_set: dict[str, Any] | None = None,
    pathways: list[dict[str, Any]] | None = None,
    theoretical_fragments: list[Any] | None = None,
    fragment_ms1_matches: list[Any] | None = None,
    optional_results: dict[str, Any] | None = None,
) -> Path:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / f"RNA_MassHunter_MVP3_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"

    reporting = config.reporting or {}
    max_excel_rows = _as_positive_int(reporting.get("max_excel_rows_per_sheet"), 100000)
    truncate_large_sheets = _as_bool(reporting.get("truncate_large_sheets"), True)
    max_charge_state_peak_rows = _as_positive_int(
        reporting.get("max_charge_state_peak_rows", config.reconstruction.get("max_charge_state_peak_rows")),
        max_excel_rows,
    )
    truncations: list[dict[str, Any]] = []

    intact_rows = []
    for item in intact_results:
        raw = asdict(item) if is_dataclass(item) else dict(item)
        intact_rows.append(
            {
                "Observed_Mass": raw.get("observed_mass"),
                "Charge_State_Count": raw.get("charge_state_count"),
                "Charge_States": ",".join(map(str, raw.get("charge_states", []))),
                "Supporting_Peak_Count": raw.get("supporting_peak_count"),
                "Total_Intensity": raw.get("total_intensity"),
                "Theoretical_Mass": raw.get("theoretical_mass"),
                "Mass_Error_Da": raw.get("mass_error_da"),
                "Mass_Error_ppm": raw.get("mass_error_ppm"),
                "Assignment": raw.get("assignment"),
                "Confidence": raw.get("confidence"),
                "Warnings": raw.get("warnings"),
            }
        )

    charge_state_peak_rows = charge_state_peaks
    if len(charge_state_peaks) > max_charge_state_peak_rows and truncate_large_sheets:
        _append_excel_warning(warnings, "Charge_state_peaks", len(charge_state_peaks), max_charge_state_peak_rows)
        truncations.append(
            {
                "sheet": "Charge_state_peaks",
                "original_rows": len(charge_state_peaks),
                "written_rows": max_charge_state_peak_rows,
            }
        )
        charge_state_peak_rows = charge_state_peaks[:max_charge_state_peak_rows]

    theoretical_fragments = theoretical_fragments or []
    fragment_ms1_matches = fragment_ms1_matches or []

    input_parameters = {
        "project": config.project,
        "input": config.input,
        "organism": config.organism,
        "sequence": config.sequence,
        "experiment": config.experiment,
        "instrument": config.instrument,
        "reconstruction": config.reconstruction,
        "digestion": config.digestion,
        "alkaline_phosphatase": config.alkaline_phosphatase,
        "fragment_mapping": config.fragment_mapping,
        "peak_filtering": config.peak_filtering,
        "performance": config.performance,
        "reporting": config.reporting,
    }

    data_sheets: dict[str, pd.DataFrame] = {
        "Input_parameters": pd.DataFrame(_flatten_dict(input_parameters)),
        "mzML_diagnostics": pd.DataFrame([diagnostics] if diagnostics else [{}]),
        "Intact_mass_reconstruction": pd.DataFrame(intact_rows, columns=INTACT_COLUMNS),
        "Charge_state_peaks": pd.DataFrame(charge_state_peak_rows, columns=CHARGE_COLUMNS),
        "Theoretical_fragments": pd.DataFrame(_fragment_rows(theoretical_fragments), columns=THEORETICAL_FRAGMENT_COLUMNS),
        "Fragment_MS1_matches": pd.DataFrame(_fragment_ms1_match_rows(fragment_ms1_matches), columns=FRAGMENT_MS1_MATCH_COLUMNS),
    }
    for sheet_name, value in (optional_results or {}).items():
        if sheet_name in {"Index", "Run_summary", "Warnings"}:
            continue
        data_sheets[sheet_name[:31]] = _coerce_to_frame(value)

    truncated_data_sheets = {
        sheet_name: _truncate_frame_if_needed(
            sheet_name,
            frame,
            max_excel_rows,
            truncate_large_sheets,
            warnings,
            truncations,
        )
        for sheet_name, frame in data_sheets.items()
    }

    summary_rows = [
        {"Item": "Project", "Value": config.project.get("name")},
        {"Item": "Generated", "Value": datetime.now().isoformat(timespec="seconds")},
        {"Item": "Modification dictionary entries", "Value": len(modifications or [])},
        {"Item": "Rule set", "Value": config.organism.get("rule_set") or (rule_set or {}).get("id") or (rule_set or {}).get("name")},
        {"Item": "Pathway files", "Value": len(pathways or [])},
        {"Item": "Intact mass candidates", "Value": len(intact_results)},
        {"Item": "Theoretical fragments", "Value": len(theoretical_fragments)},
        {"Item": "Fragment MS1 matches", "Value": len(fragment_ms1_matches)},
        {"Item": "Truncated sheets", "Value": _truncation_summary(truncations)},
        {"Item": "Warnings", "Value": len(warnings)},
    ]

    sheets: dict[str, pd.DataFrame] = {
        "Run_summary": pd.DataFrame(summary_rows),
        **truncated_data_sheets,
        "Warnings": pd.DataFrame(warnings, columns=["Timestamp", "Level", "Source", "Message", "Context"]),
    }
    sheets = {
        sheet_name: _truncate_frame_if_needed(
            sheet_name,
            frame,
            max_excel_rows,
            truncate_large_sheets,
            warnings,
            truncations,
        )
        if sheet_name == "Warnings"
        else frame
        for sheet_name, frame in sheets.items()
    }

    index_rows = [
        {
            "Sheet": sheet_name,
            "Description": SHEET_DESCRIPTIONS.get(sheet_name, "Optional result sheet."),
            "Notes": "Data starts at A3.",
        }
        for sheet_name in sheets
    ]

    with pd.ExcelWriter(report_path, engine="openpyxl") as writer:
        pd.DataFrame(index_rows, columns=["Sheet", "Description", "Notes"]).to_excel(writer, sheet_name="Index", index=False)
        for sheet_name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=sheet_name, index=False, startrow=2)
        _add_index_and_backlinks(writer, list(sheets))
        _autosize_and_freeze(writer)
    return report_path
