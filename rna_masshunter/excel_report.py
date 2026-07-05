from dataclasses import asdict, is_dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
from openpyxl.utils import get_column_letter


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

SHEET_DESCRIPTIONS = {
    "Run_summary": "Run-level summary for this RNA_MassHunter MVP-1 report.",
    "Input_parameters": "Flattened parameters loaded from config.yaml.",
    "mzML_diagnostics": "mzML scan counts, ranges, precursor metadata, and warnings.",
    "Intact_mass_reconstruction": "Reconstructed intact mass clusters and mass errors.",
    "Charge_state_peaks": "Peak and charge-state evidence supporting reconstructed masses.",
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
    optional_results: dict[str, Any] | None = None,
) -> Path:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / f"RNA_MassHunter_MVP1_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"

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

    summary_rows = [
        {"Item": "Project", "Value": config.project.get("name")},
        {"Item": "Generated", "Value": datetime.now().isoformat(timespec="seconds")},
        {"Item": "Modification dictionary entries", "Value": len(modifications or [])},
        {"Item": "Rule set", "Value": (rule_set or {}).get("id") or (rule_set or {}).get("name")},
        {"Item": "Pathway files", "Value": len(pathways or [])},
        {"Item": "Intact mass candidates", "Value": len(intact_results)},
        {"Item": "Warnings", "Value": len(warnings)},
    ]

    input_parameters = {
        "project": config.project,
        "input": config.input,
        "organism": config.organism,
        "sequence": config.sequence,
        "experiment": config.experiment,
        "instrument": config.instrument,
        "reconstruction": config.reconstruction,
        "peak_filtering": config.peak_filtering,
        "performance": config.performance,
        "reporting": config.reporting,
    }

    sheets: dict[str, pd.DataFrame] = {
        "Run_summary": pd.DataFrame(summary_rows),
        "Input_parameters": pd.DataFrame(_flatten_dict(input_parameters)),
        "mzML_diagnostics": pd.DataFrame([diagnostics] if diagnostics else [{}]),
        "Intact_mass_reconstruction": pd.DataFrame(intact_rows, columns=INTACT_COLUMNS),
        "Charge_state_peaks": pd.DataFrame(charge_state_peaks, columns=CHARGE_COLUMNS),
        "Warnings": pd.DataFrame(warnings, columns=["Timestamp", "Level", "Source", "Message", "Context"]),
    }
    for sheet_name, value in (optional_results or {}).items():
        if sheet_name == "Index":
            continue
        sheets[sheet_name[:31]] = _coerce_to_frame(value)

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
