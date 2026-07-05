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
        worksheet.freeze_panes = "A2"
        for column_cells in worksheet.columns:
            max_length = 0
            column = get_column_letter(column_cells[0].column)
            for cell in column_cells:
                value = "" if cell.value is None else str(cell.value)
                max_length = max(max_length, min(len(value), 60))
            worksheet.column_dimensions[column].width = max(10, max_length + 2)


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

    with pd.ExcelWriter(report_path, engine="openpyxl") as writer:
        pd.DataFrame(summary_rows).to_excel(writer, sheet_name="Run_summary", index=False)
        pd.DataFrame(_flatten_dict(input_parameters)).to_excel(writer, sheet_name="Input_parameters", index=False)
        pd.DataFrame([diagnostics] if diagnostics else [{}]).to_excel(writer, sheet_name="mzML_diagnostics", index=False)
        pd.DataFrame(intact_rows, columns=INTACT_COLUMNS).to_excel(writer, sheet_name="Intact_mass_reconstruction", index=False)
        pd.DataFrame(charge_state_peaks, columns=CHARGE_COLUMNS).to_excel(writer, sheet_name="Charge_state_peaks", index=False)
        pd.DataFrame(warnings, columns=["Timestamp", "Level", "Source", "Message", "Context"]).to_excel(writer, sheet_name="Warnings", index=False)
        _autosize_and_freeze(writer)
    return report_path
