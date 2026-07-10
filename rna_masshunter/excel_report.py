from dataclasses import asdict, is_dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
from openpyxl.utils import get_column_letter

from rna_masshunter.ms2_annotation import (
    MS2_FRAGMENT_EVIDENCE_COLUMNS,
    MS2_ION_MATCH_COLUMNS,
    MS2_MODIFIED_PRECURSOR_COLUMNS,
    MS2_MODIFIED_THEORETICAL_ION_COLUMNS,
    MS2_MODIFIED_ION_MATCH_COLUMNS,
    MS2_LOCALIZATION_EVIDENCE_COLUMNS,
    MS2_PARENT_CANDIDATE_COLUMNS,
    MS2_SPECTRA_COLUMNS,
    MS2_SUMMARY_COLUMNS,
    MS2_THEORETICAL_ION_COLUMNS,
    MS2_UNMATCHED_COLUMNS,
)
from rna_masshunter.evidence_ranking import RANKING_COLUMNS, SUMMARY_COLUMNS
from rna_masshunter.p1_annotation import (
    P1_ANNOTATION_COLUMNS,
    P1_SUMMARY_COLUMNS,
    P1_THEORETICAL_COLUMNS,
    P1_UNMATCHED_COLUMNS,
)


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

FRAGMENT_MS1_FILTERED_COLUMNS = [
    "Match_ID",
    "Fragment_ID",
    "Target_ID",
    "Sequence",
    "Length",
    "Start",
    "End",
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

FRAGMENT_MS1_SUMMARY_COLUMNS = [
    "Fragment_ID",
    "Target_ID",
    "Sequence",
    "Length",
    "Start",
    "End",
    "Enzyme",
    "Missed_Cleavages",
    "Terminal_Form",
    "Best_Charge",
    "Best_Theoretical_mz",
    "Best_Observed_mz",
    "Best_Mass_Error_ppm",
    "Best_Intensity",
    "Best_RT",
    "Best_Peak_Tier",
    "Best_Confidence",
    "Match_Count",
    "Major_Count",
    "Minor_Count",
    "Trace_Count",
    "High_Count",
    "Medium_Count",
    "Low_Count",
]

KNOWN_MODIFICATION_CANDIDATE_COLUMNS = [
    "candidate_id",
    "source_type",
    "source_id",
    "target_id",
    "sequence",
    "start",
    "end",
    "observed_mz",
    "theoretical_mz",
    "observed_mass",
    "unmodified_mass",
    "mass_error_unmodified_da",
    "mass_error_unmodified_ppm",
    "modification_id",
    "modification_symbol",
    "modification_name",
    "target_base",
    "modification_mass_shift",
    "modified_mass",
    "mass_error_modified_da",
    "mass_error_modified_ppm",
    "charge",
    "intensity",
    "rt",
    "peak_tier",
    "confidence",
    "priority_score",
    "notes",
    "warnings",
]

KNOWN_MODIFICATION_SUMMARY_COLUMNS = [
    "Modification_ID",
    "Modification_Name",
    "Symbol",
    "Target_Base",
    "Candidate_Count",
    "Best_Source_ID",
    "Best_Sequence",
    "Best_Mass_Error_Modified_ppm",
    "Best_Intensity",
    "Best_Peak_Tier",
    "Best_Confidence",
    "Best_Priority_Score",
]

SHEET_DESCRIPTIONS = {
    "Run_summary": "Run-level summary for this RNA_MassHunter MVP-3 report.",
    "Input_parameters": "Flattened parameters loaded from config.yaml.",
    "mzML_diagnostics": "mzML scan counts, ranges, precursor metadata, and warnings.",
    "Intact_mass_reconstruction": "Reconstructed intact mass clusters and mass errors.",
    "Charge_state_peaks": "Peak and charge-state evidence supporting reconstructed masses.",
    "Theoretical_fragments": "Theoretical RNase digestion fragments and terminal forms.",
    "Fragment_MS1_matches": "MS1 peak matches for unmodified theoretical fragments.",
    "Fragment_MS1_filtered": "Filtered MS1 fragment matches for practical review.",
    "Fragment_MS1_summary": "Best MS1 match per fragment with match counts.",
    "Known_Modification_Candidates": "Known modification candidates explaining fragment or intact mass shifts.",
    "Known_Modification_Summary": "Grouped summary of known modification candidates.",
    "Modification_Evidence_Summary": "Run-level counts for integrated modification evidence ranking.",
    "Modification_Evidence_Ranking": "Integrated evidence scores for prioritizing modification candidates.",
    "P1_Summary": "Summary of P1 observed peak annotation results.",
    "P1_Theoretical_Structures": "P1 monomer and short oligonucleotide theoretical structure candidates.",
    "P1_Peak_Annotations": "Observed P1 peaks matched to theoretical structure candidates, retaining unmatched peaks.",
    "P1_Unmatched_Peaks": "Observed P1 peaks outside tolerance retained for unknown/adduct/phosphate review.",
    "MS2_Summary": "Run-level summary of MS2 c/y ion annotation.",
    "MS2_Spectra": "MS2 spectrum metadata, peak counts, and annotation status.",
    "MS2_Parent_Candidates": "Precursor m/z matches between MS2 spectra and theoretical digestion fragments.",
    "MS2_Theoretical_Ions": "Theoretical c/y RNA fragment ions generated from digestion fragments.",
    "MS2_Ion_Matches": "Matched observed MS2 peaks only; unmatched peaks are reported separately.",
    "MS2_Unmatched_Peaks": "Observed MS2 peaks outside tolerance retained for review.",
    "MS2_Fragment_Evidence": "Spectrum-parent fragment evidence summary from matched MS2 ions.",
    "MS2_Peak_Annotations": "Optional all-peak MS2 annotation sheet, disabled by default.",
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


def _analysis_mode(config) -> str:
    reconstruction_enabled = _as_bool((config.reconstruction or {}).get("enabled"), True)
    digestion_enabled = _as_bool((config.digestion or {}).get("enabled"), True)
    if reconstruction_enabled and digestion_enabled:
        return "Intact + digested fragment analysis"
    if reconstruction_enabled:
        return "Intact reconstruction only"
    if digestion_enabled:
        return "Digested fragment MS1 mapping"
    return "No active mass analysis"


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
                "Enzyme": raw.get("enzyme"),
                "Missed_Cleavages": raw.get("missed_cleavages"),
                "Terminal_Form": raw.get("terminal_form"),
                "Unmodified_Mass": raw.get("unmodified_mass"),
                "Warnings": fragment_warnings,
            }
        )
    return rows


def _match_raw(item: Any) -> dict[str, Any]:
    return asdict(item) if is_dataclass(item) else dict(item)


def _normalize_filter_values(values: Any) -> set[str]:
    if values is None:
        return set()
    if isinstance(values, str):
        return {values.lower()}
    return {str(value).lower() for value in values}


def _fragment_ms1_match_rows(fragment_ms1_matches: list[Any], include_length: bool = False) -> list[dict[str, Any]]:
    rows = []
    for item in fragment_ms1_matches:
        raw = _match_raw(item)
        match_warnings = raw.get("warnings", [])
        if isinstance(match_warnings, list):
            match_warnings = "; ".join(map(str, match_warnings))
        row = {
            "Match_ID": raw.get("match_id"),
            "Fragment_ID": raw.get("fragment_id"),
            "Target_ID": raw.get("target_id"),
            "Sequence": raw.get("sequence"),
            "Start": raw.get("start"),
            "End": raw.get("end"),
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
        if include_length:
            row["Length"] = len(raw.get("sequence") or "")
        rows.append(row)
    return rows


def _filter_fragment_ms1_matches(fragment_ms1_matches: list[Any], mapping_config: dict[str, Any]) -> list[Any]:
    min_length = _as_positive_int(mapping_config.get("min_fragment_length_for_filtered"), 3)
    allowed_tiers = _normalize_filter_values(mapping_config.get("filtered_peak_tiers", ["Major", "Minor"]))
    allowed_confidence = _normalize_filter_values(mapping_config.get("filtered_confidence", ["High", "Medium"]))
    filtered = []
    for item in fragment_ms1_matches:
        raw = _match_raw(item)
        if len(raw.get("sequence") or "") < min_length:
            continue
        if allowed_tiers and str(raw.get("peak_tier") or "").lower() not in allowed_tiers:
            continue
        if allowed_confidence and str(raw.get("confidence") or "").lower() not in allowed_confidence:
            continue
        filtered.append(item)
    return filtered


def _confidence_rank(value: Any) -> int:
    return {"high": 3, "medium": 2, "low": 1}.get(str(value or "").lower(), 0)


def _peak_tier_rank(value: Any) -> int:
    return {"major": 3, "minor": 2, "trace": 1}.get(str(value or "").lower(), 0)


def _best_match_sort_key(item: Any) -> tuple[int, int, float, float]:
    raw = _match_raw(item)
    return (
        -_confidence_rank(raw.get("confidence")),
        -_peak_tier_rank(raw.get("peak_tier")),
        abs(float(raw.get("mass_error_ppm") or 0.0)),
        -float(raw.get("intensity") or 0.0),
    )


def _fragment_ms1_summary_rows(fragment_ms1_matches: list[Any], mapping_config: dict[str, Any]) -> list[dict[str, Any]]:
    group_key = str(mapping_config.get("summary_best_match_by", "fragment_id") or "fragment_id")
    if group_key != "fragment_id":
        group_key = "fragment_id"

    grouped: dict[str, list[Any]] = {}
    for item in fragment_ms1_matches:
        raw = _match_raw(item)
        key = str(raw.get(group_key) or "")
        if not key:
            continue
        grouped.setdefault(key, []).append(item)

    rows = []
    for fragment_id, matches in grouped.items():
        best = min(matches, key=_best_match_sort_key)
        best_raw = _match_raw(best)
        tier_counts = {"Major": 0, "Minor": 0, "Trace": 0}
        confidence_counts = {"High": 0, "Medium": 0, "Low": 0}
        for item in matches:
            raw = _match_raw(item)
            tier = str(raw.get("peak_tier") or "")
            confidence = str(raw.get("confidence") or "")
            if tier in tier_counts:
                tier_counts[tier] += 1
            if confidence in confidence_counts:
                confidence_counts[confidence] += 1
        rows.append(
            {
                "Fragment_ID": fragment_id,
                "Target_ID": best_raw.get("target_id"),
                "Sequence": best_raw.get("sequence"),
                "Length": len(best_raw.get("sequence") or ""),
                "Start": best_raw.get("start"),
                "End": best_raw.get("end"),
                "Enzyme": best_raw.get("enzyme"),
                "Missed_Cleavages": best_raw.get("missed_cleavages"),
                "Terminal_Form": best_raw.get("terminal_form"),
                "Best_Charge": best_raw.get("charge"),
                "Best_Theoretical_mz": best_raw.get("theoretical_mz"),
                "Best_Observed_mz": best_raw.get("observed_mz"),
                "Best_Mass_Error_ppm": best_raw.get("mass_error_ppm"),
                "Best_Intensity": best_raw.get("intensity"),
                "Best_RT": best_raw.get("rt"),
                "Best_Peak_Tier": best_raw.get("peak_tier"),
                "Best_Confidence": best_raw.get("confidence"),
                "Match_Count": len(matches),
                "Major_Count": tier_counts["Major"],
                "Minor_Count": tier_counts["Minor"],
                "Trace_Count": tier_counts["Trace"],
                "High_Count": confidence_counts["High"],
                "Medium_Count": confidence_counts["Medium"],
                "Low_Count": confidence_counts["Low"],
            }
        )
    return sorted(rows, key=lambda row: (row["Start"] or 0, row["End"] or 0, row["Fragment_ID"]))


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
    known_modification_candidates: list[dict[str, Any]] | None = None,
    known_modification_summary: list[dict[str, Any]] | None = None,
    optional_results: dict[str, Any] | None = None,
) -> Path:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / f"RNA_MassHunter_MVP5_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"

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
    known_modification_candidates = known_modification_candidates or []
    known_modification_summary = known_modification_summary or []
    fragment_ms1_filtered = _filter_fragment_ms1_matches(fragment_ms1_matches, config.fragment_mapping or {})
    fragment_ms1_summary_rows = _fragment_ms1_summary_rows(fragment_ms1_matches, config.fragment_mapping or {})
    reconstruction_enabled = _as_bool(config.reconstruction.get("enabled"), True)

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
        "modification_search": config.modification_search,
        "peak_filtering": config.peak_filtering,
        "p1_annotation": config.p1_annotation,
        "ms2_annotation": config.ms2_annotation,
        "modification_evidence_ranking": config.modification_evidence_ranking,
        "performance": config.performance,
        "reporting": config.reporting,
    }

    data_sheets: dict[str, pd.DataFrame] = {
        "Input_parameters": pd.DataFrame(_flatten_dict(input_parameters)),
        "mzML_diagnostics": pd.DataFrame([diagnostics] if diagnostics else [{}]),
        "Theoretical_fragments": pd.DataFrame(_fragment_rows(theoretical_fragments), columns=THEORETICAL_FRAGMENT_COLUMNS),
        "Fragment_MS1_matches": pd.DataFrame(_fragment_ms1_match_rows(fragment_ms1_matches), columns=FRAGMENT_MS1_MATCH_COLUMNS),
        "Fragment_MS1_filtered": pd.DataFrame(_fragment_ms1_match_rows(fragment_ms1_filtered, include_length=True), columns=FRAGMENT_MS1_FILTERED_COLUMNS),
        "Fragment_MS1_summary": pd.DataFrame(fragment_ms1_summary_rows, columns=FRAGMENT_MS1_SUMMARY_COLUMNS),
        "Known_Modification_Candidates": pd.DataFrame(known_modification_candidates, columns=KNOWN_MODIFICATION_CANDIDATE_COLUMNS),
        "Known_Modification_Summary": pd.DataFrame(known_modification_summary, columns=KNOWN_MODIFICATION_SUMMARY_COLUMNS),
    }
    if reconstruction_enabled:
        data_sheets = {
            "Input_parameters": data_sheets["Input_parameters"],
            "mzML_diagnostics": data_sheets["mzML_diagnostics"],
            "Intact_mass_reconstruction": pd.DataFrame(intact_rows, columns=INTACT_COLUMNS),
            "Charge_state_peaks": pd.DataFrame(charge_state_peak_rows, columns=CHARGE_COLUMNS),
            **{key: value for key, value in data_sheets.items() if key not in {"Input_parameters", "mzML_diagnostics"}},
        }
    optional_columns = {
        "P1_Summary": P1_SUMMARY_COLUMNS,
        "P1_Theoretical_Structures": P1_THEORETICAL_COLUMNS,
        "P1_Peak_Annotations": P1_ANNOTATION_COLUMNS,
        "P1_Unmatched_Peaks": P1_UNMATCHED_COLUMNS,
        "MS2_Summary": MS2_SUMMARY_COLUMNS,
        "MS2_Spectra": MS2_SPECTRA_COLUMNS,
        "MS2_Parent_Candidates": MS2_PARENT_CANDIDATE_COLUMNS,
        "MS2_Modified_Precursor_Candidates": MS2_MODIFIED_PRECURSOR_COLUMNS,
        "MS2_Modified_Theoretical_Ions": MS2_MODIFIED_THEORETICAL_ION_COLUMNS,
        "MS2_Modified_Ion_Matches": MS2_MODIFIED_ION_MATCH_COLUMNS,
        "MS2_Modification_Localization_Evidence": MS2_LOCALIZATION_EVIDENCE_COLUMNS,
        "Modification_Evidence_Summary": SUMMARY_COLUMNS,
        "Modification_Evidence_Ranking": RANKING_COLUMNS,
        "MS2_Theoretical_Ions": MS2_THEORETICAL_ION_COLUMNS,
        "MS2_Ion_Matches": MS2_ION_MATCH_COLUMNS,
        "MS2_Unmatched_Peaks": MS2_UNMATCHED_COLUMNS,
        "MS2_Fragment_Evidence": MS2_FRAGMENT_EVIDENCE_COLUMNS,
        "MS2_Peak_Annotations": MS2_ION_MATCH_COLUMNS,
    }
    for sheet_name, value in (optional_results or {}).items():
        if sheet_name in {"Index", "Run_summary", "Warnings"}:
            continue
        frame = _coerce_to_frame(value)
        columns = optional_columns.get(sheet_name)
        if columns:
            frame = pd.DataFrame(frame, columns=columns)
        data_sheets[sheet_name[:31]] = frame

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
        {"Item": "Analysis mode", "Value": _analysis_mode(config)},
        {"Item": "Modification dictionary entries", "Value": len(modifications or [])},
        {"Item": "Rule set", "Value": config.organism.get("rule_set") or (rule_set or {}).get("id") or (rule_set or {}).get("name")},
        {"Item": "Pathway files", "Value": len(pathways or [])},
        {"Item": "Intact mass candidates", "Value": len(intact_results)},
        {"Item": "Theoretical fragments", "Value": len(theoretical_fragments)},
        {"Item": "Fragment MS1 matches", "Value": len(fragment_ms1_matches)},
        {"Item": "Fragment MS1 filtered", "Value": len(fragment_ms1_filtered)},
        {"Item": "Fragment MS1 summary", "Value": len(fragment_ms1_summary_rows)},
        {"Item": "Known modification candidates", "Value": len(known_modification_candidates)},
        {"Item": "Known modification summary", "Value": len(known_modification_summary)},
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
