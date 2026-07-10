#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
from openpyxl.utils import get_column_letter

COMPARISON_NOTE = (
    "Comparison report highlights condition-dependent review patterns; it does not confirm modification identity or causality."
)

SUMMARY_COLUMNS = [
    "Total_Conditions",
    "Conditions",
    "Total_Unique_Candidates",
    "Candidates_Present_In_All_Conditions",
    "Candidates_Present_In_One_Condition",
    "Candidates_With_Review_Priority_Change",
    "Candidates_With_Confidence_Change",
    "Candidates_With_Ambiguity_Status_Change",
    "Top_Condition_Specific_Candidates",
    "Notes",
]

CANDIDATE_COLUMNS = [
    "Comparison_Key",
    "Modification_ID",
    "Modification_Name",
    "Parent_Fragment_ID",
    "Parent_Sequence",
    "Conditions_Present",
    "Num_Conditions_Present",
    "Best_Review_Priority",
    "Best_Final_Confidence",
    "Best_Final_Score",
    "Review_Priority_By_Condition",
    "Final_Confidence_By_Condition",
    "Final_Score_By_Condition",
    "Candidate_Positions_By_Condition",
    "Ambiguity_Status_By_Condition",
    "Evidence_Summary_By_Condition",
    "Key_Warnings_By_Condition",
    "Suggested_Comparison_Interpretation",
]

PRIORITY_CHANGE_COLUMNS = [
    "Comparison_Key",
    "Modification_ID",
    "Parent_Fragment_ID",
    "Review_Priority_By_Condition",
    "Final_Confidence_By_Condition",
    "Final_Score_By_Condition",
    "Priority_Change",
    "Confidence_Change",
    "Score_Range",
    "Notes",
]

AMBIGUITY_COLUMNS = [
    "Comparison_Key",
    "Modification_ID",
    "Parent_Fragment_ID",
    "Ambiguity_Status_By_Condition",
    "Position_Resolution_Basis_By_Condition",
    "Candidate_Positions_By_Condition",
    "Ambiguity_Change",
    "Notes",
]

DELTA_COLUMNS = [
    "Condition",
    "Num_Candidates",
    "Num_Unique_To_Condition",
    "Num_A_Strong_Review",
    "Num_B_Medium_Review",
    "Num_C_Ambiguous_Review",
    "Top_Unique_Candidates",
    "Top_Scoring_Candidates",
    "Notes",
]

SHEET_PRIORITY = [
    "Top_Modification_Candidates",
    "Modification_Evidence_Ranking",
    "Candidate_Decision_Summary",
    "Evidence_Checklist",
]

PRIORITY_ORDER = {
    "A_strong_review": 1,
    "B_medium_review": 2,
    "C_ambiguous_review": 3,
    "D_weak_review": 4,
    "E_low_information": 5,
}

CONFIDENCE_ORDER = {"Very_High": 5, "High": 4, "Medium": 3, "Low": 2, "Very_Low": 1}


@dataclass
class InputReport:
    condition: str
    path: Path
    candidates: pd.DataFrame


def parse_input(value: str) -> tuple[str | None, Path]:
    if "=" in value:
        condition, raw_path = value.split("=", 1)
        condition = condition.strip() or None
    else:
        condition, raw_path = None, value
    path = Path(raw_path.strip()).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    return condition, path


def normalize_condition(condition: str | None, path: Path, used: set[str]) -> str:
    base = condition or path.stem
    safe = "".join(ch if ch.isalnum() or ch in "_-" else "_" for ch in base).strip("_") or "condition"
    name = safe
    index = 2
    while name in used:
        name = f"{safe}_{index}"
        index += 1
    used.add(name)
    return name


def read_sheet(path: Path, sheet_name: str) -> pd.DataFrame:
    candidate_headers = [2, 0]
    for header in candidate_headers:
        try:
            frame = pd.read_excel(path, sheet_name=sheet_name, header=header)
        except ValueError:
            return pd.DataFrame()
        if frame.empty:
            continue
        columns = {str(column).lower() for column in frame.columns}
        if {"modification_id", "parent_fragment_id"} & columns:
            return frame
    return frame if "frame" in locals() else pd.DataFrame()


def first_existing(row: pd.Series, names: list[str], default: Any = "") -> Any:
    lowered = {str(key).lower(): key for key in row.index}
    for name in names:
        key = lowered.get(name.lower())
        if key is None:
            continue
        value = row.get(key)
        if pd.notna(value) and value != "":
            return value
    return default


def first_column(frame: pd.DataFrame, names: list[str]) -> str | None:
    lowered = {str(key).lower(): key for key in frame.columns}
    for name in names:
        key = lowered.get(name.lower())
        if key is not None:
            return key
    return None


def text(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value)


def number(value: Any, default: float = 0.0) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def join_values(values: list[Any]) -> str:
    seen: list[str] = []
    for value in values:
        item = text(value).strip()
        if item and item not in seen:
            seen.append(item)
    return "; ".join(seen)


def by_condition(items: dict[str, Any]) -> str:
    return "; ".join(f"{condition}={text(value)}" for condition, value in items.items() if text(value) != "")


def priority_rank(priority: Any) -> int:
    return PRIORITY_ORDER.get(text(priority), 99)


def confidence_rank(confidence: Any) -> int:
    raw = text(confidence).replace("-", "_").replace(" ", "_")
    normalized = "_".join(part.capitalize() for part in raw.split("_") if part)
    return CONFIDENCE_ORDER.get(normalized, 0)


def comparison_key(row: pd.Series, key_columns: list[str]) -> str:
    parts = []
    for column in key_columns:
        parts.append(text(first_existing(row, [column])).strip())
    return "|".join(parts)


def normalize_candidates(path: Path, include_low_confidence: bool, min_score: float, key_columns: list[str]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    source_sheet = ""
    for sheet_name in SHEET_PRIORITY:
        frame = read_sheet(path, sheet_name)
        if not frame.empty:
            source_sheet = sheet_name
            frames.append(frame)
            break

    if not frames:
        return pd.DataFrame(columns=["Comparison_Key", *key_columns])

    source = frames[0].copy()
    rows: list[dict[str, Any]] = []
    for _, row in source.iterrows():
        mod_id = first_existing(row, ["Modification_ID", "Modification", "modification_id"])
        parent_id = first_existing(row, ["Parent_Fragment_ID", "Fragment_ID", "parent_fragment_id"])
        if not text(mod_id).strip() and not text(parent_id).strip():
            continue
        final_score = number(first_existing(row, ["Best_Final_Score", "Final_Score", "Score", "Ranking_Score"], 0))
        confidence = first_existing(row, ["Best_Final_Confidence", "Final_Confidence", "Confidence"])
        if not include_low_confidence and confidence_rank(confidence) <= CONFIDENCE_ORDER["Low"]:
            continue
        if final_score < min_score:
            continue
        normalized = {
            "Modification_ID": mod_id,
            "Modification_Name": first_existing(row, ["Modification_Name", "Name", "modification_name"]),
            "Parent_Fragment_ID": parent_id,
            "Parent_Sequence": first_existing(row, ["Parent_Sequence", "Sequence", "parent_sequence"]),
            "Ambiguity_Group_ID": first_existing(row, ["Ambiguity_Group_ID", "Group_ID", "ambiguity_group_id"]),
            "Candidate_Positions_In_tRNA": first_existing(row, ["Candidate_Positions_In_tRNA", "Positions_In_tRNA", "tRNA_Position"]),
            "Review_Priority": first_existing(row, ["Review_Priority"], "E_low_information"),
            "Final_Confidence": confidence,
            "Final_Score": final_score,
            "Position_Ambiguity_Status": first_existing(row, ["Position_Ambiguity_Status", "Ambiguity_Status", "Position_Status"]),
            "Position_Resolution_Basis": first_existing(row, ["Position_Resolution_Basis", "Resolution_Basis"]),
            "Evidence_Summary": first_existing(row, ["Evidence_Summary", "Evidence_For", "Summary", "Evidence"]),
            "Key_Warnings": join_values(
                [
                    first_existing(row, ["Key_Warnings", "Warnings"]),
                    first_existing(row, ["Evidence_Against"]),
                    first_existing(row, ["Confidence_Limiting_Factors", "Limiting_Factors"]),
                    first_existing(row, ["Near_Isobaric_Warning"]),
                ]
            ),
            "Source_Sheet": source_sheet,
        }
        normalized["Comparison_Key"] = comparison_key(pd.Series(normalized), key_columns)
        rows.append(normalized)

    if not rows:
        return pd.DataFrame(columns=["Comparison_Key", *key_columns])
    frame = pd.DataFrame(rows)
    frame["_priority_rank"] = frame["Review_Priority"].map(priority_rank).fillna(99)
    frame = frame.sort_values(["Comparison_Key", "_priority_rank", "Final_Score"], ascending=[True, True, False])
    frame = frame.drop_duplicates("Comparison_Key", keep="first").drop(columns=["_priority_rank"])
    return frame


def load_reports(inputs: list[str], include_low_confidence: bool, min_score: float, key_columns: list[str]) -> list[InputReport]:
    reports = []
    used: set[str] = set()
    for value in inputs:
        condition, path = parse_input(value)
        name = normalize_condition(condition, path, used)
        if not path.exists():
            raise FileNotFoundError(f"Input Excel not found for {name}: {path}")
        reports.append(
            InputReport(
                condition=name,
                path=path,
                candidates=normalize_candidates(path, include_low_confidence, min_score, key_columns),
            )
        )
    return reports


def candidate_groups(reports: list[InputReport]) -> dict[str, dict[str, pd.Series]]:
    grouped: dict[str, dict[str, pd.Series]] = defaultdict(dict)
    for report in reports:
        if report.candidates.empty:
            continue
        for _, row in report.candidates.iterrows():
            key = text(row.get("Comparison_Key"))
            if key:
                grouped[key][report.condition] = row
    return dict(grouped)


def best_value(rows: list[pd.Series], column: str, ranker) -> Any:
    values = [row.get(column, "") for row in rows if text(row.get(column, ""))]
    if not values:
        return ""
    return sorted(values, key=ranker)[0] if ranker is priority_rank else max(values, key=ranker)


def interpretation(condition_rows: dict[str, pd.Series], reports: list[InputReport]) -> str:
    present = len(condition_rows)
    total = len(reports)
    statuses = {text(row.get("Position_Ambiguity_Status")) for row in condition_rows.values() if text(row.get("Position_Ambiguity_Status"))}
    if not condition_rows:
        return "insufficient-comparison-evidence"
    if len(statuses) > 1:
        return "ambiguity-changes-between-conditions"
    if present == 1:
        return "condition-specific-review-candidate"
    if present == total:
        return "shared-across-conditions"
    scores = {condition: number(row.get("Final_Score")) for condition, row in condition_rows.items()}
    if scores and max(scores.values()) > min(scores.values()):
        return "stronger-in-treatment"
    return "shared-across-conditions"


def build_candidate_comparison(reports: list[InputReport]) -> pd.DataFrame:
    grouped = candidate_groups(reports)
    rows = []
    for key, condition_rows in grouped.items():
        row_list = list(condition_rows.values())
        conditions = [report.condition for report in reports if report.condition in condition_rows]
        scores = {condition: condition_rows[condition].get("Final_Score", "") for condition in conditions}
        priorities = {condition: condition_rows[condition].get("Review_Priority", "") for condition in conditions}
        confidences = {condition: condition_rows[condition].get("Final_Confidence", "") for condition in conditions}
        positions = {condition: condition_rows[condition].get("Candidate_Positions_In_tRNA", "") for condition in conditions}
        ambiguities = {condition: condition_rows[condition].get("Position_Ambiguity_Status", "") for condition in conditions}
        evidence = {condition: condition_rows[condition].get("Evidence_Summary", "") for condition in conditions}
        warnings = {condition: condition_rows[condition].get("Key_Warnings", "") for condition in conditions}
        rows.append(
            {
                "Comparison_Key": key,
                "Modification_ID": join_values([row.get("Modification_ID") for row in row_list]),
                "Modification_Name": join_values([row.get("Modification_Name") for row in row_list]),
                "Parent_Fragment_ID": join_values([row.get("Parent_Fragment_ID") for row in row_list]),
                "Parent_Sequence": join_values([row.get("Parent_Sequence") for row in row_list]),
                "Conditions_Present": "; ".join(conditions),
                "Num_Conditions_Present": len(conditions),
                "Best_Review_Priority": best_value(row_list, "Review_Priority", priority_rank),
                "Best_Final_Confidence": best_value(row_list, "Final_Confidence", confidence_rank),
                "Best_Final_Score": max(number(row.get("Final_Score")) for row in row_list),
                "Review_Priority_By_Condition": by_condition(priorities),
                "Final_Confidence_By_Condition": by_condition(confidences),
                "Final_Score_By_Condition": by_condition(scores),
                "Candidate_Positions_By_Condition": by_condition(positions),
                "Ambiguity_Status_By_Condition": by_condition(ambiguities),
                "Evidence_Summary_By_Condition": by_condition(evidence),
                "Key_Warnings_By_Condition": by_condition(warnings),
                "Suggested_Comparison_Interpretation": interpretation(condition_rows, reports),
            }
        )
    frame = pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)
    if frame.empty:
        return frame
    frame["_priority_rank"] = frame["Best_Review_Priority"].map(priority_rank).fillna(99)
    frame = frame.sort_values(["_priority_rank", "Best_Final_Score", "Num_Conditions_Present"], ascending=[True, False, False])
    return frame.drop(columns=["_priority_rank"])


def changed(values: list[Any]) -> bool:
    cleaned = {text(value) for value in values if text(value)}
    return len(cleaned) > 1


def build_presence_matrix(reports: list[InputReport], comparison: pd.DataFrame) -> pd.DataFrame:
    rows = []
    lookup = {report.condition: report.candidates.set_index("Comparison_Key") if not report.candidates.empty else pd.DataFrame() for report in reports}
    for _, candidate in comparison.iterrows():
        key = candidate["Comparison_Key"]
        row = {
            "Comparison_Key": key,
            "Modification_ID": candidate["Modification_ID"],
            "Parent_Fragment_ID": candidate["Parent_Fragment_ID"],
        }
        for report in reports:
            table = lookup[report.condition]
            present = not table.empty and key in table.index
            row[report.condition] = "present" if present else "absent"
        for report in reports:
            table = lookup[report.condition]
            row[f"Score_{report.condition}"] = number(table.loc[key].get("Final_Score")) if not table.empty and key in table.index else ""
        rows.append(row)
    base_columns = ["Comparison_Key", "Modification_ID", "Parent_Fragment_ID"]
    return pd.DataFrame(rows, columns=base_columns + [r.condition for r in reports] + [f"Score_{r.condition}" for r in reports])


def build_priority_changes(comparison: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in comparison.iterrows():
        priorities = [part.split("=", 1)[1] for part in text(row["Review_Priority_By_Condition"]).split("; ") if "=" in part]
        confidences = [part.split("=", 1)[1] for part in text(row["Final_Confidence_By_Condition"]).split("; ") if "=" in part]
        scores = [number(part.split("=", 1)[1]) for part in text(row["Final_Score_By_Condition"]).split("; ") if "=" in part]
        priority_change = changed(priorities)
        confidence_change = changed(confidences)
        if not priority_change and not confidence_change:
            continue
        rows.append(
            {
                "Comparison_Key": row["Comparison_Key"],
                "Modification_ID": row["Modification_ID"],
                "Parent_Fragment_ID": row["Parent_Fragment_ID"],
                "Review_Priority_By_Condition": row["Review_Priority_By_Condition"],
                "Final_Confidence_By_Condition": row["Final_Confidence_By_Condition"],
                "Final_Score_By_Condition": row["Final_Score_By_Condition"],
                "Priority_Change": priority_change,
                "Confidence_Change": confidence_change,
                "Score_Range": max(scores) - min(scores) if scores else 0,
                "Notes": COMPARISON_NOTE,
            }
        )
    return pd.DataFrame(rows, columns=PRIORITY_CHANGE_COLUMNS)


def build_ambiguity_comparison(reports: list[InputReport], comparison: pd.DataFrame) -> pd.DataFrame:
    grouped = candidate_groups(reports)
    rows = []
    for _, candidate in comparison.iterrows():
        key = candidate["Comparison_Key"]
        condition_rows = grouped.get(key, {})
        statuses = {condition: row.get("Position_Ambiguity_Status", "") for condition, row in condition_rows.items()}
        if not changed(list(statuses.values())):
            continue
        basis = {condition: row.get("Position_Resolution_Basis", "") for condition, row in condition_rows.items()}
        positions = {condition: row.get("Candidate_Positions_In_tRNA", "") for condition, row in condition_rows.items()}
        rows.append(
            {
                "Comparison_Key": key,
                "Modification_ID": candidate["Modification_ID"],
                "Parent_Fragment_ID": candidate["Parent_Fragment_ID"],
                "Ambiguity_Status_By_Condition": by_condition(statuses),
                "Position_Resolution_Basis_By_Condition": by_condition(basis),
                "Candidate_Positions_By_Condition": by_condition(positions),
                "Ambiguity_Change": True,
                "Notes": COMPARISON_NOTE,
            }
        )
    return pd.DataFrame(rows, columns=AMBIGUITY_COLUMNS)


def build_delta_summary(reports: list[InputReport], comparison: pd.DataFrame, top_n: int) -> pd.DataFrame:
    grouped = candidate_groups(reports)
    unique_keys_by_condition: dict[str, list[str]] = defaultdict(list)
    for key, condition_rows in grouped.items():
        if len(condition_rows) == 1:
            condition = next(iter(condition_rows))
            unique_keys_by_condition[condition].append(key)

    rows = []
    for report in reports:
        frame = report.candidates
        priorities = Counter(frame.get("Review_Priority", pd.Series(dtype=str)).fillna("").astype(str)) if not frame.empty else Counter()
        unique_keys = unique_keys_by_condition.get(report.condition, [])
        unique_labels = []
        for key in unique_keys[:top_n]:
            match = comparison[comparison["Comparison_Key"] == key]
            unique_labels.append(key if match.empty else text(match.iloc[0].get("Modification_ID")) or key)
        top_scoring = []
        if not frame.empty:
            for _, row in frame.sort_values("Final_Score", ascending=False).head(top_n).iterrows():
                label = text(row.get("Modification_ID")) or text(row.get("Comparison_Key"))
                score = number(row.get("Final_Score"))
                top_scoring.append(f"{label}:{score:g}")
        rows.append(
            {
                "Condition": report.condition,
                "Num_Candidates": len(frame),
                "Num_Unique_To_Condition": len(unique_keys),
                "Num_A_Strong_Review": priorities.get("A_strong_review", 0),
                "Num_B_Medium_Review": priorities.get("B_medium_review", 0),
                "Num_C_Ambiguous_Review": priorities.get("C_ambiguous_review", 0),
                "Top_Unique_Candidates": "; ".join(unique_labels),
                "Top_Scoring_Candidates": "; ".join(top_scoring),
                "Notes": COMPARISON_NOTE,
            }
        )
    return pd.DataFrame(rows, columns=DELTA_COLUMNS)


def build_summary(reports: list[InputReport], comparison: pd.DataFrame, priority_changes: pd.DataFrame, ambiguity: pd.DataFrame) -> pd.DataFrame:
    total = len(reports)
    present_counts = comparison.get("Num_Conditions_Present", pd.Series(dtype=int))
    condition_specific = comparison[comparison.get("Num_Conditions_Present", pd.Series(dtype=int)) == 1] if not comparison.empty else pd.DataFrame()
    top_specific = "; ".join(
        text(row.Modification_ID) or text(row.Comparison_Key)
        for row in condition_specific[["Comparison_Key", "Modification_ID"]].head(10).itertuples(index=False)
    ) if not condition_specific.empty else ""
    confidence_changes = 0
    if not priority_changes.empty:
        confidence_changes = int(priority_changes["Confidence_Change"].fillna(False).sum())
    row = {
        "Total_Conditions": total,
        "Conditions": "; ".join(report.condition for report in reports),
        "Total_Unique_Candidates": len(comparison),
        "Candidates_Present_In_All_Conditions": int((present_counts == total).sum()) if not comparison.empty else 0,
        "Candidates_Present_In_One_Condition": int((present_counts == 1).sum()) if not comparison.empty else 0,
        "Candidates_With_Review_Priority_Change": int(priority_changes["Priority_Change"].fillna(False).sum()) if not priority_changes.empty else 0,
        "Candidates_With_Confidence_Change": confidence_changes,
        "Candidates_With_Ambiguity_Status_Change": len(ambiguity),
        "Top_Condition_Specific_Candidates": top_specific,
        "Notes": COMPARISON_NOTE,
    }
    return pd.DataFrame([row], columns=SUMMARY_COLUMNS)


def autosize(writer: pd.ExcelWriter) -> None:
    for worksheet in writer.book.worksheets:
        worksheet.freeze_panes = "A2"
        for column_cells in worksheet.columns:
            width = max(10, min(60, max(len("" if cell.value is None else str(cell.value)) for cell in column_cells) + 2))
            worksheet.column_dimensions[get_column_letter(column_cells[0].column)].width = width


def write_report(output: Path, sheets: dict[str, pd.DataFrame]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        for sheet_name, frame in sheets.items():
            frame.to_excel(writer, sheet_name=sheet_name[:31], index=False)
        autosize(writer)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare RNA_MassHunter Excel review reports across conditions.")
    parser.add_argument("--input", action="append", required=True, help="condition=path.xlsx, or path.xlsx to use file stem as condition")
    parser.add_argument("--output", required=True, help="Output comparison workbook path")
    parser.add_argument("--key-columns", nargs="+", default=["Modification_ID", "Parent_Fragment_ID"], help="Candidate key columns")
    parser.add_argument("--include-low-confidence", action="store_true", default=True, help="Include low-confidence candidates")
    parser.add_argument("--exclude-low-confidence", dest="include_low_confidence", action="store_false", help="Exclude Low/Very_Low candidates")
    parser.add_argument("--min-score", type=float, default=0.0, help="Minimum final score to include")
    parser.add_argument("--top-n", type=int, default=50, help="Maximum candidates to keep in ranked comparison outputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if len(args.input) < 2:
        raise SystemExit("At least two --input workbooks are required.")
    output = Path(args.output).expanduser()
    if not output.is_absolute():
        output = Path.cwd() / output
    reports = load_reports(args.input, args.include_low_confidence, args.min_score, args.key_columns)
    comparison = build_candidate_comparison(reports)
    presence = build_presence_matrix(reports, comparison)
    priority_changes = build_priority_changes(comparison)
    ambiguity = build_ambiguity_comparison(reports, comparison)
    delta = build_delta_summary(reports, comparison, args.top_n)
    summary = build_summary(reports, comparison, priority_changes, ambiguity)
    write_report(
        output,
        {
            "Comparison_Summary": summary,
            "Candidate_Comparison": comparison,
            "Condition_Presence_Matrix": presence,
            "Review_Priority_Changes": priority_changes,
            "Ambiguity_Comparison": ambiguity,
            "Candidate_Delta_Summary": delta,
        },
    )
    print(f"Comparison report written: {output}")


if __name__ == "__main__":
    main()
