"""Safe common parser and diagnostics for SCIEX profile text exports."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Any
import math
import re
import unicodedata

NEUTRAL_MASS_PROFILE = "NEUTRAL_MASS_PROFILE"
MZ_PROFILE = "MZ_PROFILE"
UNKNOWN_PROFILE = "UNKNOWN"

FORMAL_FALSE = {
    "Applied_To_Formal_Result": False,
    "Formal_Change_Ready": False,
    "Formal_Result_Changed": False,
}

DIAGNOSTIC_COLUMNS = [
    "Source_File", "Source_File_Name", "Header_Line_Number", "Raw_Headers",
    "Normalized_Headers", "Coordinate_Header", "Intensity_Header", "Profile_Type",
    "Expected_Profile_Type", "Input_Status", "Parsed_Row_Count", "Invalid_Row_Count",
    "Invalid_Row_Numbers", "Ignored_Blank_Line_Count", "Ignored_Comment_Line_Count",
    "Coordinate_Min", "Coordinate_Max", "Intensity_Min", "Intensity_Max",
    "Zero_Intensity_Count", "Duplicate_Coordinate_Count",
    "Descending_Transition_Count", "Strictly_Increasing", "Step_Count",
    "Step_Min", "Step_Max", "Step_Median", "Uniform_Step",
    "Eligible_For_Neutral_Mass_Analysis", "Neutral_Mass_Analysis_Skip_Reason",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]

INPUT_COLUMNS = [
    "Source_File", "Source_Row_Number", "Profile_Type", "Neutral_Mass", "MZ",
    "Intensity", "Coordinate_Duplicate", "Coordinate_Descending_From_Previous",
    "Eligible_For_Neutral_Mass_Analysis", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]


@dataclass(frozen=True)
class SCIEXProfileParseResult:
    profile_type: str
    expected_profile_type: str
    input_status: str
    diagnostic_rows: list[dict[str, Any]]
    input_rows: list[dict[str, Any]]

    @property
    def neutral_mass_analysis_eligible(self) -> bool:
        return self.input_status == "SUPPORTED_INPUT" and self.profile_type == NEUTRAL_MASS_PROFILE

    def sheets(self, audit_level: str = "full") -> dict[str, list[dict[str, Any]]]:
        if str(audit_level).lower() == "standard":
            return {}
        output = {"SCIEX_Profile_Diagnostics": self.diagnostic_rows}
        if str(audit_level).lower() == "full":
            output["SCIEX_Profile_Input"] = self.input_rows
        return output


def _normalize_header(value: str) -> str:
    normalized = unicodedata.normalize("NFKC", value).replace("\ufeff", "")
    return " ".join(normalized.strip().casefold().split())


def _canonical_header(value: str) -> str:
    normalized = _normalize_header(value)
    if normalized == "mass":
        return "MASS"
    if normalized == "intensity":
        return "INTENSITY"
    compact = normalized.replace(" ", "")
    if compact in {"mass/charge", "m/z", "mz", "mass-to-charge"}:
        return "MASS_TO_CHARGE"
    return "UNRECOGNIZED"


def _split_header(line: str) -> list[str]:
    stripped = line.strip()
    if "\t" in stripped:
        return [item.strip() for item in stripped.split("\t") if item.strip()]
    return [item.strip() for item in re.split(r"\s{2,}", stripped) if item.strip()]


def _split_data(line: str) -> list[str]:
    return [item for item in re.split(r"\s+", line.strip()) if item]


def _expected_profile(path: Path) -> str:
    stem = path.stem.casefold()
    if stem.endswith("-full"):
        return NEUTRAL_MASS_PROFILE
    if stem.endswith("-t1"):
        return MZ_PROFILE
    return UNKNOWN_PROFILE


def _status(profile_type: str, expected: str, invalid: int, parsed: int,
            recognized: bool, has_header: bool) -> str:
    if not has_header:
        return "EMPTY_INPUT"
    if not recognized:
        return "UNRECOGNIZED_COLUMNS"
    if parsed == 0 and invalid == 0:
        return "EMPTY_INPUT"
    if invalid:
        return "INVALID_NUMERIC_DATA"
    if expected != UNKNOWN_PROFILE and expected != profile_type:
        return "HEADER_FILENAME_CONFLICT"
    if profile_type == NEUTRAL_MASS_PROFILE:
        return "SUPPORTED_INPUT"
    return "UNSUPPORTED_PROFILE_TYPE"


def parse_sciex_profile(path: str | Path) -> SCIEXProfileParseResult:
    """Parse a SCIEX text profile without converting m/z to neutral mass."""
    source = Path(path)
    lines = source.read_text(encoding="utf-8-sig").splitlines()
    blank_count = 0
    comment_count = 0
    header_line_number = 0
    header_line = ""
    data_lines: list[tuple[int, str]] = []
    for number, line in enumerate(lines, 1):
        if not line.strip():
            blank_count += 1
            continue
        if line.lstrip().startswith("#"):
            comment_count += 1
            continue
        if not header_line:
            header_line = line
            header_line_number = number
        else:
            data_lines.append((number, line))

    raw_headers = _split_header(header_line) if header_line else []
    canonical = [_canonical_header(value) for value in raw_headers]
    neutral_header = canonical.count("MASS") == 1
    mz_header = canonical.count("MASS_TO_CHARGE") == 1
    intensity_header = canonical.count("INTENSITY") == 1
    recognized = (
        len(raw_headers) == 2 and intensity_header and neutral_header != mz_header
        and "UNRECOGNIZED" not in canonical
    )
    profile_type = (
        NEUTRAL_MASS_PROFILE if recognized and neutral_header else
        MZ_PROFILE if recognized and mz_header else UNKNOWN_PROFILE
    )
    expected = _expected_profile(source)
    coordinate_index = canonical.index("MASS") if neutral_header else canonical.index("MASS_TO_CHARGE") if mz_header else -1
    intensity_index = canonical.index("INTENSITY") if intensity_header else -1

    parsed: list[tuple[int, float, float]] = []
    invalid_rows: list[int] = []
    if recognized:
        for row_number, line in data_lines:
            fields = _split_data(line)
            if len(fields) != 2:
                invalid_rows.append(row_number)
                continue
            try:
                coordinate = float(fields[coordinate_index])
                intensity = float(fields[intensity_index])
            except (TypeError, ValueError):
                invalid_rows.append(row_number)
                continue
            if not math.isfinite(coordinate) or not math.isfinite(intensity) or coordinate < 0 or intensity < 0:
                invalid_rows.append(row_number)
                continue
            parsed.append((row_number, coordinate, intensity))

    status = _status(profile_type, expected, len(invalid_rows), len(parsed), recognized, bool(header_line))
    eligible = status == "SUPPORTED_INPUT" and profile_type == NEUTRAL_MASS_PROFILE
    coordinates = [item[1] for item in parsed]
    intensities = [item[2] for item in parsed]
    seen: set[float] = set()
    duplicate_flags: list[bool] = []
    descending_flags: list[bool] = []
    for index, coordinate in enumerate(coordinates):
        duplicate_flags.append(coordinate in seen)
        seen.add(coordinate)
        descending_flags.append(index > 0 and coordinate < coordinates[index - 1])
    steps = [coordinates[index] - coordinates[index - 1] for index in range(1, len(coordinates))]
    positive_steps = [step for step in steps if step > 0]
    uniform_step = False
    if positive_steps and not any(duplicate_flags) and not any(descending_flags):
        reference = median(positive_steps)
        tolerance = max(1e-12, abs(reference) * 1e-9)
        uniform_step = all(abs(step - reference) <= tolerance for step in positive_steps)
    skip_reason = "" if eligible else {
        "UNSUPPORTED_PROFILE_TYPE": "mz_profile_not_eligible_for_neutral_mass_analysis",
        "HEADER_FILENAME_CONFLICT": "header_filename_profile_type_conflict",
        "UNRECOGNIZED_COLUMNS": "profile_columns_unrecognized",
        "INVALID_NUMERIC_DATA": "invalid_numeric_profile_rows",
        "EMPTY_INPUT": "profile_has_no_data_rows",
    }.get(status, "input_not_supported")

    input_rows = []
    for (row_number, coordinate, intensity), duplicate, descending in zip(
        parsed, duplicate_flags, descending_flags, strict=False,
    ):
        input_rows.append({
            "Source_File": str(source), "Source_Row_Number": row_number,
            "Profile_Type": profile_type,
            "Neutral_Mass": coordinate if profile_type == NEUTRAL_MASS_PROFILE else "",
            "MZ": coordinate if profile_type == MZ_PROFILE else "",
            "Intensity": intensity, "Coordinate_Duplicate": duplicate,
            "Coordinate_Descending_From_Previous": descending,
            "Eligible_For_Neutral_Mass_Analysis": eligible, **FORMAL_FALSE,
        })
    diagnostics = [{
        "Source_File": str(source), "Source_File_Name": source.name,
        "Header_Line_Number": header_line_number or "", "Raw_Headers": ";".join(raw_headers),
        "Normalized_Headers": ";".join(canonical),
        "Coordinate_Header": raw_headers[coordinate_index] if coordinate_index >= 0 and coordinate_index < len(raw_headers) else "",
        "Intensity_Header": raw_headers[intensity_index] if intensity_index >= 0 and intensity_index < len(raw_headers) else "",
        "Profile_Type": profile_type, "Expected_Profile_Type": expected,
        "Input_Status": status, "Parsed_Row_Count": len(parsed),
        "Invalid_Row_Count": len(invalid_rows),
        "Invalid_Row_Numbers": ";".join(map(str, invalid_rows)),
        "Ignored_Blank_Line_Count": blank_count,
        "Ignored_Comment_Line_Count": comment_count,
        "Coordinate_Min": min(coordinates) if coordinates else "",
        "Coordinate_Max": max(coordinates) if coordinates else "",
        "Intensity_Min": min(intensities) if intensities else "",
        "Intensity_Max": max(intensities) if intensities else "",
        "Zero_Intensity_Count": sum(value == 0 for value in intensities),
        "Duplicate_Coordinate_Count": sum(duplicate_flags),
        "Descending_Transition_Count": sum(descending_flags),
        "Strictly_Increasing": bool(coordinates) and not any(duplicate_flags) and not any(descending_flags),
        "Step_Count": len(steps), "Step_Min": min(positive_steps) if positive_steps else "",
        "Step_Max": max(positive_steps) if positive_steps else "",
        "Step_Median": median(positive_steps) if positive_steps else "",
        "Uniform_Step": uniform_step,
        "Eligible_For_Neutral_Mass_Analysis": eligible,
        "Neutral_Mass_Analysis_Skip_Reason": skip_reason, **FORMAL_FALSE,
    }]
    return SCIEXProfileParseResult(profile_type, expected, status, diagnostics, input_rows)
