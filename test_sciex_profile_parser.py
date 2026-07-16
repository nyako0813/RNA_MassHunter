from pathlib import Path

import pytest

from rna_masshunter.audit_policy import AuditPolicy, included_sheet_names
from rna_masshunter.sciex_profile_parser import (
    MZ_PROFILE, NEUTRAL_MASS_PROFILE, parse_sciex_profile,
)


def write(tmp_path: Path, name: str, content: str, *, bom: bool = False) -> Path:
    path = tmp_path / name
    encoding = "utf-8-sig" if bom else "utf-8"
    path.write_text(content, encoding=encoding)
    return path


def test_neutral_mass_tab_supported_and_fields_separated(tmp_path):
    result = parse_sciex_profile(write(tmp_path, "sample-full.txt", "Mass\tIntensity\n25000\t100\n25001\t0\n"))
    assert result.profile_type == NEUTRAL_MASS_PROFILE
    assert result.input_status == "SUPPORTED_INPUT"
    assert result.neutral_mass_analysis_eligible is True
    assert result.input_rows[0]["Neutral_Mass"] == 25000
    assert result.input_rows[0]["MZ"] == ""
    assert result.diagnostic_rows[0]["Zero_Intensity_Count"] == 1


@pytest.mark.parametrize("header", ["Mass/Charge", "Mass / Charge", "m/z", "MZ", "Mass-to-Charge"])
def test_mz_aliases_are_not_neutral_mass(tmp_path, header):
    result = parse_sciex_profile(write(tmp_path, "sample-T1.txt", f"{header}\tIntensity\n500\t10\n"))
    assert result.profile_type == MZ_PROFILE
    assert result.input_status == "UNSUPPORTED_PROFILE_TYPE"
    assert result.neutral_mass_analysis_eligible is False
    assert result.input_rows[0]["Neutral_Mass"] == ""
    assert result.input_rows[0]["MZ"] == 500
    assert result.input_rows[0]["Eligible_For_Neutral_Mass_Analysis"] is False


def test_bom_case_whitespace_comments_and_blank_lines(tmp_path):
    path = write(tmp_path, "mixed-full.txt", "# comment\n\n  MASS  \t intensity  \n25000  10\n25001  20\n", bom=True)
    result = parse_sciex_profile(path)
    assert result.input_status == "SUPPORTED_INPUT"
    diagnostic = result.diagnostic_rows[0]
    assert diagnostic["Ignored_Comment_Line_Count"] == 1
    assert diagnostic["Ignored_Blank_Line_Count"] == 1
    assert diagnostic["Parsed_Row_Count"] == 2


def test_consecutive_space_header_and_data(tmp_path):
    result = parse_sciex_profile(write(tmp_path, "sample-full.txt", "Mass    Intensity\n25000    10\n"))
    assert result.input_status == "SUPPORTED_INPUT"


@pytest.mark.parametrize("name,header", [
    ("sample-full.txt", "Mass/Charge"),
    ("sample-T1.txt", "Mass"),
])
def test_filename_header_conflict(tmp_path, name, header):
    result = parse_sciex_profile(write(tmp_path, name, f"{header}\tIntensity\n500\t10\n"))
    assert result.input_status == "HEADER_FILENAME_CONFLICT"
    assert result.neutral_mass_analysis_eligible is False


@pytest.mark.parametrize("value", ["abc", "NaN", "inf", "-1"])
def test_invalid_coordinate_data(tmp_path, value):
    result = parse_sciex_profile(write(tmp_path, "sample-full.txt", f"Mass\tIntensity\n{value}\t10\n"))
    assert result.input_status == "INVALID_NUMERIC_DATA"
    assert result.diagnostic_rows[0]["Invalid_Row_Count"] == 1


@pytest.mark.parametrize("value", ["abc", "NaN", "-inf", "-1"])
def test_invalid_intensity_data(tmp_path, value):
    result = parse_sciex_profile(write(tmp_path, "sample-full.txt", f"Mass\tIntensity\n25000\t{value}\n"))
    assert result.input_status == "INVALID_NUMERIC_DATA"


def test_unrecognized_and_empty(tmp_path):
    unknown = parse_sciex_profile(write(tmp_path, "unknown.txt", "Time\tSignal\n1\t2\n"))
    assert unknown.input_status == "UNRECOGNIZED_COLUMNS"
    empty = parse_sciex_profile(write(tmp_path, "empty.txt", "# only comment\n\n"))
    assert empty.input_status == "EMPTY_INPUT"


def test_duplicate_and_descending_diagnostics(tmp_path):
    text = "Mass\tIntensity\n10\t1\n10\t2\n9\t3\n"
    result = parse_sciex_profile(write(tmp_path, "sample-full.txt", text))
    diagnostic = result.diagnostic_rows[0]
    assert diagnostic["Duplicate_Coordinate_Count"] == 1
    assert diagnostic["Descending_Transition_Count"] == 1
    assert diagnostic["Strictly_Increasing"] is False
    assert result.input_rows[1]["Coordinate_Duplicate"] is True
    assert result.input_rows[2]["Coordinate_Descending_From_Previous"] is True


def test_uniform_and_nonuniform_step_diagnostics(tmp_path):
    uniform = parse_sciex_profile(write(tmp_path, "u-full.txt", "Mass\tIntensity\n1\t1\n2\t1\n3\t1\n"))
    assert uniform.diagnostic_rows[0]["Uniform_Step"] is True
    assert uniform.diagnostic_rows[0]["Step_Median"] == 1
    nonuniform = parse_sciex_profile(write(tmp_path, "n-full.txt", "Mass\tIntensity\n1\t1\n2\t1\n4\t1\n"))
    assert nonuniform.diagnostic_rows[0]["Uniform_Step"] is False
    assert nonuniform.diagnostic_rows[0]["Step_Min"] == 1
    assert nonuniform.diagnostic_rows[0]["Step_Max"] == 2


def test_all_formal_flags_false(tmp_path):
    result = parse_sciex_profile(write(tmp_path, "sample-full.txt", "Mass\tIntensity\n1\t1\n"))
    for rows in (result.diagnostic_rows, result.input_rows):
        for row in rows:
            assert row["Applied_To_Formal_Result"] is False
            assert row["Formal_Change_Ready"] is False
            assert row["Formal_Result_Changed"] is False


def test_result_sheet_levels_and_policy():
    names = ["SCIEX_Profile_Diagnostics", "SCIEX_Profile_Input"]
    standard, _ = included_sheet_names(names, AuditPolicy.from_level("standard"))
    audit, _ = included_sheet_names(names, AuditPolicy.from_level("audit"))
    full, _ = included_sheet_names(names, AuditPolicy.from_level("full"))
    assert standard == [] and audit == [names[0]] and full == names
