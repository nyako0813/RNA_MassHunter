from types import SimpleNamespace

from openpyxl import load_workbook

from rna_masshunter.excel_report import write_excel_report
from rna_masshunter.intact_reconstruction import build_intact_reconstruction_qc
from rna_masshunter.models import IntactMassCandidate


BASE_QC_CONFIG = {
    "intact_reconstruction": {
        "min_charge_states_for_reliable": 3,
        "min_charge_states_for_review": 2,
        "require_contiguous_charge_states": True,
        "max_neutral_mass_sd_da": 0.5,
        "max_neutral_mass_range_da": 1.5,
        "max_mass_error_ppm": 20,
        "min_relative_intensity_percent": 0.5,
        "max_competing_envelopes": 3,
        "comparison_ready_statuses": ["Reliable", "Review"],
    }
}


def _candidate(cluster_id, charges, mass=1000.0, intensity=10000.0, mass_error_ppm=1.0):
    return IntactMassCandidate(
        observed_mass=mass,
        charge_state_count=len(charges),
        charge_states=charges,
        supporting_peak_count=len(charges),
        total_intensity=intensity,
        theoretical_mass=mass,
        mass_error_da=0.0,
        mass_error_ppm=mass_error_ppm,
        cluster_id=cluster_id,
    )


def _peaks(cluster_id, charges, masses=None, intensity=1000.0):
    masses = masses or [1000.0 for _ in charges]
    return [
        {
            "Cluster_ID": cluster_id,
            "Charge": charge,
            "Neutral_Mass": mass,
            "Intensity": intensity,
        }
        for charge, mass in zip(charges, masses)
    ]


def _qc(candidate, peaks):
    rows, diagnostics = build_intact_reconstruction_qc([candidate], peaks, BASE_QC_CONFIG, reconstruction_enabled=True)
    return rows[0], diagnostics[0], candidate


def test_contiguous_three_charge_states_are_reliable_and_comparison_ready():
    row, diagnostics, candidate = _qc(_candidate("C1", [10, 11, 12]), _peaks("C1", [10, 11, 12]))
    assert row["Reconstruction_Status"] == "Reliable"
    assert row["Reconstruction_Confidence"] == "High"
    assert row["Charge_State_Continuity"] == "contiguous"
    assert row["Comparison_Ready"] is True
    assert candidate.comparison_ready is True
    assert diagnostics["Reliable_Count"] == 1


def test_two_charge_states_are_review_and_comparison_ready():
    row, _, _ = _qc(_candidate("C2", [10, 11]), _peaks("C2", [10, 11]))
    assert row["Reconstruction_Status"] == "Review"
    assert row["Reconstruction_Confidence"] == "Medium"
    assert row["Comparison_Ready"] is True


def test_single_charge_state_is_insufficient_not_comparison_ready():
    row, _, _ = _qc(_candidate("C3", [10]), _peaks("C3", [10]))
    assert row["Reconstruction_Status"] == "Insufficient"
    assert row["Primary_Limiting_Factor"] == "insufficient_charge_states"
    assert row["Comparison_Ready"] is False


def test_non_contiguous_charge_series_warns():
    row, _, _ = _qc(_candidate("C4", [10, 12, 13]), _peaks("C4", [10, 12, 13]))
    assert row["Reconstruction_Status"] == "Review"
    assert row["Charge_State_Continuity"] == "non_contiguous"
    assert row["Primary_Limiting_Factor"] == "non_contiguous_charge_states"


def test_large_mass_spread_adds_limiting_reason():
    row, _, _ = _qc(_candidate("C5", [10, 11, 12]), _peaks("C5", [10, 11, 12], masses=[1000.0, 1001.2, 1002.4]))
    assert row["Reconstruction_Status"] == "Review"
    assert row["Neutral_Mass_Range"] > 1.5
    assert row["Primary_Limiting_Factor"] == "mass_spread_too_large"


def _excel_config(reconstruction):
    return SimpleNamespace(
        project={"name": "test"},
        input={},
        organism={},
        sequence={},
        experiment={},
        instrument={},
        reconstruction=reconstruction,
        digestion={"enabled": True},
        alkaline_phosphatase={},
        fragment_mapping={},
        modification_search={},
        peak_filtering={},
        p1_annotation={},
        ms2_annotation={},
        modification_evidence_ranking={},
        biological_context={},
        performance={},
        reporting={"max_excel_rows_per_sheet": 1000, "truncate_large_sheets": True},
    )


def _write_empty_excel(tmp_path, reconstruction):
    return write_excel_report(
        output_dir=tmp_path,
        config=_excel_config(reconstruction),
        diagnostics={},
        intact_results=[],
        charge_state_peaks=[],
        warnings=[],
        modifications=[],
        rule_set={},
        pathways=[],
        theoretical_fragments=[],
        fragment_ms1_matches=[],
        known_modification_candidates=[],
        known_modification_summary=[],
        optional_results={},
    )


def _diagnostic_row(workbook):
    diagnostics = workbook["Intact_Reconstruction_Diag"]
    headers = [cell.value for cell in next(diagnostics.iter_rows(min_row=3, max_row=3))]
    values = [cell.value for cell in next(diagnostics.iter_rows(min_row=4, max_row=4))]
    return dict(zip(headers, values))


def _assert_intact_qc_sheet_names(workbook):
    assert "Intact_Reconstruction_QC" in workbook.sheetnames
    assert "Intact_Reconstruction_Diag" in workbook.sheetnames
    old_sheet_name = "Intact_Reconstruction_" + "Diagnostics"
    assert old_sheet_name not in workbook.sheetnames
    assert all(len(name) <= 31 for name in workbook.sheetnames)


def test_zero_candidates_excel_output_succeeds_with_qc_sheets(tmp_path, recwarn):
    report = _write_empty_excel(tmp_path, {"enabled": True, **BASE_QC_CONFIG})
    workbook = load_workbook(report, read_only=True, data_only=True)
    try:
        assert "Intact_mass_reconstruction" in workbook.sheetnames
        _assert_intact_qc_sheet_names(workbook)
        row = _diagnostic_row(workbook)
        assert row["Total_Reconstruction_Candidates"] == 0
        assert "no_charge_state_candidates" in row["Failure_Reason_Counts"]
        assert not [warning for warning in recwarn if "Title is more than 31 characters" in str(warning.message)]
    finally:
        workbook.close()


def test_zero_candidates_excel_output_uses_python_defaults_without_qc_config(tmp_path):
    report = _write_empty_excel(tmp_path, {"enabled": True})
    workbook = load_workbook(report, read_only=True, data_only=True)
    try:
        _assert_intact_qc_sheet_names(workbook)
        row = _diagnostic_row(workbook)
        assert row["Min_Charge_States_For_Reliable"] == 3
        assert row["Min_Charge_States_For_Review"] == 2
        assert row["Comparison_Ready_Statuses"] == "Reliable; Review"
    finally:
        workbook.close()


def test_reconstruction_disabled_still_writes_diagnostic_sheet(tmp_path):
    report = _write_empty_excel(tmp_path, {"enabled": False, **BASE_QC_CONFIG})
    workbook = load_workbook(report, read_only=True, data_only=True)
    try:
        assert "Intact_mass_reconstruction" not in workbook.sheetnames
        _assert_intact_qc_sheet_names(workbook)
        row = _diagnostic_row(workbook)
        assert row["Reconstruction_Enabled"] is False
        assert "reconstruction_disabled" in row["Failure_Reason_Counts"]
    finally:
        workbook.close()
