from types import SimpleNamespace

from openpyxl import load_workbook

from rna_masshunter.excel_report import write_excel_report
from rna_masshunter.intact_reconstruction import build_intact_reconstruction_qc
from rna_masshunter.models import IntactMassCandidate


UNMODIFIED_MASS = 25082.3
MODIFIED_MASS = 25325.5
ALT_MODIFIED_MASS = 25342.0
LOW_FULL_LENGTH_MASS = 25082.0
FRAGMENT_LIKE_MASS = 4818.0

BASE_QC_CONFIG = {
    "intact_reconstruction": {
        "min_charge_states_for_reliable": 3,
        "min_charge_states_for_review": 2,
        "require_contiguous_charge_states": True,
        "max_neutral_mass_sd_da": 0.5,
        "max_neutral_mass_range_da": 1.5,
        "max_mass_error_ppm": 20,
        "max_envelope_internal_error_ppm": 20,
        "min_relative_envelope_intensity_percent_for_reliable": 1.0,
        "min_relative_envelope_intensity_percent_for_review": 0.1,
        "max_competing_envelopes": 3,
        "comparison_ready_statuses": ["Reliable", "Review"],
        "max_rt_range_min_for_reliable": 0.15,
        "max_rt_range_min_for_review": 0.30,
        "allow_trace_only_reliable": False,
        "search_mode": "untargeted",
        "reference_mass_tolerance_ppm": 20,
    }
}


def _candidate(cluster_id, charges, mass=1000.0, intensity=10000.0, theoretical_mass=1000.0):
    delta = mass - theoretical_mass if theoretical_mass is not None else None
    ppm = delta / theoretical_mass * 1_000_000 if theoretical_mass else None
    return IntactMassCandidate(
        observed_mass=mass,
        charge_state_count=len(charges),
        charge_states=charges,
        supporting_peak_count=len(charges),
        total_intensity=intensity,
        theoretical_mass=theoretical_mass,
        mass_error_da=delta,
        mass_error_ppm=ppm,
        cluster_id=cluster_id,
    )


def _peaks(cluster_id, charges, masses=None, rts=None, intensity=1000.0, tier="Major"):
    masses = masses or [1000.0 for _ in charges]
    rts = rts if rts is not None else [5.0 + index * 0.02 for index, _ in enumerate(charges)]
    return [
        {
            "Cluster_ID": cluster_id,
            "Charge": charge,
            "Neutral_Mass": mass,
            "Intensity": intensity,
            "RT": rt,
            "Peak_Tier": tier,
        }
        for charge, mass, rt in zip(charges, masses, rts)
    ]


def _qc(candidates, peaks, config=None, enabled=True):
    rows, diagnostics = build_intact_reconstruction_qc(candidates, peaks, config or BASE_QC_CONFIG, reconstruction_enabled=enabled)
    return rows, diagnostics[0]


def test_modified_envelope_far_from_unmodified_theory_can_be_reliable():
    candidate = _candidate("C1", [10, 11, 12], mass=MODIFIED_MASS, intensity=50000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("C1", [10, 11, 12], masses=[25325.45, 25325.50, 25325.55], intensity=10000)
    rows, _ = _qc([candidate], peaks)
    row = rows[0]
    assert row["Reconstruction_Status"] == "Reliable"
    assert row["Comparison_Ready_Strict"] is True
    assert row["Unmodified_Theory_Delta_Da"] == candidate.mass_error_da
    assert round(row["Unmodified_Theory_Delta_Da"], 1) == 243.2
    assert abs(row["Unmodified_Theory_Delta_ppm"]) > 9000


def test_large_unmodified_delta_alone_does_not_lower_reliable_status():
    candidate = _candidate("C2", [20, 21, 22], mass=MODIFIED_MASS, intensity=40000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("C2", [20, 21, 22], masses=[25325.49, 25325.50, 25325.51], intensity=9000)
    rows, _ = _qc([candidate], peaks)
    assert rows[0]["Reconstruction_Status"] == "Reliable"
    assert "mass_error_too_large" not in rows[0]["Limiting_Factors"]


def test_envelope_internal_mass_error_prevents_reliable():
    candidate = _candidate("C3", [10, 11, 12], mass=1000.0, intensity=30000)
    peaks = _peaks("C3", [10, 11, 12], masses=[1000.0, 1001.0, 1002.0], intensity=8000)
    rows, _ = _qc([candidate], peaks)
    assert rows[0]["Reconstruction_Status"] == "Review"
    assert "internal_mass_error_too_large" in rows[0]["Limiting_Factors"]
    assert rows[0]["Comparison_Ready"] is False


def test_same_rt_three_contiguous_charge_states_are_reliable():
    candidate = _candidate("C4", [7, 8, 9], mass=1500.0, intensity=25000)
    peaks = _peaks("C4", [7, 8, 9], masses=[1499.99, 1500.0, 1500.01], rts=[3.00, 3.04, 3.08], intensity=7000)
    rows, _ = _qc([candidate], peaks)
    assert rows[0]["Charge_State_Continuity"] == "contiguous"
    assert rows[0]["RT_Consistency"] == "consistent"
    assert rows[0]["Reconstruction_Status"] == "Reliable"


def test_rt_range_too_large_adds_rt_inconsistent():
    candidate = _candidate("C5", [10, 11, 12], mass=1000.0, intensity=30000)
    peaks = _peaks("C5", [10, 11, 12], masses=[999.99, 1000.0, 1000.01], rts=[1.0, 1.3, 1.7], intensity=9000)
    rows, _ = _qc([candidate], peaks)
    assert rows[0]["RT_Consistency"] == "inconsistent"
    assert "rt_inconsistent" in rows[0]["Limiting_Factors"]
    assert rows[0]["Reconstruction_Status"] == "Review"


def test_trace_only_envelope_is_not_reliable_by_default():
    candidate = _candidate("C6", [10, 11, 12], mass=1000.0, intensity=30000)
    peaks = _peaks("C6", [10, 11, 12], masses=[999.99, 1000.0, 1000.01], intensity=9000, tier="Trace")
    rows, _ = _qc([candidate], peaks)
    assert rows[0]["Trace_Only_Envelope"] is True
    assert rows[0]["Reconstruction_Status"] == "Review"
    assert "trace_only_envelope" in rows[0]["Limiting_Factors"]


def test_highest_intensity_envelope_is_dominant():
    low = _candidate("LOW", [10, 11, 12], mass=1000.0, intensity=10000)
    high = _candidate("HIGH", [13, 14, 15], mass=2000.0, intensity=90000)
    peaks = _peaks("LOW", [10, 11, 12], masses=[999.99, 1000.0, 1000.01], intensity=3000) + _peaks(
        "HIGH", [13, 14, 15], masses=[1999.99, 2000.0, 2000.01], intensity=30000
    )
    _, diag = _qc([low, high], peaks)
    assert diag["Dominant_Envelope_Mass"] == 2000.0
    assert diag["Dominant_Envelope_Intensity"] == 90000


def test_comparison_ready_review_allows_minor_noncontiguous_case():
    candidate = _candidate("C7", [10, 12], mass=MODIFIED_MASS, intensity=30000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("C7", [10, 12], masses=[25325.49, 25325.51], rts=[2.0, 2.1], intensity=10000)
    rows, _ = _qc([candidate], peaks)
    assert rows[0]["Reconstruction_Status"] == "Review"
    assert rows[0]["Comparison_Ready_Strict"] is False
    assert rows[0]["Comparison_Ready_Review"] is True
    assert rows[0]["Comparison_Ready"] is True
    assert rows[0]["In_Neutral_Mass_Search_Range"] is True


def test_multiple_limiting_factors_are_preserved():
    candidate = _candidate("C8", [10, 12, 14], mass=1000.0, intensity=30000)
    peaks = _peaks("C8", [10, 12, 14], masses=[1000.0, 1002.0, 1004.0], rts=[1.0, 1.4, 1.8], intensity=8000, tier="Trace")
    rows, _ = _qc([candidate], peaks)
    factors = rows[0]["Limiting_Factors"]
    assert "non_contiguous_charge_states" in factors
    assert "mass_spread_too_large" in factors
    assert "rt_inconsistent" in factors
    assert "trace_only_envelope" in factors
    assert rows[0]["Num_Limiting_Factors"] >= 4


def test_reference_mass_match_is_calculated():
    config = {"intact_reconstruction": {**BASE_QC_CONFIG["intact_reconstruction"], "reference_masses": [{"label": "SCiex_major_25325", "mass_da": 25325.5}]}}
    candidate = _candidate("C9", [10, 11, 12], mass=25325.52, intensity=40000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("C9", [10, 11, 12], masses=[25325.51, 25325.52, 25325.53], intensity=10000)
    rows, diag = _qc([candidate], peaks, config=config)
    assert rows[0]["Best_Reference_Label"] == "SCiex_major_25325"
    assert abs(rows[0]["Reference_Mass_Error_Da"] - 0.02) < 1e-6
    assert rows[0]["Reference_Mass_Matched"] is True
    assert diag["Reference_Match_Count"] == 1


def test_reference_not_configured_is_safe():
    candidate = _candidate("C10", [10, 11, 12], mass=1000.0, intensity=30000)
    peaks = _peaks("C10", [10, 11, 12], masses=[999.99, 1000.0, 1000.01], intensity=9000)
    rows, diag = _qc([candidate], peaks, config={"intact_reconstruction": {}})
    assert rows[0]["Best_Reference_Label"] == "not_configured"
    assert rows[0]["Reference_Mass_Matched"] is False
    assert diag["Reference_Masses_Used"] == "not_configured"


def test_python_defaults_work_without_new_qc_settings():
    candidate = _candidate("C11", [10, 11, 12], mass=MODIFIED_MASS, intensity=30000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("C11", [10, 11, 12], masses=[25325.49, 25325.5, 25325.51], intensity=9000)
    rows, diag = _qc([candidate], peaks, config={"intact_reconstruction": {}})
    assert rows[0]["Reconstruction_Status"] == "Reliable"
    assert rows[0]["In_Neutral_Mass_Search_Range"] is True
    assert rows[0]["Neutral_Mass_Search_Min_Da"] == 20000.0
    assert rows[0]["Neutral_Mass_Search_Max_Da"] == 30000.0
    assert diag["Search_Mode"] == "untargeted"
    assert diag["Neutral_Mass_Search_Min_Da"] == 20000.0
    assert diag["Neutral_Mass_Search_Max_Da"] == 30000.0
    assert diag["Min_Relative_Envelope_Intensity_Percent_For_Reliable"] == 1.0


def test_fragment_like_4818_da_is_outside_default_neutral_mass_range():
    candidate = _candidate("FRAG", [10, 11, 12], mass=FRAGMENT_LIKE_MASS, intensity=50000)
    peaks = _peaks("FRAG", [10, 11, 12], masses=[4817.99, 4818.0, 4818.01], intensity=10000)
    rows, diag = _qc([candidate], peaks)
    assert rows[0]["In_Neutral_Mass_Search_Range"] is False
    assert rows[0]["Neutral_Mass_Range_Status"] == "outside_range"
    assert rows[0]["Comparison_Ready"] is False
    assert "outside_neutral_mass_search_range" in rows[0]["Limiting_Factors"]
    assert diag["Total_Candidates_Outside_Mass_Range"] == 1
    assert diag["Total_Candidates_In_Mass_Range"] == 0


def test_full_length_masses_are_inside_default_neutral_mass_range():
    candidates = [
        _candidate("M25082", [10, 11, 12], mass=LOW_FULL_LENGTH_MASS, intensity=30000, theoretical_mass=UNMODIFIED_MASS),
        _candidate("M25325", [10, 11, 12], mass=MODIFIED_MASS, intensity=40000, theoretical_mass=UNMODIFIED_MASS),
        _candidate("M25342", [10, 11, 12], mass=ALT_MODIFIED_MASS, intensity=35000, theoretical_mass=UNMODIFIED_MASS),
    ]
    peaks = (
        _peaks("M25082", [10, 11, 12], masses=[25081.99, 25082.0, 25082.01], intensity=7000)
        + _peaks("M25325", [10, 11, 12], masses=[25325.49, 25325.5, 25325.51], intensity=8000)
        + _peaks("M25342", [10, 11, 12], masses=[25341.99, 25342.0, 25342.01], intensity=7500)
    )
    rows, diag = _qc(candidates, peaks)
    assert {row["Reconstructed_Mass"]: row["In_Neutral_Mass_Search_Range"] for row in rows} == {
        LOW_FULL_LENGTH_MASS: True,
        MODIFIED_MASS: True,
        ALT_MODIFIED_MASS: True,
    }
    assert diag["Total_Candidates_In_Mass_Range"] == 3
    assert diag["Total_Candidates_Outside_Mass_Range"] == 0


def test_neutral_mass_range_can_be_changed_and_changes_readiness():
    config = {
        "intact_reconstruction": {
            **BASE_QC_CONFIG["intact_reconstruction"],
            "neutral_mass_range": {"enabled": True, "min_da": 24000, "max_da": 26000},
        }
    }
    in_candidate = _candidate("IN", [10, 11, 12], mass=MODIFIED_MASS, intensity=30000, theoretical_mass=UNMODIFIED_MASS)
    out_candidate = _candidate("OUT", [10, 11, 12], mass=23000.0, intensity=30000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("IN", [10, 11, 12], masses=[25325.49, 25325.5, 25325.51], intensity=9000) + _peaks(
        "OUT", [10, 11, 12], masses=[22999.99, 23000.0, 23000.01], intensity=9000
    )
    rows, diag = _qc([in_candidate, out_candidate], peaks, config=config)
    by_cluster = {row["Cluster_ID"]: row for row in rows}
    assert by_cluster["IN"]["In_Neutral_Mass_Search_Range"] is True
    assert by_cluster["IN"]["Comparison_Ready"] is True
    assert by_cluster["OUT"]["In_Neutral_Mass_Search_Range"] is False
    assert by_cluster["OUT"]["Comparison_Ready"] is False
    assert diag["Neutral_Mass_Search_Min_Da"] == 24000.0
    assert diag["Neutral_Mass_Search_Max_Da"] == 26000.0
    assert diag["Total_Candidates_In_Mass_Range"] == 1
    assert diag["Total_Candidates_Outside_Mass_Range"] == 1


def test_neutral_mass_range_change_can_exclude_previous_default_candidate():
    default_candidate = _candidate("DEFAULT", [10, 11, 12], mass=MODIFIED_MASS, intensity=30000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("DEFAULT", [10, 11, 12], masses=[25325.49, 25325.5, 25325.51], intensity=9000)
    default_rows, _ = _qc([default_candidate], peaks)
    assert default_rows[0]["Comparison_Ready"] is True

    narrowed = {
        "intact_reconstruction": {
            **BASE_QC_CONFIG["intact_reconstruction"],
            "neutral_mass_range": {"enabled": True, "min_da": 26000, "max_da": 30000},
        }
    }
    narrowed_rows, narrowed_diag = _qc([default_candidate], peaks, config=narrowed)
    assert narrowed_rows[0]["In_Neutral_Mass_Search_Range"] is False
    assert narrowed_rows[0]["Comparison_Ready"] is False
    assert narrowed_diag["Total_Candidates_In_Mass_Range"] == 0


def test_dominant_envelope_in_search_range_ignores_out_of_range_overall_dominant():
    fragment = _candidate("FRAG_DOM", [10, 11, 12], mass=FRAGMENT_LIKE_MASS, intensity=90000)
    in_range_low = _candidate("IN_LOW", [10, 11, 12], mass=LOW_FULL_LENGTH_MASS, intensity=30000, theoretical_mass=UNMODIFIED_MASS)
    in_range_high = _candidate("IN_HIGH", [10, 11, 12], mass=MODIFIED_MASS, intensity=50000, theoretical_mass=UNMODIFIED_MASS)
    peaks = (
        _peaks("FRAG_DOM", [10, 11, 12], masses=[4817.99, 4818.0, 4818.01], intensity=30000)
        + _peaks("IN_LOW", [10, 11, 12], masses=[25081.99, 25082.0, 25082.01], intensity=9000)
        + _peaks("IN_HIGH", [10, 11, 12], masses=[25325.49, 25325.5, 25325.51], intensity=16000)
    )
    _, diag = _qc([fragment, in_range_low, in_range_high], peaks)
    assert diag["Dominant_Envelope_Overall_Mass"] == FRAGMENT_LIKE_MASS
    assert diag["Dominant_Envelope_Overall_Intensity"] == 90000
    assert diag["Dominant_Envelope_In_Mass_Range_Mass"] == MODIFIED_MASS
    assert diag["Dominant_Envelope_In_Mass_Range_Intensity"] == 50000


def test_bad_high_intensity_in_range_candidate_is_raw_dominant_not_intact_dominant():
    bad = _candidate("BAD_RAW", [10, 11, 12], mass=22286.0, intensity=90000, theoretical_mass=UNMODIFIED_MASS)
    good = _candidate("GOOD", [10, 11, 12], mass=MODIFIED_MASS, intensity=30000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("BAD_RAW", [10, 11, 12], masses=[22284.0, 22286.0, 22288.0], rts=[1.0, 1.4, 1.8], intensity=30000) + _peaks(
        "GOOD", [10, 11, 12], masses=[25325.49, 25325.50, 25325.51], intensity=10000
    )
    rows, diag = _qc([bad, good], peaks)
    by_cluster = {row["Cluster_ID"]: row for row in rows}
    assert diag["Dominant_Envelope_In_Search_Range_Raw_Mass"] == 22286.0
    assert by_cluster["BAD_RAW"]["Envelope_QC_Eligible"] is False
    assert by_cluster["BAD_RAW"]["Comparison_Ready"] is False
    assert diag["Dominant_Intact_Eligible_Envelope_Mass"] == MODIFIED_MASS
    assert by_cluster["GOOD"]["Dominant_Intact_Envelope_Flag"] is True


def test_good_qc_candidate_can_be_intact_eligible_despite_lower_intensity():
    bad = _candidate("BAD", [10, 11, 12], mass=22286.0, intensity=90000, theoretical_mass=UNMODIFIED_MASS)
    good = _candidate("LOW_GOOD", [10, 11, 12], mass=LOW_FULL_LENGTH_MASS, intensity=15000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("BAD", [10, 11, 12], masses=[22284.0, 22286.0, 22288.0], rts=[1.0, 1.4, 1.8], intensity=30000) + _peaks(
        "LOW_GOOD", [10, 11, 12], masses=[25081.99, 25082.0, 25082.01], intensity=5000
    )
    rows, _ = _qc([bad, good], peaks)
    by_cluster = {row["Cluster_ID"]: row for row in rows}
    assert by_cluster["BAD"]["Intact_Strict_Eligible"] is False
    assert by_cluster["LOW_GOOD"]["Intact_Strict_Eligible"] is True
    assert by_cluster["LOW_GOOD"]["Comparison_Ready_Strict"] is True


def test_strict_candidate_is_preferred_over_review_candidate():
    strict = _candidate("STRICT", [10, 11, 12], mass=MODIFIED_MASS, intensity=20000, theoretical_mass=UNMODIFIED_MASS)
    review = _candidate("REVIEW", [10, 12], mass=LOW_FULL_LENGTH_MASS, intensity=80000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("STRICT", [10, 11, 12], masses=[25325.49, 25325.50, 25325.51], intensity=7000) + _peaks(
        "REVIEW", [10, 12], masses=[25081.99, 25082.01], rts=[3.0, 3.1], intensity=40000
    )
    rows, diag = _qc([strict, review], peaks)
    by_cluster = {row["Cluster_ID"]: row for row in rows}
    assert by_cluster["STRICT"]["Intact_Strict_Eligible"] is True
    assert by_cluster["REVIEW"]["Intact_Review_Eligible"] is True
    assert diag["Dominant_Intact_Eligible_Envelope_Mass"] == MODIFIED_MASS
    assert by_cluster["STRICT"]["Dominant_Intact_Envelope_Flag"] is True


def test_review_top_is_dominant_when_no_strict_candidate_exists():
    better_review = _candidate("REVIEW_A", [10, 12, 13], mass=MODIFIED_MASS, intensity=25000, theoretical_mass=UNMODIFIED_MASS)
    lower_review = _candidate("REVIEW_B", [10, 12], mass=LOW_FULL_LENGTH_MASS, intensity=50000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("REVIEW_A", [10, 12, 13], masses=[25325.49, 25325.50, 25325.51], rts=[2.0, 2.1, 2.2], intensity=8000) + _peaks(
        "REVIEW_B", [10, 12], masses=[25081.99, 25082.01], rts=[3.0, 3.1], intensity=25000
    )
    rows, diag = _qc([better_review, lower_review], peaks)
    by_cluster = {row["Cluster_ID"]: row for row in rows}
    assert diag["Intact_Strict_Eligible_Count"] == 0
    assert by_cluster["REVIEW_A"]["Intact_Review_Eligible"] is True
    assert by_cluster["REVIEW_B"]["Intact_Review_Eligible"] is True
    assert diag["Dominant_Intact_Review_Envelope_Mass"] == MODIFIED_MASS
    assert diag["Dominant_Intact_Eligible_Envelope_Mass"] == MODIFIED_MASS


def test_trace_only_alone_does_not_remove_review_eligibility():
    candidate = _candidate("TRACE_REVIEW", [10, 12], mass=MODIFIED_MASS, intensity=30000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("TRACE_REVIEW", [10, 12], masses=[25325.49, 25325.51], rts=[2.0, 2.1], intensity=10000, tier="Trace")
    rows, _ = _qc([candidate], peaks)
    row = rows[0]
    assert row["Trace_Only_Envelope"] is True
    assert row["Intact_Review_Eligible"] is True
    assert row["Comparison_Ready_Review"] is True
    assert "trace_only_envelope" in row["Comparison_Readiness_Reason"] or row["Comparison_Readiness_Reason"] == "review"


def test_relative_intact_eligible_intensity_is_normalized_within_eligible_candidates():
    bad = _candidate("BAD", [10, 11, 12], mass=22286.0, intensity=90000, theoretical_mass=UNMODIFIED_MASS)
    eligible_max = _candidate("ELIGIBLE_MAX", [10, 11, 12], mass=MODIFIED_MASS, intensity=30000, theoretical_mass=UNMODIFIED_MASS)
    eligible_half = _candidate("ELIGIBLE_HALF", [10, 12], mass=LOW_FULL_LENGTH_MASS, intensity=15000, theoretical_mass=UNMODIFIED_MASS)
    peaks = (
        _peaks("BAD", [10, 11, 12], masses=[22284.0, 22286.0, 22288.0], rts=[1.0, 1.4, 1.8], intensity=30000)
        + _peaks("ELIGIBLE_MAX", [10, 11, 12], masses=[25325.49, 25325.5, 25325.51], intensity=10000)
        + _peaks("ELIGIBLE_HALF", [10, 12], masses=[25081.99, 25082.01], intensity=7500)
    )
    rows, _ = _qc([bad, eligible_max, eligible_half], peaks)
    by_cluster = {row["Cluster_ID"]: row for row in rows}
    assert by_cluster["BAD"]["Relative_Intact_Eligible_Intensity_Percent"] == 0.0
    assert by_cluster["ELIGIBLE_MAX"]["Relative_Intact_Eligible_Intensity_Percent"] == 100.0
    assert by_cluster["ELIGIBLE_HALF"]["Relative_Intact_Eligible_Intensity_Percent"] == 50.0


def test_target_review_range_enabled_disabled_and_counts():
    candidates = [
        _candidate("M22286", [10, 11, 12], mass=22286.0, intensity=10000, theoretical_mass=UNMODIFIED_MASS),
        _candidate("M25082", [10, 11, 12], mass=LOW_FULL_LENGTH_MASS, intensity=30000, theoretical_mass=UNMODIFIED_MASS),
        _candidate("M25325", [10, 11, 12], mass=MODIFIED_MASS, intensity=40000, theoretical_mass=UNMODIFIED_MASS),
        _candidate("M25343", [10, 11, 12], mass=25343.0, intensity=35000, theoretical_mass=UNMODIFIED_MASS),
    ]
    peaks = (
        _peaks("M22286", [10, 11, 12], masses=[22285.99, 22286.0, 22286.01], intensity=3000)
        + _peaks("M25082", [10, 11, 12], masses=[25081.99, 25082.0, 25082.01], intensity=7000)
        + _peaks("M25325", [10, 11, 12], masses=[25325.49, 25325.5, 25325.51], intensity=8000)
        + _peaks("M25343", [10, 11, 12], masses=[25342.99, 25343.0, 25343.01], intensity=7500)
    )
    default_rows, default_diag = _qc(candidates, peaks)
    assert {row["Target_Review_Mass_Range_Status"] for row in default_rows} == {"not_configured"}
    assert default_diag["Target_Review_Mass_Range_Settings"] == "disabled"

    config = {
        "intact_reconstruction": {
            **BASE_QC_CONFIG["intact_reconstruction"],
            "target_review_mass_range": {"enabled": True, "min_da": 24000, "max_da": 26000},
        }
    }
    rows, diag = _qc(candidates, peaks, config=config)
    by_cluster = {row["Cluster_ID"]: row for row in rows}
    assert by_cluster["M22286"]["In_Target_Review_Mass_Range"] is False
    assert by_cluster["M22286"]["Target_Review_Mass_Range_Status"] == "outside_range"
    assert by_cluster["M25082"]["In_Target_Review_Mass_Range"] is True
    assert by_cluster["M25325"]["In_Target_Review_Mass_Range"] is True
    assert by_cluster["M25343"]["In_Target_Review_Mass_Range"] is True
    assert diag["Target_Review_Candidate_Count"] == 3


def test_sciex_reference_mass_error_for_25325_581():
    config = {
        "intact_reconstruction": {
            **BASE_QC_CONFIG["intact_reconstruction"],
            "reference_masses": [
                {"label": "SCiex_25325_5", "mass_da": 25325.5},
                {"label": "SCiex_25342", "mass_da": 25342.0},
                {"label": "SCiex_25343", "mass_da": 25343.0},
            ],
            "reference_mass_tolerance_ppm": 20,
        }
    }
    candidate = _candidate("REF", [10, 11, 12], mass=25325.581, intensity=40000, theoretical_mass=UNMODIFIED_MASS)
    peaks = _peaks("REF", [10, 11, 12], masses=[25325.571, 25325.581, 25325.591], intensity=10000)
    rows, _ = _qc([candidate], peaks, config=config)
    row = rows[0]
    assert row["Reference_Mass_Matched"] is True
    assert row["Best_Reference_Label"] == "SCiex_25325_5"
    assert round(row["Reference_Mass_Error_Da"], 3) == 0.081
    assert round(row["Reference_Mass_Error_ppm"], 1) == 3.2
    assert row["Reconstruction_Status"] == "Reliable"


def test_generic_logic_works_without_reference_masses_or_trna_specific_values():
    config = {
        "intact_reconstruction": {
            **BASE_QC_CONFIG["intact_reconstruction"],
            "neutral_mass_range": {"enabled": True, "min_da": 15000, "max_da": 16000},
        }
    }
    good = _candidate("GENERIC_GOOD", [6, 7, 8], mass=15555.5, intensity=42000, theoretical_mass=12345.6)
    out = _candidate("GENERIC_OUT", [6, 7, 8], mass=17000.0, intensity=84000, theoretical_mass=12345.6)
    peaks = _peaks("GENERIC_GOOD", [6, 7, 8], masses=[15555.49, 15555.50, 15555.51], intensity=14000) + _peaks(
        "GENERIC_OUT", [6, 7, 8], masses=[16999.99, 17000.00, 17000.01], intensity=28000
    )
    rows, diag = _qc([good, out], peaks, config=config)
    by_cluster = {row["Cluster_ID"]: row for row in rows}
    assert by_cluster["GENERIC_GOOD"]["Best_Reference_Label"] == "not_configured"
    assert by_cluster["GENERIC_GOOD"]["Intact_Strict_Eligible"] is True
    assert by_cluster["GENERIC_GOOD"]["Comparison_Ready"] is True
    assert by_cluster["GENERIC_OUT"]["In_Neutral_Mass_Search_Range"] is False
    assert by_cluster["GENERIC_OUT"]["Comparison_Ready"] is False
    assert diag["Dominant_Intact_Eligible_Envelope_Mass"] == 15555.5


def test_reference_annotations_do_not_change_qc_score_rank_or_dominant_selection():
    candidates_without_ref = [
        _candidate("A", [6, 7, 8], mass=15555.5, intensity=40000, theoretical_mass=12345.6),
        _candidate("B", [6, 8], mass=15620.0, intensity=50000, theoretical_mass=12345.6),
    ]
    peaks = _peaks("A", [6, 7, 8], masses=[15555.49, 15555.50, 15555.51], intensity=13000) + _peaks(
        "B", [6, 8], masses=[15619.99, 15620.01], rts=[3.0, 3.1], intensity=25000
    )
    base_config = {
        "intact_reconstruction": {
            **BASE_QC_CONFIG["intact_reconstruction"],
            "neutral_mass_range": {"enabled": True, "min_da": 15000, "max_da": 16000},
        }
    }
    ref_config = {
        "intact_reconstruction": {
            **base_config["intact_reconstruction"],
            "reference_masses": [{"label": "external_check_A", "mass_da": 15555.5}],
            "reference_mass_tolerance_ppm": 20,
        }
    }
    rows_no_ref, diag_no_ref = _qc(candidates_without_ref, peaks, config=base_config)
    candidates_with_ref = [
        _candidate("A", [6, 7, 8], mass=15555.5, intensity=40000, theoretical_mass=12345.6),
        _candidate("B", [6, 8], mass=15620.0, intensity=50000, theoretical_mass=12345.6),
    ]
    rows_with_ref, diag_with_ref = _qc(candidates_with_ref, peaks, config=ref_config)
    no_ref = {row["Cluster_ID"]: row for row in rows_no_ref}
    with_ref = {row["Cluster_ID"]: row for row in rows_with_ref}
    for cluster_id in ["A", "B"]:
        for field in [
            "Reconstruction_Status",
            "Envelope_QC_Eligible",
            "Intact_Strict_Eligible",
            "Intact_Review_Eligible",
            "Comparison_Ready",
            "Intact_Envelope_QC_Score",
            "Intact_Envelope_QC_Rank",
            "Dominant_Intact_Envelope_Flag",
        ]:
            assert with_ref[cluster_id][field] == no_ref[cluster_id][field]
    assert diag_with_ref["Dominant_Intact_Eligible_Envelope_Mass"] == diag_no_ref["Dominant_Intact_Eligible_Envelope_Mass"]
    assert with_ref["A"]["Reference_Mass_Matched"] is True
    assert no_ref["A"]["Reference_Mass_Matched"] is False


def test_target_review_range_enabled_without_bounds_is_not_configured_and_non_filtering():
    config = {
        "intact_reconstruction": {
            **BASE_QC_CONFIG["intact_reconstruction"],
            "neutral_mass_range": {"enabled": True, "min_da": 15000, "max_da": 16000},
            "target_review_mass_range": {"enabled": True},
        }
    }
    candidate = _candidate("NO_BOUNDS", [6, 7, 8], mass=15555.5, intensity=40000, theoretical_mass=12345.6)
    peaks = _peaks("NO_BOUNDS", [6, 7, 8], masses=[15555.49, 15555.50, 15555.51], intensity=13000)
    rows, diag = _qc([candidate], peaks, config=config)
    assert rows[0]["Target_Review_Mass_Range_Status"] == "not_configured"
    assert rows[0]["Comparison_Ready"] is True
    assert diag["Target_Review_Mass_Range_Settings"] == "not_configured"
    assert diag["Target_Review_Candidate_Count"] == 0


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


def test_existing_excel_sheets_and_new_columns_are_present(tmp_path):
    report = _write_empty_excel(tmp_path, {"enabled": True, **BASE_QC_CONFIG})
    workbook = load_workbook(report, read_only=True, data_only=True)
    try:
        for sheet in ["Run_summary", "Input_parameters", "mzML_diagnostics", "Intact_mass_reconstruction", "Charge_state_peaks"]:
            assert sheet in workbook.sheetnames
        qc_headers = [cell.value for cell in next(workbook["Intact_Reconstruction_QC"].iter_rows(min_row=3, max_row=3))]
        diag_headers = [cell.value for cell in next(workbook["Intact_Reconstruction_Diag"].iter_rows(min_row=3, max_row=3))]
        assert "Comparison_Ready_Strict" in qc_headers
        assert "Comparison_Ready_Review" in qc_headers
        assert "In_Neutral_Mass_Search_Range" in qc_headers
        assert "Neutral_Mass_Search_Min_Da" in qc_headers
        assert "Neutral_Mass_Search_Max_Da" in qc_headers
        assert "Neutral_Mass_Range_Status" in qc_headers
        assert "Limiting_Factors" in qc_headers
        assert "Envelope_QC_Eligible" in qc_headers
        assert "Intact_Review_Eligible" in qc_headers
        assert "Intact_Strict_Eligible" in qc_headers
        assert "Intact_Envelope_QC_Score" in qc_headers
        assert "Relative_Intact_Eligible_Intensity_Percent" in qc_headers
        assert "In_Target_Review_Mass_Range" in qc_headers
        assert "Dominant_Envelope_Mass" in diag_headers
        assert "Dominant_Envelope_Overall_Mass" in diag_headers
        assert "Dominant_Envelope_In_Mass_Range_Mass" in diag_headers
        assert "Dominant_Envelope_In_Search_Range_Raw_Mass" in diag_headers
        assert "Dominant_Intact_Eligible_Envelope_Mass" in diag_headers
        assert "Envelope_QC_Eligible_Count" in diag_headers
        assert "Intact_Strict_Eligible_Count" in diag_headers
        assert "Intact_Review_Eligible_Count" in diag_headers
        assert "Target_Review_Candidate_Count" in diag_headers
        assert "Total_Candidates_In_Mass_Range" in diag_headers
        assert "Total_Candidates_Outside_Mass_Range" in diag_headers
    finally:
        workbook.close()
