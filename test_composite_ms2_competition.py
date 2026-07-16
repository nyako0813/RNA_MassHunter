from copy import deepcopy

import pytest

from rna_masshunter.audit_policy import AuditPolicy, included_sheet_names
from rna_masshunter.composite_ms2_matcher import match_composite_ms2
from rna_masshunter.models import MS2SpectrumInfo, RunConfig
from rna_masshunter.ms2_identity_evidence import physical_observed_peak_key


def config():
    return RunConfig(
        instrument={"polarity": "negative"},
        ms2_annotation={
            "mz_tolerance_ppm": 20,
            "precursor_match_tolerance_ppm": 20,
            "constrain_by_precursor": False,
        },
    )


def ion(candidate="C1", structure="S1", ion_id="I1", mz=500.0,
        positions="1", bonds=""):
    return {
        "Candidate_ID": candidate, "Complete_Structure_ID": structure,
        "Parent_Fragment_ID": "F1", "Parent_Neutral_Mass": 1000.0,
        "Ion_ID": ion_id, "Ion_Series": "c", "Ion_Number": 2,
        "Cleavage_Position": 2, "Included_Positions": positions,
        "Included_Modified_Positions": positions,
        "Included_Backbone_Bonds": bonds, "Theoretical_Neutral_Mass": 501.0,
        "Charge": 1, "Theoretical_mz": mz,
        "Position_Informative": bool(positions),
        "Backbone_Informative": bool(bonds),
    }


def spectrum(mz=500.001, intensity=100.0):
    return MS2SpectrumInfo(
        "SP1", 7, 1.25, 1000.0, 1, 1000.0, 1, mz, intensity,
        intensity, peaks=[(mz, intensity)],
    )


def run(ions, observed=500.001, intensity=100.0):
    return match_composite_ms2(
        ions, [spectrum(observed, intensity)], config(),
        return_competition=True,
    )


def test_one_physical_peak_one_assignment_and_standard_key():
    best, detail = run([ion()])
    assert len(best) == len(detail) == 1
    row = detail[0]
    expected = physical_observed_peak_key({
        "Spectrum_ID": "SP1", "Observed_mz": 500.001,
        "Observed_Intensity": 100.0, "RT": 1.25,
    })
    assert row["Physical_Observed_Peak_Key"] == expected
    assert row["Observed_Peak_Index"] == 0
    assert row["Raw_Peak_Index"] == ""
    assert row["Raw_Peak_Index_Missing_Reason"]
    assert row["Assignment_Rank"] == 1 and row["Best_Assignment"] is True
    assert row["Within_Tolerance_Assignment_Count"] == 1
    assert row["Candidate_Specific"] is True
    assert row["Complete_Structure_Specific"] is True
    assert row["Theoretical_Ion_Specific"] is True


def test_multiple_theoretical_ions_and_nonbest_retained():
    best, detail = run([
        ion(ion_id="I2", mz=500.002),
        ion(ion_id="I1", mz=500.0005),
    ])
    assert len(best) == 1 and len(detail) == 2
    assert best[0]["Theoretical_mz"] == 500.0005
    assert [row["Assignment_Rank"] for row in detail] == [1, 2]
    assert [row["Best_Assignment"] for row in detail] == [True, False]
    assert all(row["Within_Tolerance_Assignment_Count"] == 2 for row in detail)
    assert all(row["Competing_Theoretical_Ion_Count"] == 1 for row in detail)
    assert detail[0]["Competing_Ion_IDs"] == "I2"
    assert detail[1]["Competing_Ion_IDs"] == "I1"


def test_multiple_candidates_and_complete_structures():
    _, detail = run([
        ion(candidate="C2", structure="S2", ion_id="I2", mz=500.002),
        ion(candidate="C1", structure="S1", ion_id="I1", mz=500.0005),
    ])
    first = detail[0]
    assert first["Competing_Candidate_Count"] == 1
    assert first["Competing_Candidate_IDs"] == "C2"
    assert first["Competing_Complete_Structure_Count"] == 1
    assert first["Competing_Complete_Structure_IDs"] == "S2"
    assert first["Candidate_Specific"] is False
    assert first["Complete_Structure_Specific"] is False


def test_same_candidate_multiple_complete_structures():
    _, detail = run([
        ion(candidate="C1", structure="S2", ion_id="I2", mz=500.002),
        ion(candidate="C1", structure="S1", ion_id="I1", mz=500.0005),
    ])
    assert all(row["Candidate_Specific"] is True for row in detail)
    assert all(row["Complete_Structure_Specific"] is False for row in detail)
    assert detail[0]["Competing_Candidate_Count"] == 0
    assert detail[0]["Competing_Complete_Structure_Count"] == 1


def test_best_assignment_matches_legacy_rule_and_margin():
    ions = [
        ion(candidate="B", ion_id="I0", mz=500.0005),
        ion(candidate="A", ion_id="I9", mz=500.0015),
    ]
    best, detail = run(ions)
    # Equal absolute errors are resolved by Candidate_ID, as in the legacy prefix.
    assert best[0]["Candidate_ID"] == "A"
    assert detail[0]["Candidate_ID"] == "A"
    assert detail[0]["Best_Error_ppm"] == pytest.approx(
        abs(best[0]["Mass_Error_ppm"])
    )
    assert detail[0]["Second_Best_Error_ppm"] == pytest.approx(
        abs(detail[1]["Mass_Error_ppm"])
    )
    assert detail[0]["Best_vs_Second_Error_Margin_ppm"] == pytest.approx(
        detail[0]["Second_Best_Error_ppm"] - detail[0]["Best_Error_ppm"]
    )


def test_exact_duplicate_only_is_deduplicated():
    duplicate = ion()
    _, detail = run([
        duplicate, deepcopy(duplicate),
        ion(candidate="C1", structure="S1", ion_id="I2", mz=500.002),
        ion(candidate="C2", structure="S2", ion_id="I3", mz=500.003),
    ])
    assert len(detail) == 3
    assert {row["Ion_ID"] for row in detail} == {"I1", "I2", "I3"}


def test_input_order_is_deterministic():
    ions = [
        ion(candidate="C2", structure="S2", ion_id="I3", mz=500.003),
        ion(candidate="C1", structure="S1", ion_id="I2", mz=500.002),
        ion(candidate="C1", structure="S1", ion_id="I1", mz=500.0005),
    ]
    first = run(ions)
    second = run(list(reversed(deepcopy(ions))))
    assert first == second


def test_zero_intensity_state_and_formal_flags():
    best, detail = run([ion()], intensity=0.0)
    assert best[0]["Observed_Intensity"] == 0.0
    assert detail[0]["Observed_Intensity_State"] == "zero"
    for row in detail:
        assert row["Applied_To_Formal_Result"] is False
        assert row["Formal_Change_Ready"] is False
        assert row["Formal_Result_Changed"] is False


def test_default_return_is_legacy_best_list():
    ions = [ion(ion_id="I2", mz=500.002), ion(ion_id="I1", mz=500.0005)]
    legacy = match_composite_ms2(ions, [spectrum()], config())
    best, _detail = run(ions)
    assert legacy == best
    assert set(legacy[0]) == {
        "Candidate_ID", "Complete_Structure_ID", "Spectrum_ID", "Precursor_mz",
        "Precursor_Charge", "Ion_Series", "Ion_Number", "Cleavage_Position",
        "Included_Positions", "Included_Modified_Positions",
        "Included_Backbone_Bonds", "Theoretical_Neutral_Mass",
        "Theoretical_mz", "Observed_mz", "Mass_Error_Da", "Mass_Error_ppm",
        "Observed_Intensity", "Position_Informative", "Backbone_Informative",
        "Candidate_Discriminating", "Isomer_Discriminating",
        "Legacy_Competition_Class", "Audit_Level",
        "Applied_To_Formal_Result", "Formal_Change_Ready",
    }


def test_sheet_policy_standard_audit_full():
    logical = "Composite_MS2_Assignment_Competition"
    alias = "Composite_MS2_Assignment_Compe"
    for name in (logical, alias):
        standard, _ = included_sheet_names([name], AuditPolicy.from_level("standard"))
        audit, _ = included_sheet_names([name], AuditPolicy.from_level("audit"))
        full, _ = included_sheet_names([name], AuditPolicy.from_level("full"))
        assert standard == [] and audit == [] and full == [name]
