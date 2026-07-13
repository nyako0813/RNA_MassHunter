from dataclasses import asdict

import pandas as pd
import pytest

from rna_masshunter.models import Fragment, Peak, RunConfig
from rna_masshunter.ms1_mapping import map_fragments_to_ms1_peaks, theoretical_mz_from_mass
from rna_masshunter.ms1_match_truncation_audit import (
    AUDIT_COLUMNS, DETAIL_COLUMNS, SUMMARY_COLUMNS, TOP_COLUMNS,
    append_diagnostic_shadow_columns, append_top_shadow_columns,
    build_ms1_truncation_audit,
)


def config(limit=20, report_limit=100000, max_charge=1):
    return RunConfig(
        instrument={"polarity": "negative"},
        fragment_mapping={
            "enabled": True, "polarity": "negative", "min_charge": 1,
            "max_charge": max_charge, "mz_tolerance_ppm": 10,
            "max_matches_per_fragment": limit, "use_peak_tiers": True,
            "include_trace_peaks": True, "min_fragment_length_for_filtered": 3,
            "filtered_peak_tiers": ["Major", "Minor"],
            "filtered_confidence": ["High", "Medium"],
        },
        modification_search={"enabled": False},
        modification_evidence_ranking={"enabled": False},
        reporting={"max_excel_rows_per_sheet": report_limit},
    )


def fragment():
    return Fragment("F1", "target", "ACGU", 1, 4, 1, 4, "RNase_T1", 0, "default", 1000.0)


def peaks(count, tiers=None, intensities=None):
    mz = theoretical_mz_from_mass(1000.0, 1, "negative")
    values = []
    for i in range(count):
        ppm = ((i // 2) + 1) * (1 if i % 2 else -1) * 0.15
        values.append(Peak(
            mz=mz * (1 + ppm / 1_000_000),
            intensity=(intensities[i] if intensities else 1000 - i),
            rt=1.0 + i / 1000,
            scan_id=f"S{i}",
            tier=(tiers[i] if tiers else "Major"),
        ))
    return values


def run(count, limit=20, report_limit=100000, tiers=None, intensities=None):
    cfg = config(limit, report_limit)
    ctx = {}
    formal = map_fragments_to_ms1_peaks([fragment()], peaks(count, tiers, intensities), cfg, audit_context=ctx)
    audit = build_ms1_truncation_audit(ctx, cfg, [], [], formal, [], [], {})
    return formal, ctx, audit


@pytest.mark.parametrize("count", [0, 1, 20, 21, 100])
def test_boundary_counts(count):
    formal, ctx, audit = run(count)
    assert len(formal) == min(count, 20)
    assert audit["summary"]["Total_Pre_Truncation_Matches"] == count
    assert audit["summary"]["Total_Discarded_Matches"] == max(0, count - 20)
    assert len(ctx["fragments"]) == 1


def test_unlimited_formal_configuration_is_normalized_without_shadow_mutation():
    cfg = config(limit=0)
    ctx = {}
    formal = map_fragments_to_ms1_peaks([fragment()], peaks(21), cfg, audit_context=ctx)
    assert len(formal) == 20
    assert ctx["configured_max_matches"] == 20


def test_formal_result_is_identical_with_and_without_capture():
    cfg = config()
    baseline = map_fragments_to_ms1_peaks([fragment()], peaks(100), cfg)
    captured = {}
    audited = map_fragments_to_ms1_peaks([fragment()], peaks(100), cfg, audit_context=captured)
    assert [asdict(x) for x in baseline] == [asdict(x) for x in audited]


def test_deterministic_sort_and_ties_use_peak_input_order():
    mz = theoretical_mz_from_mass(1000.0, 1, "negative")
    tied = [Peak(mz, 100, 1.0, f"S{i}", tier="Major") for i in range(25)]
    cfg = config()
    ctx = {}
    first = map_fragments_to_ms1_peaks([fragment()], tied, cfg, audit_context=ctx)
    assert [getattr(x, "_audit_peak_index") for x in first] == list(range(1, 21))
    ctx2 = {}
    second = map_fragments_to_ms1_peaks([fragment()], tied, cfg, audit_context=ctx2)
    assert [asdict(x) for x in first] == [asdict(x) for x in second]


def test_error_precedes_intensity_in_formal_sort():
    mz = theoretical_mz_from_mass(1000.0, 1, "negative")
    ps = [Peak(mz * (1 + (i + 1) / 1_000_000), 1, scan_id=str(i), tier="Major") for i in range(20)]
    ps.append(Peak(mz * (1 + 0.01 / 1_000_000), 1_000_000, scan_id="best", tier="Major"))
    cfg = config(); ctx = {}
    result = map_fragments_to_ms1_peaks([fragment()], ps, cfg, audit_context=ctx)
    assert any(x.scan_id == "best" for x in result)


def test_intensity_breaks_equal_error_tie():
    mz = theoretical_mz_from_mass(1000.0, 1, "negative")
    ps = [Peak(mz * (1 + 1 / 1_000_000), i, scan_id=str(i), tier="Major") for i in range(21)]
    cfg = config(); ctx = {}
    result = map_fragments_to_ms1_peaks([fragment()], ps, cfg, audit_context=ctx)
    assert "0" not in {x.scan_id for x in result}
    assert "20" in {x.scan_id for x in result}


def test_retained_discarded_and_quality_flags():
    formal, _, audit = run(21)
    row = audit["audit_rows"][0]
    assert row["Retained_Match_Count"] == 20
    assert row["Discarded_Match_Count"] == 1
    assert sum(x["Retained_Or_Discarded"] == "retained" for x in audit["detail_rows"]) == 20
    assert sum(x["Retained_Or_Discarded"] == "discarded" for x in audit["detail_rows"]) == 1


def test_filter_first_selection_detects_trace_slot_occupation():
    tiers = ["Trace"] * 20 + ["Major"]
    intensities = list(range(100, 79, -1))
    _, _, audit = run(21, tiers=tiers, intensities=intensities)
    summary = audit["summary"]
    assert summary["Filter_First_Filter_Passing_Retained"] > summary["Current_Sort_Filter_Passing_Retained"]
    assert audit["audit_rows"][0]["Discarded_Filter_Passing_Count"] == 1


def test_physical_peak_and_charge_bias_are_reported():
    mz1 = theoretical_mz_from_mass(1000.0, 1, "negative")
    mz2 = theoretical_mz_from_mass(1000.0, 2, "negative")
    ps = [Peak(mz1 * (1 + i * 0.05 / 1_000_000), 100-i, scan_id=f"a{i}", tier="Major") for i in range(20)]
    ps += [Peak(mz2 * (1 + (5 + i * 0.05) / 1_000_000), 100-i, scan_id=f"b{i}", tier="Major") for i in range(5)]
    cfg = config(max_charge=2); ctx = {}
    formal = map_fragments_to_ms1_peaks([fragment()], ps, cfg, audit_context=ctx)
    audit = build_ms1_truncation_audit(ctx, cfg, [], [], formal, [], [], {})
    assert audit["audit_rows"][0]["New_Charge_State_Only_In_Discarded"]
    assert audit["audit_rows"][0]["New_Physical_Peak_Only_In_Discarded"]


def test_detail_truncation_is_deterministic_and_summarized():
    _, _, first = run(100, report_limit=17)
    _, _, second = run(100, report_limit=17)
    assert first["detail_rows"] == second["detail_rows"]
    assert first["summary"]["Detail_Original_Row_Count"] == 100
    assert first["summary"]["Detail_Written_Row_Count"] == 17
    assert first["summary"]["Detail_Truncated"]


def test_required_columns_and_sheet_names():
    assert len("MS1_Truncation_Audit") <= 31
    assert len("MS1_Truncation_Detail") <= 31
    assert len("MS1_Truncation_Summary") <= 31
    assert len(AUDIT_COLUMNS) == len(set(AUDIT_COLUMNS))
    assert len(DETAIL_COLUMNS) == len(set(DETAIL_COLUMNS))
    assert len(SUMMARY_COLUMNS) == len(set(SUMMARY_COLUMNS))
    assert "Applied_To_Final_Score" in AUDIT_COLUMNS
    assert "Applied_To_Final_Score" in SUMMARY_COLUMNS


def test_top_columns_append_at_right_without_changing_existing_values():
    _, _, audit = run(21)
    top = pd.DataFrame([{"Review_Rank": 1, "Modification_ID": "m1", "Parent_Fragment_ID": "F1", "Candidate_Positions_In_tRNA": 36, "Best_Final_Score": 7.0}])
    result = append_top_shadow_columns(top, audit)
    assert list(result.columns[:len(top.columns)]) == list(top.columns)
    assert list(result.columns[-len(TOP_COLUMNS):]) == TOP_COLUMNS
    assert result.iloc[0]["Best_Final_Score"] == 7.0
    assert not bool(result.iloc[0]["MS1_Truncation_Applied_To_Final_Score"])


def test_diagnostics_append_and_formal_false():
    _, _, audit = run(21)
    result = append_diagnostic_shadow_columns([{"Existing": 1}], audit)
    assert result[0]["Existing"] == 1
    assert result[0]["MS1_Discarded_Match_Count"] == 1
    assert result[0]["MS1_Truncation_Applied_To_Final_Score"] is False


def test_all_applied_flags_are_false():
    _, _, audit = run(100)
    assert all(row["Applied_To_Final_Score"] is False for row in audit["audit_rows"])
    assert audit["summary"]["Applied_To_Final_Score"] is False
