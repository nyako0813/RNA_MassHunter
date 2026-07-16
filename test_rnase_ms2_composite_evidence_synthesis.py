from copy import deepcopy

from rna_masshunter.audit_policy import AuditPolicy, included_sheet_names
from rna_masshunter.rnase_ms2_composite_evidence_synthesis import (
    build_rnase_ms2_composite_evidence_synthesis,
)


def assignment(key="P1", candidate="C1", structure="S1", ion="I1", intensity=100.0,
               candidate_specific=True, structure_specific=True, ion_specific=True,
               position=False, backbone=False, rank=1, candidate_comp=0,
               structure_comp=0, ion_comp=0):
    return {
        "Composite_Match_ID": f"M-{key}-{candidate}-{structure}-{ion}",
        "Physical_Observed_Peak_Key": key, "Spectrum_ID": "SP1",
        "Observed_Peak_Index": 0, "Raw_Peak_Index": "", "RT": 1.0,
        "Observed_mz": 500.0, "Observed_Intensity": intensity,
        "Observed_Intensity_State": "positive" if intensity > 0 else "zero",
        "Candidate_ID": candidate, "Complete_Structure_ID": structure,
        "Ion_ID": ion, "Parent_Fragment_ID": "F1", "Ion_Series": "c",
        "Ion_Number": 2, "Cleavage_Position": 2,
        "Included_Modified_Positions": "2" if position else "",
        "Included_Backbone_Bonds": "2_3" if backbone else "",
        "Position_Informative": position, "Backbone_Informative": backbone,
        "Candidate_Specific": candidate_specific,
        "Complete_Structure_Specific": structure_specific,
        "Theoretical_Ion_Specific": ion_specific,
        "Position_Specific": position and candidate_specific and structure_specific,
        "Backbone_Bond_Specific": backbone and candidate_specific and structure_specific,
        "Assignment_Rank": rank, "Best_Assignment": rank == 1,
        "Within_Tolerance_Assignment_Count": 1 + max(candidate_comp, structure_comp, ion_comp),
        "Competing_Candidate_Count": candidate_comp,
        "Competing_Candidate_IDs": "C2" if candidate_comp else "",
        "Competing_Complete_Structure_Count": structure_comp,
        "Competing_Complete_Structure_IDs": "S2" if structure_comp else "",
        "Competing_Theoretical_Ion_Count": ion_comp,
        "Competing_Ion_IDs": "I2" if ion_comp else "", "Mass_Error_ppm": 1.0,
    }


def build(rows=(), *, ions=(), support=(), compare=(), matches=()):
    return build_rnase_ms2_composite_evidence_synthesis(
        list(ions), list(matches), list(rows), list(support), list(compare),
    )


def candidate(result, candidate="C1", structure="S1"):
    return next(row for row in result.evidence_rows
                if row["Candidate_ID"] == candidate and row["Complete_Structure_ID"] == structure)


def test_empty_input():
    result = build()
    assert result.evidence_rows == [] and result.peak_rows == []
    assert result.summary_rows[0]["Composite_Candidate_Count"] == 0


def test_positive_match_is_only_provisional_and_structure_unresolved():
    result = build([assignment()])
    row = candidate(result)
    assert row["Composite_Identity_Status"] == "PROVISIONAL_CANDIDATE_SUPPORT"
    assert row["Composite_Structure_Status"] == "UNRESOLVED"
    assert result.summary_rows[0]["Structure_SUPPORTED_Count"] == 0


def test_candidate_structure_and_ion_competition_are_ambiguous():
    for field in ("candidate_comp", "structure_comp", "ion_comp"):
        kwargs = {field: 1, "candidate_specific": False,
                  "structure_specific": False, "ion_specific": False}
        row = candidate(build([assignment(**kwargs)]))
        assert row["Composite_Identity_Status"] == "AMBIGUOUS"
        assert row["Composite_Localization_Status"] == "AMBIGUOUS"
        assert row["Composite_Backbone_Status"] == "AMBIGUOUS"
    row = candidate(build([assignment(structure_comp=1, structure_specific=False)]))
    assert row["Composite_Structure_Status"] == "AMBIGUOUS"


def test_position_support_one_and_multiple_independent_peaks():
    one = candidate(build([assignment(position=True)]))
    assert one["Composite_Localization_Status"] == "POSITION_COMPATIBLE"
    two = candidate(build([assignment(position=True), assignment(key="P2", ion="I2", position=True)]))
    assert two["Composite_Localization_Status"] == "PARTIALLY_SUPPORTED"
    assert two["Position_Informative_Support_Peak_Count"] == 2


def test_backbone_support_one_and_multiple_independent_peaks():
    one = candidate(build([assignment(backbone=True)]))
    assert one["Composite_Backbone_Status"] == "BACKBONE_COMPATIBLE"
    two = candidate(build([assignment(backbone=True), assignment(key="P2", ion="I2", backbone=True)]))
    assert two["Composite_Backbone_Status"] == "PARTIALLY_SUPPORTED"


def test_shared_assignment_is_excluded_from_individual_support():
    result = build([assignment(candidate_specific=False)])
    row = candidate(result)
    assert row["Individual_Support_Physical_Peak_Count"] == 0
    assert row["Shared_Physical_Peak_Count"] == 1
    assert result.peak_rows[0]["Counts_For_Individual_Support"] is False


def test_zero_intensity_is_not_support():
    result = build([assignment(intensity=0.0)])
    row = candidate(result)
    assert row["Composite_Identity_Status"] == "UNSUPPORTED"
    assert result.peak_rows[0]["Positive_Assignment"] is False


def test_structural_isomer_is_ambiguous_and_never_supported():
    result = build([assignment()], compare=[{"Candidate_ID": "C1", "Comparison_Class": "BOTH_ISOMERIC"}])
    row = candidate(result)
    assert row["Composite_Structure_Status"] == "AMBIGUOUS"
    assert all(item["Composite_Structure_Status"] != "SUPPORTED" for item in result.evidence_rows)


def test_input_order_is_deterministic():
    rows = [assignment(key="P2", candidate="C2", structure="S2", ion="I2"), assignment()]
    assert build(rows) == build(list(reversed(deepcopy(rows))))


def test_all_formal_flags_false():
    result = build([assignment()])
    for rows in (result.summary_rows, result.evidence_rows, result.peak_rows):
        for row in rows:
            assert row["Applied_To_Formal_Result"] is False
            assert row["Formal_Change_Ready"] is False
            assert row["Formal_Result_Changed"] is False


def test_sheet_inclusion():
    summary = "RNase_MS2_Composite_Summary"
    evidence = "RNase_MS2_Composite_Evidence"
    peak_alias = "RNase_MS2_Composite_Peak_Eviden"
    names = [summary, evidence, peak_alias]
    standard, _ = included_sheet_names(names, AuditPolicy.from_level("standard"))
    audit, _ = included_sheet_names(names, AuditPolicy.from_level("audit"))
    full, _ = included_sheet_names(names, AuditPolicy.from_level("full"))
    assert standard == []
    assert audit == [summary]
    assert full == names
