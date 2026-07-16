from copy import deepcopy

from rna_masshunter.rnase_ms2_evidence_synthesis import (
    CANDIDATE_EVIDENCE_COLUMNS,
    PEAK_EVIDENCE_COLUMNS,
    SUMMARY_COLUMNS,
    build_rnase_ms2_evidence_synthesis,
)


def build(**overrides):
    values = {
        "ranking_rows": [], "ambiguity_groups": [], "modified_precursors": [],
        "modified_theoretical_ions": [], "modified_ion_matches": [],
        "localization_rows": [], "identity_rows": [],
        "identity_peak_assignments": [], "ambiguous_clusters": [],
        "ambiguous_peak_details": [], "effective_ambiguity_rows": [],
        "effective_ambiguity_details": [],
    }
    values.update(overrides)
    return build_rnase_ms2_evidence_synthesis(**values)


def ranking(mod="mod", parent="F1", position=2, trna=12):
    return {
        "Modification_ID": mod, "Modification_Name": mod,
        "Parent_Fragment_ID": parent, "Candidate_Position_In_Parent": position,
        "Candidate_tRNA_Position": trna, "Candidate_Base": "A",
    }


def match(mod="mod", parent="F1", position=2, *, intensity=100.0,
          informative=False, discriminating=False, ion="ION1", observed=500.0,
          theoretical=500.001):
    return {
        "Spectrum_ID": "S1", "Scan_Index": 7, "Modification_ID": mod,
        "Parent_Fragment_ID": parent,
        "Candidate_Modification_Position_In_Parent": position,
        "Ion_Contains_Modification": True, "Informative_Ion": informative,
        "Position_Discriminating_Ion": discriminating,
        "Observed_mz": observed, "Observed_Intensity": intensity,
        "Ion_ID": ion, "Theoretical_mz": theoretical,
    }


def match_id(row):
    values = (
        row.get("Spectrum_ID"), row.get("Scan_Index"), row.get("Observed_mz"),
        row.get("Ion_ID"), row.get("Theoretical_mz"),
    )
    return ":".join(str(value if value not in (None, "") else "NA") for value in values)


def assignment(row, *, scope="candidate_specific", counts=True,
               theories=1, isomer=False, candidates=1):
    return {
        "Physical_Observed_Peak_Key": f"{row['Spectrum_ID']}|peak={row['Ion_ID']}",
        "Match_ID": match_id(row), "Spectrum_ID": row["Spectrum_ID"],
        "Observed_mz": row["Observed_mz"],
        "Observed_Intensity": row["Observed_Intensity"],
        "Modification_ID": row["Modification_ID"],
        "Parent_Fragment_ID": row["Parent_Fragment_ID"],
        "Candidate_Position_In_Parent": row["Candidate_Modification_Position_In_Parent"],
        "Candidate_tRNA_Position": 12, "Theoretical_Ion_ID": row["Ion_ID"],
        "Theoretical_mz": row["Theoretical_mz"],
        "Physical_Peak_Assignment_Count": max(theories, candidates),
        "Physical_Peak_Candidate_Count": candidates,
        "Physical_Peak_Theoretical_Ion_Count": theories,
        "Physical_Peak_Shared_Across_Candidates": candidates > 1,
        "Physical_Peak_Shared_Across_Isomers": isomer,
        "Physical_Peak_Assignment_Status": scope,
        "Evidence_Scope": scope, "Counts_For_Individual_Identity": counts,
    }


def candidate_row(result):
    assert len(result.candidate_rows) == 1
    return result.candidate_rows[0]


def test_empty_input():
    result = build()
    assert list(result.summary_rows[0]) == SUMMARY_COLUMNS
    assert result.summary_rows[0]["Candidate_Count"] == 0
    assert result.candidate_rows == [] and result.peak_rows == []
    assert result.summary_rows[0]["Applied_To_Formal_Result"] is False


def test_precursor_only():
    row = candidate_row(build(modified_precursors=[{
        "Spectrum_ID": "S1", "Modification_ID": "mod",
        "Parent_Fragment_ID": "F1",
    }]))
    assert row["Modification_Identity_Status"] == "PRECURSOR_COMPATIBLE"
    assert row["Localization_Status"] == "UNRESOLVED"
    assert row["Structure_Status"] == "NOT_EVALUATED"


def test_candidate_specific_modified_peak_supports_identity():
    ion = match()
    row = candidate_row(build(
        ranking_rows=[ranking()], modified_ion_matches=[ion],
        identity_peak_assignments=[assignment(ion)],
    ))
    assert row["Modification_Identity_Status"] == "FRAGMENT_SUPPORTED"
    assert row["Candidate_Specific_Physical_Peak_Count"] == 1
    assert row["Localization_Status"] == "UNRESOLVED"
    assert row["Structure_Status"] == "UNRESOLVED"
    assert row["Ambiguity_Status"] == "NONE"


def test_shared_peak_does_not_support_individual_identity_or_localization():
    ion = match(informative=True, discriminating=True)
    shared = assignment(
        ion, scope="position_group_level", counts=False, candidates=2,
    )
    result = build(
        ranking_rows=[ranking()], modified_ion_matches=[ion],
        identity_peak_assignments=[shared],
    )
    row = candidate_row(result)
    assert row["Modification_Identity_Status"] == "AMBIGUOUS"
    assert row["Localization_Status"] == "AMBIGUOUS"
    assert row["Candidate_Specific_Physical_Peak_Count"] == 0
    assert row["Shared_Physical_Peak_Count"] == 1
    assert row["Ambiguity_Status"] == "POSITION_GROUP"
    assert list(row) == CANDIDATE_EVIDENCE_COLUMNS
    assert list(result.peak_rows[0]) == PEAK_EVIDENCE_COLUMNS


def test_informative_candidate_specific_position_discriminating_ion_localizes():
    ion = match(informative=True, discriminating=True)
    row = candidate_row(build(
        ranking_rows=[ranking()], modified_ion_matches=[ion],
        identity_peak_assignments=[assignment(ion)],
    ))
    assert row["Modification_Identity_Status"] == "FRAGMENT_SUPPORTED"
    assert row["Localization_Status"] == "LOCALIZED"
    assert row["Position_Discriminating_Ion_Count"] == 1


def test_single_candidate_position_alone_is_not_localized():
    ion = match()
    row = candidate_row(build(
        ranking_rows=[ranking()], modified_ion_matches=[ion],
        identity_peak_assignments=[assignment(ion)],
        ambiguity_groups=[{
            "Modification_ID": "mod", "Parent_Fragment_ID": "F1",
            "Position_Ambiguity_Status": "single_candidate_position",
            "Candidate_Positions_In_Parent": "2",
        }],
    ))
    assert row["Localization_Status"] == "UNRESOLVED"
    assert "SINGLE_CANDIDATE_POSITION_DOES_NOT_ESTABLISH_LOCALIZATION" in row["Limiting_Reasons"]


def test_competing_theoretical_ions_limit_identity_and_localization():
    ion = match(informative=True, discriminating=True)
    competing = assignment(ion, theories=2)
    row = candidate_row(build(
        ranking_rows=[ranking()], modified_ion_matches=[ion],
        identity_peak_assignments=[competing],
        ambiguous_peak_details=[{
            "Physical_Observed_Peak_Key": competing["Physical_Observed_Peak_Key"],
            "Theoretical_Ion_Count": 2,
            "Competing_Theoretical_Ion_IDs": "ION1;ION2",
        }],
    ))
    assert row["Modification_Identity_Status"] == "AMBIGUOUS"
    assert row["Localization_Status"] == "AMBIGUOUS"
    assert row["Ambiguity_Status"] == "THEORETICAL_ION_COMPETITION"
    assert row["Competing_Theoretical_Ion_Count"] == 2


def test_structural_isomer_sharing_keeps_structure_ambiguous():
    ion = match()
    shared = assignment(
        ion, scope="structural_isomer_group_level", counts=False,
        isomer=True, candidates=2,
    )
    row = candidate_row(build(
        ranking_rows=[ranking()], modified_ion_matches=[ion],
        identity_peak_assignments=[shared],
    ))
    assert row["Structure_Status"] == "AMBIGUOUS"
    assert row["Ambiguity_Status"] == "STRUCTURAL_ISOMER"


def test_cross_candidate_ambiguity():
    ion = match()
    shared = assignment(
        ion, scope="cross_candidate_ambiguous", counts=False, candidates=2,
    )
    row = candidate_row(build(
        ranking_rows=[ranking()], modified_ion_matches=[ion],
        identity_peak_assignments=[shared],
    ))
    assert row["Modification_Identity_Status"] == "AMBIGUOUS"
    assert row["Ambiguity_Status"] == "CROSS_CANDIDATE"


def test_raw_only_ambiguity_is_not_formal_ambiguity():
    ion = match()
    specific = assignment(ion)
    peak = specific["Physical_Observed_Peak_Key"]
    row = candidate_row(build(
        ranking_rows=[ranking()], modified_ion_matches=[ion],
        identity_peak_assignments=[specific],
        effective_ambiguity_details=[{
            "Physical_Peak_ID": peak, "Effective_Ambiguity_Level": "raw_only",
            "Formal_Match_Ambiguous": False, "Cross_Candidate_Shared": True,
        }],
    ))
    assert row["Modification_Identity_Status"] == "FRAGMENT_SUPPORTED"
    assert row["Ambiguity_Status"] == "NONE"


def test_zero_intensity_peak_is_excluded_from_support():
    ion = match(intensity=0)
    row = candidate_row(build(
        ranking_rows=[ranking()], modified_ion_matches=[ion],
        identity_peak_assignments=[assignment(ion)],
    ))
    assert row["Modification_Identity_Status"] == "UNSUPPORTED"
    assert row["Localization_Status"] == "UNRESOLVED"
    assert row["Candidate_Specific_Physical_Peak_Count"] == 0


def test_structure_supported_requires_existing_structure_resolution():
    ion = match()
    identity = {
        "Modification_ID": "mod", "Parent_Fragment_ID": "F1",
        "Candidate_tRNA_Position": 12,
        "Structure_Resolution_Status": "structure_resolved",
    }
    row = candidate_row(build(
        ranking_rows=[ranking()], identity_rows=[identity],
        modified_ion_matches=[ion],
        identity_peak_assignments=[assignment(ion)],
    ))
    assert row["Structure_Status"] == "SUPPORTED"


def test_input_order_determinism_and_formal_flags():
    ion1 = match(ion="ION2", observed=501.0, theoretical=501.001)
    ion2 = match(ion="ION1")
    values = {
        "ranking_rows": [ranking(position=3, trna=13), ranking()],
        "modified_ion_matches": [
            dict(ion1, Candidate_Modification_Position_In_Parent=3), ion2,
        ],
        "identity_peak_assignments": [
            assignment(dict(ion1, Candidate_Modification_Position_In_Parent=3)),
            assignment(ion2),
        ],
        "modified_precursors": [{
            "Spectrum_ID": "S1", "Modification_ID": "mod",
            "Parent_Fragment_ID": "F1",
        }],
    }
    first = build(**values)
    second = build(**{
        key: list(reversed(deepcopy(rows))) for key, rows in values.items()
    })
    assert first.to_jsonable() == second.to_jsonable()
    for rows in (first.summary_rows, first.candidate_rows, first.peak_rows):
        assert all(row["Applied_To_Formal_Result"] is False for row in rows)
        assert all(row["Formal_Change_Ready"] is False for row in rows)
        assert all(row["Formal_Result_Changed"] is False for row in rows)
