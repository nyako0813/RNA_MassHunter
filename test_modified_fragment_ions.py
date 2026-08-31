from pathlib import Path
from types import SimpleNamespace

from rna_masshunter.masses import load_base_masses
from rna_masshunter.models import MS2SpectrumInfo
from rna_masshunter.modified_fragment_ions import (
    annotate_position_discrimination,
    build_localization_evidence,
    generate_modified_theoretical_ions,
    match_modified_ions,
)


def _config():
    return SimpleNamespace(
        instrument={"polarity": "negative"},
        ms2_annotation={
            "include_modified_fragment_ions": True,
            "modified_fragment_require_target_base": True,
            "modified_fragment_include_unmodified_counterparts": True,
            "modified_fragment_max_positions_per_candidate": 20,
            "modified_fragment_min_ion_length": 1,
            "modified_fragment_min_ion_length_for_localization": 2,
            "modified_fragment_max_rows": 100000,
            "mz_tolerance_ppm": 20,
        },
    )


def _parent():
    return {
        "Spectrum_ID": "S1", "Candidate_Parent_Fragment_ID": "F1",
        "Candidate_Parent_Sequence": "ACGU", "Candidate_Parent_Start": 10,
        "Candidate_Parent_End": 13, "Candidate_Type": "modified",
        "Modification_ID": "mC", "Modification_Name": "C modification",
        "Modification_Target_Base": "C", "Modification_Mass_Shift": 14.0156,
        "Parent_Charge": 1,
    }


def test_modified_ion_generation_matching_and_localization():
    base_masses = load_base_masses(Path(__file__).parent / "data" / "base_masses.yaml")
    ions = generate_modified_theoretical_ions([_parent()], _config(), base_masses)
    assert ions
    assert {ion["Candidate_Modification_Position_In_Parent"] for ion in ions} == {2}
    modified = [ion for ion in ions if ion["Ion_Contains_Modification"] and ion["Informative_Ion"]]
    selected = []
    for ion_type in ("d", "w"):
        selected.extend([ion for ion in modified if ion["Ion_Type"] == ion_type][:2])
    assert len(selected) >= 3 and {ion["Ion_Type"] for ion in selected} == {"d", "w"}
    peaks = [(ion["Theoretical_mz"], 1000.0 - index * 10) for index, ion in enumerate(selected)]
    spectrum = MS2SpectrumInfo("S1", 1, 1.5, 1000.0, 1, None, len(peaks), peaks[0][0], 1000.0, 3000.0, peaks)
    matches = match_modified_ions([spectrum], ions, _config())
    assert any(row["Match_Status"] == "matched_modified_ion" for row in matches)
    evidence = build_localization_evidence(ions, matches)
    assert evidence[0]["Candidate_Modification_Position_In_Parent"] == 2
    assert evidence[0]["Candidate_Modification_Position_In_tRNA"] == 11
    assert evidence[0]["Localization_Level"] == "Strong"


def test_one_nt_support_is_not_strong_and_multiple_positions_are_ambiguous():
    ion_template = {
        "Spectrum_ID": "S2", "Parent_Fragment_ID": "F2", "Modification_ID": "mC",
        "Modification_Name": "C modification", "Parent_Sequence": "CCG", "Parent_Start": 20,
        "Parent_End": 22, "Candidate_Modification_Base": "C",
    }
    ions = []
    matches = []
    for position in (1, 2):
        for index, ion_type in enumerate(("d", "d", "w"), start=1):
            ion = {**ion_template, "Candidate_Modification_Position_In_Parent": position,
                   "Ion_ID": f"I{position}{index}", "Ion_Type": ion_type}
            ions.append(ion)
            matches.append({
                "Spectrum_ID": "S2", "Parent_Fragment_ID": "F2", "Modification_ID": "mC",
                "Candidate_Modification_Position_In_Parent": position, "Ion_Contains_Modification": True,
                "Informative_Ion": True, "Ion_Type": ion_type, "Mass_Error_ppm": 1.0,
                "Observed_Intensity": 100.0, "RT": 1.0, "Precursor_mz": 500.0,
            })
    evidence = build_localization_evidence(ions, matches)
    assert all(row["Localization_Level"] == "Strong" for row in evidence)
    assert all(row["Localization_Interpretation"] == "ambiguous-multiple-positions" for row in evidence)

    weak_ion = {**ion_template, "Candidate_Modification_Position_In_Parent": 1, "Ion_ID": "weak", "Ion_Type": "d"}
    weak_match = {
        "Spectrum_ID": "S2", "Parent_Fragment_ID": "F2", "Modification_ID": "mC",
        "Candidate_Modification_Position_In_Parent": 1, "Ion_Contains_Modification": True,
        "Informative_Ion": False, "Ion_Type": "d", "Mass_Error_ppm": 1.0,
        "Observed_Intensity": 100.0, "RT": 1.0, "Precursor_mz": 500.0,
    }
    weak = build_localization_evidence([weak_ion], [weak_match])
    assert weak[0]["Localization_Level"] == "Weak"


def test_generates_a_and_z_ions_with_correct_mass_offsets():
    from rna_masshunter.masses import calculate_unmodified_rna_mass, fragment_ion_series_offsets

    base_masses = load_base_masses(Path(__file__).parent / "data" / "base_masses.yaml")
    ions = generate_modified_theoretical_ions([_parent()], _config(), base_masses)
    assert {ion["Ion_Type"] for ion in ions} == {"d", "w", "a", "z"}

    offsets = fragment_ion_series_offsets(base_masses)
    unmodified_counterparts = [ion for ion in ions if not ion["Ion_Contains_Modification"]]
    for ion in unmodified_counterparts:
        base_mass = calculate_unmodified_rna_mass(ion["Ion_Sequence"], base_masses, terminal_form="default")
        assert abs(ion["Theoretical_Mass"] - (base_mass + offsets[ion["Ion_Type"]])) < 1e-9


def test_ion_series_config_restricts_generated_modified_ion_types():
    config = _config()
    config.ms2_annotation["ion_series"] = ["d", "w"]
    base_masses = load_base_masses(Path(__file__).parent / "data" / "base_masses.yaml")
    ions = generate_modified_theoretical_ions([_parent()], config, base_masses)
    assert ions
    assert {ion["Ion_Type"] for ion in ions} == {"d", "w"}


def test_position_discriminating_ion_excludes_alternative_position():
    ions = [
        {"Spectrum_ID": "S3", "Parent_Fragment_ID": "F3", "Modification_ID": "mA",
         "Candidate_Modification_Position_In_Parent": position, "Parent_Start": 10}
        for position in (1, 3)
    ]
    matches = [{
        "Spectrum_ID": "S3", "Parent_Fragment_ID": "F3", "Modification_ID": "mA",
        "Candidate_Modification_Position_In_Parent": 1, "Ion_Start": 1, "Ion_End": 1,
        "Informative_Ion": True, "Ion_Contains_Modification": True,
    }]
    annotated = annotate_position_discrimination(matches, ions, _config())
    assert annotated[0]["Position_Discriminating_Ion"] is True
    assert annotated[0]["Discriminates_Position"] == 10
    assert annotated[0]["Also_Explains_Positions"] == ""


if __name__ == "__main__":
    test_modified_ion_generation_matching_and_localization()
    test_one_nt_support_is_not_strong_and_multiple_positions_are_ambiguous()
    test_generates_a_and_z_ions_with_correct_mass_offsets()
    test_ion_series_config_restricts_generated_modified_ion_types()
    test_position_discriminating_ion_excludes_alternative_position()
    print("synthetic modified fragment ion tests: OK")
