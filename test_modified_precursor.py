from types import SimpleNamespace

from rna_masshunter.models import Fragment, Modification, MS2SpectrumInfo
from rna_masshunter.ms1_mapping import theoretical_mz_from_mass
from rna_masshunter.ms2_annotation import find_parent_candidates, parent_candidate_rows


def test_modified_precursor_rescue_and_isobaric_exclusion():
    fragment = Fragment("F1", "T", "AUC", 1, 3, None, None, "RNase_T1", 0, "default", 1000.0)
    modification = Modification("mX", "mX", 14.0, "test", ["C"], raw={"name": "C mod"})
    config = SimpleNamespace(
        ms2_annotation={
            "precursor_match_tolerance_ppm": 20,
            "include_modified_precursor_candidates": True,
            "modified_precursor_require_base_compatibility": True,
            "modified_precursor_include_isobaric": False,
            "modified_precursor_max_candidates_per_spectrum": 20,
            "ion_charge_states": [1],
        },
        instrument={"polarity": "negative"},
    )
    mz = theoretical_mz_from_mass(1014.0, 1, "negative")
    spectrum = MS2SpectrumInfo("S1", 1, 1.0, mz, 1, None, 0, None, None, 0.0, [])
    candidates = find_parent_candidates(spectrum, [fragment], config, [modification])
    rows = parent_candidate_rows(spectrum, candidates)
    assert len(candidates) == 1
    assert candidates[0]["candidate_type"] == "modified"
    assert candidates[0]["modified_rescue"] is True
    assert rows[0]["Parent_Match_Status"] == "matched_modified"
    assert rows[0]["Modified_Precursor_Rescue"] is True

    isobaric = Modification("psi", "psi", 0.0, "test", ["U"])
    assert find_parent_candidates(spectrum, [fragment], config, [isobaric]) == []


if __name__ == "__main__":
    test_modified_precursor_rescue_and_isobaric_exclusion()
    print("synthetic modified precursor test: OK")
