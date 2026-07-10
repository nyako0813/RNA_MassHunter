from pathlib import Path

import yaml

from rna_masshunter.modification_search import _candidate_policy_allows_mass_search
from rna_masshunter.modifications import load_modifications
from tools.import_curated_modifications import curated_record


def test_curated_policy_and_loader_metadata(tmp_path: Path | None = None):
    pseudouridine = curated_record({
        "id": "Y", "symbol": "Y", "common_name": "pseudouridine", "base": "U", "target_bases": "U",
        "modified_nucleoside_mass_mono": 244.0695, "mass_shift_from_unmodified": 0.0,
        "mass_basis": "nucleoside_delta", "category": "biological", "detectability_ms1": True,
        "detectability_ms2": "limited", "source_priority": "user_pdf_for_mass_shift",
        "curation_status": "manually_checked",
    })
    assert pseudouridine["mass_shift_from_unmodified"] == 0.0
    assert pseudouridine["detectability"]["ms1"] is False
    assert pseudouridine["candidate_policy"]["include_by_mass_search"] is False
    assert pseudouridine["candidate_policy"]["include_if_position_rule_exists"] is True

    trimethyl = curated_record({
        "id": "tri", "symbol": "tri", "common_name": "trimethyl nucleoside", "target_bases": "A",
        "modified_nucleoside_mass_mono": 300.0, "mass_shift_from_unmodified": 42.0,
        "category": "biological", "chemical_group": "trimethylation",
    })
    acetyl = curated_record({
        "id": "ac", "symbol": "ac", "common_name": "acetyl nucleoside", "target_bases": "C",
        "modified_nucleoside_mass_mono": 300.0, "mass_shift_from_unmodified": 42.0,
        "category": "biological", "chemical_group": "acetylation",
    })
    assert trimethyl["isobaric_group"] == "trimethylation_group"
    assert acetyl["isobaric_group"] == "acetylation_group"
    assert trimethyl["near_isobaric_group"] == acetyl["near_isobaric_group"] == "near_isobaric_42Da_group"

    path = Path("output/test_curated_loader.yaml")
    path.parent.mkdir(exist_ok=True)
    path.write_text(yaml.safe_dump({"modifications": [pseudouridine]}, sort_keys=False), encoding="utf-8")
    loaded = load_modifications(path)
    assert loaded[0].curation_status == "manually_checked"
    assert loaded[0].source_priority == "user_pdf_for_mass_shift"
    assert loaded[0].candidate_policy["include_by_mass_search"] is False
    assert _candidate_policy_allows_mass_search(loaded[0]) is False


if __name__ == "__main__":
    test_curated_policy_and_loader_metadata()
    print("synthetic curated modification tests: OK")
