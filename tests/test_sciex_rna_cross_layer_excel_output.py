import pytest
import pandas as pd
from types import SimpleNamespace as NS
from rna_masshunter.models import RunConfig
from main import build_sciex_cross_layer_evidence_optional_results
from rna_masshunter.sciex_rna_cross_layer_evidence_reconciliation import OPTIONAL_RESULT_KEY
from rna_masshunter.excel_report import _sciex_cross_layer_excel_sheets

def full(pattern=(0,18,34,50),rna="TRNA_TEST"):
    return {"source_id":"FULL","rna_identity":rna,"digest_type":"FULL_LENGTH",
        "input_path":"full.mzML","observed_masses":tuple(100+x for x in pattern)}

def source(layer,path,rna="TRNA_TEST"):
    return NS(source_id=layer,rna_identity=rna,digest_type=layer,input_path=path)

def fragment(match_id="M1",sequence="ACG",start=1,end=3,amb="UNAMBIGUOUS"):
    return NS(match_id=match_id,fragment_id=f"F{start}",peak_id=f"P{start}",fragment_sequence=sequence,
        start_position=start,end_position=end,ion_mode="NEGATIVE_DEPROTONATED",observed_apex_mz=500+start,
        fragment_ambiguity_status=amb,candidate_count_for_peak=1 if amb=="UNAMBIGUOUS" else 2,match_block_reasons=())

def t1(*,families=(),matches=None,rna="TRNA_TEST"):
    matches=(fragment(),) if matches is None else tuple(matches)
    return NS(run_summary=source("T1_DIGEST","t1.mzML",rna),fragment_matches=matches,
        state_families=tuple(families),summary=NS(state_family_count=len(families)))

def get_config(enabled=True):
    return RunConfig(
        sequence={"name": "TRNA_TEST"},
        sciex_profile={
            "cross_layer_evidence_reconciliation": {"enabled": enabled}
        },
        input={},
    )

from rna_masshunter.sciex_rna_cross_layer_evidence_reconciliation import OPTIONAL_RESULT_KEY
from rna_masshunter.excel_report import _sciex_cross_layer_excel_sheets

def test_excel_sheets_generated_and_deterministic():
    config = get_config()
    optional_results = {
        "sciex_intact_oxygen_water_state_audit": full(),
        "sciex_t1_fragment_state_series_audit": t1(),
    }
    warnings = []
    res = build_sciex_cross_layer_evidence_optional_results(config, optional_results, warnings)
    cross_layer_result = res[OPTIONAL_RESULT_KEY]

    sheets = _sciex_cross_layer_excel_sheets(cross_layer_result)

    expected_sheets = {
        "XL_Nodes", "XL_Edges", "XL_Hypotheses",
        "XL_Layer_Summary", "XL_Consensus", "XL_Next_Evidence"
    }
    assert set(sheets.keys()) == expected_sheets

    # Check max 31 characters
    for name in sheets:
        assert len(name) <= 31

    # Deterministic output (alphabetical keys logic inside _format)
    # The columns should be sorted alphabetically. We check if columns of a DataFrame are sorted.
    df = sheets["XL_Nodes"]
    if not df.empty:
        columns = list(df.columns)
        assert columns == sorted(columns)

    df_consensus = sheets["XL_Consensus"]
    assert not df_consensus.empty

    # Check formal propagation safeguard inside the output DataFrame
    assert "formal_propagation" not in df_consensus.columns or df_consensus["formal_propagation"].iloc[0] == False
    assert "chemical_identity_assigned" not in df_consensus.columns or df_consensus["chemical_identity_assigned"].iloc[0] == False

def test_excel_disabled_returns_empty():
    sheets = _sciex_cross_layer_excel_sheets(None)
    assert not sheets

def test_safeguards_in_records():
    config = get_config()
    res = build_sciex_cross_layer_evidence_optional_results(
        config, {"sciex_intact_oxygen_water_state_audit": full()}, []
    )
    sheets = _sciex_cross_layer_excel_sheets(res[OPTIONAL_RESULT_KEY])

    for sheet_name, df in sheets.items():
        if df.empty:
            continue

        # Ensure that exact identity/localization attributes are never marked True
        safeguard_cols = [
            "formal_propagation", "chemical_identity_assigned", "modification_assigned",
            "exact_candidate_identity_confirmed", "exact_isomer_identity_confirmed",
            "exact_nucleotide_localization", "exact_atom_localization",
            "reaction_order_assigned", "applied_to_formal_score",
            "applied_to_ranking", "applied_to_candidate_filtering",
            "applied_to_final_consensus"
        ]

        for col in safeguard_cols:
            if col in df.columns:
                # Value must be False or equivalent (not True)
                assert all(val is not True for val in df[col].values)
