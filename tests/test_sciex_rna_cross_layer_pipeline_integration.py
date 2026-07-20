import pytest
from types import SimpleNamespace as NS
from copy import deepcopy

from rna_masshunter.models import RunConfig
from main import build_sciex_cross_layer_evidence_optional_results
from rna_masshunter.sciex_rna_cross_layer_evidence_reconciliation import OPTIONAL_RESULT_KEY

# Fixtures borrowed from unit tests
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

def t1family(pattern=(0,16),family_id="TF1"):
    return NS(state_family_id=family_id,fragment_id="F1",fragment_sequence="ACG",start_position=1,end_position=3,
        localization_level="FRAGMENT_RANGE_ONLY",charge=1,ion_mode="NEGATIVE_DEPROTONATED",base_observed_mz=500,
        observed_neutral_deltas=pattern,series_ambiguity_status="UNAMBIGUOUS",series_confidence="HIGH",series_block_reasons=())

def p1match(name="A",candidate_id=None,cls="NEUTRAL_NUCLEOSIDE",ambiguity="UNAMBIGUOUS"):
    candidate_id=candidate_id or name
    return NS(match_id=f"PM_{candidate_id}",candidate_id=candidate_id,candidate_name=name,
        candidate_class=cls,identity_ambiguity_status=ambiguity,ion_mode="POSITIVE",observed_apex_mz=268,
        candidate_count_for_peak=1,match_block_reasons=())

def p1family(pattern=(0,16),ambiguity="IDENTITY_AMBIGUOUS"):
    return NS(state_family_id="PF1",base_candidate_id="A",base_candidate_name="A",ion_mode="POSITIVE",
        base_observed_mz=268,observed_neutral_deltas=pattern,identity_ambiguity_status=ambiguity,
        series_confidence="LOW",alternative_explanation_count=3,series_block_reasons=())

def p1(*,families=None,matches=None,rna="TRNA_TEST"):
    families=(p1family(),) if families is None else tuple(families)
    matches=(p1match(),p1match("massX","X","MASS_ONLY_MODIFIED_NUCLEOSIDE","IDENTITY_AMBIGUOUS")) if matches is None else tuple(matches)
    return NS(run_summary=source("P1_AP_DIGEST","p1.mzML",rna),matches=matches,state_families=families,
        summary=NS(state_family_count=len(families)))

def ms2summary(cid,name,status):
    return NS(candidate_id=cid,candidate_name=name,ms2_identity_evidence_status=status,
        ms2_identity_confidence="MEDIUM" if "SUPPORTS_CANONICAL" in status else "LOW",
        identity_ambiguity_after_ms2="CLASS_SUPPORTED_EXACT_IDENTITY_UNCONFIRMED",ms2_block_reasons=())

def ms2(rows=None):
    rows=(ms2summary("A","A","MS2_SUPPORTS_CANONICAL_NUCLEOSIDE_CLASS"),
          ms2summary("G","G","MS2_SUPPORTS_CANONICAL_NUCLEOSIDE_CLASS"),
          ms2summary("X","massX","MS2_PRECURSOR_COMPATIBLE_ONLY")) if rows is None else tuple(rows)
    return NS(candidate_summary_records=rows,summary=NS(source_id="P1_AP_DIGEST",
        overall_evidence_status="P1AP_MS2_NUCLEOSIDE_CLASS_SUPPORT_NONDISCRIMINATING",
        overall_confidence="LOW",overall_block_reasons=("MS2_NONDISCRIMINATING",)))


def get_config(enabled=True):
    return RunConfig(
        sequence={"name": "TRNA_TEST"},
        sciex_profile={
            "cross_layer_evidence_reconciliation": {"enabled": enabled}
        },
        input={},
    )


def test_integration_enabled_four_layers():
    config = get_config()
    optional_results = {
        "sciex_intact_oxygen_water_state_audit": full(),
        "sciex_t1_fragment_state_series_audit": t1(),
        "sciex_p1ap_nucleoside_state_audit": p1(),
        "sciex_p1ap_nucleoside_ms2_identity_audit": ms2(),
    }
    warnings = []
    res = build_sciex_cross_layer_evidence_optional_results(config, optional_results, warnings)
    assert OPTIONAL_RESULT_KEY in res, f"Failed with warnings: {warnings}"
    assert not warnings
    consensus = res[OPTIONAL_RESULT_KEY].consensus
    assert not getattr(consensus, "formal_propagation", False)
    assert not getattr(consensus, "chemical_identity_assigned", False)

def test_integration_disabled():
    config = get_config(enabled=False)
    optional_results = {
        "sciex_intact_oxygen_water_state_audit": full(),
    }
    warnings = []
    res = build_sciex_cross_layer_evidence_optional_results(config, optional_results, warnings)
    assert not res
    assert not warnings

def test_integration_section_missing():
    config = RunConfig(sequence={"name": "TRNA_TEST"}, sciex_profile={}, input={})
    optional_results = {"sciex_intact_oxygen_water_state_audit": full()}
    res = build_sciex_cross_layer_evidence_optional_results(config, optional_results, [])
    assert not res

def test_integration_missing_layers():
    config = get_config()
    # Missing everything
    res = build_sciex_cross_layer_evidence_optional_results(config, {}, [])
    assert not res

    # Only full-length
    warnings = []
    res2 = build_sciex_cross_layer_evidence_optional_results(config, {"sciex_intact_oxygen_water_state_audit": full()}, warnings)
    assert OPTIONAL_RESULT_KEY in res2, f"Failed with warnings: {warnings}"

def test_malformed_layer_object():
    config = get_config()
    warnings = []
    # Pass string instead of expected objects
    res = build_sciex_cross_layer_evidence_optional_results(
        config, {"sciex_intact_oxygen_water_state_audit": "malformed_string"}, warnings
    )
    # The pipeline should not crash. It returns a result, but full_length layer won't be parsed correctly
    assert OPTIONAL_RESULT_KEY in res
    nodes = res[OPTIONAL_RESULT_KEY].nodes
    assert any(n.evidence_status == "NO_FULL_LENGTH_STATE_PATTERN" for n in nodes)

def test_raw_decode_not_called(monkeypatch):
    import rna_masshunter.sciex_profile_parser
    import rna_masshunter.sciex_intact_peak_detection

    def raise_error(*args, **kwargs):
        raise RuntimeError("Raw decode should not be called!")

    monkeypatch.setattr(rna_masshunter.sciex_profile_parser, "parse_sciex_profile", raise_error)
    monkeypatch.setattr(rna_masshunter.sciex_intact_peak_detection, "detect_sciex_intact_peaks", raise_error, raising=False)

    config = get_config()
    optional_results = {
        "sciex_intact_oxygen_water_state_audit": full(),
    }
    warnings = []
    res = build_sciex_cross_layer_evidence_optional_results(config, optional_results, warnings)
    assert OPTIONAL_RESULT_KEY in res

def test_ab_formal_non_propagation():
    # Comparing workflow output with and without cross layer
    config_a = get_config(enabled=False)
    config_b = get_config(enabled=True)

    opt_a = {"sciex_intact_oxygen_water_state_audit": full()}
    opt_b = {"sciex_intact_oxygen_water_state_audit": full()}

    res_a = build_sciex_cross_layer_evidence_optional_results(config_a, opt_a, [])
    warnings_b = []
    res_b = build_sciex_cross_layer_evidence_optional_results(config_b, opt_b, warnings_b)

    assert OPTIONAL_RESULT_KEY not in res_a
    assert OPTIONAL_RESULT_KEY in res_b, f"Failed with warnings: {warnings_b}"

    # Assert formal objects in opt_a and opt_b are identical
    assert opt_a == opt_b

def test_integration_incompatible_rna_context():
    config = get_config(enabled=True)

    # Provide a full-length layer with a mismatched RNA identity
    opt = {
        "sciex_intact_oxygen_water_state_audit": full(rna="COMPLETELY_DIFFERENT_RNA"),
    }

    warnings = []
    res = build_sciex_cross_layer_evidence_optional_results(config, opt, warnings)
    assert OPTIONAL_RESULT_KEY in res

    nodes = res[OPTIONAL_RESULT_KEY].nodes
    full_length_node = next(n for n in nodes if n.layer == "FULL_LENGTH" and n.evidence_type == "SOURCE_METADATA")

    # It should be flagged as CONTEXT_MISMATCH due to internal reconciliation logic
    assert full_length_node.evidence_status == "CONTEXT_MISMATCH"
    assert "RNA_CONTEXT_MISMATCH" in full_length_node.block_reasons
