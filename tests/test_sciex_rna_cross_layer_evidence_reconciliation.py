from dataclasses import replace
from types import SimpleNamespace as NS

from rna_masshunter.sciex_rna_cross_layer_evidence_reconciliation import (
    CrossLayerEvidenceNode, audit_optional_result,
    audit_rna_cross_layer_evidence_reconciliation,
    build_cross_layer_evidence_edges, generate_cross_layer_interpretation_hypotheses,
    normalize_cross_layer_evidence_nodes,
)

CTX={"RNA_Identity":"TRNA_TEST","Context_Source":"USER_PROVIDED_RUNTIME_MANIFEST","Context_Confidence":"USER_CONFIRMED"}


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


def audit(**overrides):
    values={"full_length_result":full(),"t1_result":t1(),"p1ap_ms1_result":p1(),"p1ap_ms2_result":ms2(),"runtime_context":CTX}
    values.update(overrides); return audit_rna_cross_layer_evidence_reconciliation(**values)


def test_full_length_only_never_confirms_cross_layer_identity():
    result=audit(t1_result=None,p1ap_ms1_result=None,p1ap_ms2_result=None)
    assert result.consensus.observed_full_length_multi_state_pattern
    assert result.consensus.cross_layer_evidence_status!="CROSS_LAYER_SUPPORT_STRONG"
    assert not result.consensus.chemical_identity_assigned and not result.consensus.reaction_order_assigned


def test_full_length_and_t1_exact_compatible_series():
    result=audit(full_length_result=full((0,16)),t1_result=t1(families=(t1family(),)),p1ap_ms1_result=None,p1ap_ms2_result=None)
    assert any(x.compatibility_status=="EXACT_DELTA_PATTERN_MATCH" for x in result.edges)


def test_t1_state_absent_is_missing_support_not_conflict():
    result=audit()
    edges=[x for x in result.edges if x.target_node_id in {n.evidence_node_id for n in result.nodes if n.evidence_status=="T1_STATE_NOT_OBSERVED"}]
    assert any(x.edge_type=="MISSING_EXPECTED_SUPPORT" for x in edges)
    assert all(x.compatibility_status!="T1_STATE_CHEMICALLY_ABSENT" for x in edges)


def test_p1ap_plus16_family_preserves_ambiguity():
    node=next(x for x in audit().nodes if x.evidence_type=="P1AP_MS1_STATE_FAMILY")
    assert node.ambiguity_status=="IDENTITY_AMBIGUOUS"
    assert "P1AP_STATE_IDENTITY_AMBIGUOUS" in node.block_reasons


def test_ms2_canonical_alternative_does_not_resolve_state():
    result=audit()
    assert result.consensus.p1ap_ms2_supports_base_class
    assert result.consensus.p1ap_ms2_supports_canonical_alternative
    assert not result.consensus.p1ap_ms2_resolves_state_interpretation
    assert any(x.edge_type=="ALTERNATIVE_EXPLANATION" for x in result.edges)


def test_p1ap_ms1_ms2_shared_source_is_not_independent():
    result=audit()
    shared=[x for x in result.edges if x.edge_type=="SHARED_SOURCE_DEPENDENCE"]
    assert shared and all(x.shared_data_source and not x.independent_evidence for x in shared)


def test_distinct_full_t1_p1_sources_make_three_independence_groups():
    assert audit().consensus.independence_group_count==3


def test_exact_delta_pattern_status():
    result=audit(full_length_result=full((0,16)),p1ap_ms1_result=p1(families=(p1family((0,16)),)))
    assert result.consensus.cross_layer_delta_compatibility=="EXACT_DELTA_PATTERN_MATCH"


def test_partial_delta_uses_single_step_compatibility():
    assert audit().consensus.cross_layer_delta_compatibility=="SINGLE_STEP_COMPATIBLE"


def test_incompatible_delta_patterns_report_conflict():
    result=audit(full_length_result=full((0,10)),p1ap_ms1_result=p1(families=(p1family((0,16)),)))
    assert result.consensus.cross_layer_delta_compatibility=="DELTA_PATTERN_CONFLICT"
    assert result.consensus.cross_layer_evidence_status=="CROSS_LAYER_CONFLICTING"
    assert any(x.conflicting_node_ids for x in result.hypotheses if x.delta_compatibility_status=="DELTA_PATTERN_CONFLICT")


def test_no_t1_state_family_never_localizes_state():
    result=audit()
    absent=next(x for x in result.nodes if x.evidence_status=="T1_STATE_NOT_OBSERVED")
    assert absent.localization_level=="NO_LOCALIZATION"
    assert not result.consensus.exact_nucleotide_localization


def test_canonical_coexistence_hypothesis_is_retained_as_alternative():
    hyp=next(x for x in audit().hypotheses if x.hypothesis_class=="H_CANONICAL_A_AND_G_COEXISTENCE")
    assert hyp.evidence_status=="SUPPORTED_AS_ALTERNATIVE" and len(hyp.supporting_node_ids)==2
    assert not hyp.chemical_identity_assigned


def test_modified_state_hypothesis_is_not_confirmed():
    hyp=next(x for x in audit().hypotheses if x.hypothesis_class=="H_MODIFIED_NUCLEOSIDE_STATE_RELATION")
    assert hyp.evidence_status=="MASS_COMPATIBLE_PRECURSOR_ONLY"
    assert "EXACT_IDENTITY_UNSUPPORTED" in hyp.deterministic_block_reasons


def test_missing_background_controls_create_blocks_and_recommendation():
    result=audit()
    assert "BACKGROUND_CONTROL_MISSING" in result.consensus.overall_block_reasons
    assert "BUFFER_CONTROL_MISSING" in result.consensus.overall_block_reasons
    assert result.next_evidence[0].priority=="HIGH"
    assert any(x.evidence_type=="BLANK_AND_BUFFER_CONTROL" for x in result.next_evidence)


def test_missing_layers_return_safe_summary():
    result=audit_rna_cross_layer_evidence_reconciliation(full_length_result=None,t1_result=None,
        p1ap_ms1_result=None,p1ap_ms2_result=None,runtime_context=CTX)
    assert result.consensus.evidence_node_count==1
    assert result.consensus.cross_layer_evidence_status=="INSUFFICIENT_CROSS_LAYER_EVIDENCE"
    assert "FULL_LENGTH_RESULT_MISSING" in result.consensus.overall_block_reasons
    assert "P1AP_MS2_RESULT_MISSING" in result.consensus.overall_block_reasons
    assert not result.consensus.chemical_identity_assigned



def test_equivalent_rna_identity_formatting_is_not_a_context_mismatch():
    context={**CTX,"RNA_Identity":"tRNA^Test"}
    result=audit(runtime_context=context)
    assert not any("RNA_CONTEXT_MISMATCH" in x.block_reasons for x in result.nodes)

def test_rna_context_mismatch_is_blocked_and_not_integrated():
    result=audit(t1_result=t1(rna="OTHER_RNA"))
    mismatched=[x for x in result.nodes if "RNA_CONTEXT_MISMATCH" in x.block_reasons]
    assert mismatched and all(x.evidence_status=="CONTEXT_MISMATCH" for x in mismatched)
    ids={x.evidence_node_id for x in mismatched}
    assert not any(x.source_node_id in ids or x.target_node_id in ids for x in result.edges)


def test_deterministic_under_input_record_reordering():
    left=audit()
    right=audit(t1_result=t1(matches=(fragment("M2",start=4),fragment())),
        p1ap_ms1_result=p1(matches=tuple(reversed(p1().matches))),p1ap_ms2_result=ms2(tuple(reversed(ms2().candidate_summary_records))))
    # Compare a run with the same two T1 matches in both orders.
    base=audit(t1_result=t1(matches=(fragment(),fragment("M2",start=4))))
    assert [x.evidence_node_id for x in base.nodes]==[x.evidence_node_id for x in right.nodes]
    assert [x.evidence_edge_id for x in base.edges]==[x.evidence_edge_id for x in right.edges]
    assert [x.hypothesis_class for x in left.hypotheses]==[x.hypothesis_class for x in right.hypotheses]


def test_optional_records_preserve_formal_nonpropagation():
    result=audit(); payload=audit_optional_result(result)
    for rows in payload.values():
        for row in rows:
            assert row["formal_propagation"] is False
            assert row["chemical_identity_assigned"] is False
            assert row["applied_to_formal_score"] is False
            assert row["applied_to_final_consensus"] is False
