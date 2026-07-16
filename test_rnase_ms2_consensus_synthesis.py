from copy import deepcopy

from rna_masshunter.audit_policy import AuditPolicy, included_sheet_names
from rna_masshunter.rnase_ms2_consensus_synthesis import build_rnase_ms2_consensus_synthesis


def standard(mod="M1", parent="F1", position=3, trna=12, identity="FRAGMENT_SUPPORTED",
             localization="LOCALIZED", structure="UNRESOLVED", ambiguity="NONE"):
    return {
        "Modification_ID": mod, "Parent_Fragment_ID": parent,
        "Candidate_Position_In_Parent": position, "Candidate_tRNA_Position": trna,
        "Modification_Identity_Status": identity, "Localization_Status": localization,
        "Structure_Status": structure, "Ambiguity_Status": ambiguity,
    }


def composite(candidate="C1", structure_id="S1", identity="PROVISIONAL_CANDIDATE_SUPPORT",
              localization="POSITION_COMPATIBLE", backbone="NOT_EVALUATED",
              structure="UNRESOLVED", ambiguity="NONE"):
    return {
        "Candidate_ID": candidate, "Complete_Structure_ID": structure_id,
        "Composite_Identity_Status": identity,
        "Composite_Localization_Status": localization,
        "Composite_Backbone_Status": backbone,
        "Composite_Structure_Status": structure,
        "Composite_Ambiguity_Status": ambiguity,
    }


def link(status="EXACT_MATCH", candidate="C1", structure="S1", cardinality="ONE_TO_ONE",
         mod="M1", parent="F1", position=3, trna=12):
    return {
        "Modification_ID": mod, "Parent_Fragment_ID": parent,
        "Candidate_Position_In_Parent": position, "Candidate_tRNA_Position": trna,
        "Candidate_ID": candidate, "Complete_Structure_ID": structure,
        "Crosswalk_Status": status, "Crosswalk_Cardinality": cardinality,
        "Position_Match_Status": "MATCH" if status != "POSITION_CONFLICT" else "CONFLICT",
        "Parent_Fragment_Match_Status": "CONFLICT" if status == "PARENT_FRAGMENT_CONFLICT" else "MATCH",
    }


def build(standards=(), composites=(), links=()):
    return build_rnase_ms2_consensus_synthesis(list(standards), list(composites), list(links))


def row(result, mod="M1"):
    return next(item for item in result.evidence_rows if item["Modification_ID"] == mod)


def test_exact_match_fragment_support_gives_supported_identity():
    item = row(build([standard()], [composite()], [link()]))
    assert item["Standard_Composite_Consistency_Status"] == "CONSISTENT"
    assert item["Modification_Identity_Consensus"] == "SUPPORTED"
    assert item["Localization_Consensus"] == "LOCALIZED"


def test_precursor_only_with_composite_unsupported_is_provisional_at_most():
    item = row(build(
        [standard(identity="PRECURSOR_COMPATIBLE", localization="UNRESOLVED")],
        [composite(identity="UNSUPPORTED", localization="UNRESOLVED")], [link()],
    ))
    assert item["Modification_Identity_Consensus"] == "PROVISIONAL"
    assert item["Structure_Consensus"] != "SUPPORTED"


def test_many_to_many_caps_consensus_and_is_high_priority():
    composites = [composite(), composite("C2", "S2")]
    links = [link(cardinality="MANY_TO_MANY"), link(candidate="C2", structure="S2", cardinality="MANY_TO_MANY")]
    item = row(build([standard()], composites, links))
    assert item["Modification_Identity_Consensus"] == "AMBIGUOUS"
    assert item["Consensus_Ambiguity_Status"] == "CROSSWALK_MULTIPLICITY"
    assert item["Review_Priority"] == "HIGH"


def test_position_identity_and_parent_conflicts_remain_explicit():
    cases = [
        ("POSITION_CONFLICT", "POSITION_CONFLICT"),
        ("MODIFICATION_IDENTITY_CONFLICT", "IDENTITY_CONFLICT"),
        ("PARENT_FRAGMENT_CONFLICT", "PARENT_FRAGMENT_CONFLICT"),
    ]
    for crosswalk_status, expected in cases:
        item = row(build([standard()], [composite()], [link(crosswalk_status)]))
        assert item["Standard_Composite_Consistency_Status"] == expected
        assert item["Consensus_Ambiguity_Status"] in {"EVIDENCE_CONFLICT", "MULTIPLE"}
        assert item["Review_Priority"] == "HIGH"


def test_mass_only_is_not_identity_consensus_support():
    item = row(build([standard()], [composite()], [link("MASS_EQUIVALENT_ONLY", cardinality="ONE_TO_ONE")]))
    assert item["Exact_Crosswalk_Count"] == 0
    assert item["Standard_Composite_Consistency_Status"] == "INSUFFICIENT_PROVENANCE"
    assert item["Modification_Identity_Consensus"] != "SUPPORTED"


def test_composite_adds_backbone_context_independently():
    item = row(build([standard()], [composite(backbone="BACKBONE_COMPATIBLE")], [link()]))
    assert item["Standard_Composite_Consistency_Status"] == "COMPOSITE_ADDS_CONTEXT"
    assert item["Backbone_Consensus"] == "BOND_COMPATIBLE"


def test_composite_structure_unresolved_never_becomes_supported():
    item = row(build([standard(structure="SUPPORTED")], [composite(structure="UNRESOLVED")], [link()]))
    assert item["Structure_Consensus"] == "UNRESOLVED"


def test_structure_ambiguity_is_high_priority():
    item = row(build([standard()], [composite(structure="AMBIGUOUS")], [link()]))
    assert item["Structure_Consensus"] == "AMBIGUOUS"
    assert item["Review_Priority"] == "HIGH"


def test_standard_only_composite_only_and_not_linked():
    standard_only = row(build([standard()], [], []))
    assert standard_only["Consensus_Ambiguity_Status"] == "STANDARD_ONLY"
    assert standard_only["Standard_Composite_Consistency_Status"] == "NOT_LINKED"
    composite_only = build([], [composite()], []).evidence_rows[0]
    assert composite_only["Consensus_Ambiguity_Status"] == "COMPOSITE_ONLY"
    not_linked = row(build([standard()], [composite()], [link("NOT_MAPPABLE", cardinality="")]))
    assert not_linked["Standard_Composite_Consistency_Status"] == "NOT_LINKED"


def test_input_order_is_deterministic():
    standards = [standard("M2"), standard("M1")]
    composites = [composite("C2", "S2"), composite()]
    links = [link(mod="M2", candidate="C2", structure="S2"), link()]
    first = build(standards, composites, links)
    second = build(list(reversed(deepcopy(standards))), list(reversed(deepcopy(composites))), list(reversed(deepcopy(links))))
    assert first == second


def test_all_formal_flags_false_and_inputs_unchanged():
    standards=[standard()]; composites=[composite()]; links=[link()]
    before=deepcopy((standards,composites,links)); result=build(standards,composites,links)
    for rows in (result.summary_rows,result.evidence_rows):
        for item in rows:
            assert item["Applied_To_Formal_Result"] is False
            assert item["Formal_Change_Ready"] is False
            assert item["Formal_Result_Changed"] is False
    assert (standards,composites,links)==before


def test_sheet_inclusion():
    names=["RNase_MS2_Consensus_Summary","RNase_MS2_Consensus_Evidence"]
    standard_names,_=included_sheet_names(names,AuditPolicy.from_level("standard"))
    audit_names,_=included_sheet_names(names,AuditPolicy.from_level("audit"))
    full_names,_=included_sheet_names(names,AuditPolicy.from_level("full"))
    assert standard_names==[]
    assert audit_names==[names[0]]
    assert full_names==names
