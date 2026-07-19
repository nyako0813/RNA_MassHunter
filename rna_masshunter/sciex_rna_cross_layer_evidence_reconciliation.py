"""Deterministic shadow reconciliation across intact, T1, and P1/AP evidence.

This module consumes existing audit results.  It does not parse raw data, alter an
input result, assign chemistry/localization/order, or propagate into formal output.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, replace
from enum import Enum
from hashlib import sha256
from typing import Any, Iterable, Mapping, Sequence

from rna_masshunter.sciex_t1_fragment_state_series_audit import (
    StateSeriesAuditParameters, build_default_state_delta_definitions,
)

OPTIONAL_RESULT_KEY = "sciex_rna_cross_layer_evidence_reconciliation"
ALGORITHM_VERSION = "sciex-rna-cross-layer-evidence-reconciliation-v1"
DELTA_TOLERANCE_DA = StateSeriesAuditParameters().apex_centroid_neutral_disagreement_da
KNOWN_STATE_DELTAS = tuple(x.target_neutral_delta for x in build_default_state_delta_definitions())

_BLOCK_ORDER = (
    "FULL_LENGTH_RESULT_MISSING", "T1_RESULT_MISSING", "P1AP_MS1_RESULT_MISSING",
    "P1AP_MS2_RESULT_MISSING", "SOURCE_PROVENANCE_MISSING", "RNA_CONTEXT_MISMATCH",
    "DIGEST_CONTEXT_MISMATCH", "ION_MODE_CONTEXT_MISMATCH", "NO_FULL_LENGTH_STATE_PATTERN",
    "NO_T1_FRAGMENT_MATCHES", "NO_T1_STATE_SERIES", "P1AP_STATE_IDENTITY_AMBIGUOUS",
    "P1AP_MS2_NONDISCRIMINATING", "SHARED_SOURCE_DEPENDENCY",
    "MISSING_INDEPENDENT_SUPPORT", "DELTA_PATTERN_PARTIAL_ONLY", "DELTA_PATTERN_CONFLICT",
    "LOCALIZATION_UNSUPPORTED", "EXACT_IDENTITY_UNSUPPORTED",
    "EXACT_ISOMER_IDENTITY_UNSUPPORTED", "REACTION_ORDER_UNSUPPORTED",
    "BACKGROUND_CONTROL_MISSING", "BUFFER_CONTROL_MISSING",
    "BACKGROUND_EXPLANATION_NOT_EXCLUDED", "ALTERNATIVE_CANONICAL_SPECIES_EXPLANATION",
    "INSUFFICIENT_CROSS_LAYER_EVIDENCE",
)


def _blocks(values: Iterable[str]) -> tuple[str, ...]:
    found = {str(x) for x in values if x}
    order = {value: index for index, value in enumerate(_BLOCK_ORDER)}
    return tuple(sorted(found, key=lambda x: (order.get(x, len(order)), x)))


def _id(prefix: str, *values: Any) -> str:
    return prefix + "__" + sha256("|".join(map(str, values)).encode()).hexdigest()[:20].upper()


def _value(obj: Any, *names: str, default: Any = None) -> Any:
    for name in names:
        if isinstance(obj, Mapping) and name in obj:
            return obj[name]
        if obj is not None and hasattr(obj, name):
            return getattr(obj, name)
    return default


def _text(value: Any, default: str = "UNKNOWN") -> str:
    if isinstance(value, Enum): value = value.value
    text = str(value or "").strip()
    return text or default


def _rna_key(value: Any) -> str:
    return "".join(character for character in _text(value, "").upper() if character.isalnum())


def _tuple_float(value: Any) -> tuple[float, ...]:
    if value is None or isinstance(value, (str, bytes)):
        return ()
    try:
        values = tuple(float(x) for x in value)
    except (TypeError, ValueError, OverflowError):
        return ()
    return values


def _normalized(values: Sequence[float]) -> tuple[float, ...]:
    values = tuple(float(x) for x in values)
    return tuple(x - values[0] for x in values) if values else ()


def _enum_values(values: Any) -> tuple[str, ...]:
    return tuple(_text(x) for x in (values or ()))


@dataclass(frozen=True, kw_only=True)
class CrossLayerSafeguards:
    shadow_analysis_only: bool = True
    cross_layer_reconciliation_only: bool = True
    formal_propagation: bool = False
    chemical_identity_assigned: bool = False
    modification_assigned: bool = False
    exact_candidate_identity_confirmed: bool = False
    exact_isomer_identity_confirmed: bool = False
    exact_nucleotide_localization: bool = False
    exact_atom_localization: bool = False
    reaction_order_assigned: bool = False
    applied_to_formal_score: bool = False
    applied_to_ranking: bool = False
    applied_to_candidate_filtering: bool = False
    applied_to_final_consensus: bool = False


@dataclass(frozen=True, kw_only=True)
class CrossLayerEvidenceNode(CrossLayerSafeguards):
    evidence_node_id: str
    layer: str
    evidence_type: str
    source_audit: str
    source_record_id: str
    rna_identity: str
    digest_type: str
    ion_mode: str
    observed_mz_or_mass: tuple[float, ...]
    neutral_delta_pattern: tuple[float, ...]
    candidate_name: str
    candidate_class: str
    fragment_id: str
    fragment_sequence: str
    fragment_start: int | None
    fragment_end: int | None
    identity_level: str
    localization_level: str
    evidence_status: str
    evidence_confidence: str
    ambiguity_status: str
    alternative_explanation_count: int
    block_reasons: tuple[str, ...]
    source_file_id: str
    source_run_id: str
    source_digest: str
    evidence_independence_group: str
    shared_source_dependency: bool


@dataclass(frozen=True, kw_only=True)
class CrossLayerEvidenceEdge(CrossLayerSafeguards):
    evidence_edge_id: str
    source_node_id: str
    target_node_id: str
    edge_type: str
    compatibility_status: str
    compatibility_confidence: str
    shared_data_source: bool
    independent_evidence: bool
    support_direction: str
    conflict_direction: str
    reason: str
    block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class CrossLayerInterpretationHypothesis(CrossLayerSafeguards):
    hypothesis_id: str
    hypothesis_class: str
    description: str
    supporting_node_ids: tuple[str, ...]
    conflicting_node_ids: tuple[str, ...]
    missing_expected_node_types: tuple[str, ...]
    independent_support_group_count: int
    shared_source_support_count: int
    alternative_explanation_count: int
    delta_compatibility_status: str
    identity_compatibility_status: str
    localization_status: str
    evidence_status: str
    evidence_confidence: str
    deterministic_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class CrossLayerSummary(CrossLayerSafeguards):
    observed_full_length_multi_state_pattern: bool
    t1_fragment_matches_present: bool
    t1_state_series_present: bool
    p1ap_ms1_state_relation_present: bool
    p1ap_ms1_state_relation_ambiguous: bool
    p1ap_ms2_class_support_present: bool
    p1ap_ms2_supports_base_class: bool
    p1ap_ms2_supports_canonical_alternative: bool
    p1ap_ms2_resolves_state_interpretation: bool
    cross_layer_delta_compatibility: str
    cross_layer_identity_compatibility: str
    cross_layer_localization_status: str
    best_supported_interpretation_class: str
    alternative_interpretation_classes: tuple[str, ...]
    cross_layer_evidence_status: str
    cross_layer_confidence: str
    evidence_node_count: int
    evidence_edge_count: int
    hypothesis_count: int
    independence_group_count: int
    overall_block_reasons: tuple[str, ...]
    recommended_next_evidence: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class CrossLayerLayerSummary(CrossLayerSafeguards):
    layer: str
    node_count: int
    supported_node_count: int
    ambiguous_node_count: int
    independence_groups: tuple[str, ...]
    layer_status: str
    layer_block_reasons: tuple[str, ...]


@dataclass(frozen=True, kw_only=True)
class RecommendedNextEvidence(CrossLayerSafeguards):
    evidence_type: str
    priority: str
    reason: str
    ambiguity_addressed: str


@dataclass(frozen=True)
class CrossLayerEvidenceAuditResult:
    nodes: tuple[CrossLayerEvidenceNode, ...]
    edges: tuple[CrossLayerEvidenceEdge, ...]
    hypotheses: tuple[CrossLayerInterpretationHypothesis, ...]
    layer_summaries: tuple[CrossLayerLayerSummary, ...]
    consensus: CrossLayerSummary
    next_evidence: tuple[RecommendedNextEvidence, ...]
    algorithm_version: str = ALGORITHM_VERSION
    formal_propagation: bool = False


def _context(runtime_context: Mapping[str, Any] | None) -> dict[str, Any]:
    context = dict(runtime_context or {})
    return {
        "rna": _text(context.get("RNA_Identity", context.get("rna_identity", "UNKNOWN"))),
        "source": _text(context.get("Context_Source", context.get("context_source", "UNKNOWN"))),
        "confidence": _text(context.get("Context_Confidence", context.get("context_confidence", "UNKNOWN"))),
        "background": bool(context.get("Background_Control_Available", context.get("background_control_available", False))),
        "buffer": bool(context.get("Buffer_Control_Available", context.get("buffer_control_available", False))),
    }


def _full_length_info(result: Any) -> dict[str, Any]:
    if result is None: return {}
    if isinstance(result, Sequence) and not isinstance(result, (str, bytes, Mapping)):
        masses = _tuple_float(result); source = {}
    else:
        source = result
        masses = _tuple_float(_value(result, "observed_masses", "member_apex_masses", "masses", "full_length_series"))
        if not masses:
            series = _value(result, "series", default=()) or ()
            if series:
                best = max(series, key=lambda x: (_value(x, "member_count", default=0), _value(x, "mass_span_da", default=0)))
                masses = _tuple_float(_value(best, "member_apex_masses", default=()))
    pattern = _tuple_float(_value(source, "normalized_delta_pattern", "neutral_delta_pattern")) or _normalized(masses)
    return {
        "masses": masses, "pattern": pattern,
        "source_id": _text(_value(source, "source_id", default="FULL_LENGTH_RUNTIME")),
        "rna": _text(_value(source, "rna_identity", "RNA_Identity")),
        "digest": _text(_value(source, "digest_type", default="FULL_LENGTH")),
        "path": _text(_value(source, "input_path", "source_file_id", default="UNKNOWN")),
        "run": _text(_value(source, "source_run_id", "run_id", default=_value(source, "source_id", default="FULL_LENGTH_RUNTIME"))),
    }


def _result_source(result: Any, layer: str) -> dict[str, str]:
    summary = _value(result, "run_summary", default=_value(result, "summary"))
    source_id = _text(_value(summary, "source_id", default=f"{layer}_RUNTIME"))
    path = _text(_value(summary, "input_path", default="UNKNOWN"))
    return {
        "source_id": source_id, "path": path,
        "rna": _text(_value(summary, "rna_identity", default="UNKNOWN")),
        "digest": _text(_value(summary, "digest_type", default=layer)),
        "run": source_id,
    }


def _node(*, layer: str, evidence_type: str, source_audit: str, record_id: str,
          source: Mapping[str, str], context_rna: str, identity_level: str,
          localization_level: str, evidence_status: str, evidence_confidence: str,
          ion_mode: str = "UNKNOWN", observed: Sequence[float] = (), deltas: Sequence[float] = (),
          candidate_name: str = "", candidate_class: str = "", fragment_id: str = "",
          fragment_sequence: str = "", fragment_start: int | None = None,
          fragment_end: int | None = None, ambiguity: str = "UNAMBIGUOUS",
          alternatives: int = 0, blocks: Iterable[str] = (), shared: bool = False) -> CrossLayerEvidenceNode:
    rna = _text(source.get("rna"), context_rna)
    block_list = list(blocks)
    if context_rna != "UNKNOWN" and rna != "UNKNOWN" and _rna_key(rna) != _rna_key(context_rna):
        block_list.append("RNA_CONTEXT_MISMATCH"); evidence_status = "CONTEXT_MISMATCH"; evidence_confidence = "NONE"
    path = _text(source.get("path")); run = _text(source.get("run")); digest = _text(source.get("digest"))
    group = _id("INDEPENDENCE", path, run, digest)
    return CrossLayerEvidenceNode(
        evidence_node_id=_id("XLNODE", layer, evidence_type, record_id, path), layer=layer,
        evidence_type=evidence_type, source_audit=source_audit, source_record_id=record_id,
        rna_identity=rna, digest_type=digest, ion_mode=ion_mode,
        observed_mz_or_mass=tuple(float(x) for x in observed),
        neutral_delta_pattern=tuple(float(x) for x in deltas), candidate_name=candidate_name,
        candidate_class=candidate_class, fragment_id=fragment_id, fragment_sequence=fragment_sequence,
        fragment_start=fragment_start, fragment_end=fragment_end, identity_level=identity_level,
        localization_level=localization_level, evidence_status=evidence_status,
        evidence_confidence=evidence_confidence, ambiguity_status=ambiguity,
        alternative_explanation_count=int(alternatives), block_reasons=_blocks(block_list),
        source_file_id=path, source_run_id=run, source_digest=digest,
        evidence_independence_group=group, shared_source_dependency=shared,
    )


def normalize_cross_layer_evidence_nodes(*, full_length_result: Any | None,
        t1_result: Any | None, p1ap_ms1_result: Any | None, p1ap_ms2_result: Any | None,
        runtime_context: Mapping[str, Any] | None = None) -> list[CrossLayerEvidenceNode]:
    ctx = _context(runtime_context); nodes: list[CrossLayerEvidenceNode] = []
    context_source = {"path": ctx["source"], "run": ctx["source"], "digest": "CONTEXT", "rna": ctx["rna"]}
    nodes.append(_node(layer="CONTEXT", evidence_type="USER_MANIFEST_CONTEXT",
        source_audit="runtime_context", record_id=ctx["source"], source=context_source,
        context_rna=ctx["rna"], identity_level="RNA_CONTEXT", localization_level="RNA_LEVEL",
        evidence_status="USER_CONFIRMED" if ctx["confidence"] == "USER_CONFIRMED" else "CONTEXT_RECORDED",
        evidence_confidence="HIGH" if ctx["confidence"] == "USER_CONFIRMED" else "LOW"))

    full = _full_length_info(full_length_result)
    if full:
        source = {"path": full["path"], "run": full["run"], "digest": full["digest"], "rna": full["rna"]}
        nodes.append(_node(layer="FULL_LENGTH", evidence_type="SOURCE_METADATA", source_audit="full_length_result",
            record_id=full["source_id"], source=source, context_rna=ctx["rna"], identity_level="RNA_CONTEXT",
            localization_level="RNA_LEVEL", evidence_status="SOURCE_PROVENANCE_RECORDED",
            evidence_confidence="HIGH" if full["path"] != "UNKNOWN" else "LOW",
            blocks=() if full["path"] != "UNKNOWN" else ("SOURCE_PROVENANCE_MISSING",)))
        nodes.append(_node(layer="FULL_LENGTH", evidence_type="FULL_LENGTH_STATE_SERIES",
            source_audit="sciex_intact_state_series", record_id=full["source_id"] + "__SERIES",
            source=source, context_rna=ctx["rna"], identity_level="STATE_DELTA_ONLY",
            localization_level="RNA_LEVEL", evidence_status="FULL_LENGTH_MULTI_STATE_MASS_PATTERN_OBSERVED" if len(full["pattern"]) > 1 else "NO_FULL_LENGTH_STATE_PATTERN",
            evidence_confidence="MEDIUM" if len(full["pattern"]) > 1 else "NONE",
            observed=full["masses"], deltas=full["pattern"], ambiguity="CHEMICAL_SEQUENCE_UNASSIGNED",
            blocks=("EXACT_IDENTITY_UNSUPPORTED", "REACTION_ORDER_UNSUPPORTED") if len(full["pattern"]) > 1 else ("NO_FULL_LENGTH_STATE_PATTERN",)))

    if t1_result is not None:
        source = _result_source(t1_result, "T1_DIGEST")
        nodes.append(_node(layer="T1", evidence_type="SOURCE_METADATA", source_audit=type(t1_result).__name__,
            record_id=source["source_id"], source=source, context_rna=ctx["rna"], identity_level="RNA_CONTEXT",
            localization_level="DIGEST_LAYER_ONLY", evidence_status="SOURCE_PROVENANCE_RECORDED",
            evidence_confidence="HIGH" if source["path"] != "UNKNOWN" else "LOW",
            blocks=() if source["path"] != "UNKNOWN" else ("SOURCE_PROVENANCE_MISSING",)))
        for match in sorted(_value(t1_result, "fragment_matches", default=()) or (), key=lambda x: _text(_value(x, "match_id"))):
            ambiguous = _text(_value(match, "fragment_ambiguity_status"))
            nodes.append(_node(layer="T1", evidence_type="T1_FRAGMENT_MATCH", source_audit=type(t1_result).__name__,
                record_id=_text(_value(match, "match_id")), source=source, context_rna=ctx["rna"],
                identity_level="FRAGMENT_SEQUENCE_CANDIDATE", localization_level="FRAGMENT_RANGE_ONLY",
                evidence_status="FRAGMENT_MASS_MATCH", evidence_confidence="MEDIUM" if ambiguous == "UNAMBIGUOUS" else "LOW",
                ion_mode=_text(_value(match, "ion_mode")), observed=(_value(match, "observed_apex_mz", default=0.0),),
                candidate_name=_text(_value(match, "fragment_sequence"), ""), candidate_class="T1_FRAGMENT",
                fragment_id=_text(_value(match, "fragment_id"), ""), fragment_sequence=_text(_value(match, "fragment_sequence"), ""),
                fragment_start=_value(match, "start_position"), fragment_end=_value(match, "end_position"),
                ambiguity=ambiguous, alternatives=max(0, int(_value(match, "candidate_count_for_peak", default=1)) - 1),
                blocks=_value(match, "match_block_reasons", default=())))
        families = tuple(_value(t1_result, "state_families", default=()) or ())
        if families:
            for family in sorted(families, key=lambda x: _text(_value(x, "state_family_id"))):
                nodes.append(_node(layer="T1", evidence_type="T1_FRAGMENT_STATE_FAMILY", source_audit=type(t1_result).__name__,
                    record_id=_text(_value(family, "state_family_id")), source=source, context_rna=ctx["rna"],
                    identity_level="STATE_DELTA_ONLY", localization_level="FRAGMENT_RANGE_ONLY",
                    evidence_status="T1_STATE_SERIES_OBSERVED", evidence_confidence=_text(_value(family, "series_confidence"), "LOW"),
                    ion_mode=_text(_value(family, "ion_mode")), observed=(_value(family, "base_observed_mz", default=0.0),),
                    deltas=_value(family, "observed_neutral_deltas", default=()), fragment_id=_text(_value(family, "fragment_id"), ""),
                    fragment_sequence=_text(_value(family, "fragment_sequence"), ""), fragment_start=_value(family, "start_position"),
                    fragment_end=_value(family, "end_position"), ambiguity=_text(_value(family, "series_ambiguity_status")),
                    blocks=_value(family, "series_block_reasons", default=())))
        else:
            nodes.append(_node(layer="T1", evidence_type="T1_FRAGMENT_STATE_FAMILY", source_audit=type(t1_result).__name__,
                record_id="NO_T1_STATE_FAMILY", source=source, context_rna=ctx["rna"], identity_level="STATE_DELTA_ONLY",
                localization_level="NO_LOCALIZATION", evidence_status="T1_STATE_NOT_OBSERVED", evidence_confidence="LOW",
                ambiguity="NOT_APPLICABLE", blocks=("NO_T1_STATE_SERIES", "LOCALIZATION_UNSUPPORTED")))

    p1_source = _result_source(p1ap_ms1_result, "P1_AP_DIGEST") if p1ap_ms1_result is not None else {}
    if p1ap_ms1_result is not None:
        source = p1_source
        nodes.append(_node(layer="P1AP_MS1", evidence_type="SOURCE_METADATA", source_audit=type(p1ap_ms1_result).__name__,
            record_id=source["source_id"], source=source, context_rna=ctx["rna"], identity_level="RNA_CONTEXT",
            localization_level="DIGEST_LAYER_ONLY", evidence_status="SOURCE_PROVENANCE_RECORDED",
            evidence_confidence="HIGH" if source["path"] != "UNKNOWN" else "LOW",
            blocks=() if source["path"] != "UNKNOWN" else ("SOURCE_PROVENANCE_MISSING",), shared=True))
        for match in sorted(_value(p1ap_ms1_result, "matches", default=()) or (), key=lambda x: _text(_value(x, "match_id"))):
            ambiguity = _text(_value(match, "identity_ambiguity_status"))
            cls = _text(_value(match, "candidate_class"))
            identity = "CANONICAL_NUCLEOSIDE_CLASS" if cls == "NEUTRAL_NUCLEOSIDE" else "MODIFIED_NUCLEOSIDE_MASS_CLASS"
            nodes.append(_node(layer="P1AP_MS1", evidence_type="P1AP_MS1_CANDIDATE_MATCH", source_audit=type(p1ap_ms1_result).__name__,
                record_id=_text(_value(match, "match_id")), source=source, context_rna=ctx["rna"], identity_level=identity,
                localization_level="NUCLEOSIDE_CLASS_ONLY", evidence_status="MASS_COMPATIBLE_ONLY",
                evidence_confidence="MEDIUM" if ambiguity == "UNAMBIGUOUS" else "LOW", ion_mode=_text(_value(match, "ion_mode")),
                observed=(_value(match, "observed_apex_mz", default=0.0),), candidate_name=_text(_value(match, "candidate_name"), ""),
                candidate_class=cls, ambiguity=ambiguity, alternatives=max(0, int(_value(match, "candidate_count_for_peak", default=1)) - 1),
                blocks=_value(match, "match_block_reasons", default=()), shared=True))
        for family in sorted(_value(p1ap_ms1_result, "state_families", default=()) or (), key=lambda x: _text(_value(x, "state_family_id"))):
            alternatives = int(_value(family, "alternative_explanation_count", default=0))
            ambiguity = _text(_value(family, "identity_ambiguity_status"))
            blocks = list(_value(family, "series_block_reasons", default=()))
            if ambiguity != "UNAMBIGUOUS": blocks.append("P1AP_STATE_IDENTITY_AMBIGUOUS")
            nodes.append(_node(layer="P1AP_MS1", evidence_type="P1AP_MS1_STATE_FAMILY", source_audit=type(p1ap_ms1_result).__name__,
                record_id=_text(_value(family, "state_family_id")), source=source, context_rna=ctx["rna"], identity_level="STATE_DELTA_ONLY",
                localization_level="NUCLEOSIDE_CLASS_ONLY", evidence_status="P1AP_STATE_RELATION_OBSERVED_AMBIGUOUS" if ambiguity != "UNAMBIGUOUS" else "P1AP_STATE_RELATION_OBSERVED",
                evidence_confidence=_text(_value(family, "series_confidence"), "LOW"), ion_mode=_text(_value(family, "ion_mode")),
                observed=(_value(family, "base_observed_mz", default=0.0),), deltas=_value(family, "observed_neutral_deltas", default=()),
                candidate_name=_text(_value(family, "base_candidate_name"), ""), candidate_class="NUCLEOSIDE_STATE_FAMILY",
                ambiguity=ambiguity, alternatives=alternatives, blocks=blocks, shared=True))

    if p1ap_ms2_result is not None:
        source = p1_source or {
            "source_id": _text(_value(_value(p1ap_ms2_result, "summary"), "source_id", default="P1AP_MS2_RUNTIME")),
            "path": "UNKNOWN", "rna": ctx["rna"], "digest": "P1_AP_DIGEST", "run": "P1AP_MS2_RUNTIME",
        }
        for summary in sorted(_value(p1ap_ms2_result, "candidate_summary_records", default=()) or (), key=lambda x: _text(_value(x, "candidate_id"))):
            status = _text(_value(summary, "ms2_identity_evidence_status"))
            canonical = status == "MS2_SUPPORTS_CANONICAL_NUCLEOSIDE_CLASS"
            nodes.append(_node(layer="P1AP_MS2", evidence_type="P1AP_MS2_CANDIDATE_CLASS_SUPPORT",
                source_audit=type(p1ap_ms2_result).__name__, record_id=_text(_value(summary, "candidate_id")),
                source=source, context_rna=ctx["rna"], identity_level="CANONICAL_NUCLEOSIDE_CLASS" if canonical else "MODIFIED_NUCLEOSIDE_MASS_CLASS",
                localization_level="NUCLEOSIDE_CLASS_ONLY", evidence_status=status,
                evidence_confidence=_text(_value(summary, "ms2_identity_confidence"), "LOW"),
                candidate_name=_text(_value(summary, "candidate_name"), ""), candidate_class="MS2_CLASS_SUPPORT",
                ambiguity=_text(_value(summary, "identity_ambiguity_after_ms2")),
                alternatives=1 if "UNRESOLVED" in _text(_value(summary, "identity_ambiguity_after_ms2")) else 0,
                blocks=_value(summary, "ms2_block_reasons", default=()), shared=True))
        ms2_summary = _value(p1ap_ms2_result, "summary")
        nodes.append(_node(layer="P1AP_MS2", evidence_type="P1AP_MS2_NONDISCRIMINATING_RESULT",
            source_audit=type(p1ap_ms2_result).__name__, record_id="P1AP_MS2_OVERALL", source=source,
            context_rna=ctx["rna"], identity_level="NUCLEOSIDE_CLASS", localization_level="NO_LOCALIZATION",
            evidence_status=_text(_value(ms2_summary, "overall_evidence_status"), "P1AP_MS2_NONDISCRIMINATING"),
            evidence_confidence=_text(_value(ms2_summary, "overall_confidence"), "LOW"), ambiguity="STATE_INTERPRETATION_UNRESOLVED",
            alternatives=1, blocks=tuple(_value(ms2_summary, "overall_block_reasons", default=())) + ("P1AP_MS2_NONDISCRIMINATING",), shared=True))
    return sorted(nodes, key=lambda x: (x.layer, x.evidence_type, x.source_record_id, x.evidence_node_id))


def _pattern_status(left: Sequence[float], right: Sequence[float], tolerance: float) -> str:
    left, right = tuple(left), tuple(right)
    if len(left) < 2 or len(right) < 2: return "INSUFFICIENT_DELTA_EVIDENCE"
    if len(left) == len(right) and all(abs(a - b) <= tolerance for a, b in zip(left, right)):
        return "EXACT_DELTA_PATTERN_MATCH"
    left_steps = tuple(b - a for a, b in zip(left, left[1:])); right_steps = tuple(b - a for a, b in zip(right, right[1:]))
    if any(abs(a - b) <= tolerance for a in left_steps for b in right_steps): return "SINGLE_STEP_COMPATIBLE"
    if any(abs(a - b) <= tolerance for a in left for b in right if abs(a) > tolerance and abs(b) > tolerance): return "PARTIAL_DELTA_PATTERN_MATCH"
    return "DELTA_PATTERN_CONFLICT"


def _edge(source: CrossLayerEvidenceNode, target: CrossLayerEvidenceNode, edge_type: str,
          status: str, confidence: str, reason: str, blocks: Iterable[str] = (),
          support: str = "SOURCE_TO_TARGET", conflict: str = "NONE") -> CrossLayerEvidenceEdge:
    shared = source.evidence_independence_group == target.evidence_independence_group
    return CrossLayerEvidenceEdge(
        evidence_edge_id=_id("XLEDGE", source.evidence_node_id, target.evidence_node_id, edge_type),
        source_node_id=source.evidence_node_id, target_node_id=target.evidence_node_id,
        edge_type=edge_type, compatibility_status=status, compatibility_confidence=confidence,
        shared_data_source=shared, independent_evidence=not shared,
        support_direction=support, conflict_direction=conflict, reason=reason,
        block_reasons=_blocks(tuple(blocks) + (("SHARED_SOURCE_DEPENDENCY",) if shared else ())),
    )


def build_cross_layer_evidence_edges(nodes: Sequence[CrossLayerEvidenceNode], *,
        reconciliation_config: Mapping[str, Any] | None = None) -> list[CrossLayerEvidenceEdge]:
    config = dict(reconciliation_config or {}); tolerance = float(config.pop("delta_tolerance_da", DELTA_TOLERANCE_DA))
    if config: raise ValueError(f"unsupported reconciliation_config keys: {sorted(config)}")
    usable = [x for x in nodes if "RNA_CONTEXT_MISMATCH" not in x.block_reasons]
    by_type = defaultdict(list); edges = []
    for node in usable: by_type[node.evidence_type].append(node)
    contexts = by_type["USER_MANIFEST_CONTEXT"]
    for context in contexts:
        for source in by_type["SOURCE_METADATA"]:
            edges.append(_edge(context, source, "IDENTITY_CLASS_COMPATIBLE", "RNA_CONTEXT_COMPATIBLE", "HIGH", "User-confirmed RNA context is compatible with source provenance."))
    sources = by_type["SOURCE_METADATA"]
    for i, left in enumerate(sources):
        for right in sources[i + 1:]:
            if left.evidence_independence_group != right.evidence_independence_group:
                edges.append(_edge(left, right, "INDEPENDENT_SUPPORT", "INDEPENDENT_RUN_OR_DIGEST", "HIGH", "Distinct source/run/digest independence groups."))
    full = by_type["FULL_LENGTH_STATE_SERIES"]
    p1families = by_type["P1AP_MS1_STATE_FAMILY"]
    t1families = by_type["T1_FRAGMENT_STATE_FAMILY"]
    for left in full:
        for right in p1families:
            status = _pattern_status(left.neutral_delta_pattern, right.neutral_delta_pattern, tolerance)
            edge_type = "DELTA_PATTERN_COMPATIBLE" if status == "EXACT_DELTA_PATTERN_MATCH" else "PARTIALLY_DELTA_COMPATIBLE" if status in {"SINGLE_STEP_COMPATIBLE", "PARTIAL_DELTA_PATTERN_MATCH"} else "CONFLICTING"
            blocks = ("DELTA_PATTERN_PARTIAL_ONLY",) if edge_type == "PARTIALLY_DELTA_COMPATIBLE" else ("DELTA_PATTERN_CONFLICT",) if edge_type == "CONFLICTING" else ()
            edges.append(_edge(left, right, edge_type, status, "MEDIUM" if edge_type != "CONFLICTING" else "LOW", "Neutral-delta patterns compared using the established state-series tolerance.", blocks, conflict="BIDIRECTIONAL" if edge_type == "CONFLICTING" else "NONE"))
        for right in t1families:
            if right.evidence_status == "T1_STATE_NOT_OBSERVED":
                edges.append(_edge(left, right, "MISSING_EXPECTED_SUPPORT", "T1_STATE_NOT_OBSERVED_NOT_CHEMICAL_CONFLICT", "LOW", "T1 state-series localization support is missing; chemical absence is not inferred.", ("NO_T1_STATE_SERIES", "LOCALIZATION_UNSUPPORTED"), support="NONE"))
            else:
                status = _pattern_status(left.neutral_delta_pattern, right.neutral_delta_pattern, tolerance)
                edges.append(_edge(left, right, "DELTA_PATTERN_COMPATIBLE" if status == "EXACT_DELTA_PATTERN_MATCH" else "PARTIALLY_DELTA_COMPATIBLE", status, "MEDIUM", "Full-length and T1 neutral-delta patterns compared."))
    for match in by_type["T1_FRAGMENT_MATCH"]:
        for family in t1families:
            if family.evidence_status == "T1_STATE_NOT_OBSERVED":
                edges.append(_edge(match, family, "LOCALIZATION_UNSUPPORTED", "UNMODIFIED_FRAGMENT_MATCH_DOES_NOT_LOCALIZE_STATE", "HIGH", "A fragment mass match does not create a state-family localization.", ("LOCALIZATION_UNSUPPORTED",), support="NONE"))
    p1matches = by_type["P1AP_MS1_CANDIDATE_MATCH"]
    ms2supports = by_type["P1AP_MS2_CANDIDATE_CLASS_SUPPORT"]
    for match in p1matches:
        for support in ms2supports:
            if match.candidate_name and match.candidate_name == support.candidate_name:
                edges.append(_edge(match, support, "PRECURSOR_PRODUCT_COMPATIBLE", "MS1_PRECURSOR_MS2_CLASS_COMPATIBLE", "MEDIUM", "MS1 precursor candidate and MS2 class support share a candidate record and raw source."))
                edges.append(_edge(match, support, "SHARED_SOURCE_DEPENDENCE", "NOT_INDEPENDENT_EVIDENCE", "HIGH", "P1/AP MS1 and MS2 derive from the same raw run.", ("SHARED_SOURCE_DEPENDENCY",), support="NONE"))
    for family in p1families:
        for support in ms2supports:
            if support.evidence_status == "MS2_SUPPORTS_CANONICAL_NUCLEOSIDE_CLASS":
                if support.candidate_name == family.candidate_name:
                    edges.append(_edge(family, support, "IDENTITY_CLASS_COMPATIBLE", "BASE_CLASS_SUPPORTED_STATE_UNRESOLVED", "MEDIUM", "MS2 supports the base canonical class without resolving the state relation."))
                else:
                    edges.append(_edge(support, family, "ALTERNATIVE_EXPLANATION", "CANONICAL_SPECIES_ALTERNATIVE", "MEDIUM", "A distinct canonical nucleoside class can explain the state-region mass feature.", ("ALTERNATIVE_CANONICAL_SPECIES_EXPLANATION",), support="ALTERNATIVE_TO_TARGET"))
        for nondiscriminating in by_type["P1AP_MS2_NONDISCRIMINATING_RESULT"]:
            edges.append(_edge(family, nondiscriminating, "NONDISCRIMINATING", "STATE_INTERPRETATION_UNRESOLVED", "HIGH", "MS2 class evidence does not resolve the P1/AP state interpretation.", ("P1AP_MS2_NONDISCRIMINATING",), support="NONE"))
        for t1family in t1families:
            if t1family.evidence_status == "T1_STATE_NOT_OBSERVED":
                edges.append(_edge(family, t1family, "MISSING_EXPECTED_SUPPORT", "T1_LOCALIZATION_SUPPORT_MISSING", "LOW", "P1/AP state relation lacks T1 state-family localization support.", ("NO_T1_STATE_SERIES", "LOCALIZATION_UNSUPPORTED"), support="NONE"))
    return sorted({x.evidence_edge_id: x for x in edges}.values(), key=lambda x: (x.edge_type, x.source_node_id, x.target_node_id))


def _hypothesis(nodes: Sequence[CrossLayerEvidenceNode], hypothesis_class: str, description: str,
                support_types: Sequence[str], *, status: str, confidence: str,
                missing: Sequence[str] = (), blocks: Sequence[str] = (),
                predicate=lambda x: True, delta_status: str = "INSUFFICIENT_DELTA_EVIDENCE",
                identity_status: str = "IDENTITY_UNRESOLVED") -> CrossLayerInterpretationHypothesis:
    supporting = tuple(sorted(x.evidence_node_id for x in nodes if x.evidence_type in support_types and predicate(x)))
    groups = {x.evidence_independence_group for x in nodes if x.evidence_node_id in supporting and x.evidence_type not in {"USER_MANIFEST_CONTEXT", "SOURCE_METADATA"}}
    shared = Counter(x.evidence_independence_group for x in nodes if x.evidence_node_id in supporting)
    alternatives = sum(x.alternative_explanation_count for x in nodes if x.evidence_node_id in supporting)
    missing_types = tuple(sorted(x for x in missing if not any(n.evidence_type == x and n.evidence_status not in {"T1_STATE_NOT_OBSERVED", "CONTEXT_MISMATCH"} for n in nodes)))
    return CrossLayerInterpretationHypothesis(
        hypothesis_id=_id("XLHYP", hypothesis_class), hypothesis_class=hypothesis_class,
        description=description, supporting_node_ids=supporting, conflicting_node_ids=(),
        missing_expected_node_types=missing_types, independent_support_group_count=len(groups),
        shared_source_support_count=sum(count - 1 for count in shared.values() if count > 1),
        alternative_explanation_count=alternatives, delta_compatibility_status=delta_status,
        identity_compatibility_status=identity_status,
        localization_status="FRAGMENT_RANGE_ONLY" if any(x.localization_level == "FRAGMENT_RANGE_ONLY" for x in nodes if x.evidence_node_id in supporting) else "NO_LOCALIZATION",
        evidence_status=status if supporting else "INSUFFICIENT_CROSS_LAYER_EVIDENCE",
        evidence_confidence=confidence if supporting else "NONE",
        deterministic_block_reasons=_blocks(tuple(blocks) + (("MISSING_INDEPENDENT_SUPPORT",) if len(groups) < 2 else ()) + (("INSUFFICIENT_CROSS_LAYER_EVIDENCE",) if not supporting else ())),
    )


def generate_cross_layer_interpretation_hypotheses(nodes: Sequence[CrossLayerEvidenceNode], edges: Sequence[CrossLayerEvidenceEdge]) -> list[CrossLayerInterpretationHypothesis]:
    delta_edges = [x for x in edges if x.edge_type in {"DELTA_PATTERN_COMPATIBLE", "PARTIALLY_DELTA_COMPATIBLE", "CONFLICTING"}]
    delta_status = delta_edges[0].compatibility_status if delta_edges else "INSUFFICIENT_DELTA_EVIDENCE"
    specs = [
        _hypothesis(nodes, "H_FULL_LENGTH_MULTI_STATE_PROCESS", "Multiple full-length mass states are observed; chemistry and order remain unassigned.", ("FULL_LENGTH_STATE_SERIES",), status="MASS_PATTERN_OBSERVED_CHEMISTRY_UNRESOLVED", confidence="MEDIUM", blocks=("EXACT_IDENTITY_UNSUPPORTED", "REACTION_ORDER_UNSUPPORTED"), delta_status=delta_status, identity_status="MASS_PATTERN_ONLY"),
        _hypothesis(nodes, "H_SINGLE_NUCLEOSIDE_PLUS16_RELATION", "A single-nucleoside +16-equivalent mass relation is possible but identity remains ambiguous.", ("P1AP_MS1_STATE_FAMILY",), status="MASS_COMPATIBLE_UNRESOLVED", confidence="LOW", missing=("T1_FRAGMENT_STATE_FAMILY",), blocks=("P1AP_STATE_IDENTITY_AMBIGUOUS", "P1AP_MS2_NONDISCRIMINATING", "LOCALIZATION_UNSUPPORTED", "EXACT_IDENTITY_UNSUPPORTED"), delta_status=delta_status),
        _hypothesis(nodes, "H_CANONICAL_A_AND_G_COEXISTENCE", "Distinct canonical nucleoside classes can coexist and explain the two P1/AP regions.", ("P1AP_MS2_CANDIDATE_CLASS_SUPPORT",), status="SUPPORTED_AS_ALTERNATIVE", confidence="MEDIUM", blocks=("SHARED_SOURCE_DEPENDENCY", "ALTERNATIVE_CANONICAL_SPECIES_EXPLANATION", "EXACT_IDENTITY_UNSUPPORTED"), predicate=lambda x: x.evidence_status == "MS2_SUPPORTS_CANONICAL_NUCLEOSIDE_CLASS", identity_status="CANONICAL_CLASSES_SUPPORTED_EXACT_SPECIES_UNCONFIRMED"),
        _hypothesis(nodes, "H_MODIFIED_NUCLEOSIDE_STATE_RELATION", "Modified nucleoside mass classes remain possible without structure-specific MS2 support.", ("P1AP_MS1_CANDIDATE_MATCH", "P1AP_MS2_CANDIDATE_CLASS_SUPPORT"), status="MASS_COMPATIBLE_PRECURSOR_ONLY", confidence="LOW", blocks=("P1AP_MS2_NONDISCRIMINATING", "EXACT_IDENTITY_UNSUPPORTED", "EXACT_ISOMER_IDENTITY_UNSUPPORTED"), predicate=lambda x: x.identity_level == "MODIFIED_NUCLEOSIDE_MASS_CLASS"),
        _hypothesis(nodes, "H_MULTIPLE_INDEPENDENT_SPECIES", "Multiple independently observed species may contribute across runs and digest layers.", ("FULL_LENGTH_STATE_SERIES", "T1_FRAGMENT_MATCH", "P1AP_MS1_CANDIDATE_MATCH"), status="CROSS_LAYER_SUPPORT_LIMITED", confidence="LOW", blocks=("EXACT_IDENTITY_UNSUPPORTED",), identity_status="MULTIPLE_MASS_OR_CLASS_FEATURES"),
        _hypothesis(nodes, "H_BACKGROUND_OR_ISOTOPE_EXPLANATION", "Background, buffer, or isotope/envelope alternatives have not been excluded.", ("P1AP_MS1_STATE_FAMILY", "P1AP_MS2_NONDISCRIMINATING_RESULT"), status="ALTERNATIVE_NOT_EXCLUDED", confidence="MEDIUM", blocks=("BACKGROUND_CONTROL_MISSING", "BUFFER_CONTROL_MISSING", "BACKGROUND_EXPLANATION_NOT_EXCLUDED")),
        _hypothesis(nodes, "H_T1_STATE_BELOW_DETECTION", "A T1 state may be below detection; non-observation is not chemical absence.", ("T1_FRAGMENT_STATE_FAMILY",), status="POSSIBLE_NOT_CONFIRMED", confidence="LOW", blocks=("NO_T1_STATE_SERIES", "LOCALIZATION_UNSUPPORTED"), predicate=lambda x: x.evidence_status == "T1_STATE_NOT_OBSERVED"),
        _hypothesis(nodes, "H_T1_STATE_DESTABILIZED_OR_NOT_RETAINED", "A state may not be retained through T1 preparation or detection.", ("T1_FRAGMENT_STATE_FAMILY",), status="POSSIBLE_NOT_CONFIRMED", confidence="LOW", blocks=("NO_T1_STATE_SERIES", "LOCALIZATION_UNSUPPORTED"), predicate=lambda x: x.evidence_status == "T1_STATE_NOT_OBSERVED"),
        _hypothesis(nodes, "H_INSUFFICIENT_CROSS_LAYER_EVIDENCE", "Current layers do not uniquely connect mass patterns, identities, and localization.", ("FULL_LENGTH_STATE_SERIES", "T1_FRAGMENT_STATE_FAMILY", "P1AP_MS1_STATE_FAMILY", "P1AP_MS2_NONDISCRIMINATING_RESULT"), status="SUPPORTED_LIMITATION", confidence="HIGH", blocks=("P1AP_MS2_NONDISCRIMINATING", "DELTA_PATTERN_PARTIAL_ONLY", "LOCALIZATION_UNSUPPORTED", "EXACT_IDENTITY_UNSUPPORTED", "REACTION_ORDER_UNSUPPORTED", "INSUFFICIENT_CROSS_LAYER_EVIDENCE"), delta_status=delta_status),
    ]
    conflict_ids = tuple(sorted({node_id for edge in edges if edge.edge_type == "CONFLICTING" for node_id in (edge.source_node_id, edge.target_node_id)}))
    if conflict_ids:
        specs = [replace(x, conflicting_node_ids=conflict_ids,
            evidence_status="CROSS_LAYER_CONFLICTING", evidence_confidence="LOW",
            deterministic_block_reasons=_blocks(x.deterministic_block_reasons + ("DELTA_PATTERN_CONFLICT",)))
            if x.delta_compatibility_status == "DELTA_PATTERN_CONFLICT" else x for x in specs]
    return sorted(specs, key=lambda x: x.hypothesis_class)


def _recommendations(nodes: Sequence[CrossLayerEvidenceNode], context: Mapping[str, Any]) -> tuple[RecommendedNextEvidence, ...]:
    types = {x.evidence_type for x in nodes}; ambiguous_p1 = any(x.evidence_type == "P1AP_MS1_STATE_FAMILY" and x.ambiguity_status != "UNAMBIGUOUS" for x in nodes)
    mass_only = any(x.identity_level == "MODIFIED_NUCLEOSIDE_MASS_CLASS" for x in nodes)
    no_t1 = any(x.evidence_type == "T1_FRAGMENT_STATE_FAMILY" and x.evidence_status == "T1_STATE_NOT_OBSERVED" for x in nodes)
    rows = []
    def add(kind, priority, reason, ambiguity): rows.append(RecommendedNextEvidence(evidence_type=kind, priority=priority, reason=reason, ambiguity_addressed=ambiguity))
    if ambiguous_p1:
        add("TARGETED_MS2_FOR_284_REGION", "HIGH", "Acquire targeted product evidence across the ambiguous state/canonical region.", "P1AP_PLUS16_VS_CANONICAL_ALTERNATIVE")
        add("AUTHENTIC_A_AND_G_STANDARDS", "HIGH", "Compare retention and product-ion behavior with canonical standards.", "CANONICAL_CLASS_VS_PURE_SPECIES")
        add("TARGETED_LC_SEPARATION", "HIGH", "Resolve potentially co-observed canonical and state-related species chromatographically.", "COELUTION_AND_IDENTITY_AMBIGUITY")
    if mass_only: add("AUTHENTIC_MODIFIED_NUCLEOSIDE_STANDARD", "MEDIUM", "Mass-only candidates lack structure-specific product rules.", "MODIFIED_NUCLEOSIDE_EXACT_IDENTITY")
    if not context["background"] or not context["buffer"]: add("BLANK_AND_BUFFER_CONTROL", "HIGH", "Background-like broad traces and exogenous canonical species are not excluded.", "BACKGROUND_OR_BUFFER_ORIGIN")
    if "P1AP_MS1_STATE_FAMILY" in types: add("REPLICATE_P1AP_RUN", "MEDIUM", "An independent P1/AP run would test recurrence without reusing the same MS1/MS2 source.", "SHARED_SOURCE_DEPENDENCY_AND_RECURRENCE")
    if no_t1:
        add("T1_MS2_LOCALIZATION", "HIGH", "Direct fragment product evidence is needed for state localization.", "T1_STATE_LOCALIZATION")
        add("IMPROVED_T1_SENSITIVITY", "MEDIUM", "Test whether the missing T1 state family is below detection.", "T1_STATE_NONOBSERVATION")
    order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
    return tuple(sorted({x.evidence_type: x for x in rows}.values(), key=lambda x: (order[x.priority], x.evidence_type)))


def summarize_cross_layer_evidence(nodes: Sequence[CrossLayerEvidenceNode], edges: Sequence[CrossLayerEvidenceEdge],
        hypotheses: Sequence[CrossLayerInterpretationHypothesis], *, runtime_context: Mapping[str, Any] | None = None) -> CrossLayerSummary:
    ctx = _context(runtime_context)
    full = any(x.evidence_type == "FULL_LENGTH_STATE_SERIES" and len(x.neutral_delta_pattern) > 1 for x in nodes)
    t1matches = any(x.evidence_type == "T1_FRAGMENT_MATCH" for x in nodes)
    t1series = any(x.evidence_type == "T1_FRAGMENT_STATE_FAMILY" and x.evidence_status == "T1_STATE_SERIES_OBSERVED" for x in nodes)
    p1families = [x for x in nodes if x.evidence_type == "P1AP_MS1_STATE_FAMILY"]
    p1relation = bool(p1families); p1ambiguous = any(x.ambiguity_status != "UNAMBIGUOUS" for x in p1families)
    canonical = [x for x in nodes if x.evidence_type == "P1AP_MS2_CANDIDATE_CLASS_SUPPORT" and x.evidence_status == "MS2_SUPPORTS_CANONICAL_NUCLEOSIDE_CLASS"]
    base_support = any(any(f.candidate_name == x.candidate_name for f in p1families) for x in canonical)
    alternative_support = len({x.candidate_name for x in canonical}) >= 2
    delta_edges = [x for x in edges if x.edge_type in {"DELTA_PATTERN_COMPATIBLE", "PARTIALLY_DELTA_COMPATIBLE", "CONFLICTING"}]
    delta = delta_edges[0].compatibility_status if delta_edges else "INSUFFICIENT_DELTA_EVIDENCE"
    context_mismatch = any("RNA_CONTEXT_MISMATCH" in x.block_reasons for x in nodes)
    blocks = []
    missing_layers = []
    if not any(x.layer == "FULL_LENGTH" for x in nodes): missing_layers.append("FULL_LENGTH_RESULT_MISSING")
    if not any(x.layer == "T1" for x in nodes): missing_layers.append("T1_RESULT_MISSING")
    if not any(x.layer == "P1AP_MS1" for x in nodes): missing_layers.append("P1AP_MS1_RESULT_MISSING")
    if not any(x.layer == "P1AP_MS2" for x in nodes): missing_layers.append("P1AP_MS2_RESULT_MISSING")
    blocks.extend(missing_layers)
    if not full: blocks.append("NO_FULL_LENGTH_STATE_PATTERN")
    if not t1matches: blocks.append("NO_T1_FRAGMENT_MATCHES")
    if not t1series: blocks += ["NO_T1_STATE_SERIES", "LOCALIZATION_UNSUPPORTED"]
    if p1ambiguous: blocks.append("P1AP_STATE_IDENTITY_AMBIGUOUS")
    if any(x.evidence_type == "P1AP_MS2_NONDISCRIMINATING_RESULT" for x in nodes): blocks.append("P1AP_MS2_NONDISCRIMINATING")
    if any(x.shared_data_source for x in edges): blocks.append("SHARED_SOURCE_DEPENDENCY")
    if delta in {"SINGLE_STEP_COMPATIBLE", "PARTIAL_DELTA_PATTERN_MATCH"}: blocks.append("DELTA_PATTERN_PARTIAL_ONLY")
    if delta == "DELTA_PATTERN_CONFLICT": blocks.append("DELTA_PATTERN_CONFLICT")
    blocks += ["EXACT_IDENTITY_UNSUPPORTED", "EXACT_ISOMER_IDENTITY_UNSUPPORTED", "REACTION_ORDER_UNSUPPORTED"]
    if not ctx["background"]: blocks.append("BACKGROUND_CONTROL_MISSING")
    if not ctx["buffer"]: blocks.append("BUFFER_CONTROL_MISSING")
    if not ctx["background"] or not ctx["buffer"]: blocks.append("BACKGROUND_EXPLANATION_NOT_EXCLUDED")
    if alternative_support: blocks.append("ALTERNATIVE_CANONICAL_SPECIES_EXPLANATION")
    if context_mismatch: blocks.append("RNA_CONTEXT_MISMATCH")
    if missing_layers or context_mismatch: status, confidence = "INSUFFICIENT_CROSS_LAYER_EVIDENCE", "LOW"
    elif delta == "DELTA_PATTERN_CONFLICT": status, confidence = "CROSS_LAYER_CONFLICTING", "LOW"
    elif p1ambiguous or not t1series: status, confidence = "CROSS_LAYER_SUPPORT_AMBIGUOUS", "LOW"
    elif delta == "EXACT_DELTA_PATTERN_MATCH": status, confidence = "CROSS_LAYER_SUPPORT_MODERATE", "MEDIUM"
    else: status, confidence = "CROSS_LAYER_SUPPORT_LIMITED", "LOW"
    best = "H_CANONICAL_A_AND_G_COEXISTENCE" if alternative_support else "H_FULL_LENGTH_MULTI_STATE_PROCESS" if full else "H_INSUFFICIENT_CROSS_LAYER_EVIDENCE"
    alternatives = tuple(sorted(x.hypothesis_class for x in hypotheses if x.hypothesis_class != best and x.supporting_node_ids))
    recommendations = _recommendations(nodes, ctx)
    groups = {x.evidence_independence_group for x in nodes if x.evidence_type == "SOURCE_METADATA" and "RNA_CONTEXT_MISMATCH" not in x.block_reasons}
    return CrossLayerSummary(
        observed_full_length_multi_state_pattern=full, t1_fragment_matches_present=t1matches,
        t1_state_series_present=t1series, p1ap_ms1_state_relation_present=p1relation,
        p1ap_ms1_state_relation_ambiguous=p1ambiguous, p1ap_ms2_class_support_present=bool(canonical),
        p1ap_ms2_supports_base_class=base_support, p1ap_ms2_supports_canonical_alternative=alternative_support,
        p1ap_ms2_resolves_state_interpretation=False,
        cross_layer_delta_compatibility=delta,
        cross_layer_identity_compatibility="CANONICAL_NUCLEOSIDE_CLASSES_SUPPORTED_STATE_IDENTITY_UNRESOLVED" if canonical else "IDENTITY_UNRESOLVED",
        cross_layer_localization_status="FRAGMENT_RANGE_ONLY_FOR_T1_MATCHES_STATE_LOCALIZATION_UNSUPPORTED" if t1matches else "NO_LOCALIZATION",
        best_supported_interpretation_class=best, alternative_interpretation_classes=alternatives,
        cross_layer_evidence_status=status, cross_layer_confidence=confidence,
        evidence_node_count=len(nodes), evidence_edge_count=len(edges), hypothesis_count=len(hypotheses),
        independence_group_count=len(groups), overall_block_reasons=_blocks(blocks),
        recommended_next_evidence=tuple(x.evidence_type for x in recommendations),
    )


def _layer_summaries(nodes: Sequence[CrossLayerEvidenceNode]) -> tuple[CrossLayerLayerSummary, ...]:
    output = []
    by_layer = defaultdict(list)
    for node in nodes: by_layer[node.layer].append(node)
    for layer, rows in sorted(by_layer.items()):
        ambiguous = sum(x.ambiguity_status not in {"UNAMBIGUOUS", "NOT_APPLICABLE"} for x in rows)
        blocks = _blocks(b for x in rows for b in x.block_reasons)
        output.append(CrossLayerLayerSummary(layer=layer, node_count=len(rows),
            supported_node_count=sum(x.evidence_confidence not in {"NONE", "UNKNOWN"} for x in rows),
            ambiguous_node_count=ambiguous, independence_groups=tuple(sorted({x.evidence_independence_group for x in rows})),
            layer_status="AMBIGUOUS" if ambiguous else "SUPPORTED" if rows else "MISSING", layer_block_reasons=blocks))
    return tuple(output)


def audit_rna_cross_layer_evidence_reconciliation(*, full_length_result: Any | None,
        t1_result: Any | None, p1ap_ms1_result: Any | None, p1ap_ms2_result: Any | None,
        runtime_context: Mapping[str, Any] | None = None,
        reconciliation_config: Mapping[str, Any] | None = None) -> CrossLayerEvidenceAuditResult:
    nodes = tuple(normalize_cross_layer_evidence_nodes(full_length_result=full_length_result,
        t1_result=t1_result, p1ap_ms1_result=p1ap_ms1_result, p1ap_ms2_result=p1ap_ms2_result,
        runtime_context=runtime_context))
    edges = tuple(build_cross_layer_evidence_edges(nodes, reconciliation_config=reconciliation_config))
    hypotheses = tuple(generate_cross_layer_interpretation_hypotheses(nodes, edges))
    consensus = summarize_cross_layer_evidence(nodes, edges, hypotheses, runtime_context=runtime_context)
    next_evidence = _recommendations(nodes, _context(runtime_context))
    return CrossLayerEvidenceAuditResult(nodes, edges, hypotheses, _layer_summaries(nodes), consensus, next_evidence)


def _record(value: Any) -> dict[str, Any]:
    row = asdict(value)
    def normalize(item: Any) -> Any:
        if isinstance(item, Enum): return item.value
        if isinstance(item, dict): return {key: normalize(val) for key, val in item.items()}
        if isinstance(item, (tuple, list)): return [normalize(x) for x in item]
        return item
    return normalize(row)


def audit_optional_result(result: CrossLayerEvidenceAuditResult) -> dict[str, Any]:
    return {
        "node_records": [_record(x) for x in result.nodes],
        "edge_records": [_record(x) for x in result.edges],
        "hypothesis_records": [_record(x) for x in result.hypotheses],
        "layer_summary_records": [_record(x) for x in result.layer_summaries],
        "consensus_records": [_record(result.consensus)],
        "next_evidence_records": [_record(x) for x in result.next_evidence],
    }
