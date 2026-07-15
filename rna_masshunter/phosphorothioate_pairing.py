"""Generic normal-phosphate/phosphorothioate pair construction."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping
import hashlib
import yaml
from rna_masshunter.backbone_state import BackboneTransformation, normal_bond
from rna_masshunter.cleavage_site_discovery import CleavageBondCandidate, discover_candidate_cleavage_bonds
from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.modification_composer import apply_transform_ids, compose_modifications
from rna_masshunter.structure_fragment import CompleteStructureState, StructureFragment, extract_fragment_from_structure

@dataclass(frozen=True)
class PositionStateSpec:
    position: int
    parent_base: str
    transform_ids: tuple[str, ...]

@dataclass(frozen=True)
class PTPairSpec:
    candidate_id: str
    hypothesis_id: str
    search_mode: str
    sequence_id: str
    sequence: str
    enzyme: str
    bond_id: str
    position_states: tuple[PositionStateSpec, ...]
    fragment_start: int
    fragment_end: int
    terminal_form: str = "default"

@dataclass(frozen=True)
class PTPair:
    spec: PTPairSpec
    cleavage_candidate: CleavageBondCandidate
    normal_structure: CompleteStructureState
    modified_structure: CompleteStructureState
    normal_fragment: StructureFragment
    modified_fragment: StructureFragment
    shared_modification_composition: ElementalComposition
    composition_delta: ElementalComposition
    expected_backbone_delta: ElementalComposition
    delta_consistency_error: float
    block_status: str
    block_evidence_status: str

@dataclass(frozen=True)
class PTHypothesisLoadResult:
    schema_version: int
    enabled: bool
    specs: tuple[PTPairSpec, ...]
    invalid_rows: tuple[dict[str, Any], ...]


def _invalid(candidate_id: str, reason: str, detail: str = "") -> dict[str, Any]:
    return {"Candidate_ID": candidate_id, "Valid": False, "Invalid_Reason": reason,
        "Invalid_Detail": detail, "Applied_To_Formal_Result": False,
        "Formal_Change_Ready": False, "Formal_Result_Changed": False}


def composition_difference(left: ElementalComposition, right: ElementalComposition) -> ElementalComposition:
    values = left.to_dict()
    for element, count in right.to_dict().items():
        values[element] = values.get(element, 0) - count
    return ElementalComposition.delta(values)


def load_pt_pair_hypotheses(path: str | Path, *, sequence: str, sequence_id: str,
    organism: str, rule_set: str, transformations: list[Any]) -> PTHypothesisLoadResult:
    source = Path(path)
    if not source.exists():
        return PTHypothesisLoadResult(1, False, (), (_invalid("", "fixture_not_found", str(source)),))
    raw = yaml.safe_load(source.read_text(encoding="utf-8")) or {}
    schema_version = int(raw.get("schema_version", 1) or 1); enabled = bool(raw.get("enabled", True))
    seq = str(sequence or "").upper().replace("T", "U")
    required = raw.get("target") or {}
    actual = {"name": sequence_id, "length": len(seq),
        "sequence_sha256": hashlib.sha256(seq.encode("utf-8")).hexdigest() if seq else "",
        "organism": organism, "rule_set": rule_set}
    failures = [f"{key}:{actual.get(key, '')}!={value}" for key, value in required.items()
        if str(actual.get(key, "")) != str(value)]
    if failures:
        return PTHypothesisLoadResult(schema_version, enabled, (),
            (_invalid("", "target_identity_mismatch", ";".join(failures)),))
    transform_map = {item.id: item for item in transformations}; specs = []; invalid = []
    for index, row in enumerate(raw.get("hypotheses") or (), 1):
        hid = str(row.get("hypothesis_id") or f"PT_HYP_{index:04d}"); reasons = []; details = []
        context = row.get("enzyme_context") or {}; enzyme = str(context.get("enzyme") or "")
        backbone = list(row.get("backbone") or ()); bond_id = str((backbone[0] if backbone else {}).get("bond_id") or "")
        try: left, right = (int(x) for x in bond_id.split("_", 1))
        except (ValueError, TypeError): left = right = -1; reasons.append("invalid_bond_id")
        candidates = {c.bond_id: c for c in discover_candidate_cleavage_bonds(seq, sequence_id, enzyme)} if enzyme else {}
        candidate = candidates.get(bond_id)
        if candidate is None: reasons.append("not_normal_cleavage_site"); details.append(bond_id)
        expected = context.get("expected_cleavage_position")
        if expected not in (None, "") and int(expected) != left: reasons.append("cleavage_position_mismatch")
        positions = []
        for pos_raw, payload in (row.get("positions") or {}).items():
            pos = int(pos_raw); parent = str((payload or {}).get("parent_base") or "").upper().replace("T", "U")
            ids = tuple(str(x) for x in ((payload or {}).get("transformations") or ()))
            if pos < 1 or pos > len(seq) or seq[pos - 1] != parent:
                reasons.append("parent_base_mismatch"); details.append(f"{pos}:{parent}")
            if any(x not in transform_map for x in ids): reasons.append("unknown_transform_id")
            elif any(parent not in transform_map[x].parent_bases for x in ids): reasons.append("transform_parent_base_incompatible")
            positions.append(PositionStateSpec(pos, parent, ids))
        if reasons:
            invalid.append(_invalid(hid, ";".join(dict.fromkeys(reasons)), ";".join(details))); continue
        specs.append(PTPairSpec(hid, hid, "hypothesis_driven", sequence_id, seq, enzyme, bond_id,
            tuple(positions), candidate.fragment_start, candidate.fragment_end,
            str(row.get("terminal_form") or "default")))
    return PTHypothesisLoadResult(schema_version, enabled, tuple(specs) if enabled else (), tuple(invalid))


def build_pt_pair(spec: PTPairSpec, transformations: list[Any], slot_schema_path: str | Path,
    backbone_transform: BackboneTransformation) -> PTPair:
    candidates = {c.bond_id: c for c in discover_candidate_cleavage_bonds(spec.sequence, spec.sequence_id, spec.enzyme)}
    cleavage = candidates.get(spec.bond_id)
    if cleavage is None: raise ValueError("not_normal_cleavage_site")
    left, right = (int(x) for x in spec.bond_id.split("_", 1)); states = {}
    shared = ElementalComposition.delta()
    for item in spec.position_states:
        state, result, _ = apply_transform_ids(item.parent_base, item.position, item.transform_ids,
            transformations, slot_schema_path)
        if not result.valid: raise ValueError(result.reason_code)
        if item.transform_ids:
            states[item.position] = state; shared = shared + state.elemental_composition_delta
    normal = normal_bond(left, right); modified = normal.apply(backbone_transform)
    normal_structure = CompleteStructureState(spec.candidate_id + "|normal", states, {spec.bond_id: normal}, "inherited", "inherited")
    modified_structure = CompleteStructureState(spec.candidate_id + "|pt", states, {spec.bond_id: modified}, "inherited", "inherited")
    normal_fragment = extract_fragment_from_structure(normal_structure, spec.sequence, spec.fragment_start,
        spec.fragment_end, fragment_id=spec.candidate_id + "|normal", fragment_type=spec.enzyme,
        terminal_form=spec.terminal_form)
    modified_fragment = extract_fragment_from_structure(modified_structure, spec.sequence, spec.fragment_start,
        spec.fragment_end, fragment_id=spec.candidate_id + "|pt", fragment_type=spec.enzyme,
        terminal_form=spec.terminal_form)
    delta = composition_difference(modified_fragment.elemental_composition, normal_fragment.elemental_composition)
    error = (modified_fragment.neutral_exact_mass - normal_fragment.neutral_exact_mass
        - backbone_transform.exact_mass_delta)
    block_status, block_evidence = backbone_transform.rule_for(spec.enzyme)
    return PTPair(spec, cleavage, normal_structure, modified_structure, normal_fragment, modified_fragment,
        shared, delta, backbone_transform.composition_delta, error, block_status, block_evidence)


def discovery_pair_specs(sequence: str, sequence_id: str, enzyme: str, transformations: list[Any],
    slot_schema_path: str | Path, *, organism: str = "", pathway_context: str | None = None,
    max_nucleoside_modifications: int = 1) -> tuple[list[PTPairSpec], list[dict[str, Any]]]:
    specs = []; metadata = []
    for bond in discover_candidate_cleavage_bonds(sequence, sequence_id, enzyme):
        eligible = [t for t in transformations if t.target_scope == "nucleoside"
            and bond.left_base in t.parent_bases and (not t.allowed_enzymes or enzyme in t.allowed_enzymes)]
        result = compose_modifications(bond.left_base, bond.left_position, eligible, slot_schema_path,
            max_components=max_nucleoside_modifications, pathway_context=pathway_context,
            organism_context=organism)
        state_sets = [()] + [tuple(c.state.applied_transform_ids) for c in result.valid_candidates
            if c.position_compatible and c.pathway_compatible]
        for state_index, transform_ids in enumerate(state_sets):
            state_name = "unmodified" if not transform_ids else "+".join(transform_ids)
            cid = f"DISC|{enzyme}|{bond.bond_id}|{state_name}"
            positions = () if not transform_ids else (PositionStateSpec(bond.left_position, bond.left_base, transform_ids),)
            specs.append(PTPairSpec(cid, "", "discovery", sequence_id, sequence, enzyme, bond.bond_id,
                positions, bond.fragment_start, bond.fragment_end, "default"))
        metadata.append({"candidate": bond, "state_count": len(state_sets),
            "invalid_count": len(result.invalid_attempts), "pair_count": len(state_sets)})
    return specs, metadata
