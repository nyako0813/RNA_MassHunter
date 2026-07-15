"""Schema and validation for explicit sample structure hypotheses."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping
import hashlib
import yaml

ALLOWED_TERMINAL_STATES = {"default", "inherited", "dephosphorylated", "residual_phosphate", "cyclic_phosphate"}

@dataclass(frozen=True)
class PositionHypothesis:
    position: int
    parent_base: str
    transformation_ids: tuple[str, ...]

@dataclass(frozen=True)
class BackboneHypothesis:
    bond_id: str
    state: str

@dataclass(frozen=True)
class SampleStructureHypothesis:
    hypothesis_id: str
    positions: tuple[PositionHypothesis, ...]
    backbone: tuple[BackboneHypothesis, ...]
    five_prime: str = "inherited"
    three_prime: str = "inherited"

@dataclass(frozen=True)
class SampleStructureLoadResult:
    schema_version: int
    enabled: bool
    hypotheses: tuple[SampleStructureHypothesis, ...]
    invalid_rows: tuple[dict[str, Any], ...]

def _invalid(hypothesis_id: str, reason: str, detail: Any = "") -> dict[str, Any]:
    return {"Candidate_ID": hypothesis_id, "Valid": False, "Invalid_Reason": reason,
            "Invalid_Detail": detail, "Applied_To_Formal_Result": False, "Formal_Change_Ready": False}

def load_sample_structure_hypotheses(path: str | Path, *, sequence: str = "",
    transformations: list[Any] | tuple[Any, ...] = (),
    backbone_bond_ids: set[str] | None = None,
    target_identity: Mapping[str, Any] | None = None) -> SampleStructureLoadResult:
    """Load a fixture, validating identifiers and sample coordinates."""
    source = Path(path)
    if not source.exists():
        return SampleStructureLoadResult(1, False, (), (_invalid("", "fixture_not_found", str(source)),))
    with source.open(encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    schema_version = int(raw.get("schema_version", 1) or 1)
    enabled = bool(raw.get("enabled", True))
    transform_map = {str(item.id): item for item in transformations}
    seq = str(sequence or "").upper().replace("T", "U")
    required_target = raw.get("target") or {}
    actual_target = dict(target_identity or {})
    target_failures: list[str] = []
    
    if required_target:
        checks = {
            "name": str(actual_target.get("name") or ""),
            "length": len(seq),
            "sequence_sha256": (
                hashlib.sha256(seq.encode("utf-8")).hexdigest()
                if seq
                else ""
            ),
            "organism": str(actual_target.get("organism") or ""),
            "rule_set": str(actual_target.get("rule_set") or ""),
        }

        raw_aliases = required_target.get("name_aliases") or ()
        if isinstance(raw_aliases, str):
            aliases = {raw_aliases}
        else:
            aliases = {str(value) for value in raw_aliases}

        for field, required in required_target.items():
            if field == "name_aliases":
                continue

            actual = str(checks.get(field, ""))
            matches = actual == str(required)

            if field == "name":
                matches = matches or actual in aliases

            if field not in checks or not matches:
                target_failures.append(
                    f"{field}:{actual}!={required}"
                )

    if target_failures:
        return SampleStructureLoadResult(schema_version, enabled, (),
            (_invalid("", "target_identity_mismatch", ";".join(target_failures)),))
    seen: set[str] = set()
    valid: list[SampleStructureHypothesis] = []
    invalid: list[dict[str, Any]] = []
    for index, payload in enumerate(raw.get("hypotheses") or (), 1):
        hid = str((payload or {}).get("hypothesis_id") or f"HYP_{index:04d}")
        reasons: list[str] = []
        details: list[str] = []
        if hid in seen:
            reasons.append("duplicate_hypothesis_id")
        seen.add(hid)
        position_states: list[PositionHypothesis] = []
        for position_raw, state_raw in ((payload or {}).get("positions") or {}).items():
            try:
                position = int(position_raw)
            except (TypeError, ValueError):
                reasons.append("non_integer_position"); details.append(str(position_raw)); continue
            state_raw = state_raw if isinstance(state_raw, Mapping) else {}
            parent = str(state_raw.get("parent_base") or "").upper().replace("T", "U")
            ids = tuple(str(item) for item in (state_raw.get("transformations") or ()))
            if not seq or position < 1 or position > len(seq):
                reasons.append("position_out_of_range"); details.append(str(position))
            elif parent != seq[position - 1]:
                reasons.append("parent_base_mismatch"); details.append(f"{position}:{parent}!={seq[position - 1]}")
            if not ids:
                reasons.append("empty_transform_list"); details.append(str(position))
            for transform_id in ids:
                transform = transform_map.get(transform_id)
                if transform is None:
                    reasons.append("unknown_transform_id"); details.append(transform_id)
                elif parent not in tuple(transform.parent_bases):
                    reasons.append("transform_parent_base_incompatible"); details.append(transform_id)
            position_states.append(PositionHypothesis(position, parent, ids))
        backbone_states: list[BackboneHypothesis] = []
        for item in (payload or {}).get("backbone") or ():
            item = item if isinstance(item, Mapping) else {"bond_id": item}
            bond_id = str(item.get("bond_id") or "")
            state = str(item.get("state") or "")
            if not bond_id or (backbone_bond_ids is not None and bond_id not in backbone_bond_ids):
                reasons.append("unknown_backbone_bond"); details.append(bond_id)
            if state not in {"normal_phosphate", "phosphorothioate"}:
                reasons.append("unsupported_backbone_state"); details.append(state)
            backbone_states.append(BackboneHypothesis(bond_id, state))
        terminal = (payload or {}).get("terminal_state") or {}
        if not isinstance(terminal, Mapping):
            reasons.append("invalid_terminal_state_format"); terminal = {}
        five = str(terminal.get("five_prime") or "inherited")
        three = str(terminal.get("three_prime") or "inherited")
        if five not in ALLOWED_TERMINAL_STATES or three not in ALLOWED_TERMINAL_STATES:
            reasons.append("unsupported_terminal_state"); details.extend((five, three))
        if reasons:
            invalid.append(_invalid(hid, ";".join(dict.fromkeys(reasons)), ";".join(details)))
        else:
            valid.append(SampleStructureHypothesis(hid, tuple(position_states), tuple(backbone_states), five, three))
    return SampleStructureLoadResult(schema_version, enabled, tuple(valid) if enabled else (), tuple(invalid))