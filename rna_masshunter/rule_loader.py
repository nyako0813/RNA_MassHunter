from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

from rna_masshunter.warnings_manager import add_warning


def _merge_position_rules(parent: dict[str, Any], child: dict[str, Any]) -> dict[str, Any]:
    merged = deepcopy(parent)
    merged.update({k: v for k, v in child.items() if k != "position_rules"})
    parent_rules = {rule.get("id"): rule for rule in parent.get("position_rules", []) if isinstance(rule, dict)}
    for rule in child.get("position_rules", []) or []:
        if isinstance(rule, dict):
            parent_rules[rule.get("id")] = rule
    merged["position_rules"] = list(parent_rules.values())
    return merged


def _inherit_names(inherits: Any) -> list[str]:
    if not inherits:
        return []
    if isinstance(inherits, str):
        return [inherits]
    if isinstance(inherits, list):
        return [str(item) for item in inherits if item]
    return [str(inherits)]


def resolve_rule_inheritance(
    rule_dir: str | Path,
    rule_data: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
    seen: set[str] | None = None,
) -> dict[str, Any]:
    parent_names = _inherit_names(rule_data.get("inherits"))
    if not parent_names:
        return rule_data
    seen = seen or set()
    merged_parent: dict[str, Any] = {}
    for parent_name in parent_names:
        parent_path = Path(rule_dir) / f"{parent_name}.yaml"
        if parent_name in seen:
            if warnings is not None:
                add_warning(warnings, "ERROR", "rule_loader", "Circular rule_set inheritance detected.", parent_name)
            continue
        if not parent_path.exists():
            if warnings is not None:
                add_warning(warnings, "ERROR", "rule_loader", "Inherited rule_set was not found.", str(parent_path))
            continue
        with parent_path.open("r", encoding="utf-8") as handle:
            parent_data = yaml.safe_load(handle) or {}
        parent_resolved = resolve_rule_inheritance(rule_dir, parent_data, warnings, seen | {parent_name})
        merged_parent = _merge_position_rules(merged_parent, parent_resolved)
    return _merge_position_rules(merged_parent, rule_data)


def load_rule_set(rule_dir: str | Path, rule_set_name: str, warnings: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    path = Path(rule_dir) / f"{rule_set_name}.yaml"
    if not path.exists():
        if warnings is not None:
            add_warning(warnings, "ERROR", "rule_loader", "Rule set file was not found.", str(path))
        return {}
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return resolve_rule_inheritance(rule_dir, data, warnings)


def validate_rule_set(rule_set: dict[str, Any], warnings: list[dict[str, Any]] | None = None) -> None:
    if not rule_set and warnings is not None:
        add_warning(warnings, "ERROR", "rule_loader", "Rule set is empty.")
        return
    if "position_rules" in rule_set and not isinstance(rule_set["position_rules"], list) and warnings is not None:
        add_warning(warnings, "ERROR", "rule_loader", "position_rules must be a list.")
