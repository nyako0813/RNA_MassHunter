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


def resolve_rule_inheritance(rule_dir: str | Path, rule_data: dict[str, Any], warnings: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    parent_name = rule_data.get("inherits")
    if not parent_name:
        return rule_data
    parent_path = Path(rule_dir) / f"{parent_name}.yaml"
    if not parent_path.exists():
        if warnings is not None:
            add_warning(warnings, "ERROR", "rule_loader", "Inherited rule_set was not found.", str(parent_path))
        return rule_data
    with parent_path.open("r", encoding="utf-8") as handle:
        parent_data = yaml.safe_load(handle) or {}
    parent_resolved = resolve_rule_inheritance(rule_dir, parent_data, warnings)
    return _merge_position_rules(parent_resolved, rule_data)


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
