"""Generic, user-configured biological context prioritization for MVP-5.6."""

from typing import Any


CONTEXT_PRIORITY_COLUMNS = ["Context_Field", "Value", "Source", "Notes"]


def _bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return list(value)
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return [value]


def _float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _contains_token(value: Any, token: str) -> bool:
    token = token.casefold()
    if isinstance(value, dict):
        return any(_contains_token(key, token) or _contains_token(child, token) for key, child in value.items())
    if isinstance(value, (list, tuple, set)):
        return any(_contains_token(child, token) for child in value)
    return token in str(value or "").casefold()


def _pathway_support(modification_id: str, pathways: list[dict[str, Any]] | None) -> bool:
    return any(_contains_token(pathway, modification_id) for pathway in pathways or []) if modification_id else False


def biological_context_priority_rows(config: Any) -> list[dict[str, Any]]:
    context = getattr(config, "biological_context", {}) or {}
    ranking = getattr(config, "modification_evidence_ranking", {}) or {}
    fields = [
        "enabled", "organism_group", "organism_species", "trna_name", "trna_type", "anticodon",
        "focus_positions", "focus_position_window", "priority_modifications", "priority_keywords",
    ]
    rows = [{
        "Context_Field": field, "Value": ";".join(map(str, _list(context.get(field)))) if isinstance(context.get(field), list) else context.get(field, ""),
        "Source": "config.biological_context", "Notes": "User-configured; empty values do not boost candidates.",
    } for field in fields]
    rows.extend([
        {"Context_Field": "boost settings", "Value": str(context.get("boost", {})), "Source": "config.biological_context", "Notes": "Additive ranking weights."},
        {"Context_Field": "penalty settings", "Value": str(context.get("penalties", {})), "Source": "config.biological_context", "Notes": "Context conflict/unrelated-family penalties."},
        {"Context_Field": "confidence gating", "Value": str({
            "cap_context_only_confidence": ranking.get("cap_context_only_confidence", "Medium"),
            "require_ms_evidence_for_context_boosted_high": ranking.get("require_ms_evidence_for_context_boosted_high", True),
        }), "Source": "config.modification_evidence_ranking", "Notes": "Context alone cannot establish high confidence."},
    ])
    return rows


def score_biological_context(
    candidate: dict[str, Any],
    modification: Any,
    config: Any,
    rule_set: dict[str, Any] | None = None,
    pathways: list[dict[str, Any]] | None = None,
    rule_supported: bool = False,
) -> dict[str, Any]:
    context = getattr(config, "biological_context", {}) or {}
    if not _bool(context.get("enabled"), True):
        return _empty_result()
    boosts = context.get("boost", {}) or {}
    penalties = context.get("penalties", {}) or {}
    score = 0.0
    notes = []
    mod_raw = getattr(modification, "raw", {}) or {}
    identities = {
        str(candidate.get("Modification_ID") or "").casefold(),
        str(getattr(modification, "symbol", "") or "").casefold(),
        str(mod_raw.get("short_name") or "").casefold(),
    }
    priorities = {str(value).casefold() for value in _list(context.get("priority_modifications"))}
    matched_priority = bool(identities & priorities) if priorities else False
    if matched_priority:
        score += _float(boosts.get("priority_modification"), 1.5)
        notes.append("user priority modification matched")

    keywords = [str(value) for value in _list(context.get("priority_keywords"))]
    searchable = " ".join(str(value or "") for value in (
        candidate.get("Modification_ID"), candidate.get("Modification_Name"), candidate.get("Chemical_Group"),
        candidate.get("Near_Isobaric_Group"), candidate.get("Source_Priority"), candidate.get("Notes"),
        mod_raw.get("name"), mod_raw.get("chemical_group"), mod_raw.get("near_isobaric_group"),
    )).casefold()
    matched_keywords = [keyword for keyword in keywords if keyword.casefold() in searchable]
    if matched_keywords:
        score += _float(boosts.get("priority_keyword_match"), 0.75)
        notes.append("user priority keyword matched")

    focus_positions = []
    for value in _list(context.get("focus_positions")):
        try:
            focus_positions.append(int(value))
        except (TypeError, ValueError):
            continue
    candidate_position = candidate.get("Candidate_tRNA_Position")
    try:
        position = int(candidate_position)
    except (TypeError, ValueError):
        position = None
    distance = min((abs(position - focus) for focus in focus_positions), default=None) if position is not None else None
    focus_match = ""
    if distance == 0:
        score += _float(boosts.get("focus_position_match"), 1.0)
        focus_match = "exact"
        notes.append("user focus position matched")
    elif distance is not None and distance <= int(context.get("focus_position_window", 2) or 2):
        score += _float(boosts.get("focus_position_nearby"), 0.5)
        focus_match = "nearby"
        notes.append("candidate is near a user focus position")

    configured_organism = getattr(config, "organism", {}) or {}
    requested_group = str(context.get("organism_group") or "").strip()
    requested_species = str(context.get("organism_species") or "").strip()
    conflict = False
    if requested_group and configured_organism.get("group") and requested_group.casefold() != str(configured_organism.get("group")).casefold():
        conflict = True
    if requested_species and configured_organism.get("species") and requested_species.casefold() != str(configured_organism.get("species")).casefold():
        conflict = True
    organism_supported = bool((requested_group or requested_species) and rule_supported and not conflict)
    if organism_supported:
        score += _float(boosts.get("organism_rule_supported"), 1.0)
        notes.append("explicit organism rule supports modification")

    pathway_supported = _pathway_support(str(candidate.get("Modification_ID") or ""), pathways)
    if pathway_supported:
        score += _float(boosts.get("pathway_supported"), 1.0)
        notes.append("loaded pathway explicitly references modification")

    requested_trna_type = str(context.get("trna_type") or "").strip()
    requested_anticodon = str(context.get("anticodon") or "").strip()
    trna_explicit = bool(requested_trna_type or requested_anticodon)
    trna_supported = bool(trna_explicit and (
        (requested_trna_type and _contains_token(rule_set or {}, requested_trna_type)) or
        (requested_anticodon and _contains_token(rule_set or {}, requested_anticodon))
    ))
    if trna_supported:
        score += _float(boosts.get("trna_context_supported"), 1.0)
        notes.append("rule set explicitly supports configured tRNA context")

    if conflict:
        score += _float(penalties.get("organism_context_conflict"), -2.0)
        notes.append("configured biological context conflicts with run organism")
    if conflict:
        level = "Conflict"
    elif score >= 3:
        level = "High"
    elif score >= 1.5:
        level = "Medium"
    elif score > 0:
        level = "Low"
    else:
        level = "None"
    return {
        "Biological_Context_Score": score, "Biological_Context_Level": level,
        "Biological_Context_Notes": "; ".join(notes),
        "Context_Matched_Priority_Modification": matched_priority,
        "Context_Matched_Keywords": ";".join(matched_keywords),
        "Context_Focus_Position_Match": focus_match,
        "Context_Focus_Position_Distance": distance if distance is not None else "",
        "Context_Pathway_Supported": pathway_supported, "Context_Organism_Supported": organism_supported,
        "Context_TRNA_Supported": trna_supported, "Context_Conflict": conflict,
    }


def _empty_result() -> dict[str, Any]:
    return {
        "Biological_Context_Score": 0.0, "Biological_Context_Level": "None", "Biological_Context_Notes": "",
        "Context_Matched_Priority_Modification": False, "Context_Matched_Keywords": "",
        "Context_Focus_Position_Match": "", "Context_Focus_Position_Distance": "",
        "Context_Pathway_Supported": False, "Context_Organism_Supported": False,
        "Context_TRNA_Supported": False, "Context_Conflict": False,
    }
