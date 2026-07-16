"""Generic position-constrained P1+AP dinucleotide candidate groups.

Sequence positions are candidate-construction constraints only.  Observations
from P1+AP preparation never localize a sequence position or source bond.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, Iterable
import math
import yaml

from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.modification_composer import compose_modifications
from rna_masshunter.modification_constraints import Transformation, load_transformations

DINUCLEOTIDE_MODEL_VERSION = "P1_SAP_DINUCLEOTIDE_CHEMICAL_STATE_v2.0"
MODEL_NOT_DEFINED = "MODEL_NOT_DEFINED"
MODEL_NOT_APPLICABLE = "MODEL_NOT_APPLICABLE"
MS2_MODEL_REASON = "P1_SAP_DINUCLEOTIDE_REQUIRES_VALIDATED_FRAGMENT_MODEL"

FORMAL_FALSE = {
    "Applied_To_Formal_Result": False,
    "Formal_Change_Ready": False,
    "Formal_Result_Changed": False,
}
LOCALIZATION_FALSE = {
    "Sequence_Position_Localized": False,
    "Original_Bond_Localized": False,
    "Position_Localized": False,
    "Position_Resolvable": False,
    "Bond_Resolvable": False,
    "Localization_Status": "POSITION_NOT_RETAINED_BY_PREPARATION",
    "Position_Localization_Status": "POSITION_NOT_RETAINED_BY_PREPARATION",
    "Source_Bond_Resolution_Status": "SOURCE_BOND_UNRESOLVED",
}
NUCLEOSIDE_COMPOSITIONS = {
    "A": ElementalComposition({"C": 10, "H": 13, "N": 5, "O": 4}),
    "C": ElementalComposition({"C": 9, "H": 13, "N": 3, "O": 5}),
    "G": ElementalComposition({"C": 10, "H": 13, "N": 5, "O": 5}),
    "U": ElementalComposition({"C": 9, "H": 12, "N": 2, "O": 6}),
}
LINKAGE_COMPOSITIONS = {
    "NORMAL_PHOSPHATE": ElementalComposition({"H": 1, "O": 3, "P": 1}),
    "PHOSPHOROTHIOATE": ElementalComposition({"H": 1, "O": 2, "P": 1, "S": 1}),
}
CONDENSATION_ADJUSTMENT = ElementalComposition.delta({"H": -2, "O": -1})
PT_DELTA = ElementalComposition.delta({"O": -1, "S": 1})

GROUP_COLUMNS = [
    "Dinucleotide_Group_ID", "Chemical_State_ID", "Model_Version",
    "Final_Elemental_Composition", "Elemental_Composition", "Neutral_Mass",
    "Theoretical_mz", "Charge", "Polarity", "Linkage_State", "Chemical_Family",
    "Product_Type", "Observable", "Search_Enabled", "Not_Observable_Reason",
    "Observable_In_Acquisition", "Observable_In_MS1_Extraction",
    "Observable_In_Dinucleotide_Search", "Search_Executed", "Structural_Assignment_Count",
    "Possible_Source_Bond_Count", "Possible_Source_Bonds", "Possible_Left_States",
    "Possible_Right_States", "Possible_Structural_Isomers", "Possible_Position_Assignments",
    "Representative_Assignment", "Representative_Is_Confirmed", "Raw_Profile_Point_Count",
    "Unique_MS1_Spectrum_Count", "Observed_Min_mz", "Observed_Max_mz", "Raw_RT_Start",
    "Raw_RT_End", "Feature_Count", "Qualified_Feature_Count", "Group_Interpretation",
    "Position_Constraint_Summary",
    "Chemical_Constraint_Summary", "Composition_Resolution_Status",
    "Linkage_Resolution_Status", "Structure_Resolution_Status",
    "Source_Bond_Resolution_Status", "Candidate_Generation_Truncated",
    "Candidate_Generation_Reason", "Sequence_Position_Localized",
    "Original_Bond_Localized", "Position_Localized", "Position_Localization_Status",
    "Localization_Status", "Applied_To_Formal_Result", "Formal_Change_Ready",
    "Formal_Result_Changed",
]
ASSIGNMENT_COLUMNS = [
    "Structural_Assignment_ID", "Dinucleotide_Group_ID", "Left_Position", "Right_Position",
    "Possible_Source_Bond", "Sequence_Direction", "Left_Base", "Right_Base",
    "Left_State_ID", "Right_State_ID", "Left_Modifications", "Right_Modifications",
    "Linkage_State", "Left_Composition", "Right_Composition", "Linkage_Composition",
    "Condensation_Adjustment", "Final_Elemental_Composition", "Neutral_Mass",
    "Theoretical_mz", "Charge", "Polarity", "Base_Compatible", "Position_Compatible",
    "Organism_Compatible", "tRNA_Compatible", "Composite_Compatible",
    "Chemical_Constraint_Status", "Biological_Plausibility", "Sequence_Position_Localized",
    "Original_Bond_Localized", "Position_Localization_Status", "Source_Bond_Resolution_Status",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]
SUMMARY_COLUMNS = [
    "Model_Version", "Input_Sequence_Length", "Adjacent_Bond_Count", "Position_State_Count",
    "Raw_Candidate_Count", "Constraint_Rejected_Count", "Grouped_Candidate_Count",
    "Observable_Group_Count", "Searched_Group_Count", "Raw_Matched_Group_Count",
    "Qualified_Group_Count", "Qualified_Feature_Count", "Normal_Phosphate_Group_Count",
    "PT_Group_Count", "Qualified_Normal_Phosphate_Group_Count", "Qualified_PT_Group_Count",
    "Competition_Unresolved_Group_Count", "Isotope_Assessed_Feature_Count",
    "Isotope_Compatible_Feature_Count", "Isotope_Incompatible_Feature_Count",
    "Precursor_Compatible_MS2_Count", "Candidate_Generation_Runtime", "Raw_Matching_Runtime",
    "Spectrum_Level_Grouping_Runtime", "Chromatographic_Grouping_Runtime", "Isotope_Runtime",
    "Competition_Runtime", "MS2_Provenance_Runtime", "Interpretation_Runtime",
    "Total_Shadow_Audit_Runtime", "Maximum_RSS_MiB", "Tracemalloc_Peak_MiB",
    "Candidate_Generation_Truncated",
    "Candidate_Generation_Reason", "Sequence_Position_Localized", "Original_Bond_Localized",
    "Position_Localization_Status", "Source_Bond_Resolution_Status",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]
# Backward-compatible names for report integrations; values are generic schemas.
CANDIDATE_COLUMNS = GROUP_COLUMNS


@dataclass(frozen=True)
class PositionState:
    position: int
    base: str
    state_id: str
    transform_ids: tuple[str, ...]
    composition: ElementalComposition
    base_compatible: bool = True
    position_compatible: bool = True
    organism_compatible: bool = True
    trna_compatible: bool = True
    composite_compatible: bool = True
    structural_id: str = ""


@dataclass
class DinucleotideCandidateResult:
    candidates: list[dict[str, Any]]
    assignments: list[dict[str, Any]]
    position_states: dict[int, list[PositionState]]
    rejected: list[dict[str, Any]]
    summary: dict[str, Any]
    truncated: bool = False
    truncation_reason: str = ""

    @property
    def raw_candidates(self) -> list[dict[str, Any]]:
        """Compatibility alias; raw candidates are structural assignments in v2."""
        return self.assignments

    def to_jsonable(self) -> dict[str, Any]:
        return {
            "model_version": DINUCLEOTIDE_MODEL_VERSION,
            "candidate_generation": dict(self.summary),
            "groups": self.candidates,
            "assignments": self.assignments,
            "constraint_rejections": self.rejected,
        }


def adjacent_bonds(sequence: str) -> list[tuple[int, int, str, str]]:
    normalized = str(sequence).upper().replace("T", "U")
    return [(i, i + 1, normalized[i - 1], normalized[i]) for i in range(1, len(normalized))]


def dinucleotide_settings(config: Any | None) -> dict[str, Any]:
    """Normalize nested v2 configuration while accepting v1 flat keys."""
    section = dict((getattr(config, "p1_sap_dinucleotide", {}) or {})) if config is not None else {}
    generation = dict(section.get("candidate_generation") or {})
    search = dict(section.get("search") or {})
    mass_accuracy = dict(section.get("mass_accuracy") or {})
    feature_quality = dict(section.get("feature_quality") or {})
    isotope = dict(section.get("isotope") or {})
    ms2 = dict(section.get("ms2_provenance") or {})
    for key in ("max_modifications_per_side", "max_composite_states_per_position", "max_candidate_count"):
        if key in section and key not in generation:
            generation[key] = section[key]
    if "candidate_mass_tolerance_ppm" in section and "tolerance_ppm" not in search:
        search["tolerance_ppm"] = section["candidate_mass_tolerance_ppm"]
    instrument = (getattr(config, "instrument", {}) or {}) if config is not None else {}
    p1 = (getattr(config, "p1_annotation", {}) or {}) if config is not None else {}
    search_min = float(search.get("mz_min", 100.0))
    search_max = float(search.get("mz_max", 1000.0))
    acquisition_min = float(search.get("acquisition_mz_min", instrument.get("ms1_mz_min", search_min)))
    acquisition_max = float(search.get("acquisition_mz_max", instrument.get("ms1_mz_max", search_max)))
    extraction_min = float(search.get("ms1_extraction_mz_min", p1.get("mz_min", search_min)))
    extraction_max = float(search.get("ms1_extraction_mz_max", p1.get("mz_max", search_max)))
    return {
        "enabled": bool(section.get("enabled", True)),
        "max_modifications_per_side": max(0, int(generation.get("max_modifications_per_side", 3))),
        "max_composite_states_per_position": max(1, int(generation.get("max_composite_states_per_position", 64))),
        "max_candidate_count": max(1, int(generation.get("max_candidate_count", 100000))),
        "include_normal_phosphate": bool(generation.get("include_normal_phosphate", True)),
        "include_phosphorothioate": bool(generation.get("include_phosphorothioate", True)),
        "charges": tuple(int(value) for value in generation.get("charges", section.get("charges", [1]))),
        "polarity": str(generation.get("polarity", section.get("polarity", "positive"))).lower(),
        "strict_positions": {str(k): {int(x) for x in v} for k, v in (generation.get("strict_positions") or {}).items()},
        "search_mz_min": search_min, "search_mz_max": search_max,
        "acquisition_mz_min": acquisition_min, "acquisition_mz_max": acquisition_max,
        "ms1_extraction_mz_min": extraction_min, "ms1_extraction_mz_max": extraction_max,
        "search_tolerance_ppm": float(search.get("tolerance_ppm", 10.0)),
        "strong_mass_accuracy_ppm": float(mass_accuracy.get("strong_ppm", 2.0)),
        "moderate_mass_accuracy_ppm": float(mass_accuracy.get("moderate_ppm", 5.0)),
        "search_mass_accuracy_ppm": float(mass_accuracy.get("search_ppm", search.get("tolerance_ppm", 10.0))),
        "feature_quality": feature_quality,
        "isotope_enabled": bool(isotope.get("enabled", True)),
        "isotope_tolerance_ppm": float(isotope.get("tolerance_ppm", 20.0)),
        "isotope_require_same_scan": bool(isotope.get("require_same_scan", True)),
        "ms2_provenance_enabled": bool(ms2.get("enabled", True)),
        "targets": list(section.get("targets") or []),
    }


def _matching_position_hypotheses(root: Path, sequence: str, config: Any | None) -> dict[str, set[int]]:
    path = root / "data/modification_position_hypotheses.yaml"
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    normalized = str(sequence).upper().replace("T", "U")
    digest = sha256(normalized.encode("utf-8")).hexdigest()
    sequence_cfg = (getattr(config, "sequence", {}) or {}) if config is not None else {}
    organism_cfg = (getattr(config, "organism", {}) or {}) if config is not None else {}
    sequence_name = str(sequence_cfg.get("name", ""))
    organism = str(organism_cfg.get("species", ""))
    positions: dict[str, set[int]] = defaultdict(set)
    for target in raw.get("targets", []):
        hash_value = str(target.get("sequence_sha256", ""))
        names = {str(target.get(key, "")) for key in ("sequence_id", "sequence_name")}
        names.update(str(x) for key in ("sequence_id_aliases", "sequence_name_aliases") for x in target.get(key, []))
        target_organism = str(target.get("organism", ""))
        identified = bool(hash_value or any(names) or target_organism)
        if hash_value:
            matches = hash_value == digest
        elif any(names) and sequence_name:
            matches = sequence_name in names
        elif target_organism and organism:
            matches = target_organism == organism and (not target.get("sequence_length") or int(target["sequence_length"]) == len(normalized))
        else:
            matches = not identified
        if identified and not matches:
            continue
        for row in target.get("hypotheses", []):
            if row.get("position") is None:
                continue
            ids = row.get("modification_ids") or [row.get("modification_id")]
            for mod_id in ids:
                if mod_id and str(mod_id) != "phosphorothioate":
                    positions[str(mod_id)].add(int(row["position"]))
    return dict(positions)


def _family_landmarks(root: Path, config: Any | None) -> dict[str, set[int]]:
    path = root / "data/modification_position_priors.yaml"
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as handle:
        raw = yaml.safe_load(handle) or {}
    sequence_cfg = (getattr(config, "sequence", {}) or {}) if config is not None else {}
    configured = {str(k): int(v) for k, v in (sequence_cfg.get("canonical_landmarks") or {}).items()}
    if sequence_cfg.get("wobble_position") is not None:
        configured.setdefault("wobble", int(sequence_cfg["wobble_position"]))
    result: dict[str, set[int]] = {}
    for family in raw.get("families", []):
        position = configured.get(str(family.get("canonical_landmark")))
        if position is None:
            continue
        for mod_id in family.get("modification_ids", []):
            result.setdefault(str(mod_id), set()).add(position)
    return result


def _allowed_positions(transform: Transformation, exact: dict[str, set[int]], landmarks: dict[str, set[int]], strict: dict[str, set[int]]) -> set[int]:
    if transform.id in strict:
        return strict[transform.id]
    if transform.allowed_positions:
        return set(transform.allowed_positions)
    if transform.id in exact:
        return exact[transform.id]
    return landmarks.get(transform.id, set())


def build_position_states(
    sequence: str,
    project_root: str | Path,
    *,
    config: Any | None = None,
    organism: str | None = None,
    max_modifications_per_side: int | None = None,
    max_composite_states_per_position: int | None = None,
) -> tuple[dict[int, list[PositionState]], list[dict[str, Any]]]:
    root = Path(project_root)
    settings = dinucleotide_settings(config)
    max_components = settings["max_modifications_per_side"] if max_modifications_per_side is None else int(max_modifications_per_side)
    max_states = settings["max_composite_states_per_position"] if max_composite_states_per_position is None else int(max_composite_states_per_position)
    transforms = load_transformations(root / "data/modification_transforms_v2.yaml")
    exact = _matching_position_hypotheses(root, sequence, config)
    landmarks = _family_landmarks(root, config)
    schema_path = root / "data/nucleoside_slots.yaml"
    normalized = str(sequence).upper().replace("T", "U")
    organism_value = organism or str(((getattr(config, "organism", {}) or {}).get("species", "") if config is not None else ""))
    states_by_position: dict[int, list[PositionState]] = {}
    rejected: list[dict[str, Any]] = []
    for position, base in enumerate(normalized, 1):
        if base not in NUCLEOSIDE_COMPOSITIONS:
            states_by_position[position] = []
            rejected.append({"Position": position, "Base": base, "Reason": "unsupported_parent_base"})
            continue
        states = [PositionState(position, base, base, (), NUCLEOSIDE_COMPOSITIONS[base], structural_id=f"{base}@{position}|unmodified")]
        if max_components:
            composed = compose_modifications(base, position, transforms, schema_path, max_components=max_components, organism_context=organism_value or None)
            for invalid in composed.invalid_attempts:
                rejected.append({"Position": position, "Base": base, "Transform_IDs": ";".join(invalid.transform_ids), "Reason": invalid.result.reason_code})
            for candidate in composed.valid_candidates:
                ids = tuple(transform.id for transform in candidate.transforms)
                disallowed = [transform.id for transform in candidate.transforms if _allowed_positions(transform, exact, landmarks, settings["strict_positions"]) and position not in _allowed_positions(transform, exact, landmarks, settings["strict_positions"])]
                if disallowed or not candidate.position_compatible:
                    rejected.append({"Position": position, "Base": base, "Transform_IDs": ";".join(ids), "Reason": "position_disallowed"})
                    continue
                if not candidate.pathway_compatible:
                    rejected.append({"Position": position, "Base": base, "Transform_IDs": ";".join(ids), "Reason": "organism_or_pathway_disallowed"})
                    continue
                try:
                    composition = NUCLEOSIDE_COMPOSITIONS[base] + candidate.state.elemental_composition_delta
                except ValueError:
                    rejected.append({"Position": position, "Base": base, "Transform_IDs": ";".join(ids), "Reason": MODEL_NOT_DEFINED})
                    continue
                states.append(PositionState(position, base, "+".join(ids), ids, composition, organism_compatible=candidate.pathway_compatible, structural_id=candidate.state.canonical_structure_id))
        states.sort(key=lambda state: (len(state.transform_ids), state.state_id, state.structural_id))
        if len(states) > max_states:
            rejected.extend({"Position": position, "Base": base, "Transform_IDs": ";".join(state.transform_ids), "Reason": "max_composite_states_per_position"} for state in states[max_states:])
            states = states[:max_states]
        states_by_position[position] = states
    return states_by_position, rejected


def _combine(*parts: ElementalComposition) -> ElementalComposition | None:
    try:
        result = ElementalComposition()
        for part in parts:
            result = result + part
        return result if all(value >= 0 for value in result.to_dict().values()) else None
    except (TypeError, ValueError):
        return None


def _assignment(left: PositionState, right: PositionState, linkage_state: str, charge: int, polarity: str) -> dict[str, Any] | None:
    linkage = LINKAGE_COMPOSITIONS.get(linkage_state)
    if linkage is None:
        return None
    final = _combine(left.composition, right.composition, linkage, CONDENSATION_ADJUSTMENT)
    if final is None:
        return None
    neutral = final.exact_mass
    theoretical = mz_from_neutral_mass(neutral, charge, polarity)
    return {
        "Left_Position": left.position, "Right_Position": right.position,
        "Possible_Source_Bond": f"{left.position}-{right.position}",
        "Sequence_Direction": "5prime_to_3prime", "Left_Base": left.base, "Right_Base": right.base,
        "Left_State_ID": left.state_id, "Right_State_ID": right.state_id,
        "Left_Modifications": ";".join(left.transform_ids), "Right_Modifications": ";".join(right.transform_ids),
        "Linkage_State": linkage_state, "Left_Composition": left.composition.canonical_string(),
        "Right_Composition": right.composition.canonical_string(), "Linkage_Composition": linkage.canonical_string(),
        "Condensation_Adjustment": CONDENSATION_ADJUSTMENT.canonical_string(),
        "Final_Elemental_Composition": final.canonical_string(), "Neutral_Mass": neutral,
        "Theoretical_mz": theoretical, "Charge": charge, "Polarity": polarity,
        "Base_Compatible": left.base_compatible and right.base_compatible,
        "Position_Compatible": left.position_compatible and right.position_compatible,
        "Organism_Compatible": left.organism_compatible and right.organism_compatible,
        "tRNA_Compatible": left.trna_compatible and right.trna_compatible,
        "Composite_Compatible": left.composite_compatible and right.composite_compatible,
        "Chemical_Constraint_Status": "COMPATIBLE", "Biological_Plausibility": "CONSTRAINT_COMPATIBLE",
        "Structural_Isomer": f"{left.structural_id}|{linkage_state}|{right.structural_id}",
        **LOCALIZATION_FALSE, **FORMAL_FALSE,
    }


def _group_assignments(assignments: list[dict[str, Any]], settings: dict[str, Any], truncated: bool, reason: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in assignments:
        grouped[(row["Final_Elemental_Composition"], row["Charge"], row["Polarity"], row["Linkage_State"])].append(row)
    groups: list[dict[str, Any]] = []
    output_assignments: list[dict[str, Any]] = []
    for serial, (_key, members) in enumerate(sorted(grouped.items(), key=lambda item: (float(item[1][0]["Theoretical_mz"]), item[0][3], item[0][1])), 1):
        group_id = f"P1SAP_DG_{serial:06d}"
        members.sort(key=lambda row: (row["Left_Position"], row["Right_Position"], row["Left_State_ID"], row["Right_State_ID"]))
        for assignment_number, member in enumerate(members, 1):
            member["Dinucleotide_Group_ID"] = group_id
            member["Structural_Assignment_ID"] = f"{group_id}_A{assignment_number:05d}"
            output_assignments.append(member)
        first = members[0]
        mz = float(first["Theoretical_mz"])
        acquisition = settings["acquisition_mz_min"] <= mz <= settings["acquisition_mz_max"]
        extraction = settings["ms1_extraction_mz_min"] <= mz <= settings["ms1_extraction_mz_max"]
        search = settings["search_mz_min"] <= mz <= settings["search_mz_max"]
        observable = acquisition and extraction and search
        reasons = []
        if not acquisition: reasons.append("OUTSIDE_ACQUISITION_RANGE")
        if not extraction: reasons.append("OUTSIDE_MS1_EXTRACTION_RANGE")
        if not search: reasons.append("OUTSIDE_DINUCLEOTIDE_SEARCH_RANGE")
        bonds = sorted({row["Possible_Source_Bond"] for row in members}, key=lambda value: tuple(map(int, value.split("-"))))
        left_states = sorted({row["Left_State_ID"] for row in members})
        right_states = sorted({row["Right_State_ID"] for row in members})
        isomers = sorted({row["Structural_Isomer"] for row in members})
        assignments_text = [f"{row['Left_Position']}:{row['Left_State_ID']}-{row['Linkage_State']}-{row['Right_Position']}:{row['Right_State_ID']}" for row in members]
        representative = assignments_text[0]
        family = "P1_RESISTANT_PT_DINUCLEOTIDE" if first["Linkage_State"] == "PHOSPHOROTHIOATE" else "NORMAL_PHOSPHATE_DINUCLEOTIDE"
        groups.append({
            "Dinucleotide_Group_ID": group_id, "Chemical_State_ID": group_id,
            "Model_Version": DINUCLEOTIDE_MODEL_VERSION,
            "Final_Elemental_Composition": first["Final_Elemental_Composition"],
            "Elemental_Composition": first["Final_Elemental_Composition"],
            "Neutral_Mass": first["Neutral_Mass"], "Theoretical_mz": mz,
            "Charge": first["Charge"], "Polarity": first["Polarity"],
            "Linkage_State": first["Linkage_State"], "Chemical_Family": family,
            "Product_Type": "dinucleotide", "Observable": observable, "Search_Enabled": observable,
            "Not_Observable_Reason": ";".join(reasons), "Observable_In_Acquisition": acquisition,
            "Observable_In_MS1_Extraction": extraction, "Observable_In_Dinucleotide_Search": search,
            "Search_Executed": False, "Structural_Assignment_Count": len(members),
            "Possible_Source_Bond_Count": len(bonds), "Possible_Source_Bonds": ";".join(bonds),
            "Possible_Left_States": ";".join(left_states), "Possible_Right_States": ";".join(right_states),
            "Possible_Structural_Isomers": ";".join(isomers), "Possible_Position_Assignments": ";".join(assignments_text),
            "Representative_Assignment": representative, "Representative_Is_Confirmed": False,
            "Position_Constraint_Summary": "ALL_ASSIGNMENTS_POSITION_COMPATIBLE",
            "Chemical_Constraint_Summary": "ALL_ASSIGNMENTS_CHEMICALLY_COMPATIBLE",
            "Composition_Resolution_Status": "COMPOSITION_UNRESOLVED_WITHOUT_OBSERVATION",
            "Linkage_Resolution_Status": "LINKAGE_UNRESOLVED_WITHOUT_OBSERVATION",
            "Structure_Resolution_Status": "STRUCTURE_UNRESOLVED",
            "Source_Bond_Resolution_Status": "SOURCE_BOND_UNRESOLVED",
            "Candidate_Generation_Truncated": truncated, "Candidate_Generation_Reason": reason,
            **LOCALIZATION_FALSE, **FORMAL_FALSE,
        })
    return groups, output_assignments


def generate_dinucleotide_candidates(
    sequence: str,
    project_root: str | Path,
    *,
    config: Any | None = None,
    charges: tuple[int, ...] | None = None,
    polarity: str | None = None,
    linkage_states: tuple[str, ...] | None = None,
) -> DinucleotideCandidateResult:
    settings = dinucleotide_settings(config)
    states, rejected = build_position_states(sequence, project_root, config=config)
    selected_charges = tuple(charges or settings["charges"])
    selected_polarity = str(polarity or settings["polarity"]).lower()
    if linkage_states is None:
        selected_linkages = []
        if settings["include_normal_phosphate"]: selected_linkages.append("NORMAL_PHOSPHATE")
        if settings["include_phosphorothioate"]: selected_linkages.append("PHOSPHOROTHIOATE")
    else:
        selected_linkages = list(linkage_states)
    assignments: list[dict[str, Any]] = []
    limit_hit = False
    limit = settings["max_candidate_count"]
    for left_position, right_position, _left_base, _right_base in adjacent_bonds(sequence):
        for left in states.get(left_position, []):
            for right in states.get(right_position, []):
                for linkage in selected_linkages:
                    for charge in selected_charges:
                        row = _assignment(left, right, linkage, charge, selected_polarity)
                        if row is None:
                            rejected.append({"Position": f"{left_position}-{right_position}", "Reason": MODEL_NOT_DEFINED})
                            continue
                        assignments.append(row)
                        if len(assignments) >= limit:
                            limit_hit = True; break
                    if limit_hit: break
                if limit_hit: break
            if limit_hit: break
        if limit_hit: break
    state_limit_hit = any(row.get("Reason") == "max_composite_states_per_position" for row in rejected)
    reasons = []
    if limit_hit: reasons.append(f"max_candidate_count={limit}")
    if state_limit_hit: reasons.append("max_composite_states_per_position")
    truncated = bool(reasons)
    reason = ";".join(reasons)
    groups, assignments = _group_assignments(assignments, settings, truncated, reason)
    summary = {
        "Model_Version": DINUCLEOTIDE_MODEL_VERSION, "Input_Sequence_Length": len(str(sequence)),
        "Adjacent_Bond_Count": len(adjacent_bonds(sequence)), "Position_State_Count": sum(len(value) for value in states.values()),
        "Raw_Candidate_Count": len(assignments), "Constraint_Rejected_Count": len(rejected),
        "Grouped_Candidate_Count": len(groups), "Observable_Group_Count": sum(bool(row["Observable"]) for row in groups),
        "Searched_Group_Count": 0, "Raw_Matched_Group_Count": 0, "Qualified_Group_Count": 0,
        "Qualified_Feature_Count": 0,
        "Normal_Phosphate_Group_Count": sum(row["Linkage_State"] == "NORMAL_PHOSPHATE" for row in groups),
        "PT_Group_Count": sum(row["Linkage_State"] == "PHOSPHOROTHIOATE" for row in groups),
        "Qualified_Normal_Phosphate_Group_Count": 0, "Qualified_PT_Group_Count": 0,
        "Competition_Unresolved_Group_Count": 0, "Isotope_Assessed_Feature_Count": 0,
        "Isotope_Compatible_Feature_Count": 0, "Isotope_Incompatible_Feature_Count": 0,
        "Precursor_Compatible_MS2_Count": 0, "Candidate_Generation_Truncated": truncated,
        "Candidate_Generation_Reason": reason, **LOCALIZATION_FALSE, **FORMAL_FALSE,
    }
    return DinucleotideCandidateResult(groups, assignments, states, rejected, summary, truncated, reason)


def extract_target_candidates(candidates: Iterable[dict[str, Any]], low: float, high: float) -> list[dict[str, Any]]:
    """Generic range helper used only for optional display/validation filters."""
    return [dict(row) for row in candidates if row.get("Theoretical_mz") is not None and float(low) <= float(row["Theoretical_mz"]) <= float(high)]


def build_dinucleotide_audit(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Backward-compatible entry point delegated to the generic orchestrator."""
    from rna_masshunter.p1_sap_dinucleotide_interpretation import build_p1_sap_dinucleotide_audit
    return build_p1_sap_dinucleotide_audit(*args, **kwargs)
