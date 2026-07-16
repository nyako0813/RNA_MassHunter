"""Conservative shadow synthesis for complete-structure RNase MS/MS evidence.

This module integrates existing tables only. It does not inspect spectra, read
mzML, use configuration, rematch peaks, or modify formal results.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Any

from rna_masshunter.ms2_identity_evidence import physical_observed_peak_key

FORMAL_FALSE = {
    "Applied_To_Formal_Result": False,
    "Formal_Change_Ready": False,
    "Formal_Result_Changed": False,
}

SUMMARY_COLUMNS = [
    "Composite_Candidate_Count", "Composite_Peak_Evidence_Row_Count",
    "Positive_Assignment_Count", "Individual_Support_Assignment_Count",
    "Shared_Assignment_Count", "Identity_UNSUPPORTED_Count",
    "Identity_PROVISIONAL_CANDIDATE_SUPPORT_Count", "Identity_AMBIGUOUS_Count",
    "Localization_UNRESOLVED_Count", "Localization_POSITION_COMPATIBLE_Count",
    "Localization_PARTIALLY_SUPPORTED_Count", "Localization_AMBIGUOUS_Count",
    "Backbone_UNRESOLVED_Count", "Backbone_BACKBONE_COMPATIBLE_Count",
    "Backbone_PARTIALLY_SUPPORTED_Count", "Backbone_AMBIGUOUS_Count",
    "Structure_NOT_EVALUATED_Count", "Structure_UNRESOLVED_Count",
    "Structure_AMBIGUOUS_Count", "Structure_SUPPORTED_Count",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]

EVIDENCE_COLUMNS = [
    "Composite_Candidate_Key", "Candidate_ID", "Complete_Structure_ID",
    "Composite_Identity_Status", "Composite_Localization_Status",
    "Composite_Backbone_Status", "Composite_Structure_Status",
    "Composite_Ambiguity_Status", "Positive_Physical_Peak_Count",
    "Individual_Support_Physical_Peak_Count", "Shared_Physical_Peak_Count",
    "Position_Informative_Support_Peak_Count",
    "Backbone_Informative_Support_Peak_Count", "Competing_Candidate_Peak_Count",
    "Competing_Complete_Structure_Peak_Count", "Competing_Theoretical_Ion_Peak_Count",
    "Structural_Isomer_Sharing", "Composite_Support_Status",
    "Legacy_Comparison_Class", "Evidence_Basis", "Limiting_Reasons",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]

PEAK_EVIDENCE_COLUMNS = [
    "Composite_Match_ID", "Physical_Observed_Peak_Key", "Spectrum_ID",
    "Observed_Peak_Index", "Raw_Peak_Index", "RT", "Observed_mz",
    "Observed_Intensity", "Observed_Intensity_State", "Candidate_ID",
    "Complete_Structure_ID", "Ion_ID", "Parent_Fragment_ID", "Ion_Series",
    "Ion_Number", "Cleavage_Position", "Included_Modified_Positions",
    "Included_Backbone_Bonds", "Position_Informative", "Backbone_Informative",
    "Candidate_Specific", "Complete_Structure_Specific",
    "Theoretical_Ion_Specific", "Position_Specific", "Backbone_Bond_Specific",
    "Positive_Assignment", "Shared_Assignment", "Counts_For_Individual_Support",
    "Counts_For_Position_Support", "Counts_For_Backbone_Support",
    "Assignment_Rank", "Best_Assignment", "Within_Tolerance_Assignment_Count",
    "Competing_Candidate_Count", "Competing_Candidate_IDs",
    "Competing_Complete_Structure_Count", "Competing_Complete_Structure_IDs",
    "Competing_Theoretical_Ion_Count", "Competing_Ion_IDs", "Mass_Error_ppm",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]


@dataclass(frozen=True)
class RNaseMS2CompositeEvidenceSynthesisResult:
    summary_rows: list[dict[str, Any]]
    evidence_rows: list[dict[str, Any]]
    peak_rows: list[dict[str, Any]]

    @property
    def sheets(self) -> dict[str, list[dict[str, Any]]]:
        return {
            "RNase_MS2_Composite_Summary": self.summary_rows,
            "RNase_MS2_Composite_Evidence": self.evidence_rows,
            "RNase_MS2_Composite_Peak_Evidence": self.peak_rows,
        }

    def to_jsonable(self) -> dict[str, Any]:
        return {
            "summary": self.summary_rows[0] if self.summary_rows else {},
            "candidates": self.evidence_rows,
            "peaks": self.peak_rows,
        }


def _bool(value: Any) -> bool:
    return value is True or str(value or "").strip().lower() in {"1", "true", "yes", "on"}


def _float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _positive(row: dict[str, Any]) -> bool:
    intensity = _float(row.get("Observed_Intensity"))
    return intensity is not None and intensity > 0


def _key(row: dict[str, Any]) -> tuple[str, str]:
    return str(row.get("Candidate_ID") or ""), str(row.get("Complete_Structure_ID") or "")


def _text_key(key: tuple[str, str]) -> str:
    return "|".join(value or "NA" for value in key)


def _row_sort(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        str(row.get("Physical_Observed_Peak_Key") or ""),
        int(row.get("Assignment_Rank") or 0), str(row.get("Candidate_ID") or ""),
        str(row.get("Complete_Structure_ID") or ""), str(row.get("Ion_ID") or ""),
    )


def _ambiguity(types: set[str]) -> str:
    return "NONE" if not types else next(iter(types)) if len(types) == 1 else "MULTIPLE"


def _fallback_peak(row: dict[str, Any], index: int) -> dict[str, Any]:
    spectrum = str(row.get("Spectrum_ID") or "")
    observed = row.get("Observed_mz", "")
    physical = physical_observed_peak_key({
        "Spectrum_ID": spectrum, "Observed_mz": observed,
        "Observed_Intensity": row.get("Observed_Intensity"), "RT": row.get("RT"),
    })
    return {
        **row, "Composite_Match_ID": f"LEGACY_COMPOSITE_MATCH_{index:06d}",
        "Physical_Observed_Peak_Key": physical, "Ion_ID": row.get("Ion_ID", ""),
        "Observed_Peak_Index": "", "Raw_Peak_Index": "", "RT": row.get("RT", ""),
        "Observed_Intensity_State": "positive" if _positive(row) else "zero",
        "Candidate_Specific": row.get("Candidate_Discriminating", False),
        "Complete_Structure_Specific": row.get("Candidate_Discriminating", False),
        "Theoretical_Ion_Specific": row.get("Candidate_Discriminating", False),
        "Position_Specific": row.get("Candidate_Discriminating", False),
        "Backbone_Bond_Specific": row.get("Candidate_Discriminating", False),
        "Assignment_Rank": 1, "Best_Assignment": True,
        "Within_Tolerance_Assignment_Count": 1,
        "Competing_Candidate_Count": 0, "Competing_Candidate_IDs": "",
        "Competing_Complete_Structure_Count": 0,
        "Competing_Complete_Structure_IDs": "",
        "Competing_Theoretical_Ion_Count": 0, "Competing_Ion_IDs": "",
    }


def build_rnase_ms2_composite_evidence_synthesis(
    composite_ions: list[dict[str, Any]],
    composite_matches: list[dict[str, Any]],
    composite_assignment_competition: list[dict[str, Any]],
    composite_support: list[dict[str, Any]],
    legacy_composite_compare: list[dict[str, Any]],
) -> RNaseMS2CompositeEvidenceSynthesisResult:
    """Integrate complete-structure evidence without changing formal results."""
    ions = list(composite_ions or [])
    matches = list(composite_matches or [])
    competition = list(composite_assignment_competition or [])
    support = list(composite_support or [])
    comparisons = list(legacy_composite_compare or [])

    universe: set[tuple[str, str]] = set()
    for rows in (ions, matches, competition, support):
        universe.update(_key(row) for row in rows if any(_key(row)))

    source_peaks = competition or [_fallback_peak(row, i) for i, row in enumerate(matches, 1)]
    peak_rows: list[dict[str, Any]] = []
    by_candidate: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for source in sorted(source_peaks, key=_row_sort):
        positive = _positive(source)
        candidate_specific = _bool(source.get("Candidate_Specific"))
        structure_specific = _bool(source.get("Complete_Structure_Specific"))
        ion_specific = _bool(source.get("Theoretical_Ion_Specific"))
        individual = positive and candidate_specific and structure_specific and ion_specific
        position_support = (
            individual and _bool(source.get("Position_Informative"))
            and _bool(source.get("Position_Specific"))
        )
        backbone_support = (
            individual and _bool(source.get("Backbone_Informative"))
            and _bool(source.get("Backbone_Bond_Specific"))
        )
        row = {column: source.get(column, "") for column in PEAK_EVIDENCE_COLUMNS}
        row.update({
            "Positive_Assignment": positive, "Shared_Assignment": positive and not individual,
            "Counts_For_Individual_Support": individual,
            "Counts_For_Position_Support": position_support,
            "Counts_For_Backbone_Support": backbone_support,
            **FORMAL_FALSE,
        })
        peak_rows.append(row)
        by_candidate[_key(source)].append(row)

    support_by_key = {_key(row): row for row in support}
    compare_by_candidate: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in comparisons:
        compare_by_candidate[str(row.get("Candidate_ID") or "")].append(row)

    evidence_rows: list[dict[str, Any]] = []
    for key in sorted(universe):
        rows = by_candidate.get(key, [])
        positive = [row for row in rows if row["Positive_Assignment"]]
        individual = [row for row in rows if row["Counts_For_Individual_Support"]]
        shared = [row for row in positive if row["Shared_Assignment"]]
        physical = lambda selected: {str(row.get("Physical_Observed_Peak_Key") or "") for row in selected}
        position_peaks = physical([row for row in rows if row["Counts_For_Position_Support"]])
        backbone_peaks = physical([row for row in rows if row["Counts_For_Backbone_Support"]])
        candidate_comp = physical([row for row in positive if int(row.get("Competing_Candidate_Count") or 0) > 0])
        structure_comp = physical([row for row in positive if int(row.get("Competing_Complete_Structure_Count") or 0) > 0])
        ion_comp = physical([row for row in positive if int(row.get("Competing_Theoretical_Ion_Count") or 0) > 0])
        ambiguity_types = set()
        if candidate_comp: ambiguity_types.add("CANDIDATE_COMPETITION")
        if structure_comp: ambiguity_types.add("COMPLETE_STRUCTURE_COMPETITION")
        if ion_comp: ambiguity_types.add("THEORETICAL_ION_COMPETITION")
        comparison_classes = {
            str(row.get("Comparison_Class") or "") for row in compare_by_candidate.get(key[0], [])
        }
        isomer = any("ISOMER" in value.upper() for value in comparison_classes)
        if isomer: ambiguity_types.add("STRUCTURAL_ISOMER")
        has_competition = bool(candidate_comp or structure_comp or ion_comp)

        if not positive:
            identity = "UNSUPPORTED"
        elif has_competition:
            identity = "AMBIGUOUS"
        elif individual:
            identity = "PROVISIONAL_CANDIDATE_SUPPORT"
        else:
            identity = "AMBIGUOUS"

        def axis_status(peaks: set[str], compatible: str) -> str:
            if has_competition:
                return "AMBIGUOUS"
            if len(peaks) >= 2:
                return "PARTIALLY_SUPPORTED"
            if len(peaks) == 1:
                return compatible
            return "UNRESOLVED"

        localization = axis_status(position_peaks, "POSITION_COMPATIBLE")
        backbone = axis_status(backbone_peaks, "BACKBONE_COMPATIBLE")
        structure = (
            "NOT_EVALUATED" if not positive else
            "AMBIGUOUS" if isomer or structure_comp else "UNRESOLVED"
        )
        reasons = []
        if not positive: reasons.append("no_positive_composite_ms2_match")
        if shared: reasons.append("shared_assignment_excluded_from_individual_support")
        if not position_peaks: reasons.append("no_noncompeting_position_informative_support")
        if not backbone_peaks: reasons.append("no_noncompeting_backbone_informative_support")
        if has_competition: reasons.append("assignment_competition")
        if isomer: reasons.append("structural_isomer_sharing")
        reasons.append("complete_structure_support_not_enabled")
        basis = []
        if individual: basis.append("candidate_specific_positive_composite_fragment_match")
        if position_peaks: basis.append("position_informative_independent_physical_peak")
        if backbone_peaks: basis.append("backbone_informative_independent_physical_peak")
        evidence_rows.append({
            "Composite_Candidate_Key": _text_key(key), "Candidate_ID": key[0],
            "Complete_Structure_ID": key[1], "Composite_Identity_Status": identity,
            "Composite_Localization_Status": localization,
            "Composite_Backbone_Status": backbone,
            "Composite_Structure_Status": structure,
            "Composite_Ambiguity_Status": _ambiguity(ambiguity_types),
            "Positive_Physical_Peak_Count": len(physical(positive)),
            "Individual_Support_Physical_Peak_Count": len(physical(individual)),
            "Shared_Physical_Peak_Count": len(physical(shared)),
            "Position_Informative_Support_Peak_Count": len(position_peaks),
            "Backbone_Informative_Support_Peak_Count": len(backbone_peaks),
            "Competing_Candidate_Peak_Count": len(candidate_comp),
            "Competing_Complete_Structure_Peak_Count": len(structure_comp),
            "Competing_Theoretical_Ion_Peak_Count": len(ion_comp),
            "Structural_Isomer_Sharing": isomer,
            "Composite_Support_Status": support_by_key.get(key, {}).get("Support_Status", ""),
            "Legacy_Comparison_Class": ";".join(sorted(value for value in comparison_classes if value)),
            "Evidence_Basis": ";".join(basis), "Limiting_Reasons": ";".join(reasons),
            **FORMAL_FALSE,
        })

    identity_counts = Counter(row["Composite_Identity_Status"] for row in evidence_rows)
    localization_counts = Counter(row["Composite_Localization_Status"] for row in evidence_rows)
    backbone_counts = Counter(row["Composite_Backbone_Status"] for row in evidence_rows)
    structure_counts = Counter(row["Composite_Structure_Status"] for row in evidence_rows)
    summary = [{
        "Composite_Candidate_Count": len(evidence_rows),
        "Composite_Peak_Evidence_Row_Count": len(peak_rows),
        "Positive_Assignment_Count": sum(row["Positive_Assignment"] for row in peak_rows),
        "Individual_Support_Assignment_Count": sum(row["Counts_For_Individual_Support"] for row in peak_rows),
        "Shared_Assignment_Count": sum(row["Shared_Assignment"] for row in peak_rows),
        **{f"Identity_{status}_Count": identity_counts[status] for status in ("UNSUPPORTED", "PROVISIONAL_CANDIDATE_SUPPORT", "AMBIGUOUS")},
        **{f"Localization_{status}_Count": localization_counts[status] for status in ("UNRESOLVED", "POSITION_COMPATIBLE", "PARTIALLY_SUPPORTED", "AMBIGUOUS")},
        **{f"Backbone_{status}_Count": backbone_counts[status] for status in ("UNRESOLVED", "BACKBONE_COMPATIBLE", "PARTIALLY_SUPPORTED", "AMBIGUOUS")},
        **{f"Structure_{status}_Count": structure_counts[status] for status in ("NOT_EVALUATED", "UNRESOLVED", "AMBIGUOUS")},
        "Structure_SUPPORTED_Count": 0, **FORMAL_FALSE,
    }]
    return RNaseMS2CompositeEvidenceSynthesisResult(summary, evidence_rows, peak_rows)
