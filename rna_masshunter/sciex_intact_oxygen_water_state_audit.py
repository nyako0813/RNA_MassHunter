"""Shadow-only oxygen/water-equivalent state-series audit for intact profiles.

The audit consumes an existing :class:`SciexIntactPeakFamilyResult`; it never
redetects peaks and does not assign chemistry, structure, or reaction direction.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from hashlib import sha256
from typing import Iterable

from rna_masshunter.intact_rna_average_mass import TheoreticalMassDefinition
from rna_masshunter.sciex_intact_peak_family import (
    CandidateBaselineRelation,
    DeltaMassDefinition,
    DeltaMatchClass,
    DeltaReference,
    PeakDeltaMatch,
    PeakDeltaPair,
    PeakFamilyPeak,
    SciexIntactPeakFamilyResult,
)

ALGORITHM_VERSION = "sciex-intact-oxygen-water-state-audit-v1"
_O_AVERAGE_REFERENCE_ID = "O_AVERAGE_DELTA"
_WATER_AVERAGE_REFERENCE_ID = "WATER_AVERAGE_DELTA_HYDRATION"
_O_MONO_REFERENCE_ID = "O_MONOISOTOPIC_DELTA"
_WATER_MONO_REFERENCE_ID = "WATER_MONOISOTOPIC_DELTA_HYDRATION"


class StateRelationClass(str, Enum):
    O_EQUIVALENT_STRICT = "O_EQUIVALENT_STRICT"
    O_EQUIVALENT_EXPLORATORY = "O_EQUIVALENT_EXPLORATORY"
    H2O_EQUIVALENT_STRICT = "H2O_EQUIVALENT_STRICT"
    H2O_EQUIVALENT_EXPLORATORY = "H2O_EQUIVALENT_EXPLORATORY"
    OTHER_AVERAGE_COMPATIBLE = "OTHER_AVERAGE_COMPATIBLE"
    MASS_DEFINITION_MISMATCH_DIAGNOSTIC = "MASS_DEFINITION_MISMATCH_DIAGNOSTIC"
    NO_MATCH = "NO_MATCH"


class StateEdgeType(str, Enum):
    PLUS_O_EQUIVALENT = "PLUS_O_EQUIVALENT"
    PLUS_H2O_EQUIVALENT = "PLUS_H2O_EQUIVALENT"


class StateSeriesPattern(str, Enum):
    SINGLE_O_STEP = "SINGLE_O_STEP"
    MULTIPLE_SEQUENTIAL_O_STEPS = "MULTIPLE_SEQUENTIAL_O_STEPS"
    SINGLE_H2O_STEP = "SINGLE_H2O_STEP"
    MULTIPLE_SEQUENTIAL_H2O_STEPS = "MULTIPLE_SEQUENTIAL_H2O_STEPS"
    MIXED_H2O_AND_O_STEPS = "MIXED_H2O_AND_O_STEPS"
    BRANCHED_STATE_SERIES = "BRANCHED_STATE_SERIES"
    UNRESOLVED_STATE_SERIES = "UNRESOLVED_STATE_SERIES"


class PeakShapeSupportClass(str, Enum):
    STRONG_DISTINCT_PEAK_SUPPORT = "STRONG_DISTINCT_PEAK_SUPPORT"
    MIXED_PEAK_SUPPORT = "MIXED_PEAK_SUPPORT"
    WEAK_OR_OVERLAPPING_SUPPORT = "WEAK_OR_OVERLAPPING_SUPPORT"


@dataclass(frozen=True, kw_only=True)
class ChemicalStateSafeguards:
    shadow_analysis_only: bool = True
    mass_evidence_only: bool = True
    oxidation_assigned: bool = False
    hydration_assigned: bool = False
    dehydration_assigned: bool = False
    thioamide_assigned: bool = False
    thioamide_oxidation_state_assigned: bool = False
    modification_assigned: bool = False
    modification_composition_assigned: bool = False
    position_assigned: bool = False
    structure_assigned: bool = False
    reaction_direction_assigned: bool = False
    precursor_product_assigned: bool = False
    biological_cause_assigned: bool = False
    applied_to_formal_score: bool = False
    applied_to_ranking: bool = False
    applied_to_candidate_filtering: bool = False
    applied_to_final_consensus: bool = False


@dataclass(frozen=True, kw_only=True)
class StateRelation(ChemicalStateSafeguards):
    state_relation_id: str
    source_id: str
    lower_peak_id: str
    higher_peak_id: str
    lower_apex_mass: float
    higher_apex_mass: float
    observed_apex_delta_da: float
    observed_centroid_delta_da: float | None
    reference_name: str
    reference_id: str
    reference_delta_da: float | None
    apex_error_da: float | None
    centroid_error_da: float | None
    state_relation_class: StateRelationClass
    delta_match_class: DeltaMatchClass
    oxygen_equivalent: bool
    water_equivalent: bool
    oxygen_equivalent_mass_difference_detected: bool
    water_equivalent_mass_difference_detected: bool
    observed_delta_mass_definition: DeltaMassDefinition
    reference_delta_mass_definition: DeltaMassDefinition
    mass_definition_compatible: bool
    eligible_for_state_series: bool


@dataclass(frozen=True, kw_only=True)
class DirectedStateEdge(ChemicalStateSafeguards):
    state_edge_id: str
    state_relation_id: str
    source_id: str
    lower_peak_id: str
    higher_peak_id: str
    lower_apex_mass: float
    higher_apex_mass: float
    edge_type: StateEdgeType
    delta_match_class: DeltaMatchClass
    oxygen_equivalent_mass_difference_detected: bool
    water_equivalent_mass_difference_detected: bool
    observed_delta_mass_definition: DeltaMassDefinition
    reference_delta_mass_definition: DeltaMassDefinition
    mass_definition_compatible: bool
    eligible_for_state_series: bool
    oxidation_direction_assigned: bool = False


@dataclass(frozen=True)
class StateCandidateBaselineLink:
    peak_id: str
    nearest_candidate_id: str
    nearest_reference_mode: TheoreticalMassDefinition
    nearest_theoretical_mass: float
    nearest_delta_da: float
    nearest_tolerance_class: str


@dataclass(frozen=True, kw_only=True)
class StateSeries(ChemicalStateSafeguards):
    state_series_id: str
    source_id: str
    member_peak_ids: tuple[str, ...]
    member_apex_masses: tuple[float, ...]
    member_centroid_masses: tuple[float | None, ...]
    member_count: int
    lowest_mass: float
    highest_mass: float
    mass_span_da: float
    o_equivalent_edge_count: int
    h2o_equivalent_edge_count: int
    strict_edge_count: int
    exploratory_edge_count: int
    highest_intensity_peak_id: str
    highest_intensity_apex: float
    highest_relative_apex_intensity: float
    sequential_o_step_count: int
    sequential_h2o_step_count: int
    mixed_o_h2o_series: bool
    branched_series: bool
    series_pattern: StateSeriesPattern
    member_quality_classes: tuple[str, ...]
    member_relative_apex_intensities: tuple[float, ...]
    member_relative_integrated_intensities: tuple[float | None, ...]
    member_prominences: tuple[float, ...]
    member_fwhms: tuple[float | None, ...]
    member_peak_widths: tuple[float | None, ...]
    all_members_independent_peaks: bool
    any_shoulder: bool
    any_duplicate: bool
    peak_shape_support_class: PeakShapeSupportClass
    mass_ordered_apex_masses: tuple[float, ...]
    mass_ordered_relative_apex_intensities: tuple[float, ...]
    mass_ordered_relative_integrated_intensities: tuple[float | None, ...]
    candidate_baseline_links: tuple[StateCandidateBaselineLink, ...]
    oxygen_equivalent_mass_difference_detected: bool
    water_equivalent_mass_difference_detected: bool
    kinetic_order_assigned: bool = False
    pathway_order_assigned: bool = False
    oxidation_order_assigned: bool = False


@dataclass(frozen=True)
class OxygenWaterStateAuditResult:
    source_id: str
    status: str
    references: tuple[DeltaReference, ...]
    relations: tuple[StateRelation, ...]
    edges: tuple[DirectedStateEdge, ...]
    series: tuple[StateSeries, ...]
    algorithm_version: str = ALGORITHM_VERSION


def oxygen_water_reference_provenance(
    references: Iterable[DeltaReference],
) -> tuple[DeltaReference, ...]:
    """Return the existing average and monoisotopic O/H2O references unchanged."""
    wanted = {
        _O_AVERAGE_REFERENCE_ID, _WATER_AVERAGE_REFERENCE_ID,
        _O_MONO_REFERENCE_ID, _WATER_MONO_REFERENCE_ID,
    }
    selected = tuple(reference for reference in references if reference.reference_id in wanted)
    found = {reference.reference_id for reference in selected}
    missing = wanted - found
    if missing:
        raise ValueError(f"peak-family result lacks O/H2O references: {sorted(missing)}")
    return tuple(sorted(selected, key=lambda reference: reference.reference_id))


def _best_match(matches: tuple[PeakDeltaMatch, ...]) -> PeakDeltaMatch | None:
    matches = tuple(match for match in matches if match.delta_match_class is not DeltaMatchClass.NO_MATCH)
    priority = {
        _O_AVERAGE_REFERENCE_ID: 0,
        _WATER_AVERAGE_REFERENCE_ID: 1,
    }
    average_state = [match for match in matches if match.reference_id in priority]
    if average_state:
        return min(average_state, key=lambda match: (
            match.delta_match_class is DeltaMatchClass.EXPLORATORY,
            match.absolute_apex_delta_error_da or 0.0,
            priority[match.reference_id],
        ))
    compatible = [match for match in matches if match.mass_definition_compatible]
    if compatible:
        return min(compatible, key=lambda match: (
            match.delta_match_class is DeltaMatchClass.EXPLORATORY,
            match.absolute_apex_delta_error_da or 0.0,
            match.reference_id,
        ))
    if matches:
        return min(matches, key=lambda match: (
            match.delta_match_class is DeltaMatchClass.EXPLORATORY,
            match.absolute_apex_delta_error_da or 0.0,
            match.reference_id,
        ))
    return None


def _relation(pair: PeakDeltaPair, matches: tuple[PeakDeltaMatch, ...], source_id: str) -> StateRelation:
    match = _best_match(matches)
    if match is None:
        return StateRelation(
            state_relation_id="STATE_REL__" + pair.peak_pair_id,
            source_id=source_id, lower_peak_id=pair.lower_peak_id,
            higher_peak_id=pair.higher_peak_id, lower_apex_mass=pair.lower_apex_mass,
            higher_apex_mass=pair.higher_apex_mass,
            observed_apex_delta_da=pair.delta_mass_da,
            observed_centroid_delta_da=pair.observed_centroid_delta_da,
            reference_name="NO_MATCH", reference_id="NO_MATCH", reference_delta_da=None,
            apex_error_da=None, centroid_error_da=None,
            state_relation_class=StateRelationClass.NO_MATCH,
            delta_match_class=DeltaMatchClass.NO_MATCH, oxygen_equivalent=False,
            water_equivalent=False, oxygen_equivalent_mass_difference_detected=False,
            water_equivalent_mass_difference_detected=False,
            observed_delta_mass_definition=pair.observed_delta_mass_definition,
            reference_delta_mass_definition=DeltaMassDefinition.UNKNOWN,
            mass_definition_compatible=False, eligible_for_state_series=False,
        )
    oxygen = match.reference_id == _O_AVERAGE_REFERENCE_ID
    water = match.reference_id == _WATER_AVERAGE_REFERENCE_ID
    if oxygen:
        klass = (StateRelationClass.O_EQUIVALENT_STRICT if match.delta_match_class is DeltaMatchClass.STRICT
                 else StateRelationClass.O_EQUIVALENT_EXPLORATORY)
    elif water:
        klass = (StateRelationClass.H2O_EQUIVALENT_STRICT if match.delta_match_class is DeltaMatchClass.STRICT
                 else StateRelationClass.H2O_EQUIVALENT_EXPLORATORY)
    elif match.mass_definition_compatible:
        klass = StateRelationClass.OTHER_AVERAGE_COMPATIBLE
    else:
        klass = StateRelationClass.MASS_DEFINITION_MISMATCH_DIAGNOSTIC
    eligible = bool(
        (oxygen or water) and match.mass_definition_compatible
        and match.delta_match_class in {DeltaMatchClass.STRICT, DeltaMatchClass.EXPLORATORY}
        and not pair.possible_shoulder and not pair.possible_duplicate_peak
    )
    return StateRelation(
        state_relation_id="STATE_REL__" + match.delta_match_id,
        source_id=source_id, lower_peak_id=pair.lower_peak_id,
        higher_peak_id=pair.higher_peak_id, lower_apex_mass=pair.lower_apex_mass,
        higher_apex_mass=pair.higher_apex_mass, observed_apex_delta_da=pair.delta_mass_da,
        observed_centroid_delta_da=pair.observed_centroid_delta_da,
        reference_name=match.reference_name, reference_id=match.reference_id,
        reference_delta_da=match.reference_delta_da, apex_error_da=match.apex_delta_error_da,
        centroid_error_da=match.centroid_delta_error_da, state_relation_class=klass,
        delta_match_class=match.delta_match_class, oxygen_equivalent=oxygen,
        water_equivalent=water, oxygen_equivalent_mass_difference_detected=oxygen,
        water_equivalent_mass_difference_detected=water,
        observed_delta_mass_definition=match.observed_delta_mass_definition,
        reference_delta_mass_definition=match.reference_delta_mass_definition,
        mass_definition_compatible=match.mass_definition_compatible,
        eligible_for_state_series=eligible,
    )


def _pattern(edges: tuple[DirectedStateEdge, ...], branched: bool) -> StateSeriesPattern:
    if branched:
        return StateSeriesPattern.BRANCHED_STATE_SERIES
    o_count = sum(edge.edge_type is StateEdgeType.PLUS_O_EQUIVALENT for edge in edges)
    water_count = len(edges) - o_count
    if o_count and water_count:
        return StateSeriesPattern.MIXED_H2O_AND_O_STEPS
    if o_count == 1:
        return StateSeriesPattern.SINGLE_O_STEP
    if o_count > 1:
        return StateSeriesPattern.MULTIPLE_SEQUENTIAL_O_STEPS
    if water_count == 1:
        return StateSeriesPattern.SINGLE_H2O_STEP
    if water_count > 1:
        return StateSeriesPattern.MULTIPLE_SEQUENTIAL_H2O_STEPS
    return StateSeriesPattern.UNRESOLVED_STATE_SERIES


def _candidate_links(
    members: tuple[PeakFamilyPeak, ...], relations: tuple[CandidateBaselineRelation, ...],
) -> tuple[StateCandidateBaselineLink, ...]:
    result = []
    modes = tuple(TheoreticalMassDefinition)
    for peak in members:
        for mode in modes:
            available = [relation for relation in relations
                         if relation.peak_id == peak.peak_id and relation.nearest_reference_mode is mode]
            if not available:
                continue
            nearest = min(available, key=lambda relation: (
                abs(relation.nearest_delta_da), relation.nearest_candidate_id, relation.relation_id,
            ))
            result.append(StateCandidateBaselineLink(
                peak.peak_id, nearest.nearest_candidate_id, mode,
                nearest.nearest_theoretical_mass, nearest.nearest_delta_da,
                nearest.nearest_tolerance_class,
            ))
    return tuple(result)


def _build_series(
    peaks: tuple[PeakFamilyPeak, ...], edges: tuple[DirectedStateEdge, ...],
    candidate_relations: tuple[CandidateBaselineRelation, ...], source_id: str,
) -> tuple[StateSeries, ...]:
    by_id = {peak.peak_id: peak for peak in peaks}
    neighbors = {peak.peak_id: set() for peak in peaks}
    for edge in edges:
        neighbors[edge.lower_peak_id].add(edge.higher_peak_id)
        neighbors[edge.higher_peak_id].add(edge.lower_peak_id)
    components = []
    unseen = set(neighbors)
    while unseen:
        start = min(unseen, key=lambda peak_id: (by_id[peak_id].apex_mass, peak_id))
        stack, component = [start], set()
        while stack:
            current = stack.pop()
            if current in component:
                continue
            component.add(current); unseen.discard(current)
            stack.extend(sorted(neighbors[current] - component, reverse=True))
        components.append(tuple(sorted(component, key=lambda peak_id: (by_id[peak_id].apex_mass, peak_id))))
    result = []
    for component in components:
        members = tuple(by_id[peak_id] for peak_id in component)
        member_edges = tuple(edge for edge in edges
                             if edge.lower_peak_id in component and edge.higher_peak_id in component)
        degrees = {peak_id: 0 for peak_id in component}
        for edge in member_edges:
            degrees[edge.lower_peak_id] += 1; degrees[edge.higher_peak_id] += 1
        branched = any(degree > 2 for degree in degrees.values())
        o_count = sum(edge.edge_type is StateEdgeType.PLUS_O_EQUIVALENT for edge in member_edges)
        water_count = len(member_edges) - o_count
        strongest = max(members, key=lambda peak: (peak.apex_intensity, -peak.apex_mass))
        shoulder = any(peak.possible_shoulder for peak in members)
        duplicate = any(peak.possible_duplicate_peak for peak in members)
        independent = not shoulder and not duplicate
        if independent and all(peak.relative_prominence >= 0.02 for peak in members):
            support = PeakShapeSupportClass.STRONG_DISTINCT_PEAK_SUPPORT
        elif independent:
            support = PeakShapeSupportClass.MIXED_PEAK_SUPPORT
        else:
            support = PeakShapeSupportClass.WEAK_OR_OVERLAPPING_SUPPORT
        series_id = "STATE_SERIES__" + sha256("|".join(component).encode()).hexdigest()[:16].upper()
        result.append(StateSeries(
            state_series_id=series_id, source_id=source_id, member_peak_ids=component,
            member_apex_masses=tuple(peak.apex_mass for peak in members),
            member_centroid_masses=tuple(peak.centroid_mass for peak in members),
            member_count=len(members), lowest_mass=members[0].apex_mass,
            highest_mass=members[-1].apex_mass, mass_span_da=members[-1].apex_mass-members[0].apex_mass,
            o_equivalent_edge_count=o_count, h2o_equivalent_edge_count=water_count,
            strict_edge_count=sum(edge.delta_match_class is DeltaMatchClass.STRICT for edge in member_edges),
            exploratory_edge_count=sum(edge.delta_match_class is DeltaMatchClass.EXPLORATORY for edge in member_edges),
            highest_intensity_peak_id=strongest.peak_id, highest_intensity_apex=strongest.apex_mass,
            highest_relative_apex_intensity=strongest.relative_apex_intensity,
            sequential_o_step_count=o_count if not branched else 0,
            sequential_h2o_step_count=water_count if not branched else 0,
            mixed_o_h2o_series=bool(o_count and water_count), branched_series=branched,
            series_pattern=_pattern(member_edges, branched),
            member_quality_classes=tuple(peak.peak_quality_class.value for peak in members),
            member_relative_apex_intensities=tuple(peak.relative_apex_intensity for peak in members),
            member_relative_integrated_intensities=tuple(peak.relative_integrated_intensity for peak in members),
            member_prominences=tuple(peak.prominence for peak in members),
            member_fwhms=tuple(peak.fwhm_da for peak in members),
            member_peak_widths=tuple(peak.peak_width_da for peak in members),
            all_members_independent_peaks=independent, any_shoulder=shoulder, any_duplicate=duplicate,
            peak_shape_support_class=support,
            mass_ordered_apex_masses=tuple(peak.apex_mass for peak in members),
            mass_ordered_relative_apex_intensities=tuple(peak.relative_apex_intensity for peak in members),
            mass_ordered_relative_integrated_intensities=tuple(peak.relative_integrated_intensity for peak in members),
            candidate_baseline_links=_candidate_links(members, candidate_relations),
            oxygen_equivalent_mass_difference_detected=bool(o_count),
            water_equivalent_mass_difference_detected=bool(water_count),
        ))
    return tuple(sorted(result, key=lambda series: (series.lowest_mass, series.state_series_id)))


def audit_oxygen_water_state_series(
    peak_family_result: SciexIntactPeakFamilyResult,
) -> OxygenWaterStateAuditResult:
    """Audit O/H2O-equivalent relations using an existing peak-family result."""
    if peak_family_result.status != "COMPLETED":
        return OxygenWaterStateAuditResult("UNKNOWN", "SKIPPED", (), (), (), ())
    references = oxygen_water_reference_provenance(peak_family_result.delta_references)
    source_id = peak_family_result.selected_peaks[0].source_id if peak_family_result.selected_peaks else "UNKNOWN"
    matches_by_pair: dict[str, list[PeakDeltaMatch]] = {}
    for match in peak_family_result.delta_matches:
        matches_by_pair.setdefault(match.peak_pair_id, []).append(match)
    relations = tuple(
        _relation(pair, tuple(matches_by_pair.get(pair.peak_pair_id, ())), source_id)
        for pair in peak_family_result.delta_pairs
    )
    edges = tuple(DirectedStateEdge(
        state_edge_id="STATE_EDGE__" + relation.state_relation_id,
        state_relation_id=relation.state_relation_id, source_id=source_id,
        lower_peak_id=relation.lower_peak_id, higher_peak_id=relation.higher_peak_id,
        lower_apex_mass=relation.lower_apex_mass, higher_apex_mass=relation.higher_apex_mass,
        edge_type=(StateEdgeType.PLUS_O_EQUIVALENT if relation.oxygen_equivalent
                   else StateEdgeType.PLUS_H2O_EQUIVALENT),
        delta_match_class=relation.delta_match_class,
        oxygen_equivalent_mass_difference_detected=relation.oxygen_equivalent,
        water_equivalent_mass_difference_detected=relation.water_equivalent,
        observed_delta_mass_definition=relation.observed_delta_mass_definition,
        reference_delta_mass_definition=relation.reference_delta_mass_definition,
        mass_definition_compatible=relation.mass_definition_compatible,
        eligible_for_state_series=True,
    ) for relation in relations if relation.eligible_for_state_series)
    series = _build_series(
        peak_family_result.selected_peaks, edges,
        peak_family_result.candidate_relations, source_id,
    )
    return OxygenWaterStateAuditResult(source_id, "COMPLETED", references, relations, edges, series)
