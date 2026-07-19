"""Shadow-only peak-family and delta-network analysis for intact SCIEX profiles.

This layer consumes the established intact peak detector result.  It does not
redetect peaks, calibrate masses, assign modifications, or affect formal output.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from hashlib import sha256
from math import isfinite
from types import MappingProxyType
from typing import Any, Iterable, Mapping

from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.intact_rna_average_mass import (
    AVERAGE_ATOMIC_MASSES,
    PROTON_MASS_DA,
    TheoreticalMassDefinition,
    calculate_average_neutral_mass_from_composition,
    calculate_intact_rna_average_mass,
)
from rna_masshunter.intact_rna_mass import (
    FivePrimeState,
    IntactRnaMassParameters,
)
from rna_masshunter.masses import MONOISOTOPIC_ATOMIC_MASSES
from rna_masshunter.structure_fragment import RNA_RESIDUE_COMPOSITIONS

ALGORITHM_VERSION = "sciex-intact-peak-family-v2"
SODIUM_MONOISOTOPIC_MASS_DA = 22.9897692820
POTASSIUM_MONOISOTOPIC_MASS_DA = 38.9637064864


class PeakQualityClass(str, Enum):
    MAJOR_SHARP = "MAJOR_SHARP"
    MAJOR_BROAD = "MAJOR_BROAD"
    MINOR_SHARP = "MINOR_SHARP"
    MINOR_BROAD = "MINOR_BROAD"
    SHOULDER_OR_OVERLAP = "SHOULDER_OR_OVERLAP"
    LOW_SUPPORT = "LOW_SUPPORT"


class DeltaMatchClass(str, Enum):
    STRICT = "STRICT"
    EXPLORATORY = "EXPLORATORY"
    NO_MATCH = "NO_MATCH"


class DeltaMassDefinition(str, Enum):
    AVERAGE_DELTA = "AVERAGE_DELTA"
    MONOISOTOPIC_DELTA = "MONOISOTOPIC_DELTA"
    EXACT_ION_DELTA = "EXACT_ION_DELTA"
    UNKNOWN = "UNKNOWN"


class DeltaComparisonRole(str, Enum):
    FAMILY_EDGE_REFERENCE = "FAMILY_EDGE_REFERENCE"
    OUTPUT_CONVENTION_DIAGNOSTIC_ONLY = "OUTPUT_CONVENTION_DIAGNOSTIC_ONLY"
    MASS_DEFINITION_MISMATCH_DIAGNOSTIC_ONLY = "MASS_DEFINITION_MISMATCH_DIAGNOSTIC_ONLY"
    UNKNOWN_MASS_DEFINITION_DIAGNOSTIC_ONLY = "UNKNOWN_MASS_DEFINITION_DIAGNOSTIC_ONLY"


class DeltaReferenceCategory(str, Enum):
    CCA_OR_TERMINAL_STATE = "CCA_OR_TERMINAL_STATE"
    OUTPUT_CONVENTION_DIAGNOSTIC = "OUTPUT_CONVENTION_DIAGNOSTIC"
    ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY = "ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY"
    KNOWN_RNA_MODIFICATION_DIAGNOSTIC_ONLY = "KNOWN_RNA_MODIFICATION_DIAGNOSTIC_ONLY"
    UNRESOLVED = "UNRESOLVED"


@dataclass(frozen=True)
class PeakFamilyParameters:
    maximum_major_peaks_per_profile: int = 50
    minimum_relative_apex_intensity: float = 0.01
    minimum_relative_integrated_intensity: float = 0.01
    major_relative_intensity: float = 0.05
    minimum_relative_prominence: float = 0.005
    minor_sharp_minimum_relative_prominence: float = 0.01
    major_minimum_relative_prominence: float = 0.02
    sharp_fwhm_max_da: float = 6.0
    minimum_peak_separation_da: float = 1.0
    maximum_peak_overlap_fraction: float = 0.50
    mass_range_coverage_bins: int = 10
    strict_delta_tolerance_da: float = 0.5
    exploratory_delta_tolerance_da: float = 1.0
    candidate_strict_tolerance_da: float = 1.0
    candidate_exploratory_tolerance_da: float = 5.0

    def validate(self) -> None:
        if self.maximum_major_peaks_per_profile < 1:
            raise ValueError("maximum_major_peaks_per_profile must be positive")
        for name in (
            "minimum_relative_apex_intensity",
            "minimum_relative_integrated_intensity",
            "major_relative_intensity",
            "minimum_relative_prominence",
            "minor_sharp_minimum_relative_prominence",
            "major_minimum_relative_prominence",
            "maximum_peak_overlap_fraction",
        ):
            value = float(getattr(self, name))
            if not 0 <= value <= 1:
                raise ValueError(f"{name} must be between zero and one")
        for name in (
            "sharp_fwhm_max_da",
            "minimum_peak_separation_da",
            "strict_delta_tolerance_da",
            "exploratory_delta_tolerance_da",
            "candidate_strict_tolerance_da",
            "candidate_exploratory_tolerance_da",
        ):
            if float(getattr(self, name)) <= 0:
                raise ValueError(f"{name} must be positive")
        if self.exploratory_delta_tolerance_da < self.strict_delta_tolerance_da:
            raise ValueError("exploratory delta tolerance must be >= strict tolerance")
        if self.candidate_exploratory_tolerance_da < self.candidate_strict_tolerance_da:
            raise ValueError("candidate exploratory tolerance must be >= strict tolerance")
        if self.mass_range_coverage_bins < 1:
            raise ValueError("mass_range_coverage_bins must be positive")


@dataclass(frozen=True, kw_only=True)
class ShadowSafeguards:
    shadow_analysis_only: bool = True
    mass_evidence_only: bool = True
    rna_identity_confirmed: bool = False
    target_rna_identity_confirmed_by_mass: bool = False
    co_captured_rna_excluded: bool = False
    modification_assigned: bool = False
    modification_composition_assigned: bool = False
    position_assigned: bool = False
    structure_assigned: bool = False
    cca_state_confirmed: bool = False
    terminal_state_confirmed: bool = False
    biological_cause_assigned: bool = False
    rnase_t_assigned: bool = False
    applied_to_formal_score: bool = False
    applied_to_ranking: bool = False
    applied_to_candidate_filtering: bool = False
    applied_to_final_consensus: bool = False


@dataclass(frozen=True, kw_only=True)
class PeakFamilyPeak(ShadowSafeguards):
    peak_id: str
    source_id: str
    measurement_id: str
    rna_identity: str
    apex_mass: float
    centroid_mass: float | None
    apex_intensity: float
    integrated_intensity: float | None
    relative_apex_intensity: float
    relative_integrated_intensity: float | None
    left_boundary_mass: float | None
    right_boundary_mass: float | None
    peak_width_da: float | None
    fwhm_da: float | None
    prominence: float
    relative_prominence: float
    sharpness_score: float | None
    nearest_peak_separation_da: float | None
    peak_overlap_fraction: float | None
    peak_detection_status: str
    peak_quality_class: PeakQualityClass
    selected_as_major_peak: bool
    possible_isotope_or_reconstruction_artifact: bool
    possible_shoulder: bool
    possible_duplicate_peak: bool
    possible_adduct: bool
    possible_output_convention_offset: bool


@dataclass(frozen=True, kw_only=True)
class PeakDeltaPair(ShadowSafeguards):
    peak_pair_id: str
    lower_peak_id: str
    higher_peak_id: str
    lower_apex_mass: float
    higher_apex_mass: float
    delta_mass_da: float
    observed_centroid_delta_da: float | None
    observed_delta_mass_definition: DeltaMassDefinition
    lower_relative_intensity: float
    higher_relative_intensity: float
    pair_min_relative_intensity: float
    lower_quality_class: PeakQualityClass
    higher_quality_class: PeakQualityClass
    pair_quality_class: str
    possible_isotope_or_reconstruction_artifact: bool
    possible_shoulder: bool
    possible_duplicate_peak: bool
    possible_adduct: bool
    possible_output_convention_offset: bool


@dataclass(frozen=True)
class DeltaReference:
    reference_id: str
    reference_name: str
    reference_category: DeltaReferenceCategory
    reference_delta_da: float
    signed_delta_da: float
    elemental_difference: str
    mass_constant_set: str
    diagnostic_only: bool
    delta_mass_definition: DeltaMassDefinition
    mass_definition_compatible: bool
    eligible_for_family_edge: bool
    comparison_role: DeltaComparisonRole
    modification_assigned: bool = False
    position_assigned: bool = False
    structure_assigned: bool = False


@dataclass(frozen=True, kw_only=True)
class PeakDeltaMatch(ShadowSafeguards):
    delta_match_id: str
    peak_pair_id: str
    observed_apex_delta_da: float
    observed_centroid_delta_da: float | None
    reference_delta_da: float | None
    apex_delta_error_da: float | None
    centroid_delta_error_da: float | None
    absolute_apex_delta_error_da: float | None
    absolute_centroid_delta_error_da: float | None
    delta_match_class: DeltaMatchClass
    reference_category: DeltaReferenceCategory
    reference_name: str
    reference_id: str
    observed_delta_mass_definition: DeltaMassDefinition
    reference_delta_mass_definition: DeltaMassDefinition
    mass_definition_compatible: bool
    eligible_for_family_edge: bool
    comparison_role: DeltaComparisonRole


@dataclass(frozen=True, kw_only=True)
class CandidateBaselineRelation(ShadowSafeguards):
    relation_id: str
    peak_id: str
    nearest_candidate_id: str
    nearest_reference_mode: TheoreticalMassDefinition
    nearest_theoretical_mass: float
    nearest_delta_da: float
    nearest_tolerance_class: str
    observed_output_species: str = "UNKNOWN"
    observed_output_species_confirmed: bool = False


@dataclass(frozen=True, kw_only=True)
class PeakFamily(ShadowSafeguards):
    peak_family_id: str
    member_peak_ids: tuple[str, ...]
    member_count: int
    lowest_mass: float
    highest_mass: float
    mass_span_da: float
    highest_intensity_peak_id: str
    highest_quality_peak_id: str
    supported_delta_reference_count: int
    strict_edge_count: int
    exploratory_edge_count: int
    unresolved_edge_count: int
    candidate_relation_ids: tuple[str, ...]
    hypotheses: tuple[str, ...]
    primary_biological_hypothesis: str
    hypothesis_confirmed: bool
    target_modification_isoform_possible: bool
    co_captured_rna_possible: bool
    cca_or_terminal_state_possible: bool
    adduct_or_reconstruction_possible: bool
    native_modifications_expected: bool


@dataclass(frozen=True)
class SciexIntactPeakFamilyResult:
    status: str
    reason: str
    parameters: PeakFamilyParameters
    peaks: tuple[PeakFamilyPeak, ...]
    selected_peaks: tuple[PeakFamilyPeak, ...]
    delta_pairs: tuple[PeakDeltaPair, ...]
    delta_references: tuple[DeltaReference, ...]
    delta_matches: tuple[PeakDeltaMatch, ...]
    candidate_relations: tuple[CandidateBaselineRelation, ...]
    families: tuple[PeakFamily, ...]
    observed_delta_mass_definition: DeltaMassDefinition = DeltaMassDefinition.AVERAGE_DELTA
    observed_output_species: str = "UNKNOWN"
    observed_output_species_confirmed: bool = False
    output_species_assigned: bool = False
    algorithm_version: str = ALGORITHM_VERSION


_DIAGNOSTIC_ATOMIC_MASSES = MappingProxyType({
    "H": MONOISOTOPIC_ATOMIC_MASSES["H"],
    "O": MONOISOTOPIC_ATOMIC_MASSES["O"],
    "Na": SODIUM_MONOISOTOPIC_MASS_DA,
    "K": POTASSIUM_MONOISOTOPIC_MASS_DA,
})
_AVERAGE_METAL_ATOMIC_MASSES = MappingProxyType({
    "Na": 22.98976928,
    "K": 39.0983,
})
_QUALITY_PRIORITY = {
    PeakQualityClass.MAJOR_SHARP: 0,
    PeakQualityClass.MAJOR_BROAD: 1,
    PeakQualityClass.MINOR_SHARP: 2,
    PeakQualityClass.MINOR_BROAD: 3,
    PeakQualityClass.SHOULDER_OR_OVERLAP: 4,
    PeakQualityClass.LOW_SUPPORT: 5,
}


def _float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if isfinite(result) else None


def _peak_rows(detection_result: Any) -> tuple[dict[str, Any], ...]:
    values = (
        detection_result.peak_rows()
        if hasattr(detection_result, "peak_rows")
        else getattr(detection_result, "peaks", ())
    )
    rows = []
    for value in values or ():
        if isinstance(value, Mapping):
            rows.append(dict(value))
        elif hasattr(value, "to_dict"):
            rows.append(dict(value.to_dict()))
    return tuple(rows)


def _detection_status(detection_result: Any) -> str:
    value = (
        detection_result.diagnostics_row()
        if hasattr(detection_result, "diagnostics_row")
        else getattr(detection_result, "diagnostics", {})
    )
    if hasattr(value, "to_dict"):
        value = value.to_dict()
    return str(value.get("Detection_Status") or "UNKNOWN") if isinstance(value, Mapping) else "UNKNOWN"


def _nearest_geometry(rows: tuple[dict[str, Any], ...], index: int) -> tuple[float | None, float | None]:
    current = rows[index]
    current_mass = float(current["Apex_Mass"])
    neighbors = []
    if index:
        neighbors.append(rows[index - 1])
    if index + 1 < len(rows):
        neighbors.append(rows[index + 1])
    if not neighbors:
        return None, None
    neighbor = min(neighbors, key=lambda row: abs(float(row["Apex_Mass"]) - current_mass))
    separation = abs(float(neighbor["Apex_Mass"]) - current_mass)
    current_fwhm = _float(current.get("FWHM_Da"))
    neighbor_fwhm = _float(neighbor.get("FWHM_Da"))
    if current_fwhm is None or neighbor_fwhm is None or min(current_fwhm, neighbor_fwhm) <= 0:
        return separation, None
    overlap_da = max(0.0, (current_fwhm + neighbor_fwhm) / 2.0 - separation)
    return separation, min(1.0, overlap_da / min(current_fwhm, neighbor_fwhm))


def _classify_peak(
    *,
    relative_apex: float,
    relative_area: float | None,
    relative_prominence: float,
    fwhm_da: float | None,
    possible_shoulder: bool,
    possible_duplicate: bool,
    overlap_fraction: float | None,
    parameters: PeakFamilyParameters,
) -> PeakQualityClass:
    if possible_shoulder or possible_duplicate or (
        overlap_fraction is not None and overlap_fraction >= parameters.maximum_peak_overlap_fraction
    ):
        return PeakQualityClass.SHOULDER_OR_OVERLAP
    area = relative_area if relative_area is not None else 0.0
    if (
        relative_apex < parameters.minimum_relative_apex_intensity
        or area < parameters.minimum_relative_integrated_intensity
        or relative_prominence < parameters.minimum_relative_prominence
    ):
        return PeakQualityClass.LOW_SUPPORT
    broad = fwhm_da is None or fwhm_da > parameters.sharp_fwhm_max_da
    major = (
        max(relative_apex, area) >= parameters.major_relative_intensity
        and relative_prominence >= parameters.major_minimum_relative_prominence
    )
    if major:
        return PeakQualityClass.MAJOR_BROAD if broad else PeakQualityClass.MAJOR_SHARP
    return PeakQualityClass.MINOR_BROAD if broad else PeakQualityClass.MINOR_SHARP


def calculate_peak_metrics(
    detection_result: Any,
    *,
    source_id: str,
    measurement_id: str,
    rna_identity: str,
    parameters: PeakFamilyParameters | None = None,
) -> tuple[PeakFamilyPeak, ...]:
    """Normalize established detector peaks and assign deterministic quality classes."""
    params = parameters or PeakFamilyParameters()
    params.validate()
    original_rows = _peak_rows(detection_result)
    rows = tuple(sorted(original_rows, key=lambda row: (float(row["Apex_Mass"]), str(row.get("Peak_ID")))))
    if not rows:
        return ()
    max_apex = max(float(row.get("Apex_Intensity_Raw") or 0.0) for row in rows)
    areas = [_float(row.get("Peak_Area_Baseline_Corrected")) for row in rows]
    max_area = max((value for value in areas if value is not None), default=0.0)
    max_prominence = max(float(row.get("Prominence") or 0.0) for row in rows)
    status = _detection_status(detection_result)
    interim = []
    for index, row in enumerate(rows):
        apex = float(row["Apex_Mass"])
        centroid = _float(row.get("Centroid_Mass"))
        intensity = float(row.get("Apex_Intensity_Raw") or 0.0)
        area = areas[index]
        prominence = float(row.get("Prominence") or 0.0)
        relative_apex = intensity / max_apex if max_apex > 0 else 0.0
        relative_area = area / max_area if area is not None and max_area > 0 else None
        relative_prominence = prominence / max_prominence if max_prominence > 0 else 0.0
        left = _float(row.get("Left_Boundary_Mass"))
        right = _float(row.get("Right_Boundary_Mass"))
        width = right - left if left is not None and right is not None else None
        fwhm = _float(row.get("FWHM_Da"))
        separation, overlap = _nearest_geometry(rows, index)
        possible_shoulder = bool(row.get("Possible_Shoulder", False))
        possible_duplicate = bool(
            separation is not None and separation < params.minimum_peak_separation_da
        ) or bool(overlap is not None and overlap >= params.maximum_peak_overlap_fraction)
        isotope_artifact = bool(
            separation is not None
            and separation <= 1.5
            and abs(separation * 2.0 - round(separation * 2.0)) <= 1e-9
        )
        quality = _classify_peak(
            relative_apex=relative_apex,
            relative_area=relative_area,
            relative_prominence=relative_prominence,
            fwhm_da=fwhm,
            possible_shoulder=possible_shoulder,
            possible_duplicate=possible_duplicate,
            overlap_fraction=overlap,
            parameters=params,
        )
        sharpness = relative_prominence / fwhm if fwhm is not None and fwhm > 0 else None
        interim.append(PeakFamilyPeak(
            peak_id=str(row.get("Peak_ID") or f"PEAK_{index + 1:05d}"),
            source_id=source_id,
            measurement_id=measurement_id,
            rna_identity=rna_identity,
            apex_mass=apex,
            centroid_mass=centroid,
            apex_intensity=intensity,
            integrated_intensity=area,
            relative_apex_intensity=relative_apex,
            relative_integrated_intensity=relative_area,
            left_boundary_mass=left,
            right_boundary_mass=right,
            peak_width_da=width,
            fwhm_da=fwhm,
            prominence=prominence,
            relative_prominence=relative_prominence,
            sharpness_score=sharpness,
            nearest_peak_separation_da=separation,
            peak_overlap_fraction=overlap,
            peak_detection_status=status,
            peak_quality_class=quality,
            selected_as_major_peak=False,
            possible_isotope_or_reconstruction_artifact=isotope_artifact,
            possible_shoulder=possible_shoulder,
            possible_duplicate_peak=possible_duplicate,
            possible_adduct=False,
            possible_output_convention_offset=False,
        ))
    selected_ids = {peak.peak_id for peak in select_major_peaks(tuple(interim), parameters=params)}
    return tuple(
        PeakFamilyPeak(**{
            **peak.__dict__,
            "selected_as_major_peak": peak.peak_id in selected_ids,
        })
        for peak in interim
    )


def _selection_rank(peak: PeakFamilyPeak) -> tuple[Any, ...]:
    return (
        _QUALITY_PRIORITY[peak.peak_quality_class],
        -peak.relative_apex_intensity,
        -(peak.relative_integrated_intensity or 0.0),
        -peak.relative_prominence,
        peak.apex_mass,
        peak.peak_id,
    )


def select_major_peaks(
    peaks: tuple[PeakFamilyPeak, ...],
    *,
    parameters: PeakFamilyParameters | None = None,
) -> tuple[PeakFamilyPeak, ...]:
    """Select quality-supported peaks with a deterministic coverage-aware bound."""
    params = parameters or PeakFamilyParameters()
    params.validate()
    eligible = [
        peak for peak in peaks
        if peak.peak_quality_class in {PeakQualityClass.MAJOR_SHARP, PeakQualityClass.MAJOR_BROAD}
        or (
            peak.peak_quality_class is PeakQualityClass.MINOR_SHARP
            and peak.relative_prominence >= params.minor_sharp_minimum_relative_prominence
        )
    ]
    eligible.sort(key=_selection_rank)
    maximum = params.maximum_major_peaks_per_profile
    if len(eligible) <= maximum:
        return tuple(sorted(eligible, key=lambda peak: (peak.apex_mass, peak.peak_id)))
    low_mass = min(peak.apex_mass for peak in eligible)
    high_mass = max(peak.apex_mass for peak in eligible)
    span = high_mass - low_mass
    covered: list[PeakFamilyPeak] = []
    if span > 0:
        for bin_index in range(params.mass_range_coverage_bins):
            left = low_mass + span * bin_index / params.mass_range_coverage_bins
            right = low_mass + span * (bin_index + 1) / params.mass_range_coverage_bins
            members = [
                peak for peak in eligible
                if left <= peak.apex_mass <= right
                and (bin_index + 1 == params.mass_range_coverage_bins or peak.apex_mass < right)
            ]
            if members:
                covered.append(min(members, key=_selection_rank))
    class_representatives = []
    for quality in (
        PeakQualityClass.MAJOR_SHARP,
        PeakQualityClass.MAJOR_BROAD,
        PeakQualityClass.MINOR_SHARP,
    ):
        members = [peak for peak in eligible if peak.peak_quality_class is quality]
        if members:
            class_representatives.append(min(members, key=_selection_rank))
    selected = list(dict.fromkeys((*class_representatives, *covered)))[:maximum]
    for peak in eligible:
        if peak not in selected:
            selected.append(peak)
        if len(selected) == maximum:
            break
    return tuple(sorted(selected, key=lambda peak: (peak.apex_mass, peak.peak_id)))


def generate_delta_pairs(peaks: tuple[PeakFamilyPeak, ...]) -> tuple[PeakDeltaPair, ...]:
    ordered = tuple(sorted(peaks, key=lambda peak: (peak.apex_mass, peak.peak_id)))
    pairs = []
    for lower_index, lower in enumerate(ordered):
        for higher in ordered[lower_index + 1:]:
            if lower.peak_id == higher.peak_id:
                raise ValueError("same peak cannot form a delta pair")
            centroid_delta = (
                higher.centroid_mass - lower.centroid_mass
                if lower.centroid_mass is not None and higher.centroid_mass is not None
                else None
            )
            shoulder = lower.possible_shoulder or higher.possible_shoulder
            duplicate = lower.possible_duplicate_peak or higher.possible_duplicate_peak
            artifact = (
                lower.possible_isotope_or_reconstruction_artifact
                or higher.possible_isotope_or_reconstruction_artifact
            )
            pair_quality = (
                "SHOULDER_OR_DUPLICATE_RISK" if shoulder or duplicate
                else "BROAD_MEMBER" if (
                    lower.peak_quality_class in {PeakQualityClass.MAJOR_BROAD, PeakQualityClass.MINOR_BROAD}
                    or higher.peak_quality_class in {PeakQualityClass.MAJOR_BROAD, PeakQualityClass.MINOR_BROAD}
                )
                else "SHARP_PAIR"
            )
            pair_id = "PAIR__" + sha256(
                f"{lower.peak_id}|{higher.peak_id}".encode("utf-8")
            ).hexdigest()[:16].upper()
            pairs.append(PeakDeltaPair(
                peak_pair_id=pair_id,
                lower_peak_id=lower.peak_id,
                higher_peak_id=higher.peak_id,
                lower_apex_mass=lower.apex_mass,
                higher_apex_mass=higher.apex_mass,
                delta_mass_da=higher.apex_mass - lower.apex_mass,
                observed_centroid_delta_da=centroid_delta,
                observed_delta_mass_definition=DeltaMassDefinition.AVERAGE_DELTA,
                lower_relative_intensity=lower.relative_apex_intensity,
                higher_relative_intensity=higher.relative_apex_intensity,
                pair_min_relative_intensity=min(
                    lower.relative_apex_intensity, higher.relative_apex_intensity
                ),
                lower_quality_class=lower.peak_quality_class,
                higher_quality_class=higher.peak_quality_class,
                pair_quality_class=pair_quality,
                possible_isotope_or_reconstruction_artifact=artifact,
                possible_shoulder=shoulder,
                possible_duplicate_peak=duplicate,
                possible_adduct=False,
                possible_output_convention_offset=False,
            ))
    return tuple(pairs)


def _average_transition(sequence_before: str, sequence_after: str) -> float:
    parameters = IntactRnaMassParameters(five_prime_state=FivePrimeState.OH)
    before = calculate_intact_rna_average_mass(sequence_before, parameters=parameters)
    after = calculate_intact_rna_average_mass(sequence_after, parameters=parameters)
    return after.average_neutral_molecular_mass_m - before.average_neutral_molecular_mass_m


def _reference(
    reference_id: str,
    name: str,
    category: DeltaReferenceCategory,
    signed_delta: float,
    elemental_difference: str,
    mass_constant_set: str,
    mass_definition: DeltaMassDefinition,
    compatible: bool,
    family_edge_eligible: bool,
    comparison_role: DeltaComparisonRole,
) -> DeltaReference:
    return DeltaReference(
        reference_id,
        name,
        category,
        abs(float(signed_delta)),
        float(signed_delta),
        elemental_difference,
        mass_constant_set,
        True,
        mass_definition,
        compatible,
        family_edge_eligible,
        comparison_role,
    )


def _average_reference(
    reference_id: str,
    name: str,
    category: DeltaReferenceCategory,
    signed_delta: float,
    elemental_difference: str,
    mass_constant_set: str = "AVERAGE_ATOMIC_MASSES",
) -> DeltaReference:
    return _reference(
        reference_id,
        name,
        category,
        signed_delta,
        elemental_difference,
        mass_constant_set,
        DeltaMassDefinition.AVERAGE_DELTA,
        True,
        True,
        DeltaComparisonRole.FAMILY_EDGE_REFERENCE,
    )


def _diagnostic_reference(
    reference_id: str,
    name: str,
    category: DeltaReferenceCategory,
    signed_delta: float,
    elemental_difference: str,
    mass_constant_set: str,
    mass_definition: DeltaMassDefinition,
    role: DeltaComparisonRole = DeltaComparisonRole.MASS_DEFINITION_MISMATCH_DIAGNOSTIC_ONLY,
) -> DeltaReference:
    return _reference(
        reference_id,
        name,
        category,
        signed_delta,
        elemental_difference,
        mass_constant_set,
        mass_definition,
        False,
        False,
        role,
    )


def _modification_mass_definition(modification: Any) -> DeltaMassDefinition:
    raw = getattr(modification, "raw", {})
    if not isinstance(raw, Mapping):
        return DeltaMassDefinition.UNKNOWN
    explicit = " ".join(
        str(raw.get(key) or "")
        for key in ("delta_mass_definition", "mass_definition", "mass_type")
    ).lower()
    if "average" in explicit:
        return DeltaMassDefinition.AVERAGE_DELTA
    if "mono" in explicit or "modified_nucleoside_mass_mono" in raw:
        return DeltaMassDefinition.MONOISOTOPIC_DELTA
    sources = raw.get("sources", ())
    if isinstance(sources, list) and any(
        "mono" in str(item.get("source_field") or "").lower()
        for item in sources if isinstance(item, Mapping)
    ):
        return DeltaMassDefinition.MONOISOTOPIC_DELTA
    return DeltaMassDefinition.UNKNOWN


def _modification_composition_delta(modification: Any) -> Mapping[str, int] | None:
    raw = getattr(modification, "raw", {})
    if not isinstance(raw, Mapping):
        return None
    for key in ("elemental_composition_delta", "composition_delta"):
        value = raw.get(key)
        if not isinstance(value, Mapping):
            continue
        counts: dict[str, int] = {}
        for element, count in value.items():
            symbol = str(element)
            if symbol not in AVERAGE_ATOMIC_MASSES or isinstance(count, bool):
                return None
            try:
                integer = int(count)
            except (TypeError, ValueError):
                return None
            if integer != count:
                return None
            counts[symbol] = integer
        return counts
    return None


def _average_signed_composition_mass(composition: Mapping[str, int]) -> float:
    return sum(AVERAGE_ATOMIC_MASSES[element] * count for element, count in composition.items())


def build_delta_reference_registry(
    known_modifications: Iterable[Any] = (),
) -> tuple[DeltaReference, ...]:
    """Build mass-definition-explicit references for an average-delta profile."""
    c_addition = _average_transition("A", "AC")
    a_addition = _average_transition("C", "CA")
    oh = calculate_intact_rna_average_mass(
        "A", parameters=IntactRnaMassParameters(five_prime_state=FivePrimeState.OH),
    )
    phosphate = calculate_intact_rna_average_mass(
        "A", parameters=IntactRnaMassParameters(five_prime_state=FivePrimeState.MONOPHOSPHATE),
    )
    c_formula = ElementalComposition(RNA_RESIDUE_COMPOSITIONS["C"]).canonical_string()
    a_formula = ElementalComposition(RNA_RESIDUE_COMPOSITIONS["A"]).canonical_string()
    average_h = AVERAGE_ATOMIC_MASSES["H"]
    average_o = AVERAGE_ATOMIC_MASSES["O"]
    average_water = calculate_average_neutral_mass_from_composition({"H": 2, "O": 1})
    mono_h = _DIAGNOSTIC_ATOMIC_MASSES["H"]
    mono_water = ElementalComposition({"H": 2, "O": 1}).exact_mass
    mono_oxygen = ElementalComposition({"O": 1}).exact_mass
    refs = [
        _average_reference("CCA_NONE_TO_C", "NONE_TO_C", DeltaReferenceCategory.CCA_OR_TERMINAL_STATE,
                           c_addition, c_formula),
        _average_reference("CCA_C_TO_CC", "C_TO_CC", DeltaReferenceCategory.CCA_OR_TERMINAL_STATE,
                           c_addition, c_formula),
        _average_reference("CCA_CC_TO_CCA", "CC_TO_CCA", DeltaReferenceCategory.CCA_OR_TERMINAL_STATE,
                           a_addition, a_formula),
        _average_reference("RESIDUE_C_ADDITION", "C_RESIDUE_ADDITION", DeltaReferenceCategory.CCA_OR_TERMINAL_STATE,
                           c_addition, c_formula),
        _average_reference("RESIDUE_A_ADDITION", "A_RESIDUE_ADDITION", DeltaReferenceCategory.CCA_OR_TERMINAL_STATE,
                           a_addition, a_formula),
        _average_reference(
            "TERMINAL_5P_VS_5OH", "5_PRIME_OH_TO_MONOPHOSPHATE",
            DeltaReferenceCategory.CCA_OR_TERMINAL_STATE,
            phosphate.average_neutral_molecular_mass_m - oh.average_neutral_molecular_mass_m,
            "H1O3P1",
        ),
        _diagnostic_reference(
            "OUTPUT_PLUS_PROTON", "+PROTON_MASS_DA",
            DeltaReferenceCategory.OUTPUT_CONVENTION_DIAGNOSTIC,
            PROTON_MASS_DA, "H+", "PROTON_MASS_DA", DeltaMassDefinition.EXACT_ION_DELTA,
            DeltaComparisonRole.OUTPUT_CONVENTION_DIAGNOSTIC_ONLY,
        ),
        _diagnostic_reference(
            "OUTPUT_MINUS_PROTON", "-PROTON_MASS_DA",
            DeltaReferenceCategory.OUTPUT_CONVENTION_DIAGNOSTIC,
            -PROTON_MASS_DA, "-H+", "PROTON_MASS_DA", DeltaMassDefinition.EXACT_ION_DELTA,
            DeltaComparisonRole.OUTPUT_CONVENTION_DIAGNOSTIC_ONLY,
        ),
        _average_reference(
            "NA_H_AVERAGE_DELTA", "NA_H_EXCHANGE_AVERAGE",
            DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY,
            _AVERAGE_METAL_ATOMIC_MASSES["Na"] - average_h, "Na1H-1",
            "AVERAGE_ATOMIC_MASSES_WITH_NA_K",
        ),
        _diagnostic_reference(
            "NA_H_EXACT_ION_DELTA", "NA_H_EXCHANGE_EXACT_ION",
            DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY,
            _DIAGNOSTIC_ATOMIC_MASSES["Na"] - mono_h, "Na1H-1",
            "MONOISOTOPIC_ATOMIC_MASSES", DeltaMassDefinition.EXACT_ION_DELTA,
        ),
        _average_reference(
            "K_H_AVERAGE_DELTA", "K_H_EXCHANGE_AVERAGE",
            DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY,
            _AVERAGE_METAL_ATOMIC_MASSES["K"] - average_h, "K1H-1",
            "AVERAGE_ATOMIC_MASSES_WITH_NA_K",
        ),
        _diagnostic_reference(
            "K_H_EXACT_ION_DELTA", "K_H_EXCHANGE_EXACT_ION",
            DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY,
            _DIAGNOSTIC_ATOMIC_MASSES["K"] - mono_h, "K1H-1",
            "MONOISOTOPIC_ATOMIC_MASSES", DeltaMassDefinition.EXACT_ION_DELTA,
        ),
        _average_reference(
            "WATER_AVERAGE_DELTA_HYDRATION", "HYDRATION_H2O_AVERAGE",
            DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY,
            average_water, "H2O1",
        ),
        _average_reference(
            "WATER_AVERAGE_DELTA_DEHYDRATION", "DEHYDRATION_H2O_AVERAGE",
            DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY,
            -average_water, "H-2O-1",
        ),
        _diagnostic_reference(
            "WATER_MONOISOTOPIC_DELTA_HYDRATION", "HYDRATION_H2O_MONOISOTOPIC",
            DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY,
            mono_water, "H2O1", "MONOISOTOPIC_ATOMIC_MASSES", DeltaMassDefinition.MONOISOTOPIC_DELTA,
        ),
        _diagnostic_reference(
            "WATER_MONOISOTOPIC_DELTA_DEHYDRATION", "DEHYDRATION_H2O_MONOISOTOPIC",
            DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY,
            -mono_water, "H-2O-1", "MONOISOTOPIC_ATOMIC_MASSES", DeltaMassDefinition.MONOISOTOPIC_DELTA,
        ),
        _average_reference(
            "O_AVERAGE_DELTA", "OXIDATION_O_AVERAGE",
            DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY,
            average_o, "O1",
        ),
        _diagnostic_reference(
            "O_MONOISOTOPIC_DELTA", "OXIDATION_O_MONOISOTOPIC",
            DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY,
            mono_oxygen, "O1", "MONOISOTOPIC_ATOMIC_MASSES", DeltaMassDefinition.MONOISOTOPIC_DELTA,
        ),
    ]
    for modification in known_modifications:
        mass = _float(getattr(modification, "mass_shift_from_unmodified", None))
        name = str(getattr(modification, "id", "") or getattr(modification, "symbol", ""))
        if mass is None or mass == 0 or not name:
            continue
        definition = _modification_mass_definition(modification)
        if definition is DeltaMassDefinition.AVERAGE_DELTA:
            refs.append(_average_reference(
                "KNOWN_MOD_AVERAGE__" + sha256(name.encode("utf-8")).hexdigest()[:16].upper(),
                name,
                DeltaReferenceCategory.KNOWN_RNA_MODIFICATION_DIAGNOSTIC_ONLY,
                mass,
                "NOT_ASSIGNED",
                "CURATED_AVERAGE_MODIFICATION_MASS_SOURCE",
            ))
        else:
            role = (
                DeltaComparisonRole.MASS_DEFINITION_MISMATCH_DIAGNOSTIC_ONLY
                if definition in {DeltaMassDefinition.MONOISOTOPIC_DELTA, DeltaMassDefinition.EXACT_ION_DELTA}
                else DeltaComparisonRole.UNKNOWN_MASS_DEFINITION_DIAGNOSTIC_ONLY
            )
            refs.append(_diagnostic_reference(
                "KNOWN_MOD__" + sha256(name.encode("utf-8")).hexdigest()[:16].upper(),
                name,
                DeltaReferenceCategory.KNOWN_RNA_MODIFICATION_DIAGNOSTIC_ONLY,
                mass,
                "NOT_ASSIGNED",
                "CURATED_MODIFICATION_MASS_SOURCE",
                definition,
                role,
            ))
        composition = _modification_composition_delta(modification)
        if composition is not None:
            average_mass = _average_signed_composition_mass(composition)
            refs.append(_average_reference(
                "KNOWN_MOD_COMPOSITION_AVERAGE__" + sha256(name.encode("utf-8")).hexdigest()[:16].upper(),
                name + "__COMPOSITION_AVERAGE",
                DeltaReferenceCategory.KNOWN_RNA_MODIFICATION_DIAGNOSTIC_ONLY,
                average_mass,
                ElementalComposition.delta(composition).canonical_string(),
                "AVERAGE_ATOMIC_MASSES_FROM_CURATED_COMPOSITION",
            ))
    return tuple(refs)

def match_delta_pairs(
    pairs: tuple[PeakDeltaPair, ...],
    references: tuple[DeltaReference, ...],
    *,
    parameters: PeakFamilyParameters | None = None,
) -> tuple[PeakDeltaMatch, ...]:
    params = parameters or PeakFamilyParameters()
    params.validate()
    matches = []
    for pair in pairs:
        retained = []
        for reference in references:
            compatible = (
                reference.mass_definition_compatible
                and reference.delta_mass_definition is pair.observed_delta_mass_definition
            )
            family_edge_eligible = reference.eligible_for_family_edge and compatible
            apex_error = pair.delta_mass_da - reference.reference_delta_da
            centroid_error = (
                pair.observed_centroid_delta_da - reference.reference_delta_da
                if pair.observed_centroid_delta_da is not None else None
            )
            absolute = abs(apex_error)
            if absolute <= params.strict_delta_tolerance_da:
                match_class = DeltaMatchClass.STRICT
            elif absolute <= params.exploratory_delta_tolerance_da:
                match_class = DeltaMatchClass.EXPLORATORY
            else:
                continue
            retained.append(PeakDeltaMatch(
                delta_match_id=f"{pair.peak_pair_id}__{reference.reference_id}",
                peak_pair_id=pair.peak_pair_id,
                observed_apex_delta_da=pair.delta_mass_da,
                observed_centroid_delta_da=pair.observed_centroid_delta_da,
                reference_delta_da=reference.reference_delta_da,
                apex_delta_error_da=apex_error,
                centroid_delta_error_da=centroid_error,
                absolute_apex_delta_error_da=absolute,
                absolute_centroid_delta_error_da=(
                    abs(centroid_error) if centroid_error is not None else None
                ),
                delta_match_class=match_class,
                reference_category=reference.reference_category,
                reference_name=reference.reference_name,
                reference_id=reference.reference_id,
                observed_delta_mass_definition=pair.observed_delta_mass_definition,
                reference_delta_mass_definition=reference.delta_mass_definition,
                mass_definition_compatible=compatible,
                eligible_for_family_edge=family_edge_eligible,
                comparison_role=reference.comparison_role,
            ))
        if retained:
            matches.extend(sorted(retained, key=lambda item: (
                item.absolute_apex_delta_error_da,
                item.reference_category.value,
                item.reference_name,
            )))
        else:
            matches.append(PeakDeltaMatch(
                delta_match_id=f"{pair.peak_pair_id}__NO_MATCH",
                peak_pair_id=pair.peak_pair_id,
                observed_apex_delta_da=pair.delta_mass_da,
                observed_centroid_delta_da=pair.observed_centroid_delta_da,
                reference_delta_da=None,
                apex_delta_error_da=None,
                centroid_delta_error_da=None,
                absolute_apex_delta_error_da=None,
                absolute_centroid_delta_error_da=None,
                delta_match_class=DeltaMatchClass.NO_MATCH,
                reference_category=DeltaReferenceCategory.UNRESOLVED,
                reference_name="UNRESOLVED",
                reference_id="",
                observed_delta_mass_definition=pair.observed_delta_mass_definition,
                reference_delta_mass_definition=DeltaMassDefinition.UNKNOWN,
                mass_definition_compatible=False,
                eligible_for_family_edge=False,
                comparison_role=DeltaComparisonRole.UNKNOWN_MASS_DEFINITION_DIAGNOSTIC_ONLY,
            ))
    return tuple(matches)



def _annotate_pair_diagnostics(
    pairs: tuple[PeakDeltaPair, ...],
    matches: tuple[PeakDeltaMatch, ...],
) -> tuple[PeakDeltaPair, ...]:
    categories_by_pair: dict[str, set[DeltaReferenceCategory]] = {}
    for match in matches:
        if match.delta_match_class is not DeltaMatchClass.NO_MATCH:
            categories_by_pair.setdefault(match.peak_pair_id, set()).add(match.reference_category)
    return tuple(
        replace(
            pair,
            possible_adduct=(
                DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY
                in categories_by_pair.get(pair.peak_pair_id, set())
            ),
            possible_output_convention_offset=(
                DeltaReferenceCategory.OUTPUT_CONVENTION_DIAGNOSTIC
                in categories_by_pair.get(pair.peak_pair_id, set())
            ),
        )
        for pair in pairs
    )


def connect_candidate_baselines(
    peaks: tuple[PeakFamilyPeak, ...],
    candidates: Iterable[Any],
    *,
    parameters: PeakFamilyParameters | None = None,
) -> tuple[CandidateBaselineRelation, ...]:
    params = parameters or PeakFamilyParameters()
    candidate_items = tuple(candidates)
    if not candidate_items:
        return ()
    modes = (
        (TheoreticalMassDefinition.AVERAGE_NEUTRAL_M, "theoretical_average_neutral_molecular_mass_m"),
        (TheoreticalMassDefinition.AVERAGE_M_PLUS_H, "theoretical_average_m_plus_h"),
        (TheoreticalMassDefinition.AVERAGE_M_MINUS_H, "theoretical_average_m_minus_h"),
        (TheoreticalMassDefinition.MONOISOTOPIC_NEUTRAL_M, "theoretical_monoisotopic_neutral_mass"),
    )
    result = []
    for peak in peaks:
        for mode, attribute in modes:
            nearest = min(
                candidate_items,
                key=lambda candidate: (
                    abs(peak.apex_mass - float(getattr(candidate, attribute))),
                    str(candidate.candidate_id),
                ),
            )
            theoretical = float(getattr(nearest, attribute))
            delta = peak.apex_mass - theoretical
            absolute = abs(delta)
            tolerance = (
                "STRICT" if absolute <= params.candidate_strict_tolerance_da
                else "EXPLORATORY" if absolute <= params.candidate_exploratory_tolerance_da
                else "NO_MATCH"
            )
            result.append(CandidateBaselineRelation(
                relation_id=f"{peak.peak_id}__{mode.value}",
                peak_id=peak.peak_id,
                nearest_candidate_id=str(nearest.candidate_id),
                nearest_reference_mode=mode,
                nearest_theoretical_mass=theoretical,
                nearest_delta_da=delta,
                nearest_tolerance_class=tolerance,
            ))
    return tuple(result)


def _family_hypotheses(
    rna_identity: str,
    categories: set[DeltaReferenceCategory],
) -> tuple[tuple[str, ...], str, bool, bool, bool, bool]:
    glu = rna_identity == "TRNA_GLU_UUC"
    hypotheses = ["TARGET_TRNA_MODIFICATION_ISOFORM"]
    cca = DeltaReferenceCategory.CCA_OR_TERMINAL_STATE in categories
    adduct = bool(categories & {
        DeltaReferenceCategory.ADDUCT_OR_CHEMICAL_STATE_DIAGNOSTIC_ONLY,
        DeltaReferenceCategory.OUTPUT_CONVENTION_DIAGNOSTIC,
    })
    if cca:
        hypotheses.append("CCA_OR_TERMINAL_STATE_SERIES")
    if adduct:
        hypotheses.append("ADDUCT_OR_OUTPUT_CONVENTION_SERIES")
    hypotheses.append("CO_CAPTURED_RNA_POSSIBLE")
    hypotheses.append("MIXED_OR_UNRESOLVED")
    primary = "TARGET_TRNA_MODIFICATION_ISOFORM" if glu else "MIXED_OR_UNRESOLVED"
    return tuple(hypotheses), primary, True, True, cca, adduct


def build_peak_families(
    peaks: tuple[PeakFamilyPeak, ...],
    pairs: tuple[PeakDeltaPair, ...],
    matches: tuple[PeakDeltaMatch, ...],
    candidate_relations: tuple[CandidateBaselineRelation, ...],
    *,
    rna_identity: str,
) -> tuple[PeakFamily, ...]:
    peak_by_id = {peak.peak_id: peak for peak in peaks}
    pair_by_id = {pair.peak_pair_id: pair for pair in pairs}
    supported_by_pair: dict[str, list[PeakDeltaMatch]] = {}
    for match in matches:
        if (
            match.delta_match_class in {DeltaMatchClass.STRICT, DeltaMatchClass.EXPLORATORY}
            and match.mass_definition_compatible
            and match.eligible_for_family_edge
        ):
            pair = pair_by_id[match.peak_pair_id]
            if pair.possible_duplicate_peak or pair.possible_shoulder:
                continue
            supported_by_pair.setdefault(match.peak_pair_id, []).append(match)
    adjacency = {peak.peak_id: set() for peak in peaks}
    for pair_id in supported_by_pair:
        pair = pair_by_id[pair_id]
        adjacency[pair.lower_peak_id].add(pair.higher_peak_id)
        adjacency[pair.higher_peak_id].add(pair.lower_peak_id)
    components = []
    unseen = set(adjacency)
    while unseen:
        start = min(unseen, key=lambda peak_id: (peak_by_id[peak_id].apex_mass, peak_id))
        stack = [start]
        member_ids = set()
        while stack:
            current = stack.pop()
            if current in member_ids:
                continue
            member_ids.add(current)
            unseen.discard(current)
            stack.extend(sorted(adjacency[current] - member_ids, reverse=True))
        components.append(tuple(sorted(member_ids, key=lambda peak_id: (
            peak_by_id[peak_id].apex_mass, peak_id
        ))))
    families = []
    relation_by_peak: dict[str, list[CandidateBaselineRelation]] = {}
    for relation in candidate_relations:
        relation_by_peak.setdefault(relation.peak_id, []).append(relation)
    for component in components:
        member_set = set(component)
        member_peaks = [peak_by_id[peak_id] for peak_id in component]
        internal_pairs = [
            pair for pair in pairs
            if pair.lower_peak_id in member_set and pair.higher_peak_id in member_set
        ]
        internal_pair_ids = {pair.peak_pair_id for pair in internal_pairs}
        internal_matches = [
            match for match in matches
            if match.peak_pair_id in internal_pair_ids
            and match.delta_match_class in {DeltaMatchClass.STRICT, DeltaMatchClass.EXPLORATORY}
            and match.mass_definition_compatible
            and match.eligible_for_family_edge
        ]
        categories = {match.reference_category for match in internal_matches}
        hypotheses, primary, target_possible, co_possible, cca_possible, adduct_possible = (
            _family_hypotheses(rna_identity, categories)
        )
        lowest = member_peaks[0].apex_mass
        highest = member_peaks[-1].apex_mass
        strongest = min(member_peaks, key=lambda peak: (-peak.apex_intensity, peak.apex_mass))
        quality_best = min(member_peaks, key=_selection_rank)
        strict_edges = len({
            match.peak_pair_id for match in internal_matches
            if match.delta_match_class is DeltaMatchClass.STRICT
        })
        exploratory_edges = len({
            match.peak_pair_id for match in internal_matches
            if match.delta_match_class is DeltaMatchClass.EXPLORATORY
            and not any(
                other.peak_pair_id == match.peak_pair_id
                and other.delta_match_class is DeltaMatchClass.STRICT
                for other in internal_matches
            )
        })
        unresolved = sum(
            not any(match.peak_pair_id == pair.peak_pair_id for match in internal_matches)
            for pair in internal_pairs
        )
        family_id = "FAMILY__" + sha256("|".join(component).encode("utf-8")).hexdigest()[:16].upper()
        relation_ids = tuple(
            relation.relation_id
            for peak_id in component
            for relation in relation_by_peak.get(peak_id, ())
        )
        families.append(PeakFamily(
            peak_family_id=family_id,
            member_peak_ids=component,
            member_count=len(component),
            lowest_mass=lowest,
            highest_mass=highest,
            mass_span_da=highest - lowest,
            highest_intensity_peak_id=strongest.peak_id,
            highest_quality_peak_id=quality_best.peak_id,
            supported_delta_reference_count=len({match.reference_id for match in internal_matches}),
            strict_edge_count=strict_edges,
            exploratory_edge_count=exploratory_edges,
            unresolved_edge_count=unresolved,
            candidate_relation_ids=relation_ids,
            hypotheses=hypotheses,
            primary_biological_hypothesis=primary,
            hypothesis_confirmed=False,
            target_modification_isoform_possible=target_possible,
            co_captured_rna_possible=co_possible,
            cca_or_terminal_state_possible=cca_possible,
            adduct_or_reconstruction_possible=adduct_possible,
            native_modifications_expected=True,
        ))
    return tuple(sorted(families, key=lambda family: (family.lowest_mass, family.peak_family_id)))


def analyze_sciex_intact_peak_families(
    detection_result: Any,
    *,
    source_id: str,
    measurement_id: str,
    rna_identity: str,
    candidates: Iterable[Any] = (),
    known_modifications: Iterable[Any] = (),
    parameters: PeakFamilyParameters | None = None,
) -> SciexIntactPeakFamilyResult:
    params = parameters or PeakFamilyParameters()
    params.validate()
    peaks = calculate_peak_metrics(
        detection_result,
        source_id=source_id,
        measurement_id=measurement_id,
        rna_identity=rna_identity,
        parameters=params,
    )
    selected = tuple(peak for peak in peaks if peak.selected_as_major_peak)
    pairs = generate_delta_pairs(selected)
    references = build_delta_reference_registry(known_modifications)
    matches = match_delta_pairs(pairs, references, parameters=params)
    pairs = _annotate_pair_diagnostics(pairs, matches)
    relations = connect_candidate_baselines(selected, candidates, parameters=params)
    families = build_peak_families(
        selected, pairs, matches, relations, rna_identity=rna_identity,
    )
    return SciexIntactPeakFamilyResult(
        status="COMPLETED",
        reason="SHADOW_PEAK_FAMILY_ANALYSIS",
        parameters=params,
        peaks=peaks,
        selected_peaks=selected,
        delta_pairs=pairs,
        delta_references=references,
        delta_matches=matches,
        candidate_relations=relations,
        families=families,
    )


def skipped_peak_family_result(
    reason: str,
    *,
    parameters: PeakFamilyParameters | None = None,
) -> SciexIntactPeakFamilyResult:
    params = parameters or PeakFamilyParameters()
    params.validate()
    return SciexIntactPeakFamilyResult(
        status="SKIPPED",
        reason=reason,
        parameters=params,
        peaks=(),
        selected_peaks=(),
        delta_pairs=(),
        delta_references=(),
        delta_matches=(),
        candidate_relations=(),
        families=(),
    )
