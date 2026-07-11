from dataclasses import dataclass, field
from typing import Any


@dataclass
class RunConfig:
    project: dict[str, Any] = field(default_factory=dict)
    input: dict[str, Any] = field(default_factory=dict)
    organism: dict[str, Any] = field(default_factory=dict)
    sequence: dict[str, Any] = field(default_factory=dict)
    experiment: dict[str, Any] = field(default_factory=dict)
    instrument: dict[str, Any] = field(default_factory=dict)
    reconstruction: dict[str, Any] = field(default_factory=dict)
    digestion: dict[str, Any] = field(default_factory=dict)
    alkaline_phosphatase: dict[str, Any] = field(default_factory=dict)
    fragment_mapping: dict[str, Any] = field(default_factory=dict)
    modification_search: dict[str, Any] = field(default_factory=dict)
    peak_filtering: dict[str, Any] = field(default_factory=dict)
    p1_annotation: dict[str, Any] = field(default_factory=dict)
    ms2_annotation: dict[str, Any] = field(default_factory=dict)
    modification_evidence_ranking: dict[str, Any] = field(default_factory=dict)
    biological_context: dict[str, Any] = field(default_factory=dict)
    performance: dict[str, Any] = field(default_factory=dict)
    reporting: dict[str, Any] = field(default_factory=dict)
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass
class Modification:
    id: str
    symbol: str | None
    mass_shift_from_unmodified: float
    category: str
    target_bases: list[str]
    detectability: Any = None
    curation: Any = None
    sources: Any = None
    source: Any = None
    source_priority: Any = None
    curation_status: str = ""
    candidate_policy: dict[str, Any] = field(default_factory=dict)
    chemical_group: str = ""
    near_isobaric_group: str = ""
    raw: dict[str, Any] = field(default_factory=dict)


@dataclass
class Peak:
    mz: float
    intensity: float
    rt: float | None = None
    scan_id: str | None = None
    ms_level: int = 1
    tier: str | None = None


@dataclass
class PeakTierResult:
    major: list[Peak] = field(default_factory=list)
    minor: list[Peak] = field(default_factory=list)
    trace: list[Peak] = field(default_factory=list)
    below_threshold: list[Peak] = field(default_factory=list)

    @property
    def usable_peaks(self) -> list[Peak]:
        return self.major + self.minor + self.trace


@dataclass
class IntactMassCandidate:
    observed_mass: float
    charge_state_count: int
    charge_states: list[int]
    supporting_peak_count: int
    total_intensity: float
    theoretical_mass: float | None = None
    mass_error_da: float | None = None
    mass_error_ppm: float | None = None
    assignment: str = "Unassigned intact mass"
    confidence: str = "Low"
    warnings: str = ""
    cluster_id: str = ""
    reconstruction_status: str = ""
    reconstruction_confidence: str = ""
    num_supporting_charge_states: int = 0
    charge_state_range: str = ""
    charge_state_continuity: str = ""
    neutral_mass_sd: float | None = None
    neutral_mass_range: float | None = None
    max_mass_error_ppm: float | None = None
    envelope_internal_error_ppm: float | None = None
    unmodified_theory_delta_da: float | None = None
    unmodified_theory_delta_ppm: float | None = None
    best_reference_label: str = ""
    best_reference_mass_da: float | None = None
    reference_mass_error_da: float | None = None
    reference_mass_error_ppm: float | None = None
    reference_mass_matched: bool = False
    in_neutral_mass_search_range: bool = False
    neutral_mass_search_min_da: float | None = None
    neutral_mass_search_max_da: float | None = None
    neutral_mass_range_status: str = ""
    in_target_review_mass_range: bool = False
    target_review_mass_range_status: str = ""
    target_review_priority: str = ""
    envelope_qc_eligible: bool = False
    intact_review_eligible: bool = False
    intact_strict_eligible: bool = False
    intact_envelope_qc_score: float | None = None
    intact_envelope_qc_rank: int | None = None
    strict_eligible_rank: int | None = None
    review_eligible_rank: int | None = None
    dominant_intact_envelope_flag: bool = False
    supporting_peak_ids: str = ""
    supporting_scan_ids: str = ""
    supporting_rt_values: str = ""
    supporting_charge_states: str = ""
    exact_peak_set_key: str = ""
    exact_duplicate_group_id: str = ""
    exact_duplicate_count: int = 1
    is_exact_duplicate_representative: bool = True
    intact_envelope_group_id: str = ""
    envelope_group_size: int = 1
    group_representative: bool = True
    group_ambiguity_status: str = "unique"
    comparison_representative: bool = False
    comparison_representative_reason: str = ""
    comparison_representative_rank: int | None = None
    excluded_from_comparison_reason: str = ""
    target_review_group_representative: bool = False
    target_review_rank: int | None = None
    dominant_target_review_eligible_flag: bool = False
    rt_min: float | None = None
    rt_max: float | None = None
    rt_mean: float | None = None
    rt_range_min: float | None = None
    max_rt_difference_min: float | None = None
    rt_consistency: str = ""
    total_supporting_intensity: float = 0.0
    mean_supporting_intensity: float = 0.0
    max_supporting_intensity: float = 0.0
    relative_envelope_intensity_percent: float = 0.0
    relative_overall_envelope_intensity_percent: float = 0.0
    relative_in_range_raw_intensity_percent: float = 0.0
    relative_intact_eligible_intensity_percent: float = 0.0
    supporting_peak_classes: str = ""
    trace_only_envelope: bool = False
    competing_envelope_count: int = 0
    limiting_factors: str = ""
    severe_limiting_factors: str = ""
    num_limiting_factors: int = 0
    primary_limiting_factor: str = ""
    comparison_ready_strict: bool = False
    comparison_ready_review: bool = False
    comparison_ready: bool = False
    comparison_readiness_reason: str = ""


@dataclass
class Fragment:
    fragment_id: str
    target_id: str
    sequence: str
    start: int
    end: int
    standard_start: int | None
    standard_end: int | None
    enzyme: str
    missed_cleavages: int
    terminal_form: str
    unmodified_mass: float
    warnings: list[str] = field(default_factory=list)


@dataclass
class FragmentMS1Match:
    match_id: str
    fragment_id: str
    target_id: str
    sequence: str
    start: int
    end: int
    standard_start: int | None
    standard_end: int | None
    enzyme: str
    missed_cleavages: int
    terminal_form: str
    fragment_mass: float
    charge: int
    theoretical_mz: float
    observed_mz: float
    mass_error_da: float
    mass_error_ppm: float
    intensity: float
    rt: float | None
    scan_id: str | None
    peak_tier: str | None
    confidence: str
    warnings: list[str] = field(default_factory=list)


@dataclass
class MS2SpectrumInfo:
    spectrum_id: str
    scan_index: int
    rt: float | None
    precursor_mz: float | None
    precursor_charge: int | None
    precursor_intensity: float | None
    num_peaks: int
    base_peak_mz: float | None
    base_peak_intensity: float | None
    total_ion_current: float
    peaks: list[tuple[float, float]] = field(default_factory=list)


@dataclass
class TheoreticalMS2Ion:
    ion_id: str
    parent_fragment_id: str
    parent_sequence: str
    ion_type: str
    ion_sequence: str
    ion_start: int
    ion_end: int
    charge: int
    theoretical_mz: float
    theoretical_mass: float
    neutral_loss: str = ""
    modification_id: str = ""
    modification_name: str = ""
    comment: str = ""


@dataclass
class MS2IonMatch:
    spectrum_id: str
    scan_index: int
    rt: float | None
    precursor_mz: float | None
    precursor_charge: int | None
    observed_mz: float
    observed_intensity: float
    best_ion_id: str
    best_ion_type: str
    best_ion_sequence: str
    parent_fragment_id: str
    parent_fragment_sequence: str
    ion_charge: int
    theoretical_mz: float
    mass_error_da: float
    mass_error_ppm: float
    match_status: str
    confidence: str
    alternative_candidates: str = ""
    comment: str = ""


@dataclass
class KnownModificationCandidate:
    candidate_id: str
    source_type: str
    source_id: str
    target_id: str
    sequence: str
    start: int | None
    end: int | None
    standard_start: int | None
    standard_end: int | None
    observed_mz: float | None
    theoretical_mz: float | None
    observed_mass: float
    unmodified_mass: float
    mass_error_unmodified_da: float
    mass_error_unmodified_ppm: float
    modification_id: str
    modification_symbol: str | None
    modification_name: str
    target_base: str
    modification_mass_shift: float
    modified_mass: float
    mass_error_modified_da: float
    mass_error_modified_ppm: float
    charge: int | None
    intensity: float
    rt: float | None
    peak_tier: str | None
    confidence: str
    priority_score: float
    notes: str = ""
    warnings: list[str] = field(default_factory=list)
