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
