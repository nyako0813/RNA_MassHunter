"""Generic enzyme-rule-driven cleavage-bond discovery for shadow audits."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any
from rna_masshunter.enzymes import find_cleavage_sites, get_enzyme_rule, normalize_enzyme_name

@dataclass(frozen=True)
class CleavageBondCandidate:
    sequence_id: str
    enzyme: str
    bond_id: str
    left_position: int
    right_position: int
    left_base: str
    right_base: str
    is_normal_cleavage_site: bool
    fragment_start: int
    fragment_end: int
    fragment_context: str


def missed_cleavage_fragment_range(sequence: str, enzyme: str, cleavage_position: int) -> tuple[int, int]:
    """Return the one-site missed-cleavage fragment spanning a normal cleavage bond."""
    seq = str(sequence or "").upper().replace("T", "U")
    sites = sorted(x for x in set(find_cleavage_sites(seq, enzyme)) if 0 < x <= len(seq))
    if cleavage_position not in sites or cleavage_position >= len(seq):
        raise ValueError("not_normal_cleavage_site")
    previous = max((x for x in sites if x < cleavage_position), default=0)
    following = min((x for x in sites if x > cleavage_position), default=len(seq))
    return previous + 1, following


def discover_candidate_cleavage_bonds(sequence: str, sequence_id: str, enzyme: str) -> list[CleavageBondCandidate]:
    seq = str(sequence or "").upper().replace("T", "U")
    enzyme_name = normalize_enzyme_name(enzyme)
    get_enzyme_rule(enzyme_name)
    candidates = []
    for position in sorted(set(find_cleavage_sites(seq, enzyme_name))):
        if position < 1 or position >= len(seq):
            continue
        start, end = missed_cleavage_fragment_range(seq, enzyme_name, position)
        candidates.append(CleavageBondCandidate(
            sequence_id, enzyme_name, f"{position}_{position + 1}", position, position + 1,
            seq[position - 1], seq[position], True, start, end, seq[start - 1:end],
        ))
    return candidates


def discovery_candidate_row(candidate: CleavageBondCandidate, *, nucleoside_state_count: int,
    backbone_state_count: int, pair_count: int, position_compatible: bool = True,
    pathway_compatible: bool = True, generated: bool = True, invalid_reason: str = "") -> dict[str, Any]:
    return {
        "Sequence_ID": candidate.sequence_id, "Enzyme": candidate.enzyme, "Bond_ID": candidate.bond_id,
        "Left_Position": candidate.left_position, "Right_Position": candidate.right_position,
        "Left_Base": candidate.left_base, "Right_Base": candidate.right_base,
        "Is_Normal_Cleavage_Site": candidate.is_normal_cleavage_site,
        "Fragment_Context": candidate.fragment_context,
        "Candidate_Nucleoside_State_Count": nucleoside_state_count,
        "Candidate_Backbone_State_Count": backbone_state_count, "Pair_Count": pair_count,
        "Position_Compatible": position_compatible, "Pathway_Compatible": pathway_compatible,
        "Generated_For_Evaluation": generated, "Invalid_Reason": invalid_reason,
        "Applied_To_Formal_Result": False, "Formal_Change_Ready": False, "Formal_Result_Changed": False,
    }
