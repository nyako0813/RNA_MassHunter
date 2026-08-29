"""Explicit 3-prime CCA maturation processing for registered RNA sequences."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from rna_masshunter.cca_tail_state import RegisteredSequenceCCAMode


class CCAMaturationState(str, Enum):
    NONE = "NONE"
    C = "C"
    CC = "CC"
    CCA = "CCA"


@dataclass(frozen=True)
class CCAProcessingResult:
    original_sequence: str
    processed_sequence: str
    original_tail_state: CCAMaturationState
    added_suffix: str
    added_nucleotide_count: int
    cca_processing_applied: bool


def _normalize_sequence(sequence: str) -> str:
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")

    normalized = "".join(sequence.split()).upper().replace("T", "U")

    if not normalized:
        raise ValueError("sequence must not be empty")

    allowed = set("ACGU")
    invalid = next((base for base in normalized if base not in allowed), None)

    if invalid is not None:
        raise ValueError(f"invalid RNA base: {invalid!r}")

    return normalized


def _tail_state(sequence: str) -> CCAMaturationState:
    if sequence.endswith("CCA"):
        return CCAMaturationState.CCA
    if sequence.endswith("CC"):
        return CCAMaturationState.CC
    if sequence.endswith("C"):
        return CCAMaturationState.C
    return CCAMaturationState.NONE


def process_cca_tail(
    sequence: str,
    registered_sequence_cca_mode: RegisteredSequenceCCAMode | str,
    *,
    enabled: bool = True,
) -> CCAProcessingResult:
    """Apply theoretical 3-prime CCA maturation without mutating the input."""
    sequence = _normalize_sequence(sequence)

    try:
        mode = (
            registered_sequence_cca_mode
            if isinstance(registered_sequence_cca_mode, RegisteredSequenceCCAMode)
            else RegisteredSequenceCCAMode(registered_sequence_cca_mode)
        )
    except (TypeError, ValueError) as exc:
        allowed = ", ".join(item.value for item in RegisteredSequenceCCAMode)
        raise ValueError(
            f"unknown registered_sequence_cca_mode: "
            f"{registered_sequence_cca_mode!r}; "
            f"expected one of {allowed}"
        ) from exc

    if mode is RegisteredSequenceCCAMode.UNKNOWN:
        raise ValueError("UNKNOWN registered CCA mode cannot be processed")

    state = _tail_state(sequence)

    # 登録配列そのものに完全なCCAが含まれる場合、
    # CCA処理による追加は行わない。
    if mode is RegisteredSequenceCCAMode.INCLUDES_COMPLETE_CCA:
        if state is not CCAMaturationState.CCA:
            raise ValueError(
                "INCLUDES_COMPLETE_CCA registered sequence must end with CCA"
            )

        return CCAProcessingResult(
            original_sequence=sequence,
            processed_sequence=sequence,
            original_tail_state=state,
            added_suffix="",
            added_nucleotide_count=0,
            cca_processing_applied=False,
        )

    if not enabled:
        return CCAProcessingResult(
            original_sequence=sequence,
            processed_sequence=sequence,
            original_tail_state=state,
            added_suffix="",
            added_nucleotide_count=0,
            cca_processing_applied=False,
        )

    suffix_by_state = {
        CCAMaturationState.CCA: "",
        CCAMaturationState.CC: "A",
        CCAMaturationState.C: "CA",
        CCAMaturationState.NONE: "CCA",
    }

    suffix = suffix_by_state[state]

    return CCAProcessingResult(
        original_sequence=sequence,
        processed_sequence=sequence + suffix,
        original_tail_state=state,
        added_suffix=suffix,
        added_nucleotide_count=len(suffix),
        cca_processing_applied=bool(suffix),
    )
