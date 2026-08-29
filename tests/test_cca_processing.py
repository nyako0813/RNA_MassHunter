from __future__ import annotations

import pytest

from rna_masshunter.cca_processing import (
    CCAProcessingResult,
    CCAMaturationState,
    process_cca_tail,
)
from rna_masshunter.cca_tail_state import RegisteredSequenceCCAMode


def test_cca_is_unchanged():
    result = process_cca_tail(
        "ACGUCCA",
        RegisteredSequenceCCAMode.EXCLUDES_CCA,
    )

    assert result.processed_sequence == "ACGUCCA"
    assert result.original_tail_state is CCAMaturationState.CCA
    assert result.added_suffix == ""
    assert result.added_nucleotide_count == 0
    assert result.cca_processing_applied is False


def test_cc_gets_a():
    result = process_cca_tail(
        "ACGUCC",
        RegisteredSequenceCCAMode.EXCLUDES_CCA,
    )

    assert result.processed_sequence == "ACGUCCA"
    assert result.original_tail_state is CCAMaturationState.CC
    assert result.added_suffix == "A"
    assert result.added_nucleotide_count == 1
    assert result.cca_processing_applied is True


def test_c_gets_ca():
    result = process_cca_tail(
        "ACGUC",
        RegisteredSequenceCCAMode.EXCLUDES_CCA,
    )

    assert result.processed_sequence == "ACGUCCA"
    assert result.original_tail_state is CCAMaturationState.C
    assert result.added_suffix == "CA"
    assert result.added_nucleotide_count == 2
    assert result.cca_processing_applied is True


def test_non_c_terminal_gets_complete_cca():
    result = process_cca_tail(
        "ACGU",
        RegisteredSequenceCCAMode.EXCLUDES_CCA,
    )

    assert result.processed_sequence == "ACGUCCA"
    assert result.original_tail_state is CCAMaturationState.NONE
    assert result.added_suffix == "CCA"
    assert result.added_nucleotide_count == 3
    assert result.cca_processing_applied is True


def test_disabled_processing_preserves_sequence():
    result = process_cca_tail(
        "ACGUCC",
        RegisteredSequenceCCAMode.EXCLUDES_CCA,
        enabled=False,
    )

    assert result.processed_sequence == "ACGUCC"
    assert result.added_suffix == ""
    assert result.added_nucleotide_count == 0
    assert result.cca_processing_applied is False


def test_includes_complete_cca_does_not_append():
    result = process_cca_tail(
        "ACGUCCA",
        RegisteredSequenceCCAMode.INCLUDES_COMPLETE_CCA,
    )

    assert result.processed_sequence == "ACGUCCA"
    assert result.added_suffix == ""
    assert result.cca_processing_applied is False


def test_includes_complete_cca_without_cca_is_rejected():
    with pytest.raises(ValueError, match="must end with CCA"):
        process_cca_tail(
            "ACGUCC",
            RegisteredSequenceCCAMode.INCLUDES_COMPLETE_CCA,
        )


def test_unknown_mode_is_rejected():
    with pytest.raises(ValueError, match="UNKNOWN"):
        process_cca_tail(
            "ACGU",
            RegisteredSequenceCCAMode.UNKNOWN,
        )


def test_invalid_sequence_is_rejected():
    with pytest.raises(ValueError, match="invalid RNA base"):
        process_cca_tail(
            "ACGX",
            RegisteredSequenceCCAMode.EXCLUDES_CCA,
        )


def test_sequence_is_normalized_without_mutating_input():
    original = " acgucc "

    result = process_cca_tail(
        original,
        RegisteredSequenceCCAMode.EXCLUDES_CCA,
    )

    assert original == " acgucc "
    assert result.original_sequence == "ACGUCC"
    assert result.processed_sequence == "ACGUCCA"
    assert isinstance(result, CCAProcessingResult)
