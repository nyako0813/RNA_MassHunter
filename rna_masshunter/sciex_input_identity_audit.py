"""Shadow audit of SCIEX source-name identity against configured RNA metadata."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping
import re
import unicodedata

AUDIT_RESULT_KEY = "sciex_input_identity_audit"
SHEET_NAME = "SCIEX_Input_Identity_Audit"
WARNING_CODE = "SCIEX_INPUT_IDENTITY_CONFLICT"
ERROR_CODE = "SCIEX_INPUT_IDENTITY_AUDIT_ERROR"
AMINO_ACID_CODES = frozenset({
    "ala", "arg", "asn", "asp", "cys", "gln", "glu", "gly", "his", "ile",
    "leu", "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val",
})
STOP_TOKENS = frozenset({
    "full", "profile", "spectrum", "mass", "intact", "wt", "wildtype",
    "wild", "type", "trna", "rna", "sample", "run", "rep", "replicate", "r",
})
REPLICATE_TOKEN = re.compile(r"^(?:rep(?:licate)?|run|r)\d+$")
ANTICODON_TOKEN = re.compile(r"^[acgu]{3}$")
OUTPUT_COLUMNS = [
    "Audit_Status", "Audit_Eligible", "SCIEX_Source_File", "SCIEX_Source_Basename",
    "SCIEX_File_Stem", "SCIEX_Parent_Directory", "Parser_Sample_Name",
    "SCIEX_Filename_Tokens", "Configured_Sequence_Name", "Configured_Sequence",
    "Configured_Anticodon", "Configured_Organism_Group", "Configured_Species",
    "Configured_Condition_Name", "SCIEX_Amino_Acid_Tokens", "SCIEX_Anticodon_Tokens",
    "SCIEX_Combined_Identity_Tokens", "Configured_Amino_Acid_Tokens",
    "Configured_Anticodon_Tokens", "Configured_Combined_Identity_Tokens",
    "Amino_Acid_Match", "Anticodon_Match", "Combined_Identity_Match",
    "Identity_Conflict", "Identity_Evidence_Level", "Warning_Code", "Warning_Message",
    "Biological_Interpretation_Eligible", "Shadow_Only", "Applied_To_Formal_Score",
    "Applied_To_Ranking", "Applied_To_Candidate_Filtering",
    "Molecular_Identity_Assigned", "Notes",
]


@dataclass(frozen=True)
class IdentityTokens:
    filename_tokens: tuple[str, ...]
    amino_acids: tuple[str, ...]
    anticodons: tuple[str, ...]
    combined: tuple[str, ...]


@dataclass(frozen=True)
class SciexInputIdentityAuditResult:
    values: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "values", MappingProxyType(dict(self.values)))

    def row(self) -> dict[str, Any]:
        return dict(self.values)


def _normalized_words(value: Any) -> list[str]:
    normalized = unicodedata.normalize("NFKC", str(value or "")).casefold()
    return [token for token in re.split(r"[^a-z0-9]+", normalized) if token]


def _is_noise_token(token: str) -> bool:
    return (
        token in STOP_TOKENS
        or token.isdigit()
        or bool(REPLICATE_TOKEN.fullmatch(token))
    )


def _identity_from_words(words: list[str], *, filename: bool) -> IdentityTokens:
    meaningful = [token for token in words if not _is_noise_token(token)]
    amino: set[str] = set()
    anticodons: set[str] = set()
    combined: set[str] = set()

    for token in meaningful:
        if token in AMINO_ACID_CODES:
            amino.add(token)
        if ANTICODON_TOKEN.fullmatch(token):
            anticodons.add(token)
        if len(token) == 6 and token[:3] in AMINO_ACID_CODES and ANTICODON_TOKEN.fullmatch(token[3:]):
            amino.add(token[:3])
            anticodons.add(token[3:])
            combined.add(token)

    for left, right in zip(meaningful, meaningful[1:]):
        if left in AMINO_ACID_CODES and ANTICODON_TOKEN.fullmatch(right):
            amino.add(left)
            anticodons.add(right)
            combined.add(left + right)

    # Config metadata may provide the amino acid in the name and anticodon separately.
    if not filename and amino and anticodons:
        combined.update(left + right for left in amino for right in anticodons)

    candidate_tokens = set(meaningful)
    candidate_tokens.update(amino)
    candidate_tokens.update(anticodons)
    candidate_tokens.update(combined)
    return IdentityTokens(
        tuple(sorted(candidate_tokens)), tuple(sorted(amino)),
        tuple(sorted(anticodons)), tuple(sorted(combined)),
    )


def extract_sciex_identity_tokens(source_path: str | Path) -> IdentityTokens:
    source = Path(source_path)
    return _identity_from_words(_normalized_words(source.stem), filename=True)


def extract_configured_identity_tokens(sequence_name: Any, anticodon: Any) -> IdentityTokens:
    name_tokens = _identity_from_words(_normalized_words(sequence_name), filename=False)
    anticodon_words = _normalized_words(anticodon)
    explicit_anticodons = {
        token for token in anticodon_words if ANTICODON_TOKEN.fullmatch(token)
    }
    amino = set(name_tokens.amino_acids)
    anticodons = set(name_tokens.anticodons) | explicit_anticodons
    combined = set(name_tokens.combined)
    combined.update(left + right for left in amino for right in anticodons)
    filename_tokens = set(name_tokens.filename_tokens) | set(anticodon_words)
    return IdentityTokens(
        tuple(sorted(filename_tokens)), tuple(sorted(amino)),
        tuple(sorted(anticodons)), tuple(sorted(combined)),
    )


def _match(left: tuple[str, ...], right: tuple[str, ...]) -> bool | None:
    if not left or not right:
        return None
    return bool(set(left) & set(right))


def _joined(values: tuple[str, ...]) -> str:
    return ";".join(values)


def _joined_amino(values: tuple[str, ...]) -> str:
    return ";".join(value.title() for value in values)


def _joined_anticodons(values: tuple[str, ...]) -> str:
    return ";".join(value.upper() for value in values)


def _joined_combined(values: tuple[str, ...]) -> str:
    return ";".join(value[:3].title() + value[3:].upper() for value in values)


def _display_identity(tokens: IdentityTokens) -> str:
    if tokens.combined:
        return "/".join(value[:3].title() + "-" + value[3:].upper() for value in tokens.combined)
    parts = [*(value.title() for value in tokens.amino_acids), *(value.upper() for value in tokens.anticodons)]
    return "/".join(parts) or "unresolved identity"


def audit_sciex_input_identity(
    source_path: str | Path,
    *,
    sequence_name: Any = "",
    sequence: Any = "",
    anticodon: Any = "",
    organism_group: Any = "",
    species: Any = "",
    condition_name: Any = "",
    parser_sample_name: Any = "",
) -> SciexInputIdentityAuditResult:
    """Return a deterministic one-row identity audit without reading or changing the input."""
    source = Path(source_path)
    sciex = extract_sciex_identity_tokens(source)
    configured = extract_configured_identity_tokens(sequence_name, anticodon)
    eligible = bool(str(sequence_name or "").strip() or str(sequence or "").strip())

    amino_match = _match(sciex.amino_acids, configured.amino_acids)
    anticodon_match = _match(sciex.anticodons, configured.anticodons)
    combined_match = _match(sciex.combined, configured.combined)
    comparable_count = sum(value is not None for value in (amino_match, anticodon_match))
    conflict = any(value is False for value in (amino_match, anticodon_match, combined_match))

    if not eligible:
        status = "NOT_ELIGIBLE"
    elif conflict:
        status = "CONFLICT"
    elif combined_match is True or (amino_match is True and anticodon_match is True):
        status = "MATCH"
    elif amino_match is True or anticodon_match is True:
        status = "PARTIAL_MATCH"
    else:
        status = "INSUFFICIENT_INFORMATION"

    if comparable_count == 2:
        evidence = "HIGH"
    elif comparable_count == 1 or combined_match is not None:
        evidence = "MEDIUM"
    elif sciex.amino_acids or sciex.anticodons or configured.amino_acids or configured.anticodons:
        evidence = "LOW"
    else:
        evidence = "NONE"

    warning_code = WARNING_CODE if status == "CONFLICT" else ""
    warning_message = ""
    if warning_code:
        warning_message = (
            f"SCIEX source filename suggests tRNA-{_display_identity(sciex)}, but the "
            f"configured sequence name/anticodon suggests tRNA-{_display_identity(configured)}. "
            "SCIEX comparison results are shadow diagnostics only and must not be interpreted "
            "as a biological identity match."
        )
    values = {
        "Audit_Status": status,
        "Audit_Eligible": eligible,
        "SCIEX_Source_File": str(source),
        "SCIEX_Source_Basename": source.name,
        "SCIEX_File_Stem": source.stem,
        "SCIEX_Parent_Directory": source.parent.name,
        "Parser_Sample_Name": str(parser_sample_name or ""),
        "SCIEX_Filename_Tokens": _joined(sciex.filename_tokens),
        "Configured_Sequence_Name": str(sequence_name or ""),
        "Configured_Sequence": str(sequence or ""),
        "Configured_Anticodon": str(anticodon or ""),
        "Configured_Organism_Group": str(organism_group or ""),
        "Configured_Species": str(species or ""),
        "Configured_Condition_Name": str(condition_name or ""),
        "SCIEX_Amino_Acid_Tokens": _joined_amino(sciex.amino_acids),
        "SCIEX_Anticodon_Tokens": _joined_anticodons(sciex.anticodons),
        "SCIEX_Combined_Identity_Tokens": _joined_combined(sciex.combined),
        "Configured_Amino_Acid_Tokens": _joined_amino(configured.amino_acids),
        "Configured_Anticodon_Tokens": _joined_anticodons(configured.anticodons),
        "Configured_Combined_Identity_Tokens": _joined_combined(configured.combined),
        "Amino_Acid_Match": amino_match,
        "Anticodon_Match": anticodon_match,
        "Combined_Identity_Match": combined_match,
        "Identity_Conflict": status == "CONFLICT",
        "Identity_Evidence_Level": evidence,
        "Warning_Code": warning_code,
        "Warning_Message": warning_message,
        "Biological_Interpretation_Eligible": status in {"MATCH", "PARTIAL_MATCH"},
        "Shadow_Only": True,
        "Applied_To_Formal_Score": False,
        "Applied_To_Ranking": False,
        "Applied_To_Candidate_Filtering": False,
        "Molecular_Identity_Assigned": False,
        "Notes": (
            "Filename/config metadata consistency only; no sequence inference, modification "
            "assignment, or molecular identity assignment."
        ),
    }
    return SciexInputIdentityAuditResult(values)
