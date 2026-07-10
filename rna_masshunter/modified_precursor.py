"""Modified precursor candidate generation for MVP-5.2 MS2 annotation."""

from typing import Any

from rna_masshunter.models import Fragment, Modification, MS2SpectrumInfo
from rna_masshunter.ms1_mapping import ppm_error, theoretical_mz_from_mass


def _as_bool(value: Any, default: bool) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _target_bases(modification: Modification) -> list[str]:
    values = modification.target_bases or modification.raw.get("target_bases") or modification.raw.get("target_base") or []
    if isinstance(values, str):
        values = [values]
    return [str(value).upper().replace("T", "U") for value in values if value is not None]


def _compatible(sequence: str, modification: Modification, required: bool) -> tuple[bool, str]:
    bases = _target_bases(modification)
    if not required or not bases or any(base in {"", "ANY", "UNKNOWN", "N", "*"} for base in bases):
        return True, "target base unknown/any; compatibility accepted" if required else ""
    normalized = sequence.upper().replace("T", "U")
    return any(base in normalized for base in bases), ""


def modification_name(modification: Modification) -> str:
    return str(modification.raw.get("name") or modification.raw.get("full_name") or modification.symbol or modification.id)


def find_modified_parent_candidates(
    spectrum: MS2SpectrumInfo,
    fragments: list[Fragment],
    modifications: list[Modification],
    config: Any,
) -> list[dict[str, Any]]:
    """Return one-known-modification precursor candidates compatible with a spectrum."""
    if spectrum.precursor_mz is None:
        return []
    ms2 = getattr(config, "ms2_annotation", {}) or {}
    if not _as_bool(ms2.get("include_modified_precursor_candidates"), True):
        return []
    if int(ms2.get("modified_precursor_max_mods_per_fragment", 1) or 1) < 1:
        return []
    tolerance_ppm = float(ms2.get("precursor_match_tolerance_ppm", 20) or 20)
    zero_tolerance = float(ms2.get("modified_precursor_mass_shift_tolerance_da", 1e-6) or 1e-6)
    include_isobaric = _as_bool(ms2.get("modified_precursor_include_isobaric"), False)
    require_base = _as_bool(ms2.get("modified_precursor_require_base_compatibility"), True)
    polarity = str(getattr(config, "instrument", {}).get("polarity", "negative") or "negative").lower()
    charges = [abs(int(spectrum.precursor_charge))] if spectrum.precursor_charge else _charges(ms2)
    candidates: list[dict[str, Any]] = []
    for fragment in fragments or []:
        for modification in modifications or []:
            policy = getattr(modification, "candidate_policy", None) or (getattr(modification, "raw", {}) or {}).get("candidate_policy", {})
            if not _as_bool(policy.get("include_by_mass_search"), True):
                continue
            shift = modification.mass_shift_from_unmodified
            if shift != shift or (not include_isobaric and abs(shift) <= zero_tolerance):
                continue
            compatible, note = _compatible(fragment.sequence, modification, require_base)
            if not compatible:
                continue
            modified_mass = float(fragment.unmodified_mass) + float(shift)
            if modified_mass <= 0:
                continue
            bases = _target_bases(modification)
            for charge in charges:
                mz = theoretical_mz_from_mass(modified_mass, charge, polarity)
                error_ppm = ppm_error(float(spectrum.precursor_mz), mz)
                if abs(error_ppm) <= tolerance_ppm:
                    candidates.append({
                        "fragment": fragment, "charge": charge, "theoretical_mz": mz,
                        "error_da": float(spectrum.precursor_mz) - mz, "error_ppm": error_ppm,
                        "candidate_type": "modified", "modification_id": modification.id,
                        "modification_name": modification_name(modification),
                        "modification_target_base": ",".join(bases),
                        "modification_mass_shift": float(shift),
                        "unmodified_mass": float(fragment.unmodified_mass), "modified_mass": modified_mass,
                        "comment": note,
                    })
    return candidates


def _charges(ms2_config: dict[str, Any]) -> list[int]:
    values = ms2_config.get("precursor_charge_states") or ms2_config.get("ion_charge_states") or [1]
    if isinstance(values, (int, float, str)):
        values = [values]
    charges = []
    for value in values:
        try:
            charge = abs(int(value))
        except (TypeError, ValueError):
            continue
        if charge and charge not in charges:
            charges.append(charge)
    return charges or [1]
