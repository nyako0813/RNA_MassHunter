"""Base-loss theoretical ion generation and matching for the unmodified
d/w/a/z + precursor MS2 pathway (Feature G, Part 2/2).

Consumes the per-position base-loss mass lookup table built once at startup
by ``rna_masshunter.base_loss_masses.build_base_loss_masses`` and expands
each already-generated backbone ion (d/w/a/z, plus a synthetic full-length
"precursor" ion built here) into one candidate per base position in that
ion, representing the neutral loss of that base's free nucleobase.
"""

from typing import Any

from rna_masshunter.masses import calculate_unmodified_rna_mass
from rna_masshunter.models import Fragment, MS2IonMatch, MS2SpectrumInfo, TheoreticalMS2Ion
from rna_masshunter.ms1_mapping import ppm_error, theoretical_mz_from_mass

BASE_LOSS_THEORETICAL_ION_COLUMNS = [
    "Ion_ID",
    "Parent_Fragment_ID",
    "Parent_Sequence",
    "Ion_Type",
    "Ion_Sequence",
    "Ion_Start",
    "Ion_End",
    "Ion_Length",
    "Charge",
    "Theoretical_Mass",
    "Theoretical_mz",
    "Base_Loss_Base",
    "Base_Loss_Position",
    "Comment",
]

BASE_LOSS_ION_MATCH_COLUMNS = [
    "Spectrum_ID",
    "Scan_Index",
    "RT",
    "Precursor_mz",
    "Precursor_Charge",
    "Observed_mz",
    "Observed_Intensity",
    "Best_Ion_ID",
    "Best_Ion_Type",
    "Best_Ion_Sequence",
    "Ion_Length",
    "Parent_Fragment_ID",
    "Parent_Fragment_Sequence",
    "Ion_Charge",
    "Theoretical_mz",
    "Mass_Error_Da",
    "Mass_Error_ppm",
    "Match_Status",
    "Confidence",
    "Alternative_Candidates",
    "Base_Loss_Base",
    "Base_Loss_Position",
    "Comment",
]


def _as_bool(value: Any, default: bool) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _positive_int(value: Any, default: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed > 0 else default


def _ion_charges(ms2_config: dict[str, Any]) -> list[int]:
    raw = ms2_config.get("charge_states") or [1]
    values = raw if isinstance(raw, list) else [raw]
    charges = []
    for value in values:
        try:
            charge = abs(int(value))
        except (TypeError, ValueError):
            continue
        if charge > 0 and charge not in charges:
            charges.append(charge)
    return charges or [1]


def generate_precursor_ions_for_base_loss(
    theoretical_fragments: list[Fragment],
    config: Any,
    base_masses: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
) -> list[TheoreticalMS2Ion]:
    """Build one synthetic full-length "precursor" ion per (fragment, charge).

    These rows exist only to feed ``generate_base_loss_ions``; unlike d/w/a/z
    ions they are never added to the existing ``MS2_Theoretical_Ions`` sheet
    or ``ions`` list used elsewhere, so existing backbone-ion output is
    untouched. Unlike the "default"-terminal-form sub-fragment ions, the
    precursor mass must use the fragment's actual terminal form (e.g. a CCA
    tail), since it represents the intact fragment rather than a cut piece.
    """
    ms2 = getattr(config, "ms2_annotation", {}) or {}
    if not _as_bool(ms2.get("include_base_loss"), False):
        return []
    if not _as_bool(ms2.get("base_loss_include_precursor"), True):
        return []
    polarity = str(getattr(config, "instrument", {}).get("polarity", "negative") or "negative").lower()
    charges = _ion_charges(ms2)
    ions: list[TheoreticalMS2Ion] = []
    for fragment in theoretical_fragments or []:
        sequence = (fragment.sequence or "").upper().replace("T", "U")
        if not sequence:
            continue
        mass = calculate_unmodified_rna_mass(
            sequence, base_masses, warnings=warnings, terminal_form=fragment.terminal_form,
        )
        if mass is None:
            continue
        for charge in charges:
            ions.append(TheoreticalMS2Ion(
                ion_id=f"PRECION_{len(ions) + 1:06d}",
                parent_fragment_id=fragment.fragment_id,
                parent_sequence=sequence,
                ion_type="precursor",
                ion_sequence=sequence,
                ion_start=1,
                ion_end=len(sequence),
                charge=charge,
                theoretical_mass=float(mass),
                theoretical_mz=theoretical_mz_from_mass(float(mass), charge, polarity),
                comment="full-length precursor ion; base-loss generation input only",
            ))
    return ions


def generate_base_loss_ions(
    ions: list[TheoreticalMS2Ion],
    config: Any,
    base_loss_masses: dict[str, float],
) -> list[TheoreticalMS2Ion]:
    """Expand each backbone/precursor ion into one row per base-loss position."""
    ms2 = getattr(config, "ms2_annotation", {}) or {}
    if not _as_bool(ms2.get("include_base_loss"), False):
        return []
    min_length = _positive_int(ms2.get("base_loss_min_ion_length"), 2)
    max_rows = _positive_int(ms2.get("base_loss_max_rows"), 50000)
    polarity = str(getattr(config, "instrument", {}).get("polarity", "negative") or "negative").lower()

    rows: list[TheoreticalMS2Ion] = []
    for ion in ions:
        sequence = ion.ion_sequence or ""
        if len(sequence) < min_length:
            continue
        for offset, base in enumerate(sequence):
            loss_mass = base_loss_masses.get(base)
            if loss_mass is None:
                continue
            new_mass = ion.theoretical_mass - loss_mass
            position = ion.ion_start + offset
            rows.append(TheoreticalMS2Ion(
                ion_id=f"BL_{len(rows) + 1:08d}",
                parent_fragment_id=ion.parent_fragment_id,
                parent_sequence=ion.parent_sequence,
                ion_type=ion.ion_type,
                ion_sequence=sequence,
                ion_start=ion.ion_start,
                ion_end=ion.ion_end,
                charge=ion.charge,
                theoretical_mass=new_mass,
                theoretical_mz=theoretical_mz_from_mass(new_mass, ion.charge, polarity),
                neutral_loss=base,
                base_loss_position=position,
                base_loss_base=base,
                comment=f"base loss: {base} at parent position {position}",
            ))
            if len(rows) >= max_rows:
                return rows
    return rows


def match_base_loss_ions(
    spectra: list[MS2SpectrumInfo],
    base_loss_ions: list[TheoreticalMS2Ion],
    config: Any,
) -> list[MS2IonMatch]:
    ms2 = getattr(config, "ms2_annotation", {}) or {}
    tolerance_ppm = float(ms2.get("mz_tolerance_ppm", 20) or 20)
    matches: list[MS2IonMatch] = []
    if not base_loss_ions:
        return matches

    for spectrum in spectra:
        for observed_mz, observed_intensity in spectrum.peaks:
            ranked = sorted(
                (
                    (ion, observed_mz - ion.theoretical_mz, ppm_error(observed_mz, ion.theoretical_mz))
                    for ion in base_loss_ions
                ),
                key=lambda item: abs(item[2]),
            )
            within = [item for item in ranked if abs(item[2]) <= tolerance_ppm]
            if not within:
                continue
            ion, error_da, error_ppm = within[0]
            alternatives = within[1:]
            status = "multiple_candidates" if alternatives else "matched"
            relative = observed_intensity / spectrum.base_peak_intensity if spectrum.base_peak_intensity else 1.0
            if status == "multiple_candidates":
                confidence = "Low"
            elif abs(error_ppm) <= tolerance_ppm / 3 and relative >= 0.1:
                confidence = "High"
            elif abs(error_ppm) <= tolerance_ppm:
                confidence = "Medium"
            else:
                confidence = "Low"
            matches.append(MS2IonMatch(
                spectrum_id=spectrum.spectrum_id,
                scan_index=spectrum.scan_index,
                rt=spectrum.rt,
                precursor_mz=spectrum.precursor_mz,
                precursor_charge=spectrum.precursor_charge,
                observed_mz=observed_mz,
                observed_intensity=observed_intensity,
                best_ion_id=ion.ion_id,
                best_ion_type=ion.ion_type,
                best_ion_sequence=ion.ion_sequence,
                parent_fragment_id=ion.parent_fragment_id,
                parent_fragment_sequence=ion.parent_sequence,
                ion_charge=ion.charge,
                theoretical_mz=ion.theoretical_mz,
                mass_error_da=error_da,
                mass_error_ppm=error_ppm,
                match_status=status,
                confidence=confidence,
                alternative_candidates="; ".join(candidate[0].ion_id for candidate in alternatives[:5]),
                base_loss_position=ion.base_loss_position,
                base_loss_base=ion.base_loss_base,
                comment=f"base loss: {ion.base_loss_base} at parent position {ion.base_loss_position}",
            ))
    return matches


def base_loss_ion_row(ion: TheoreticalMS2Ion, config: Any) -> dict[str, Any]:
    return {
        "Ion_ID": ion.ion_id,
        "Parent_Fragment_ID": ion.parent_fragment_id,
        "Parent_Sequence": ion.parent_sequence,
        "Ion_Type": ion.ion_type,
        "Ion_Sequence": ion.ion_sequence,
        "Ion_Start": ion.ion_start,
        "Ion_End": ion.ion_end,
        "Ion_Length": len(ion.ion_sequence or ""),
        "Charge": ion.charge,
        "Theoretical_Mass": ion.theoretical_mass,
        "Theoretical_mz": ion.theoretical_mz,
        "Base_Loss_Base": ion.base_loss_base,
        "Base_Loss_Position": ion.base_loss_position,
        "Comment": ion.comment,
    }


def base_loss_match_row(match: MS2IonMatch, config: Any) -> dict[str, Any]:
    return {
        "Spectrum_ID": match.spectrum_id,
        "Scan_Index": match.scan_index,
        "RT": match.rt,
        "Precursor_mz": match.precursor_mz,
        "Precursor_Charge": match.precursor_charge,
        "Observed_mz": match.observed_mz,
        "Observed_Intensity": match.observed_intensity,
        "Best_Ion_ID": match.best_ion_id,
        "Best_Ion_Type": match.best_ion_type,
        "Best_Ion_Sequence": match.best_ion_sequence,
        "Ion_Length": len(match.best_ion_sequence or ""),
        "Parent_Fragment_ID": match.parent_fragment_id,
        "Parent_Fragment_Sequence": match.parent_fragment_sequence,
        "Ion_Charge": match.ion_charge,
        "Theoretical_mz": match.theoretical_mz,
        "Mass_Error_Da": match.mass_error_da,
        "Mass_Error_ppm": match.mass_error_ppm,
        "Match_Status": match.match_status,
        "Confidence": match.confidence,
        "Alternative_Candidates": match.alternative_candidates,
        "Base_Loss_Base": match.base_loss_base,
        "Base_Loss_Position": match.base_loss_position,
        "Comment": match.comment,
    }
