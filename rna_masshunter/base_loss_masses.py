"""Base-loss mass lookup table (Feature G: base-loss ions).

Computes, once at startup, the neutral base-loss mass for every entry in
the modification catalog plus the four canonical (unmodified) RNA bases.
"Base loss" is the mass of the free nucleobase (BH form) released when the
N-glycosidic bond breaks during CID/HCD/EAD fragmentation. Downstream ion
generation code should look up this precomputed dict rather than
recomputing per-ion, since the same base-loss mass is reused across many
ions/positions/spectra.

Decomposition rule for compound base+sugar modifications
----------------------------------------------------------
A modification's own curated ``mass_shift_from_unmodified`` mixes the
base-side and sugar-side contribution and cannot be used directly for
base-loss when the modification also carries a 2'-O-sugar component
(e.g. m1Am = m1A base modification + Am sugar modification). This catalog
follows a consistent naming convention: a compound entry's id/symbol is
the base-only counterpart's id with a trailing lowercase "m" appended
(e.g. "m1A" + "m" -> "m1Am"; "s2U" + "m" -> "s2Um"). This module exploits
that convention:

- id ends in "m" and the id with "m" stripped equals a canonical base
  letter (A/U/G/C): pure 2'-O-methyl sugar modification, base itself is
  unmodified -> base-loss shift = 0.
- id ends in "m" and the stripped id matches another catalog entry:
  compound modification -> reuse that base-only entry's own curated
  ``mass_shift_from_unmodified`` as the base-side contribution.
- id ends in "m" and the stripped id matches neither: the base-side
  component cannot be resolved from the catalog. The modification is
  excluded from base-loss generation and flagged via ``warnings`` rather
  than guessed.
- ids in ``PURE_SUGAR_SPECIAL_IDS`` (2'-O-ribosyl/phosphate additions that
  do not follow the trailing-"m" convention): treated as pure sugar
  modifications, same as the first case above.
- any other id: treated as a pure base-side modification; the full
  curated shift applies to the lost base.

This rule was validated against all 27 catalog entries with a 2'-O-related
name: 26 resolved to totals matching the PDF-curated compound shift to
within 0.005 Da when re-added to the sugar-side (+14.0157 Da) shift, and
the 1 exception (m4_4Cm) was covered by adding a derived base-only
counterpart entry (m4_4C) to the catalog.
"""

from typing import Any

from rna_masshunter.models import Modification
from rna_masshunter.warnings_manager import add_warning

CANONICAL_BASES = {"A", "U", "G", "C"}

# Ids that are pure 2'-O-sugar modifications (base itself unmodified) but
# whose symbol does not follow the "<base_component_id>m" suffix
# convention, e.g. 2'-O-ribosylation/phosphate additions.
PURE_SUGAR_SPECIAL_IDS = {"Arp", "Grp"}


def build_base_loss_masses(
    modifications: list[Modification],
    base_masses: dict[str, Any],
    warnings: list[dict[str, Any]] | None = None,
) -> dict[str, float]:
    """Return a dict mapping a lookup key to a neutral base-loss mass.

    Keys are the four canonical base letters ("A", "U", "G", "C") for
    unmodified positions, plus every resolvable modification id for
    modified positions. Values are the monoisotopic mass of the free
    nucleobase released at that position upon N-glycosidic bond cleavage.
    """
    nucleobase_loss = base_masses.get("nucleobase_loss_masses") or {}
    missing_canonical = CANONICAL_BASES - set(nucleobase_loss)
    if missing_canonical and warnings is not None:
        add_warning(
            warnings, "ERROR", "base_loss_masses",
            "base_masses.yaml is missing nucleobase_loss_masses for some canonical bases; base-loss ions cannot be generated for them.",
            sorted(missing_canonical),
        )

    table: dict[str, float] = {
        base: float(mass) for base, mass in nucleobase_loss.items() if base in CANONICAL_BASES
    }

    by_id = {mod.id: mod for mod in modifications if mod.id}

    for mod in modifications:
        mod_id = mod.id
        if not mod_id:
            continue
        target = mod.target_bases[0] if mod.target_bases else None
        canonical_loss = nucleobase_loss.get(target) if target else None
        if canonical_loss is None:
            if warnings is not None:
                add_warning(
                    warnings, "WARNING", "base_loss_masses",
                    "Modification has no resolvable canonical target base for base-loss; excluded from base-loss ion generation.",
                    mod_id,
                )
            continue

        if mod_id in PURE_SUGAR_SPECIAL_IDS:
            table[mod_id] = float(canonical_loss)
            continue

        if mod_id.endswith("m") and len(mod_id) > 1:
            stem = mod_id[:-1]
            if stem in CANONICAL_BASES:
                table[mod_id] = float(canonical_loss)
                continue
            if stem in by_id:
                base_component = by_id[stem]
                table[mod_id] = float(canonical_loss) + float(base_component.mass_shift_from_unmodified)
                continue
            if warnings is not None:
                add_warning(
                    warnings, "WARNING", "base_loss_masses",
                    "Compound base+sugar modification could not be decomposed (base-only counterpart not found in catalog); excluded from base-loss ion generation.",
                    {"modification_id": mod_id, "expected_base_component_id": stem},
                )
            continue

        # Pure base-side modification: the full curated shift applies to the lost base.
        table[mod_id] = float(canonical_loss) + float(mod.mass_shift_from_unmodified)

    return table
