"""Position-independent P1 + SAP chemical-state shadow audit.

This module deliberately does not propagate source positions or original bonds.  It
uses canonical product chemistry, groups profile points into features, and keeps
all formal-isolation flags false.
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, is_dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import subprocess
import time
import tracemalloc
from typing import Any

from rna_masshunter.elemental_composition import ElementalComposition
from rna_masshunter.enzymes import normalize_enzyme_name
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.modification_constraints import load_transformations
from rna_masshunter.resource_utils import get_maximum_rss_mib
from rna_masshunter.p1_sap_feature_quality import build_p1_sap_feature_quality
from rna_masshunter.p1_sap_dinucleotide_interpretation import build_p1_sap_dinucleotide_audit
from rna_masshunter.mzml_diagnostics import _rt_minutes
from rna_masshunter.mzml_reader import iter_spectra

FALSE_FLAGS = {
    "Applied_To_Formal_Result": False,
    "Formal_Change_Ready": False,
    "Formal_Result_Changed": False,
}
LOCALIZATION = {
    "Chemical_State_Supported": False,
    "Sequence_Position_Localized": False,
    "Original_Bond_Localized": False,
    "Position_Resolvable": False,
    "Bond_Resolvable": False,
    "Localization_Status": "POSITION_NOT_RETAINED_BY_PREPARATION",
}
NUCLEOSIDE_FORMULAS = {
    "A": ElementalComposition({"C": 10, "H": 13, "N": 5, "O": 4}),
    "C": ElementalComposition({"C": 9, "H": 13, "N": 3, "O": 5}),
    "G": ElementalComposition({"C": 10, "H": 13, "N": 5, "O": 5}),
    "U": ElementalComposition({"C": 9, "H": 12, "N": 2, "O": 6}),
}
PHOSPHATE_MONOESTER = ElementalComposition.delta({"H": 1, "O": 3, "P": 1})
PHOSPHODIESTER = ElementalComposition.delta({"H": -1, "O": 2, "P": 1})
O_TO_S = ElementalComposition.delta({"O": -1, "S": 1})
ISOTOPE_SPACING = 1.00335483507
MODEL_NOT_DEFINED = "MODEL_NOT_DEFINED"

CHEMICAL_STATE_COLUMNS = [
    "Chemical_State_ID", "Chemical_Family", "Product_Type", "Base_or_Oligomer_Composition",
    "Nucleoside_Modification_State", "Phosphate_State", "Phosphorothioate_State",
    "Sulfur_Count", "Oxidation_State", "Terminal_State", "Elemental_Composition",
    "Neutral_Mass", "Charge", "Theoretical_mz", "Observable", "Not_Observable_Reason",
    "Raw_Profile_Point_Match_Count", "Grouped_Feature_Count", "Independent_Feature_Count",
    "Candidate_Specific_Feature_Count", "Best_Feature_ID", "Best_Observed_mz",
    "Best_Mass_Error_ppm", "Best_RT_Apex", "Best_Integrated_Intensity",
    "Possible_Source_Position_Count", "Possible_Source_Positions", "Possible_Source_Bond_Count",
    "Possible_Source_Bonds", "Possible_Source_Hypothesis_IDs", "Chemical_State_Supported",
    "PT_Chemical_State_Interpretation", "Sequence_Position_Localized", "Original_Bond_Localized",
    "Position_Resolvable", "Bond_Resolvable", "Localization_Status", "Model_Status",
    "Search_Enabled", "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]
FEATURE_COLUMNS = [
    "Feature_ID", "Physical_Feature_ID", "Chemical_State_ID", "Chemical_Family", "Product_Type",
    "Elemental_Composition", "Charge", "Theoretical_mz", "RT_Start", "RT_End", "RT_Apex",
    "RT_Span", "Apex_mz", "mz_Centroid", "mz_SD", "Apex_Intensity", "Integrated_Intensity",
    "Profile_Point_Count", "Spectrum_Count", "Mass_Error_ppm_at_Apex",
    "Mass_Error_ppm_at_Centroid", "Expected_Isotope_Spacing", "Observed_Isotope_Spacing",
    "Estimated_Charge", "Monoisotopic_Candidate", "Minus1_Isotope_Candidate",
    "Plus1_Isotope_Candidate", "Envelope_Assessed", "Isotope_Status",
    "Feature_Continuity_Status", "Feature_Eligible_For_Support", "Feature_Exclusion_Reason",
    "Normal_Phosphate_Competition_Count", "PT_Competition_Count",
    "Thiophosphate_Competition_Count", "Oxidized_PT_Competition_Count",
    "Sulfur_Non_PT_Competition_Count", "Isotope_Competition_Count", "Charge_Competition_Count",
    "Composition_Competition_Count", "Competition_Count", "Candidate_Specific",
    "Final_Interpretation", "Chemical_State_Supported", "Sequence_Position_Localized",
    "Original_Bond_Localized", "Localization_Status", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]
PT_FAMILY_COLUMNS = [
    "Family_ID", "Chemical_Base_State", "Normal_State_ID", "Dephosphorylated_State_ID",
    "PT_State_ID", "Thiophosphate_State_ID", "Oxidized_PT_State_ID",
    "Normal_Composition", "PT_Composition", "Normal_Theoretical_mz", "PT_Theoretical_mz",
    "Expected_Delta_Da", "Expected_Delta_mz", "Normal_Feature_Count",
    "Dephosphorylated_Feature_Count", "PT_Feature_Count", "Thiophosphate_Feature_Count",
    "Oxidized_PT_Feature_Count", "P1_Resistant_PT_Oligomer_Feature_Count",
    "Ambiguous_Sulfur_Feature_Count", "Family_Interpretation", "Sequence_Position_Localized",
    "Original_Bond_Localized", "Localization_Status", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]
TERMINAL_COLUMNS = [
    "Chemical_Base_State", "SAP_Substrate_State", "Dephosphorylated_Feature_Present",
    "Residual_Phosphate_Feature_Present", "PT_Like_Feature_Present", "SAP_Interpretation",
    "SAP_Removal_Expected", "SAP_Removal_Confirmed", "SAP_Removal_Unknown", "SAP_Model_Reason",
    "Observable", "Feature_Present", "Feature_ID", "Observed_mz", "Mass_Error_ppm", "RT_Apex",
    "Integrated_Intensity", "Candidate_Specific", "Sequence_Position_Localized",
    "Original_Bond_Localized", "Localization_Status", "Applied_To_Formal_Result",
    "Formal_Change_Ready", "Formal_Result_Changed",
]
CROSS_ENZYME_COLUMNS = [
    "Chemical_State_ID", "Elemental_Composition", "T1_Support_Status", "P1_SAP_Support_Status",
    "Cross_Enzyme_Chemical_Status", "Position_Confirmed", "Reason",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]
COMPETITION_COLUMNS = [
    "Physical_Feature_ID", "Feature_ID", "Chemical_State_ID", "Competing_State_ID",
    "Candidate_Chemical_Family", "Competing_Chemical_Family", "Same_Elemental_Composition",
    "Same_Charge", "Competition_Type", "Applied_To_Formal_Result", "Formal_Change_Ready",
    "Formal_Result_Changed",
]
MS2_PROVENANCE_COLUMNS = [
    "Chemical_State_ID", "Precursor_Compatible_MS2_Count", "MS2_Spectrum_ID", "Spectrum_RT",
    "Isolation_Window", "Collision_Energy", "Precursor_mz", "Precursor_Charge",
    "Precursor_Mass_Error_ppm", "MS2_Model_Applicable", "MS2_Model_Not_Applicable_Reason",
    "Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed",
]
SUMMARY_COLUMNS = [
    "Configured_Enzyme", "Rule_ID", "Cleavage_Specificity", "Cleavage_Side",
    "Expected_Product_Type", "Expected_5prime_Terminal_State", "Expected_3prime_Terminal_State",
    "Expected_Phosphate_Product", "Missed_Cleavage_Model", "P1_Cleaves_Normal_Phosphodiester",
    "P1_PT_Cleavage_Behavior", "Candidate_Count", "Searchable_Candidate_Count",
    "Raw_Profile_Point_Match_Count", "Grouped_Feature_Count", "Independent_Feature_Count",
    "Feature_Apex_Match_Count", "Dephosphorylated_Feature_Count",
    "Residual_Normal_Phosphate_Feature_Count", "PT_Feature_Count",
    "Thiophosphate_Like_Feature_Count", "Oxidized_PT_Like_Feature_Count",
    "P1_Resistant_PT_Oligomer_Feature_Count", "Ambiguous_Sulfur_Feature_Count",
    "PT_Chemical_State_Final_Interpretation", "Position_Localization_Disabled",
    "Bond_Localization_Disabled", "Audit_Runtime", "Maximum_RSS_MiB", "Tracemalloc_Peak_MiB",
    "Feature_Detail_Row_Count", "Applied_To_Formal_Result", "Formal_Change_Ready",
    "Formal_Result_Changed",
]

@dataclass
class P1SAPAuditResult:
    sheets: dict[str, list[dict[str, Any]]]
    metrics: dict[str, Any]
    summary_payload: dict[str, Any]


def _safe_id(value: Any) -> str:
    return "".join(ch if str(ch).isalnum() else "_" for ch in str(value)).strip("_")


def _positions(sequence: str, base: str) -> tuple[int, ...]:
    return tuple(i for i, value in enumerate(sequence, 1) if value == base)


def _bonds(sequence: str, base: str) -> tuple[str, ...]:
    return tuple(f"{i}_{i+1}" for i, value in enumerate(sequence[:-1], 1) if value == base)


def _composition_add(base: ElementalComposition | None, delta: ElementalComposition) -> ElementalComposition | None:
    if base is None:
        return None
    try:
        return base + delta
    except ValueError:
        return None


def _candidate_row(*, state_id: str, family: str, product_type: str, base_key: str,
                   modification: str, phosphate: str, pt_state: str, composition: ElementalComposition | None,
                   neutral_mass: float | None, charge: int, oxidation: str, terminal: str,
                   positions: tuple[int, ...] = (), bonds: tuple[str, ...] = (), source_hypotheses: str = "",
                   search_enabled: bool = True, model_status: str = "defined") -> dict[str, Any]:
    theoretical = mz_from_neutral_mass(neutral_mass, charge, "positive") if neutral_mass is not None else None
    observable = theoretical is not None and 0 < theoretical < math.inf
    sulfur = composition.to_dict().get("S", 0) if composition is not None else "unknown"
    return {
        "Chemical_State_ID": state_id, "Chemical_Family": family, "Product_Type": product_type,
        "Base_or_Oligomer_Composition": base_key, "Nucleoside_Modification_State": modification,
        "Phosphate_State": phosphate, "Phosphorothioate_State": pt_state,
        "Sulfur_Count": sulfur, "Oxidation_State": oxidation, "Terminal_State": terminal,
        "Elemental_Composition": composition.canonical_string() if composition is not None else MODEL_NOT_DEFINED,
        "Neutral_Mass": neutral_mass, "Charge": charge, "Theoretical_mz": theoretical,
        "Observable": observable, "Not_Observable_Reason": "" if observable else MODEL_NOT_DEFINED,
        "Possible_Source_Position_Count": len(positions), "Possible_Source_Positions": ";".join(map(str, positions)),
        "Possible_Source_Bond_Count": len(bonds), "Possible_Source_Bonds": ";".join(bonds),
        "Possible_Source_Hypothesis_IDs": source_hypotheses, "Model_Status": model_status,
        "Search_Enabled": bool(search_enabled and theoretical is not None), **LOCALIZATION, **FALSE_FLAGS,
    }


def generate_chemical_state_candidates(sequence: str, modifications: list[Any], project_root: Path,
                                       charges: tuple[int, ...] = (1,)) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Generate canonical chemistry; source positions/bonds never enter state identity."""
    transforms = {item.id: item for item in load_transformations(project_root / "data/modification_transforms_v2.yaml")}
    chemical_bases: list[dict[str, Any]] = []
    for base, composition in NUCLEOSIDE_FORMULAS.items():
        chemical_bases.append({"key": base, "base": base, "mod": "unmodified", "composition": composition,
            "mass": composition.exact_mass, "hypotheses": ""})
    seen = {(x["base"], x["mod"]) for x in chemical_bases}
    for item in modifications or ():
        for base in tuple(getattr(item, "target_bases", ()) or ()):
            key = (base, str(item.id))
            if base not in NUCLEOSIDE_FORMULAS or key in seen:
                continue
            transform = transforms.get(str(item.id)); composition = None
            if transform is not None:
                composition = _composition_add(NUCLEOSIDE_FORMULAS[base], transform.composition_delta)
            mass = NUCLEOSIDE_FORMULAS[base].exact_mass + float(item.mass_shift_from_unmodified)
            chemical_bases.append({"key": f"{base}:{item.id}", "base": base, "mod": str(item.id),
                "composition": composition, "mass": composition.exact_mass if composition is not None else mass,
                "hypotheses": "Mac_m22G10" if str(item.id) == "m22G" else ""})
            seen.add(key)
    # Canonical m22G is transform-defined even when the broad modification dictionary uses an alias ID.
    if ("G", "m22G") not in seen:
        m22_composition = NUCLEOSIDE_FORMULAS["G"] + transforms["m22G"].composition_delta
        chemical_bases.append({"key":"G:m22G","base":"G","mod":"m22G",
            "composition":m22_composition,"mass":m22_composition.exact_mass,"hypotheses":"Mac_m22G10"})
        seen.add(("G", "m22G"))
    # Explicit schema-defined U side-chain composite states.
    for mod_ids, label, hypotheses in (
        (("s2U", "cnm5U", "side_chain_thioamide"), "U_side_chain_thioamide", "U37_side_chain_thioamide"),
        (("s2U", "cnm5U", "side_chain_thioamide_oxo1"), "U_side_chain_thioamide_oxo1", "U37_side_chain_thioamide_oxo1"),
    ):
        delta = ElementalComposition.delta()
        for mod_id in mod_ids:
            delta = delta + transforms[mod_id].composition_delta
        composition = NUCLEOSIDE_FORMULAS["U"] + delta
        chemical_bases.append({"key": f"U:{label}", "base": "U", "mod": label,
            "composition": composition, "mass": composition.exact_mass, "hypotheses": hypotheses})

    candidates: list[dict[str, Any]] = []
    family_specs: list[dict[str, Any]] = []
    for chemical in chemical_bases:
        base = chemical["base"]; key = chemical["key"]; slug = _safe_id(key)
        positions = _positions(sequence, base); bonds = _bonds(sequence, base)
        comp = chemical["composition"]; mass = float(chemical["mass"])
        sulfur_non_pt = comp is not None and comp.to_dict().get("S", 0) > 0
        for charge in charges:
            suffix = f"z{charge}"
            dephos_family = "SULFUR_CONTAINING_NON_PT_ALTERNATIVE" if sulfur_non_pt else "DEPHOSPHORYLATED"
            dephos = _candidate_row(state_id=f"P1SAP_{slug}_dephos_{suffix}", family=dephos_family,
                product_type="monomer", base_key=key, modification=chemical["mod"], phosphate="none",
                pt_state="none", composition=comp, neutral_mass=mass, charge=charge, oxidation="none",
                terminal="dephosphorylated", positions=positions, bonds=bonds,
                source_hypotheses=chemical["hypotheses"])
            normal_comp = _composition_add(comp, PHOSPHATE_MONOESTER)
            normal_mass = mass + PHOSPHATE_MONOESTER.exact_mass
            normal = _candidate_row(state_id=f"P1SAP_{slug}_normal_phosphate_reference_{suffix}",
                family="NORMAL_PHOSPHATE", product_type="monomer", base_key=key,
                modification=chemical["mod"], phosphate="normal_phosphate_reference", pt_state="none",
                composition=normal_comp, neutral_mass=normal_mass, charge=charge, oxidation="none",
                terminal="normal_phosphate", positions=positions, bonds=bonds,
                source_hypotheses=chemical["hypotheses"], search_enabled=False)
            residual = _candidate_row(state_id=f"P1SAP_{slug}_residual_phosphate_{suffix}",
                family="RESIDUAL_NORMAL_PHOSPHATE", product_type="monomer", base_key=key,
                modification=chemical["mod"], phosphate="residual_normal_phosphate", pt_state="none",
                composition=normal_comp, neutral_mass=normal_mass, charge=charge, oxidation="none",
                terminal="residual_phosphate", positions=positions, bonds=bonds,
                source_hypotheses=chemical["hypotheses"])
            pt_comp = _composition_add(normal_comp, O_TO_S)
            pt_mass = normal_mass + O_TO_S.exact_mass
            pt = _candidate_row(state_id=f"P1SAP_{slug}_phosphorothioate_{suffix}",
                family="PHOSPHOROTHIOATE", product_type="monomer", base_key=key,
                modification=chemical["mod"], phosphate="phosphate_monoester", pt_state="O_to_S",
                composition=pt_comp, neutral_mass=pt_mass, charge=charge, oxidation="unoxidized",
                terminal="phosphorothioate_monoester", positions=positions, bonds=bonds,
                source_hypotheses=chemical["hypotheses"])
            thiophosphate = _candidate_row(state_id=f"P1SAP_{slug}_thiophosphate_like_{suffix}",
                family="THIOPHOSPHATE_LIKE", product_type="monomer", base_key=key,
                modification=chemical["mod"], phosphate="unknown", pt_state="thiophosphate_like",
                composition=None, neutral_mass=None, charge=charge, oxidation="unknown",
                terminal="thiophosphate_like", positions=positions, bonds=bonds,
                source_hypotheses=chemical["hypotheses"], search_enabled=False, model_status=MODEL_NOT_DEFINED)
            oxidized = _candidate_row(state_id=f"P1SAP_{slug}_oxidized_pt_{suffix}",
                family="OXIDIZED_PT_DERIVATIVE", product_type="monomer", base_key=key,
                modification=chemical["mod"], phosphate="unknown", pt_state="oxidized_PT",
                composition=None, neutral_mass=None, charge=charge, oxidation="unknown",
                terminal="oxidized_PT", positions=positions, bonds=bonds,
                source_hypotheses=chemical["hypotheses"], search_enabled=False, model_status=MODEL_NOT_DEFINED)
            candidates.extend((dephos, normal, residual, pt, thiophosphate, oxidized))
            family_specs.append({"Family_ID": f"P1SAP_FAMILY_{slug}_{suffix}", "Chemical_Base_State": key,
                "Normal_State_ID": normal["Chemical_State_ID"], "Dephosphorylated_State_ID": dephos["Chemical_State_ID"],
                "Residual_State_ID": residual["Chemical_State_ID"], "PT_State_ID": pt["Chemical_State_ID"],
                "Thiophosphate_State_ID": thiophosphate["Chemical_State_ID"], "Oxidized_PT_State_ID": oxidized["Chemical_State_ID"],
                "Normal_Composition": normal["Elemental_Composition"], "PT_Composition": pt["Elemental_Composition"],
                "Normal_Theoretical_mz": normal["Theoretical_mz"], "PT_Theoretical_mz": pt["Theoretical_mz"],
                "Expected_Delta_Da": O_TO_S.exact_mass, "Expected_Delta_mz": O_TO_S.exact_mass / charge})

    # Bounded, composition-deduplicated hypothetical P1-resistant PT dimers.
    dimer_sources: dict[str, dict[str, set[Any]]] = defaultdict(lambda: {"positions": set(), "bonds": set()})
    for i in range(len(sequence) - 1):
        key = "".join(sorted(sequence[i:i + 2]))
        dimer_sources[key]["positions"].update((i + 1, i + 2)); dimer_sources[key]["bonds"].add(f"{i+1}_{i+2}")
    for key, sources in sorted(dimer_sources.items()):
        comp = NUCLEOSIDE_FORMULAS[key[0]] + NUCLEOSIDE_FORMULAS[key[1]] + PHOSPHODIESTER + O_TO_S
        for charge in charges:
            candidates.append(_candidate_row(state_id=f"P1SAP_{key}_P1_resistant_internal_PT_z{charge}",
                family="P1_RESISTANT_PT_OLIGOMER", product_type="dimer", base_key=key,
                modification="unmodified", phosphate="internal_phosphodiester", pt_state="internal_O_to_S",
                composition=comp, neutral_mass=comp.exact_mass, charge=charge, oxidation="unoxidized",
                terminal="dephosphorylated_ends", positions=tuple(sorted(sources["positions"])),
                bonds=tuple(sorted(sources["bonds"])), search_enabled=True,
                model_status="hypothetical_P1_resistance;P1_PT_cleavage_behavior_unknown"))
    # Canonical state ID is the dedup key; source coordinates are metadata only.
    unique = {row["Chemical_State_ID"]: row for row in candidates}
    return list(unique.values()), family_specs


def _raw_peak(peak: Any, index: int) -> dict[str, Any]:
    raw = asdict(peak) if is_dataclass(peak) else dict(peak)
    return {"index": index, "mz": float(raw.get("mz") or 0), "intensity": float(raw.get("intensity") or 0),
        "rt": float(raw.get("rt")) if raw.get("rt") is not None else None, "scan_id": str(raw.get("scan_id") or "")}


def _group_candidate_points(candidate: dict[str, Any], points: list[dict[str, Any]], max_rt_gap: float) -> list[dict[str, Any]]:
    if not points:
        return []
    ordered = sorted(points, key=lambda x: ((x["rt"] if x["rt"] is not None else -1), x["scan_id"], x["mz"]))
    groups: list[list[dict[str, Any]]] = [[ordered[0]]]
    for point in ordered[1:]:
        previous = groups[-1][-1]
        gap = abs(float(point["rt"] or 0) - float(previous["rt"] or 0))
        if gap <= max_rt_gap:
            groups[-1].append(point)
        else:
            groups.append([point])
    result=[]
    for number, group in enumerate(groups, 1):
        apex=max(group,key=lambda x:x["intensity"]); total=sum(x["intensity"] for x in group)
        centroid=sum(x["mz"]*x["intensity"] for x in group)/total if total else sum(x["mz"] for x in group)/len(group)
        variance=sum(x["intensity"]*(x["mz"]-centroid)**2 for x in group)/total if total else 0.0
        rts=[x["rt"] for x in group if x["rt"] is not None]; theoretical=float(candidate["Theoretical_mz"])
        result.append({"Feature_ID":f"F_{_safe_id(candidate['Chemical_State_ID'])}_{number}",
            "Chemical_State_ID":candidate["Chemical_State_ID"], "Chemical_Family":candidate["Chemical_Family"],
            "Product_Type":candidate["Product_Type"], "Elemental_Composition":candidate["Elemental_Composition"],
            "Charge":candidate["Charge"], "Theoretical_mz":theoretical,
            "RT_Start":min(rts) if rts else None, "RT_End":max(rts) if rts else None,
            "RT_Apex":apex["rt"], "RT_Span":max(rts)-min(rts) if rts else 0.0,
            "Apex_mz":apex["mz"], "mz_Centroid":centroid, "mz_SD":math.sqrt(max(variance,0.0)),
            "Apex_Intensity":apex["intensity"], "Integrated_Intensity":total,
            "Profile_Point_Count":len(group), "Spectrum_Count":len({x["scan_id"] for x in group}),
            "Mass_Error_ppm_at_Apex":(apex["mz"]-theoretical)/theoretical*1e6,
            "Mass_Error_ppm_at_Centroid":(centroid-theoretical)/theoretical*1e6,
            "_point_ids":{x["index"] for x in group}, "_scan_ids":{x["scan_id"] for x in group}, **FALSE_FLAGS})
    return result


def _assign_physical_features(features: list[dict[str, Any]]) -> None:
    groups: list[dict[str, Any]]=[]
    for feature in sorted(features,key=lambda x:(float(x.get("RT_Start") or 0),float(x.get("Apex_mz") or 0))):
        chosen=None
        for group in groups:
            if feature["Charge"] != group["charge"]: continue
            overlap=bool(feature["_point_ids"] & group["points"])
            rt_overlap=float(feature.get("RT_Start") or 0)<=group["rt_end"]+0.08 and float(feature.get("RT_End") or 0)>=group["rt_start"]-0.08
            mz_close=abs(float(feature["Apex_mz"])-group["mz"])<=max(float(feature["Apex_mz"]),group["mz"])*10e-6
            if overlap or (rt_overlap and mz_close): chosen=group;break
        if chosen is None:
            chosen={"id":f"PF_{len(groups)+1:05d}","charge":feature["Charge"],"points":set(),
                "rt_start":float(feature.get("RT_Start") or 0),"rt_end":float(feature.get("RT_End") or 0),"mz":float(feature["Apex_mz"])}
            groups.append(chosen)
        chosen["points"].update(feature["_point_ids"]);chosen["rt_start"]=min(chosen["rt_start"],float(feature.get("RT_Start") or 0));chosen["rt_end"]=max(chosen["rt_end"],float(feature.get("RT_End") or 0))
        feature["Physical_Feature_ID"]=chosen["id"]


def _competition_and_isotopes(features: list[dict[str, Any]], peaks: list[dict[str, Any]], candidates: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    by_physical=defaultdict(list)
    for feature in features:by_physical[feature["Physical_Feature_ID"]].append(feature)
    rt_peaks=sorted(((float(p["rt"]),p) for p in peaks if p["rt"] is not None), key=lambda item:item[0])
    rt_values=[item[0] for item in rt_peaks]
    competition=[]
    for feature in features:
        candidate=candidates[feature["Chemical_State_ID"]]; competitors=[x for x in by_physical[feature["Physical_Feature_ID"]] if x["Chemical_State_ID"]!=feature["Chemical_State_ID"]]
        family_counts=Counter(x["Chemical_Family"] for x in competitors)
        same_comp=sum(x["Elemental_Composition"]==feature["Elemental_Composition"] for x in competitors)
        charge_comp=sum(x["Charge"]!=feature["Charge"] for x in competitors)
        rt=float(feature.get("RT_Apex") or 0);spacing=ISOTOPE_SPACING/abs(int(feature["Charge"]));mz=float(feature["Apex_mz"])
        nearby=[item[1] for item in rt_peaks[bisect_left(rt_values,rt-0.08):bisect_right(rt_values,rt+0.08)]]
        minus=any(abs(p["mz"]-(mz-spacing))<=max(mz-spacing,1)*10e-6 for p in nearby)
        plus=any(abs(p["mz"]-(mz+spacing))<=max(mz+spacing,1)*10e-6 for p in nearby)
        if feature["Profile_Point_Count"] == 1:
            continuity="single_isolated_profile_point";eligible=False;exclusion="single_profile_point_is_not_independent_feature_support"
        elif feature["Spectrum_Count"] < 2:
            continuity="single_spectrum_profile_cluster";eligible=False;exclusion="single_spectrum_cluster_is_not_independent_feature_support"
        elif float(feature.get("RT_Span") or 0) > 1.0:
            continuity="continuous_background_trace";eligible=False;exclusion="RT_span_exceeds_1_minute_chromatographic_feature_limit"
        else:
            continuity="continuous_profile_feature";eligible=True;exclusion=""
        feature.update({"Expected_Isotope_Spacing":spacing,"Observed_Isotope_Spacing":spacing if minus or plus else "",
            "Estimated_Charge":feature["Charge"],"Monoisotopic_Candidate":True,"Minus1_Isotope_Candidate":minus,
            "Plus1_Isotope_Candidate":plus,"Envelope_Assessed":False,"Isotope_Status":"provisional",
            "Feature_Continuity_Status":continuity,"Feature_Eligible_For_Support":eligible,"Feature_Exclusion_Reason":exclusion,
            "Normal_Phosphate_Competition_Count":family_counts["NORMAL_PHOSPHATE"]+family_counts["RESIDUAL_NORMAL_PHOSPHATE"],
            "PT_Competition_Count":family_counts["PHOSPHOROTHIOATE"]+family_counts["P1_RESISTANT_PT_OLIGOMER"],
            "Thiophosphate_Competition_Count":family_counts["THIOPHOSPHATE_LIKE"],
            "Oxidized_PT_Competition_Count":family_counts["OXIDIZED_PT_DERIVATIVE"],
            "Sulfur_Non_PT_Competition_Count":family_counts["SULFUR_CONTAINING_NON_PT_ALTERNATIVE"],
            "Isotope_Competition_Count":int(minus)+int(plus),"Charge_Competition_Count":charge_comp,
            "Composition_Competition_Count":same_comp,"Competition_Count":len(competitors),
            "Candidate_Specific":not competitors})
        pt_like=candidate["Chemical_Family"] in {"PHOSPHOROTHIOATE","THIOPHOSPHATE_LIKE","OXIDIZED_PT_DERIVATIVE","P1_RESISTANT_PT_OLIGOMER"}
        if pt_like:
            interpretation=("PT_CHEMICAL_STATE_SUPPORTED" if not competitors else "PT_LIKE_STATE_AMBIGUOUS") if eligible else "NOT_EVALUABLE"
        elif candidate["Chemical_Family"]=="RESIDUAL_NORMAL_PHOSPHATE":interpretation="NORMAL_PHOSPHATE_ONLY"
        elif candidate["Chemical_Family"] in {"DEPHOSPHORYLATED","SULFUR_CONTAINING_NON_PT_ALTERNATIVE"}:interpretation="DEPHOSPHORYLATED_ONLY"
        else:interpretation="NOT_EVALUABLE"
        feature.update({"Final_Interpretation":interpretation,"Chemical_State_Supported":eligible,
            "Sequence_Position_Localized":False,"Original_Bond_Localized":False,
            "Localization_Status":"POSITION_NOT_RETAINED_BY_PREPARATION"})
        for other in competitors:
            competition.append({"Physical_Feature_ID":feature["Physical_Feature_ID"],"Feature_ID":feature["Feature_ID"],
                "Chemical_State_ID":feature["Chemical_State_ID"],"Competing_State_ID":other["Chemical_State_ID"],
                "Candidate_Chemical_Family":feature["Chemical_Family"],"Competing_Chemical_Family":other["Chemical_Family"],
                "Same_Elemental_Composition":feature["Elemental_Composition"]==other["Elemental_Composition"],
                "Same_Charge":feature["Charge"]==other["Charge"],"Competition_Type":"same_physical_feature",**FALSE_FLAGS})
    return competition


def match_and_group_features(candidates: list[dict[str, Any]], peaks: list[Any], tolerance_ppm: float = 10.0,
                             max_rt_gap: float = 0.08) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, int]]:
    raw=[_raw_peak(p,i) for i,p in enumerate(peaks)];ordered=sorted(raw,key=lambda x:x["mz"]);mzs=[p["mz"] for p in ordered]
    features=[];raw_counts={}
    for candidate in candidates:
        if not candidate.get("Search_Enabled") or candidate.get("Theoretical_mz") is None:
            raw_counts[candidate["Chemical_State_ID"]]=0;continue
        target=float(candidate["Theoretical_mz"]);delta=target*tolerance_ppm/1e6
        points=ordered[bisect_left(mzs,target-delta):bisect_right(mzs,target+delta)]
        raw_counts[candidate["Chemical_State_ID"]]=len(points)
        features.extend(_group_candidate_points(candidate,points,max_rt_gap))
    _assign_physical_features(features)
    competition=_competition_and_isotopes(features,raw,{x["Chemical_State_ID"]:x for x in candidates})
    return features,competition,raw_counts


def _update_candidates(candidates: list[dict[str, Any]], features: list[dict[str, Any]], raw_counts: dict[str,int]) -> None:
    by_state=defaultdict(list)
    for feature in features:by_state[feature["Chemical_State_ID"]].append(feature)
    for candidate in candidates:
        rows=by_state[candidate["Chemical_State_ID"]];eligible=[x for x in rows if x.get("Feature_Eligible_For_Support",False)];best=max(eligible,key=lambda x:x["Integrated_Intensity"],default=None)
        interpretation="NO_PT_LIKE_SUPPORT_IN_CURRENT_DATA" if candidate["Chemical_Family"] in {"PHOSPHOROTHIOATE","P1_RESISTANT_PT_OLIGOMER"} and candidate["Observable"] else "NOT_EVALUABLE"
        if eligible:
            values={x["Final_Interpretation"] for x in eligible};interpretation="PT_LIKE_STATE_AMBIGUOUS" if "PT_LIKE_STATE_AMBIGUOUS" in values else next(iter(values))
        candidate.update({"Raw_Profile_Point_Match_Count":raw_counts.get(candidate["Chemical_State_ID"],0),
            "Grouped_Feature_Count":len({x["Physical_Feature_ID"] for x in rows}),
            "Independent_Feature_Count":len({x["Physical_Feature_ID"] for x in eligible}),
            "Candidate_Specific_Feature_Count":sum(bool(x["Candidate_Specific"]) for x in eligible),
            "Best_Feature_ID":best["Feature_ID"] if best else "","Best_Observed_mz":best["Apex_mz"] if best else "",
            "Best_Mass_Error_ppm":best["Mass_Error_ppm_at_Apex"] if best else "","Best_RT_Apex":best["RT_Apex"] if best else "",
            "Best_Integrated_Intensity":best["Integrated_Intensity"] if best else "",
            "Chemical_State_Supported":bool(eligible),"PT_Chemical_State_Interpretation":interpretation})


def _family_rows(specs: list[dict[str, Any]], candidates: list[dict[str, Any]], features: list[dict[str, Any]]) -> list[dict[str, Any]]:
    cmap={x["Chemical_State_ID"]:x for x in candidates};fmap=defaultdict(list)
    for f in features:
        if f.get("Feature_Eligible_For_Support",True):fmap[f["Chemical_State_ID"]].append(f)
    rows=[]
    for spec in specs:
        counts={key:len({x["Physical_Feature_ID"] for x in fmap.get(spec[state_key],[])}) for key,state_key in (
            ("dephos","Dephosphorylated_State_ID"),("normal","Residual_State_ID"),("pt","PT_State_ID"),
            ("thio","Thiophosphate_State_ID"),("oxidized","Oxidized_PT_State_ID"))}
        pt_features=fmap.get(spec["PT_State_ID"],[]);ambiguous=sum(not x["Candidate_Specific"] for x in pt_features)
        if counts["pt"]:interpretation="PT_LIKE_STATE_AMBIGUOUS" if ambiguous else "PT_RETAINED_AFTER_SAP"
        elif counts["normal"] and counts["dephos"]:interpretation="MIXED_TERMINAL_STATES"
        elif counts["normal"]:interpretation="RESIDUAL_NORMAL_PHOSPHATE_PRESENT"
        elif counts["dephos"]:interpretation="DEPHOSPHORYLATED_DOMINANT"
        else:interpretation="NO_TERMINAL_STATE_SUPPORT"
        row=dict(spec);row.update({"Normal_Feature_Count":counts["normal"],"Dephosphorylated_Feature_Count":counts["dephos"],
            "PT_Feature_Count":counts["pt"],"Thiophosphate_Feature_Count":counts["thio"],
            "Oxidized_PT_Feature_Count":counts["oxidized"],"P1_Resistant_PT_Oligomer_Feature_Count":0,
            "Ambiguous_Sulfur_Feature_Count":ambiguous,"Family_Interpretation":interpretation,
            "Sequence_Position_Localized":False,"Original_Bond_Localized":False,
            "Localization_Status":"POSITION_NOT_RETAINED_BY_PREPARATION",**FALSE_FLAGS})
        rows.append(row)
    return rows


def _terminal_rows(families: list[dict[str, Any]], candidates: list[dict[str, Any]], features: list[dict[str, Any]]) -> list[dict[str, Any]]:
    cmap={x["Chemical_State_ID"]:x for x in candidates};fmap=defaultdict(list)
    for f in features:
        if f.get("Feature_Eligible_For_Support",True):fmap[f["Chemical_State_ID"]].append(f)
    rows=[]
    for family in families:
        state_ids=[family["Dephosphorylated_State_ID"],family["Residual_State_ID"],family["PT_State_ID"]]
        state_features=[x for state in state_ids for x in fmap.get(state,[])];best=max(state_features,key=lambda x:x["Integrated_Intensity"],default=None)
        pt_present=bool(fmap.get(family["PT_State_ID"]));normal_present=bool(fmap.get(family["Residual_State_ID"]));dephos_present=bool(fmap.get(family["Dephosphorylated_State_ID"]))
        if pt_present:interpretation="PT_LIKE_STATE_RETAINED"
        elif normal_present and dephos_present:interpretation="MIXED_TERMINAL_STATES"
        elif normal_present:interpretation="RESIDUAL_NORMAL_PHOSPHATE_PRESENT"
        elif dephos_present:interpretation="DEPHOSPHORYLATED_DOMINANT"
        else:interpretation="NO_TERMINAL_STATE_SUPPORT"
        for substrate,expected,unknown,reason in (
            ("normal phosphate monoester",True,False,"SAP is expected to remove a normal phosphate monoester; no untreated control is available."),
            ("residual normal phosphate",True,False,"Residual normal phosphate is compatible with incomplete SAP reaction or another phosphate-bearing species."),
            ("thiophosphate-like terminal state",False,True,"SAP reactivity is not established for the undefined thiophosphate-like state."),
            ("phosphorothioate monoester",False,True,"SAP removal of a phosphorothioate monoester is unknown in this model."),
            ("internal phosphorothioate",False,True,"An internal phosphorothioate is not modeled as a confirmed SAP substrate."),
            ("dephosphorylated counterpart",False,False,"No phosphate remains for SAP removal."),
        ):
            rows.append({"Chemical_Base_State":family["Chemical_Base_State"],"SAP_Substrate_State":substrate,
                "Dephosphorylated_Feature_Present":dephos_present,"Residual_Phosphate_Feature_Present":normal_present,
                "PT_Like_Feature_Present":pt_present,"SAP_Interpretation":interpretation,
                "SAP_Removal_Expected":expected,"SAP_Removal_Confirmed":False,"SAP_Removal_Unknown":unknown,
                "SAP_Model_Reason":reason,"Observable":any(cmap[x]["Observable"] for x in state_ids),
                "Feature_Present":bool(state_features),"Feature_ID":best["Feature_ID"] if best else "",
                "Observed_mz":best["Apex_mz"] if best else "","Mass_Error_ppm":best["Mass_Error_ppm_at_Apex"] if best else "",
                "RT_Apex":best["RT_Apex"] if best else "","Integrated_Intensity":best["Integrated_Intensity"] if best else "",
                "Candidate_Specific":best["Candidate_Specific"] if best else False,
                "Sequence_Position_Localized":False,"Original_Bond_Localized":False,
                "Localization_Status":"POSITION_NOT_RETAINED_BY_PREPARATION",**FALSE_FLAGS})
    return rows


def _cross_enzyme_rows(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows=[]
    for candidate in candidates:
        if candidate["Chemical_Family"] not in {"PHOSPHOROTHIOATE","P1_RESISTANT_PT_OLIGOMER"}:continue
        p1="FEATURE_LEVEL_SUPPORT" if candidate.get("Chemical_State_Supported") else "NO_FEATURE_LEVEL_SUPPORT"
        rows.append({"Chemical_State_ID":candidate["Chemical_State_ID"],"Elemental_Composition":candidate["Elemental_Composition"],
            "T1_Support_Status":"NOT_COMPARABLE_FROM_POSITION_PRIOR_ONLY","P1_SAP_Support_Status":p1,
            "Cross_Enzyme_Chemical_Status":"NOT_COMPARABLE","Position_Confirmed":False,
            "Reason":"T1 manifest evidence is position/fragment based or absent; P1+SAP chemistry does not retain position.",**FALSE_FLAGS})
    return rows


def _precursor(spectrum: dict[str, Any]) -> tuple[float|None,int|None,str,str]:
    precursors=((spectrum.get("precursorList") or {}).get("precursor") or [])
    if not precursors:return None,None,"",""
    precursor=precursors[0];ions=((precursor.get("selectedIonList") or {}).get("selectedIon") or []);ion=ions[0] if ions else {}
    mz=ion.get("selected ion m/z");charge=ion.get("charge state")
    isolation=precursor.get("isolationWindow") or {};activation=precursor.get("activation") or {}
    window=";".join(f"{k}={v}" for k,v in isolation.items() if "isolation window" in str(k))
    energy=activation.get("collision energy","")
    try:mz=float(mz) if mz is not None else None
    except (TypeError,ValueError):mz=None
    try:charge=int(charge) if charge is not None else None
    except (TypeError,ValueError):charge=None
    return mz,charge,window,str(energy)


def compatible_ms2_provenance(mzml_path: str|Path, candidates: list[dict[str, Any]], tolerance_ppm: float=20.0) -> list[dict[str, Any]]:
    searchable=[x for x in candidates if x.get("Search_Enabled") and x.get("Theoretical_mz") is not None]
    targets=sorted(((float(x["Theoretical_mz"]),x) for x in searchable), key=lambda item:item[0]);mzs=[x[0] for x in targets];rows=[];counts=Counter()
    for spectrum in iter_spectra(mzml_path):
        if int(spectrum.get("ms level",0) or 0)!=2:continue
        observed,charge,window,energy=_precursor(spectrum)
        if observed is None:continue
        delta=observed*tolerance_ppm/1e6
        for theoretical,candidate in targets[bisect_left(mzs,observed-delta):bisect_right(mzs,observed+delta)]:
            if charge not in (None,0,candidate["Charge"]):continue
            counts[candidate["Chemical_State_ID"]]+=1
            rows.append({"Chemical_State_ID":candidate["Chemical_State_ID"],"Precursor_Compatible_MS2_Count":"",
                "MS2_Spectrum_ID":str(spectrum.get("id") or ""),"Spectrum_RT":_rt_minutes(spectrum),
                "Isolation_Window":window,"Collision_Energy":energy,"Precursor_mz":observed,
                "Precursor_Charge":charge,"Precursor_Mass_Error_ppm":(observed-theoretical)/theoretical*1e6,
                "MS2_Model_Applicable":False,"MS2_Model_Not_Applicable_Reason":"P1_SAP_small_product_requires_dedicated_fragment_model",**FALSE_FLAGS})
    for row in rows:row["Precursor_Compatible_MS2_Count"]=counts[row["Chemical_State_ID"]]
    return rows


def build_p1_sap_chemical_state_audit(project_root: str|Path, sequence: str, peaks: list[Any], config: Any,
                                       modifications: list[Any], *, audit_level: str="audit",
                                       mzml_path: str|Path|None=None) -> P1SAPAuditResult:
    started=time.perf_counter();tracemalloc.start();root=Path(project_root)
    candidates,family_specs=generate_chemical_state_candidates(sequence,modifications,root,charges=(1,))
    mz_min=float((getattr(config,"reconstruction",{}) or {}).get("mz_min",0));mz_max=float((getattr(config,"reconstruction",{}) or {}).get("mz_max",float("inf")))
    for candidate in candidates:
        if candidate["Theoretical_mz"] is not None:
            candidate["Observable"]=mz_min<=float(candidate["Theoretical_mz"])<=mz_max
            candidate["Not_Observable_Reason"]="" if candidate["Observable"] else "outside_acquisition_range"
            candidate["Search_Enabled"]=candidate["Search_Enabled"] and candidate["Observable"]
    tolerance=float((getattr(config,"p1_annotation",{}) or {}).get("mz_tolerance_ppm") or (getattr(config,"instrument",{}) or {}).get("ms1_tolerance_ppm",10) or 10)
    raw_peaks=[_raw_peak(p,i) for i,p in enumerate(peaks)]
    features,competition,raw_counts=match_and_group_features(candidates,peaks,tolerance_ppm=tolerance)
    quality_results=build_p1_sap_feature_quality(candidates, features, raw_peaks, config, tolerance)
    dinucleotide_audit=build_p1_sap_dinucleotide_audit(root, sequence, peaks, config, audit_level=audit_level, mzml_path=mzml_path)
    _update_candidates(candidates,features,raw_counts);families=_family_rows(family_specs,candidates,features);terminal=_terminal_rows(family_specs,candidates,features);cross=_cross_enzyme_rows(candidates)
    for feature in features:
        feature.pop("_point_ids",None); feature.pop("_scan_ids",None)
    ms2=compatible_ms2_provenance(mzml_path,candidates,20.0) if audit_level=="full" and mzml_path else []
    pt_features=[x for x in features if x["Chemical_Family"] in {"PHOSPHOROTHIOATE","P1_RESISTANT_PT_OLIGOMER"}]
    supported_pt_features=[x for x in pt_features if x.get("Feature_Eligible_For_Support",False)]
    if supported_pt_features:
        final="PT_LIKE_STATE_AMBIGUOUS" if any(not x["Candidate_Specific"] for x in supported_pt_features) else "PT_RETAINED_AFTER_SAP"
    else:final="NO_PT_LIKE_SUPPORT_IN_CURRENT_DATA" if any(x["Chemical_Family"]=="PHOSPHOROTHIOATE" and x["Observable"] for x in candidates) else "NOT_EVALUABLE"
    eligible_features=[x for x in features if x.get("Feature_Eligible_For_Support",False)]
    physical={x["Physical_Feature_ID"] for x in eligible_features};raw_total=sum(raw_counts.values());runtime=time.perf_counter()-started;_,peak_bytes=tracemalloc.get_traced_memory();tracemalloc.stop()
    maximum_rss_mib=get_maximum_rss_mib()
    count_family=lambda family:len({x["Physical_Feature_ID"] for x in eligible_features if x["Chemical_Family"]==family})
    summary={"Configured_Enzyme":normalize_enzyme_name((getattr(config,"digestion",{}) or {}).get("enzyme","")),
        "Rule_ID":"Nuclease_P1","Cleavage_Specificity":"nonspecific_all_standard_RNA_bonds","Cleavage_Side":"3prime_of_each_residue",
        "Expected_Product_Type":"nucleoside_or_nucleotide_monomer_after_SAP","Expected_5prime_Terminal_State":"preparation_dependent",
        "Expected_3prime_Terminal_State":"hydroxyl_or_preparation_dependent","Expected_Phosphate_Product":"5prime_nucleotide_before_SAP;dephosphorylated_counterpart_after_SAP",
        "Missed_Cleavage_Model":"formal_complete;P1-resistant_PT_oligomer_shadow_max_length_2",
        "P1_Cleaves_Normal_Phosphodiester":True,"P1_PT_Cleavage_Behavior":"unknown",
        "Candidate_Count":len(candidates),"Searchable_Candidate_Count":sum(bool(x["Search_Enabled"]) for x in candidates),
        "Raw_Profile_Point_Match_Count":raw_total,"Grouped_Feature_Count":len(features),"Independent_Feature_Count":len(physical),
        "Feature_Apex_Match_Count":len(physical),"Dephosphorylated_Feature_Count":count_family("DEPHOSPHORYLATED")+count_family("SULFUR_CONTAINING_NON_PT_ALTERNATIVE"),
        "Residual_Normal_Phosphate_Feature_Count":count_family("RESIDUAL_NORMAL_PHOSPHATE"),"PT_Feature_Count":count_family("PHOSPHOROTHIOATE"),
        "Thiophosphate_Like_Feature_Count":count_family("THIOPHOSPHATE_LIKE"),"Oxidized_PT_Like_Feature_Count":count_family("OXIDIZED_PT_DERIVATIVE"),
        "P1_Resistant_PT_Oligomer_Feature_Count":count_family("P1_RESISTANT_PT_OLIGOMER"),
        "Ambiguous_Sulfur_Feature_Count":len({x["Physical_Feature_ID"] for x in eligible_features if not x["Candidate_Specific"] and (x["Chemical_Family"] in {"PHOSPHOROTHIOATE","P1_RESISTANT_PT_OLIGOMER","SULFUR_CONTAINING_NON_PT_ALTERNATIVE"})}),
        "PT_Chemical_State_Final_Interpretation":final,"Position_Localization_Disabled":True,"Bond_Localization_Disabled":True,
        "Audit_Runtime":runtime,"Maximum_RSS_MiB":maximum_rss_mib,"Tracemalloc_Peak_MiB":peak_bytes/(1024*1024),
        "Feature_Detail_Row_Count":len(features) if audit_level=="full" else 0,**FALSE_FLAGS}
    sheets={"P1_SAP_Chemical_State":candidates,"P1_SAP_PT_Family":families,"P1_SAP_Terminal_Audit":terminal,"P1_SAP_Chemistry_Summary":[summary],
        "P1_SAP_Spectrum_Peaks":quality_results["spectrum_peaks"],
        "P1_SAP_Refined_Features":quality_results["refined_features"],
        "P1_SAP_Feature_Quality":quality_results["quality_rows"],
        "P1_SAP_Isotope_Audit":quality_results["isotope_rows"],
        "P1_SAP_Quality_Summary":[quality_results["summary_row"]]}
    if audit_level=="full":sheets.update({"P1_SAP_Features":features,"P1_SAP_Competition":competition,"Cross_Enzyme_Chemistry":cross,"P1_SAP_MS2_Provenance":ms2})
    sheets.update(dinucleotide_audit["sheets"])
    metrics=dict(summary)
    payload={"experiment_treatment":["Nuclease P1 digestion","SAP alkaline phosphatase","SCIEX ZenoTOF positive profile"],
        "p1_cleavage_model":{k:summary[k] for k in ("Configured_Enzyme","Rule_ID","Cleavage_Specificity","Cleavage_Side","Expected_Product_Type","Expected_5prime_Terminal_State","Expected_3prime_Terminal_State","Expected_Phosphate_Product","Missed_Cleavage_Model","P1_Cleaves_Normal_Phosphodiester","P1_PT_Cleavage_Behavior")},
        "sap_reaction_model":{"normal_phosphate_removal_expected":True,"normal_phosphate_removal_confirmed":False,"PT_removal_unknown":True,"thiophosphate_removal_unknown":True,"complete_reaction_asserted":False},
        "candidate_counts":{"total":len(candidates),"searchable":summary["Searchable_Candidate_Count"]},
        "profile_match_counts":{"raw":raw_total,"grouped_candidate_features":len(features),"independent_physical_features":len(physical)},
        "feature_family_counts":{k:summary[k] for k in summary if k.endswith("Feature_Count")},
        "PT_like_features":[{k:x.get(k) for k in FEATURE_COLUMNS if not k.startswith("Applied_") and not k.startswith("Formal_")} for x in pt_features],
        "terminal_state_audit":terminal,"position_localization_disabled":True,"bond_localization_disabled":True,
        "feature_quality":{"summary":quality_results["summary_row"],"features":quality_results["quality_rows"],"isotope_audit":quality_results["isotope_rows"]},
        "cross_enzyme_comparison":cross,"formal_flags":FALSE_FLAGS,"performance":{"audit_runtime_seconds":runtime,"tracemalloc_peak_mib":peak_bytes/(1024*1024),"feature_detail_rows":summary["Feature_Detail_Row_Count"]},
        "final_interpretation":final}
    payload.update(dinucleotide_audit["payload"])
    payload["performance"]["maximum_rss_mib"]=maximum_rss_mib
    return P1SAPAuditResult(sheets,metrics,payload)


def _hash_file(path: str|Path) -> str:
    h=sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()


def write_p1_sap_summary_json(result: P1SAPAuditResult, output_dir: str|Path, *, mzml_path: str|Path,
                              config_path: str|Path, original_config_path: str|Path,
                              audit_level: str, excel_path: str|Path|None=None) -> Path:
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=True);stamp=time.strftime("%Y%m%d_%H%M%S")
    path=out/f"Nsd01_P1_SAP_PT_chemical_state_summary_{stamp}.json"
    try:commit=subprocess.run(["git","rev-parse","HEAD"],capture_output=True,text=True,check=True).stdout.strip()
    except Exception:commit="unknown"
    payload={"schema_version":1,"audit_level":audit_level,"input_mzML_identity":{"path":str(mzml_path),"size_bytes":Path(mzml_path).stat().st_size},
        "git_commit":commit,"config_sha256":_hash_file(original_config_path),"alternate_config_sha256":_hash_file(config_path),
        "excel_path":str(excel_path or ""),**result.summary_payload}
    path.write_text(json.dumps(payload,ensure_ascii=False,indent=2,default=str)+"\n",encoding="utf-8")
    return path
