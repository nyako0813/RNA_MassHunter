"""Orchestration and Excel rows for shadow-only PT paired evidence audits."""
from __future__ import annotations
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any
import yaml
from rna_masshunter.backbone_state import load_backbone_transformations
from rna_masshunter.cleavage_site_discovery import discovery_candidate_row
from rna_masshunter.composite_ms2_matcher import match_composite_ms2
from rna_masshunter.composite_ms2_propagation import generate_composite_theoretical_ions
from rna_masshunter.enzymes import normalize_enzyme_name
from rna_masshunter.masses import mz_from_neutral_mass
from rna_masshunter.models import Fragment
from rna_masshunter.modification_constraints import load_transformations
from rna_masshunter.phosphorothioate_evidence import build_pt_evidence
from rna_masshunter.phosphorothioate_pairing import (
    PTPairSpec, build_pt_pair, discovery_pair_specs, load_pt_pair_hypotheses,
)

STATUS_COLUMNS = ["Applied_To_Formal_Result", "Formal_Change_Ready", "Formal_Result_Changed"]
PT_PAIRED_EVIDENCE_COLUMNS = [
    "Audit_Level","Candidate_ID","Hypothesis_ID","Search_Mode","Sequence_ID","Enzyme","Bond_ID",
    "Left_Position","Right_Position","Left_Base","Right_Base","Is_Normal_Cleavage_Site","Fragment_ID",
    "Fragment_Start","Fragment_End","Terminal_Form","Charge","Shared_Nucleoside_States",
    "Shared_Modified_Positions","Applied_Transformations","Normal_Backbone_State","Modified_Backbone_State",
    "Normal_Composition","Modified_Composition","Composition_Delta","Shared_Modification_Composition_Delta",
    "Backbone_Composition_Delta","Normal_Neutral_Mass","Modified_Neutral_Mass","Neutral_Mass_Delta",
    "Expected_O_to_S_Delta","Delta_Consistency_Error","Normal_Theoretical_mz","Modified_Theoretical_mz",
    "Theoretical_mz_Delta","Normal_Observed_mz","Modified_Observed_mz","Normal_Mass_Error_Da",
    "Modified_Mass_Error_Da","Normal_Mass_Error_ppm","Modified_Mass_Error_ppm","Normal_Intensity",
    "Modified_Intensity","Normal_Scan","Modified_Scan","Normal_RT","Modified_RT","Normal_Physical_Peak_ID",
    "Modified_Physical_Peak_ID","Same_Physical_Peak","Normal_Competition_Count","Modified_Competition_Count",
    "Normal_Competing_Candidate_IDs","Modified_Competing_Candidate_IDs","Peak_Ambiguity_Status",
    "Normal_Cleavage_Mechanism","Modified_Cleavage_Mechanism","Nucleoside_Blocking_Status","Cleavage_Status",
    "Evidence_Class","Evidence_Reason","Mechanism_Discriminating","Position_Localizing","Candidate_Specific",
    "Observable", *STATUS_COLUMNS,
]
PT_STATE_SEARCH_COLUMNS = [
    "Audit_Level","Candidate_ID","Hypothesis_ID","Search_Mode","Candidate_State_ID","Sequence_ID","Enzyme",
    "Bond_ID","Shared_Nucleoside_States","Shared_Modified_Positions","Applied_Transformations","Backbone_State",
    "Fragment_Start","Fragment_End","Terminal_Form","Charge","Elemental_Composition","Neutral_Mass",
    "Theoretical_mz","Nearest_Observed_mz","Mass_Error_Da","Mass_Error_ppm","Intensity","Scan","RT",
    "Physical_Peak_ID","Competition_Count","Competing_Candidate_IDs","Observable","Matched", *STATUS_COLUMNS,
]
PT_DISCOVERY_COLUMNS = [
    "Sequence_ID","Enzyme","Bond_ID","Left_Position","Right_Position","Left_Base","Right_Base",
    "Is_Normal_Cleavage_Site","Fragment_Context","Candidate_Nucleoside_State_Count",
    "Candidate_Backbone_State_Count","Pair_Count","Position_Compatible","Pathway_Compatible",
    "Generated_For_Evaluation","Invalid_Reason", *STATUS_COLUMNS,
]
PT_SUMMARY_COLUMNS = [
    "Audit_Level","Search_Mode","Target_Hypothesis_Count","Valid_Hypothesis_Count","Invalid_Hypothesis_Count",
    "Evaluated_Bond_Count","Evaluated_Nucleoside_State_Count","Neutral_Pair_Count","PT_Pair_Count",
    "Observable_Pair_Count","PT_Strong_Paired_Support_Count","PT_Candidate_Specific_MS1_Support_Count",
    "PT_Only_Support_Count",
    "Normal_Only_Support_Count","Both_Present_Count","Neither_Present_Count",
    "Mass_Shift_Inconsistent_Count","Ambiguous_Count","Not_Observable_Count","Candidate_Specific_Count",
    "Candidate_Specific_Physical_Peak_Count","Physical_Peak_Competition_Count","Filtered_Invalid_Combination_Count",
    "Total_MS2_Spectra","Precursor_Compatible_Spectra","PT_Ion_Matched_Spectra","PT_MS2_Matched_Rows",
    "Position_Localizing_Ion_Count","Backbone_Localizing_Ion_Count","MS2_Reason_For_Zero",
    "PT_Audit_Runtime_Seconds","PT_Audit_Peak_Memory_MiB", *STATUS_COLUMNS,
]

@dataclass(frozen=True)
class PTPairedAuditResult:
    sheets: dict[str, list[dict[str, Any]]]
    pairs: tuple[Any, ...]
    invalid_rows: tuple[dict[str, Any], ...]


def _precursor_compatible_count(ions: list[dict[str, Any]], spectra: list[Any], config: Any) -> int:
    ms2 = getattr(config, "ms2_annotation", {}) or {}; tolerance = float(ms2.get("precursor_match_tolerance_ppm", 20) or 20)
    polarity = str((getattr(config, "instrument", {}) or {}).get("polarity") or "negative").lower(); compatible = set()
    for spectrum in spectra or ():
        mz = getattr(spectrum, "precursor_mz", None); charge = getattr(spectrum, "precursor_charge", None)
        if mz in (None, "") or charge in (None, "", 0): continue
        z = abs(int(charge))
        for ion in ions:
            theoretical = mz_from_neutral_mass(float(ion["Parent_Neutral_Mass"]), z, polarity)
            error = (float(mz) - theoretical) / theoretical * 1_000_000 if theoretical else 0.0
            if abs(error) <= tolerance:
                compatible.add(getattr(spectrum, "spectrum_id", "")); break
    return len(compatible)


def _summaries(evidence: list[dict[str, Any]], pairs: list[Any], invalid_count: int,
    hypothesis_count: int, discovery_invalid: int, spectra_count: int, ms2: dict[str, Any], audit_level: str) -> list[dict[str, Any]]:
    rows = []
    for mode in ("hypothesis_driven", "discovery"):
        group = [r for r in evidence if r["Search_Mode"] == mode]; mode_pairs = [p for p in pairs if p.spec.search_mode == mode]
        classes = [r["Evidence_Class"] for r in group]
        physical = {r["Modified_Physical_Peak_ID"] for r in group if r.get("Candidate_Specific") and r.get("Modified_Physical_Peak_ID")}
        rows.append({"Audit_Level": audit_level, "Search_Mode": mode,
            "Target_Hypothesis_Count": hypothesis_count if mode == "hypothesis_driven" else 0,
            "Valid_Hypothesis_Count": len({p.spec.hypothesis_id for p in mode_pairs if p.spec.hypothesis_id}) if mode == "hypothesis_driven" else 0,
            "Invalid_Hypothesis_Count": invalid_count if mode == "hypothesis_driven" else 0,
            "Evaluated_Bond_Count": len({p.spec.bond_id for p in mode_pairs}),
            "Evaluated_Nucleoside_State_Count": len({(p.spec.bond_id, tuple((x.position,x.transform_ids) for x in p.spec.position_states)) for p in mode_pairs}),
            "Neutral_Pair_Count": len(mode_pairs), "PT_Pair_Count": len(group),
            "Observable_Pair_Count": sum(bool(r["Observable"]) for r in group),
            "PT_Strong_Paired_Support_Count": classes.count("PT_STRONG_PAIRED_SUPPORT"),
            "PT_Candidate_Specific_MS1_Support_Count": classes.count("PT_CANDIDATE_SPECIFIC_MS1_SUPPORT"),
            "PT_Only_Support_Count": classes.count("PT_ONLY_SUPPORT"),
            "Normal_Only_Support_Count": classes.count("NORMAL_ONLY_SUPPORT"),
            "Both_Present_Count": classes.count("BOTH_PRESENT"), "Neither_Present_Count": classes.count("NEITHER_PRESENT"),
            "Mass_Shift_Inconsistent_Count": classes.count("MASS_SHIFT_INCONSISTENT"),
            "Ambiguous_Count": classes.count("AMBIGUOUS_PEAK_ASSIGNMENT"),
            "Not_Observable_Count": classes.count("NOT_OBSERVABLE"),
            "Candidate_Specific_Count": sum(bool(r["Candidate_Specific"]) for r in group),
            "Candidate_Specific_Physical_Peak_Count": len(physical),
            "Physical_Peak_Competition_Count": sum(bool(r["Normal_Competition_Count"] or r["Modified_Competition_Count"]) for r in group),
            "Filtered_Invalid_Combination_Count": discovery_invalid if mode == "discovery" else invalid_count,
            "Total_MS2_Spectra": spectra_count, "Precursor_Compatible_Spectra": ms2["compatible"] if mode == "hypothesis_driven" else 0,
            "PT_Ion_Matched_Spectra": ms2["matched_spectra"] if mode == "hypothesis_driven" else 0,
            "PT_MS2_Matched_Rows": ms2["matched_rows"] if mode == "hypothesis_driven" else 0,
            "Position_Localizing_Ion_Count": ms2["position"] if mode == "hypothesis_driven" else 0,
            "Backbone_Localizing_Ion_Count": ms2["backbone"] if mode == "hypothesis_driven" else 0,
            "MS2_Reason_For_Zero": ms2["reason"] if mode == "hypothesis_driven" else "not_evaluated_for_discovery_summary",
            "PT_Audit_Runtime_Seconds": "", "PT_Audit_Peak_Memory_MiB": "",
            "Applied_To_Formal_Result": False, "Formal_Change_Ready": False, "Formal_Result_Changed": False})
    return rows


def build_pt_paired_audit(project_root: str | Path, sequence: str, sequence_id: str, peaks: list[Any],
    spectra: list[Any], config: Any, *, legacy_matches: list[Any] = (), other_composite_matches: list[Any] = (),
    audit_level: str = "full", fixture_path: str | Path | None = None) -> PTPairedAuditResult:
    configured_enzyme = normalize_enzyme_name((getattr(config, "digestion", {}) or {}).get("enzyme", ""))
    if configured_enzyme and configured_enzyme != "RNase_T1":
        ms2 = {"compatible": 0, "matched_spectra": 0, "matched_rows": 0, "position": 0, "backbone": 0,
            "reason": "RNase_T1_hypothesis_on_non_RNase_T1_run"}
        summary = _summaries([], [], 0, 1, 0, len(spectra or ()), ms2, audit_level)
        for row in summary:
            row.update({"Configured_Enzyme": configured_enzyme, "Enzyme_Context_Applicable": False,
                "Evidence_Not_Applicable_Reason": "RNase_T1_hypothesis_on_RNase_P1_run" if configured_enzyme == "Nuclease_P1" else "RNase_T1_hypothesis_on_other_enzyme_run"})
        return PTPairedAuditResult({"PT_Paired_Summary": summary, "PT_Discovery_Candidates": []}, (), ())
    root = Path(project_root); transforms = load_transformations(root / "data/modification_transforms_v2.yaml")
    backbone_transform = load_backbone_transformations(root / "data/backbone_modifications.yaml")[0]
    fixture = Path(fixture_path) if fixture_path else root / "data/sample_pt_pair_hypotheses.yaml"
    organism = str((getattr(config, "organism", {}) or {}).get("species", "")); rule_set = str((getattr(config, "organism", {}) or {}).get("rule_set", ""))
    loaded = load_pt_pair_hypotheses(fixture, sequence=sequence, sequence_id=sequence_id,
        organism=organism, rule_set=rule_set, transformations=transforms)
    hypothesis_specs = []
    for primary in loaded.specs:
        hypothesis_specs.append(replace(primary, candidate_id=primary.candidate_id + "|unmodified", position_states=()))
        hypothesis_specs.append(replace(primary, candidate_id=primary.candidate_id + "|modified"))
    discovery_specs, discovery_meta = discovery_pair_specs(sequence, sequence_id, "RNase_T1", transforms,
        root / "data/nucleoside_slots.yaml", organism=organism, pathway_context=rule_set,
        max_nucleoside_modifications=1)
    pairs = []; invalid = list(loaded.invalid_rows)
    for spec in hypothesis_specs + discovery_specs:
        try: pairs.append(build_pt_pair(spec, transforms, root / "data/nucleoside_slots.yaml", backbone_transform))
        except (KeyError, ValueError) as exc:
            invalid.append({"Candidate_ID": spec.candidate_id, "Valid": False, "Invalid_Reason": str(exc),
                "Invalid_Detail": "pair_construction", "Applied_To_Formal_Result": False,
                "Formal_Change_Ready": False, "Formal_Result_Changed": False})
    evidence, states = build_pt_evidence(pairs, peaks, config, legacy_matches=legacy_matches,
        other_composite_matches=other_composite_matches, audit_level=audit_level, include_detail=audit_level == "full")
    discovery_rows = [discovery_candidate_row(item["candidate"], nucleoside_state_count=item["state_count"],
        backbone_state_count=2, pair_count=item["pair_count"], generated=True) for item in discovery_meta]
    # Targeted modified/PT structures only; lack of a compatible spectrum is reported, not penalized.
    primary_pairs = [p for p in pairs if p.spec.search_mode == "hypothesis_driven" and p.spec.position_states]
    ions = []; ms2_matches = []
    for pair in primary_pairs:
        parent = Fragment(pair.spec.candidate_id + "|parent", sequence_id,
            sequence[pair.spec.fragment_start - 1:pair.spec.fragment_end], pair.spec.fragment_start,
            pair.spec.fragment_end, pair.spec.fragment_start, pair.spec.fragment_end, pair.spec.enzyme, 1,
            pair.spec.terminal_form, pair.normal_fragment.neutral_exact_mass)
        ions.extend(generate_composite_theoretical_ions([pair.modified_structure], [parent], sequence, config, audit_level=audit_level))
    if ions: ms2_matches = match_composite_ms2(ions, spectra, config, audit_level=audit_level)
    compatible = _precursor_compatible_count(ions, spectra, config) if ions else 0
    ms2 = {"compatible": compatible, "matched_spectra": len({r["Spectrum_ID"] for r in ms2_matches}),
        "matched_rows": len(ms2_matches), "position": sum(bool(r.get("Position_Informative")) for r in ms2_matches),
        "backbone": sum(bool(r.get("Backbone_Informative")) for r in ms2_matches),
        "reason": "no_precursor_compatible_spectra" if not compatible else ("no_fragment_peak_match" if not ms2_matches else "matched")}
    summary = _summaries(evidence, pairs, len(loaded.invalid_rows), len((yaml.safe_load(fixture.read_text(encoding='utf-8')) or {}).get('hypotheses') or ()),
        sum(int(x["invalid_count"]) for x in discovery_meta), len(spectra or ()), ms2, audit_level)
    sheets = {"PT_Paired_Summary": summary, "PT_Discovery_Candidates": discovery_rows}
    if audit_level == "full":
        sheets.update({"PT_Paired_Evidence": evidence, "PT_State_Search": states})
    return PTPairedAuditResult(sheets, tuple(pairs), tuple(invalid))
