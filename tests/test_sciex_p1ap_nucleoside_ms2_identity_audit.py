from dataclasses import replace
from pathlib import Path

import pytest

from rna_masshunter.sciex_p1ap_nucleoside_ms2_identity_audit import (
    MS2CandidateRecord, NucleosideProductIonHypothesis, ProcessedMS2Spectrum,
    PrecursorCompatibleMS2Record, audit_optional_result,
    compare_candidate_product_ions, decode_selected_ms2_spectra,
    generate_nucleoside_product_ion_hypotheses,
    match_product_ions_to_ms2_peaks, reconcile_p1ap_ms1_and_ms2_evidence,
    select_precursor_compatible_ms2_spectra, summarize_candidate_ms2_identity_evidence,
    _recurrence_and_energy, _spectrum_summaries, P1APNucleosideMS2AuditResult,
    P1APNucleosideMS2Summary,
)


def candidate(candidate_id="A", mz=268.104, formula="C10H13N5O4", candidate_class="NEUTRAL_NUCLEOSIDE", ambiguity="UNAMBIGUOUS"):
    return MS2CandidateRecord(
        candidate_id=candidate_id, candidate_name=candidate_id, parent_base="A",
        candidate_class=candidate_class, molecular_formula=formula,
        theoretical_neutral_mass=mz - 1.007276466812, theoretical_precursor_mz=mz,
        observed_ms1_peak_id=f"P_{candidate_id}", observed_ms1_mz=mz + 0.0002,
        ms1_mass_error=0.0002, ms1_identity_ambiguity_status=ambiguity,
        ms1_state_family_id="F1",
    )


def spectrum(spectrum_id="s1", selected=268.104, target=268.104, lower=0.5, upper=0.5, energy=25, arrays=None, polarity="positive"):
    ion = {} if selected is None else {"selected ion m/z": selected, "charge state": 1}
    isolation = {}
    if target is not None: isolation["isolation window target m/z"] = target
    if lower is not None: isolation["isolation window lower offset"] = lower
    if upper is not None: isolation["isolation window upper offset"] = upper
    activation = {} if energy is None else {"collision energy": energy, "unitName": "eV"}
    row = {
        "id": spectrum_id, "ms level": 2, f"{polarity} scan": True,
        "centroid spectrum": True, "defaultArrayLength": 3,
        "precursorList": {"precursor": [{"selectedIonList": {"selectedIon": [ion]},
            "isolationWindow": isolation, "activation": activation}]},
        "scanList": {"scan": [{"scan start time": 1.0, "unitName": "minute"}]},
    }
    if arrays:
        row["m/z array"], row["intensity array"] = arrays
    return row


def precursor(candidate_id="A", spectrum_id="s1", energy=25):
    return PrecursorCompatibleMS2Record(
        candidate_id=candidate_id, candidate_name=candidate_id, ms2_spectrum_id=spectrum_id,
        ms2_scan_time=1.0, selected_ion_mz=268.104, isolation_target_mz=268.104,
        isolation_lower_offset=0.5, isolation_upper_offset=0.5,
        isolation_lower_bound=267.604, isolation_upper_bound=268.604,
        candidate_theoretical_mz=268.104, candidate_observed_ms1_mz=268.1042,
        selected_ion_delta=0, isolation_contains_theoretical_mz=True,
        isolation_contains_observed_ms1_mz=True, precursor_charge=1,
        collision_energy=energy, collision_energy_unit="eV", default_array_length=3,
        precursor_compatibility_status="BOTH_SELECTED_ION_AND_ISOLATION_SUPPORT",
        precursor_compatibility_confidence="HIGH", precursor_block_reasons=(),
    )


def processed(spectrum_id="s1", peaks=((136.06177, 100.0), (200.0, 10.0))):
    total = sum(x[1] for x in peaks)
    return ProcessedMS2Spectrum(
        ms2_spectrum_id=spectrum_id, scan_time=1.0, raw_peak_count=len(peaks),
        positive_intensity_peak_count=len(peaks), zero_intensity_peak_count=0,
        negative_intensity_peak_count=0, filtered_peak_count=len(peaks),
        base_peak_mz=max(peaks, key=lambda x:x[1])[0], base_peak_intensity=max(x[1] for x in peaks),
        tic=total, mz_min=min(x[0] for x in peaks), mz_max=max(x[0] for x in peaks),
        profile_or_centroid_metadata="CENTROID", ms2_preprocessing_status="COMPLETED",
        ms2_preprocessing_block_reasons=(), peaks=tuple(peaks),
    )


def test_selected_ion_within_tolerance():
    rows=select_precursor_compatible_ms2_spectra(Path("x"),[candidate()],metadata_spectra=[spectrum(target=300)])
    assert rows[0].precursor_compatibility_status=="SELECTED_ION_MZ_WITHIN_TOLERANCE"


def test_isolation_window_containment_without_selected_support():
    rows=select_precursor_compatible_ms2_spectra(Path("x"),[candidate()],metadata_spectra=[spectrum(selected=270,target=268.2)])
    assert rows[0].precursor_compatibility_status=="ISOLATION_WINDOW_CONTAINS_CANDIDATE"


def test_missing_isolation_metadata_lowers_confidence():
    rows=select_precursor_compatible_ms2_spectra(Path("x"),[candidate()],metadata_spectra=[spectrum(target=None,lower=None,upper=None)])
    assert rows[0].precursor_compatibility_confidence=="LOW"
    assert "MISSING_ISOLATION_WINDOW" in rows[0].precursor_block_reasons


def test_candidate_incompatible_precursor_is_excluded():
    assert not select_precursor_compatible_ms2_spectra(Path("x"),[candidate()],metadata_spectra=[spectrum(selected=500,target=500)])


def test_negative_ms2_is_excluded():
    assert not select_precursor_compatible_ms2_spectra(Path("x"),[candidate()],metadata_spectra=[spectrum(polarity="negative")])


class LazyArray:
    def __init__(self, values, allowed=True): self.values=values; self.allowed=allowed; self.calls=0
    def decode(self):
        self.calls+=1
        if not self.allowed: raise AssertionError("unselected binary decoded")
        return self.values


def test_selected_ms2_only_decode_and_unique_reuse():
    a=LazyArray([100,136.06177]); b=LazyArray([1,100]); trap=LazyArray([],False)
    source=[spectrum("s1",arrays=(a,b)),spectrum("s2",arrays=(trap,trap))]
    rows=decode_selected_ms2_spectra(Path("x"),["s1","s1"],spectrum_source=source)
    assert len(rows)==1 and a.calls==1 and b.calls==1 and trap.calls==0


def test_zero_and_negative_intensity_excluded_but_counted():
    row=decode_selected_ms2_spectra(Path("x"),["s1"],spectrum_source=[spectrum(arrays=([100,136,150],[0,-2,4]))])[0]
    assert row.zero_intensity_peak_count==1 and row.negative_intensity_peak_count==1
    assert row.peaks==((150.0,4.0),)


def test_canonical_formula_products_and_match():
    products=generate_nucleoside_product_ion_hypotheses([candidate()])
    assert {x.product_ion_label for x in products}=={"BASE_MOLECULAR_ION","PROTONATED_BASE","DEHYDRATED_PRECURSOR"}
    matches=match_product_ions_to_ms2_peaks(products,[processed()],precursor_records=[precursor()])
    assert any(x.product_ion_label=="PROTONATED_BASE" for x in matches)


def test_shared_isobaric_product_is_nondiscriminating():
    c1,c2=candidate("X",300),candidate("Y",300)
    registry={x.candidate_id:[{"theoretical_product_mz":150,"candidate_specific":False}] for x in (c1,c2)}
    products=generate_nucleoside_product_ion_hypotheses([c1,c2],product_ion_registry=registry)
    pairs=compare_candidate_product_ions([c1,c2],products)
    assert pairs[0].shared_product_ion_count>=1 and not pairs[0].theoretical_discrimination_possible


def test_candidate_unique_product_allows_theoretical_discrimination():
    c1,c2=candidate("X",300),candidate("Y",300)
    products=generate_nucleoside_product_ion_hypotheses([c1,c2],product_ion_registry={"X":[{"theoretical_product_mz":151,"candidate_specific":True}]})
    assert compare_candidate_product_ions([c1,c2],products)[0].theoretical_discrimination_possible


def test_mass_only_candidate_gets_no_structure_specific_products():
    c=candidate("mass",310,MODEL_NOT_DEFINED,"MASS_ONLY_MODIFIED_NUCLEOSIDE")
    assert generate_nucleoside_product_ion_hypotheses([c])==[]


MODEL_NOT_DEFINED="MODEL_NOT_DEFINED"


def test_isomer_like_precursor_only_remains_unresolved():
    cs=[candidate("ac6A",310,MODEL_NOT_DEFINED,"MASS_ONLY_MODIFIED_NUCLEOSIDE","IDENTITY_AMBIGUOUS"),candidate("m6_6Am",310,MODEL_NOT_DEFINED,"MASS_ONLY_MODIFIED_NUCLEOSIDE","IDENTITY_AMBIGUOUS")]
    ps=[replace(precursor(c.candidate_id),candidate_name=c.candidate_id) for c in cs]
    out=summarize_candidate_ms2_identity_evidence(cs,ps,[],spectra=[processed()])
    assert all(x.ms2_identity_evidence_status=="MS2_PRECURSOR_COMPATIBLE_ONLY" and not x.exact_isomer_identity_confirmed for x in out)


def test_collision_energy_groups_are_not_merged():
    c=candidate(); ps=[precursor(energy=10),replace(precursor(spectrum_id="s2"),collision_energy=20)]
    recurrence, energies=_recurrence_and_energy([c],ps,[],[],[processed(),processed("s2")])
    assert len(energies)==2 and not recurrence


def test_missing_collision_energy_is_not_inferred():
    _, energies=_recurrence_and_energy([candidate()],[precursor(energy=None)],[],[],[processed()])
    assert energies[0].collision_energy_status=="NOT_RECORDED"


def test_product_recurrence_across_two_spectra():
    c=candidate(); products=generate_nucleoside_product_ion_hypotheses([c]); spectra=[processed(),processed("s2")]
    ps=[precursor(),precursor(spectrum_id="s2")]
    matches=match_product_ions_to_ms2_peaks(products,spectra,precursor_records=ps)
    recurrence,_=_recurrence_and_energy([c],ps,products,matches,spectra)
    base=next(x for x in recurrence if x.product_ion_id.endswith("PROTONATED_BASE"))
    assert base.supporting_ms2_spectrum_count==2 and base.product_ion_recurrence_status=="RECURRENT_HIGH"


def test_single_spectrum_not_overstated():
    c=candidate(); products=generate_nucleoside_product_ion_hypotheses([c]); matches=match_product_ions_to_ms2_peaks(products,[processed()],precursor_records=[precursor()])
    recurrence,_=_recurrence_and_energy([c],[precursor()],products,matches,[processed()])
    assert next(x for x in recurrence if x.supporting_ms2_spectrum_count).product_ion_recurrence_status=="SINGLE_SPECTRUM_ONLY"


def test_explained_intensity_uses_filtered_positive_peaks():
    c=candidate(); products=generate_nucleoside_product_ion_hypotheses([c]); s=processed(peaks=((136.06177,100),(200,100)))
    matches=match_product_ions_to_ms2_peaks(products,[s],precursor_records=[precursor()])
    summary=_spectrum_summaries([c],[precursor()],[s],products,matches)[0]
    assert summary.explained_intensity_fraction==pytest.approx(.5)


def test_product_ambiguity_for_one_observed_peak():
    c1,c2=candidate("X"),candidate("Y")
    registry={"X":[{"theoretical_product_mz":150}],"Y":[{"theoretical_product_mz":150}]}
    products=generate_nucleoside_product_ion_hypotheses([c1,c2],product_ion_registry=registry)
    ps=[precursor("X"),precursor("Y")]
    matches=match_product_ions_to_ms2_peaks(products,[processed(peaks=((150,10),))],precursor_records=ps)
    assert any(x.candidate_count_for_observed_peak==2 and "PRODUCT_ION_AMBIGUITY" in x.match_block_reasons for x in matches)


def test_no_product_match_is_insufficient_not_identity_support():
    c=candidate(); products=generate_nucleoside_product_ion_hypotheses([c]); s=processed(peaks=((200,10),))
    ss=_spectrum_summaries([c],[precursor()],[s],products,[])
    out=summarize_candidate_ms2_identity_evidence([c],[precursor()],[],product_ions=products,spectrum_summaries=ss,spectra=[s])
    assert out[0].ms2_identity_evidence_status=="MS2_INSUFFICIENT"


def test_ms1_state_reconciliation_never_resolves_plus16_chemistry():
    class Family: state_family_id="F"; base_candidate_id="A"
    class Result: state_families=(Family(),)
    summaries=summarize_candidate_ms2_identity_evidence([candidate()],[precursor()],[],spectra=[processed()])
    row=reconcile_p1ap_ms1_and_ms2_evidence(Result(),summaries,t1_result=object(),full_length_series=(0,18,34,50))[0]
    assert not row.state_interpretation_resolved and not row.reaction_order_assigned


def test_safeguards_are_false_for_localization_and_formal_propagation():
    c=candidate()
    assert not c.exact_nucleotide_localization and not c.exact_atom_localization
    assert not c.applied_to_formal_score and not c.applied_to_final_consensus


def test_deterministic_selection_product_and_match_order():
    cs=[candidate("Y",300),candidate("X",300)]
    specs=[spectrum("s2",selected=300,target=300),spectrum("s1",selected=300,target=300)]
    left=select_precursor_compatible_ms2_spectra(Path("x"),cs,metadata_spectra=specs)
    right=select_precursor_compatible_ms2_spectra(Path("x"),list(reversed(cs)),metadata_spectra=list(reversed(specs)))
    assert left==right
    assert generate_nucleoside_product_ion_hypotheses(cs)==generate_nucleoside_product_ion_hypotheses(list(reversed(cs)))


def test_optional_result_only_adds_shadow_records_and_strips_peak_arrays():
    c=candidate(); summary=P1APNucleosideMS2Summary(source_id="S",status="COMPLETED",total_spectra_seen=1,total_ms2_metadata_records=1,precursor_compatible_record_count=0,decoded_unique_spectrum_count=0,candidate_count=1,theoretical_product_ion_count=0,product_match_count=0,overall_evidence_status="MS2_INSUFFICIENT",overall_confidence="LOW",overall_block_reasons=(),runtime_seconds=0)
    result=P1APNucleosideMS2AuditResult((c,),(),(processed(),),(),(),(),(),(),(),(),(),summary)
    payload=audit_optional_result(result)
    assert "peaks" not in payload["spectrum_records"][0]
    assert payload["summary_records"][0]["formal_propagation"] is False
    assert payload["summary_records"][0]["applied_to_ranking"] is False


def test_missing_ms1_result_is_deterministically_blocked(tmp_path):
    from rna_masshunter.sciex_p1ap_nucleoside_ms2_identity_audit import audit_p1ap_nucleoside_ms2_identity
    result=audit_p1ap_nucleoside_ms2_identity(tmp_path/"missing.mzML",p1ap_ms1_result=None)
    assert result.summary.status=="BLOCKED"
    assert result.summary.overall_block_reasons==("INPUT_FILE_NOT_FOUND","SOURCE_METADATA_RECORD_MISSING","P1AP_MS1_AUDIT_RESULT_MISSING")
