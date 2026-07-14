from pathlib import Path
import yaml
from rna_masshunter.audit_policy import AUDIT_DETAIL, AUDIT_SUMMARY, AuditPolicy, sheet_category
from rna_masshunter.modification_composer import apply_transform_ids
from rna_masshunter.modification_constraints import load_transformations
from rna_masshunter.composite_modification_audit import build_composite_modification_audit, append_composite_diagnostics
from rna_masshunter.masses import load_base_masses
from rna_masshunter.modifications import load_modifications
ROOT=Path(__file__).parent
def build(): return build_composite_modification_audit(ROOT,'AGCU',load_modifications(ROOT/'data/modifications.yaml'),load_base_masses(ROOT/'data/base_masses.yaml'))
def test_five_shadow_sheets_and_excel_names():
    r=build(); assert set(r.sheets)=={'Composite_Mod_Candidates','Composite_Mod_Invalid','Composite_Mod_Summary','Backbone_Mod_Candidates','Cleavage_Block_Audit'} and all(len(x)<=31 for x in r.sheets)
def test_registry_and_modes():
    assert sheet_category('Composite_Mod_Summary')==AUDIT_SUMMARY and sheet_category('Cleavage_Block_Audit')==AUDIT_SUMMARY
    assert sheet_category('Composite_Mod_Candidates')==AUDIT_DETAIL
    assert not AuditPolicy.from_level('standard').run_shadow_audits and not AuditPolicy.from_level('audit').include_detail and AuditPolicy.from_level('full').include_detail
def test_summary_nonpropagation_flags():
    row=build().sheets['Composite_Mod_Summary'][0]; assert row['Formal_Result_Changed'] is False and row['Applied_To_Formal_Result'] is False and row['Formal_Change_Ready'] is False
def test_all_rows_not_applied():
    r=build(); assert all(row['Applied_To_Formal_Result'] is False for rows in r.sheets.values() for row in rows if 'Applied_To_Formal_Result' in row)
def test_diagnostics_standard_not_run_not_zero():
    row=append_composite_diagnostics([{}],None)[0]; assert row['Composite_Mod_Valid_Count']=='not_run' and row['Composite_Mod_Applied_To_Formal_Result'] is False
def test_deterministic_rerun():
    a=build(); b=build(); assert a.sheets==b.sheets
def test_legacy_overlap_and_all_bases():
    r=build(); summary=r.sheets['Composite_Mod_Summary'][0]; assert summary['Legacy_Mapped_Transform_Count']>0 and all(summary[f'{base}_Candidate_Count']>0 for base in 'AGCU')
def test_example_fixture_is_hypothesis_and_composable():
    example=yaml.safe_load((ROOT/'data/examples/composite_structure_example.yaml').read_text()); assert example['example_only'] is True
    wanted=example['nucleoside_states'][0]['transforms']; transforms=load_transformations(ROOT/'data/modification_transforms_v2.yaml')
    state,result,order=apply_transform_ids('U',37,wanted,transforms,ROOT/'data/nucleoside_slots.yaml')
    assert result.valid and set(order)==set(wanted) and state.slot_state_dict['U_O2']=='sulfur' and state.slot_state_dict['U_C5_side_chain_carbonyl']=='oxidized_sulfur_1'
def test_backbone_and_cleavage_counts():
    r=build(); assert len(r.sheets['Backbone_Mod_Candidates'])==3 and len(r.sheets['Cleavage_Block_Audit'])==12
def test_formal_objects_are_not_arguments_or_outputs():
    r=build(); assert all('Final_Score' not in row and 'Rank' not in row for rows in r.sheets.values() for row in rows)
