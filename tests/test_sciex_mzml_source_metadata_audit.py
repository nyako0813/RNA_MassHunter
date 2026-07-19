from dataclasses import replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import os
import pytest
from rna_masshunter.sciex_mzml_source_metadata_audit import *

NS='http://psi.hupo.org/ms/mzml'
def spectrum(index,ms=1,pol='negative',rep='profile',time=.1,*,both_pol=False,both_rep=False,include=True,mz=True):
 cv=[]
 if include:cv.append(f'<cvParam accession="MS:1000511" name="ms level" value="{ms}"/>')
 if pol:cv.append(f'<cvParam accession="MS:1000129" name="negative scan" value=""/>' if pol=='negative' else '<cvParam accession="MS:1000130" name="positive scan" value=""/>')
 if both_pol:cv.append('<cvParam accession="MS:1000130" name="positive scan" value=""/>')
 if rep:cv.append('<cvParam accession="MS:1000128" name="profile spectrum" value=""/>' if rep=='profile' else '<cvParam accession="MS:1000127" name="centroid spectrum" value=""/>')
 if both_rep:cv.append('<cvParam accession="MS:1000127" name="centroid spectrum" value=""/>')
 if time is not None:cv.append(f'<scanList><scan><cvParam accession="MS:1000016" name="scan start time" value="{time}" unitAccession="UO:0000031" unitName="minute"/>{"<scanWindowList><scanWindow><cvParam accession=\"MS:1000501\" name=\"scan window lower limit\" value=\"100\"/><cvParam accession=\"MS:1000500\" name=\"scan window upper limit\" value=\"2000\"/></scanWindow></scanWindowList>" if mz else ""}</scan></scanList>')
 return f'<spectrum index="{index}" id="scan={index}">{"".join(cv)}<binaryDataArrayList><binaryDataArray><binary>QUJDREVGRw==</binary></binaryDataArray></binaryDataArrayList></spectrum>'
def document(spectra,*,run='_x0030_1',start='2026-01-01T00:00:00Z',source=True,declared=None,software='3.0',processing=True,instrument=True,indexed=True):
 count=len(spectra) if declared is None else declared;src='<sourceFileList count="1"><sourceFile id="SF1" name="01.wiff2" location="file://C:\\Data\\Run"><cvParam accession="MS:1000562" name="ABI WIFF format" value=""/></sourceFile></sourceFileList>' if source else ''
 soft=f'<softwareList count="2"><software id="Analyst" version="unknown"><cvParam accession="MS:1000551" name="Analyst" value=""/></software><software id="pwiz_Reader_ABI" version="{software}"><cvParam accession="MS:1000615" name="ProteoWizard software" value=""/></software></softwareList>'
 inst='<instrumentConfigurationList><instrumentConfiguration id="IC1"><cvParam accession="MS:1003293" name="ZenoTOF 7600" value=""/><componentList><source order="1"><cvParam accession="MS:1000073" name="electrospray ionization" value=""/></source><analyzer order="2"><cvParam accession="MS:1000084" name="time-of-flight" value=""/></analyzer><detector order="3"><cvParam accession="MS:1000253" name="electron multiplier" value=""/></detector></componentList></instrumentConfiguration></instrumentConfigurationList>' if instrument else ''
 dp='<dataProcessingList><dataProcessing id="DP1"><processingMethod order="0" softwareRef="pwiz_Reader_ABI"><cvParam accession="MS:1000544" name="Conversion to mzML" value=""/><userParam name="msconvert filter" value="peakPicking true 1-"/></processingMethod></dataProcessing></dataProcessingList>' if processing else ''
 attrs=f'id="{run}"' if run is not None else '';attrs+=f' startTimeStamp="{start}"' if start is not None else '';attrs+=' defaultSourceFileRef="SF1" defaultInstrumentConfigurationRef="IC1"'
 mz=f'<mzML xmlns="{NS}" id="DOC" version="1.1.0"><fileDescription>{src}</fileDescription>{soft}{inst}{dp}<run {attrs}><spectrumList count="{count}" defaultDataProcessingRef="DP1">{"".join(spectra)}</spectrumList></run></mzML>'
 return f'<?xml version="1.0"?><indexedmzML xmlns="{NS}">{mz}</indexedmzML>' if indexed else '<?xml version="1.0"?>'+mz
def write(tmp_path,name='a.mzML',**kwargs):
 p=tmp_path/name;p.write_text(document(**kwargs),encoding='utf-8');return p

def test_basic_metadata_parse(tmp_path):
 p=write(tmp_path,spectra=[spectrum(0,1,time=.1),spectrum(1,2,time=.2)]);r=parse_sciex_mzml_source_metadata(p);assert r.run_id_decoded=='01';assert r.run_start_time_normalized=='2026-01-01T00:00:00Z';assert r.source_file_names==('01.wiff2',);assert r.declared_spectrum_count==r.parsed_spectrum_count==2;assert (r.ms1_spectrum_count,r.ms2_spectrum_count)==(1,1);assert r.negative_spectrum_count==2;assert r.profile_spectrum_count==2;assert r.first_scan_start_time==.1 and r.last_scan_start_time==.2;assert r.acquisition_software==('Analyst',);assert r.conversion_software==('pwiz_Reader_ABI',);assert r.processing_history_status is MetadataStatus.RECORDED;assert r.default_instrument_configuration_ref=='IC1'
def test_mixed_ms_level_set(tmp_path):
 r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0,1),spectrum(1,2)]));assert r.ms_level_set==(1,2) and r.ms_level_counts==((1,1),(2,1))
def test_negative_profile_file(tmp_path):
 r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0),spectrum(1)]));assert r.polarity_status is PolarityStatus.NEGATIVE_ONLY and r.representation_status is RepresentationStatus.PROFILE_ONLY
@pytest.mark.parametrize('missing,block',[('polarity','MISSING_POLARITY_METADATA'),('representation','MISSING_SPECTRUM_REPRESENTATION'),('ms','MISSING_MS_LEVEL_METADATA'),('time','MISSING_SCAN_TIME')])
def test_missing_spectrum_metadata(tmp_path,missing,block):
 kwargs={'pol':'negative','rep':'profile','include':True,'time':.1};kwargs[{'polarity':'pol','representation':'rep','ms':'include','time':'time'}[missing]]=None if missing!='ms' else False;r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0,**kwargs)]));assert block in r.block_reasons
def test_missing_document_metadata(tmp_path):
 p=write(tmp_path,spectra=[spectrum(0)],run=None,start=None,source=False,processing=False,instrument=False);r=parse_sciex_mzml_source_metadata(p);assert {'MISSING_RUN_ID','MISSING_RUN_START_TIMESTAMP','MISSING_SOURCE_FILE','MISSING_INSTRUMENT_METADATA','MISSING_PROCESSING_METADATA'}<=set(r.block_reasons)
def test_conflicting_within_spectrum(tmp_path):
 r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0,both_pol=True,both_rep=True)]));assert r.polarity_status is PolarityStatus.CONFLICTING_WITHIN_SPECTRUM;assert r.representation_status is RepresentationStatus.CONFLICTING_WITHIN_SPECTRUM
def test_declared_parsed_mismatch(tmp_path):
 r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0),spectrum(1)],declared=3));assert r.spectrum_count_consistent is False and 'DECLARED_PARSED_SPECTRUM_COUNT_MISMATCH' in r.block_reasons
@pytest.mark.parametrize('raw,decoded,status',[('_x0030_1','01',RunIDDecodeStatus.DECODED),('_x0020_',' ',RunIDDecodeStatus.DECODED),('_x0028_2_x0029_','(2)',RunIDDecodeStatus.DECODED),('normal','normal',RunIDDecodeStatus.UNCHANGED),('_x003_','_x003_',RunIDDecodeStatus.PARTIALLY_INVALID_ESCAPE_PRESERVED),('_x00a9_','©',RunIDDecodeStatus.DECODED),('_x0030_1-01_LeuUAA_T1_2','01-01_LeuUAA_T1_2',RunIDDecodeStatus.DECODED),('_x0030_1-01_LeuUAA_T1_2_x0020__x0028_2_x0029_','01-01_LeuUAA_T1_2 (2)',RunIDDecodeStatus.DECODED)])
def test_id_decode(raw,decoded,status):assert decode_mzml_xml_safe_id(raw)==(decoded,status)
@pytest.mark.parametrize('raw,expected',[('file://C:\\Data\\Run\\','file://c:/Data/Run'),('file:///C:/Data%20Set','file://c:/Data Set'),(None,None)])
def test_location_normalization(raw,expected):assert normalize_source_location_for_comparison(raw)==expected
def test_exact_duplicate_different_paths(tmp_path):
 p=write(tmp_path,'a.mzML',spectra=[spectrum(0)]);q=tmp_path/'b.mzML';q.write_bytes(p.read_bytes());rel=compare_mzml_source_metadata(parse_sciex_mzml_source_metadata(p),parse_sciex_mzml_source_metadata(q));assert rel.duplicate_file_status is RelationshipStatus.EXACT_DUPLICATE and rel.exact_duplicate and not rel.keep_both_files
def test_same_file_object(tmp_path):
 p=write(tmp_path,'a.mzML',spectra=[spectrum(0)]);q=tmp_path/'b.mzML';os.link(p,q);rel=compare_mzml_source_metadata(parse_sciex_mzml_source_metadata(p),parse_sciex_mzml_source_metadata(q));assert rel.duplicate_file_status is RelationshipStatus.SAME_FILE_OBJECT and rel.same_inode
def test_different_run_same_source(tmp_path):
 a=parse_sciex_mzml_source_metadata(write(tmp_path,'a.mzML',spectra=[spectrum(0)],run='R1',start='2026-01-01T00:00:00Z'));b=parse_sciex_mzml_source_metadata(write(tmp_path,'b.mzML',spectra=[spectrum(0),spectrum(1,2)],run='R2',start='2026-01-01T00:01:00Z'));rel=compare_mzml_source_metadata(a,b);assert rel.duplicate_file_status is RelationshipStatus.DIFFERENT_RUN_SAME_SOURCE_WIFF and rel.keep_both_files and not rel.same_run
def test_same_run_different_export(tmp_path):
 a=parse_sciex_mzml_source_metadata(write(tmp_path,'a.mzML',spectra=[spectrum(0)],run='R',software='1'));b=parse_sciex_mzml_source_metadata(write(tmp_path,'b.mzML',spectra=[spectrum(0)],run='R',software='2'));assert compare_mzml_source_metadata(a,b).duplicate_file_status is RelationshipStatus.SAME_RUN_DIFFERENT_EXPORT
def test_possible_reexport_missing_processing(tmp_path):
 a=parse_sciex_mzml_source_metadata(write(tmp_path,'a.mzML',spectra=[spectrum(0)],run='R',processing=False));b=parse_sciex_mzml_source_metadata(write(tmp_path,'b.mzML',spectra=[spectrum(0)],run='R',software='2',processing=False));assert compare_mzml_source_metadata(a,b).duplicate_file_status is RelationshipStatus.POSSIBLE_REEXPORT_OF_SAME_RUN
def test_different_source(tmp_path):
 a=parse_sciex_mzml_source_metadata(write(tmp_path,'a.mzML',spectra=[spectrum(0)]));b=replace(parse_sciex_mzml_source_metadata(write(tmp_path,'b.mzML',spectra=[spectrum(0)])),sha256='different',source_files=(SourceFileMetadata('X','other.wiff','file://x','file://x',()),));assert compare_mzml_source_metadata(a,b).duplicate_file_status is RelationshipStatus.DIFFERENT_SOURCE_FILE
def test_invalid_and_unreadable_relationship(tmp_path):
 p=tmp_path/'bad.mzML';p.write_text('<bad>',encoding='utf-8');a=parse_sciex_mzml_source_metadata(p);b=parse_sciex_mzml_source_metadata(tmp_path/'missing.mzML');assert compare_mzml_source_metadata(a,b).duplicate_file_status is RelationshipStatus.UNREADABLE_FILE
def test_invalid_root(tmp_path):
 p=tmp_path/'bad.mzML';p.write_text('<root/>',encoding='utf-8');assert parse_sciex_mzml_source_metadata(p).read_status is ReadStatus.INVALID_MZML
def test_xml_parse_error(tmp_path):
 p=tmp_path/'bad.mzML';p.write_text('<mzML>',encoding='utf-8');assert parse_sciex_mzml_source_metadata(p).read_status is ReadStatus.XML_PARSE_ERROR
def test_deterministic_audit_order(tmp_path):
 paths=[write(tmp_path,'b.mzML',spectra=[spectrum(0)],run='B'),write(tmp_path,'a.mzML',spectra=[spectrum(0)],run='A')];x=audit_mzml_source_metadata_files(paths);y=audit_mzml_source_metadata_files(reversed(paths));assert x==y and x.relationship_records[0].left_file=='a.mzML'
def test_all_eight_pair_count_synthetic(tmp_path):
 paths=[write(tmp_path,f'{i}.mzML',spectra=[spectrum(0)],run=str(i)) for i in range(8)];assert len(audit_mzml_source_metadata_files(paths).relationship_records)==28
def test_table_rows_are_serializable_scalars(tmp_path):
 r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0)]));row=file_metadata_row(r);assert all(not isinstance(v,(Path,list,dict,tuple)) for v in row.values())
def test_binary_not_decoded_required_metadata_works(tmp_path):
 r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0)]));assert r.parsed_spectrum_count==1
def test_mz_range_metadata(tmp_path):
 r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0,mz=True)]));assert r.mz_range_status is MZRangeStatus.PARTIALLY_RECORDED and r.scan_window_lower_limit_min==100
def test_instrument_components(tmp_path):
 r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0)]));assert r.instrument_model_name==('ZenoTOF 7600',) and {x.component_type for x in r.instrument_components}=={'source','analyzer','detector'}
def test_processing_filter_history(tmp_path):
 r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0)]));assert r.msconvert_filter_history==('peakPicking true 1-',)
@pytest.mark.parametrize('field',['formal_propagation','chemical_identity_assigned','rna_identity_confirmed','applied_to_formal_score','applied_to_ranking','applied_to_candidate_filtering','applied_to_final_consensus'])
def test_formal_nonpropagation(tmp_path,field):
 r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0)]));assert getattr(r,field) is False
def test_optional_result_key_shapes(tmp_path):
 result=audit_mzml_source_metadata_files([write(tmp_path,spectra=[spectrum(0)])]);payload=audit_optional_result(result);assert set(payload)=={'file_records','relationship_records','summary'}



def test_normalized_child_metadata_properties(tmp_path):
 r=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0)]));assert r.source_file_count==1 and r.source_file_ids==('SF1',) and r.source_file_accession_metadata;assert r.software_ids==('Analyst','pwiz_Reader_ABI') and r.software_versions==('unknown','3.0');assert r.processing_method_count==1 and r.processing_method_software_refs==('pwiz_Reader_ABI',) and r.processing_cv_params

def test_user_manifest_metadata_conflict_is_retained(tmp_path):
 parsed=parse_sciex_mzml_source_metadata(write(tmp_path,spectra=[spectrum(0)]));parsed=replace(parsed,source_confirmed_context='INTERNAL_OTHER_CONTEXT');context=RuntimeSourceContext('TRNA_LEU_UAA','RNASE_T1_DIGEST','RUN_1');result=apply_runtime_source_context(parsed,context);assert result.rna_identity=='TRNA_LEU_UAA' and result.source_confirmed_context=='INTERNAL_OTHER_CONTEXT';assert result.context_conflict and 'USER_MANIFEST_METADATA_CONFLICT' in result.block_reasons

def test_user_runtime_context_is_independent_provenance(tmp_path):
 p=write(tmp_path,spectra=[spectrum(0)]);ctx=RuntimeSourceContext('TRNA_LEU_UAA','RNASE_T1_DIGEST','UAA_T1_RUN_1');r=audit_mzml_source_metadata_files([p],runtime_contexts={p.name:ctx}).file_records[0];assert r.rna_identity=='TRNA_LEU_UAA' and r.digest_type=='RNASE_T1_DIGEST' and r.technical_run_label=='UAA_T1_RUN_1';assert r.context_source=='USER_PROVIDED_RUNTIME_MANIFEST' and r.context_confidence=='USER_CONFIRMED' and not r.mzml_metadata_confirmed and r.rna_identity_confirmed;assert r.source_confirmed_context=='RNA_Identity=TRNA_LEU_UAA;Digest_Type=RNASE_T1_DIGEST;Technical_Run_Label=UAA_T1_RUN_1' and r.filename_label_only=='NOT_APPLICABLE_USER_MANIFEST_SUPPLIED'

def test_runtime_context_input_mapping_is_generic(tmp_path):
 p=write(tmp_path,'arbitrary.mzML',spectra=[spectrum(0)]);ctx=RuntimeSourceContext('ANY_RNA','ANY_DIGEST','RUN_X');assert audit_mzml_source_metadata_files([p],runtime_contexts={str(p):ctx}).file_records[0].technical_run_label=='RUN_X'

def test_context_does_not_change_relationship(tmp_path):
 a=write(tmp_path,'a.mzML',spectra=[spectrum(0)],run='A');b=write(tmp_path,'b.mzML',spectra=[spectrum(0),spectrum(1)],run='B');ctx={a.name:RuntimeSourceContext('R','T1','R1'),b.name:RuntimeSourceContext('R','T1','R2')};assert audit_mzml_source_metadata_files([a,b],runtime_contexts=ctx).relationship_records[0].duplicate_file_status is RelationshipStatus.DIFFERENT_RUN_SAME_SOURCE_WIFF
