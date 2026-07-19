"""Streaming, binary-array-free metadata reconciliation for SCIEX-derived mzML."""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass,replace
from datetime import datetime,timezone
from enum import Enum
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import re
from typing import Iterable,Mapping,Sequence
from urllib.parse import unquote
import xml.etree.ElementTree as ET

class ReadStatus(str,Enum):COMPLETED="COMPLETED";FILE_NOT_FOUND="FILE_NOT_FOUND";FILE_NOT_REGULAR="FILE_NOT_REGULAR";FILE_UNREADABLE="FILE_UNREADABLE";XML_PARSE_ERROR="XML_PARSE_ERROR";INVALID_MZML="INVALID_MZML"
class RunIDDecodeStatus(str,Enum):DECODED="DECODED";UNCHANGED="UNCHANGED";PARTIALLY_INVALID_ESCAPE_PRESERVED="PARTIALLY_INVALID_ESCAPE_PRESERVED"
class PolarityStatus(str,Enum):NEGATIVE_ONLY="NEGATIVE_ONLY";POSITIVE_ONLY="POSITIVE_ONLY";MIXED_POLARITY="MIXED_POLARITY";NOT_RECORDED="NOT_RECORDED";CONFLICTING_WITHIN_SPECTRUM="CONFLICTING_WITHIN_SPECTRUM"
class RepresentationStatus(str,Enum):PROFILE_ONLY="PROFILE_ONLY";CENTROID_ONLY="CENTROID_ONLY";MIXED_REPRESENTATION="MIXED_REPRESENTATION";NOT_RECORDED="NOT_RECORDED";CONFLICTING_WITHIN_SPECTRUM="CONFLICTING_WITHIN_SPECTRUM"
class MZRangeStatus(str,Enum):RECORDED="RECORDED";PARTIALLY_RECORDED="PARTIALLY_RECORDED";NOT_RECORDED="NOT_RECORDED"
class MetadataStatus(str,Enum):RECORDED="RECORDED";PARTIALLY_RECORDED="PARTIALLY_RECORDED";NOT_RECORDED="NOT_RECORDED"
SCIEX_MZML_SOURCE_METADATA_AUDIT_OPTIONAL_RESULT_KEY="sciex_mzml_source_metadata_audit"
class RelationshipStatus(str,Enum):EXACT_DUPLICATE="EXACT_DUPLICATE";SAME_FILE_OBJECT="SAME_FILE_OBJECT";DIFFERENT_RUN_SAME_SOURCE_WIFF="DIFFERENT_RUN_SAME_SOURCE_WIFF";POSSIBLE_DIFFERENT_RUN_SAME_SOURCE="POSSIBLE_DIFFERENT_RUN_SAME_SOURCE";SAME_RUN_DIFFERENT_EXPORT="SAME_RUN_DIFFERENT_EXPORT";POSSIBLE_REEXPORT_OF_SAME_RUN="POSSIBLE_REEXPORT_OF_SAME_RUN";DIFFERENT_SOURCE_FILE="DIFFERENT_SOURCE_FILE";INSUFFICIENT_SOURCE_METADATA="INSUFFICIENT_SOURCE_METADATA";UNREADABLE_FILE="UNREADABLE_FILE";INVALID_MZML="INVALID_MZML"

_BLOCK_ORDER=("FILE_NOT_FOUND","FILE_NOT_REGULAR","FILE_UNREADABLE","XML_PARSE_ERROR","INVALID_MZML_ROOT","MISSING_RUN_ID","MISSING_RUN_START_TIMESTAMP","MISSING_SOURCE_FILE","MULTIPLE_SOURCE_FILES","DEFAULT_SOURCE_REF_UNRESOLVED","MISSING_SPECTRUM_LIST","DECLARED_PARSED_SPECTRUM_COUNT_MISMATCH","MISSING_MS_LEVEL_METADATA","MIXED_MS_LEVEL_METADATA","MISSING_POLARITY_METADATA","MIXED_POLARITY_METADATA","CONFLICTING_POLARITY_WITHIN_SPECTRUM","MISSING_SPECTRUM_REPRESENTATION","MIXED_SPECTRUM_REPRESENTATION","CONFLICTING_REPRESENTATION_WITHIN_SPECTRUM","MISSING_SCAN_TIME","MIXED_SCAN_TIME_UNITS","MISSING_MZ_RANGE_METADATA","MISSING_INSTRUMENT_METADATA","MISSING_PROCESSING_METADATA","SOURCE_FILE_REFERENCE_UNRESOLVED","INSUFFICIENT_RUN_IDENTITY","INSUFFICIENT_SOURCE_IDENTITY","INSUFFICIENT_TXT_MZML_LINKAGE","USER_MANIFEST_METADATA_CONFLICT")
def _ordered_blocks(values:Iterable[str]):
 found=set(values);return tuple(x for x in _BLOCK_ORDER if x in found)+tuple(sorted(found-set(_BLOCK_ORDER)))
def _local(tag:str):return tag.rsplit("}",1)[-1]
def _cv(elem):return tuple((str(x.get("accession","")),str(x.get("name","")),str(x.get("value","")),str(x.get("unitAccession","")),str(x.get("unitName",""))) for x in elem.iter() if _local(x.tag)=="cvParam")
def _users(elem):return tuple((str(x.get("name","")),str(x.get("value",""))) for x in elem.iter() if _local(x.tag)=="userParam")
def _float(value):
 try:return float(value)
 except (TypeError,ValueError):return None

def decode_mzml_xml_safe_id(value:str):
 changed=False
 def sub(match):
  nonlocal changed
  try:result=chr(int(match.group(1),16));changed=True;return result
  except (ValueError,OverflowError):return match.group(0)
 decoded=re.sub(r"_x([0-9A-Fa-f]{4})_",sub,str(value));invalid=bool(re.search(r"_x[^_]*$|_x[0-9A-Fa-f]{0,3}_",decoded));status=RunIDDecodeStatus.PARTIALLY_INVALID_ESCAPE_PRESERVED if invalid else RunIDDecodeStatus.DECODED if changed else RunIDDecodeStatus.UNCHANGED
 return decoded,status
def calculate_sha256(path:Path,chunk_size:int=1024*1024):
 digest=sha256()
 with Path(path).open("rb") as handle:
  for block in iter(lambda:handle.read(chunk_size),b""):digest.update(block)
 return digest.hexdigest()
def normalize_source_location_for_comparison(value:str|None):
 if not value:return None
 text=unquote(str(value).strip()).replace("\\","/");text=re.sub(r"^file:(?:/{0,3})", "file://",text,flags=re.I);text=re.sub(r"/+","/",text.replace("file://","file§§",1)).replace("file§§","file://",1);text=text.rstrip("/")
 m=re.match(r"^(file://)([A-Za-z]):",text)
 if m:text=m.group(1)+m.group(2).lower()+text[m.end(2):]
 return text

def _normalize_timestamp(value:str|None):
 if not value:return "NOT_RECORDED"
 try:
  dt=datetime.fromisoformat(value.replace("Z","+00:00"));return dt.astimezone(timezone.utc).isoformat().replace("+00:00","Z")
 except ValueError:return str(value)
@dataclass(frozen=True)
class RuntimeSourceContext:
 rna_identity:str;digest_type:str;technical_run_label:str;context_source:str="USER_PROVIDED_RUNTIME_MANIFEST";context_confidence:str="USER_CONFIRMED";mzml_metadata_confirmed:bool=False
@dataclass(frozen=True)
class SourceFileMetadata:
 source_file_id:str;name:str;location_raw:str;location_normalized_for_comparison:str|None;accession_metadata:tuple[tuple[str,str,str],...]
@dataclass(frozen=True)
class SoftwareMetadata:
 software_id:str;name:str;version:str;cv_params:tuple[tuple[str,str,str],...]
@dataclass(frozen=True)
class ProcessingMethodMetadata:
 data_processing_id:str;order:str;software_ref:str;cv_params:tuple[tuple[str,str,str],...];user_params:tuple[tuple[str,str],...]
@dataclass(frozen=True)
class InstrumentComponentMetadata:
 instrument_configuration_id:str;component_type:str;order:str;cv_params:tuple[tuple[str,str,str],...]
@dataclass(frozen=True,kw_only=True)
class MzMLSourceMetadataRecord:
 input_path:str;file_name:str;file_size_bytes:int|str="NOT_AVAILABLE";sha256:str="NOT_AVAILABLE";file_exists:bool=False;is_regular_file:bool=False;is_symlink:bool=False;device_id:int|str="NOT_AVAILABLE";inode:int|str="NOT_AVAILABLE";read_status:ReadStatus=ReadStatus.FILE_NOT_FOUND
 mzml_version:str="NOT_RECORDED";mzml_id:str="NOT_RECORDED";run_id_raw:str="NOT_RECORDED";run_id_decoded:str="NOT_RECORDED";run_id_decode_status:RunIDDecodeStatus=RunIDDecodeStatus.UNCHANGED;run_start_time_raw:str="NOT_RECORDED";run_start_time_normalized:str="NOT_RECORDED";default_source_file_ref:str="NOT_RECORDED";default_instrument_configuration_ref:str="NOT_RECORDED";default_data_processing_ref:str="NOT_RECORDED";declared_spectrum_count:int|str="NOT_RECORDED";parsed_spectrum_count:int=0;spectrum_count_consistent:bool|str="NOT_AVAILABLE"
 source_files:tuple[SourceFileMetadata,...]=();default_source_file_resolved:bool=False;ms_level_set:tuple[int,...]=();ms_level_counts:tuple[tuple[int,int],...]=();ms1_spectrum_count:int=0;ms2_spectrum_count:int=0;other_ms_level_count:int=0;missing_ms_level_count:int=0;positive_spectrum_count:int=0;negative_spectrum_count:int=0;mixed_polarity_spectrum_count:int=0;missing_polarity_count:int=0;polarity_status:PolarityStatus=PolarityStatus.NOT_RECORDED;profile_spectrum_count:int=0;centroid_spectrum_count:int=0;mixed_representation_spectrum_count:int=0;missing_representation_count:int=0;representation_status:RepresentationStatus=RepresentationStatus.NOT_RECORDED
 first_scan_start_time:float|str="NOT_RECORDED";last_scan_start_time:float|str="NOT_RECORDED";minimum_scan_start_time:float|str="NOT_RECORDED";maximum_scan_start_time:float|str="NOT_RECORDED";scan_time_unit_set:tuple[str,...]=();approximate_run_duration:float|str="NOT_RECORDED"
 lowest_observed_mz_min:float|str="NOT_RECORDED";lowest_observed_mz_max:float|str="NOT_RECORDED";highest_observed_mz_min:float|str="NOT_RECORDED";highest_observed_mz_max:float|str="NOT_RECORDED";scan_window_lower_limit_min:float|str="NOT_RECORDED";scan_window_lower_limit_max:float|str="NOT_RECORDED";scan_window_upper_limit_min:float|str="NOT_RECORDED";scan_window_upper_limit_max:float|str="NOT_RECORDED";mz_range_source:tuple[str,...]=();mz_range_status:MZRangeStatus=MZRangeStatus.NOT_RECORDED
 software:tuple[SoftwareMetadata,...]=();acquisition_software:tuple[str,...]=();conversion_software:tuple[str,...]=();conversion_software_version:tuple[str,...]=();data_processing_ids:tuple[str,...]=();processing_methods:tuple[ProcessingMethodMetadata,...]=();msconvert_filter_history:tuple[str,...]=();processing_history_status:MetadataStatus=MetadataStatus.NOT_RECORDED
 instrument_configuration_ids:tuple[str,...]=();instrument_model_cv:tuple[str,...]=();instrument_model_name:tuple[str,...]=();ion_source_cv:tuple[str,...]=();analyzer_cv:tuple[str,...]=();detector_cv:tuple[str,...]=();instrument_components:tuple[InstrumentComponentMetadata,...]=();default_instrument_configuration_resolved:bool=False;instrument_metadata_status:MetadataStatus=MetadataStatus.NOT_RECORDED
 metadata_confidence:str="LOW";block_reasons:tuple[str,...]=();filename_label_only:str="UNKNOWN";source_confirmed_context:str="UNKNOWN";rna_identity:str="UNKNOWN";digest_type:str="UNKNOWN";technical_run_label:str="UNKNOWN";context_source:str="NOT_AVAILABLE";context_confidence:str="UNVERIFIED";mzml_metadata_confirmed:bool=False;context_conflict:bool=False;formal_propagation:bool=False;shadow_analysis_only:bool=True;source_metadata_audit_only:bool=True;chemical_identity_assigned:bool=False;rna_identity_confirmed:bool=False;applied_to_formal_score:bool=False;applied_to_ranking:bool=False;applied_to_candidate_filtering:bool=False;applied_to_final_consensus:bool=False
 @property
 def source_file_count(self):return len(self.source_files)
 @property
 def source_file_ids(self):return tuple(x.source_file_id for x in self.source_files)
 @property
 def source_file_names(self):return tuple(x.name for x in self.source_files)
 @property
 def source_file_locations(self):return tuple(x.location_raw for x in self.source_files)
 @property
 def source_file_accession_metadata(self):return tuple(x.accession_metadata for x in self.source_files)
 @property
 def software_ids(self):return tuple(x.software_id for x in self.software)
 @property
 def software_names(self):return tuple(x.name for x in self.software)
 @property
 def software_versions(self):return tuple(x.version for x in self.software)
 @property
 def processing_method_count(self):return len(self.processing_methods)
 @property
 def processing_method_orders(self):return tuple(x.order for x in self.processing_methods)
 @property
 def processing_method_software_refs(self):return tuple(x.software_ref for x in self.processing_methods)
 @property
 def processing_cv_params(self):return tuple(x.cv_params for x in self.processing_methods)
 @property
 def processing_user_params(self):return tuple(x.user_params for x in self.processing_methods)
@dataclass(frozen=True,kw_only=True)
class MzMLFileRelationship:
 left_file:str;right_file:str;exact_duplicate:bool;same_file_size:bool;same_sha256:bool;same_inode:bool;same_source_file:bool;same_source_location:bool;same_default_source_ref:bool;same_run_id_raw:bool;same_run_id_decoded:bool;same_run_start_time:bool;same_declared_spectrum_count:bool;same_parsed_spectrum_count:bool;same_ms_level_distribution:bool;same_polarity_distribution:bool;same_representation_distribution:bool;same_conversion_software:bool;same_instrument_configuration:bool;duplicate_file_status:RelationshipStatus;keep_both_files:bool;relationship_confidence:str;relationship_block_reasons:tuple[str,...];relationship_notes:tuple[str,...];formal_propagation:bool=False
 @property
 def same_run(self):return self.same_run_id_raw and self.same_run_start_time
@dataclass(frozen=True)
class MzMLSourceMetadataAuditResult:
 file_records:tuple[MzMLSourceMetadataRecord,...];relationship_records:tuple[MzMLFileRelationship,...];summary:tuple[tuple[str,int],...];formal_propagation:bool=False

def _empty_record(path:Path,status:ReadStatus,blocks):
 exists=path.exists();regular=path.is_file() if exists else False
 try:st=path.stat() if exists else None
 except OSError:st=None
 return MzMLSourceMetadataRecord(input_path=str(path),file_name=path.name,file_size_bytes=st.st_size if st else "NOT_AVAILABLE",file_exists=exists,is_regular_file=regular,is_symlink=path.is_symlink(),device_id=st.st_dev if st else "NOT_AVAILABLE",inode=st.st_ino if st else "NOT_AVAILABLE",read_status=status,block_reasons=_ordered_blocks(blocks))
def parse_sciex_mzml_source_metadata(path:Path,*,calculate_hash:bool=True):
 path=Path(path)
 if not path.exists():return _empty_record(path,ReadStatus.FILE_NOT_FOUND,("FILE_NOT_FOUND",))
 if not path.is_file():return _empty_record(path,ReadStatus.FILE_NOT_REGULAR,("FILE_NOT_REGULAR",))
 try:st=path.stat();digest=calculate_sha256(path) if calculate_hash else "NOT_CALCULATED"
 except OSError:return _empty_record(path,ReadStatus.FILE_UNREADABLE,("FILE_UNREADABLE",))
 mzml_version=mzml_id="NOT_RECORDED";run_raw=run_start=default_source=default_ic=default_dp="NOT_RECORDED";declared="NOT_RECORDED";root_name=None;mzml_seen=False;spectrum_list_seen=False;sources=[];software=[];methods=[];instrument_ids=[];instrument_models=[];components=[];ms_counts=Counter();missing_ms=positive=negative=mixed_pol=missing_pol=profile=centroid=mixed_rep=missing_rep=parsed=0;scan_times=[];scan_units=set();low=[];high=[];window_low=[];window_high=[];range_sources=set();blocks=[]
 try:
  for event,elem in ET.iterparse(path,events=("start","end")):
   tag=_local(elem.tag)
   if root_name is None and event=="start":root_name=tag
   if event=="start":
    if tag=="mzML":mzml_seen=True;mzml_version=str(elem.get("version") or "NOT_RECORDED");mzml_id=str(elem.get("id") or "NOT_RECORDED")
    elif tag=="run":run_raw=str(elem.get("id") or "NOT_RECORDED");run_start=str(elem.get("startTimeStamp") or "NOT_RECORDED");default_source=str(elem.get("defaultSourceFileRef") or "NOT_RECORDED");default_ic=str(elem.get("defaultInstrumentConfigurationRef") or "NOT_RECORDED")
    elif tag=="spectrumList":spectrum_list_seen=True;raw=elem.get("count");declared=int(raw) if raw and raw.isdigit() else "NOT_RECORDED";default_dp=str(elem.get("defaultDataProcessingRef") or "NOT_RECORDED")
    continue
   if tag=="binary":elem.text=None;elem.clear();continue
   if tag=="sourceFile":sources.append(SourceFileMetadata(str(elem.get("id") or "NOT_RECORDED"),str(elem.get("name") or "NOT_RECORDED"),str(elem.get("location") or "NOT_RECORDED"),normalize_source_location_for_comparison(elem.get("location")),tuple((a,n,v) for a,n,v,_,_ in _cv(elem))));elem.clear()
   elif tag=="software":
    cvs=_cv(elem);name=next((n for _,n,_,_,_ in cvs if n),"NOT_RECORDED");software.append(SoftwareMetadata(str(elem.get("id") or "NOT_RECORDED"),name,str(elem.get("version") or "NOT_RECORDED"),tuple((a,n,v) for a,n,v,_,_ in cvs)));elem.clear()
   elif tag=="instrumentConfiguration":
    iid=str(elem.get("id") or "NOT_RECORDED");instrument_ids.append(iid);direct=[x for x in elem if _local(x.tag)=="cvParam"]
    instrument_models.extend((str(x.get("accession") or ""),str(x.get("name") or "")) for x in direct if x.get("name") and x.get("accession") not in {"MS:1000529"})
    for child in elem.iter():
     kind=_local(child.tag)
     if kind in {"source","analyzer","detector"}:components.append(InstrumentComponentMetadata(iid,kind,str(child.get("order") or "NOT_RECORDED"),tuple((a,n,v) for a,n,v,_,_ in _cv(child))))
    elem.clear()
   elif tag=="dataProcessing":
    did=str(elem.get("id") or "NOT_RECORDED")
    for child in elem.iter():
     if _local(child.tag)=="processingMethod":methods.append(ProcessingMethodMetadata(did,str(child.get("order") or "NOT_RECORDED"),str(child.get("softwareRef") or "NOT_RECORDED"),tuple((a,n,v) for a,n,v,_,_ in _cv(child)),_users(child)))
    elem.clear()
   elif tag=="spectrum":
    parsed+=1;cvs=_cv(elem);acc={a for a,_,_,_,_ in cvs};ms=[int(float(v)) for a,_,v,_,_ in cvs if a=="MS:1000511" and _float(v) is not None]
    if ms:ms_counts[ms[0]]+=1
    else:missing_ms+=1
    pos="MS:1000130" in acc;neg="MS:1000129" in acc
    if pos and neg:mixed_pol+=1
    elif pos:positive+=1
    elif neg:negative+=1
    else:missing_pol+=1
    pro="MS:1000128" in acc;cen="MS:1000127" in acc
    if pro and cen:mixed_rep+=1
    elif pro:profile+=1
    elif cen:centroid+=1
    else:missing_rep+=1
    for a,_,v,uacc,uname in cvs:
     val=_float(v)
     if a=="MS:1000016" and val is not None:
      unit=uacc or uname or "NOT_RECORDED";scan_units.add(unit);minutes=val/60 if uacc=="UO:0000010" or "second" in uname.lower() else val;scan_times.append(minutes)
     elif a=="MS:1000528" and val is not None:low.append(val);range_sources.add(a)
     elif a=="MS:1000527" and val is not None:high.append(val);range_sources.add(a)
     elif a=="MS:1000501" and val is not None:window_low.append(val);range_sources.add(a)
     elif a=="MS:1000500" and val is not None:window_high.append(val);range_sources.add(a)
    elem.clear()
  if not mzml_seen or root_name not in {"mzML","indexedmzML"}:return MzMLSourceMetadataRecord(input_path=str(path),file_name=path.name,file_size_bytes=st.st_size,sha256=digest,file_exists=True,is_regular_file=True,is_symlink=path.is_symlink(),device_id=st.st_dev,inode=st.st_ino,read_status=ReadStatus.INVALID_MZML,block_reasons=("INVALID_MZML_ROOT",))
 except ET.ParseError:return MzMLSourceMetadataRecord(input_path=str(path),file_name=path.name,file_size_bytes=st.st_size,sha256=digest,file_exists=True,is_regular_file=True,is_symlink=path.is_symlink(),device_id=st.st_dev,inode=st.st_ino,read_status=ReadStatus.XML_PARSE_ERROR,block_reasons=("XML_PARSE_ERROR",))
 decoded,decode_status=decode_mzml_xml_safe_id(run_raw);source_ids={x.source_file_id for x in sources};default_resolved=default_source in source_ids if default_source!="NOT_RECORDED" else False;ic_resolved=default_ic in set(instrument_ids) if default_ic!="NOT_RECORDED" else False
 if run_raw=="NOT_RECORDED":blocks.append("MISSING_RUN_ID")
 if run_start=="NOT_RECORDED":blocks.append("MISSING_RUN_START_TIMESTAMP")
 if not sources:blocks.append("MISSING_SOURCE_FILE")
 if len(sources)>1:blocks.append("MULTIPLE_SOURCE_FILES")
 if default_source!="NOT_RECORDED" and not default_resolved:blocks.extend(("DEFAULT_SOURCE_REF_UNRESOLVED","SOURCE_FILE_REFERENCE_UNRESOLVED"))
 if not spectrum_list_seen:blocks.append("MISSING_SPECTRUM_LIST")
 consistent=declared==parsed if isinstance(declared,int) else "NOT_AVAILABLE"
 if consistent is False:blocks.append("DECLARED_PARSED_SPECTRUM_COUNT_MISMATCH")
 if missing_ms:blocks.append("MISSING_MS_LEVEL_METADATA")
 if len(ms_counts)>1:blocks.append("MIXED_MS_LEVEL_METADATA")
 if missing_pol:blocks.append("MISSING_POLARITY_METADATA")
 if positive and negative:blocks.append("MIXED_POLARITY_METADATA")
 if mixed_pol:blocks.append("CONFLICTING_POLARITY_WITHIN_SPECTRUM")
 if missing_rep:blocks.append("MISSING_SPECTRUM_REPRESENTATION")
 if profile and centroid:blocks.append("MIXED_SPECTRUM_REPRESENTATION")
 if mixed_rep:blocks.append("CONFLICTING_REPRESENTATION_WITHIN_SPECTRUM")
 if not scan_times:blocks.append("MISSING_SCAN_TIME")
 if len(scan_units)>1:blocks.append("MIXED_SCAN_TIME_UNITS")
 if not range_sources:blocks.append("MISSING_MZ_RANGE_METADATA")
 if not instrument_ids:blocks.append("MISSING_INSTRUMENT_METADATA")
 if not methods:blocks.append("MISSING_PROCESSING_METADATA")
 pol_status=PolarityStatus.CONFLICTING_WITHIN_SPECTRUM if mixed_pol else PolarityStatus.MIXED_POLARITY if positive and negative else PolarityStatus.NEGATIVE_ONLY if negative and not positive and not missing_pol else PolarityStatus.POSITIVE_ONLY if positive and not negative and not missing_pol else PolarityStatus.NOT_RECORDED
 rep_status=RepresentationStatus.CONFLICTING_WITHIN_SPECTRUM if mixed_rep else RepresentationStatus.MIXED_REPRESENTATION if profile and centroid else RepresentationStatus.PROFILE_ONLY if profile and not centroid and not missing_rep else RepresentationStatus.CENTROID_ONLY if centroid and not profile and not missing_rep else RepresentationStatus.NOT_RECORDED
 other=sum(v for k,v in ms_counts.items() if k not in {1,2});all_range=[bool(low),bool(high),bool(window_low),bool(window_high)];range_status=MZRangeStatus.RECORDED if all(all_range) else MZRangeStatus.PARTIALLY_RECORDED if any(all_range) else MZRangeStatus.NOT_RECORDED;softmap={x.software_id:x for x in software};conversion_ids=sorted({m.software_ref for m in methods if any(a=="MS:1000544" or "conversion" in n.lower() for a,n,_ in m.cv_params)});conversion_names=tuple(softmap[x].software_id if x in softmap else x for x in conversion_ids);conversion_versions=tuple(softmap[x].version if x in softmap else "NOT_RECORDED" for x in conversion_ids);acquisition=tuple(sorted(x.software_id for x in software if any("Analyst" in n for _,n,_ in x.cv_params) or x.software_id.lower()=="analyst"));filters=tuple(sorted(v for m in methods for n,v in m.user_params if "filter" in n.lower()));instrument_names=tuple(sorted({n for _,n in instrument_models}));confidence="HIGH" if not any(x in blocks for x in ("MISSING_RUN_ID","MISSING_SOURCE_FILE","MISSING_SPECTRUM_LIST")) else "MODERATE" if parsed else "LOW"
 return MzMLSourceMetadataRecord(input_path=str(path),file_name=path.name,file_size_bytes=st.st_size,sha256=digest,file_exists=True,is_regular_file=True,is_symlink=path.is_symlink(),device_id=st.st_dev,inode=st.st_ino,read_status=ReadStatus.COMPLETED,mzml_version=mzml_version,mzml_id=mzml_id,run_id_raw=run_raw,run_id_decoded=decoded,run_id_decode_status=decode_status,run_start_time_raw=run_start,run_start_time_normalized=_normalize_timestamp(None if run_start=="NOT_RECORDED" else run_start),default_source_file_ref=default_source,default_instrument_configuration_ref=default_ic,default_data_processing_ref=default_dp,declared_spectrum_count=declared,parsed_spectrum_count=parsed,spectrum_count_consistent=consistent,source_files=tuple(sorted(sources,key=lambda x:x.source_file_id)),default_source_file_resolved=default_resolved,ms_level_set=tuple(sorted(ms_counts)),ms_level_counts=tuple(sorted(ms_counts.items())),ms1_spectrum_count=ms_counts[1],ms2_spectrum_count=ms_counts[2],other_ms_level_count=other,missing_ms_level_count=missing_ms,positive_spectrum_count=positive,negative_spectrum_count=negative,mixed_polarity_spectrum_count=mixed_pol,missing_polarity_count=missing_pol,polarity_status=pol_status,profile_spectrum_count=profile,centroid_spectrum_count=centroid,mixed_representation_spectrum_count=mixed_rep,missing_representation_count=missing_rep,representation_status=rep_status,first_scan_start_time=scan_times[0] if scan_times else "NOT_RECORDED",last_scan_start_time=scan_times[-1] if scan_times else "NOT_RECORDED",minimum_scan_start_time=min(scan_times) if scan_times else "NOT_RECORDED",maximum_scan_start_time=max(scan_times) if scan_times else "NOT_RECORDED",scan_time_unit_set=tuple(sorted(scan_units)),approximate_run_duration=max(scan_times)-min(scan_times) if scan_times else "NOT_RECORDED",lowest_observed_mz_min=min(low) if low else "NOT_RECORDED",lowest_observed_mz_max=max(low) if low else "NOT_RECORDED",highest_observed_mz_min=min(high) if high else "NOT_RECORDED",highest_observed_mz_max=max(high) if high else "NOT_RECORDED",scan_window_lower_limit_min=min(window_low) if window_low else "NOT_RECORDED",scan_window_lower_limit_max=max(window_low) if window_low else "NOT_RECORDED",scan_window_upper_limit_min=min(window_high) if window_high else "NOT_RECORDED",scan_window_upper_limit_max=max(window_high) if window_high else "NOT_RECORDED",mz_range_source=tuple(sorted(range_sources)),mz_range_status=range_status,software=tuple(sorted(software,key=lambda x:x.software_id)),acquisition_software=acquisition,conversion_software=conversion_names,conversion_software_version=conversion_versions,data_processing_ids=tuple(sorted({m.data_processing_id for m in methods})),processing_methods=tuple(sorted(methods,key=lambda x:(x.data_processing_id,x.order,x.software_ref))),msconvert_filter_history=filters,processing_history_status=MetadataStatus.RECORDED if methods else MetadataStatus.NOT_RECORDED,instrument_configuration_ids=tuple(sorted(instrument_ids)),instrument_model_cv=tuple(sorted({a for a,_ in instrument_models})),instrument_model_name=instrument_names,ion_source_cv=tuple(sorted({a for c in components if c.component_type=="source" for a,_,_ in c.cv_params})),analyzer_cv=tuple(sorted({a for c in components if c.component_type=="analyzer" for a,_,_ in c.cv_params})),detector_cv=tuple(sorted({a for c in components if c.component_type=="detector" for a,_,_ in c.cv_params})),instrument_components=tuple(sorted(components,key=lambda x:(x.instrument_configuration_id,int(x.order) if x.order.isdigit() else 999,x.component_type))),default_instrument_configuration_resolved=ic_resolved,instrument_metadata_status=MetadataStatus.RECORDED if instrument_ids and instrument_names else MetadataStatus.PARTIALLY_RECORDED if instrument_ids else MetadataStatus.NOT_RECORDED,metadata_confidence=confidence,block_reasons=_ordered_blocks(blocks))

def _known(value):return value not in {None,"","NOT_RECORDED","NOT_AVAILABLE","NOT_CALCULATED"}
def _same_known(a,b):return _known(a) and _known(b) and a==b
def compare_mzml_source_metadata(left:MzMLSourceMetadataRecord,right:MzMLSourceMetadataRecord):
 a,b=sorted((left,right),key=lambda x:(x.file_name,x.input_path));unread={ReadStatus.FILE_NOT_FOUND,ReadStatus.FILE_NOT_REGULAR,ReadStatus.FILE_UNREADABLE};invalid={ReadStatus.XML_PARSE_ERROR,ReadStatus.INVALID_MZML};same_sha=_same_known(a.sha256,b.sha256);same_inode=_known(a.inode) and _known(b.inode) and a.inode==b.inode and a.device_id==b.device_id;same_size=_same_known(a.file_size_bytes,b.file_size_bytes);an={x.casefold() for x in a.source_file_names if _known(x)};bn={x.casefold() for x in b.source_file_names if _known(x)};al={x.location_normalized_for_comparison.casefold() for x in a.source_files if x.location_normalized_for_comparison};bl={x.location_normalized_for_comparison.casefold() for x in b.source_files if x.location_normalized_for_comparison};same_source=bool(an and bn and an==bn);same_location=bool(al and bl and al==bl);same_default=_same_known(a.default_source_file_ref,b.default_source_file_ref);same_raw=_same_known(a.run_id_raw,b.run_id_raw);same_decoded=_same_known(a.run_id_decoded,b.run_id_decoded);same_start=_same_known(a.run_start_time_normalized,b.run_start_time_normalized);same_declared=_same_known(a.declared_spectrum_count,b.declared_spectrum_count);same_parsed=a.parsed_spectrum_count==b.parsed_spectrum_count;same_ms=a.ms_level_counts==b.ms_level_counts;same_pol=(a.positive_spectrum_count,a.negative_spectrum_count,a.mixed_polarity_spectrum_count,a.missing_polarity_count)==(b.positive_spectrum_count,b.negative_spectrum_count,b.mixed_polarity_spectrum_count,b.missing_polarity_count);same_rep=(a.profile_spectrum_count,a.centroid_spectrum_count,a.mixed_representation_spectrum_count,a.missing_representation_count)==(b.profile_spectrum_count,b.centroid_spectrum_count,b.mixed_representation_spectrum_count,b.missing_representation_count);same_conv=(a.conversion_software,a.conversion_software_version,a.processing_methods)==(b.conversion_software,b.conversion_software_version,b.processing_methods);same_ic=(a.default_instrument_configuration_ref,a.instrument_configuration_ids,a.instrument_model_cv)==(b.default_instrument_configuration_ref,b.instrument_configuration_ids,b.instrument_model_cv);run_diff=(_known(a.run_id_raw) and _known(b.run_id_raw) and not same_raw) or (_known(a.run_start_time_normalized) and _known(b.run_start_time_normalized) and not same_start);notes=[];blocks=[]
 if a.read_status in unread or b.read_status in unread:status=RelationshipStatus.UNREADABLE_FILE;confidence="HIGH";blocks.append("FILE_UNREADABLE")
 elif a.read_status in invalid or b.read_status in invalid:status=RelationshipStatus.INVALID_MZML;confidence="HIGH";blocks.append("XML_PARSE_ERROR" if ReadStatus.XML_PARSE_ERROR in {a.read_status,b.read_status} else "INVALID_MZML_ROOT")
 elif same_inode:status=RelationshipStatus.SAME_FILE_OBJECT;confidence="HIGH";notes.append("SAME_DEVICE_AND_INODE")
 elif same_sha:status=RelationshipStatus.EXACT_DUPLICATE;confidence="HIGH";notes.append("IDENTICAL_SHA256")
 elif same_source and same_location and run_diff:status=RelationshipStatus.DIFFERENT_RUN_SAME_SOURCE_WIFF;confidence="HIGH" if (not same_parsed or not same_ms) else "MODERATE";notes.append("SOURCE_FILE_AND_LOCATION_MATCH_RUN_IDENTITY_DIFFERS")
 elif same_source and run_diff:status=RelationshipStatus.POSSIBLE_DIFFERENT_RUN_SAME_SOURCE;confidence="MODERATE";notes.append("SOURCE_FILE_MATCH_LOCATION_INCOMPLETE_OR_DIFFERENT")
 elif same_source and (same_raw or same_decoded) and same_start and same_ms and same_pol and same_rep:
  status=RelationshipStatus.SAME_RUN_DIFFERENT_EXPORT if not same_conv and a.processing_history_status is MetadataStatus.RECORDED and b.processing_history_status is MetadataStatus.RECORDED else RelationshipStatus.POSSIBLE_REEXPORT_OF_SAME_RUN;confidence="HIGH" if status is RelationshipStatus.SAME_RUN_DIFFERENT_EXPORT else "MODERATE";notes.append("RUN_ID_START_AND_SPECTRUM_DISTRIBUTIONS_MATCH")
 elif an and bn and an!=bn:status=RelationshipStatus.DIFFERENT_SOURCE_FILE;confidence="HIGH";notes.append("SOURCE_FILE_NAMES_DIFFER")
 else:status=RelationshipStatus.INSUFFICIENT_SOURCE_METADATA;confidence="LOW";blocks.extend(("INSUFFICIENT_RUN_IDENTITY","INSUFFICIENT_SOURCE_IDENTITY"))
 keep=status not in {RelationshipStatus.EXACT_DUPLICATE,RelationshipStatus.SAME_FILE_OBJECT}
 return MzMLFileRelationship(left_file=a.file_name,right_file=b.file_name,exact_duplicate=same_sha,same_file_size=same_size,same_sha256=same_sha,same_inode=same_inode,same_source_file=same_source,same_source_location=same_location,same_default_source_ref=same_default,same_run_id_raw=same_raw,same_run_id_decoded=same_decoded,same_run_start_time=same_start,same_declared_spectrum_count=same_declared,same_parsed_spectrum_count=same_parsed,same_ms_level_distribution=same_ms,same_polarity_distribution=same_pol,same_representation_distribution=same_rep,same_conversion_software=same_conv,same_instrument_configuration=same_ic,duplicate_file_status=status,keep_both_files=keep,relationship_confidence=confidence,relationship_block_reasons=_ordered_blocks(blocks),relationship_notes=tuple(sorted(notes)))
def apply_runtime_source_context(record:MzMLSourceMetadataRecord,context:RuntimeSourceContext):
 existing=record.source_confirmed_context;has_internal=_known(existing) and existing not in {"UNKNOWN","UNCONFIRMED_FROM_MZML_METADATA"};manifest_values={context.rna_identity,context.digest_type,context.technical_run_label};conflict=has_internal and existing not in manifest_values;blocks=_ordered_blocks(record.block_reasons+(("USER_MANIFEST_METADATA_CONFLICT",) if conflict else ()));manifest_context=f"RNA_Identity={context.rna_identity};Digest_Type={context.digest_type};Technical_Run_Label={context.technical_run_label}"
 return replace(record,rna_identity=context.rna_identity,digest_type=context.digest_type,technical_run_label=context.technical_run_label,context_source=context.context_source,context_confidence=context.context_confidence,mzml_metadata_confirmed=context.mzml_metadata_confirmed,rna_identity_confirmed=True,source_confirmed_context=existing if conflict else manifest_context,filename_label_only="NOT_APPLICABLE_USER_MANIFEST_SUPPLIED",context_conflict=conflict,block_reasons=blocks)
def audit_mzml_source_metadata_files(paths:Sequence[Path],*,calculate_hash:bool=True,runtime_contexts:Mapping[str,RuntimeSourceContext]|None=None):
 contexts=runtime_contexts or {};parsed=[]
 for p in paths:
  record=parse_sciex_mzml_source_metadata(Path(p),calculate_hash=calculate_hash);context=contexts.get(record.input_path) or contexts.get(record.file_name);parsed.append(apply_runtime_source_context(record,context) if context else record)
 records=tuple(sorted(parsed,key=lambda x:(x.file_name,x.input_path)));relationships=tuple(compare_mzml_source_metadata(a,b) for a,b in combinations(records,2));counts=Counter(x.duplicate_file_status.value for x in relationships);summary=(("File_Count",len(records)),("Relationship_Count",len(relationships)))+tuple(sorted(counts.items()))
 return MzMLSourceMetadataAuditResult(records,relationships,summary)
def file_metadata_row(record:MzMLSourceMetadataRecord):
 return {"File_Name":record.file_name,"File_Size_Bytes":record.file_size_bytes,"SHA256":record.sha256,"Read_Status":record.read_status.value,"Run_ID_Raw":record.run_id_raw,"Run_ID_Decoded":record.run_id_decoded,"Run_Start_Time":record.run_start_time_normalized,"Source_File_Count":record.source_file_count,"Source_File_IDs":";".join(record.source_file_ids),"Source_File_Names":";".join(record.source_file_names),"Source_File_Locations":";".join(record.source_file_locations),"Source_File_Accession_Metadata":repr(record.source_file_accession_metadata),"Declared_Spectrum_Count":record.declared_spectrum_count,"Parsed_Spectrum_Count":record.parsed_spectrum_count,"MS_Level_Set":";".join(map(str,record.ms_level_set)),"MS1_Spectrum_Count":record.ms1_spectrum_count,"MS2_Spectrum_Count":record.ms2_spectrum_count,"Other_MS_Level_Count":record.other_ms_level_count,"Polarity_Status":record.polarity_status.value,"Positive_Spectrum_Count":record.positive_spectrum_count,"Negative_Spectrum_Count":record.negative_spectrum_count,"Representation_Status":record.representation_status.value,"Profile_Spectrum_Count":record.profile_spectrum_count,"Centroid_Spectrum_Count":record.centroid_spectrum_count,"First_Scan_Start_Time":record.first_scan_start_time,"Last_Scan_Start_Time":record.last_scan_start_time,"MZ_Range_Status":record.mz_range_status.value,"Instrument_Model_Name":";".join(record.instrument_model_name),"Software_IDs":";".join(record.software_ids),"Software_Names":";".join(record.software_names),"Software_Versions":";".join(record.software_versions),"Conversion_Software":";".join(record.conversion_software),"Conversion_Software_Version":";".join(record.conversion_software_version),"Processing_Method_Count":record.processing_method_count,"Processing_Method_Orders":";".join(record.processing_method_orders),"Processing_Method_Software_Refs":";".join(record.processing_method_software_refs),"Processing_CV_Params":repr(record.processing_cv_params),"Processing_User_Params":repr(record.processing_user_params),"Processing_History_Status":record.processing_history_status.value,"Metadata_Confidence":record.metadata_confidence,"RNA_Identity":record.rna_identity,"Digest_Type":record.digest_type,"Technical_Run_Label":record.technical_run_label,"Context_Source":record.context_source,"Context_Confidence":record.context_confidence,"MzML_Metadata_Confirmed":record.mzml_metadata_confirmed,"Block_Reasons":";".join(record.block_reasons),"Formal_Propagation":False}
def relationship_row(record:MzMLFileRelationship):
 return {"Left_File":record.left_file,"Right_File":record.right_file,"Same_SHA256":record.same_sha256,"Same_Source_File":record.same_source_file,"Same_Source_Location":record.same_source_location,"Same_Run_ID":record.same_run_id_raw,"Same_Run_Start_Time":record.same_run_start_time,"Same_Spectrum_Count":record.same_parsed_spectrum_count,"Same_MS_Level_Distribution":record.same_ms_level_distribution,"Same_Polarity_Distribution":record.same_polarity_distribution,"Duplicate_File_Status":record.duplicate_file_status.value,"Keep_Both_Files":record.keep_both_files,"Relationship_Confidence":record.relationship_confidence,"Block_Reasons":";".join(record.relationship_block_reasons),"Formal_Propagation":False}
def audit_optional_result(result:MzMLSourceMetadataAuditResult):
 return {"file_records":[file_metadata_row(x) for x in result.file_records],"relationship_records":[relationship_row(x) for x in result.relationship_records],"summary":dict(result.summary)}
