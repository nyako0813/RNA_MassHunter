"""Resolution-aware shadow audit for SCIEX integer/isotope-like spacing flags."""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
from math import isfinite
from statistics import mean, median
from types import MappingProxyType
from typing import Any, Mapping, Sequence
AUDIT_RESULT_KEY = 'sciex_spacing_resolution_audit'
SUMMARY_SHEET = 'SCIEX_Spacing_Resolution'
DETAIL_SHEET = 'SCIEX_Spacing_Resolution_Detail'
WARNING_CODE = 'SCIEX_SPACING_NOT_RESOLUTION_DISTINGUISHABLE'
ERROR_CODE = 'SCIEX_SPACING_RESOLUTION_AUDIT_ERROR'
ALGORITHM_VERSION = 'sciex-spacing-resolution-v1'
FORMAL_FALSE = {'Shadow_Only': True, 'Applied_To_Formal_Score': False, 'Applied_To_Ranking': False, 'Applied_To_Candidate_Filtering': False, 'Molecular_Identity_Assigned': False}
DEFAULT_PARAMETERS = {'enabled': True, 'minimum_spacing_sample_count': 20, 'quantization_tolerance_da': 0.02, 'distinguishability_margin_factor': 2.0, 'maximum_spacing_multiple': 10}
SUMMARY_COLUMNS = ['Audit_Status', 'Audit_Eligible', 'SCIEX_Source_File', 'Input_Mass_Point_Count', 'Input_Adjacent_Spacing_Count', 'Input_Grid_Min_Da', 'Input_Grid_Max_Da', 'Input_Grid_Mean_Da', 'Input_Grid_Median_Da', 'Input_Grid_Mode_Da', 'Input_Grid_Unique_Spacing_Count', 'Input_Grid_Is_Uniform', 'Input_Grid_Uniformity_Fraction', 'Detected_Apex_Count', 'Apex_Adjacent_Spacing_Min_Da', 'Apex_Adjacent_Spacing_Median_Da', 'Apex_Adjacent_Spacing_Mode_Da', 'Apex_Decimal_Pattern_Count', 'Apex_Quantization_Step_Da', 'Apex_Quantization_Confidence', 'Estimated_Effective_Grid_Da', 'Grid_Estimation_Source', 'Grid_Estimation_Method', 'Grid_Residual_Median_Da', 'Grid_Residual_Max_Da', 'Grid_Confidence', 'Integer_Spacing_Da', 'Isotope_Spacing_Da', 'Single_Step_Target_Separation_Da', 'Integer_Tolerance_Da', 'Isotope_Tolerance_Da', 'Tolerance_Windows_Overlap', 'Single_Step_Grid_Steps_Per_Separation', 'Resolution_Status', 'Theoretically_Distinguishable', 'Total_Relation_Count', 'Integer_Candidate_Count', 'Isotope_Candidate_Count', 'Dual_Integer_And_Isotope_Count', 'Integer_Only_Count', 'Isotope_Only_Count', 'Neither_Count', 'Dual_Flag_Fraction', 'Observed_Relations_Resolution_Compatible', 'Observed_Dual_Flag_Expected_From_Resolution', 'Isotope_Interpretation_Eligible', 'Numerical_Spacing_Interpretation_Eligible', 'Chemical_Interpretation_Eligible', 'Input_Identity_Audit_Status', 'Input_Identity_Conflict', 'Biological_Interpretation_Eligible', 'Warning_Code', 'Warning_Message', 'Minimum_Spacing_Sample_Count', 'Quantization_Tolerance_Da', 'Distinguishability_Margin_Factor', 'Maximum_Spacing_Multiple', 'Algorithm_Version', 'Shadow_Only', 'Applied_To_Formal_Score', 'Applied_To_Ranking', 'Applied_To_Candidate_Filtering', 'Molecular_Identity_Assigned', 'Notes']
DETAIL_COLUMNS = ['Spacing_Multiple', 'Integer_Target_Da', 'Isotope_Target_Da', 'Target_Separation_Da', 'Integer_Window_Min_Da', 'Integer_Window_Max_Da', 'Isotope_Window_Min_Da', 'Isotope_Window_Max_Da', 'Integer_Tolerance_Da', 'Isotope_Tolerance_Da', 'Tolerance_Window_Overlap_Da', 'Tolerance_Windows_Overlap', 'Estimated_Effective_Grid_Da', 'Grid_Steps_Per_Target_Separation', 'Resolution_Status', 'Theoretically_Distinguishable', 'Resolution_Limited', 'Algorithm_Version', 'Shadow_Only', 'Applied_To_Formal_Score', 'Applied_To_Ranking', 'Applied_To_Candidate_Filtering', 'Molecular_Identity_Assigned', 'Notes']

@dataclass(frozen=True)
class SpacingResolutionParameters:
    minimum_spacing_sample_count: int = 20
    quantization_tolerance_da: float = 0.02
    distinguishability_margin_factor: float = 2.0
    maximum_spacing_multiple: int = 10

    @classmethod
    def from_mapping(cls, value):
        source = dict(value or {})
        source.pop('enabled', None)
        allowed = set(cls.__dataclass_fields__)
        result = cls(**{k: v for k, v in source.items() if k in allowed})
        result.validate()
        return result

    def validate(self):
        if isinstance(self.minimum_spacing_sample_count, bool) or not isinstance(self.minimum_spacing_sample_count, int) or self.minimum_spacing_sample_count < 2:
            raise ValueError('minimum_spacing_sample_count must be an integer >= 2')
        if isinstance(self.maximum_spacing_multiple, bool) or not isinstance(self.maximum_spacing_multiple, int) or self.maximum_spacing_multiple < 1:
            raise ValueError('maximum_spacing_multiple must be an integer >= 1')
        for name in ('quantization_tolerance_da', 'distinguishability_margin_factor'):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or (not isfinite(float(value))) or (value <= 0):
                raise ValueError(f'{name} must be finite and positive')

@dataclass(frozen=True)
class SciexSpacingResolutionAuditResult:
    summary_rows: tuple[Mapping[str, Any], ...]
    detail_rows: tuple[Mapping[str, Any], ...]

    def __post_init__(self):
        object.__setattr__(self, 'summary_rows', tuple((MappingProxyType(dict(x)) for x in self.summary_rows)))
        object.__setattr__(self, 'detail_rows', tuple((MappingProxyType(dict(x)) for x in self.detail_rows)))

    def summaries(self):
        return [dict(x) for x in self.summary_rows]

    def details(self):
        return [dict(x) for x in self.detail_rows]

def _finite_values(values):
    out = []
    for value in values or ():
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if isfinite(number):
            out.append(number)
    return out

def _positive_spacings(values):
    ordered = sorted(set(_finite_values(values)))
    return [b - a for a, b in zip(ordered, ordered[1:]) if b - a > 0]

def _mode(values):
    if not values:
        return None
    rounded = [round(float(x), 6) for x in values]
    counts = Counter(rounded)
    return min(counts, key=lambda x: (-counts[x], x))

def _grid_stats(masses, tolerance):
    points = _finite_values(masses)
    spacings = _positive_spacings(points)
    mode = _mode(spacings)
    fraction = sum((abs(x - mode) <= tolerance for x in spacings)) / len(spacings) if spacings and mode is not None else 0.0
    return {'points': points, 'spacings': spacings, 'min': min(spacings) if spacings else None, 'max': max(spacings) if spacings else None, 'mean': mean(spacings) if spacings else None, 'median': median(spacings) if spacings else None, 'mode': mode, 'unique': len({round(x, 6) for x in spacings}), 'fraction': fraction, 'uniform': bool(spacings and fraction >= 0.99)}

def _residuals(values, step):
    return [abs(x - round(x / step) * step) for x in values]

def _quantization(values, tolerance, minimum_count):
    points = _finite_values(values)
    if len(points) < 2:
        return (None, 'NONE', None, None)
    candidates = (0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0)
    chosen = None
    chosen_res = []
    for step in candidates:
        residuals = _residuals(points, step)
        residual_limit = min(tolerance, step * 0.1)
        fraction = sum((r <= residual_limit for r in residuals)) / len(residuals)
        if fraction >= 0.99:
            chosen = step
            chosen_res = residuals
    if chosen is None:
        return (None, 'NONE', None, None)
    confidence = 'HIGH' if len(points) >= minimum_count and max(chosen_res) <= tolerance else 'MEDIUM' if len(points) >= minimum_count else 'LOW'
    return (chosen, confidence, median(chosen_res), max(chosen_res))

def _records(value, method, attribute):
    source = getattr(value, method)() if hasattr(value, method) else getattr(value, attribute, ())
    return [dict(x) for x in source or () if isinstance(x, Mapping)]

def _resolution_status(overlap, grid, separation, margin):
    if grid is None:
        return 'INSUFFICIENT_INFORMATION'
    ratio = separation / grid
    grid_limited = separation < margin * grid
    if overlap and grid_limited:
        return 'NOT_DISTINGUISHABLE_BOTH'
    if overlap:
        return 'NOT_DISTINGUISHABLE_TOLERANCE_OVERLAP'
    if not grid_limited:
        return 'DISTINGUISHABLE'
    if ratio >= 1:
        return 'MARGINALLY_DISTINGUISHABLE'
    return 'NOT_DISTINGUISHABLE_GRID_LIMITED'

def audit_sciex_spacing_resolution(input_masses: Sequence[Any], apex_masses: Sequence[Any], cluster_result: Any, cluster_parameters: Mapping[str, Any], parameters: SpacingResolutionParameters | Mapping[str, Any] | None=None, *, source_file: str=''):
    params = parameters if isinstance(parameters, SpacingResolutionParameters) else SpacingResolutionParameters.from_mapping(parameters)
    params.validate()
    isotope = float(cluster_parameters.get('isotope_spacing_da', 1.003355))
    iso_tol = float(cluster_parameters.get('isotope_spacing_tolerance_da', 0.15))
    int_tol = float(cluster_parameters.get('integer_spacing_tolerance_da', 0.15))
    for value in (isotope, iso_tol, int_tol):
        if not isfinite(value) or value <= 0:
            raise ValueError('cluster spacing parameters must be finite and positive')
    grid = _grid_stats(input_masses, params.quantization_tolerance_da)
    apex = _grid_stats(apex_masses, params.quantization_tolerance_da)
    apex_step, apex_conf, _, _ = _quantization(apex['points'], params.quantization_tolerance_da, params.minimum_spacing_sample_count)
    effective = None
    source = ''
    method = ''
    confidence = 'NONE'
    if len(grid['spacings']) >= params.minimum_spacing_sample_count and grid['fraction'] >= 0.9:
        effective = grid['mode']
        source = 'INPUT_PROFILE_GRID'
        method = 'ROUNDED_ADJACENT_SPACING_MODE'
        confidence = 'HIGH' if grid['fraction'] >= 0.99 else 'MEDIUM'
    elif apex_step is not None and len(apex['points']) >= params.minimum_spacing_sample_count:
        effective = apex_step
        source = 'DETECTED_APEX_QUANTIZATION'
        method = 'LARGEST_CANDIDATE_STEP_WITH_RESIDUAL_COVERAGE'
        confidence = apex_conf
    elif grid['median'] is not None:
        effective = grid['median']
        source = 'INPUT_PROFILE_GRID'
        method = 'ADJACENT_SPACING_MEDIAN_FALLBACK'
        confidence = 'LOW'
    elif apex_step is not None:
        effective = apex_step
        source = 'DETECTED_APEX_QUANTIZATION'
        method = 'CANDIDATE_STEP_FALLBACK'
        confidence = 'LOW'
    residuals = _residuals(grid['points'] if source == 'INPUT_PROFILE_GRID' else apex['points'], effective) if effective else []
    relation_rows = _records(cluster_result, 'relations', 'relation_rows')
    cluster_summaries = _records(cluster_result, 'summaries', 'summary_rows')
    cluster_summary = cluster_summaries[0] if cluster_summaries else {}
    detail = []
    for n in range(1, params.maximum_spacing_multiple + 1):
        integer = float(n)
        iso = n * isotope
        separation = abs(iso - integer)
        imin, imax = (integer - int_tol, integer + int_tol)
        smin, smax = (iso - iso_tol, iso + iso_tol)
        overlap_amount = max(0.0, min(imax, smax) - max(imin, smin))
        overlap = max(imin, smin) <= min(imax, smax)
        status = _resolution_status(overlap, effective, separation, params.distinguishability_margin_factor)
        dist = status in {'DISTINGUISHABLE', 'MARGINALLY_DISTINGUISHABLE'}
        detail.append({'Spacing_Multiple': n, 'Integer_Target_Da': integer, 'Isotope_Target_Da': iso, 'Target_Separation_Da': separation, 'Integer_Window_Min_Da': imin, 'Integer_Window_Max_Da': imax, 'Isotope_Window_Min_Da': smin, 'Isotope_Window_Max_Da': smax, 'Integer_Tolerance_Da': int_tol, 'Isotope_Tolerance_Da': iso_tol, 'Tolerance_Window_Overlap_Da': overlap_amount, 'Tolerance_Windows_Overlap': overlap, 'Estimated_Effective_Grid_Da': effective, 'Grid_Steps_Per_Target_Separation': separation / effective if effective else None, 'Resolution_Status': status, 'Theoretically_Distinguishable': dist, 'Resolution_Limited': not dist, 'Algorithm_Version': ALGORITHM_VERSION, **FORMAL_FALSE, 'Notes': 'Target-window and grid-resolution diagnostic only; no isotope or chemical assignment.'})
    first = detail[0]
    integer_count = sum((bool(x.get('Integer_Spacing_Candidate')) for x in relation_rows))
    isotope_count = sum((bool(x.get('Isotope_Spacing_Candidate')) for x in relation_rows))
    dual = sum((bool(x.get('Integer_Spacing_Candidate')) and bool(x.get('Isotope_Spacing_Candidate')) for x in relation_rows))
    integer_only = integer_count - dual
    isotope_only = isotope_count - dual
    neither = len(relation_rows) - integer_only - isotope_only - dual
    not_dist = first['Resolution_Status'].startswith('NOT_DISTINGUISHABLE')
    warning = bool(not_dist and isotope_count and dual)
    message = 'The current SCIEX mass grid and spacing tolerances cannot distinguish integer-Da spacing from 1.003355-Da isotope-like spacing. Isotope-like relation flags are numerical proximity diagnostics only and must not be interpreted as isotope assignments.' if warning else ''
    summary = {'Audit_Status': 'AUDIT_COMPLETED', 'Audit_Eligible': True, 'SCIEX_Source_File': source_file, 'Input_Mass_Point_Count': len(grid['points']), 'Input_Adjacent_Spacing_Count': len(grid['spacings']), 'Input_Grid_Min_Da': grid['min'], 'Input_Grid_Max_Da': grid['max'], 'Input_Grid_Mean_Da': grid['mean'], 'Input_Grid_Median_Da': grid['median'], 'Input_Grid_Mode_Da': grid['mode'], 'Input_Grid_Unique_Spacing_Count': grid['unique'], 'Input_Grid_Is_Uniform': grid['uniform'], 'Input_Grid_Uniformity_Fraction': grid['fraction'], 'Detected_Apex_Count': len(apex['points']), 'Apex_Adjacent_Spacing_Min_Da': apex['min'], 'Apex_Adjacent_Spacing_Median_Da': apex['median'], 'Apex_Adjacent_Spacing_Mode_Da': apex['mode'], 'Apex_Decimal_Pattern_Count': len({round(x % 1, 6) for x in apex['points']}), 'Apex_Quantization_Step_Da': apex_step, 'Apex_Quantization_Confidence': apex_conf, 'Estimated_Effective_Grid_Da': effective, 'Grid_Estimation_Source': source, 'Grid_Estimation_Method': method, 'Grid_Residual_Median_Da': median(residuals) if residuals else None, 'Grid_Residual_Max_Da': max(residuals) if residuals else None, 'Grid_Confidence': confidence, 'Integer_Spacing_Da': 1.0, 'Isotope_Spacing_Da': isotope, 'Single_Step_Target_Separation_Da': first['Target_Separation_Da'], 'Integer_Tolerance_Da': int_tol, 'Isotope_Tolerance_Da': iso_tol, 'Tolerance_Windows_Overlap': first['Tolerance_Windows_Overlap'], 'Single_Step_Grid_Steps_Per_Separation': first['Grid_Steps_Per_Target_Separation'], 'Resolution_Status': first['Resolution_Status'], 'Theoretically_Distinguishable': first['Theoretically_Distinguishable'], 'Total_Relation_Count': len(relation_rows), 'Integer_Candidate_Count': integer_count, 'Isotope_Candidate_Count': isotope_count, 'Dual_Integer_And_Isotope_Count': dual, 'Integer_Only_Count': integer_only, 'Isotope_Only_Count': isotope_only, 'Neither_Count': neither, 'Dual_Flag_Fraction': dual / len(relation_rows) if relation_rows else 0.0, 'Observed_Relations_Resolution_Compatible': not warning or dual > 0, 'Observed_Dual_Flag_Expected_From_Resolution': bool(not_dist and dual), 'Isotope_Interpretation_Eligible': bool(first['Theoretically_Distinguishable'] and isotope_count), 'Numerical_Spacing_Interpretation_Eligible': bool(relation_rows), 'Chemical_Interpretation_Eligible': False, 'Input_Identity_Audit_Status': cluster_summary.get('Input_Identity_Audit_Status', 'NOT_RUN'), 'Input_Identity_Conflict': bool(cluster_summary.get('Input_Identity_Conflict', False)), 'Biological_Interpretation_Eligible': bool(cluster_summary.get('Biological_Interpretation_Eligible', False)), 'Warning_Code': WARNING_CODE if warning else '', 'Warning_Message': message, 'Minimum_Spacing_Sample_Count': params.minimum_spacing_sample_count, 'Quantization_Tolerance_Da': params.quantization_tolerance_da, 'Distinguishability_Margin_Factor': params.distinguishability_margin_factor, 'Maximum_Spacing_Multiple': params.maximum_spacing_multiple, 'Algorithm_Version': ALGORITHM_VERSION, **FORMAL_FALSE, 'Notes': 'Resolution-aware interpretation audit only; existing detector, comparison, cluster, and relation flags are unchanged.'}
    return SciexSpacingResolutionAuditResult((summary,), tuple(detail))

def annotate_cluster_summary(cluster_result, spacing_result):
    from rna_masshunter.sciex_delta_mass_cluster_audit import SciexDeltaMassClusterAuditResult
    summaries = _records(cluster_result, 'summaries', 'summary_rows')
    spacing = _records(spacing_result, 'summaries', 'summary_rows')
    meta = spacing[0] if spacing else {}
    updated = []
    for row in summaries:
        copy = dict(row)
        copy.update({'Spacing_Resolution_Audit_Status': meta.get('Resolution_Status', 'NOT_RUN'), 'Estimated_Effective_Grid_Da': meta.get('Estimated_Effective_Grid_Da'), 'Integer_Isotope_Distinguishable': bool(meta.get('Theoretically_Distinguishable', False)), 'Isotope_Interpretation_Eligible': bool(meta.get('Isotope_Interpretation_Eligible', False)), 'Spacing_Resolution_Warning_Code': meta.get('Warning_Code', '')})
        updated.append(copy)
    return SciexDeltaMassClusterAuditResult(tuple(_records(cluster_result, 'clusters', 'cluster_rows')), tuple(updated), tuple(_records(cluster_result, 'relations', 'relation_rows')))
