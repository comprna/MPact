"""Shared compact output schemas and bundled checkpoint selection."""
from pathlib import Path

SUPPORTED_WINDOWS = (101, 201, 501, 801, 1001)
ACCESS_SUFFIXES = ('unpaired_probability_1nt', 'accessibility_5nt',
                   'accessibility_10nt', 'accessibility_20nt')


def bundled_model_path(window_size, override=None):
    if override:
        return str(override)
    return str(Path(__file__).resolve().parent / 'models' /
               f'mpact_dtm6a_{window_size}nt_seed42.keras')


def variant_output_columns(input_columns, accessibility=False, conservation=False, full=False):
    columns = list(input_columns) + [
        'strand_gencode', 'score_type', 'a_genomic_pos1', 'snp_to_a_mRNA_offset',
        'ref_center_is_A', 'alt_center_is_A',
        'mpact_ref_score', 'mpact_alt_score', 'mpact_delta_alt_minus_ref', 'mpact_abs_delta',
        'ref_motif5', 'alt_motif5', 'ref_is_DRACH', 'alt_is_DRACH', 'ref_is_RAC', 'alt_is_RAC',
        'overlaps_AtoI_exact', 'near_AtoI_5nt', 'near_AtoI_10nt',
    ]
    if accessibility:
        columns += [f'mpact_{allele}_{suffix}' for suffix in ACCESS_SUFFIXES
                    for allele in ('ref', 'alt', 'delta')]
    if conservation:
        columns += ['variant_conservation_score', 'a_site_conservation_score']
    if full:
        columns += ['strand_gencode_source', 'strand_gencode_inference',
                    'ref_scan_seq', 'alt_scan_seq']
        if accessibility:
            columns += ['accessibility_source']
        if conservation:
            columns += ['conservation_source']
    return list(dict.fromkeys(columns))


def site_output_columns(input_columns, accessibility=False, conservation=False, full=False):
    columns = list(input_columns) + [
        'mpact_score', 'center_base_transcript', 'scoreable_center_A',
        'motif5', 'is_DRACH', 'is_RAC',
    ]
    if accessibility:
        columns += [f'mpact_ref_{suffix}' for suffix in ACCESS_SUFFIXES]
    if conservation:
        columns += ['site_conservation_score']
    if full:
        columns += ['mpact_chr', 'mpact_pos1', 'mpact_strand']
        if accessibility:
            columns += ['accessibility_source']
        if conservation:
            columns += ['conservation_source']
    return list(dict.fromkeys(columns))
