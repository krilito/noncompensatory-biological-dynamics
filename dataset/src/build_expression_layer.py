"""Build the PAIR v0.3 Level C1 technology-native quantitative expression layer.

Scope: give every expression source a numerically valid scale *within that source* on
canonical-gene coordinates, using the supplied contract in ``expression_transforms.py``
as the mathematical authority.

Explicitly out of scope (Level C2 and later): cross-cohort batch correction, ComBat,
global z-scoring, joint quantile normalization, ranks, paired deltas, feature
selection, models.  Nothing here makes two cohorts comparable; it makes each cohort
internally correct and gene-addressable.

Derived matrices are written locally and are not committed; only code, contracts,
metrics, manifests and documentation are published.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

import build_gene_space as bg
from expression_transforms import (
    ExpressionScale, SourceExpressionContract, build_canonical_quantitative_matrix)

DATASET_ROOT = bg.DATASET_ROOT
LAYER = DATASET_ROOT / 'expression_layer'
MANIFESTS = LAYER / 'manifests'
# Local only, until redistribution is authorized by the Owner.
MATRIX_DIR = bg.PROJECT / 'longitudinal-data' / 'data' / 'processed' / 'expression_v0.3'

GEO = 'GEO series-matrix !Sample_data_processing (verbatim)'
GEO_SAMPLE = 'GEO sample metadata field value (verbatim)'
FILENAME = 'file name as distributed by the submitter'
OUR_SCRIPT = 'derivation code in this repository (we produced this file)'

# One entry per expression source.  ``scale`` is adjudicated from documentary
# provenance only; numeric magnitude plays no part in the declaration (validate_scale
# merely checks that a declaration is not contradicted by the data).
CONTRACTS: dict[str, dict] = {
    # ------------------------------------------------------------------ bulk RNA-seq
    'GIDE_bulk': dict(
        modality='bulk RNA-seq', modality_class='BULK_RNA',
        scale=ExpressionScale.UNKNOWN, probe_based=False, binding='NATIVE_ID_HEADER',
        provenance_kind='no admissible provenance',
        provenance="the file is titled 'cancercell_normalized_counts_genenames' and is "
                   'non-integer, while the only local documentation '
                   '(data/raw/PRJEB23709/README.md) describes the accession as "the read '
                   'counts obtained from RNAseq data" and names neither the normalization '
                   'applied nor its divisor. Declaring RAW_COUNTS would fail the integer '
                   'contract; declaring any linear abundance would license a second, '
                   'undocumented normalization of already-normalized data'),
    'GSE91061_bulk': dict(
        modality='bulk RNA-seq', modality_class='BULK_RNA',
        scale=ExpressionScale.FPKM, probe_based=False, binding='NATIVE_ID_HEADER',
        provenance_kind=FILENAME,
        provenance="file name 'GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz' declares "
                   'FPKM on the hg19 KnownGene annotation'),
    'MORRISON_bulk': dict(
        modality='bulk RNA-seq', modality_class='BULK_RNA',
        scale=ExpressionScale.LOG2_CPM, probe_based=False, binding='NATIVE_ID_HEADER',
        provenance_kind=FILENAME,
        provenance="file name 'RNA-CancerCell-MORRISON1-no_batch_correction-logcpm-all_"
                   "samples.tsv' declares author-computed log2 CPM and explicitly no batch "
                   'correction, so the values pass through unchanged'),
    'MGH_GSE115821_bulk': dict(
        modality='bulk RNA-seq', modality_class='BULK_RNA',
        scale=ExpressionScale.RAW_COUNTS, probe_based=False, binding='MGH_LIBRARY_LABEL',
        provenance_kind=OUR_SCRIPT,
        provenance='per-gene count matrix assembled locally from the MGH GSE115821 HTSeq '
                   'pipeline; the native row values are integer read counts per gene'),
    'MGH_GSE168204_bulk': dict(
        modality='bulk RNA-seq', modality_class='BULK_RNA',
        scale=ExpressionScale.RAW_COUNTS, probe_based=False, binding='MGH_LIBRARY_LABEL',
        provenance_kind=OUR_SCRIPT,
        provenance='per-gene count matrix assembled locally from the MGH GSE168204 HTSeq '
                   'pipeline; the native row values are integer read counts per gene'),
    'GSE139533_bulk': dict(
        modality='bulk RNA-seq', modality_class='BULK_RNA',
        scale=ExpressionScale.RAW_COUNTS, probe_based=False, binding='GEO_SAMPLE_LABEL',
        provenance_kind=FILENAME + ' + ' + GEO,
        provenance="file name 'GSE139533_mergede_protein_coding_counts.tab.gz' declares "
                   'counts; GEO data processing: "read counts matrix produced with '
                   'featureCounts"'),
    'GSE310856_bulk': dict(
        modality='bulk RNA-seq', modality_class='BULK_RNA',
        scale=ExpressionScale.RAW_COUNTS, probe_based=False, binding='GEO_SAMPLE_LABEL',
        provenance_kind=FILENAME + ' + ' + GEO,
        provenance="file name 'GSE310856_raw_counts.txt.gz' declares raw counts; GEO data "
                   'processing: "Bulk RNA-Seq reads were aligned using STAR v2.7.9a to GENCODE '
                   'Human Release 39 (GRCh38.p13) and read counts matrix produced"'),
    'GSE319794_bulk': dict(
        modality='bulk RNA-seq', modality_class='BULK_RNA',
        scale=ExpressionScale.TPM, probe_based=False, binding='NATIVE_ID_HEADER',
        provenance_kind=FILENAME,
        provenance="file name 'GSE319794_NovaSeq36_ProteinCoding_TPM.txt.gz' declares "
                   'transcripts per million'),
    'GSE207422_bulk': dict(
        modality='bulk RNA-seq', modality_class='BULK_RNA',
        scale=ExpressionScale.LOG2_TPM, probe_based=False, binding='NATIVE_ID_HEADER',
        provenance_kind=FILENAME,
        provenance="file name 'GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz' declares log2 TPM, "
                   'so the values are already transformed and must not be logged again'),
    'GSE319641_bulk': dict(
        modality='bulk RNA-seq', modality_class='BULK_RNA',
        scale=ExpressionScale.LINEAR_ABUNDANCE_COMBAT_ADJUSTED, probe_based=False,
        binding='GEO_SAMPLE_LABEL',
        provenance_kind=GEO,
        provenance='GEO data processing (verbatim): "Gene-level transcripts per million (TPM) '
                   'values were obtained using RSEM. TPM values were log2-transformed after '
                   'adding an offset of 1, ComBat-corrected, and converted back to linear TPM." '
                   'The delivered matrix is therefore 2**y - 1 with y = ComBat(log2(TPM + 1)): '
                   'it is not linear TPM any more (column sums no longer equal 1e6 and values '
                   'reach -0.98), so the only legal within-source representation recovers the '
                   'author-computed log scale by log2(x + 1). We perform no batch correction '
                   'here; the matrix is flagged author_batch_corrected because the author '
                   'already did it'),
    # ------------------------------------------------------- microarray, single channel
    'GSE18728_array': dict(
        modality='microarray', modality_class='MICROARRAY',
        scale=ExpressionScale.ARRAY_NORMALIZED_LOG, probe_based=True, binding='GSM_ID_HEADER',
        provenance_kind=GEO,
        provenance='GEO data processing (verbatim): "Probe-level data were normalized and gene '
                   'expression summaries were computed for each probe set using Robust '
                   'Multichip Analysis (RMA)"; RMA output is already log2-scale'),
    'GSE55374_array': dict(
        modality='microarray', modality_class='MICROARRAY',
        scale=ExpressionScale.ARRAY_NORMALIZED_LOG, probe_based=True, binding='GSM_ID_HEADER',
        provenance_kind=GEO,
        provenance='GEO data processing (verbatim): "Illumina probe profiles were quantile '
                   'normalised using the lumi package"; lumi quantile-normalizes on the log2 '
                   'scale'),
    'GSE20181_array': dict(
        modality='microarray', modality_class='MICROARRAY',
        scale=ExpressionScale.ARRAY_NORMALIZED_LINEAR, probe_based=True,
        binding='GSM_ID_HEADER',
        provenance_kind=GEO_SAMPLE,
        provenance='GEO !Sample_data_processing value for all 176 samples is "MAS5"; Affymetrix '
                   'MAS5 reports a linear normalized intensity, not a log2 value, so the '
                   'flagged extension route log2(x + 1) applies and no library-size '
                   'renormalization is legitimate for array intensity'),
    'GSE87455_array': dict(
        modality='microarray', modality_class='MICROARRAY',
        scale=ExpressionScale.ARRAY_NORMALIZED_LINEAR, probe_based=True,
        binding='GSM_ID_HEADER',
        provenance_kind=GEO,
        provenance='GEO data processing (verbatim): "The data were normalised using quantile '
                   'normalisation with BASE"; Illumina BASE writes linear summated normalized '
                   'signal (all 275 columns share one fixed total), so the flagged extension '
                   'route log2(x + 1) applies'),
    # ----------------------------------------- microarray: semantics not established
    'GSE3578_array': dict(
        modality='microarray', modality_class='MICROARRAY',
        scale=ExpressionScale.UNKNOWN, probe_based=True, binding='GSM_ID_HEADER',
        provenance_kind=GEO + ' (insufficient)',
        provenance='GEO data processing states only "global median normalization" for a 2004 '
                   'spotted cDNA/oligo platform, and the delivered values are neither a '
                   'documented linear intensity nor a documented log scale (median 1.008, '
                   '0.47% negative, maximum 1221.5), which is consistent with a ratio-like '
                   'quantity. No legal transform exists without inventing one'),
    'GSE65303_array': dict(
        modality='microarray', modality_class='MICROARRAY',
        scale=ExpressionScale.UNKNOWN, probe_based=True, binding='GSM_ID_HEADER',
        provenance_kind=GEO + ' (insufficient)',
        provenance='GEO data processing names the software used (arrayQualityMetrics, GEOquery, '
                   'Biobase, marray, limma) but never states the normalization or the scale of '
                   'the delivered values; magnitude alone is not an admissible basis, so the '
                   'source fails closed until its submitters document it'),
    'TRIO_US_B07_array': dict(
        modality='two-channel microarray', modality_class='MICROARRAY',
        scale=ExpressionScale.UNKNOWN, probe_based=True, binding='GSM_ID_HEADER',
        provenance_kind=GEO,
        provenance='GEO data processing (verbatim): "Agilent Feature Extraction Software was '
                   'used with no background subtraction and Linear/LOWESS normalization ... dye '
                   'normalization select method = Rank Consistent probes" and "Agilent Feature '
                   'Extracted data was imported into Rosetta Resolver version 7.2 ... Multiple '
                   'probes with the same sequence code were combined into a single value using '
                   'the Rosetta Resolver error-weighted averaging"; the matrix is centred on '
                   'zero (50.5% negative, range -2.85 to 2.40), i.e. a two-channel log-ratio '
                   'quantity, so there is no single-sample abundance to place on a canonical '
                   'abundance scale'),
    # ------------------------------------------------------ scRNA-seq pseudobulk (ours)
    'GSE111014_pseudobulk': dict(
        modality='scRNA-seq pseudobulk', modality_class='SCRNA_PSEUDOBULK',
        scale=ExpressionScale.RAW_COUNTS, probe_based=False, binding='DERIVED_DONOR_DAY',
        provenance_kind=OUR_SCRIPT,
        provenance='written by longitudinal-data/src/cll_sparse_qc.py as per-donor/per-day sums '
                   'of raw per-cell counts, with the barcode set checked against the GEO donor '
                   'and timepoint manifest'),
    'GSE152469_pseudobulk': dict(
        modality='scRNA-seq pseudobulk', modality_class='SCRNA_PSEUDOBULK',
        scale=ExpressionScale.RAW_COUNTS, probe_based=False, binding='NATIVE_ID_HEADER',
        provenance_kind=OUR_SCRIPT,
        provenance='written by our single-cell QC as sums of the deposited UMI count table, one '
                   'column per donor/timepoint pseudobulk'),
    'GSE116256_pseudobulk': dict(
        modality='scRNA-seq pseudobulk', modality_class='SCRNA_PSEUDOBULK',
        scale=ExpressionScale.LOG2_CPM, probe_based=False, binding='NATIVE_ID_HEADER',
        provenance_kind=OUR_SCRIPT,
        provenance='written by longitudinal-data/src/single_cell_qc.py as np.log2(cpm + 1) of '
                   'our own pseudobulk counts, so the values are already log2 CPM'),
    'GSE123813_BCC_pseudobulk': dict(
        modality='scRNA-seq pseudobulk', modality_class='SCRNA_PSEUDOBULK',
        scale=ExpressionScale.LOG2_CPM, probe_based=False, binding='NATIVE_ID_HEADER',
        provenance_kind=OUR_SCRIPT,
        provenance='written by longitudinal-data/src/single_cell_qc.py as np.log2(cpm + 1) of '
                   'our own pseudobulk counts, so the values are already log2 CPM'),
    'GSE123813_SCC_pseudobulk': dict(
        modality='scRNA-seq pseudobulk', modality_class='SCRNA_PSEUDOBULK',
        scale=ExpressionScale.LOG2_CPM, probe_based=False, binding='NATIVE_ID_HEADER',
        provenance_kind=OUR_SCRIPT,
        provenance='written by longitudinal-data/src/single_cell_qc.py as np.log2(cpm + 1) of '
                   'our own pseudobulk counts, so the values are already log2 CPM'),
    'GSE165897_pseudobulk': dict(
        modality='scRNA-seq pseudobulk', modality_class='SCRNA_PSEUDOBULK',
        scale=ExpressionScale.LOG2_CPM, probe_based=False, binding='NATIVE_ID_HEADER',
        provenance_kind=OUR_SCRIPT,
        provenance='written by longitudinal-data/src/single_cell_qc.py as np.log2(cpm + 1) of '
                   'our own pseudobulk counts, so the values are already log2 CPM'),
}

SAMPLE_FIELD = {
    'GSM_ID_HEADER': 'repository_sample_id', 'GEO_SAMPLE_LABEL': 'repository_sample_id',
    'NATIVE_ID_HEADER': 'native_sample_id', 'MGH_LIBRARY_LABEL': 'native_sample_id',
    'DERIVED_DONOR_DAY': 'native_sample_id'}
# GEO fields that name a sample rather than describe its biology. Characteristics are
# deliberately excluded: they carry shared attributes such as "patient id: 151".
GEO_LABEL_FIELDS = ('!Sample_title', '!Sample_description')
# Accession of the locally present series matrix that documents each non-series-matrix source.
DOCUMENTING_SERIES = {
    'GSE139533_bulk': 'GSE139533', 'GSE310856_bulk': 'GSE310856', 'GSE319641_bulk': 'GSE319641',
    'MGH_GSE115821_bulk': 'GSE115821', 'MGH_GSE168204_bulk': 'GSE168204',
    'GSE207422_bulk': 'GSE207422', 'GSE111014_pseudobulk': 'GSE111014',
    'GSE152469_pseudobulk': 'GSE152469', 'GSE116256_pseudobulk': 'GSE116256',
    'GSE123813_BCC_pseudobulk': 'GSE123813', 'GSE123813_SCC_pseudobulk': 'GSE123813',
    'GSE165897_pseudobulk': 'GSE165897'}
METRIC_COLUMNS = [
    'source_expression', 'modality', 'modality_class', 'input_scale', 'probe_based',
    'gene_aggregation', 'transform', 'status', 'provenance_kind', 'normalization_provenance',
    'sample_binding_rule', 'matrix_features', 'matrix_columns', 'samples_declared',
    'samples_bound', 'samples_one_row_several_libraries', 'samples_no_matching_column',
    'matrix_columns_bound', 'matrix_columns_unbound', 'native_features',
    'mapped_native_features', 'features_without_canonical_gene', 'canonical_genes', 'samples',
    'input_nan_fraction', 'min', 'max', 'negative_fraction', 'zero_fraction',
    'integer_like_fraction', 'output_min', 'output_max', 'output_negative_fraction',
    'output_zero_fraction', 'route_check_expected', 'route_check_genes',
    'route_check_max_abs_diff', 'linear_library_sum_min', 'linear_library_sum_max',
    'cpm_invariant_holds', 'author_batch_corrected', 'failure', 'output_relpath',
    'output_size_bytes']


def open_text(path: Path):
    return (pd.io.common.gzip.open(path, 'rt') if str(path).endswith('.gz')
            else open(path, 'r', encoding='utf-8', errors='replace'))


def load_native_matrix(source: str) -> pd.DataFrame:
    """Native numeric matrix exactly as distributed: rows = native features, cols = libraries.

    Values are coerced to numbers with the same rule ``expression_transforms._numeric_frame``
    applies, so an unparseable cell becomes missing rather than silently turning a whole
    column into a non-numeric (object) column.  Nothing is rescaled or reshaped here.
    """
    spec = bg.SOURCES[source]
    path = bg.PROJECT / spec['file']
    if spec['reader'] == 'series_matrix':
        with open_text(path) as handle:
            for line in handle:
                if line.startswith('!series_matrix_table_begin'):
                    header = next(handle).rstrip('\n').split('\t')
                    index, data = [], []
                    for line in handle:
                        if line.startswith('!series_matrix_table_end'):
                            break
                        parts = line.rstrip('\n').split('\t')
                        index.append(parts[0].strip('"').strip())
                        data.append(parts[1:])
        frame = pd.DataFrame(data, index=index, columns=[c.strip('"').strip() for c in header[1:]])
    else:
        sep = ',' if str(path).endswith('.csv.gz') else '\t'
        frame = pd.read_csv(path, sep=sep, index_col=0, dtype=str, encoding='utf-8-sig')
    frame.index = [str(v) for v in frame.index]
    frame.columns = [str(c).strip() for c in frame.columns]
    return frame.apply(pd.to_numeric, errors='coerce')


def geo_sample_labels(source: str) -> dict[str, list[str]]:
    """{GSM: labels GEO itself uses to name that sample} from the documenting series matrix.

    Only ``!Sample_title`` and ``!Sample_description`` rows are read, and every repeated row
    contributes its own value, so a GEO record that names a library in a second description
    row is usable without anyone having to choose which row "means" the sample.
    """
    accession = DOCUMENTING_SERIES.get(source)
    if accession is None:
        return {}
    path = next((candidate for candidate in (
        bg.RAW / accession / f'{accession}_series_matrix.txt.gz',
        bg.PROJECT / 'data' / 'raw' / f'{accession}_series_matrix.txt.gz')
        if candidate.exists()), None)
    if path is None:
        return {}
    gsms: list[str] = []
    label_rows: list[list[str]] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('!series_matrix_table_begin'):
                break
            if not line.startswith('!Sample'):
                continue
            field, *values = [p.strip().strip('"') for p in line.rstrip('\n').split('\t')]
            if field == '!Sample_geo_accession' and not gsms:
                gsms = values
            elif field in GEO_LABEL_FIELDS:
                label_rows.append(values)
    labels: dict[str, list[str]] = {}
    for index, gsm in enumerate(gsms):
        found = [row[index].strip() for row in label_rows if len(row) > index
                 and row[index].strip()]
        prefixed = [strip_library_prefix(value) for value in found]
        labels[gsm] = sorted({v for v in found + prefixed if v})
    return labels


def strip_library_prefix(value: str) -> str:
    """'Library name: TNBC00001_Baseline' -> 'TNBC00001_Baseline'; anything else unchanged."""
    head, separator, tail = value.partition(':')
    if separator and head.strip().lower() == 'library name':
        return tail.strip()
    return value.strip()


def matrix_column_key(rule: str, column: str) -> str:
    if rule == 'MGH_LIBRARY_LABEL':
        return re.sub(r'(?i)[._-]bam$', '', column.strip()).lower().replace('-', '_')
    return column.strip()


def sample_column_keys(rule: str, token: str, labels: dict[str, list[str]]) -> list[str]:
    """Matrix-column keys one corpus identifier token can denote (possibly several, possibly none)."""
    token = token.strip()
    if rule == 'GSM_ID_HEADER' or rule == 'NATIVE_ID_HEADER':
        return [token]
    if rule == 'GEO_SAMPLE_LABEL':
        return labels.get(token, [])
    if rule == 'MGH_LIBRARY_LABEL':
        return [matrix_column_key(rule, token)]
    if rule == 'DERIVED_DONOR_DAY':
        match = re.search(r'from patient (\S+) after (\d+) days', token)
        return [f'{match.group(1)}_D{match.group(2)}'] if match else []
    raise ValueError(f'unsupported sample binding rule: {rule}')


def bind_samples(rule: str, samples: pd.DataFrame, columns: list[str],
                 labels: dict[str, list[str]]) -> pd.DataFrame:
    """Bind each corpus sample of one source to exactly one matrix column.

    A binding is accepted only when it is unique in both directions.  Anything else is
    reported as an explicit non-binding; nothing is resolved by approximate matching.
    """
    by_key: dict[str, list[str]] = {}
    for column in columns:
        by_key.setdefault(matrix_column_key(rule, column), []).append(column)

    rows = []
    for sample_id, raw in zip(samples.sample_id, samples[SAMPLE_FIELD[rule]]):
        tokens = [t.strip() for t in str(raw).split('|') if t.strip()] if pd.notna(raw) else []
        keys = {key: token for token in tokens for key in sample_column_keys(rule, token, labels)}
        candidates = sorted({c for key in keys for c in by_key.get(key, [])})
        if len(candidates) == 1:
            status, column = 'BOUND', candidates[0]
        elif candidates:
            status, column = 'AMBIGUOUS_MULTIPLE_MATRIX_COLUMNS', '|'.join(candidates)
        else:
            status, column = 'NO_UNIQUE_MATRIX_COLUMN', ''
        matched = sorted({key for key in keys if by_key.get(key)})
        rows.append(dict(sample_id=sample_id, native_tokens='|'.join(tokens),
                         matched_label=(matched[0] if len(matched) == 1 else ' / '.join(matched)),
                         matrix_column=column, binding_status=status))
    frame = pd.DataFrame(rows, columns=['sample_id', 'native_tokens', 'matched_label',
                                         'matrix_column', 'binding_status'])
    claims = frame.loc[frame.binding_status.eq('BOUND'), 'matrix_column'].value_counts()
    claimed = set(claims[claims > 1].index)
    if claimed:
        duplicate = frame.matrix_column.isin(claimed) & frame.binding_status.eq('BOUND')
        frame.loc[duplicate, ['binding_status', 'matrix_column']] = [
            'AMBIGUOUS_MATRIX_COLUMN_CLAIMED_TWICE', '']
    frame.insert(0, 'binding_rule', rule)
    frame.insert(0, 'source_expression', '')
    return frame


def contract_for(source: str) -> SourceExpressionContract:
    entry = CONTRACTS[source]
    return SourceExpressionContract(
        source_expression=source, modality=entry['modality'], scale=entry['scale'],
        probe_based=entry['probe_based'], normalization_provenance=entry['provenance'],
        gene_aggregation='SUM_BEFORE_TRANSFORM'
        if entry['scale'] is ExpressionScale.RAW_COUNTS else 'MEDIAN_AFTER_TRANSFORM')


def route_expectation(scale: ExpressionScale) -> str:
    """The declared route as an independent recomputation, not a call into the module."""
    if scale is ExpressionScale.RAW_COUNTS:
        return 'LOG2_OF_CPM_PLUS_1'
    if scale in (ExpressionScale.TPM, ExpressionScale.FPKM,
                 ExpressionScale.ARRAY_NORMALIZED_LINEAR,
                 ExpressionScale.LINEAR_ABUNDANCE_COMBAT_ADJUSTED):
        return 'LOG2_OF_X_PLUS_1'
    return 'IDENTITY'


def verify_route(source: str, native: pd.DataFrame, columns: list[str], feature_map: pd.DataFrame,
                 output: pd.DataFrame) -> dict:
    """Re-derive the declared route from the native matrix and compare it to what was written.

    Only canonical genes served by exactly one mapped native feature are checked, so the
    median/sum collapse cannot mask a wrong transform.  A nonzero difference means a value
    was transformed twice, not at all, or by the wrong rule.
    """
    scale = CONTRACTS[source]['scale']
    mapped = feature_map[feature_map.mapping_status.isin(
        ['EXACT', 'ALIAS', 'PREVIOUS_SYMBOL', 'PLATFORM_ANNOTATION'])]
    mapped = mapped.drop_duplicates('original_feature_id', keep=False)
    per_gene = mapped.hgnc_symbol.value_counts()
    single = [(row.original_feature_id, row.hgnc_symbol) for row in mapped.itertuples()
              if per_gene[row.hgnc_symbol] == 1 and row.original_feature_id in native.index
              and row.hgnc_symbol in output.index]
    if not single:
        return {'route_check_genes': 0, 'route_check_max_abs_diff': None,
                'route_check_expected': route_expectation(scale)}
    features = [f for f, _ in single]
    genes = [g for _, g in single]
    difference = 0.0
    for column in columns[:3]:
        values = native.loc[features, column].to_numpy(dtype=float)
        produced = output.loc[genes, column].to_numpy(dtype=float)
        if scale is ExpressionScale.RAW_COUNTS:
            mapped_rows = mapped.loc[mapped.original_feature_id.isin(native.index),
                                     'original_feature_id']
            total = native.loc[mapped_rows, column].sum()
            expected = np.log2(values / total * 1_000_000.0 + 1.0)
        elif route_expectation(scale) == 'LOG2_OF_X_PLUS_1':
            expected = np.log2(values + 1.0)
        else:
            expected = values
        usable = np.isfinite(expected) & np.isfinite(produced)
        difference = max(difference, float(np.abs(expected[usable] - produced[usable]).max()))
    linear = np.power(2.0, output.to_numpy(dtype=float)) - 1.0
    library = pd.Series(linear.sum(axis=0), index=columns).to_numpy()
    return {'route_check_genes': len(single), 'route_check_max_abs_diff': difference,
            'route_check_expected': route_expectation(scale),
            'linear_library_sum_min': float(np.nanmin(library)),
            'linear_library_sum_max': float(np.nanmax(library)),
            'cpm_invariant_holds': bool(np.allclose(library, 1e6, rtol=1e-6))
            if scale is ExpressionScale.RAW_COUNTS else None}


def build_source(source: str, samples: pd.DataFrame) -> tuple[dict, pd.DataFrame, dict[str, str]]:
    entry = CONTRACTS[source]
    native = load_native_matrix(source)
    cohort_samples = samples[samples.relpath == bg.SOURCES[source]['file']]
    binding = bind_samples(entry['binding'], cohort_samples, list(native.columns),
                           geo_sample_labels(source))
    binding['source_expression'] = source
    bound = binding[binding.binding_status.eq('BOUND')]
    columns = [c for c in native.columns if c in set(bound.matrix_column)]
    record = dict(
        modality=entry['modality'], modality_class=entry['modality_class'],
        input_scale=entry['scale'].value, probe_based=entry['probe_based'],
        gene_aggregation=contract_for(source).gene_aggregation,
        provenance_kind=entry['provenance_kind'],
        normalization_provenance=entry['provenance'], sample_binding_rule=entry['binding'],
        matrix_features=int(native.shape[0]), matrix_columns=int(native.shape[1]),
        samples_declared=int(binding.shape[0]), samples_bound=int(len(bound)),
        matrix_columns_bound=len(columns),
        matrix_columns_unbound=int(native.shape[1] - len(columns)),
        native_features=int(native.shape[0]), author_batch_corrected=False,
        failure='', transform='NOT_APPLIED', canonical_genes=0, samples=0,
        mapped_native_features=0, features_without_canonical_gene=int(native.shape[0]),
        min=None, max=None, negative_fraction=None, zero_fraction=None,
        integer_like_fraction=None, input_nan_fraction=None, output_min=None, output_max=None,
        output_negative_fraction=None, output_zero_fraction=None,
        route_check_expected=None, route_check_genes=None, route_check_max_abs_diff=None,
        linear_library_sum_min=None, linear_library_sum_max=None, cpm_invariant_holds=None,
        output_relpath='', output_size_bytes=0)
    columns_by_sample = dict(zip(bound.matrix_column, bound.sample_id))
    if entry['scale'] is ExpressionScale.UNKNOWN:
        record['status'] = 'SEMANTICS_NOT_ESTABLISHED'
        record['failure'] = entry['provenance']
        return record, binding, {}
    if not columns:
        record['status'] = 'SAMPLE_BINDING_FAILED'
        record['failure'] = (f'no sample of {source} binds uniquely to a matrix column '
                             f'under rule {entry["binding"]}')
        return record, binding, {}
    try:
        feature_map = pd.read_csv(bg.FEATURE_MAP / f'{source}.feature_map.csv.gz',
                                  dtype={'original_feature_id': str}, low_memory=False)
        output, report = build_canonical_quantitative_matrix(
            native[columns], feature_map, contract_for(source))
    except ValueError as error:
        record['status'] = 'SCALE_CONTRACT_FAILED'
        record['failure'] = str(error)
        return record, binding, {}
    if output.shape[1] != len(columns):
        record['status'] = 'SCALE_CONTRACT_FAILED'
        record['failure'] = f'column count changed during aggregation: {output.shape[1]}'
        return record, binding, {}
    produced = output.to_numpy(dtype=float).ravel()
    produced = produced[np.isfinite(produced)]
    record.update(input_nan_fraction=float(native[columns].isna().to_numpy().mean()),
                  output_min=float(produced.min()), output_max=float(produced.max()),
                  output_negative_fraction=float(np.mean(produced < 0)),
                  output_zero_fraction=float(np.mean(produced == 0)))
    path = MATRIX_DIR / f'{source}.canonical_log_expression.parquet'
    output.to_parquet(path)
    record.update(report)
    verification = verify_route(source, native[columns], columns, feature_map, output)
    record.update(verification)
    diverged = (verification['route_check_max_abs_diff'] is not None
                and verification['route_check_max_abs_diff'] > 1e-9)
    record.update(
        status='ROUTE_CHECK_FAILED' if diverged else 'QUANTITATIVE_READY',
        failure=(f'written values differ from the declared {verification["route_check_expected"]} '
                 f'route by up to {verification["route_check_max_abs_diff"]}'
                 if diverged else record['failure']),
        features_without_canonical_gene=int(native.shape[0] - report['mapped_native_features']),
        author_batch_corrected=bool(report.get('author_batch_corrected', False)),
        output_relpath=path.relative_to(bg.PROJECT).as_posix(),
        output_size_bytes=int(path.stat().st_size))
    record.pop('source_expression', None)
    present = set(output.columns)
    return record, binding, {(source, columns_by_sample[column]): column
                             for column in columns if column in present}


def main() -> int:
    for directory in (LAYER, MANIFESTS, MATRIX_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    symmetric_difference = set(bg.SOURCES) ^ set(CONTRACTS)
    if symmetric_difference:
        raise SystemExit(f'contract/source inventory mismatch: {sorted(symmetric_difference)}')

    samples = pd.read_parquet(bg.RELEASE / 'samples.parquet')
    samples['relpath'] = samples.expression_file.map(
        lambda v: (bg.PROJECT / str(v)).relative_to(bg.PROJECT).as_posix()
        if pd.notna(v) and str(v) not in ('', 'NOT_AVAILABLE') else '')
    master = pd.read_parquet(bg.RELEASE / 'PAIR_LONGITUDINAL_MASTER.parquet')[
        ['pair_uid', 'sample_t0', 'sample_t1']]

    records, bindings, pair_columns = [], [], {}
    for source in sorted(bg.SOURCES):
        record, binding, columns_by_sample = build_source(source, samples)
        record['source_expression'] = source
        records.append(record)
        bindings.append(binding)
        pair_columns.update(columns_by_sample)

    metric_frame = pd.DataFrame(records)
    for column in METRIC_COLUMNS:
        if column not in metric_frame:
            metric_frame[column] = None
    binding_frame = pd.concat(bindings, ignore_index=True)[
        ['source_expression', 'sample_id', 'binding_rule', 'native_tokens', 'matched_label',
         'matrix_column', 'binding_status']].sort_values(['source_expression', 'sample_id'])
    binding_frame.to_csv(MANIFESTS / 'SAMPLE_COLUMN_BINDING.csv.gz', index=False,
                         compression='gzip')
    metric_frame[METRIC_COLUMNS].to_csv(LAYER / 'SOURCE_EXPRESSION_METRICS.csv', index=False)
    metric_frame[['source_expression', 'modality', 'modality_class', 'input_scale', 'probe_based',
                  'gene_aggregation', 'transform', 'sample_binding_rule', 'provenance_kind',
                  'normalization_provenance', 'status']].to_csv(
        LAYER / 'SOURCE_EXPRESSION_CONTRACT.csv', index=False)

    pairs = pd.read_csv(bg.GENE_SPACE / 'PAIR_GENE_SPACE.csv.gz').merge(
        master, on='pair_uid', how='left')
    pairs = pairs[pairs.gene_space_status != 'EXPRESSION_FILE_NOT_LOCAL'].copy()
    ready_sources = set(metric_frame.loc[
        metric_frame.status.eq('QUANTITATIVE_READY'), 'source_expression'])

    def classify(row):
        if pd.isna(row.source_expression):
            return 'SOURCE_NOT_RESOLVED_IN_GENE_SPACE', str(row.gene_space_status)
        if row.source_expression not in ready_sources:
            status = metric_frame.loc[
                metric_frame.source_expression.eq(row.source_expression), 'status']
            return str(status.iloc[0]) if len(status) else 'SOURCE_NOT_DECLARED', ''
        t0 = pair_columns.get((row.source_expression, row.sample_t0), '')
        t1 = pair_columns.get((row.source_expression, row.sample_t1), '')
        if t0 and t1:
            return ('QUANTITATIVE_READY', 'SAME_MATRIX_COLUMN_FOR_BOTH_ENDPOINTS'
                    if t0 == t1 else '')
        missing = [end for end, column in (('t0', t0), ('t1', t1)) if not column]
        return ('ENDPOINT_SAMPLE_NOT_BOUND_TO_MATRIX_COLUMN', 'unbound=' + ','.join(missing))

    classified = pairs.apply(classify, axis=1, result_type='expand')
    pairs['quantitative_status'], pairs['quantitative_detail'] = classified[0], classified[1]
    for end in ('t0', 't1'):
        pairs[f'matrix_column_{end}'] = pairs.apply(
            lambda r, end=end: pair_columns.get((r.source_expression, r[f'sample_{end}']), ''),
            axis=1)
    coverage = pairs[['pair_uid', 'patient_uid', 'cohort_code', 'source_expression',
                      'expression_modality', 'modality_class', 'platform', 'sample_t0',
                      'matrix_column_t0', 'sample_t1', 'matrix_column_t1',
                      'quantitative_status', 'quantitative_detail', 'gene_space_status',
                      'clinical_endpoint_available', 'strict_prcr_vs_pd_eligible',
                      'frozen_auo_eligible']].copy()
    coverage.to_csv(LAYER / 'PAIR_EXPRESSION_COVERAGE.csv.gz', index=False, compression='gzip')

    usable = coverage.quantitative_status.eq('QUANTITATIVE_READY') & (
        coverage.matrix_column_t0 != coverage.matrix_column_t1)
    ready = metric_frame[metric_frame.status.eq('QUANTITATIVE_READY')]
    total_bytes = int(metric_frame.output_size_bytes.sum())
    summary = dict(
        expression_layer_version='v0.3-C1',
        level='C1 technology-native quantitative expression layer',
        built_from='dataset/releases/v0.1.2 (MASTER + samples) + gene_space v0.2.1 feature maps '
                   '+ locally present native expression files',
        mathematical_authority='dataset/src/expression_transforms.py: the supplied Level C1 '
                               'contract, extended by two explicitly flagged additive scales',
        headline=dict(
            sources_total=int(len(metric_frame)),
            sources_quantitative_ready=int(len(ready)),
            sources_semantics_not_established=sorted(metric_frame.loc[
                metric_frame.status.eq('SEMANTICS_NOT_ESTABLISHED'), 'source_expression']),
            sources_scale_contract_failed=sorted(metric_frame.loc[
                metric_frame.status.eq('SCALE_CONTRACT_FAILED'), 'source_expression']),
            sources_sample_binding_failed=sorted(metric_frame.loc[
                metric_frame.status.eq('SAMPLE_BINDING_FAILED'), 'source_expression']),
            sources_route_check_failed=sorted(metric_frame.loc[
                metric_frame.status.eq('ROUTE_CHECK_FAILED'), 'source_expression']),
            samples_declared=int(metric_frame.samples_declared.sum()),
            samples_quantitative_ready=int(metric_frame.loc[
                metric_frame.status.eq('QUANTITATIVE_READY'), 'samples_bound'].sum()),
            samples_bound_in_unready_sources=int(metric_frame.loc[
                metric_frame.status.ne('QUANTITATIVE_READY'), 'samples_bound'].sum()),
            pairs_with_t0_t1_expression=int(len(pairs)),
            pairs_in_canonical_gene_space=int(pairs.in_canonical_gene_space.sum()),
            pairs_quantitative_ready=int(usable.sum()),
            patients_quantitative_ready=int(coverage.loc[usable, 'patient_uid'].nunique()),
            matrix_files=int((metric_frame.output_relpath != '').sum()),
            local_derived_matrix_total_bytes=total_bytes,
            local_derived_matrix_total_mebibytes=round(total_bytes / 1048576.0, 1)),
        modality_pair_classes={
            name: dict(pairs=int((coverage.modality_class.eq(name)).sum()),
                       pairs_quantitative_ready=int((
                           coverage.modality_class.eq(name) & usable).sum()),
                       patients_quantitative_ready=int(coverage.loc[
                           coverage.modality_class.eq(name) & usable, 'patient_uid'].nunique()))
            for name in sorted(coverage.modality_class.unique())},
        per_source=metric_frame[METRIC_COLUMNS].to_dict('records'),
        scale_vocabulary=dict(
            supplied=['RAW_COUNTS', 'TPM', 'FPKM', 'LOG2_CPM', 'LOG2_TPM', 'LOG2_FPKM',
                      'LOG_EXPRESSION', 'ARRAY_NORMALIZED_LOG', 'UNKNOWN'],
            added_by_the_implementation_engineer={
                'ARRAY_NORMALIZED_LINEAR':
                    'linear array intensity already normalized by the platform procedure (MAS5, '
                    'Illumina BASE quantile): validated non-negative, then log2(x + 1), then the '
                    'median of the probes resolving to one canonical gene. No library-size '
                    'renormalization, because array intensity is not counts.',
                'LINEAR_ABUNDANCE_COMBAT_ADJUSTED':
                    'linear-scale abundance the author already batch-corrected, with a documented '
                    'chain of 2 ** y - 1 where y = ComBat(log2(TPM + 1)). Validated against the '
                    '-1 floor, then log2(x + 1) recovers y exactly. Flagged '
                    'author_batch_corrected=true so it is never pooled with uncorrected sources.'}),
        transformation_verification=dict(
            method='for every ready source the declared route is recomputed from the native '
                   'matrix and compared to the written matrix on canonical genes served by '
                   'exactly one mapped native feature, so neither the median nor the sum '
                   'collapse can hide a wrong transform',
            max_abs_difference_over_all_sources=float(
                metric_frame.route_check_max_abs_diff.dropna().max()),
            sources_where_difference_is_zero=int((
                metric_frame.route_check_max_abs_diff.dropna() == 0.0).sum()),
            cpm_invariant='every RAW_COUNTS source must have each column of the written '
                          'matrix sum to exactly 1e6 after 2**x - 1; no other source may, '
                          'because a source that was not ours to renormalize would look as if '
                          'it had been',
            cpm_invariant_holds_for_counts_sources=int((
                metric_frame.cpm_invariant_holds.dropna()).sum()),
            cpm_invariant_expected=int((metric_frame.input_scale.eq('RAW_COUNTS')
                                        & metric_frame.status.eq('QUANTITATIVE_READY')).sum())),
        adjudication=dict(
            admissible_provenance=[FILENAME, GEO, GEO_SAMPLE, OUR_SCRIPT],
            forbidden=['inferring a scale from numeric magnitude',
                       'logarithming an already-transformed matrix a second time',
                       'renormalizing an already-normalized matrix',
                       'changing a declared scale after validation fails'],
            double_log_guards=[
                'MORRISON_bulk, GSE207422_bulk, GSE116256_pseudobulk, GSE123813_BCC_pseudobulk, '
                'GSE123813_SCC_pseudobulk and GSE165897_pseudobulk are declared already-'
                'transformed and take the identity route',
                'GIDE_bulk fails closed instead of being logged on an undocumented normalization',
                'GSE319641_bulk is logged exactly once, to undo the documented back-transform, '
                'and never re-normalized'],
            counts_route='native mapped rows -> sum duplicate features per canonical gene -> '
                         'CPM -> log2(CPM + 1); never CPM per feature then average',
            probe_route='probes are never summed; the median of the probes resolving to one '
                        'canonical gene is taken after the transform'),
        redistribution=dict(
            derived_matrices_committed=False,
            reason='the repository carries no per-accession licence ledger, GEO supplementary '
                   'files are distributed under their submitters terms, and a derived matrix '
                   'inherits them; the Owner must authorize redistribution and review the size '
                   'before any matrix is pushed',
            local_only_path=MATRIX_DIR.relative_to(bg.PROJECT).as_posix(),
            files=[dict(source_expression=r.source_expression, output_relpath=r.output_relpath,
                        size_bytes=int(r.output_size_bytes))
                   for r in metric_frame[metric_frame.output_relpath.ne('')].itertuples()]),
        explicit_non_goals=[
            'no cross-cohort or cross-platform correction of any kind: no ComBat, no joint '
            'quantile normalization, no global z-score, no mean centring',
            'no ranks, no paired deltas, no feature selection, no models, no endpoint prediction',
            'the matrices are NOT on one common numerical scale: a value in one source is '
            'comparable only to other values in that same source',
            'the canonical gene space (which feature) is resolved in v0.2.1; this layer resolves '
            'what the numbers mean inside a source, and is not itself cross-platform comparability'])
    (LAYER / 'EXPRESSION_LAYER_REPORT.json').write_text(json.dumps(summary, indent=2),
                                                        encoding='utf-8')

    print(json.dumps(dict(headline=summary['headline'],
                          modality_pair_classes=summary['modality_pair_classes']), indent=2))
    print(metric_frame[['source_expression', 'input_scale', 'transform', 'status',
                        'samples_bound', 'canonical_genes', 'min', 'max']].to_string(index=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
