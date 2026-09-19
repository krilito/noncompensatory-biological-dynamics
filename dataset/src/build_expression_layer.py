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
    'canonical_count_mass_fraction_min', 'canonical_count_mass_fraction_median',
    'canonical_count_mass_fraction_max', 'native_count_mass_median',
    'library_denominator_max_rel_dev', 'canonical_cpm_mass_identity_holds',
    'preprocessing_scope', 'fold_isolation',
    'author_batch_corrected', 'failure', 'output_relpath', 'output_size_bytes']

# Downstream modelling provenance.  Every transform this layer applies is computed from one
# sample column alone and fits nothing across samples, so no train/test fold can contaminate
# it.  A source the author already batch-corrected does not have that property: its scale was
# produced using the whole source, and a model must not present it as fold-internal
# preprocessing.  GSE319641 (NeoTRIP) is that case; see the Owner's ratification.
PREPROCESSING_PROVENANCE = {
    'GSE319641_bulk': ('AUTHOR_FULL_SOURCE', 'NOT_ESTABLISHED')}
DEFAULT_PREPROCESSING_PROVENANCE = ('OUR_C1_PER_SAMPLE', 'FOLD_INDEPENDENT')

# Adjudication of every AMBIGUOUS_MULTIPLE_MATRIX_COLUMNS case: one corpus sample row that
# names several matrix columns.  ``replicate_type`` is decided from the documentary evidence
# quoted below, never from label resemblance, and no case is proven to be a technical
# sequencing library of one biospecimen, so nothing is collapsed in C1.
MULTI_COLUMN_ADJUDICATION = {
    'GSE139533_bulk': dict(
        replicate_type='BIOLOGICAL_REGION_REPLICATE',
        evidence='GEO !Series_summary (verbatim): "We investigated longitudinal transcriptomic '
                 'patterns associated with glioblastoma (GB) recurrence by integrative approach '
                 'utilizing multisampling strategy ... In total, 128 tissue samples of 44 tumors '
                 '... were analyzed". Each candidate column is a separate GEO sample with its own '
                 'BioSample accession and its own !Sample_title tissue index (e.g. GSM4143223 '
                 '"_677_1", GSM4143224 "_677_2", GSM4143225 "_677_3" of tumor 677), i.e. multiple '
                 'tissue pieces of one patient, which is the object the study measures. Some rows '
                 'additionally span two specimen codes of one patient (e.g. "383_591_3" with '
                 '"383_650_1"), which are anatomically distinct lesions'),
    'GSE3578_array': dict(
        replicate_type='UNRESOLVED',
        evidence='the corpus row pools two GEO samples whose !Sample_title values differ only by '
                 'a trailing "_1"/"_2" (e.g. "1581_p_1"/"1581_p_2"); GEO documents no aliquot or '
                 'library relationship between them, and this source fails closed on expression '
                 'scale anyway, so the case cannot and need not be resolved in C1'),
    'GSE91061_bulk': dict(
        replicate_type='UNRESOLVED',
        evidence='the corpus row pools GSM2420319/GSM2420320, titled "Pt109_On_AE527955-6" and '
                 '"Pt109_On_AE527955-5"; the suffix is not documented as either an aliquot or a '
                 'sequencing library, and the source scale is author FPKM, for which no generic '
                 'sum or mean rule exists'),
    'MGH_GSE115821_bulk': dict(
        replicate_type='UNRESOLVED',
        evidence='the candidate columns come from two different GEO platforms of this series '
                 '(GPL11154, 23 RNA-Seq samples, library names such as "MGH208_031115-1.bam", and '
                 'GPL18573, 14 RNA-Seq samples, library names such as "208-3-11-15_S13"), and the '
                 'corpus records the pooled row as replicate_class '
                 'MULTIPLE_ASSAYS_OR_SAME_VISIT_SPECIMENS_NOT_EXTRA_TIME with quality flag '
                 'MULTIOMIC_OR_SAME_VISIT. No locally present document states whether same-visit '
                 'libraries such as "MGH39_082514-1" to "-5" are aliquots of one specimen or '
                 'separate pieces'),
    'MGH_GSE168204_bulk': dict(
        replicate_type='UNRESOLVED',
        evidence='same situation as GSE115821 for this donor series: libraries "MGHIPIPD1001_070715'
                 '-1/-2/-3" share one visit date, the corpus row is flagged '
                 'MULTIOMIC_OR_SAME_VISIT, and no local document states which biospecimen each '
                 'library was made from'),
}
REPLICATE_TYPE_VOCABULARY = [
    'TECHNICAL_LIBRARY_REPLICATE', 'BIOLOGICAL_REGION_REPLICATE',
    'BIOLOGICAL_ALIQUOT_REPLICATE', 'DISTINCT_BIOSPECIMENS', 'UNRESOLVED']
# What C1 is allowed to do once a case is typed.  Nothing here collapses a row: a biological
# replicate must stay a separate measurement unit, and an unresolved one stays unbound.
RECOMMENDED_ACTION = {
    'BIOLOGICAL_REGION_REPLICATE': 'PRESERVE_SEPARATELY_DO_NOT_COLLAPSE_IN_C1',
    'BIOLOGICAL_ALIQUOT_REPLICATE': 'PRESERVE_SEPARATELY_DO_NOT_COLLAPSE_IN_C1',
    'DISTINCT_BIOSPECIMENS': 'PRESERVE_SEPARATELY_DO_NOT_COLLAPSE_IN_C1',
    'UNRESOLVED': 'REMAIN_UNBOUND'}
ADJUDICATION_COLUMNS = [
    'source_expression', 'sample_id', 'patient_uid', 'timepoint_native', 'native_tokens',
    'geo_labels', 'candidate_matrix_columns', 'n_candidate_columns', 'input_scale',
    'corpus_replicate_class', 'biospecimen_identity_evidence', 'replicate_type',
    'recommended_action']


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


def adjudicate_multi_column(binding: pd.DataFrame, samples: pd.DataFrame,
                            master: pd.DataFrame) -> pd.DataFrame:
    """Type every corpus sample row that names several matrix columns, with its evidence.

    The corpus sample is coarser than the matrix column: one row can name several libraries
    because they are technical replicates, because a tumor was multisampled, or because one
    visit produced several specimens.  These are not interchangeable, and summing or averaging
    the wrong kind erases intratumoral heterogeneity by hand.  C1 therefore binds nothing here
    and publishes the decision for the Owner instead, typed from documentary evidence.
    """
    ambiguous = binding[binding.binding_status.eq('AMBIGUOUS_MULTIPLE_MATRIX_COLUMNS')]
    if ambiguous.empty:
        return pd.DataFrame(columns=ADJUDICATION_COLUMNS)
    patient_of: dict[str, str] = {}
    for row in master.itertuples():
        for sample in (row.sample_t0, row.sample_t1):
            patient_of[sample] = row.patient_uid
    corpus = samples.set_index('sample_id')
    rows = []
    for case in ambiguous.itertuples():
        verdict = MULTI_COLUMN_ADJUDICATION.get(case.source_expression)
        if verdict is None:
            raise ValueError(
                f'{case.source_expression}: multi-column case has no adjudication; type it '
                'from documentary evidence before rebuilding')
        replicate_type = verdict['replicate_type']
        if replicate_type not in RECOMMENDED_ACTION:
            raise ValueError(
                f"{case.source_expression}: replicate_type {replicate_type} is a technical "
                'library replicate, which no C1 collapse rule handles. Summing native library '
                'columns before CPM is legal only for RAW_COUNTS and must be implemented '
                'against a per-case decision, not against a vocabulary member.')
        metadata = corpus.loc[case.sample_id]
        rows.append(dict(
            source_expression=case.source_expression, sample_id=case.sample_id,
            patient_uid=patient_of.get(case.sample_id, ''),
            timepoint_native=str(metadata.timepoint_native),
            native_tokens=case.native_tokens, geo_labels=case.matched_label,
            candidate_matrix_columns=case.matrix_column,
            n_candidate_columns=len(str(case.matrix_column).split('|')),
            input_scale=CONTRACTS[case.source_expression]['scale'].value,
            corpus_replicate_class=str(metadata.replicate_class),
            biospecimen_identity_evidence=verdict['evidence'],
            replicate_type=replicate_type,
            recommended_action=RECOMMENDED_ACTION[replicate_type]))
    return pd.DataFrame(rows, columns=ADJUDICATION_COLUMNS).sort_values(
        ['source_expression', 'sample_id'])


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
        return 'LOG2_OF_NATIVE_CPM_PLUS_1'
    if scale in (ExpressionScale.TPM, ExpressionScale.FPKM,
                 ExpressionScale.ARRAY_NORMALIZED_LINEAR,
                 ExpressionScale.LINEAR_ABUNDANCE_COMBAT_ADJUSTED):
        return 'LOG2_OF_X_PLUS_1'
    return 'IDENTITY'


def verify_route_implementation(source: str, native: pd.DataFrame, columns: list[str],
                                feature_map: pd.DataFrame, output: pd.DataFrame) -> dict:
    """Re-derive the declared route from the native matrix and compare it to what was written.

    Only canonical genes served by exactly one mapped native feature are checked, so the
    median/sum collapse cannot mask a wrong transform.  A nonzero difference means a value
    was transformed twice, not at all, or by the wrong rule.  For RAW_COUNTS the denominator
    used here is the whole native column, which is what the route declares; this check says
    nothing about whether that route is scientifically right, which is what
    ``verify_library_denominator`` is for.
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
            expected = np.log2(values / native[column].sum(skipna=True) * 1_000_000.0 + 1.0)
        elif route_expectation(scale) == 'LOG2_OF_X_PLUS_1':
            expected = np.log2(values + 1.0)
        else:
            expected = values
        usable = np.isfinite(expected) & np.isfinite(produced)
        difference = max(difference, float(np.abs(expected[usable] - produced[usable]).max()))
    # Descriptive for every scale: what 2**x - 1 of the written matrix sums to per sample.  A
    # source we were not entitled to renormalize must not land on 1e6, which is how a second
    # library-size normalization would announce itself.
    linear = np.power(2.0, output.to_numpy(dtype=float)) - 1.0
    library = pd.Series(linear.sum(axis=0), index=output.columns).to_numpy(dtype=float)
    return {'route_check_genes': len(single), 'route_check_max_abs_diff': difference,
            'route_check_expected': route_expectation(scale),
            'linear_library_sum_min': float(np.nanmin(library)),
            'linear_library_sum_max': float(np.nanmax(library))}


def verify_library_denominator(output: pd.DataFrame, native_count_mass: pd.Series,
                               mapped_count_mass: pd.Series) -> dict:
    """LIBRARY_DENOMINATOR_CHECK: the retained genes must carry their true native share.

    The v0.3 invariant ``sum(2 ** output - 1) == 1e6`` is deleted: it held by construction
    for a denominator built from the mapped features only, which silently renormalized the
    canonical genes.  What must hold instead, per sample, is

        sum(2 ** output - 1) == 1e6 * mapped_count_mass / native_count_mass

    where ``native_count_mass`` counts every native feature, including those that never
    entered the canonical gene space.  A source whose unmapped features carry no counts
    legitimately reaches 1e6; one whose unmapped features carry a fifth of the reads must
    land near 200k, and any deviation from the identity below means we renormalized a
    quantity we were not entitled to renormalize.
    """
    linear = np.power(2.0, output.to_numpy(dtype=float)) - 1.0
    observed = pd.Series(linear.sum(axis=0), index=output.columns).to_numpy(dtype=float)
    expected = 1_000_000.0 * mapped_count_mass.to_numpy(dtype=float) \
        / native_count_mass.to_numpy(dtype=float)
    deviation = np.abs(observed - expected) / np.maximum(np.abs(expected), 1.0)
    fraction = mapped_count_mass / native_count_mass
    return {'native_count_mass_median': float(native_count_mass.median()),
            'canonical_count_mass_fraction_min': float(fraction.min()),
            'canonical_count_mass_fraction_median': float(fraction.median()),
            'canonical_count_mass_fraction_max': float(fraction.max()),
            'library_denominator_max_rel_dev': float(deviation.max()),
            'canonical_cpm_mass_identity_holds': bool(np.allclose(
                observed, expected, rtol=1e-6, atol=1e-6))}


def build_source(source: str, samples: pd.DataFrame) -> tuple[
        dict, pd.DataFrame, dict[str, str], list[dict]]:
    """Build one source's C1 matrix, its sample-binding rows, its column map and its
    per-sample native library sizes.

    The fourth element is empty for every non-RAW_COUNTS source: only a count matrix has a
    native count mass to account for.
    """
    entry = CONTRACTS[source]
    native = load_native_matrix(source)
    cohort_samples = samples[samples.relpath == bg.SOURCES[source]['file']]
    binding = bind_samples(entry['binding'], cohort_samples, list(native.columns),
                           geo_sample_labels(source))
    binding['source_expression'] = source
    bound = binding[binding.binding_status.eq('BOUND')]
    columns = [c for c in native.columns if c in set(bound.matrix_column)]
    preprocessing_scope, fold_isolation = PREPROCESSING_PROVENANCE.get(
        source, DEFAULT_PREPROCESSING_PROVENANCE)
    record = dict(
        modality=entry['modality'], modality_class=entry['modality_class'],
        input_scale=entry['scale'].value, probe_based=entry['probe_based'],
        gene_aggregation=contract_for(source).gene_aggregation,
        provenance_kind=entry['provenance_kind'],
        normalization_provenance=entry['provenance'], sample_binding_rule=entry['binding'],
        matrix_features=int(native.shape[0]), matrix_columns=int(native.shape[1]),
        samples_declared=int(binding.shape[0]), samples_bound=int(len(bound)),
        samples_one_row_several_libraries=int(
            binding.binding_status.eq('AMBIGUOUS_MULTIPLE_MATRIX_COLUMNS').sum()),
        matrix_columns_bound=len(columns),
        matrix_columns_unbound=int(native.shape[1] - len(columns)),
        native_features=int(native.shape[0]), author_batch_corrected=False,
        preprocessing_scope=preprocessing_scope, fold_isolation=fold_isolation,
        failure='', transform='NOT_APPLIED', canonical_genes=0, samples=0,
        mapped_native_features=0, features_without_canonical_gene=int(native.shape[0]),
        min=None, max=None, negative_fraction=None, zero_fraction=None,
        integer_like_fraction=None, input_nan_fraction=None, output_min=None, output_max=None,
        output_negative_fraction=None, output_zero_fraction=None,
        route_check_expected=None, route_check_genes=None, route_check_max_abs_diff=None,
        linear_library_sum_min=None, linear_library_sum_max=None,
        canonical_count_mass_fraction_min=None, canonical_count_mass_fraction_median=None,
        canonical_count_mass_fraction_max=None, native_count_mass_median=None,
        library_denominator_max_rel_dev=None, canonical_cpm_mass_identity_holds=None,
        output_relpath='', output_size_bytes=0)
    columns_by_sample = dict(zip(bound.matrix_column, bound.sample_id))
    if entry['scale'] is ExpressionScale.UNKNOWN:
        record['status'] = 'SEMANTICS_NOT_ESTABLISHED'
        record['failure'] = entry['provenance']
        return record, binding, {}, []
    if not columns:
        record['status'] = 'SAMPLE_BINDING_FAILED'
        record['failure'] = (f'no sample of {source} binds uniquely to a matrix column '
                             f'under rule {entry["binding"]}')
        return record, binding, {}, []
    try:
        feature_map = pd.read_csv(bg.FEATURE_MAP / f'{source}.feature_map.csv.gz',
                                  dtype={'original_feature_id': str}, low_memory=False)
        output, report = build_canonical_quantitative_matrix(
            native[columns], feature_map, contract_for(source))
    except ValueError as error:
        record['status'] = 'SCALE_CONTRACT_FAILED'
        record['failure'] = str(error)
        return record, binding, {}, []
    if output.shape[1] != len(columns):
        record['status'] = 'SCALE_CONTRACT_FAILED'
        record['failure'] = f'column count changed during aggregation: {output.shape[1]}'
        return record, binding, {}, []
    denominators = []
    count_mass = None
    if 'native_count_mass' in report:
        native_count_mass = report.pop('native_count_mass')
        mapped_count_mass = report.pop('mapped_count_mass')
        count_mass = (native_count_mass, mapped_count_mass)
        denominators = [dict(
            source_expression=source, sample_id=columns_by_sample.get(column, ''),
            matrix_column=column, input_scale=entry['scale'].value,
            native_count_mass=float(native_count_mass[column]),
            mapped_count_mass=float(mapped_count_mass[column]),
            canonical_count_mass_fraction=float(mapped_count_mass[column]
                                                / native_count_mass[column]))
            for column in output.columns]
    produced = output.to_numpy(dtype=float).ravel()
    produced = produced[np.isfinite(produced)]
    record.update(input_nan_fraction=float(native[columns].isna().to_numpy().mean()),
                  output_min=float(produced.min()), output_max=float(produced.max()),
                  output_negative_fraction=float(np.mean(produced < 0)),
                  output_zero_fraction=float(np.mean(produced == 0)))
    path = MATRIX_DIR / f'{source}.canonical_log_expression.parquet'
    output.to_parquet(path)
    record.update(report)
    verification = verify_route_implementation(
        source, native[columns], columns, feature_map, output)
    record.update(verification)
    if count_mass is not None:
        record.update(verify_library_denominator(output, *count_mass))
    diverged = (verification['route_check_max_abs_diff'] is not None
                and verification['route_check_max_abs_diff'] > 1e-9)
    denominator_broken = record['canonical_cpm_mass_identity_holds'] is False
    record.update(
        status=('ROUTE_CHECK_FAILED' if diverged else
                'LIBRARY_DENOMINATOR_FAILED' if denominator_broken else 'QUANTITATIVE_READY'),
        failure=(f'written values differ from the declared {verification["route_check_expected"]} '
                 f'route by up to {verification["route_check_max_abs_diff"]}'
                 if diverged else
                 f'canonical CPM mass deviates from 1e6 x mapped/native count mass by up to '
                 f'{record["library_denominator_max_rel_dev"]:.3e}'
                 if denominator_broken else record['failure']),
        features_without_canonical_gene=int(native.shape[0] - report['mapped_native_features']),
        author_batch_corrected=bool(report.get('author_batch_corrected', False)),
        output_relpath=path.relative_to(bg.PROJECT).as_posix(),
        output_size_bytes=int(path.stat().st_size))
    record.pop('source_expression', None)
    present = set(output.columns)
    return record, binding, {(source, columns_by_sample[column]): column
                             for column in columns if column in present}, denominators


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
        ['pair_uid', 'patient_uid', 'sample_t0', 'sample_t1']]
    records, bindings, pair_columns, denominators = [], [], {}, []
    for source in sorted(bg.SOURCES):
        record, binding, columns_by_sample, library_sizes = build_source(source, samples)
        record['source_expression'] = source
        records.append(record)
        bindings.append(binding)
        denominators.extend(library_sizes)
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
    pd.DataFrame(denominators).sort_values(['source_expression', 'matrix_column']).to_csv(
        MANIFESTS / 'RAW_COUNT_LIBRARY_DENOMINATORS.csv.gz', index=False, compression='gzip')
    adjudication = adjudicate_multi_column(binding_frame, samples, master)
    adjudication.to_csv(MANIFESTS / 'SAMPLE_MULTI_COLUMN_ADJUDICATION.csv.gz', index=False,
                        compression='gzip')
    metric_frame[METRIC_COLUMNS].to_csv(LAYER / 'SOURCE_EXPRESSION_METRICS.csv', index=False)
    metric_frame[['source_expression', 'modality', 'modality_class', 'input_scale', 'probe_based',
                  'gene_aggregation', 'transform', 'sample_binding_rule', 'provenance_kind',
                  'normalization_provenance', 'status']].to_csv(
        LAYER / 'SOURCE_EXPRESSION_CONTRACT.csv', index=False)

    pairs = pd.read_csv(bg.GENE_SPACE / 'PAIR_GENE_SPACE.csv.gz').merge(
        master[['pair_uid', 'sample_t0', 'sample_t1']], on='pair_uid', how='left')
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
        expression_layer_version='v0.3.1-C1',
        level='C1 technology-native quantitative expression layer',
        built_from='dataset/releases/v0.1.2 (MASTER + samples) + gene_space v0.2.1 feature maps '
                   '+ locally present native expression files',
        repair_of='v0.3-C1: the RAW_COUNTS CPM denominator was the count mass of the mapped '
                  'features only, which renormalized the retained canonical genes to 1e6 per '
                  'column; v0.3.1 takes the denominator from the complete native matrix',
        mathematical_authority='dataset/src/expression_transforms.py: the supplied Level C1 '
                               'contract, extended by two additive scales that the Owner '
                               'ratified for v0.3.1',
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
            sources_library_denominator_failed=sorted(metric_frame.loc[
                metric_frame.status.eq('LIBRARY_DENOMINATOR_FAILED'), 'source_expression']),
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
            added_by_the_implementation_engineer_and_ratified_by_the_owner={
                'ARRAY_NORMALIZED_LINEAR':
                    'linear array intensity already normalized by the platform procedure (MAS5, '
                    'Illumina BASE quantile): validated non-negative, then log2(x + 1), then the '
                    'median of duplicate mapped native features per canonical gene. No '
                    'library-size renormalization, because array intensity is not counts.',
                'LINEAR_ABUNDANCE_COMBAT_ADJUSTED':
                    'linear-scale abundance the author already batch-corrected, with a documented '
                    'chain of 2 ** y - 1 where y = ComBat(log2(TPM + 1)). Validated against the '
                    '-1 floor, then log2(x + 1) recovers y exactly: this is not a second '
                    'normalization by us, it is the recovery of the log-space the author analysed. '
                    'Carries author_batch_corrected=true plus the downstream provenance flags '
                    'below, and must not be modelled as an uncorrected matrix.'}),
        route_implementation_check=dict(
            method='for every ready source the declared route is recomputed from the native '
                   'matrix and compared to the written matrix on canonical genes served by '
                   'exactly one mapped native feature, so neither the median nor the sum '
                   'collapse can hide a wrong transform',
            max_abs_difference_over_all_sources=float(
                metric_frame.route_check_max_abs_diff.dropna().max()),
            sources_where_difference_is_zero=int((
                metric_frame.route_check_max_abs_diff.dropna() == 0.0).sum()),
            what_it_does_not_prove='that the declared route is scientifically correct. In v0.3 '
                                   'every counts source passed this check with a difference of '
                                   'exactly zero while the denominator itself was wrong, which '
                                   'is why the denominator is now checked separately below.'),
        library_denominator_check=dict(
            invariant='for every RAW_COUNTS source and every sample column: '
                      'sum(2 ** output - 1) == 1e6 * mapped_count_mass / native_count_mass, '
                      'where native_count_mass is the count mass of the complete native matrix '
                      'and mapped_count_mass that of the features that entered the canonical '
                      'gene space',
            deleted_invariant='sum(2 ** output - 1) == 1e6. It held by construction under the '
                              'v0.3 denominator and proved only that the code followed the '
                              'declared route; it concealed the fact that the route renormalized '
                              'the retained genes and discarded the reads of every unmapped '
                              'native feature',
            max_relative_deviation_over_counts_sources=float(
                metric_frame.library_denominator_max_rel_dev.dropna().max()),
            counts_sources_where_identity_holds=int((
                metric_frame.canonical_cpm_mass_identity_holds.dropna()).sum()),
            counts_sources_checked=int(
                metric_frame.library_denominator_max_rel_dev.notna().sum()),
            non_counts_ready_sources_that_look_per_million=int((
                metric_frame.status.eq('QUANTITATIVE_READY')
                & metric_frame.input_scale.ne('RAW_COUNTS')
                & metric_frame.linear_library_sum_min.between(
                    1e6 - 1.0, 1e6 + 1.0)).sum()),
            per_sample_values='manifests/RAW_COUNT_LIBRARY_DENOMINATORS.csv.gz (one row per '
                              'RAW_COUNTS sample column: native_count_mass, mapped_count_mass, '
                              'canonical_count_mass_fraction)'),
        canonical_count_mass_coverage=dict(
            meaning='how much of a library was carried by features that map into the HGNC '
                    'canonical gene space. This is a different quantity from feature mapping '
                    'coverage, which counts feature IDs: a source can lose most of its IDs to a '
                    'foreign namespace and still keep almost all of its reads, or lose few IDs '
                    'and lose most of the depth. The two fractions are reported side by side per '
                    'source below so the difference is visible; no expected value is hard-coded.',
            per_source=[dict(
                source_expression=r.source_expression,
                feature_mapping_fraction=round(
                    r.mapped_native_features / r.native_features, 6),
                canonical_count_mass_fraction_median=round(
                    r.canonical_count_mass_fraction_median, 6),
                canonical_count_mass_fraction_min=round(
                    r.canonical_count_mass_fraction_min, 6),
                canonical_count_mass_fraction_max=round(
                    r.canonical_count_mass_fraction_max, 6),
                native_count_mass_median=round(r.native_count_mass_median, 1),
                canonical_linear_output_sum_range=[round(r.linear_library_sum_min, 3),
                                                   round(r.linear_library_sum_max, 3)])
                for r in ready[ready.input_scale.eq('RAW_COUNTS')].itertuples()]),
        multi_column_sample_binding=dict(
            question='a corpus sample row can name several matrix columns. That is not evidence '
                     'of technical replication: GSE139533 studies intratumoral heterogeneity '
                     'with an explicit multisampling strategy (128 tissue samples from 44 '
                     'tumors), so pooling several columns and averaging would erase the '
                     'structure the source exists to measure',
            rule='C1 collapses nothing. Only same-specimen technical libraries of a RAW_COUNTS '
                 'source may be summed, and then before CPM; region, aliquot and distinct-'
                 'specimen cases stay separate measurement units, and an unprovable case stays '
                 'unbound. For an already-normalized matrix no generic mean or sum rule is '
                 'invented at all',
            table='manifests/SAMPLE_MULTI_COLUMN_ADJUDICATION.csv.gz',
            replicate_type_vocabulary=REPLICATE_TYPE_VOCABULARY,
            cases=int(adjudication.shape[0]),
            by_replicate_type={key: int(value) for key, value in
                               adjudication.replicate_type.value_counts().items()},
            by_recommended_action={key: int(value) for key, value in
                                   adjudication.recommended_action.value_counts().items()},
            by_source={key: int(value) for key, value in
                       adjudication.source_expression.value_counts().items()},
            technical_library_replicate_cases=int(
                adjudication.replicate_type.eq('TECHNICAL_LIBRARY_REPLICATE').sum()),
            pairs_lost_to_multi_column_binding=int(coverage.quantitative_status.eq(
                'ENDPOINT_SAMPLE_NOT_BOUND_TO_MATRIX_COLUMN').sum())),
        downstream_modelling_provenance=dict(
            preprocessing_scope={
                'OUR_C1_PER_SAMPLE': 'the value is a per-column function of that column alone; '
                                     'nothing was fitted across samples, cohorts or folds, so a '
                                     'train/test split cannot leak through it',
                'AUTHOR_FULL_SOURCE': 'the scale of this source was produced by the authors using '
                                      'the whole source (batch correction across its libraries), '
                                      'so it is not a preprocessing we can refit inside a fold'},
            fold_isolation={
                'FOLD_INDEPENDENT': 'no fitted parameters exist, so hold-out is exact',
                'NOT_ESTABLISHED': 'the author-side correction cannot be re-derived without the '
                                   'other samples, so a strict model must treat it as external, '
                                   'already-corrected data and never as fold-internal '
                                   'preprocessing'},
            author_batch_corrected_sources=sorted(metric_frame.loc[
                metric_frame.author_batch_corrected, 'source_expression']),
            note='a source flagged AUTHOR_FULL_SOURCE may stay in the dataset as a real source '
                 'representation; it must not be presented as an uncorrected matrix. If '
                 'pre-ComBat values are found later, the primary model analysis should prefer '
                 'them'),
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
            counts_route='complete native count matrix -> native library size (denominator); the '
                         'same matrix -> mapped features -> sum duplicate features per canonical '
                         'gene -> divide by the native library size -> CPM -> log2(CPM + 1). The '
                         'retained genes are never renormalized to 1e6, and CPM per feature '
                         'followed by an average is never used',
            probe_route='probes are never summed; a gene takes the median of duplicate mapped '
                        'native features per canonical gene after the transform',
            multi_library_rule='one corpus sample row naming several matrix columns is not '
                               'evidence of technical replication; see '
                               'multi_column_sample_binding'),
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
            'no collapsing of several matrix columns into one sample: C1 binds one corpus sample '
            'to one column and leaves every multi-column case unbound and adjudicated',
            'the matrices are NOT on one common numerical scale: a value in one source is '
            'comparable only to other values in that same source',
            'the canonical gene space (which feature) is resolved in v0.2.1; this layer resolves '
            'what the numbers mean inside a source, and is not itself cross-platform comparability'])
    (LAYER / 'EXPRESSION_LAYER_REPORT.json').write_text(json.dumps(summary, indent=2),
                                                        encoding='utf-8')

    def without(keys, mapping):
        return {key: value for key, value in mapping.items() if key not in keys}

    print(json.dumps(dict(
        headline=summary['headline'],
        modality_pair_classes=summary['modality_pair_classes'],
        route_implementation_check=without(('method',), summary['route_implementation_check']),
        library_denominator_check=without(
            ('invariant', 'deleted_invariant', 'per_sample_values'),
            summary['library_denominator_check']),
        canonical_count_mass_coverage=summary['canonical_count_mass_coverage']['per_source'],
        multi_column_sample_binding=without(
            ('question', 'rule', 'table', 'replicate_type_vocabulary'),
            summary['multi_column_sample_binding'])), indent=2))
    print(metric_frame[['source_expression', 'input_scale', 'transform', 'status',
                        'samples_bound', 'canonical_genes', 'canonical_count_mass_fraction_median',
                        'min', 'max']].to_string(index=False))
    print(adjudication.groupby(['source_expression', 'replicate_type',
                                'recommended_action']).size().to_string())
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
