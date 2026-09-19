"""Build the PAIR v0.2 canonical gene space: per-source feature maps and coverage.

Scope is deliberately limited: identify every locally present expression source, map
every feature of every source onto the HGNC approved-gene coordinate system while
preserving the native representation, and report how many samples and longitudinal
pairs can enter a shared gene space. No normalization, no batch correction, no
ranks, no deltas, no modelling.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd

DATASET_ROOT = Path(__file__).resolve().parents[1]
PROJECT = DATASET_ROOT.parents[1]
RELEASE = DATASET_ROOT / 'releases' / 'v0.1.2'
GENE_SPACE = DATASET_ROOT / 'gene_space'
FEATURE_MAP = GENE_SPACE / 'feature_map'
HGNC_PATH = GENE_SPACE / 'reference' / 'hgnc_complete_set.txt'
HGNC_ORIGIN = 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
RAW = PROJECT / 'longitudinal-data' / 'data' / 'raw'

EXACT, ALIAS, PREVIOUS = 'EXACT', 'ALIAS', 'PREVIOUS_SYMBOL'
PLATFORM_ANNOTATION, AMBIGUOUS, UNMAPPED = 'PLATFORM_ANNOTATION', 'AMBIGUOUS', 'UNMAPPED'
STATUSES = (EXACT, ALIAS, PREVIOUS, PLATFORM_ANNOTATION, AMBIGUOUS, UNMAPPED)
MAPPED_STATUSES = (EXACT, ALIAS, PREVIOUS, PLATFORM_ANNOTATION)

GENE_SYMBOL, ENSEMBL_GENE_ID, ENTREZ_GENE_ID = 'GENE_SYMBOL', 'ENSEMBL_GENE_ID', 'ENTREZ_GENE_ID'
TRANSCRIPT_ID = 'TRANSCRIPT_ID_NOT_GENE_ID'

# One entry per distinct locally present expression file. A cohort with two expression
# files (MGH: GSE115821 and GSE168204) is two sources, and main() fails closed if any
# locally present file in samples.parquet is left undeclared here.
SOURCES = {
    'GIDE_bulk': dict(
        file='data/raw/PRJEB23709/cancercell_normalized_counts_genenames.txt',
        reader='id_and_symbol', feature_type=ENSEMBL_GENE_ID),
    'GSE111014_pseudobulk': dict(
        file='longitudinal-data/data/processed/GSE111014_pseudobulk_counts.csv.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'GSE116256_pseudobulk': dict(
        file='longitudinal-data/data/processed/GSE116256_gene_log_expression.csv.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'GSE123813_BCC_pseudobulk': dict(
        file='longitudinal-data/data/processed/GSE123813_BCC_gene_log_expression.csv.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'GSE123813_SCC_pseudobulk': dict(
        file='longitudinal-data/data/processed/GSE123813_SCC_gene_log_expression.csv.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'GSE139533_bulk': dict(
        file='longitudinal-data/data/raw/GSE139533/GSE139533_mergede_protein_coding_counts.tab.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'GSE152469_pseudobulk': dict(
        file='longitudinal-data/data/processed/GSE152469_pseudobulk_counts.csv.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'GSE165897_pseudobulk': dict(
        file='longitudinal-data/data/processed/GSE165897_gene_log_expression.csv.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'GSE18728_array': dict(
        file='longitudinal-data/data/raw/GSE18728/GSE18728_series_matrix.txt.gz',
        reader='series_matrix', feature_type='AFFYMETRIX_PROBE_ID', platform='GPL570'),
    'GSE20181_array': dict(
        file='longitudinal-data/data/raw/GSE20181/GSE20181_series_matrix.txt.gz',
        reader='series_matrix', feature_type='AFFYMETRIX_PROBE_ID', platform='GPL96'),
    'GSE207422_bulk': dict(
        file='longitudinal-data/data/raw/GSE207422/GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'GSE310856_bulk': dict(
        file='longitudinal-data/data/raw/GSE310856/GSE310856_raw_counts.txt.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'GSE319641_bulk': dict(
        file='longitudinal-data/data/raw/GSE319641/GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'GSE319794_bulk': dict(
        file='longitudinal-data/data/raw/GSE319794/GSE319794_NovaSeq36_ProteinCoding_TPM.txt.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'GSE3578_array': dict(
        file='longitudinal-data/data/raw/GSE3578/GSE3578_series_matrix.txt.gz',
        reader='series_matrix', feature_type='SPOTTED_CDNA_PROBE_ID', platform='GPL2895'),
    'GSE55374_array': dict(
        file='longitudinal-data/data/raw/GSE55374/GSE55374_series_matrix.txt.gz',
        reader='series_matrix', feature_type='ILLUMINA_PROBE_ID', platform='GPL10558'),
    'GSE65303_array': dict(
        file='longitudinal-data/data/raw/GSE65303/GSE65303_series_matrix.txt.gz',
        reader='series_matrix', feature_type='AGILENT_PROBE_ID', platform='GPL16876'),
    'GSE87455_array': dict(
        file='longitudinal-data/data/raw/GSE87455/GSE87455_series_matrix.txt.gz',
        reader='series_matrix', feature_type='ILLUMINA_PROBE_ID', platform='GPL10558'),
    'GSE91061_bulk': dict(
        file='data/raw/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz',
        reader='first_column', feature_type=ENTREZ_GENE_ID),
    'MGH_GSE115821_bulk': dict(
        file='data/raw/GSE115821_MGH_counts.csv.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'MGH_GSE168204_bulk': dict(
        file='data/raw/GSE168204_MGH_counts.csv.gz',
        reader='first_column', feature_type=GENE_SYMBOL),
    'MORRISON_bulk': dict(
        file='data/raw/MORRISON-1-public/RNASeq/data/RNA-CancerCell-MORRISON1-no_batch_correction-logcpm-all_samples.tsv',
        reader='first_column', feature_type=GENE_SYMBOL),
    'TRIO_US_B07_array': dict(
        file='longitudinal-data/data/raw/GSE130786/GSE130786_series_matrix.txt.gz',
        reader='series_matrix', feature_type='AGILENT_PROBE_ID', platform='GPL6480'),
}

# GEO platform annotation files and the columns that carry a gene assignment.
PLATFORM_ANNOTATIONS = {
    'GPL570': ('GPL570.annot.gz', 'Gene symbol', 'Gene ID', None),
    'GPL96': ('GPL96.annot.gz', 'Gene symbol', 'Gene ID', None),
    'GPL10558': ('GPL10558.annot.gz', 'Gene symbol', 'Gene ID', None),
    'GPL2895': ('GPL2895.annot.gz', 'Gene symbol', 'Gene ID', None),
    'GPL6480': ('GPL6480.annot.gz', 'Gene symbol', 'Gene ID', None),
    'GPL16876': ('GPL16876_self.soft', 'GeneSymbol', 'EntrezGeneID', 'EnsemblID'),
}
NULL_TOKENS = ('', '--', '---', 'na', 'n/a', 'unknown', 'uncategorized')

# Token shapes that explain an unmapped feature without needing a crosswalk judgement.
# A class may leave the gene-addressable denominator only when its identifiers are proven
# to live outside the HGNC gene namespace; entries forced to stay addressable denote an
# intended gene (corruption, retired alias, transcript) and count against the source.
#    (class name, pattern, always_gene_addressable, published description)
IDENTIFIER_CLASSES = (
    ('legacy_noncoding_assembly_id', re.compile(r'^NONHSAG\d+$'), False,
     'NONHSAG plus digits: a legacy non-coding assembly gene naming scheme'),
    ('sequence_record_accession',
     re.compile(r'^(?:A[BCFYZ]|B[KQ]|CR|D[QRS]|EU|GC|J[QR]|KU|Z\d)[0-9]{5,7}(?:\.[0-9]+)?$'), False,
     'a GenBank-style sequence record accession with an optional version suffix'),
    ('clone_based_temporary_name',
     re.compile(r'^(?:[A-Z]{1,4}[0-9]{1,2}(?:-[A-Z0-9]+)*|FO[0-9]{6})\.[0-9]{1,2}$'), False,
     'a clone-based temporary gene name carrying a dotted version suffix'),
    ('retired_locuslink_id', re.compile(r'^LOC\d+$'), True,
     'LOC plus digits: a retired NCBI locus-link identifier that denotes a gene'),
    ('spreadsheet_corrupted_symbol',
     re.compile(r'^[0-9]{1,2}-(?:jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)$', re.IGNORECASE),
     True, 'a month-name token produced by spreadsheet date coercion of a gene symbol'),
    ('transcript_identifier', re.compile(r'^(?:ENST|ENSP|ENSD|ENSO)\d{11}'), True,
     'an Ensembl transcript or protein accession rather than a gene identifier'),
    ('unannotated_array_feature', re.compile(r'^.+$'), True,
     'an array feature whose platform annotation carries no gene assignment'),
)
NON_ADDRESSABLE_NAMESPACE_PRESENCE = 0.01  # class must be <1% present in the HGNC namespace
NON_ADDRESSABLE_MIN_TOKENS = 25
GENE_IDENTIFIER_CLASS = 'GENE_IDENTIFIER'
ENTRY_THRESHOLD = 0.5  # coverage floor for a source to enter the core gene space

FEATURE_COLUMNS = ['cohort_code', 'source_expression', 'platform', 'original_feature_id',
                   'original_feature_type', 'original_feature_annotation', 'feature_identifier_class',
                   'gene_addressable', 'canonical_gene_id',
                   'hgnc_symbol', 'ensembl_gene_id', 'entrez_gene_id', 'mapping_method',
                   'mapping_status', 'mapping_ambiguity', 'reference_version']
METRIC_COLUMNS = ['source_expression', 'status', 'cohort_codes', 'modalities', 'platform',
                  'expression_file', 'declared_feature_type', 'local_samples',
                  'raw_feature_count', 'raw_mapped_feature_count', 'raw_mapping_fraction',
                  'gene_addressable_feature_count', 'gene_addressable_mapped_count',
                  'gene_addressable_mapping_fraction', 'non_addressable_feature_count',
                  'non_addressable_classes', 'raw_coverage_passes_entry_rule',
                  'gene_addressable_coverage_passes_entry_rule',
                  'addressable_admission_evidenced',
                  'unique_canonical_genes',
                  'exact_mapping_fraction', 'alias_fraction', 'previous_symbol_fraction',
                  'platform_annotation_fraction', 'ambiguous_fraction', 'unmapped_fraction',
                  'n_exact', 'n_alias', 'n_previous_symbol', 'n_platform_annotation',
                  'n_ambiguous', 'n_unmapped', 'duplicate_features', 'features_per_canonical_gene',
                  'enters_core_gene_space', 'unmapped_methods']

ENSEMBL_GENE_RE = re.compile(r'^ENSG\d{11}(\.\d+)?$')
ENSEMBL_TRANSCRIPT_RE = re.compile(r'^(ENST|ENSP|ENSD|ENSO)\d{11}(\.\d+)?$')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relpath(value):
    """Repository-relative POSIX path: the only path form allowed in a public artifact."""
    path = Path(str(value))
    if not path.is_absolute():
        return path.as_posix()
    try:
        return path.resolve().relative_to(PROJECT).as_posix()
    except (ValueError, OSError):
        return path.as_posix()


def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8', errors='replace') \
        if str(path).endswith('.gz') else open(path, 'r', encoding='utf-8', errors='replace')


def read_series_matrix_features(path):
    features = []
    with open_text(path) as fh:
        in_table = False
        for line in fh:
            if line.startswith('!series_matrix_table_begin'):
                in_table = True
                next(fh)  # sample-column header row
                continue
            if line.startswith('!series_matrix_table_end'):
                break
            if in_table:
                features.append(line.split('\t')[0].strip('"').strip())
    return features


def read_first_column(path):
    sep = ',' if str(path).endswith('.csv.gz') else '\t'
    frame = pd.read_csv(path, sep=sep, usecols=[0], dtype=str, encoding='utf-8-sig')
    return [(value, '') for value in frame.iloc[:, 0].fillna('').tolist()]


def read_id_and_symbol(path):
    """GIDE layout: stable Ensembl gene id column plus an author gene-name column."""
    frame = pd.read_csv(path, sep='\t', usecols=[0, 1], dtype=str)
    return list(zip(frame.iloc[:, 0].fillna(''), frame.iloc[:, 1].fillna('')))


READERS = {'series_matrix': lambda p: [(f, '') for f in read_series_matrix_features(p)],
           'first_column': read_first_column, 'id_and_symbol': read_id_and_symbol}


def read_platform_annotation(platform):
    """{probe_id: {symbol: [...], entrez: [...], ensembl: [...]}} from the GPL table."""
    filename, symbol_column, entrez_column, ensembl_column = PLATFORM_ANNOTATIONS[platform]
    path = RAW / platform / filename
    if not path.exists():
        raise FileNotFoundError(f'no local platform annotation for {platform}')
    with open_text(path) as fh:
        header = None
        for line in fh:
            if line.startswith('!platform_table_begin'):
                header = next(fh).rstrip('\n').split('\t')
                break
        if header is None:
            raise RuntimeError(f'no platform table in {path}')
        index = {name: (header.index(name) if name in header else None)
                 for name in ('ID', symbol_column, entrez_column, ensembl_column)}
        if index['ID'] is None:
            raise RuntimeError(f'{path} has no ID column')
        wanted = [('symbol', index[symbol_column]), ('entrez', index[entrez_column]),
                  ('ensembl', index[ensembl_column])]
        stop = max([index['ID']] + [i for _, i in wanted if i is not None])
        mapping = {}
        for line in fh:
            if line.startswith('!platform_table_end'):
                break
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= stop:
                continue

            def tokens(column):
                if column is None:
                    return []
                return sorted({t.strip() for t in parts[column].split('///')
                               if t.strip().lower() not in NULL_TOKENS})

            mapping[parts[index['ID']].strip('"')] = {
                name: tokens(column) for name, column in wanted}
    return mapping, path


class HgncReference:
    """HGNC approved-gene table plus the only crosswalks the mapping is allowed to use."""

    def __init__(self):
        frame = pd.read_csv(HGNC_PATH, sep='\t', low_memory=False, dtype=str).fillna('')
        approved = frame.loc[frame.status.eq('Approved')]
        self.table = {row.symbol: dict(gene_id=row.hgnc_id, ensembl=row.ensembl_gene_id,
                                       entrez=row.entrez_id, locus_group=row.locus_group)
                      for row in approved.itertuples(index=False)}
        self.prev, self.alias, self.entrez, self.ensembl = {}, {}, {}, {}
        for row in approved.itertuples(index=False):
            for field, store in (('prev_symbol', self.prev), ('alias_symbol', self.alias)):
                for value in str(getattr(row, field)).split('|'):
                    if value.strip():
                        store.setdefault(value.strip(), set()).add(row.symbol)
            entrez = str(row.entrez_id).strip().split('.')[0]
            if entrez.isdigit():
                self.entrez.setdefault(entrez, set()).add(row.symbol)
            for value in str(row.ensembl_gene_id).split(';'):
                if value.strip():
                    self.ensembl.setdefault(value.strip(), set()).add(row.symbol)
        self.digest = sha256(HGNC_PATH)
        self.namespace = set(self.table) | set(self.prev) | set(self.alias)
        self.version = f'HGNC_complete_set@{self.digest[:12]}'
        self.n_rows = len(approved)
        self.n_protein_coding = sum(1 for entry in self.table.values()
                                    if entry['locus_group'] == 'protein-coding gene')

    def fields(self, symbol):
        entry = self.table.get(symbol, {})
        entrez = str(entry.get('entrez', '')).strip().split('.')[0]
        return entry.get('gene_id', ''), entry.get('ensembl', ''), entrez if entrez.isdigit() else ''

    def by_symbol(self, token):
        """(candidates, status, route) for one native gene-symbol token."""
        token = str(token).strip()
        if not token:
            return set(), '', 'NO_TOKEN'
        if token in self.table:
            return {token}, EXACT, 'HGNC_APPROVED_SYMBOL'
        if token in self.prev:
            return set(self.prev[token]), PREVIOUS, 'HGNC_PREV_SYMBOL'
        if token in self.alias:
            return set(self.alias[token]), ALIAS, 'HGNC_ALIAS_SYMBOL'
        return set(), '', 'NOT_IN_HGNC'

    def by_ensembl(self, token):
        return set(self.ensembl.get(str(token).strip().split('.')[0], set()))

    def by_entrez(self, token):
        value = str(token).strip().split('.')[0]
        return set(self.entrez[value]) if value.isdigit() and value in self.entrez else set()


def native_feature_type(declared, feature_id):
    """Identifier shape is decisive for Ensembl ids; otherwise keep the source declaration."""
    if ENSEMBL_TRANSCRIPT_RE.match(feature_id):
        return TRANSCRIPT_ID
    if declared == GENE_SYMBOL and ENSEMBL_GENE_RE.match(feature_id):
        return ENSEMBL_GENE_ID
    return declared


FORCED_ADDRESSABLE = {name: forced for name, _, forced, _ in IDENTIFIER_CLASSES}


def identifier_class(token):
    for name, pattern, _, _ in IDENTIFIER_CLASSES[:-1]:
        if pattern.match(str(token)):
            return name
    return GENE_IDENTIFIER_CLASS


def classify_features(frame, spec):
    """Identifier class of each feature plus the token set that evidence came from.

    A probe carries no identifier claim of its own, so its class is read from the gene
    assignment its platform annotation supplies; a probe with no assignment at all is an
    unannotated array feature, which still counts against the source.
    """
    probe = bool(spec.get('platform'))
    classes, tokens_by_class = [], {}
    for feature_id, annotation in zip(frame.original_feature_id, frame.original_feature_annotation):
        annotation = '' if pd.isna(annotation) else str(annotation)
        if probe and not annotation:
            classes.append('unannotated_array_feature')
            tokens_by_class.setdefault('unannotated_array_feature', set()).add(str(feature_id))
            continue
        if probe:
            tokens = [t for t in annotation.split('|') if t and not t.isdigit()]
        else:
            tokens = [str(feature_id)]
        found = {identifier_class(token) for token in tokens}
        name = found.pop() if len(found) == 1 else GENE_IDENTIFIER_CLASS
        classes.append(name)
        if name != GENE_IDENTIFIER_CLASS:
            tokens_by_class.setdefault(name, set()).update(tokens)
    return pd.Series(classes, index=frame.index), tokens_by_class


def resolve_probe(entry, ref, platform):
    """(candidates, status, method) for one probe annotation record."""
    for column, route in (('symbol', 'GENE_SYMBOL'), ('entrez', 'ENTREZ_GENE_ID'),
                          ('ensembl', 'ENSEMBL_GENE_ID')):
        found, route_status = set(), ''
        for token in entry[column]:
            if column == 'symbol':
                candidates, route_status, _ = ref.by_symbol(token)
            elif column == 'entrez':
                candidates = ref.by_entrez(token)
            else:
                candidates = ref.by_ensembl(token)
            found |= candidates
            if len(found) > 1:
                break
        if len(found) == 1:
            method = f'PLATFORM_ANNOTATION:{platform}:{route}'
            if column == 'symbol' and route_status and route_status != EXACT:
                method += f'({route_status.lower()})'
            return found, PLATFORM_ANNOTATION, method
        if len(found) > 1:
            return found, AMBIGUOUS, f'PLATFORM_ANNOTATION:{platform}:{route}'
    return set(), UNMAPPED, f'PLATFORM_ANNOTATION:{platform}:ANNOTATION_GENE_ID_NOT_IN_HGNC'


def map_features(key, spec, ref, annotations):
    """One row per native feature of one expression source."""
    path = PROJECT / spec['file']
    platform = spec.get('platform', '')
    if platform and platform not in annotations:
        annotations[platform] = read_platform_annotation(platform)
    table = annotations[platform][0] if platform else {}
    rows = []
    for raw_id, author_symbol in READERS[spec['reader']](path):
        feature_id = str(raw_id).strip()
        ftype = native_feature_type(spec['feature_type'], feature_id)
        candidates, status, method, native = set(), UNMAPPED, '', ''

        if platform:
            record = table.get(feature_id)
            if record is None:
                method = f'PLATFORM_ANNOTATION:{platform}:PROBE_NOT_IN_TABLE'
            elif not (record['symbol'] or record['entrez'] or record['ensembl']):
                method = f'PLATFORM_ANNOTATION:{platform}:ANNOTATION_WITHOUT_GENE_ASSIGNMENT'
            else:
                native = '|'.join(record['symbol'] + record['entrez'] + record['ensembl'])
                candidates, status, method = resolve_probe(record, ref, platform)
        elif ftype == TRANSCRIPT_ID:
            method = 'TRANSCRIPT_ID_NOT_MAPPED_TO_GENE'
        elif ftype == ENSEMBL_GENE_ID:
            found = ref.by_ensembl(feature_id)
            if found:
                candidates, status = found, EXACT if len(found) == 1 else AMBIGUOUS
                method = 'ENSEMBL_GENE_ID'
            elif author_symbol:
                candidates, status, method = ref.by_symbol(author_symbol)
                method = f'AUTHOR_SYMBOL_FALLBACK:{method}'
            else:
                method = 'ENSEMBL_GENE_ID_NOT_IN_HGNC'
        elif ftype == ENTREZ_GENE_ID:
            found = ref.by_entrez(feature_id)
            if found:
                candidates, status = found, EXACT if len(found) == 1 else AMBIGUOUS
                method = 'ENTREZ_GENE_ID'
            else:
                method = 'ENTREZ_GENE_ID_NOT_IN_HGNC'
        else:
            candidates, status, method = ref.by_symbol(feature_id)

        if len(candidates) == 1 and status in MAPPED_STATUSES:
            gene = next(iter(candidates))
            gene_id, ensembl, entrez = ref.fields(gene)
            rows.append((feature_id, ftype, native, gene_id, gene, ensembl, entrez,
                         method, status, ''))
        elif len(candidates) > 1:
            rows.append((feature_id, ftype, native, '', '', '', '', method, AMBIGUOUS,
                         '|'.join(sorted(candidates))))
        else:
            rows.append((feature_id, ftype, native, '', '', '', '', method, UNMAPPED, ''))
    frame = pd.DataFrame(rows, columns=['original_feature_id', 'original_feature_type',
                                        'original_feature_annotation', 'canonical_gene_id',
                                        'hgnc_symbol', 'ensembl_gene_id', 'entrez_gene_id',
                                        'mapping_method', 'mapping_status', 'mapping_ambiguity'])
    classes, tokens_by_class = classify_features(frame, spec)
    frame['feature_identifier_class'] = classes
    frame['gene_addressable'] = True
    frame.insert(0, 'cohort_code', '')
    frame.insert(1, 'source_expression', key)
    frame.insert(2, 'platform', platform or '')
    version = ref.version
    if platform:
        annot = RAW / platform / PLATFORM_ANNOTATIONS[platform][0]
        version += f'+{platform}@{sha256(annot)[:12]}'
    frame['reference_version'] = version
    return frame[FEATURE_COLUMNS], tokens_by_class


def main():
    ref = HgncReference()
    samples = pd.read_parquet(RELEASE / 'samples.parquet')
    master = pd.read_parquet(RELEASE / 'PAIR_LONGITUDINAL_MASTER.parquet')
    samples['expression_relpath'] = samples.expression_file.map(
        lambda v: relpath(v) if v and str(v) != 'NOT_FOUND' else '')

    # Fail closed: every locally present expression file must be a declared source.
    present = {p for p in samples.expression_relpath if p and (PROJECT / p).exists()}
    declared = {spec['file'] for spec in SOURCES.values()}
    if not present <= declared:
        raise RuntimeError(f'undeclared local expression sources: {sorted(present - declared)}')
    if not declared <= present:
        raise RuntimeError(f'declared sources without a local file: {sorted(declared - present)}')

    members = {}
    for row in samples.itertuples(index=False):
        if row.expression_relpath:
            slot = members.setdefault(row.expression_relpath, dict(cohorts=set(), modalities=set(),
                                                                   platforms=set(), samples=0))
            slot['cohorts'].add(row.cohort_code)
            slot['modalities'].add(row.modality or '')
            slot['platforms'].add(row.platform or '')
            slot['samples'] += 1

    FEATURE_MAP.mkdir(parents=True, exist_ok=True)
    annotations, gene_sets, metrics, unmapped_tokens = {}, {}, {}, {}
    for key, spec in sorted(SOURCES.items()):
        slot = members.get(spec['file'], dict(cohorts=set(), modalities=set(), platforms=set(),
                                              samples=0))
        unit = dict(source_expression=key, status='MAPPED',
                    cohort_codes='|'.join(sorted(slot['cohorts'])),
                    modalities='|'.join(sorted(slot['modalities'])),
                    platform=spec.get('platform') or '|'.join(sorted(slot['platforms'])),
                    expression_file=spec['file'], declared_feature_type=spec['feature_type'],
                    local_samples=slot['samples'])
        frame, tokens_by_class = map_features(key, spec, ref, annotations)
        frame['cohort_code'] = unit['cohort_codes']
        frame['platform'] = unit['platform']
        # A class leaves the gene-addressable denominator only on namespace evidence.
        evidence, excluded = {}, set()
        for name, tokens in tokens_by_class.items():
            present_in_hgnc = len(tokens & ref.namespace)
            outside = (not FORCED_ADDRESSABLE[name]
                       and len(tokens) >= NON_ADDRESSABLE_MIN_TOKENS
                       and present_in_hgnc / len(tokens) <= NON_ADDRESSABLE_NAMESPACE_PRESENCE)
            evidence[name] = dict(features=int(frame.feature_identifier_class.eq(name).sum()),
                                  distinct_tokens=len(tokens), in_hgnc_namespace=present_in_hgnc,
                                  namespace_presence_rate=round(present_in_hgnc / len(tokens), 5),
                                  outside_hgnc_gene_namespace=outside)
            if outside:
                excluded.add(name)
        frame['gene_addressable'] = ~frame.feature_identifier_class.isin(excluded)
        frame.to_csv(FEATURE_MAP / f'{key}.feature_map.csv.gz', index=False,
                     compression={'method': 'gzip', 'mtime': 0})
        counts = Counter(frame.mapping_status)
        mapped = sum(counts[s] for s in MAPPED_STATUSES)
        addressable = frame.gene_addressable
        addressable_mapped = int((frame.mapping_status.isin(MAPPED_STATUSES) & addressable).sum())
        addressable_count = int(addressable.sum())
        genes = sorted(set(frame.loc[frame.mapping_status.isin(MAPPED_STATUSES), 'hgnc_symbol']))
        # A source may enter on its raw coverage, or — only when a foreign identifier namespace
        # provably accounts for >=90% of its unmapped rows — on its gene-addressable coverage.
        outside_features = sum(ev['features'] for ev in evidence.values()
                               if ev['outside_hgnc_gene_namespace'])
        raw_passes = mapped / len(frame) >= ENTRY_THRESHOLD
        addressable_passes = bool(addressable_count) and \
            addressable_mapped / addressable_count >= ENTRY_THRESHOLD
        explained = bool(counts[UNMAPPED]) and outside_features / counts[UNMAPPED] >= 0.9
        unit.update(raw_feature_count=int(len(frame)), raw_mapped_feature_count=int(mapped),
                    raw_mapping_fraction=round(mapped / len(frame), 4),
                    gene_addressable_feature_count=addressable_count,
                    gene_addressable_mapped_count=addressable_mapped,
                    gene_addressable_mapping_fraction=round(
                        addressable_mapped / addressable_count, 4) if addressable_count else 0.0,
                    non_addressable_feature_count=int(len(frame)) - addressable_count,
                    non_addressable_classes=evidence,
                    raw_coverage_passes_entry_rule=bool(raw_passes),
                    gene_addressable_coverage_passes_entry_rule=bool(addressable_passes),
                    addressable_admission_evidenced=bool(explained),
                    unique_canonical_genes=len(genes),
                    exact_mapping_fraction=round(counts[EXACT] / len(frame), 4),
                    alias_fraction=round(counts[ALIAS] / len(frame), 4),
                    previous_symbol_fraction=round(counts[PREVIOUS] / len(frame), 4),
                    platform_annotation_fraction=round(counts[PLATFORM_ANNOTATION] / len(frame), 4),
                    ambiguous_fraction=round(counts[AMBIGUOUS] / len(frame), 4),
                    unmapped_fraction=round(counts[UNMAPPED] / len(frame), 4),
                    n_exact=counts[EXACT], n_alias=counts[ALIAS], n_previous_symbol=counts[PREVIOUS],
                    n_platform_annotation=counts[PLATFORM_ANNOTATION], n_ambiguous=counts[AMBIGUOUS],
                    n_unmapped=counts[UNMAPPED],
                    duplicate_features=int(frame.original_feature_id.duplicated().sum()),
                    features_per_canonical_gene=round(mapped / len(genes), 3) if genes else 0.0,
                    enters_core_gene_space=bool(raw_passes or (addressable_passes and explained)),
                    unmapped_methods=dict(Counter(
                        frame.loc[frame.mapping_status.eq(UNMAPPED), 'mapping_method'])))
        gene_sets[key] = set(genes)
        metrics[key] = unit
        print(f"{key}: raw {unit['raw_mapping_fraction']:.1%} | gene-addressable "
              f"{unit['gene_addressable_mapping_fraction']:.1%} of "
              f"{unit['gene_addressable_feature_count']}/{unit['raw_feature_count']} features | "
              f"{unit['unique_canonical_genes']} genes | core={unit['enters_core_gene_space']}")
    metrics_frame = pd.DataFrame(metrics.values()).reindex(columns=METRIC_COLUMNS)
    metrics_frame['unmapped_methods'] = metrics_frame.unmapped_methods.map(json.dumps)
    metrics_frame['non_addressable_classes'] = metrics_frame.non_addressable_classes.map(json.dumps)
    metrics_frame.to_csv(GENE_SPACE / 'SOURCE_MAPPING_METRICS.csv', index=False)

    # Gene-level space: union over every mapped source, core over sources above the rule.
    eligible = {k for k, g in gene_sets.items() if metrics[k]['enters_core_gene_space']}
    genes = {}
    for key in sorted(gene_sets):
        modality = metrics[key]['modalities']
        for gene in gene_sets[key]:
            entry = genes.setdefault(gene, dict(hgnc_symbol=gene, n_sources=0,
                                                n_eligible_sources=0, n_bulk_sources=0,
                                                n_array_sources=0, n_single_cell_sources=0,
                                                sources=[]))
            entry['n_sources'] += 1
            entry['n_eligible_sources'] += int(key in eligible)
            entry['n_bulk_sources'] += int('bulk' in modality)
            entry['n_array_sources'] += int('microarray' in modality or 'array' in modality)
            entry['n_single_cell_sources'] += int('scRNA' in modality or 'snRNA' in modality)
            entry['sources'].append(key)
    gene_rows = []
    for gene in sorted(genes):
        entry = genes[gene]
        gene_id, ensembl, entrez = ref.fields(gene)
        source_list = sorted(entry.pop('sources'))
        entry.update(sources='|'.join(source_list), canonical_gene_id=gene_id,
                     ensembl_gene_id=ensembl, entrez_gene_id=entrez,
                     locus_group=ref.table[gene]['locus_group'],
                     n_cohorts=len({c for key in source_list
                                    for c in metrics[key]['cohort_codes'].split('|')}))
        gene_rows.append(entry)
    genes_frame = pd.DataFrame(gene_rows)
    n_eligible = len(eligible)
    genes_frame['in_core_gene_space'] = genes_frame.n_eligible_sources.eq(n_eligible)
    genes_frame = genes_frame[['hgnc_symbol', 'canonical_gene_id', 'ensembl_gene_id',
                               'entrez_gene_id', 'locus_group', 'n_sources', 'n_eligible_sources',
                               'n_cohorts', 'n_bulk_sources', 'n_array_sources',
                               'n_single_cell_sources', 'in_core_gene_space', 'sources']]
    genes_frame.to_csv(GENE_SPACE / 'GENE_SPACE_GENES.csv.gz', index=False,
                       compression={'method': 'gzip', 'mtime': 0})

    # Source-to-source gene overlap: coverage evidence for choosing a modelling universe.
    keys = sorted(gene_sets)
    overlap = pd.DataFrame([[len(gene_sets[a] & gene_sets[b]) for b in keys] for a in keys],
                           index=keys, columns=keys)
    overlap.to_csv(GENE_SPACE / 'SOURCE_GENE_OVERLAP.csv')

    # Cohort x modality inventory.
    rows = []
    for (cohort, modality), group in samples.groupby(['cohort_code', 'modality'], sort=True):
        paths = sorted({p for p in group.expression_relpath if p})
        local = [p for p in paths if (PROJECT / p).exists()]
        found = sorted(k for k, s in SOURCES.items() if s['file'] in local)
        rows.append(dict(cohort_code=cohort, modality=modality,
                         platform='|'.join(sorted({p or '' for p in group.platform})),
                         n_samples=len(group), n_expression_files=len(paths),
                         n_files_local=len(local),
                         expression_file_status='LOCAL_SOURCE_PRESENT' if local else 'SOURCE_NOT_LOCAL',
                         source_expressions='|'.join(found),
                         mapped_features=sum(metrics[k]['raw_mapped_feature_count'] for k in found),
                         unique_canonical_genes=len(set().union(*(gene_sets[k] for k in found)))
                         if found else 0,
                         enters_core_gene_space=any(metrics[k]['enters_core_gene_space']
                                                    for k in found) if found else False))
    inventory = pd.DataFrame(rows)
    inventory.to_csv(GENE_SPACE / 'FEATURE_IDENTIFIER_INVENTORY.csv', index=False)

    # Pair level: which longitudinal intervals have a canonical gene space on both ends.
    file_to_source = {spec['file']: key for key, spec in SOURCES.items()}

    def pair_status(row):
        refs = {row.expression_t0_ref, row.expression_t1_ref}
        local = all(r and r != 'NOT_AVAILABLE' and (PROJECT / r).exists() for r in refs)
        if not local:
            return 'EXPRESSION_FILE_NOT_LOCAL'
        if len(refs) > 1:
            return 'MIXED_ENDPOINT_SOURCES'
        key = file_to_source.get(next(iter(refs)), 'UNDECLARED_SOURCE')
        return key if key in eligible else f'BELOW_CORE_ENTRY_RULE:{key}'

    pairs = master[master.paired_valid].copy()
    pairs['gene_space_status'] = [pair_status(row) for row in pairs.itertuples()]
    usable = pairs.gene_space_status.isin(eligible)
    classes = {'bulk RNA-seq': 'BULK_RNA', 'microarray': 'MICROARRAY',
               'two-channel microarray': 'MICROARRAY', 'scRNA-seq': 'SCRNA_PSEUDOBULK',
               'snRNA-seq': 'OTHER_MODALITY', 'spatial': 'OTHER_MODALITY'}
    pairs['modality_class'] = pairs.expression_modality.map(classes).fillna('OTHER_MODALITY')
    pair_table = pairs[['pair_uid', 'patient_uid', 'cohort_code', 'expression_modality',
                        'modality_class', 'platform', 'expression_t0_ref', 'expression_t1_ref',
                        'gene_space_status', 'clinical_endpoint_available',
                        'strict_prcr_vs_pd_eligible', 'frozen_auo_eligible']].copy()
    pair_table.insert(9, 'source_expression', pairs.gene_space_status.where(usable, ''))
    pair_table['in_canonical_gene_space'] = usable
    pair_table.to_csv(GENE_SPACE / 'PAIR_GENE_SPACE.csv.gz', index=False,
                      compression={'method': 'gzip', 'mtime': 0})

    with_expression = pairs[~pairs.gene_space_status.eq('EXPRESSION_FILE_NOT_LOCAL')]
    intersection_all = set.intersection(*(g for g in gene_sets.values() if g))
    excluded = {k: int(v) for k, v in pairs.gene_space_status.value_counts().items()
                if k not in eligible and k != 'EXPRESSION_FILE_NOT_LOCAL'}

    def coverage_verdict(key):
        unit = metrics[key]
        classes = unit['non_addressable_classes']
        outside = {name: ev for name, ev in classes.items() if ev['outside_hgnc_gene_namespace']}
        unannotated = classes.get('unannotated_array_feature', {}).get('features', 0)
        parts = [f"raw coverage {unit['raw_mapped_feature_count']}/{unit['raw_feature_count']} "
                 f"features ({unit['raw_mapping_fraction']:.1%}) carry an HGNC gene"]
        if outside:
            blocks = '; '.join(
                f"{ev['features']} features are named by the {name} scheme, of whose "
                f"{ev['distinct_tokens']} distinct tokens only {ev['in_hgnc_namespace']} appear "
                'anywhere in the HGNC namespace'
                for name, ev in sorted(outside.items()))
            parts.append(f"gene-addressable coverage "
                         f"{unit['gene_addressable_mapped_count']}/"
                         f"{unit['gene_addressable_feature_count']} "
                         f"({unit['gene_addressable_mapping_fraction']:.1%}) once "
                         f"non-gene-namespace classes are removed from the denominator: {blocks}")
        if unannotated:
            parts.append(f"{unannotated} features are probe ids whose platform annotation carries "
                         'no gene symbol and no gene id; they stay inside the gene-addressable '
                         'denominator because the missing information is the platform annotation, '
                         'not a foreign identifier namespace')
        return '. '.join(parts) + '.'

    admitted = sorted(k for k in metrics if metrics[k]['enters_core_gene_space']
                      and not metrics[k]['raw_coverage_passes_entry_rule'])
    below = sorted(k for k in metrics if not metrics[k]['enters_core_gene_space'])
    unevidenced = sorted(k for k in metrics
                         if metrics[k]['gene_addressable_coverage_passes_entry_rule']
                         and not metrics[k]['addressable_admission_evidenced']
                         and not metrics[k]['raw_coverage_passes_entry_rule'])

    report = dict(
        gene_space_version='v0.2.1-gene-space',
        built_from='dataset/releases/v0.1.2 (MASTER + samples) + locally present expression files '
                   '+ HGNC reference + GEO platform annotations',
        canonical_coordinate='HGNC approved gene symbol; one row per native feature, and the '
                             'native identifier is never overwritten',
        reference=dict(hgnc_origin=HGNC_ORIGIN, hgnc_sha256=ref.digest,
                       hgnc_approved_entries=ref.n_rows, hgnc_protein_coding=ref.n_protein_coding,
                       hgnc_version_token=ref.version,
                       hgnc_note='the 16.9 MB hgnc_complete_set.txt reference is not committed; it '
                                 'is re-downloaded from hgnc_origin and checked against this digest',
                       platform_annotations={
                           platform: dict(file=relpath(RAW / platform / filename),
                                          sha256=sha256(RAW / platform / filename),
                                          symbol_column=symbol_column,
                                          entrez_column=entrez_column,
                                          ensembl_column=ensembl_column)
                           for platform, (filename, symbol_column, entrez_column, ensembl_column)
                           in PLATFORM_ANNOTATIONS.items()}),
        mapping_status_vocabulary=list(STATUSES),
        core_entry_rule=f'a source enters the core gene space when >= {ENTRY_THRESHOLD:.0%} of its '
                        'raw native features resolve to exactly one canonical gene, or when >= '
                        f'{ENTRY_THRESHOLD:.0%} of its gene-addressable features do and a foreign '
                        'identifier namespace explains >=90% of its unmapped rows; both fractions '
                        'are reported for every source and no native row is ever deleted',
        headline=dict(cohorts_with_expression=int(inventory.loc[
            inventory.expression_file_status.eq('LOCAL_SOURCE_PRESENT'), 'cohort_code'].nunique()),
            cohort_modality_slots_with_expression=int(inventory.expression_file_status
                                                      .eq('LOCAL_SOURCE_PRESENT').sum()),
            expression_sources=len(SOURCES), sources_mapped=int(len(gene_sets)),
            sources_entering_core=len(eligible),
            samples_with_local_expression_file=int(sum(
                1 for p in samples.expression_relpath if p and (PROJECT / p).exists())),
            samples_with_canonical_gene_space=int(sum(
                1 for p in samples.expression_relpath if file_to_source.get(p) in eligible)),
            pairs_with_t0_t1_expression=int(len(with_expression)),
            pairs_with_canonical_gene_space=int(usable.sum()),
            patients_with_t0_t1_expression=int(with_expression.patient_uid.nunique()),
            patients_with_canonical_gene_space=int(pairs.loc[usable, 'patient_uid'].nunique()),
            genes_in_union=len(genes_frame),
            genes_in_intersection_all_sources=len(intersection_all),
            genes_in_core_all_eligible_sources=int(genes_frame.in_core_gene_space.sum()),
            core_gene_space_protein_coding=int((genes_frame.in_core_gene_space & genes_frame
                                                .locus_group.eq('protein-coding gene')).sum())),
        pair_status_breakdown=excluded,
        modality_pair_classes={
            name: dict(pairs=int((with_expression.modality_class.eq(name)).sum()),
                       pairs_in_canonical_gene_space=int(
                           (with_expression.modality_class.eq(name) & usable).sum()),
                       patients=int(with_expression.loc[with_expression.modality_class.eq(name),
                                                        'patient_uid'].nunique()))
            for name in sorted(with_expression.modality_class.unique())},
        core_gene_space_tiers={
            'in_all_eligible_sources': int(genes_frame.in_core_gene_space.sum()),
            'in_at_least_90_percent_of_eligible_sources': int((genes_frame.n_eligible_sources
                                                               >= 0.9 * n_eligible).sum()),
            'in_at_least_75_percent': int((genes_frame.n_eligible_sources >= 0.75 * n_eligible).sum()),
            'in_at_least_half': int((genes_frame.n_eligible_sources >= 0.5 * n_eligible).sum()),
            'in_all_sources_including_low_coverage': int(genes_frame.n_sources.eq(len(gene_sets)).sum())},
        coverage_denominators=dict(
            raw_mapping_fraction='mapped native features / all native features of the published '
                                 'matrix, as it stands',
            gene_addressable_mapping_fraction='mapped features / features whose native identifier '
                                              'belongs to a class represented in the HGNC gene '
                                              'namespace; a class leaves the denominator only on '
                                              'measured namespace evidence, never because it '
                                              'failed to map',
            entry_rule=f'entry needs >= {ENTRY_THRESHOLD:.0%} of raw features mapped, or >= '
                       f'{ENTRY_THRESHOLD:.0%} of gene-addressable features mapped when a foreign '
                       'identifier namespace explains >=90% of the unmapped rows',
            anti_circularity='both denominators are reported for every source, and every native '
                             'row stays in its feature map with its identifier class recorded',
            class_definitions={name: dict(identifier_shape=description,
                                          always_gene_addressable=forced)
                               for name, _, forced, description in IDENTIFIER_CLASSES}),
        adjudication=dict(
            admitted_by_gene_addressable_rule=admitted,
            below_entry_rule=below,
            feature_universe_unresolved=unevidenced,
            verdict='RESOLVED' if not unevidenced else 'UNRESOLVED',
            note='a foreign identifier namespace excuses a low raw fraction only when it '
                 'accounts for >=90% of the unmapped rows; missing platform annotation does not, '
                 'so an array whose probes were never assigned genes stays out'),
        low_coverage=dict(
            policy='no source is hidden: every source keeps its full feature map and contributes '
                   'to the union, and both denominators are reported side by side',
            verdicts={key: coverage_verdict(key) for key in sorted(set(below) | set(admitted))}),
        sources=metrics,
        explicit_non_goals=[
            'expression values are NOT placed on one numerical scale here: the canonical gene '
            'space answers which biological feature a column represents, not whether the '
            'measured values are comparable',
            'no global normalization or ComBat across cohorts',
            'no within-cohort normalization (Level C), no rank representation (Level D), no '
            'paired delta representation (Level E)',
            'ambiguous features are recorded with all candidate genes, never assigned to one by '
            'majority vote, expression correlation or genomic proximity',
            'sample-level binding of a MASTER sample id to a matrix column is not asserted here; '
            'a pair enters the canonical gene space when both endpoint references resolve to the '
            'same core-rule expression source',
            'GSE319641 features are already author-ComBat-corrected TPMs, so that source carries '
            'another laboratory batch correction before any PAIR-side representation',
            'no expression matrices are committed; only feature-to-gene mappings'],
    )
    (GENE_SPACE / 'GENE_SPACE_REPORT.json').write_text(
        json.dumps(report, indent=2, default=str) + '\n', encoding='utf-8')
    print(json.dumps(dict(headline=report['headline'], adjudication=report['adjudication'],
                          pair_status_breakdown=excluded,
                          modality_pair_classes=report['modality_pair_classes'],
                          core_gene_space_tiers=report['core_gene_space_tiers']), indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
