"""Normalization helpers for the PAIR longitudinal dataset.

All mappings are deterministic and lossless: native values are always preserved
alongside the normalized fields produced here.
"""
from __future__ import annotations

PHASE_MAP = {
    'BASELINE': 'PRE_TREATMENT',
    'PRE_TREATMENT': 'PRE_TREATMENT',
    'EARLY_ON_TREATMENT': 'ON_TREATMENT_EARLY',
    'MID_TREATMENT': 'ON_TREATMENT_LATE',
    'ON_TREATMENT': 'ON_TREATMENT_EARLY',
    'POST_TREATMENT': 'POST_TREATMENT',
    'FOLLOW_UP': 'POST_TREATMENT',
    'SURGERY': 'SURGERY',
    'PROGRESSION': 'PROGRESSION',
    'RECURRENCE': 'RECURRENCE',
    'RELAPSE': 'RECURRENCE',
    'UNKNOWN': 'UNKNOWN',
    '': 'UNKNOWN',
    'NOT_MEASURED': 'UNKNOWN',
}

# Lower-case substring -> canonical drug name, matched against native regimens.
DRUG_LEXICON = [
    ('nivolumab', 'nivolumab'),
    ('ipilimumab', 'ipilimumab'),
    ('pembrolizumab', 'pembrolizumab'),
    ('atezolizumab', 'atezolizumab'),
    ('sotigalimab', 'sotigalimab'),
    ('lenvatinib', 'lenvatinib'),
    ('sitravatinib', 'sitravatinib'),
    ('ibrutinib', 'ibrutinib'),
    ('letrozole', 'letrozole'),
    ('paclitaxel', 'paclitaxel'),
    ('docetaxel', 'docetaxel'),
    ('capecitabine', 'capecitabine'),
    ('carboplatin', 'carboplatin'),
    ('temozolomide', 'temozolomide'),
    ('decitabine', 'decitabine'),
    ('disulfiram', 'disulfiram'),
    ('trastuzumab', 'trastuzumab'),
    ('t-dm1', 'trastuzumab_emtansine'),
    ('pertuzumab', 'pertuzumab'),
    ('neratinib', 'neratinib'),
    ('bevacizumab', 'bevacizumab'),
    ('ganitumab', 'ganitumab'),
    ('ganetespib', 'ganetespib'),
    ('mk-2206', 'MK-2206'),
    ('amg 386', 'AMG-386'),
    ('amg-386', 'AMG-386'),
    ('abt 888', 'veliparib'),
    ('bo-112', 'BO-112'),
    ('ipi pd1', 'ipilimumab+nivolumab'),
    ('ipipd1', 'ipilimumab+nivolumab'),
    ('pd1', 'anti-PD1'),
    ('ctla4', 'anti-CTLA4'),
]

CLASS_LABEL = {
    'ANTI_PD1': 'anti_PD1',
    'ANTI_PDL1': 'anti_PDL1',
    'ANTI_CTLA4': 'anti_CTLA4',
    'CHEMOTHERAPY': 'chemotherapy',
    'ENDOCRINE': 'endocrine',
    'RADIOTHERAPY': 'radiotherapy',
    'CHEMORADIATION': 'chemoradiation',
    'CELL_THERAPY': 'cellular_therapy',
    'TARGETED_KINASE': 'targeted',
    'ANTI_HER2': 'anti_HER2',
    'ANTI_ANGIOGENIC': 'anti_angiogenic',
    'ICI_COMBINATION': 'ICI_combination',
    'ICI_CHEMO_COMBINATION': 'ICI_chemo_combination',
    'OTHER': 'other',
    'NOT_MEASURED': 'NOT_MEASURED',
}

ICB_CLASSES = {'ANTI_PD1', 'ANTI_PDL1', 'ANTI_CTLA4'}

ENDPOINT_SYSTEM = {
    'RECIST_RESPONSE': 'RECIST',
    'LESION_LEVEL_RECIST_RESPONSE': 'RECIST_1_1',
    'RECIST_BINARY_RESPONSE': 'author_defined_composite',
    'PCR': 'pCR',
    'CHEMOTHERAPY_RESPONSE_SCORE': 'CRS',
    'CLINICAL_RESPONSE_ULTRASOUND_VOLUME': 'ultrasound_volume_response',
    'CLINICAL_RESPONSE_COMPOSITE_NEOADJUVANT': 'author_defined_composite',
    'DFS': 'recurrence',
    'RELAPSE': 'recurrence',
    'ENDPOINT_SEMANTICS_UNCERTAIN': 'uncertain',
}

ENDPOINT_FAMILY_MASTER = {
    'RECIST_RESPONSE': 'radiologic_response',
    'LESION_LEVEL_RECIST_RESPONSE': 'radiologic_response',
    'RECIST_BINARY_RESPONSE': 'radiologic_response',
    'PCR': 'pathological_response',
    'CHEMOTHERAPY_RESPONSE_SCORE': 'pathological_response',
    'CLINICAL_RESPONSE_ULTRASOUND_VOLUME': 'clinical_composite_response',
    'CLINICAL_RESPONSE_COMPOSITE_NEOADJUVANT': 'clinical_composite_response',
    'DFS': 'recurrence',
    'RELAPSE': 'recurrence',
    'ENDPOINT_SEMANTICS_UNCERTAIN': 'other_clinical_endpoint',
}

RECIST_HARMONIZED = {'CR': 'RESPONSE', 'PR': 'RESPONSE', 'PRCR': 'RESPONSE',
                     'CR_OR_PR': 'RESPONSE', 'SD': 'STABLE', 'PD': 'PROGRESSION'}

SOURCE_PUBLICATION = {
    'GSE165897': 'https://doi.org/10.1126/sciadv.abm1831',
    'GSE179994': 'https://doi.org/10.1038/s43018-021-00292-8',
    'GSE20181': 'https://doi.org/10.1186/bcr2611',
    'GSE5462': 'https://doi.org/10.1186/bcr2611',
    'GSE18728': 'https://pubmed.ncbi.nlm.nih.gov/20012355/',
    'GSE120575': 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6641984/',
}

GENE_SPACE = {
    'bulk RNA-seq': 'GENE_LEVEL',
    'microarray': 'PROBE_BASED',
    'two-channel microarray': 'PROBE_BASED',
    'scRNA-seq': 'GENE_LEVEL_CELLULAR',
    'snRNA-seq': 'GENE_LEVEL_CELLULAR',
    'spatial': 'GENE_LEVEL_SPATIAL',
    '': 'UNKNOWN',
    'NOT_MEASURED': 'UNKNOWN',
}


def normalize_phase(timepoint_normalized: str) -> str:
    return PHASE_MAP.get(str(timepoint_normalized), 'UNKNOWN')


def parse_drugs(native: str) -> str:
    text = str(native).lower()
    if not text or text in ('not_measured',):
        return 'NOT_PARSED'
    found = []
    for token, drug in DRUG_LEXICON:
        if token in text and drug not in found:
            found.append(drug)
    return ' + '.join(found) if found else 'NOT_PARSED'


def classify_treatment(native: str, classes: str):
    """Return (drug, drug_class, treatment_family) for one regimen.

    `classes` is the corpus treatment-class pipe list (e.g. ANTI_PD1|ANTI_CTLA4|ICI_COMBINATION).
    """
    tokens = [t for t in str(classes).split('|') if t and t != 'NOT_MEASURED']
    token_set = set(tokens)
    drug = parse_drugs(native)
    drug_class = ' + '.join(CLASS_LABEL.get(t, t.lower()) for t in tokens) if tokens else 'NOT_MEASURED'
    if not token_set or token_set == {'OTHER'}:
        return drug, drug_class or 'NOT_MEASURED', 'NOT_MEASURED' if not token_set else 'OTHER'
    has_icb = bool(token_set & ICB_CLASSES)
    dual_icb = {'ANTI_PD1', 'ANTI_CTLA4'} <= token_set or {'ANTI_PDL1', 'ANTI_CTLA4'} <= token_set
    has_chemo = 'CHEMOTHERAPY' in token_set
    has_rt = 'RADIOTHERAPY' in token_set
    if 'CHEMORADIATION' in token_set or (has_rt and has_chemo):
        family = 'CHEMORADIATION'
    elif dual_icb and not has_chemo:
        family = 'DUAL_ICB'
    elif has_icb and has_chemo:
        family = 'ICB_CHEMOTHERAPY'
    elif has_icb and (token_set - ICB_CLASSES - {'ICI_COMBINATION'}):
        family = 'ICB_COMBINATION'
    elif has_icb:
        family = 'ICB'
    elif 'CELL_THERAPY' in token_set:
        family = 'CELLULAR_THERAPY'
    elif 'ENDOCRINE' in token_set:
        family = 'ENDOCRINE_THERAPY'
    elif has_chemo:
        family = 'CHEMOTHERAPY'
    elif has_rt:
        family = 'RADIOTHERAPY'
    elif token_set & {'TARGETED_KINASE', 'ANTI_HER2', 'ANTI_ANGIOGENIC'}:
        family = 'TARGETED_THERAPY'
    else:
        family = 'COMBINATION_THERAPY' if len(token_set) > 1 else 'OTHER'
    return drug, drug_class, family


def endpoint_systems(families) -> str:
    systems = []
    for family in families:
        system = ENDPOINT_SYSTEM.get(family, 'other')
        if system not in systems:
            systems.append(system)
    return '|'.join(systems) if systems else ''


def endpoint_families_master(families) -> str:
    master = []
    for family in families:
        mapped = ENDPOINT_FAMILY_MASTER.get(family, 'other_clinical_endpoint')
        if mapped not in master:
            master.append(mapped)
    return '|'.join(master) if master else ''


def harmonize_endpoint(values, families) -> str:
    """CR/PR->RESPONSE, SD->STABLE, PD->PROGRESSION for standard RECIST only.

    Every other endpoint family keeps its native semantics and is marked NOT_HARMONIZED;
    this field never replaces the native label.
    """
    out = []
    for value, family in zip(values, families):
        if family == 'RECIST_RESPONSE' and str(value).upper() in RECIST_HARMONIZED:
            tag = RECIST_HARMONIZED[str(value).upper()]
        else:
            tag = 'NOT_HARMONIZED'
        if tag not in out:
            out.append(tag)
    return '|'.join(out) if out else ''


def gene_space_for(modality: str) -> str:
    return GENE_SPACE.get(str(modality), 'UNKNOWN')
