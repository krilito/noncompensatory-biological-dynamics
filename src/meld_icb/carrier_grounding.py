"""Carrier composition summaries; no mechanistic interpretation is implied."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import gzip
import hashlib


def cell_count_shares(counts: Mapping[str, int]) -> dict[str, float]:
    """Return cell-count shares without implying signal contribution."""
    total = sum(int(v) for v in counts.values())
    if total <= 0:
        raise ValueError("carrier counts must have positive total")
    return {key: int(value) / total for key, value in counts.items()}


def weighted_signal_shares(signal: Mapping[str, float]) -> dict[str, float]:
    """Return non-negative carrier signal shares (not cell-count shares)."""
    values = {key: float(value) for key, value in signal.items()}
    if any(value < 0 for value in values.values()):
        raise ValueError("weighted signal shares require non-negative values")
    total = sum(values.values())
    if total <= 0:
        raise ValueError("weighted signal must have positive total")
    return {key: value / total for key, value in values.items()}


_MELD_AXES = {
    "T_axis": ["MKI67", "PCNA", "TOP2A", "CCNB1", "CDK1", "MCM2", "MCM5", "E2F1", "AURKB"],
    "E_axis": ["CD8A", "CD8B", "GZMB", "GZMA", "PRF1", "IFNG", "NKG7", "GNLY", "CXCL9", "CXCL10"],
    "X_axis": ["TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "COL1A1", "COL1A2", "FN1", "VIM", "ZEB1", "SNAI1", "ACTA2", "CXCL12", "CD163", "MRC1", "CSF1R", "IL10", "ARG1", "IDO1", "VEGFA", "KDR", "FLT1"],
    "C_axis": ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "TOX2"],
}
_PUBLISHED_SIGNATURES = {
    "CYT": ["GZMA", "PRF1"],
    "IFNG6": ["IFNG", "CXCL9", "CXCL10", "IDO1", "STAT1", "HLA-DRA"],
    "TIS18": ["CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7", "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB"],
    "CD274": ["CD274"],
    "MKI67": ["MKI67"],
}


def established_signature_gene_sets() -> dict[str, tuple[str, ...]]:
    """Return the project-locked comparator definitions as immutable values."""
    return {name: tuple(genes) for name, genes in _PUBLISHED_SIGNATURES.items()}


_PANEL = {**_MELD_AXES, **_PUBLISHED_SIGNATURES}
_INTENDED_REFERENT = {
    "T_axis": "tumour-intrinsic", "MKI67": "tumour-intrinsic", "E_axis": "immune", "C_axis": "immune", "CYT": "immune", "IFNG6": "immune", "TIS18": "immune", "X_axis": "mixed-stroma-weighted", "CD274": "mixed",
}
_MARKERS = {
    "CD8_T": ["CD8A", "CD8B"], "CD4_T": ["CD4", "IL7R"], "Treg": ["FOXP3", "IL2RA"], "NK": ["NKG7", "GNLY", "KLRD1", "NCAM1"], "B_cell": ["MS4A1", "CD79A", "CD79B", "CD19"], "Plasma": ["MZB1", "JCHAIN", "XBP1", "DERL3"], "Myeloid": ["LYZ", "CD14", "CD68", "CSF1R", "ITGAX", "AIF1"], "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"], "Malignant": ["MLANA", "PMEL", "TYR", "DCT", "SOX10", "MITF"], "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB"],
}
_TCELL = ["CD3D", "CD3E", "CD3G"]
_LYMPHOID = ("CD8_T", "CD4_T", "Treg", "NK", "B_cell", "Plasma")


class CarrierInputError(ValueError):
    """Carrier raw matrix or lineage mapping cannot satisfy the contract."""


def _carrier_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _carrier_matrix(path: Path):
    import numpy as np
    import pandas as pd

    if not path.exists():
        raise CarrierInputError(f"required C input is missing: {path}")
    needed = sorted({gene for genes in _PANEL.values() for gene in genes} | {gene for genes in _MARKERS.values() for gene in genes} | set(_TCELL))
    wanted = set(needed)
    rows: dict[str, Any] = {}
    with gzip.open(path, "rt", errors="replace") as handle:
        cells = handle.readline().rstrip("\n").split("\t")[1:]
        samples = handle.readline().rstrip("\n").split("\t")[1:]
        if len(cells) != len(samples):
            raise CarrierInputError("C cell/sample header lengths differ")
        for line in handle:
            tab = line.find("\t")
            if tab < 0:
                continue
            gene = line[:tab].strip().upper()
            if gene not in wanted or gene in rows:
                continue
            values = np.fromstring(line[tab + 1:], sep="\t", dtype=np.float32)
            if values.size != len(cells):
                raise CarrierInputError(f"C expression row has {values.size} values; expected {len(cells)}")
            rows[gene] = values
    if not rows:
        raise CarrierInputError("C matrix yielded no declared genes")
    expression = pd.DataFrame(rows, index=cells).T
    expression = np.log2(expression + 1.0)
    return expression, pd.Series(samples, index=cells, name="sample")


def _carrier_score(expr, genes: list[str]):
    import pandas as pd

    available = [gene for gene in genes if gene in expr.index]
    if not available:
        return pd.Series(0.0, index=expr.columns)
    return expr.loc[available].mean(axis=0)


def _assign_carrier_lineage(expr):
    import numpy as np
    import pandas as pd

    scores = {key: _carrier_score(expr, genes) for key, genes in _MARKERS.items()}
    tcell = _carrier_score(expr, _TCELL)
    non_t = ["B_cell", "Plasma", "Myeloid", "Mast", "Malignant", "Fibroblast", "NK"]
    non_t_matrix = pd.DataFrame({key: scores[key] for key in non_t})
    out = pd.Series("Unassigned", index=expr.columns, dtype=object)
    is_t = (tcell > 0.5) & (tcell >= non_t_matrix.max(axis=1))
    out[is_t & (scores["Treg"] > 0.5) & (scores["Treg"] >= scores["CD8_T"])] = "Treg"
    rest_t = is_t & (out == "Unassigned")
    # The source authority's ``>=`` tie-break assigns exact zero-signal ties to
    # CD4_T.  A non-zero tie has no declared authority and is therefore held.
    t_ties = rest_t & np.isclose(scores["CD8_T"], scores["CD4_T"]) & (scores["CD8_T"] > 0.0)
    if bool(t_ties.any()):
        raise CarrierInputError("C lineage mapping ambiguity: tied CD8/CD4 marker scores")
    out[rest_t & (scores["CD8_T"] >= scores["CD4_T"]) & (scores["CD8_T"] > 0.2)] = "CD8_T"
    out[rest_t & (out == "Unassigned")] = "CD4_T"
    rest = out == "Unassigned"
    best = non_t_matrix.loc[rest].idxmax(axis=1)
    bestv = non_t_matrix.loc[rest].max(axis=1)
    ties = non_t_matrix.loc[rest].eq(bestv, axis=0).sum(axis=1) > 1
    if bool((ties & (bestv > 0.5)).any()):
        raise CarrierInputError("C lineage mapping ambiguity: tied non-T marker scores")
    out.loc[rest] = np.where(bestv > 0.5, best, "Unassigned")
    return out


def _carrier_shares(expr, lineage, genes: list[str]):
    import pandas as pd

    signature_score = _carrier_score(expr, genes).clip(lower=0.0)
    frame = pd.DataFrame({"lineage": lineage.values, "score": signature_score.values})
    grouped = frame.groupby("lineage")["score"].agg(["sum", "mean", "count"])
    total = float(grouped["sum"].sum())
    grouped["share"] = grouped["sum"] / total if total > 0 else float("nan")
    return grouped.rename(columns={"mean": "signature_mean", "count": "n_cells"})


def run_carrier_grounding(*, truth_root: str | Path, config_path: str | Path) -> dict[str, Any]:
    """Run C E1 from raw GSE120575 TPM; return aggregate-only evidence."""
    from .config import load_json_yaml

    root = Path(truth_root).resolve()
    config = load_json_yaml(config_path)
    relative = str(config.get("inputs", {}).get("tpm_matrix", ""))
    path = root / relative
    try:
        expr, samples = _carrier_matrix(path)
        on_treatment = samples.str.startswith("Post")
        expr_on = expr.loc[:, on_treatment.values]
        lineage = _assign_carrier_lineage(expr_on)
        t_table = _carrier_shares(expr_on, lineage, _PANEL["T_axis"])
        lymphoid_share = float(t_table.reindex(_LYMPHOID)["share"].fillna(0).sum())
        malignant_share = float(t_table["share"].get("Malignant", 0.0))
        top_two = list(t_table.sort_values("share", ascending=False).index[:2])
        gates = {
            "A_lymphoid_share": {"observed": lymphoid_share, "published": 0.871, "tolerance_pp": 5.0, "pass": abs(lymphoid_share - 0.871) <= 0.05},
            "B_malignant_share": {"observed": malignant_share, "published": 0.0016, "pass": malignant_share < 0.01},
            "C_top_two_carriers": {"observed": top_two, "published": ["CD8_T", "CD4_T"], "pass": set(top_two) == {"CD8_T", "CD4_T"}},
            "n_on_treatment_cells": {"observed": int(on_treatment.sum()), "published_reference": 10363},
        }
        if not all(gate.get("pass", True) for gate in gates.values()):
            return {"analysis_id": "E1_CARRIER_GENERALIZATION", "verdict": "BLOCKED", "status_class": "DESCRIPTIVE", "reason": "lineage/gate reproduction failed; signature result not computed", "gates": gates, "provenance": {"tpm_matrix": relative, "tpm_sha256": _carrier_sha256(path)}}
        summary: list[dict[str, Any]] = []
        for signature, genes in _PANEL.items():
            table = _carrier_shares(expr_on, lineage, genes)
            # The authority's exported aggregate is read back as float64;
            # retain that order-of-operations for the signature summary while
            # leaving the gate's float32 observed value unchanged.
            lymphoid = float(table[table.index.isin(_LYMPHOID)]["share"].astype("float64").sum())
            myeloid = float(table[table.index == "Myeloid"]["share"].astype("float64").sum())
            malignant = float(table[table.index == "Malignant"]["share"].astype("float64").sum())
            dominant = str(table.sort_values("share", ascending=False).iloc[0].name)
            intended = _INTENDED_REFERENT[signature]
            immune_share = lymphoid + myeloid
            summary.append({"signature": signature, "intended_referent": intended, "dominant_carrier": dominant, "lymphoid_share": lymphoid, "myeloid_share": myeloid, "immune_share": immune_share, "malignant_share": malignant, "intent_carrier_concordant": bool(intended == "immune"), "discordant_tumour_intrinsic": bool(intended == "tumour-intrinsic" and immune_share > 0.80)})
        return {
            "analysis_id": "E1_CARRIER_GENERALIZATION",
            "verdict": "GATES_PASSED",
            "status_class": "DESCRIPTIVE",
            "may_upgrade_b5": False,
            "substrate": "GSE120575 CD45+ single cells, on-treatment",
            "analysis_unit": "pooled cell",
            "n_cells": int(on_treatment.sum()),
            "n_signatures": len(_PANEL),
            "gates": gates,
            "tumour_intrinsic_signatures_immune_carried_over_80pct": sorted(item["signature"] for item in summary if item["discordant_tumour_intrinsic"]),
            "summary": summary,
            "semantic_limits": {"patient_ci": "NOT_ESTIMABLE", "patient_inference": "NOT_PERMITTED", "b5_upgrade": False, "routine_single_column_universe": "not mixed with 10363-cell all-cell universe"},
            "provenance": {"tpm_matrix": relative, "tpm_sha256": _carrier_sha256(path), "mapping": "marker-based lineage authority with ambiguity HOLD"},
        }
    except (CarrierInputError, OSError, ValueError, KeyError) as exc:
        return {"analysis_id": "E1_CARRIER_GENERALIZATION", "verdict": "HOLD", "status_class": "DESCRIPTIVE", "hold_reason": str(exc), "provenance": {"tpm_matrix": relative}}
