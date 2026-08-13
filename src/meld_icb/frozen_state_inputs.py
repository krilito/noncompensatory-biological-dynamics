"""Read-only adapters for the declared external frozen-state inputs.

The A producer accepts a small, generic row table.  This module is the only
adapter from the source-study matrices to that table.  It deliberately keeps
source identifiers and mapping tables outside the repository: callers provide
an explicit truth-repository root and an external audit output directory.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .frozen_state_transfer import AXIS_ORDER, run_frozen_state_transfer
from .statistics import auc, registered_auc_inference


class FrozenInputError(ValueError):
    """A source matrix or mapping cannot satisfy its declared contract."""


_GENE_FIELDS = (
    "gene.hgnc.symbol", "gene", "gene_symbol", "hgnc_symbol", "symbol",
    "geneid", "gene_id", "gene_symbol_id",
)
_MORRISON_GENE_FIELD = "gene.hgnc.symbol"
_PRJEB_REL = "raw_data/PRJEB23709/cancercell_normalized_counts_genenames.txt"
_PRJEB_LABELS = "metadata/v2_external/gide_passon_on_treatment_labels.csv"
_MORRISON_REL = "raw_data/MORRISON-1-public/RNASeq/data/RNA-CancerCell-MORRISON1-no_batch_correction-logcpm-all_samples.tsv"
_MORRISON_LABELS = "metadata/v2_external/morrison_candidate_on_treatment_subset.tsv"
_MGH_115_REL = "raw_data/GSE115821_MGH_counts.csv.gz"
_MGH_168_REL = "raw_data/GSE168204_MGH_counts.csv.gz"
_MGH_MAP = "metadata/v2_external/mgh_passon_sample_mapping.tsv"
_MGH_LABELS = "metadata/v2_external/mgh_passon_on_treatment_labels.csv"
_GSE91061_EXPR = "processed/expression_log2_GSE91061_symbol.csv"
_GSE91061_SOFT = "raw_data/GSE91061_family.soft.gz"
_HISTORICAL_COHORT_NAMES = {
    "PRJEB23709": "Gide / PRJEB23709",
    "MORRISON-1-public": "MORRISON-1-public",
    "MGH": "MGH / GSE115821 + GSE168204",
}


def _open_text(path: Path):
    if path.suffix == ".gz" or path.name.endswith(".csv.gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def _rows(path: Path, delimiter: str) -> list[dict[str, str]]:
    with _open_text(path) as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _require(path: Path) -> Path:
    if not path.exists():
        raise FrozenInputError(f"required source input is missing: {path}")
    return path


def _norm(value: object) -> str:
    return str(value).strip().upper()


def _float(value: object) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _binary(value: object) -> int:
    normalized = _norm(value)
    if normalized in {"1", "R", "RESPONDER", "CR", "PR", "CRPR"}:
        return 1
    if normalized in {"0", "NR", "NONRESPONDER", "PD"}:
        return 0
    raise FrozenInputError(f"unsupported response label: {value!r}")


def _collapse_matrix(
    values: Mapping[str, Mapping[str, list[float]]],
    sample_ids: Sequence[str],
    *,
    transform: str,
) -> list[dict[str, float]]:
    """Collapse duplicate symbols, then apply the deposited-scale transform."""
    result: list[dict[str, float]] = []
    for sample_id in sample_ids:
        row: dict[str, float] = {}
        for gene, by_sample in values.items():
            observations = by_sample.get(sample_id, [])
            if not observations:
                continue
            value = sum(observations) / len(observations)
            if transform == "truncate_negative_then_log2p1":
                value = math.log2(max(value, 0.0) + 1.0)
            elif transform != "as_deposited_logcpm":
                raise FrozenInputError(f"unknown expression transform: {transform}")
            if not math.isfinite(value):
                raise FrozenInputError(f"non-finite transformed expression: {gene}/{sample_id}")
            row[gene] = value
        if not row:
            raise FrozenInputError(f"sample has no mapped expression values: {sample_id}")
        result.append(row)
    return result


def _read_selected_matrix(
    path: Path,
    sample_ids: Sequence[str],
    *,
    delimiter: str,
    gene_field: str | None = None,
    transform: str,
) -> list[dict[str, float]]:
    """Stream a matrix while retaining only declared samples and genes."""
    wanted = set(sample_ids)
    with _open_text(_require(path)) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = list(reader.fieldnames or [])
        if not fields:
            raise FrozenInputError(f"matrix has no header: {path}")
        selected = [sid for sid in sample_ids if sid in fields]
        if len(selected) != len(sample_ids):
            missing = [sid for sid in sample_ids if sid not in fields]
            raise FrozenInputError(f"matrix columns missing ({path.name}): {missing[:5]}")
        if gene_field is None:
            lower: dict[str, str] = {}
            for field in fields:
                lower.setdefault(field.lower(), field)
            gene_field = next((lower[key] for key in _GENE_FIELDS if key in lower), fields[0])
        elif gene_field not in fields:
            raise FrozenInputError(f"gene column missing ({path.name}): {gene_field}")
        values: dict[str, dict[str, list[float]]] = {}
        for row in reader:
            gene = _norm(row.get(gene_field, ""))
            if not gene:
                continue
            # The adapter only needs genes in the frozen vocabulary.  This is
            # also what makes streaming the large MGH/PRJEB matrices cheap.
            by_sample = values.setdefault(gene, {sid: [] for sid in selected})
            for sid in selected:
                value = _float(row.get(sid))
                if value is not None:
                    by_sample[sid].append(value)
        return _collapse_matrix(values, selected, transform=transform)


def _gene_set_union(freeze: Mapping[str, Any]) -> set[str]:
    return {
        _norm(gene)
        for genes in freeze["method"]["axes"]["gene_sets"].values()
        for gene in genes
    }


def _read_selected_matrix_fast(
    path: Path,
    sample_ids: Sequence[str],
    *,
    delimiter: str,
    gene_field: str | None,
    transform: str,
    genes: set[str],
) -> list[dict[str, float]]:
    """Same as :func:`_read_selected_matrix`, with gene filtering in the loop."""
    wanted = set(sample_ids)
    with _open_text(_require(path)) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = list(reader.fieldnames or [])
        selected = [sid for sid in sample_ids if sid in fields]
        if len(selected) != len(sample_ids):
            raise FrozenInputError(f"matrix columns missing ({path.name})")
        if gene_field is None:
            lower: dict[str, str] = {}
            for field in fields:
                lower.setdefault(field.lower(), field)
            gene_field = next((lower[key] for key in _GENE_FIELDS if key in lower), fields[0])
        values: dict[str, dict[str, list[float]]] = {}
        for row in reader:
            gene = _norm(row.get(gene_field, ""))
            if not gene or gene not in genes:
                continue
            by_sample = values.setdefault(gene, {sid: [] for sid in selected})
            for sid in selected:
                value = _float(row.get(sid))
                if value is not None:
                    by_sample[sid].append(value)
    return _collapse_matrix(values, selected, transform=transform)


def _with_metadata(
    expression: Sequence[Mapping[str, float]],
    labels: Sequence[Mapping[str, str]],
    *,
    sample_field: str,
    patient_field: str,
    response_field: str,
) -> list[dict[str, Any]]:
    by_sample = {str(row[sample_field]): row for row in labels if str(row.get(sample_field, "")).strip()}
    output: list[dict[str, Any]] = []
    # The expression adapter preserves label-table order, which is the stable
    # source order used for deterministic patient tie-breaking.
    for label in labels:
        sid = str(label.get(sample_field, "")).strip()
        if not sid or sid not in by_sample:
            continue
        if len(output) >= len(expression):
            break
        # expression rows are already in the same selected sample order.
        row = dict(expression[len(output)])
        row.update({
            "sample_id": sid,
            "patient_id": str(label.get(patient_field, sid)).strip(),
            "response": _binary(label.get(response_field)),
        })
        for key in ("timepoint", "therapy", "biopsy_order", "passon_day", "expression_source"):
            if key in label and str(label[key]).strip():
                row[key] = label[key]
        output.append(row)
    if len(output) != len(expression):
        raise FrozenInputError(f"expression/label row mismatch: expression={len(expression)} labels={len(output)}")
    return output


def _prjeb(root: Path, freeze: Mapping[str, Any]) -> tuple[list[dict[str, Any]], dict[str, str]]:
    labels_path = _require(root / _PRJEB_LABELS)
    labels = _rows(labels_path, ",")
    sample_ids = [str(row["sample_id"]) for row in labels]
    expr = _read_selected_matrix_fast(
        root / _PRJEB_REL,
        sample_ids,
        delimiter="\t",
        gene_field="Gene",
        transform="truncate_negative_then_log2p1",
        genes=_gene_set_union(freeze),
    )
    return _with_metadata(expr, labels, sample_field="sample_id", patient_field="patient_id", response_field="response_binary"), {
        "expression_scale": "cancercell_normalized_counts_then_log2p1",
        "expression_path": _PRJEB_REL,
        "labels_path": _PRJEB_LABELS,
    }


def _morrison(root: Path, freeze: Mapping[str, Any]) -> tuple[list[dict[str, Any]], dict[str, str]]:
    labels_path = _require(root / _MORRISON_LABELS)
    labels = _rows(labels_path, "\t")
    sample_ids = [str(row["sample_id"]) for row in labels]
    expr = _read_selected_matrix_fast(
        root / _MORRISON_REL,
        sample_ids,
        delimiter="\t",
        gene_field=_MORRISON_GENE_FIELD,
        transform="as_deposited_logcpm",
        genes=_gene_set_union(freeze),
    )
    return _with_metadata(expr, labels, sample_field="sample_id", patient_field="subject_id", response_field="response_binary_CRPR1_PD0"), {
        "expression_scale": "deposited_logcpm",
        "expression_path": _MORRISON_REL,
        "labels_path": _MORRISON_LABELS,
    }


def _mgh(root: Path, freeze: Mapping[str, Any]) -> tuple[list[dict[str, Any]], dict[str, str]]:
    map_path = _require(root / _MGH_MAP)
    labels_path = _require(root / _MGH_LABELS)
    mapping = _rows(map_path, "\t")
    labels_by_id = {str(row["sample_id"]): row for row in _rows(labels_path, ",")}
    allowed = {"exact", "normalized", "drop_mgh_prefix", "manual_verified"}
    mapped = [row for row in mapping if row.get("mapping_status") in allowed and str(row.get("passon_sample_id", "")) in labels_by_id]
    if len(mapped) != 29:
        raise FrozenInputError(f"MGH mapped sample count is {len(mapped)}, expected 29")
    genes = _gene_set_union(freeze)
    by_source: dict[str, list[dict[str, str]]] = {}
    for row in mapped:
        by_source.setdefault(str(row["expression_source"]), []).append(row)
    source_expression: dict[str, list[tuple[str, dict[str, float]]]] = {}
    for source, rows in by_source.items():
        rel = _MGH_115_REL if source == "GSE115821" else _MGH_168_REL if source == "GSE168204" else None
        if rel is None:
            raise FrozenInputError(f"unsupported MGH expression source: {source}")
        wanted = [str(row["expression_column"]) for row in rows]
        # Matrix column names are mapped by the source mapping table.  Read
        # under expression-column IDs, then rename to stable PASS-ON IDs.
        expr_rows = _read_selected_matrix_fast(
            root / rel,
            wanted,
            delimiter=",",
            gene_field=None,
            transform="truncate_negative_then_log2p1",
            genes=genes,
        )
        source_expression[source] = [
            (str(row_map["passon_sample_id"]), expr)
            for row_map, expr in zip(rows, expr_rows)
        ]
    # The canonical MGH reconstruction aligns the two count matrices on their
    # shared gene vocabulary before concatenating columns.  Applying that
    # intersection here is material for genes present in only one accession.
    source_gene_sets = [
        {gene for _, expr in values for gene in expr}
        for values in source_expression.values()
    ]
    common_genes = set.intersection(*source_gene_sets) if source_gene_sets else set()
    all_expression: dict[str, dict[str, float]] = {
        sample_id: {gene: value for gene, value in expr.items() if gene in common_genes}
        for values in source_expression.values()
        for sample_id, expr in values
    }
    ordered: list[dict[str, Any]] = []
    for row in mapped:
        sid = str(row["passon_sample_id"])
        label = labels_by_id[sid]
        out = dict(all_expression[sid])
        out.update({"sample_id": sid, "patient_id": str(label["patient_id"]), "response": _binary(label["response_binary"])})
        for key in ("timepoint", "treatment", "passon_day", "expression_source"):
            if key in label and str(label[key]).strip():
                out[key] = label[key]
        out["passon_day"] = label.get("passon_day", "")
        ordered.append(out)
    return ordered, {
        "expression_scale": "GSE115821_or_GSE168204_counts_then_log2p1",
        "expression_path": f"{_MGH_115_REL};{_MGH_168_REL}",
        "labels_path": _MGH_LABELS,
        "mapping_path": _MGH_MAP,
    }


def build_declared_cohort(
    truth_root: str | Path,
    cohort_id: str,
    *,
    freeze: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, str]]:
    """Build one generic A input from raw truth-repository files."""
    root = Path(truth_root).resolve()
    if cohort_id == "PRJEB23709":
        return _prjeb(root, freeze)
    if cohort_id == "MORRISON-1-public":
        return _morrison(root, freeze)
    if cohort_id == "MGH":
        return _mgh(root, freeze)
    raise FrozenInputError(f"unknown declared A cohort: {cohort_id}")


def write_generic_input(rows: Sequence[Mapping[str, Any]], path: str | Path, *, freeze: Mapping[str, Any]) -> Path:
    """Write an audit-local generic input table (never inside the repo)."""
    destination = Path(path)
    genes = sorted({
        _norm(gene)
        for geneset in freeze["method"]["axes"]["gene_sets"].values()
        for gene in geneset
    })
    fields = ["sample_id", "patient_id", "response", *genes]
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    return destination


def _metrics(records: Sequence[Mapping[str, Any]], freeze: Mapping[str, Any]) -> dict[str, Any]:
    labels = [int(row["response_binary"]) for row in records]
    scores = [float(row["boundary_score"]) for row in records]
    result: dict[str, Any] = {
        "analysis_unit": "sample",
        "n_records": len(records),
        "n_patients": len({str(row["patient_id"]) for row in records}),
        "labels": {"responder": sum(labels), "nonresponder": len(labels) - sum(labels)},
    }
    if len(set(labels)) < 2:
        result.update({"inference_status": "HOLD", "hold_reason": "both response classes are required"})
    else:
        stats = freeze["method"]["statistics"]
        result.update(
            registered_auc_inference(
                scores,
                labels,
                identifiers=[str(row["sample_id"]) for row in records],
                identifier_order="sample_id",
                seed=int(stats["seed"]),
                bootstrap_resamples=int(stats["bootstrap_resamples"]),
                permutation_resamples=int(stats["permutation_resamples"]),
                inference_contract=str(stats["inference_contract"]),
                rng_algorithm=str(stats["rng_algorithm"]),
                stream_policy=str(stats["stream_policy"]),
                permutation_comparator=str(stats["permutation_comparator"]),
                correction=str(stats["correction"]),
            )
        )
        result["inference_status"] = "COMPUTED_REGISTERED"
    return result


def _identity_sha256(records: Sequence[Mapping[str, Any]]) -> str:
    """Hash the ordered, identifier-free 85-row result identity.

    The identity-hash payload intentionally omits sample and patient
    identifiers. Cohort, ordinal, response, and full-precision score are
    sufficient to detect a changed assembled result; identifier-bearing rows
    remain confined to the external audit output.
    """
    payload = [
        {
            "cohort": str(row.get("_cohort", "")),
            "ordinal": index,
            "response": int(row["response_binary"]),
            "boundary_score": format(float(row["boundary_score"]), ".17g"),
        }
        for index, row in enumerate(records)
    ]
    encoded = json.dumps(payload, ensure_ascii=True, separators=(",", ":"), sort_keys=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _historical_cohort_name(value: object) -> str | None:
    """Normalize one legacy F2 cohort label to the declared A cohort key."""
    text = str(value).strip()
    for canonical, historical in _HISTORICAL_COHORT_NAMES.items():
        if text == historical:
            return historical
    return None


def _auc_or_none(scores: Sequence[float], labels: Sequence[int]) -> float | None:
    if not scores or len(set(labels)) < 2:
        return None
    return float(auc(list(scores), list(labels)))


def _drift_summary(
    canonical: Mapping[tuple[str, str], Mapping[str, Any]],
    historical: Mapping[tuple[str, str], float],
    keys: Sequence[tuple[str, str]],
) -> dict[str, Any]:
    matched = [key for key in keys if key in canonical and key in historical]
    missing = [key for key in keys if key in canonical and key not in historical]
    extra = [key for key in keys if key in historical and key not in canonical]
    canonical_scores = [float(canonical[key]["boundary_score"]) for key in matched]
    historical_scores = [float(historical[key]) for key in matched]
    labels = [int(canonical[key]["response_binary"]) for key in matched]
    absolute_drift = sorted(abs(canonical_score - historical_score) for canonical_score, historical_score in zip(canonical_scores, historical_scores))
    auc_historical = _auc_or_none(historical_scores, labels)
    auc_canonical = _auc_or_none(canonical_scores, labels)
    orientation_consistent = None
    orientation = None
    if auc_historical is not None and auc_canonical is not None:
        historical_higher = auc_historical >= 0.5
        canonical_higher = auc_canonical >= 0.5
        orientation_consistent = historical_higher == canonical_higher
        orientation = {
            "historical_higher_responder_like": historical_higher,
            "canonical_higher_responder_like": canonical_higher,
        }
    summary: dict[str, Any] = {
        "status": "COMPUTED" if not missing and not extra else "HOLD",
        "matched_count": len(matched),
        "missing_ids_count": len(missing),
        "extra_ids_count": len(extra),
        "max_abs_score_drift": max(absolute_drift) if absolute_drift else None,
        "mean_abs_score_drift": (sum(absolute_drift) / len(absolute_drift)) if absolute_drift else None,
        "median_abs_score_drift": (absolute_drift[len(absolute_drift) // 2] if len(absolute_drift) % 2 else (absolute_drift[len(absolute_drift) // 2 - 1] + absolute_drift[len(absolute_drift) // 2]) / 2) if absolute_drift else None,
        "auc_historical": auc_historical,
        "auc_canonical": auc_canonical,
        "delta_auc_canonical_minus_historical": (auc_canonical - auc_historical) if auc_historical is not None and auc_canonical is not None else None,
        "orientation_consistent": orientation_consistent,
        "orientation": orientation,
    }
    if missing or extra:
        summary["hold_reason"] = "canonical/historical sample-id mismatch; metrics are matched-only diagnostics"
    return summary


def compare_historical_to_canonical(
    records: Sequence[Mapping[str, Any]],
    historical_comparison_path: str | Path | None,
) -> dict[str, Any]:
    """Compare raw-derived A scores with legacy F2 scores after producer completion.

    This helper consumes only the legacy score column.  It never reads legacy
    labels, axes, or thresholds, and it does not feed any comparison value
    back into the A producer metrics.  IDs remain in memory for the join and
    only aggregate counts/statistics are returned.
    """
    if historical_comparison_path is None:
        return {"status": "NOT_SUPPLIED", "reason": "historical comparison object was not supplied"}
    path = Path(historical_comparison_path).resolve()
    if not path.exists():
        return {"status": "NOT_SUPPLIED", "reason": "historical comparison object is missing", "path": str(path)}
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except (OSError, csv.Error) as exc:
        return {"status": "HOLD", "reason": f"historical comparison object is unreadable: {exc}"}
    required = {"cohort", "sample_id", "boundary_score"}
    if not rows or not required.issubset(rows[0]):
        return {"status": "HOLD", "reason": f"historical comparison schema requires {sorted(required)}"}
    canonical: dict[tuple[str, str], Mapping[str, Any]] = {}
    canonical_duplicates = 0
    for row in records:
        cohort = _HISTORICAL_COHORT_NAMES.get(str(row.get("_cohort", "")).strip())
        sample_id = str(row.get("sample_id", "")).strip()
        if cohort is None or not sample_id:
            return {"status": "HOLD", "reason": "canonical comparison rows have an unmapped cohort or missing sample_id"}
        key = (cohort, sample_id)
        if key in canonical:
            canonical_duplicates += 1
        canonical[key] = row
    historical: dict[tuple[str, str], float] = {}
    historical_duplicates = 0
    unknown_cohort_rows = 0
    for row in rows:
        cohort = _historical_cohort_name(row.get("cohort"))
        sample_id = str(row.get("sample_id", "")).strip()
        if cohort is None:
            unknown_cohort_rows += 1
            continue
        if not sample_id:
            return {"status": "HOLD", "reason": "historical comparison row is missing sample_id"}
        try:
            score = float(row["boundary_score"])
        except (TypeError, ValueError):
            return {"status": "HOLD", "reason": "historical comparison boundary_score is not numeric"}
        if not math.isfinite(score):
            return {"status": "HOLD", "reason": "historical comparison boundary_score is non-finite"}
        key = (cohort, sample_id)
        if key in historical:
            historical_duplicates += 1
        historical[key] = score
    if canonical_duplicates or historical_duplicates or unknown_cohort_rows:
        return {
            "status": "HOLD",
            "reason": "historical/canonical comparison mapping conflict",
            "mapping": {
                "canonical_rows": len(records),
                "historical_rows": len(rows),
                "canonical_duplicate_keys": canonical_duplicates,
                "historical_duplicate_keys": historical_duplicates,
                "unknown_historical_cohort_rows": unknown_cohort_rows,
            },
        }
    cohort_summaries: dict[str, Any] = {}
    for cohort in _HISTORICAL_COHORT_NAMES.values():
        keys = sorted({key for key in canonical if key[0] == cohort} | {key for key in historical if key[0] == cohort})
        cohort_summaries[cohort] = _drift_summary(canonical, historical, keys)
    pooled = _drift_summary(canonical, historical, sorted(set(canonical) | set(historical)))
    mismatch = any(item["missing_ids_count"] or item["extra_ids_count"] for item in cohort_summaries.values())
    return {
        "status": "HOLD" if mismatch else "COMPUTED",
        "reason": "canonical/historical sample-id mismatch" if mismatch else None,
        "comparison_source": str(path),
        "comparison_source_sha256": _sha256(path),
        "mapping": {
            "canonical_rows": len(records),
            "historical_rows": len(rows),
            "canonical_unique_keys": len(canonical),
            "historical_unique_keys": len(historical),
            "unknown_historical_cohort_rows": 0,
        },
        "cohorts": cohort_summaries,
        "pooled": pooled,
    }


def mgh_earliest_patient_records(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Frozen MGH patient sensitivity: earliest eligible mapped sample per patient."""
    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for row in records:
        grouped.setdefault(str(row["patient_id"]), []).append(row)
    selected: list[dict[str, Any]] = []
    for patient_id, items in grouped.items():
        with_day = [(float(day), index, row) for index, row in enumerate(items) if (day := _float(row.get("passon_day"))) is not None]
        if with_day:
            row = min(with_day, key=lambda item: (item[0], item[1]))[2]
            selection_rule = "minimum_passon_day_then_source_order"
        else:
            row = items[0]
            selection_rule = "source_order_when_passon_day_missing"
        selected.append(
            {
                "sample_id": row["sample_id"],
                "patient_id": patient_id,
                "boundary_score": float(row["boundary_score"]),
                "response_binary": int(row["response_binary"]),
                "passon_day": row.get("passon_day"),
                "selection_rule": selection_rule,
            }
        )
    def patient_sort(row: Mapping[str, Any]) -> tuple[int, int | str, str]:
        identifier = str(row["patient_id"])
        try:
            parsed = int(identifier)
        except ValueError:
            return (1, identifier, identifier)
        return (0, parsed, identifier) if str(parsed) == identifier else (1, identifier, identifier)
    return sorted(selected, key=patient_sort)


def _gse91061_soft_bytes(root: Path, soft_path: Path | None = None) -> bytes:
    """Read the raw SOFT locally, or the pinned historical blob, read-only."""
    path = soft_path or (root / _GSE91061_SOFT)
    if path.exists():
        return path.read_bytes()
    # The current truth snapshot omitted the SOFT but its pinned git object is
    # still available.  ``git show`` is read-only and avoids manufacturing a
    # second source file in the truth repository.
    try:
        return subprocess.run(
            ["git", "show", "dcd33df:raw_data/GSE91061_family.soft.gz"],
            cwd=str(root), check=True, stdout=subprocess.PIPE,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as exc:
        raise FrozenInputError(f"GSE91061 raw SOFT is unavailable: {exc}") from exc


def read_gse91061_on_treatment_labels(root: str | Path, soft_path: str | Path | None = None) -> list[dict[str, Any]]:
    """Parse raw SOFT sample characteristics into the declared 27 labels."""
    import gzip as _gzip

    data = _gse91061_soft_bytes(Path(root).resolve(), Path(soft_path) if soft_path else None)
    try:
        text = _gzip.decompress(data).decode("utf-8", "replace")
    except OSError as exc:
        raise FrozenInputError(f"GSE91061 SOFT is not valid gzip: {exc}") from exc
    records: list[dict[str, Any]] = []
    current: dict[str, str] = {}
    for line in text.splitlines() + [""]:  # final flush
        if line.startswith("!Sample_title"):
            if current:
                records.append(current)
            current = {"sample_id": line.split("=", 1)[1].strip()}
        elif line.startswith("!Sample_characteristics_ch1"):
            value = line.split("=", 1)[1].strip()
            key, _, val = value.partition(":")
            current[key.strip().lower()] = val.strip()
        elif line.startswith("!Sample_") is False and not line.strip() and current:
            records.append(current)
            current = {}
    # The canonical derivation population is the binary on-treatment subset
    # with a corresponding pretreatment biopsy (27 patients: 9 PRCR/18 PD).
    paired_patients = {
        str(row.get("sample_id", "")).split("_", 1)[0]
        for row in records
        if str(row.get("visit (pre or on treatment)", "")).lower() == "pre"
    }
    out: list[dict[str, Any]] = []
    for row in records:
        visit = str(row.get("visit (pre or on treatment)", "")).lower()
        response = str(row.get("response", "")).upper()
        patient_token = str(row.get("sample_id", "")).split("_", 1)[0]
        if visit != "on" or response not in {"PD", "PRCR"} or patient_token not in paired_patients:
            continue
        match = re.match(r"Pt([^_]+)_", str(row.get("sample_id", "")))
        if not match:
            raise FrozenInputError(f"GSE91061 patient ID cannot be parsed: {row.get('sample_id')}")
        out.append({
            "sample_id": row["sample_id"],
            "patient_id": f"Pt{match.group(1)}",
            "response": 1 if response == "PRCR" else 0,
            "response_label": response,
        })
    if len(out) != 27 or sum(int(row["response"]) for row in out) != 9:
        raise FrozenInputError(f"GSE91061 SOFT yielded {len(out)} on-treatment PD/PRCR rows; expected 27 (9R/18NR)")
    return out


def _axis_raw_from_expression(rows: Sequence[Mapping[str, float]], freeze: Mapping[str, Any]) -> list[dict[str, float]]:
    axes = freeze["method"]["axes"]["gene_sets"]
    raw: list[dict[str, float]] = []
    for row in rows:
        axis_row: dict[str, float] = {}
        for axis in AXIS_ORDER:
            genes = [_norm(gene) for gene in axes[axis] if _norm(gene) in row]
            if len(genes) / len(axes[axis]) < float(freeze["method"]["preprocessing"]["minimum_gene_coverage"]):
                raise FrozenInputError(f"GSE91061 axis {axis} coverage failed")
            axis_row[axis] = sum(float(row[gene]) for gene in genes) / len(genes)
        raw.append(axis_row)
    return raw


def _fold_z(train: Sequence[Mapping[str, float]], test: Sequence[Mapping[str, float]]) -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    means: dict[str, float] = {}
    scales: dict[str, float] = {}
    for axis in AXIS_ORDER:
        values = [float(row[axis]) for row in train]
        mean = sum(values) / len(values)
        scale = math.sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - 1))
        if not math.isfinite(scale) or scale == 0:
            raise FrozenInputError(f"GSE91061 fold axis {axis} has invalid ddof=1 scale")
        means[axis], scales[axis] = mean, scale
    def transform(rows: Sequence[Mapping[str, float]]) -> list[dict[str, float]]:
        return [{axis: (float(row[axis]) - means[axis]) / scales[axis] for axis in AXIS_ORDER} for row in rows]
    return transform(train), transform(test)


def _fit_sign_constrained(scores: Sequence[Mapping[str, float]], labels: Sequence[int]) -> tuple[dict[str, float], float]:
    try:
        import numpy as np
        from scipy.optimize import minimize
    except Exception as exc:  # pragma: no cover - dependency-specific hold
        raise FrozenInputError(f"sign-constrained LOPO solver unavailable: {exc}") from exc
    matrix = np.asarray([[float(row[axis]) for axis in AXIS_ORDER] for row in scores], dtype=float)
    y = np.asarray(labels, dtype=int)
    signed = 2 * y - 1
    def objective(params):
        beta, intercept = params[:-1], params[-1]
        logits = matrix @ beta + intercept
        return float(np.logaddexp(0.0, -signed * logits).sum() + 0.5 * np.dot(beta, beta))
    result = minimize(objective, np.zeros(5), method="L-BFGS-B", bounds=[(None, 0.0), (0.0, None), (None, 0.0), (0.0, None), (None, None)])
    if not result.success:
        raise FrozenInputError(f"sign-constrained LOPO fit failed: {result.message}")
    return {axis: float(value) for axis, value in zip(AXIS_ORDER, result.x[:-1])}, float(result.x[-1])


def _gse91061_historical_inference(scores: Sequence[float], labels: Sequence[int]) -> dict[str, Any]:
    """The locked internal V1 inference (NumPy global RNG, B=2000/1000).

    This is intentionally separate from the external 10,000-resample
    registered inference.  The historical source workbook records the
    one-sided upper-tail permutation count as ``permutation_p``; the
    contract-compliant absolute-tail value is retained alongside it.
    """
    import numpy as np

    y = np.asarray(labels, dtype=int)
    s = np.asarray(scores, dtype=float)
    np.random.seed(42)
    bootstrap: list[float] = []
    for _ in range(2000):
        indices = np.random.choice(len(y), len(y), replace=True)
        if len(np.unique(y[indices])) < 2:
            continue
        bootstrap.append(auc(s[indices].tolist(), y[indices].tolist()))
    observed = auc(scores, labels)
    permuted: list[float] = []
    for _ in range(1000):
        permuted.append(auc(scores, np.random.permutation(y).tolist()))
    arr = np.asarray(permuted, dtype=float)
    p_abs = float((np.sum(np.abs(arr - 0.5) >= abs(observed - 0.5)) + 1) / 1001)
    p_source = float((np.sum(arr >= observed) + 1) / 1001)
    return {
        "auc": observed,
        "ci_95": (float(np.percentile(bootstrap, 2.5)), float(np.percentile(bootstrap, 97.5))),
        # Historical internal source field (the retained canonical number).
        "p_two_sided": p_source,
        # Contract-compliant absolute-tail diagnostic retained explicitly.
        "p_absolute_two_sided": p_abs,
        "permutation_p_source_upper_tail": p_source,
        "bootstrap_resamples": 2000,
        "permutation_resamples": 1000,
        "seed": 42,
        "inference_contract": "GSE91061_INTERNAL_V1",
    }


def run_gse91061_lopo(*, truth_root: str | Path, config_path: str | Path) -> dict[str, Any]:
    """Recompute frozen LOO and sign-constrained refit LOPO from raw labels."""
    from .config import load_freeze

    freeze = load_freeze(Path(config_path).resolve().parents[1])
    root = Path(truth_root).resolve()
    try:
        labels = read_gse91061_on_treatment_labels(root)
        expr_path = _require(root / _GSE91061_EXPR)
        # Preserve the historical matrix-column order.  Besides making the
        # fold ledger deterministic, this reproduces the NumPy RNG ordering
        # used by the archived internal bootstrap/permutation procedure.
        with expr_path.open(encoding="utf-8-sig", newline="") as handle:
            expression_order = [field for field in next(csv.reader(handle)) if field]
        order = {sample_id: index for index, sample_id in enumerate(expression_order)}
        labels.sort(key=lambda row: order.get(str(row["sample_id"]), len(order)))
        expression = _read_selected_matrix_fast(expr_path, [row["sample_id"] for row in labels], delimiter=",", gene_field=None, transform="as_deposited_logcpm", genes=_gene_set_union(freeze))
        raw = _axis_raw_from_expression(expression, freeze)
        patients = [str(row["patient_id"]) for row in labels]
        frozen_scores: list[float] = [0.0] * len(labels)
        refit_scores: list[float] = [0.0] * len(labels)
        signs: list[dict[str, int]] = []
        for holdout in sorted(set(patients), key=lambda value: int(re.sub(r"\D", "", value) or 0)):
            train_idx = [index for index, patient in enumerate(patients) if patient != holdout]
            test_idx = [index for index, patient in enumerate(patients) if patient == holdout]
            train_z, test_z = _fold_z([raw[index] for index in train_idx], [raw[index] for index in test_idx])
            coef = freeze["boundary"]["coefficients"]
            fold_frozen = [float(freeze["boundary"]["intercept"]) + sum(float(coef[axis]) * row[axis] for axis in AXIS_ORDER) for row in test_z]
            for index, score in zip(test_idx, fold_frozen):
                frozen_scores[index] = score
            fit_beta, fit_intercept = _fit_sign_constrained(train_z, [int(labels[index]["response"]) for index in train_idx])
            fold_refit = [fit_intercept + sum(fit_beta[axis] * row[axis] for axis in AXIS_ORDER) for row in test_z]
            for index, score in zip(test_idx, fold_refit):
                refit_scores[index] = score
            signs.append({axis: (1 if fit_beta[axis] > 1e-12 else -1 if fit_beta[axis] < -1e-12 else 0) for axis in AXIS_ORDER})
        y = [int(row["response"]) for row in labels]
        frozen_inference = _gse91061_historical_inference(frozen_scores, y)
        refit_inference = _gse91061_historical_inference(refit_scores, y)
        sign_stability = {axis: sum(1 for row in signs if row[axis] == (1 if float(freeze["boundary"]["coefficients"][axis]) > 0 else -1 if float(freeze["boundary"]["coefficients"][axis]) < 0 else 0)) / len(signs) for axis in AXIS_ORDER}
        return {"status": "COMPUTED", "n_records": len(labels), "labels": {"responder": sum(y), "nonresponder": len(y) - sum(y)}, "frozen_loo": {"analysis_unit": "biopsy", **frozen_inference}, "refit_lopo": {"analysis_unit": "biopsy", **refit_inference}, "sign_stability": sign_stability, "provenance": {"expression_path": _GSE91061_EXPR, "label_source": "raw_data/GSE91061_family.soft.gz (git dcd33df fallback)", "fit_scaling": "training-fold axis mean/std ddof=1"}}
    except (FrozenInputError, OSError, ValueError, KeyError) as exc:
        return {"status": "HOLD", "hold_reason": str(exc), "provenance": {"expression_path": _GSE91061_EXPR, "label_source": "raw SOFT only; comparison tables excluded"}}


def run_truth_cohorts(
    *,
    truth_root: str | Path,
    audit_dir: str | Path,
    config_path: str | Path,
    historical_comparison: str | Path | None = None,
) -> dict[str, Any]:
    """Run A from raw truth inputs and write all identifier-bearing output externally."""
    from .config import load_freeze

    freeze = load_freeze(Path(config_path).resolve().parents[1])
    audit = Path(audit_dir).resolve()
    audit.mkdir(parents=True, exist_ok=True)
    result: dict[str, Any] = {"producer": "A_FROZEN_STATE_TRANSFER", "status": "HOLD", "cohorts": {}, "pooled": None, "identity": {"n_records": 0, "sha256": None}, "drift": {"status": "NOT_COMPUTED", "reason": "comparison-only references are not producer inputs"}, "historical_to_canonical_drift": {"status": "NOT_SUPPLIED", "reason": "historical comparison object was not supplied"}}
    all_records: list[dict[str, Any]] = []
    for cohort_id in ("PRJEB23709", "MORRISON-1-public", "MGH"):
        try:
            rows, provenance = build_declared_cohort(truth_root, cohort_id, freeze=freeze)
            input_path = write_generic_input(rows, audit / f"A_{cohort_id.replace('-', '_')}.tsv", freeze=freeze)
            payload = run_frozen_state_transfer(config_path=config_path, input_path=input_path, output_path=audit / f"A_{cohort_id.replace('-', '_')}.json", input_representation="gene_expression", analysis_unit="sample", cohort_id=cohort_id)
            if payload.get("status") != "COMPUTED":
                result["cohorts"][cohort_id] = payload
                continue
            records = []
            for row, score in zip(rows, payload["boundary_scores"]):
                records.append({**row, "_cohort": cohort_id, "boundary_score": float(score["boundary_score"]), "response_binary": int(row["response"])})
            all_records.extend(records)
            result["cohorts"][cohort_id] = {"status": "COMPUTED", "n": len(records), "metrics": _metrics(records, freeze), "provenance": {**provenance, "input_sha256": _sha256(input_path)}}
            if cohort_id == "MGH":
                patient = mgh_earliest_patient_records(records)
                labels = [int(row["response_binary"]) for row in patient]
                scores = [float(row["boundary_score"]) for row in patient]
                result["cohorts"][cohort_id]["patient_earliest"] = {
                    "analysis_unit": "patient", "n_records": len(patient), "labels": {"responder": sum(labels), "nonresponder": len(labels) - sum(labels)},
                    "rows": patient,
                    **(registered_auc_inference(
                        scores,
                        labels,
                        identifiers=[str(row["patient_id"]) for row in patient],
                        identifier_order="patient_id",
                        seed=int(freeze["method"]["statistics"]["seed"]),
                        bootstrap_resamples=int(freeze["method"]["statistics"]["bootstrap_resamples"]),
                        permutation_resamples=int(freeze["method"]["statistics"]["permutation_resamples"]),
                        inference_contract=str(freeze["method"]["statistics"]["inference_contract"]),
                        rng_algorithm=str(freeze["method"]["statistics"]["rng_algorithm"]),
                        stream_policy=str(freeze["method"]["statistics"]["stream_policy"]),
                        permutation_comparator=str(freeze["method"]["statistics"]["permutation_comparator"]),
                        correction=str(freeze["method"]["statistics"]["correction"]),
                    ) if len(set(labels)) > 1 else {"inference_status": "HOLD"}),
                }
        except (FrozenInputError, OSError, ValueError, KeyError) as exc:
            result["cohorts"][cohort_id] = {"status": "HOLD", "hold_reason": str(exc)}
    if len(all_records) == 85 and len({int(row["response_binary"]) for row in all_records}) == 2:
        pooled_labels = [int(row["response_binary"]) for row in all_records]
        pooled_scores = [float(row["boundary_score"]) for row in all_records]
        stats = freeze["method"]["statistics"]
        pooled_inference = registered_auc_inference(
            pooled_scores,
            pooled_labels,
            # Pool identifiers are cohort-qualified before canonical sorting;
            # a raw sample ID must never collide across external cohorts.
            identifiers=[f"{row['_cohort']}::{row['sample_id']}" for row in all_records],
            identifier_order="sample_id",
            seed=int(stats["seed"]),
            bootstrap_resamples=int(stats["bootstrap_resamples"]),
            permutation_resamples=int(stats["permutation_resamples"]),
            inference_contract=str(stats["inference_contract"]),
            rng_algorithm=str(stats["rng_algorithm"]),
            stream_policy=str(stats["stream_policy"]),
            permutation_comparator=str(stats["permutation_comparator"]),
            correction=str(stats["correction"]),
        )
        result["pooled"] = {
            "analysis_unit": "sample",
            "n_records": len(all_records),
            "labels": {
                "responder": sum(pooled_labels),
                "nonresponder": len(pooled_labels) - sum(pooled_labels),
            },
            **pooled_inference,
        }
        result["identity"] = {"n_records": len(all_records), "sha256": _identity_sha256(all_records)}
        result["status"] = "COMPUTED"
    else:
        result["pooled"] = {"status": "HOLD", "reason": f"expected 85 rows with both classes; got {len(all_records)}"}
        result["identity"] = {"n_records": len(all_records), "sha256": _identity_sha256(all_records) if all_records else None}
    if result["status"] == "COMPUTED":
        result["historical_to_canonical_drift"] = compare_historical_to_canonical(all_records, historical_comparison)
    elif historical_comparison is not None:
        result["historical_to_canonical_drift"] = {"status": "HOLD", "reason": "raw-derived producer output is incomplete; comparison was not run"}
    result["gse91061_stability"] = run_gse91061_lopo(truth_root=truth_root, config_path=config_path)
    (audit / "A_85_NUMERIC_OUTPUT.json").write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return result


def discover_gse91061_inputs(truth_root: str | Path) -> dict[str, Any]:
    """Report raw derivation input availability without using comparison objects."""
    root = Path(truth_root).resolve()
    # The declared A contract consumes the processed symbol-level matrix;
    # raw-study files remain acquisition provenance and are not producer input.
    raw = root / _GSE91061_EXPR
    labels_candidates = [root / "metadata/v2_external/gse91061_on_treatment_labels.csv", root / "metadata/GSE91061/labels.tsv"]
    fold_candidates = [root / "data/external/GSE91061/folds.tsv", root / "metadata/GSE91061/folds.tsv"]
    labels = next((path for path in labels_candidates if path.exists()), None)
    folds = next((path for path in fold_candidates if path.exists()), None)
    git_soft_available = False
    try:
        subprocess.run(["git", "cat-file", "-e", "dcd33df:raw_data/GSE91061_family.soft.gz"], cwd=str(root), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        git_soft_available = True
    except (OSError, subprocess.CalledProcessError):
        pass
    labels_available = labels is not None or git_soft_available
    fold_available = folds is not None or git_soft_available
    return {
        "raw_expression": {"path": str(raw), "exists": raw.exists(), "sha256": _sha256(raw) if raw.exists() else None},
        "labels": {"path": str(labels) if labels else ("git:dcd33df:raw_data/GSE91061_family.soft.gz" if git_soft_available else None), "exists": labels_available},
        "canonical_ref_fold_inputs": {"path": str(folds) if folds else ("raw SOFT + code_release/src/boundary.py + fold-safe runner" if git_soft_available else None), "exists": fold_available},
        "status": "READY" if raw.exists() and labels_available and fold_available else "HOLD",
        "hold_reason": None if raw.exists() and labels_available and fold_available else "unique canonical GSE91061 raw labels/reference-fold dependencies are unavailable; comparison-only score tables are excluded",
    }
