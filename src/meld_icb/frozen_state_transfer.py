"""The sole frozen-state transfer producer and its explicit data contracts."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .config import load_freeze
from .statistics import auc, registered_auc_inference

AXIS_ORDER = ("T", "E", "X", "C")
REQUIRED_METADATA = ("sample_id", "patient_id", "response")
VALID_REPRESENTATIONS = {"gene_expression", "canonical_axis_z"}
VALID_ANALYSIS_UNITS = {"sample", "patient"}
REPO_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class AxisNormalization:
    means: dict[str, float]
    standard_deviations: dict[str, float]
    ddof: int
    fitted_sample_ids: tuple[str, ...]
    normalization_scope_id: str


def normalize_gene_id(value: object) -> str:
    return str(value).strip().upper()


def collapse_duplicate_symbols(records: Iterable[tuple[object, object]]) -> dict[str, float]:
    grouped: dict[str, list[float]] = defaultdict(list)
    for raw_gene, raw_value in records:
        gene = normalize_gene_id(raw_gene)
        if not gene:
            continue
        value = float(raw_value)
        if not math.isfinite(value):
            raise ValueError(f"non-finite expression for {gene}")
        grouped[gene].append(value)
    if not grouped:
        raise ValueError("no gene records")
    return {gene: sum(values) / len(values) for gene, values in grouped.items()}


def _axis_raw(samples: Sequence[Mapping[str, float]], gene_sets: Mapping[str, Sequence[str]], minimum_coverage: float) -> tuple[list[dict[str, float]], dict[str, dict[str, Any]]]:
    raw = [dict() for _ in samples]
    coverage: dict[str, dict[str, Any]] = {}
    for axis in AXIS_ORDER:
        genes = tuple(normalize_gene_id(gene) for gene in gene_sets[axis])
        present = tuple(gene for gene in genes if any(gene in sample for sample in samples))
        fraction = len(present) / len(genes) if genes else 0.0
        if fraction < minimum_coverage:
            raise ValueError(f"axis {axis} minimum coverage failed: {fraction:.3f}")
        coverage[axis] = {"declared_genes": len(genes), "observed_genes": len(present), "coverage": fraction, "present_genes": list(present)}
        for index, sample in enumerate(samples):
            available = [gene for gene in present if gene in sample and sample[gene] is not None]
            if not available:
                raise ValueError(f"axis {axis} has no available genes in sample {index}")
            raw[index][axis] = sum(float(sample[gene]) for gene in available) / len(available)
    return raw, coverage


def fit_axis_normalization(raw: Sequence[Mapping[str, float]], sample_ids: Sequence[str], *, scope_id: str, ddof: int = 1) -> AxisNormalization:
    if ddof != 1:
        raise ValueError("canonical axis normalization requires ddof=1")
    if not raw or len(raw) != len(sample_ids):
        raise ValueError("normalization fit requires aligned samples and IDs")
    means: dict[str, float] = {}
    standard_deviations: dict[str, float] = {}
    for axis in AXIS_ORDER:
        values = [float(row[axis]) for row in raw]
        mean = sum(values) / len(values)
        sd = math.sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - ddof))
        if not math.isfinite(mean) or not math.isfinite(sd) or sd == 0:
            raise ValueError(f"axis {axis} has invalid normalization scale")
        means[axis] = mean
        standard_deviations[axis] = sd
    return AxisNormalization(means, standard_deviations, ddof, tuple(str(item) for item in sample_ids), scope_id)


def transform_axis_normalization(raw: Sequence[Mapping[str, float]], parameters: AxisNormalization) -> list[dict[str, float]]:
    return [{axis: (float(row[axis]) - parameters.means[axis]) / parameters.standard_deviations[axis] for axis in AXIS_ORDER} for row in raw]


def pre_mean_then_axis_z(samples: Sequence[Mapping[object, object]], freeze: Mapping[str, Any], sample_ids: Sequence[str], *, scope_id: str) -> tuple[list[dict[str, float]], AxisNormalization, dict[str, Any]]:
    method = freeze["method"]
    normalized = [collapse_duplicate_symbols(sample.items()) for sample in samples]
    raw, coverage = _axis_raw(normalized, method["axes"]["gene_sets"], float(method["preprocessing"]["minimum_gene_coverage"]))
    parameters = fit_axis_normalization(raw, sample_ids, scope_id=scope_id, ddof=int(method["preprocessing"]["ddof"]))
    return transform_axis_normalization(raw, parameters), parameters, coverage


def pre_mean_then_axis_z_transfer(reference_samples: Sequence[Mapping[object, object]], target_samples: Sequence[Mapping[object, object]], freeze: Mapping[str, Any], reference_ids: Sequence[str], *, scope_id: str) -> tuple[list[dict[str, float]], AxisNormalization, dict[str, Any]]:
    method = freeze["method"]
    ref = [collapse_duplicate_symbols(sample.items()) for sample in reference_samples]
    target = [collapse_duplicate_symbols(sample.items()) for sample in target_samples]
    ref_raw, ref_coverage = _axis_raw(ref, method["axes"]["gene_sets"], float(method["preprocessing"]["minimum_gene_coverage"]))
    target_raw, target_coverage = _axis_raw(target, method["axes"]["gene_sets"], float(method["preprocessing"]["minimum_gene_coverage"]))
    parameters = fit_axis_normalization(ref_raw, reference_ids, scope_id=scope_id, ddof=int(method["preprocessing"]["ddof"]))
    return transform_axis_normalization(target_raw, parameters), parameters, {"reference": ref_coverage, "target": target_coverage}


def frozen_boundary_score(axes: Mapping[str, float], freeze: Mapping[str, Any]) -> float:
    boundary = freeze["boundary"]
    missing = [axis for axis in AXIS_ORDER if axis not in axes]
    if missing:
        raise ValueError(f"missing canonical axes: {missing}")
    return float(boundary["intercept"]) + sum(float(boundary["coefficients"][axis]) * float(axes[axis]) for axis in AXIS_ORDER)


def score_axes(axes: Mapping[str, float], boundary: Mapping[str, Any]) -> float:
    """Compatibility façade; scientific callers must obtain ``boundary`` from load_freeze."""
    missing = [axis for axis in AXIS_ORDER if axis not in axes]
    if missing:
        raise ValueError(f"missing canonical axes: {missing}")
    return float(boundary["intercept"]) + sum(float(boundary["coefficients"][axis]) * float(axes[axis]) for axis in AXIS_ORDER)


def classify_boundary(score: float, freeze: Mapping[str, Any]) -> int:
    return int(float(score) >= float(freeze["boundary"]["threshold"]))


def _encode_response(value: object) -> int:
    normalized = str(value).strip().upper()
    if normalized in {"1", "R", "RESPONDER", "PR", "CR", "PRCR"}:
        return 1
    if normalized in {"0", "NR", "NONRESPONDER", "PD"}:
        return 0
    raise ValueError(f"unsupported binary response label: {value!r}")


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        # Wide gene-expression rows can make ``Sniffer`` ambiguous (the first
        # 4 KiB may contain no complete record).  The declared adapter writes
        # TSV, so use the extension as an unambiguous contract and retain
        # Sniffer only for generic legacy inputs.
        if path.suffix.lower() in {".tsv", ".tab"}:
            dialect = csv.excel_tab
        elif path.suffix.lower() == ".csv":
            dialect = csv.excel
        else:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,") if sample.strip() else csv.excel_tab
        return [dict(row) for row in csv.DictReader(handle, dialect=dialect)]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _logical_path(path: Path) -> str:
    """Return a repository-relative identifier without leaking local roots."""
    try:
        return path.resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return path.name


def _validate_row(row: Mapping[str, str], representation: str) -> None:
    if representation not in VALID_REPRESENTATIONS:
        raise ValueError(f"input representation must be explicit: {VALID_REPRESENTATIONS}")
    missing = [field for field in REQUIRED_METADATA if not str(row.get(field, "")).strip()]
    if missing:
        raise ValueError(f"missing required metadata fields: {missing}")
    if representation == "canonical_axis_z":
        missing_axes = [axis for axis in AXIS_ORDER if not str(row.get(axis, "")).strip()]
        missing_provenance = [field for field in ("normalization_scope_id", "normalization_source") if not str(row.get(field, "")).strip()]
        if missing_axes or missing_provenance:
            raise ValueError(f"canonical_axis_z requires axes and normalization provenance: {missing_axes + missing_provenance}")


def _sample_expression(row: Mapping[str, str]) -> dict[str, float]:
    ignored = set(REQUIRED_METADATA) | {"timepoint", "therapy", "cohort", "pair_id", "analysis_unit", "normalization_scope_id", "normalization_source"}
    values: dict[str, float] = {}
    for key, value in row.items():
        if key in ignored or str(value).strip() == "":
            continue
        try:
            values[key] = float(value)
        except (TypeError, ValueError):
            continue
    if not values:
        raise ValueError(f"sample {row.get('sample_id', '<unknown>')} has no numeric expression")
    return values


def aggregate_records(records: Sequence[Mapping[str, Any]], analysis_unit: str) -> list[dict[str, Any]]:
    if analysis_unit not in VALID_ANALYSIS_UNITS:
        raise ValueError(f"analysis unit must be explicit: {VALID_ANALYSIS_UNITS}")
    if analysis_unit == "sample":
        return [dict(record) for record in records]
    grouped: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for record in records:
        grouped[str(record["patient_id"])].append(record)
    aggregated = []
    for patient_id, items in sorted(grouped.items()):
        score = sum(float(item["boundary_score"]) for item in items) / len(items)
        labels = {int(item["response_binary"]) for item in items}
        if len(labels) != 1:
            raise ValueError(f"patient {patient_id} has conflicting response labels")
        aggregated.append({"sample_id": f"patient:{patient_id}", "patient_id": patient_id, "boundary_score": score, "response_binary": labels.pop()})
    return aggregated


def _metrics(records: Sequence[Mapping[str, Any]], analysis_unit: str, freeze: Mapping[str, Any]) -> dict[str, Any]:
    aggregated = aggregate_records(records, analysis_unit)
    labels = [int(item["response_binary"]) for item in aggregated]
    scores = [float(item["boundary_score"]) for item in aggregated]
    result: dict[str, Any] = {"analysis_unit": analysis_unit, "n_records": len(aggregated), "n_patients": len({str(item["patient_id"]) for item in records}), "labels": {"responder": sum(labels), "nonresponder": len(labels) - sum(labels)}, "auc": None, "ci_95": None, "p_two_sided": None, "inference_status": "HOLD"}
    if len(set(labels)) >= 2:
        stats = freeze["method"]["statistics"]
        inference = registered_auc_inference(
            scores,
            labels,
            identifiers=[str(item["sample_id"] if analysis_unit == "sample" else item["patient_id"]) for item in aggregated],
            identifier_order=str(stats["identifier_order"][analysis_unit]),
            seed=int(stats["seed"]),
            bootstrap_resamples=int(stats["bootstrap_resamples"]),
            permutation_resamples=int(stats["permutation_resamples"]),
            inference_contract=str(stats["inference_contract"]),
            rng_algorithm=str(stats["rng_algorithm"]),
            stream_policy=str(stats["stream_policy"]),
            permutation_comparator=str(stats["permutation_comparator"]),
            correction=str(stats["correction"]),
        )
        result.update(inference)
        result["inference_status"] = "COMPUTED_REGISTERED"
    else:
        result["inference_hold_reason"] = "both response classes are required"
    return result


def _historical_comparison(path: Path | None, canonical_rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    if path is None or not path.exists():
        return {"status": "HOLD", "reason": "historical/canonical reference object is missing"}
    try:
        reference = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return {"status": "HOLD", "reason": f"reference schema unreadable: {exc}"}
    if not isinstance(reference, dict) or "canonical_axes" not in reference:
        return {"status": "HOLD", "reason": "reference lacks canonical_axes"}
    expected = len(reference["canonical_axes"])
    actual = len(canonical_rows)
    return {"status": "COMPARABLE" if expected == actual else "MISMATCH", "expected_rows": expected, "actual_rows": actual, "row_tolerance": 0}


def run_frozen_state_transfer(*, config_path: str | Path, input_path: str | Path, output_path: str | Path, input_representation: str | None, analysis_unit: str | None, cohort_id: str = "UNDECLARED", reference_input_path: str | Path | None = None, reference_output_path: str | Path | None = None) -> dict[str, Any]:
    """Execute one declared cohort; missing data, representation, or reference remains HOLD."""
    freeze = load_freeze(Path(config_path).resolve().parents[1])
    source = Path(input_path)
    destination = Path(output_path)
    metadata: dict[str, Any] = {"producer": "A_FROZEN_STATE_TRANSFER", "cohort_id": cohort_id, "config_path": _logical_path(Path(config_path)), "config_sha256": _sha256(Path(config_path)) if Path(config_path).exists() else None, "input_path": _logical_path(source), "input_representation": input_representation, "analysis_unit": analysis_unit}
    if not source.exists():
        payload = {"status": "HOLD", "hold_reason": "required cohort input file is missing", "cohort_id": cohort_id, "metadata": metadata}
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload
    try:
        if input_representation not in VALID_REPRESENTATIONS:
            raise ValueError("input representation must be explicit; automatic axis detection is disabled")
        if analysis_unit not in VALID_ANALYSIS_UNITS:
            raise ValueError("analysis unit must be explicit")
        rows = _read_rows(source)
        for row in rows:
            _validate_row(row, input_representation)
        reference_rows = _read_rows(Path(reference_input_path)) if reference_input_path else None
        if reference_input_path and not Path(reference_input_path).exists():
            raise ValueError("external transfer reference input is missing")
        for row in reference_rows or []:
            _validate_row(row, "gene_expression")
        if input_representation == "canonical_axis_z":
            axes = [{axis: float(row[axis]) for axis in AXIS_ORDER} for row in rows]
            normalization = {"scope": "declared_input", "source": rows[0]["normalization_source"], "scope_id": rows[0]["normalization_scope_id"], "ddof": 1}
            coverage = {"status": "DECLARED_CANONICAL_AXIS_INPUT"}
        else:
            samples = [_sample_expression(row) for row in rows]
            ids = [row["sample_id"] for row in rows]
            if reference_rows:
                axes, parameters, coverage = pre_mean_then_axis_z_transfer([_sample_expression(row) for row in reference_rows], samples, freeze, [row["sample_id"] for row in reference_rows], scope_id=f"{cohort_id}:reference")
            else:
                axes, parameters, coverage = pre_mean_then_axis_z(samples, freeze, ids, scope_id=f"{cohort_id}:cohort")
            normalization = asdict(parameters)
        scored = []
        for row, axis in zip(rows, axes):
            score = frozen_boundary_score(axis, freeze)
            scored.append({"sample_id": row["sample_id"], "patient_id": row["patient_id"], **axis, "boundary_score": score, "boundary_label": classify_boundary(score, freeze), "response_binary": _encode_response(row["response"])})
        metrics = _metrics(scored, analysis_unit, freeze)
        payload = {"status": "COMPUTED", "cohort_id": cohort_id, "canonical_axes": [{"sample_id": row["sample_id"], "patient_id": row["patient_id"], **{axis: row[axis] for axis in AXIS_ORDER}} for row in scored], "boundary_scores": [{"sample_id": row["sample_id"], "patient_id": row["patient_id"], "boundary_score": row["boundary_score"], "boundary_label": row["boundary_label"]} for row in scored], "cohort_metrics": metrics, "historical_to_canonical": _historical_comparison(Path(reference_output_path) if reference_output_path else None, scored), "stability": {"status": "HOLD", "reason": "canonical F1_c coefficient sign-stability component and fold inputs are not available; omitted-row AUC is not a substitute"}, "canonical_registry": {"status": "HOLD", "reason": "canonical registry reference is not supplied"}, "metadata": {**metadata, "input_sha256": _sha256(source), "normalization": normalization, "coverage": coverage}}
    except (ValueError, OSError, csv.Error, KeyError) as exc:
        payload = {"status": "HOLD", "hold_reason": str(exc), "cohort_id": cohort_id, "metadata": metadata}
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return payload


def run_declared_cohorts(*, config_path: str | Path, data_dir: str | Path, output_dir: str | Path, input_representation: str, analysis_unit: str | None = None) -> dict[str, Any]:
    """Run every declared transfer cohort; each missing local input is independently held."""
    config_file = Path(config_path)
    freeze = load_freeze(config_file.resolve().parents[1])
    outputs: dict[str, Any] = {}
    for cohort_id, contract in freeze["method"].get("cohorts", {}).items():
        required = list(contract.get("required_files", []))
        if not required:
            outputs[cohort_id] = {"status": "HOLD", "hold_reason": "cohort contract has no required input file"}
            continue
        unit = analysis_unit or contract.get("analysis_unit")
        outputs[cohort_id] = run_frozen_state_transfer(
            config_path=config_file,
            input_path=Path(data_dir) / required[0],
            output_path=Path(output_dir) / f"{cohort_id}.json",
            input_representation=input_representation,
            analysis_unit=unit,
            cohort_id=cohort_id,
        )
    status = "COMPUTED" if outputs and all(payload.get("status") == "COMPUTED" for payload in outputs.values()) else "HOLD"
    return {"producer": "A_FROZEN_STATE_TRANSFER", "status": status, "cohorts": outputs, "data_dir": _logical_path(Path(data_dir))}
