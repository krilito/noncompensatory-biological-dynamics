"""Raw-input adapter for the frozen B paired-dynamics authority.

Only the public raw expression matrix, declared on-treatment labels, and the
hash-bound GSE91061 direction snapshot are accepted.  Historical B2/B3/B4
result tables and registries are intentionally not read here.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .config import load_freeze, load_json_yaml
from .frozen_state_inputs import _read_selected_matrix_fast
from .frozen_state_transfer import AXIS_ORDER, _axis_raw, collapse_duplicate_symbols


class PairedInputError(ValueError):
    """A raw B input cannot satisfy its declared mapping contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _require(path: Path) -> Path:
    if not path.exists():
        raise PairedInputError(f"required B input is missing: {path}")
    return path


def _rows(path: Path) -> list[dict[str, str]]:
    with _require(path).open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _binary(value: object) -> int:
    normalized = str(value).strip().upper()
    if normalized in {"1", "R", "RESPONDER"}:
        return 1
    if normalized in {"0", "NR", "NONRESPONDER"}:
        return 0
    raise PairedInputError(f"unsupported B response label: {value!r}")


def _paths(truth_root: Path, config_path: Path) -> dict[str, Path]:
    contract = load_json_yaml(config_path)
    values = contract.get("inputs", {})
    return {name: truth_root / str(relative) for name, relative in values.items() if isinstance(relative, str)}


def build_paired_records(
    truth_root: str | Path,
    *,
    config_path: str | Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build 16 PRE--EDT records from raw expression columns and labels."""
    root = Path(truth_root).resolve()
    config_file = Path(config_path).resolve()
    paths = _paths(root, config_file)
    labels = _rows(paths["labels"])
    expression_path = _require(paths["expression"])
    with expression_path.open("r", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    expression_columns = set(header)
    seen_keys: set[tuple[int, str]] = set()
    records: list[dict[str, Any]] = []
    for row in labels:
        try:
            patient_id = int(str(row["patient_id"]).strip())
        except (KeyError, TypeError, ValueError) as exc:
            raise PairedInputError(f"invalid B patient mapping: {row}") from exc
        therapy = str(row.get("treatment", "")).strip()
        sample_id = str(row.get("sample_id", "")).strip()
        if not therapy or not sample_id:
            raise PairedInputError("B label mapping has missing treatment or sample_id")
        expected_edt = f"{therapy}_{patient_id}_EDT"
        pre_sample = f"{therapy}_{patient_id}_PRE"
        if sample_id != expected_edt:
            raise PairedInputError(f"B label/sample conflict: {sample_id} != {expected_edt}")
        key = (patient_id, therapy)
        if key in seen_keys:
            raise PairedInputError(f"duplicate B therapy record: {patient_id}/{therapy}")
        seen_keys.add(key)
        if pre_sample not in expression_columns or sample_id not in expression_columns:
            # The two raw on-treatment labels without PRE columns are not
            # eligible paired records; this is a declared, outcome-independent
            # pairability gate rather than a registry fallback.
            continue
        records.append({
            "record_id": f"{patient_id}_{therapy}",
            "patient_id": patient_id,
            "therapy_type": therapy,
            "pre_sample_id": pre_sample,
            "edt_sample_id": sample_id,
            "y_true": _binary(row.get("response_binary")),
            "response_label": "Responder" if _binary(row.get("response_binary")) else "Nonresponder",
        })
    records.sort(key=lambda item: (int(item["patient_id"]), str(item["therapy_type"])))
    expected_n = int(load_json_yaml(config_file).get("contract", {}).get("expected_records", 16))
    expected_patients = int(load_json_yaml(config_file).get("contract", {}).get("expected_patients", 15))
    if len(records) != expected_n or len({item["patient_id"] for item in records}) != expected_patients:
        raise PairedInputError(f"B paired scope is {len(records)} records/{len({item['patient_id'] for item in records})} patients; expected {expected_n}/{expected_patients}")
    return records, {
        "expression_path": str(paths["expression"].relative_to(root).as_posix()),
        "labels_path": str(paths["labels"].relative_to(root).as_posix()),
        "expression_sha256": _sha256(expression_path),
        "labels_sha256": _sha256(paths["labels"]),
        "n_records": len(records),
        "n_patients": len({item["patient_id"] for item in records}),
        "response_counts": {
            "responder": sum(int(item["y_true"]) for item in records),
            "nonresponder": sum(not int(item["y_true"]) for item in records),
        },
    }


def load_paired_expression(
    truth_root: str | Path,
    records: list[Mapping[str, Any]],
    *,
    config_path: str | Path,
    freeze: Mapping[str, Any] | None = None,
) -> tuple[dict[str, dict[str, float]], dict[str, dict[str, float]], dict[str, Any]]:
    """Load selected raw columns and return raw axes plus full-scope z axes."""
    root = Path(truth_root).resolve()
    config_file = Path(config_path).resolve()
    paths = _paths(root, config_file)
    frozen = freeze or load_freeze(config_file.parents[1])
    sample_ids = [sample_id for record in records for sample_id in (str(record["pre_sample_id"]), str(record["edt_sample_id"]))]
    # Preserve first occurrence while retaining source registry order.
    ordered_ids = list(dict.fromkeys(sample_ids))
    expression_rows = _read_selected_matrix_fast(
        paths["expression"], ordered_ids, delimiter="\t", gene_field="Gene",
        transform="truncate_negative_then_log2p1",
        genes={str(gene).strip().upper() for genes in frozen["method"]["axes"]["gene_sets"].values() for gene in genes},
    )
    samples = {sample_id: row for sample_id, row in zip(ordered_ids, expression_rows)}
    collapsed = {sample_id: collapse_duplicate_symbols(row.items()) for sample_id, row in samples.items()}
    raw_rows, coverage = _axis_raw(list(collapsed.values()), frozen["method"]["axes"]["gene_sets"], float(frozen["method"]["preprocessing"]["minimum_gene_coverage"]))
    raw_axes = {sample_id: row for sample_id, row in zip(ordered_ids, raw_rows)}
    # Fit the single B2/B3 normalization over the declared PRE+EDT union.
    from .frozen_state_transfer import fit_axis_normalization, transform_axis_normalization
    parameters = fit_axis_normalization(raw_rows, ordered_ids, scope_id="paired_16records_PRE_EDT_union_n32", ddof=1)
    z_rows = transform_axis_normalization(raw_rows, parameters)
    z_axes = {sample_id: row for sample_id, row in zip(ordered_ids, z_rows)}
    return raw_axes, z_axes, {"coverage": coverage, "normalization_scope_id": parameters.normalization_scope_id, "axis_ddof": parameters.ddof}


def load_reference_direction(truth_root: str | Path, *, config_path: str | Path) -> tuple[dict[str, float], dict[str, Any]]:
    root = Path(truth_root).resolve()
    paths = _paths(root, Path(config_path).resolve())
    rows = _rows(paths["reference_direction"])
    if not rows or not {"axis", "GSE91061_delta"}.issubset(rows[0]):
        raise PairedInputError("B reference direction schema is missing axis/GSE91061_delta")
    values: dict[str, float] = {}
    for row in rows:
        axis = str(row["axis"]).strip()
        if axis in AXIS_ORDER:
            values[axis] = float(row["GSE91061_delta"])
    if tuple(values) != AXIS_ORDER:
        raise PairedInputError("B reference direction does not cover T/E/X/C")
    return values, {"reference_direction_path": str(paths["reference_direction"].relative_to(root).as_posix()), "reference_direction_sha256": _sha256(paths["reference_direction"])}
