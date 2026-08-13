"""Unit and contract checks for the sole meld_icb A entry point."""

import csv
import json
from pathlib import Path

import pytest

from meld_icb.config import load_freeze
from meld_icb.frozen_state_transfer import (
    aggregate_records,
    classify_boundary,
    collapse_duplicate_symbols,
    frozen_boundary_score,
    pre_mean_then_axis_z,
    pre_mean_then_axis_z_transfer,
    run_frozen_state_transfer,
    run_declared_cohorts,
)
from meld_icb.scoring import score_axes as legacy_score_axes

ROOT = Path(__file__).resolve().parents[1]
FREEZE = load_freeze(ROOT)
METHOD = FREEZE["method"]
AXES = METHOD["axes"]["order"]
GENE_SETS = METHOD["axes"]["gene_sets"]


def _sample(seed: float) -> dict[str, float]:
    return {gene: seed * (axis_index + 1) + gene_index * 0.01 for axis_index, axis in enumerate(AXES) for gene_index, gene in enumerate(GENE_SETS[axis])}


def test_single_method_authority_and_frozen_values():
    assert FREEZE["analysis"]["method_config"] == "configs/frozen_state_transfer.yaml"
    assert FREEZE["boundary"]["coefficients"] == {"T": -0.676, "E": 0.173, "X": 0.0, "C": 1.135}
    assert FREEZE["boundary"]["intercept"] == -0.666
    assert FREEZE["boundary"]["threshold"] == 0.735
    assert frozen_boundary_score({"T": 0, "E": 0, "X": 0, "C": 0}, FREEZE) == -0.666
    assert classify_boundary(0.735, FREEZE) == 1
    assert (ROOT / "src" / "qualified_dynamics" / "operating_envelope" / "pipeline.py").exists()
    assert not (ROOT / "configs" / "frozen_boundary.yaml").exists()
    assert not (ROOT / "configs" / "preprocessing.yaml").exists()
    assert legacy_score_axes({"T": 0, "E": 0, "X": 0, "C": 0}, FREEZE["boundary"]) == frozen_boundary_score({"T": 0, "E": 0, "X": 0, "C": 0}, FREEZE)


def test_duplicate_gene_and_coverage_policy():
    assert collapse_duplicate_symbols([("A", 1), (" a ", 3)]) == {"A": 2.0}
    samples = [_sample(1), _sample(2), _sample(3)]
    axes, parameters, coverage = pre_mean_then_axis_z(samples, FREEZE, ["s1", "s2", "s3"], scope_id="fit")
    assert parameters.ddof == 1
    assert coverage["T"]["observed_genes"] == len(GENE_SETS["T"])
    assert len(axes) == 3
    incomplete = [_sample(1), _sample(2), _sample(3)]
    for sample in incomplete:
        for gene in GENE_SETS["T"][:5]:
            sample.pop(gene)
    with pytest.raises(ValueError, match="minimum coverage"):
        pre_mean_then_axis_z(incomplete, FREEZE, ["s1", "s2", "s3"], scope_id="fit")


def test_reference_fit_is_reused_for_held_out_target():
    reference = [_sample(1), _sample(2), _sample(3)]
    target = [_sample(4), _sample(5)]
    axes, parameters, coverage = pre_mean_then_axis_z_transfer(reference, target, FREEZE, ["r1", "r2", "r3"], scope_id="reference")
    assert parameters.fitted_sample_ids == ("r1", "r2", "r3")
    assert parameters.normalization_scope_id == "reference"
    assert len(axes) == 2
    assert set(coverage) == {"reference", "target"}


def test_explicit_aggregation_unit():
    records = [
        {"sample_id": "s1", "patient_id": "p1", "boundary_score": 1.0, "response_binary": 1},
        {"sample_id": "s2", "patient_id": "p1", "boundary_score": 3.0, "response_binary": 1},
        {"sample_id": "s3", "patient_id": "p2", "boundary_score": -1.0, "response_binary": 0},
    ]
    patient = aggregate_records(records, "patient")
    assert len(patient) == 2
    assert patient[0]["boundary_score"] == 2.0
    with pytest.raises(ValueError):
        aggregate_records(records, "unknown")


def _write_input(path: Path, representation: str = "gene_expression") -> None:
    fields = ["sample_id", "patient_id", "response"]
    if representation == "gene_expression":
        fields.extend(gene for axis in AXES for gene in GENE_SETS[axis])
    else:
        fields.extend([*AXES, "normalization_scope_id", "normalization_source"])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index, response in enumerate((0, 0, 1, 1), start=1):
            row = {"sample_id": f"s{index}", "patient_id": f"p{index}", "response": response}
            if representation == "gene_expression":
                row.update(_sample(float(index)))
            else:
                row.update({axis: index * (axis_index + 1) for axis_index, axis in enumerate(AXES)})
                row.update({"normalization_scope_id": "scope", "normalization_source": "local_fixture"})
            writer.writerow(row)


def test_axis_representation_must_be_explicit_and_missing_reference_holds(tmp_path):
    source = tmp_path / "input.tsv"
    _write_input(source)
    output = tmp_path / "output.json"
    payload = run_frozen_state_transfer(config_path=ROOT / "configs" / "frozen_state_transfer.yaml", input_path=source, output_path=output, input_representation=None, analysis_unit="sample", cohort_id="PRJEB23709")
    assert payload["status"] == "HOLD"
    assert "explicit" in payload["hold_reason"]
    transfer = run_frozen_state_transfer(config_path=ROOT / "configs" / "frozen_state_transfer.yaml", input_path=source, output_path=output, input_representation="gene_expression", analysis_unit="sample", cohort_id="PRJEB23709", reference_input_path=tmp_path / "absent_reference.tsv")
    assert transfer["status"] == "HOLD"
    assert "reference" in transfer["hold_reason"]


def test_gene_expression_entrypoint_computes_registered_metrics_and_holds_stability(tmp_path):
    source = tmp_path / "input.tsv"
    _write_input(source)
    payload = run_frozen_state_transfer(config_path=ROOT / "configs" / "frozen_state_transfer.yaml", input_path=source, output_path=tmp_path / "out.json", input_representation="gene_expression", analysis_unit="sample", cohort_id="PRJEB23709")
    assert payload["status"] == "COMPUTED"
    assert payload["cohort_metrics"]["analysis_unit"] == "sample"
    assert payload["cohort_metrics"]["inference_status"] == "COMPUTED_REGISTERED"
    assert payload["stability"]["status"] == "HOLD"


def test_canonical_axis_input_requires_provenance(tmp_path):
    source = tmp_path / "axes.tsv"
    _write_input(source, representation="canonical_axis_z")
    rows = list(csv.DictReader(source.open(encoding="utf-8"), delimiter="\t"))
    for row in rows:
        row.pop("normalization_source", None)
    with source.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    payload = run_frozen_state_transfer(config_path=ROOT / "configs" / "frozen_state_transfer.yaml", input_path=source, output_path=tmp_path / "out.json", input_representation="canonical_axis_z", analysis_unit="sample")
    assert payload["status"] == "HOLD"
    assert "provenance" in payload["hold_reason"]


def test_stability_is_hold_not_omitted_auc_and_missing_input_is_hold(tmp_path):
    output = tmp_path / "missing.json"
    payload = run_frozen_state_transfer(config_path=ROOT / "configs" / "frozen_state_transfer.yaml", input_path=tmp_path / "absent.tsv", output_path=output, input_representation="gene_expression", analysis_unit="sample")
    assert payload["status"] == "HOLD"
    assert "missing" in payload["hold_reason"]
    assert "boundary_scores" not in json.loads(output.read_text(encoding="utf-8"))


def test_declared_external_cohorts_hold_independently_without_local_inputs(tmp_path):
    payload = run_declared_cohorts(config_path=ROOT / "configs" / "frozen_state_transfer.yaml", data_dir=tmp_path / "external", output_dir=tmp_path / "cohort_outputs", input_representation="gene_expression")
    assert payload["status"] == "HOLD"
    assert set(payload["cohorts"]) == {"PRJEB23709", "MORRISON-1-public", "MGH"}
    assert all(item["status"] == "HOLD" for item in payload["cohorts"].values())
