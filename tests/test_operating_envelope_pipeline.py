import hashlib
import json
import gzip
from pathlib import Path

import pytest

from meld_icb.operating_envelope import classify_auc_interval
from qualified_dynamics.operating_envelope.contracts import ContractError, load_contract, validate_contract
from qualified_dynamics.operating_envelope.metrics import MetricHold, aggregate_patient_scores, select_earliest_biopsy
from qualified_dynamics.operating_envelope.pipeline import _transform_expression, run_declared_cohorts, run_operating_envelope
from qualified_dynamics.operating_envelope.verdicts import classify_verdict


def _record(sample, patient, score, response, order="1"):
    return {"sample_id": sample, "patient_id": patient, "boundary_score": score, "response_binary": response, "biopsy_order": order}


def _contract(tmp_path, *, input_name="rows.csv", declared_status="SUPPORTED"):
    authority = tmp_path / "authority.py"
    authority.write_text("pipeline authority\n", encoding="utf-8")
    digest = hashlib.sha256(authority.read_bytes()).hexdigest()
    return {
        "dataset_id": "TEST",
        "dataset_contract": {"accession": "TEST", "disease": "melanoma", "therapy_context": "ICI", "assay": "bulk", "timepoint": "pre", "observation_substrate": "biopsy", "local_input_contract": "rows"},
        "local_config_key": "datasets.TEST.root",
        "input_files": [input_name],
        "input_adapter": "tabular_rows",
        "inclusion_rules": {"biopsy_selection": "all_eligible"},
        "exclusion_rules": {"mapping_missing_or_ambiguous": "HOLD"},
        "analysis_unit": "patient",
        "response": {"column": "response", "positive_labels": ["R"], "negative_labels": ["NR"], "one_response_per_patient": True, "conflicting_labels": "HOLD"},
        "scorer": {"authority": "meld_icb.frozen_state_transfer", "preprocessing": "PRE_MEAN_THEN_AXIS_Z", "ddof": 1, "representation": "canonical_axis_z", "expression_transform": "truncate_negative_then_log2p1", "orientation": "higher_more_responder_like"},
        "statistics": {"auc": "rank_with_ties_half_credit", "ci": "registered_stratified_percentile_bootstrap", "bootstrap_resamples": 10000, "permutation_resamples": 10000, "seed": 20260731, "rng_algorithm": "python_random_shared_v1", "stream_policy": "shared_bootstrap_then_permutation_persistent_buffer", "inference_contract": "D_OPERATING_ENVELOPE_REGISTERED_V1", "permutation_comparator": "inclusive_absolute_deviation_ge", "correction": "plus_one", "minimum_class_count": 1},
        "carrier_qualification": {"status": "NOT_APPLICABLE", "required": False},
        "verdict": {"declared_status": declared_status, "primary_rule": "patient_mean_decides_verdict", "inverted_rule": "CI upper bound below 0.5", "posthoc_inversion": False},
        "provenance": {"lineage_source_path": "authority.py", "lineage_source_sha256": digest, "source_role": "test parser reference", "reference_objects_prohibited": True},
        "expected_counts": {},
    }


def test_patient_mean_is_unweighted_and_conflicting_labels_hold():
    records = [_record("a", "p1", 1.0, 1), _record("b", "p1", 3.0, 1), _record("c", "p2", 0.0, 0)]
    result = aggregate_patient_scores(records)
    assert result[0]["boundary_score"] == 2.0
    with pytest.raises(MetricHold):
        aggregate_patient_scores([_record("a", "p1", 1.0, 1), _record("b", "p1", 3.0, 0)])


def test_earliest_sensitivity_uses_declared_order_only():
    records = [_record("late", "p1", 100.0, 1, "2"), _record("early", "p1", -100.0, 1, "1")]
    assert select_earliest_biopsy(records)[0]["sample_id"] == "early"
    with pytest.raises(MetricHold):
        select_earliest_biopsy([_record("a", "p1", 1.0, 1, "1"), _record("b", "p1", 2.0, 1, "1")])
    with pytest.raises(MetricHold):
        select_earliest_biopsy([_record("a", "p1", 1.0, 1, ""), _record("b", "p1", 2.0, 1, "")])


def test_expression_transform_applies_log2_once():
    contract = {"scorer": {"expression_transform": "truncate_negative_then_log2p1"}}
    transformed = _transform_expression({"sample_id": "s", "_expression": {"GENE": 3.0}}, contract)
    assert transformed["GENE"] == pytest.approx(2.0)


def test_verdict_rules_do_not_infer_or_invert():
    assert classify_auc_interval(0.6, 0.51, 0.8, declared_status="AUTO") == "SUPPORTED"
    assert classify_auc_interval(0.4, 0.2, 0.49, declared_status="AUTO") == "INVERTED"
    assert classify_auc_interval(0.4, 0.2, 0.8, declared_status="AUTO") == "NEAR_NULL"
    assert classify_verdict({"auc": 0.4, "ci_95_low": 0.2, "ci_95_high": 0.8}, {"declared_status": "AUTO"}) == "NEAR_NULL"


def test_contract_rejects_posthoc_and_outcome_selection(tmp_path):
    contract = _contract(tmp_path)
    contract["verdict"]["posthoc_inversion"] = True
    with pytest.raises(ContractError):
        validate_contract(contract)
    contract = _contract(tmp_path)
    contract["inclusion_rules"]["biopsy_selection"] = "highest_score"
    with pytest.raises(ContractError):
        validate_contract(contract)


def test_pipeline_missing_data_and_unverified_source_hold(tmp_path):
    contract = _contract(tmp_path)
    config = tmp_path / "contract.yaml"
    config.write_text(json.dumps(contract), encoding="utf-8")
    result = run_operating_envelope(config_path=config, input_path=tmp_path / "missing.csv", output_path=tmp_path / "out.json", repo_root=Path(__file__).parents[1])
    assert result["status"] == "HOLD"
    assert "missing" in result["hold_reason"]


def test_fresh_clone_runs_without_internal_lineage_root(tmp_path):
    contract = _contract(tmp_path, input_name="rows.csv")
    contract["verdict"]["declared_status"] = "AUTO"
    config = tmp_path / "contract.yaml"
    config.write_text(json.dumps(contract), encoding="utf-8")
    source = tmp_path / "rows.csv"
    source.write_text("sample_id,patient_id,response,T,E,X,C,normalization_scope_id,normalization_source\na,p1,R,1,1,1,1,s,src\nb,p2,NR,0,0,0,0,s,src\n", encoding="utf-8")
    result = run_operating_envelope(config_path=config, input_path=source, output_path=tmp_path / "out.json", repo_root=Path(__file__).parents[1], representation="canonical_axis_z")
    assert result["status"] == "COMPUTED"
    assert result["metadata"]["lineage_source"]["status"] == "UNVERIFIED_LINEAGE"
    assert result["metadata"]["public_authority"]["status"] == "PASS"


def test_gse78220_auxiliary_resolves_from_data_root_without_lineage(tmp_path):
    """Public GSE78220 execution uses data_root for relative metadata inputs.

    ``source_root`` is reserved for optional lineage verification, so a public
    run with ``source_root=None`` must still find the family SOFT file beside
    the workbook under the declared data root and reach COMPUTED status.
    """
    openpyxl = pytest.importorskip("openpyxl")
    root = Path(__file__).parents[1]
    data_root = tmp_path / "data"
    workbook_path = data_root / "raw_data" / "GSE78220" / "GSE78220_PatientFPKM.xlsx"
    metadata_path = data_root / "raw_data" / "GSE78220_family.soft.gz"
    workbook_path.parent.mkdir(parents=True)

    sample_ids = [f"sample{index}" for index in range(1, 28)]
    genes = [
        "MKI67", "PCNA", "TOP2A", "CCNB1", "CDK1", "MCM2", "MCM5", "E2F1", "AURKB",
        "CD8A", "CD8B", "GZMB", "GZMA", "PRF1", "IFNG", "NKG7", "GNLY", "CXCL9", "CXCL10",
        "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "COL1A1", "COL1A2", "FN1", "VIM", "ZEB1", "SNAI1", "ACTA2", "CXCL12", "CD163", "MRC1", "CSF1R", "IL10", "ARG1", "IDO1", "VEGFA", "KDR", "FLT1",
        "PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "TOX2",
    ]
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = "FPKM"
    sheet.append(["Gene", *sample_ids])
    for gene_index, gene in enumerate(genes, start=1):
        # Distinct sample scaling keeps every declared axis variance non-zero.
        sheet.append([gene, *[float(gene_index * sample_index) for sample_index in range(1, 28)]])
    workbook.save(workbook_path)

    metadata_lines: list[str] = []
    for index, sample_id in enumerate(sample_ids, start=1):
        # Two biopsies for one aggregate patient satisfy the duplicate count.
        patient = f"Patient{index}" if index <= 25 else "PatientDuplicate"
        response = "Complete Response" if index <= 13 else "Progressive Disease"
        metadata_lines.extend(
            [
                f"^SAMPLE = GSM{index}",
                f"!Sample_title = {sample_id}",
                f"!Sample_characteristics_ch1 = patient id: {patient}",
                "!Sample_characteristics_ch1 = biopsy time: pre-treatment",
                f"!Sample_characteristics_ch1 = anti-pd-1 response: {response}",
            ]
        )
    with gzip.open(metadata_path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(metadata_lines) + "\n")

    result = run_declared_cohorts(
        config_dir=root / "configs" / "cohorts",
        data_root=data_root,
        output_dir=tmp_path / "cohorts",
        source_root=None,
        repo_root=root,
    )
    gse_result = result["cohorts"]["GSE78220"]
    assert gse_result["status"] == "COMPUTED", gse_result
    assert gse_result["population_summary"]["n_biopsies"] == 27
    assert gse_result["population_summary"]["n_patients"] == 26
    assert gse_result["population_summary"]["n_duplicate_patients"] == 1
    assert "duplicate_patients" not in gse_result["population_summary"]
    assert gse_result["reference_comparison"]["checks"]["n_duplicate_patients"] is True
    assert gse_result["metadata"]["lineage_source"]["status"] == "UNVERIFIED_LINEAGE"


def test_explicit_lineage_hash_mismatch_holds(tmp_path):
    contract = _contract(tmp_path, input_name="rows.csv")
    contract["provenance"]["lineage_source_sha256"] = "0" * 64
    config = tmp_path / "contract.yaml"
    config.write_text(json.dumps(contract), encoding="utf-8")
    source = tmp_path / "rows.csv"
    source.write_text("sample_id,patient_id,response,T,E,X,C,normalization_scope_id,normalization_source\na,p1,R,1,1,1,1,s,src\nb,p2,NR,0,0,0,0,s,src\n", encoding="utf-8")
    result = run_operating_envelope(config_path=config, input_path=source, output_path=tmp_path / "out.json", source_root=tmp_path, repo_root=Path(__file__).parents[1], representation="canonical_axis_z")
    assert result["status"] == "HOLD"
    assert "SHA256 mismatch" in result["hold_reason"]


def test_registry_input_is_not_a_producer(tmp_path):
    contract = _contract(tmp_path, input_name="final_numeric_registry.csv")
    config = tmp_path / "contract.yaml"
    config.write_text(json.dumps(contract), encoding="utf-8")
    source = tmp_path / "final_numeric_registry.csv"
    source.write_text("sample_id,patient_id,response,T,E,X,C,normalization_scope_id,normalization_source\na,p1,R,1,1,1,1,s,src\nb,p2,NR,0,0,0,0,s,src\n", encoding="utf-8")
    result = run_operating_envelope(config_path=config, input_path=source, output_path=tmp_path / "out.json", source_root=tmp_path, repo_root=Path(__file__).parents[1], representation="canonical_axis_z")
    assert result["status"] == "HOLD"
    assert "reference-only" in result["hold_reason"]
    assert result["metadata"]["authority"] == "qualified_dynamics.operating_envelope.pipeline"
    assert "lineage_source" in result["metadata"]
    assert "producer_source" not in result["metadata"]


def test_public_contracts_are_machine_readable_and_gse78220_primary_patient_rule():
    root = Path(__file__).parents[1]
    contracts = [load_contract(path) for path in (root / "configs" / "cohorts").glob("*.yaml")]
    assert {item["dataset_id"] for item in contracts} == {"GSE78220", "Sade-Feldman", "IMvigor210", "GSE282471"}
    gse = next(item for item in contracts if item["dataset_id"] == "GSE78220")
    assert gse["analysis_unit"] == "patient"
    assert gse["inclusion_rules"]["patient_aggregation"] == "unweighted_arithmetic_mean_of_all_eligible_biopsies"
    assert gse["verdict"]["declared_status"] == "AUTO"
    assert next(item for item in contracts if item["dataset_id"] == "GSE282471")["verdict"]["declared_status"] == "HOLD"
