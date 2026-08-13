"""Compatibility façade for the single qualified D operating-envelope authority."""

from qualified_dynamics.operating_envelope.contracts import ContractError, load_contract, validate_contract
from qualified_dynamics.operating_envelope.metrics import aggregate_patient_scores, compute_metrics, select_earliest_biopsy
from qualified_dynamics.operating_envelope.pipeline import run_declared_cohorts, run_operating_envelope
from qualified_dynamics.operating_envelope.verdicts import classify_auc_interval, classify_verdict

__all__ = [
    "ContractError",
    "aggregate_patient_scores",
    "classify_auc_interval",
    "classify_verdict",
    "compute_metrics",
    "load_contract",
    "run_declared_cohorts",
    "run_operating_envelope",
    "select_earliest_biopsy",
    "validate_contract",
]
