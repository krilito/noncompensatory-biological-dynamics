"""The single D operating-envelope authority.

The implementation lives in :mod:`qualified_dynamics.operating_envelope.pipeline`.
The submodules are intentionally small contract, input, metric, and verdict
boundaries so that no second producer can grow in a compatibility namespace.
"""

from .contracts import ContractError, load_contract, validate_contract
from .metrics import aggregate_patient_scores, compute_metrics
from .pipeline import run_declared_cohorts, run_operating_envelope
from .verdicts import classify_verdict

__all__ = [
    "ContractError",
    "aggregate_patient_scores",
    "classify_verdict",
    "compute_metrics",
    "load_contract",
    "run_declared_cohorts",
    "run_operating_envelope",
    "validate_contract",
]
