"""Cross-sectional transfer estimands; data acquisition remains external."""

from __future__ import annotations

from typing import Sequence

from .statistics import auc


def response_auc(scores: Sequence[float], responder_labels: Sequence[int]) -> float:
    """Return the response-associated AUC without label or score inversion."""
    return auc(scores, responder_labels)
