"""Declared verdict rules with no post-hoc score inversion."""

from __future__ import annotations

from typing import Any, Mapping

VALID_STATUSES = {"AUTO", "SUPPORTED", "NEAR_NULL", "INVERTED", "UNRESOLVED", "HOLD"}


def classify_verdict(metrics: Mapping[str, Any], verdict: Mapping[str, Any]) -> str:
    declared = str(verdict.get("declared_status", "")).upper()
    if declared not in VALID_STATUSES:
        raise ValueError(f"unknown declared status: {declared}")
    if bool(verdict.get("posthoc_inversion", False)):
        raise ValueError("post hoc score inversion is prohibited")
    if declared == "HOLD":
        return declared
    auc_value = metrics.get("auc")
    ci_low = metrics.get("ci_95_low")
    ci_high = metrics.get("ci_95_high")
    if auc_value is None or ci_low is None or ci_high is None:
        if declared == "UNRESOLVED":
            return declared
        raise ValueError("computed verdict requires AUC and confidence interval")
    if not (0.0 <= float(ci_low) <= float(ci_high) <= 1.0):
        raise ValueError("AUC interval must lie in [0, 1]")
    if float(ci_low) > 0.5:
        return "SUPPORTED"
    if float(ci_high) < 0.5:
        return "INVERTED"
    return "NEAR_NULL"


def classify_auc_interval(
    auc_value: float | None,
    ci_low: float | None,
    ci_high: float | None,
    *,
    declared_status: str,
    hold: bool = False,
) -> str:
    """Compatibility signature for existing public tests/imports."""
    normalized = str(declared_status).upper()
    if hold:
        if normalized not in {"HOLD", "UNRESOLVED"}:
            raise ValueError("HOLD input cannot be relabelled")
        return normalized
    if auc_value is None or ci_low is None or ci_high is None:
        if normalized in {"HOLD", "UNRESOLVED"}:
            return normalized
        raise ValueError("computed verdict requires AUC and confidence interval")
    return classify_verdict(
        {"auc": auc_value, "ci_95_low": ci_low, "ci_95_high": ci_high},
        {"declared_status": declared_status, "posthoc_inversion": False},
    )
