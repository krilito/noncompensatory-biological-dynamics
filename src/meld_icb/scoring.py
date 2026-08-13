"""Frozen boundary score computed entirely from the checked-in configuration."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

from .config import load_freeze
from .frozen_state_transfer import score_axes

__all__ = ["score_axes", "score_sample", "classify"]


def score_sample(axes: Mapping[str, float], repo_root: str | Path) -> float:
    return score_axes(axes, load_freeze(repo_root)["boundary"])


def classify(score: float, repo_root: str | Path) -> int:
    boundary = load_freeze(repo_root)["boundary"]
    return int(float(score) >= float(boundary["threshold"]))
