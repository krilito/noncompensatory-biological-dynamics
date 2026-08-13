"""Ensure human-readable FROZEN_BOUNDARY.md cannot drift from machine YAML."""

from __future__ import annotations

import json
import os
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
MACHINE = REPO / "configs" / "frozen_state_transfer.yaml"


def _project_spec() -> Path | None:
    env = os.environ.get("MELD_PROJECT_ROOT")
    if not env:
        return None
    p = Path(env) / "algorithm-core" / "specification" / "FROZEN_BOUNDARY.md"
    return p if p.is_file() else None


def test_machine_config_coefficients_are_self_consistent() -> None:
    cfg = json.loads(MACHINE.read_text(encoding="utf-8"))
    score = cfg["score"]
    assert score["T"] == pytest.approx(-0.676)
    assert score["E"] == pytest.approx(0.173)
    assert score["X"] == pytest.approx(0.0)
    assert score["C"] == pytest.approx(1.135)
    assert score["intercept"] == pytest.approx(-0.666)
    assert cfg["threshold"] == pytest.approx(0.735)
    assert cfg["preprocessing"]["method"] == "PRE_MEAN_THEN_AXIS_Z"
    assert cfg["preprocessing"]["ddof"] == 1


def test_human_spec_matches_machine_when_present() -> None:
    spec = _project_spec()
    if spec is None:
        pytest.skip("MELD FROZEN_BOUNDARY.md not found; set MELD_PROJECT_ROOT to enforce sync")
    md = spec.read_text(encoding="utf-8")
    cfg = json.loads(MACHINE.read_text(encoding="utf-8"))
    score = cfg["score"]
    # require literal coefficient strings in the human spec
    for value in (
        str(score["T"]),
        str(score["E"]),
        str(score["C"]),
        str(score["intercept"]),
        str(cfg["threshold"]),
        "0.0",  # X
        "PRE_MEAN_THEN_AXIS_Z",
    ):
        assert value in md, f"missing {value!r} in human specification"
    assert re.search(r"ddof\s*=\s*1", md), "ddof=1 missing from human specification"
    assert "Not** a machine input" in md or "Not a machine input" in md or "must **not** parse" in md
