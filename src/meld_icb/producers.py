"""Retired raw figure producers.

Figure/table objects are now emitted by :mod:`meld_icb.object_adapters` from
shared-chain aggregates.  These names remain only as a fail-closed
compatibility surface for callers of the historical API; they never import
or invoke scientific statistics.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def _retired(name: str) -> dict[str, Any]:
    return {"producer": name, "status": "HOLD_COMPONENT_NOT_AVAILABLE", "hold_reason": "raw-to-figure producer retired; use shared-chain object adapter"}


def figure1(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return _retired("Figure1")


def figure2(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return _retired("Figure2")


def figure3(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return _retired("Figure3")


def figure4(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return _retired("Figure4")


def figure5(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return _retired("Figure5")


def write_result(result: dict[str, Any], output_dir: str | Path) -> Path:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    producer = str(result.get("producer", "retired"))
    path = destination / f"{producer}.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return path
