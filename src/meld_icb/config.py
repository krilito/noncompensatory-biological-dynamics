"""Load the single machine-readable A method authority."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_json_yaml(path: str | Path) -> dict[str, Any]:
    """Load the repository's JSON-compatible YAML files with stdlib only."""
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"configuration must be an object: {path}")
    return value


def load_freeze(repo_root: str | Path) -> dict[str, Any]:
    root = Path(repo_root)
    analysis = load_json_yaml(root / "configs" / "analysis_freeze.yaml") if (root / "configs" / "analysis_freeze.yaml").exists() else {}
    method_path = root / analysis.get("method_config", "configs/frozen_state_transfer.yaml")
    method = load_json_yaml(method_path)
    preprocessing = dict(method["preprocessing"])
    preprocessing["method"] = preprocessing.pop("method")
    preprocessing["minimum_axis_coverage"] = preprocessing["minimum_gene_coverage"]
    boundary = {
        "schema_version": method["schema_version"],
        "axes": list(method["axes"]["order"]),
        "coefficients": {axis: float(method["score"][axis]) for axis in method["axes"]["order"]},
        "intercept": float(method["score"]["intercept"]),
        "threshold": float(method["threshold"]),
        "orientation": method["score"]["orientation"],
        "preprocessing": preprocessing["method"],
        "ddof": int(preprocessing["ddof"]),
        "missing_gene_rule": preprocessing["missing_gene_rule"],
    }
    return {
        "analysis": {**analysis, "method_config": method_path.relative_to(root).as_posix()},
        "method": method,
        "boundary": boundary,
        "preprocessing": preprocessing,
    }
