"""Compatibility facade for object adapters.

``produce`` adapts an existing shared-chain result.  The historical raw-input
``produce_computed`` path is always retired and fail-closed, independent of
whether another shared output is present.
"""

from __future__ import annotations

import sys
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from meld_icb.object_adapters import run_object_adapter  # noqa: E402


def repository_root() -> Path:
    return ROOT


def produce(producer: str, claim_ids: list[str], required_cohorts: list[str], *, output_dir: str = "figures/reproduced") -> int:
    del claim_ids, required_cohorts
    return run_object_adapter(producer, shared_dir=ROOT / "outputs" / "shared_chains", output_dir=ROOT / output_dir, repo_root=ROOT)


def produce_computed(
    producer: str,
    claim_ids: list[str],
    required_cohorts: list[str],
    compute_kind: str,
    input_paths: list[str],
    *,
    output_dir: str = "figures/reproduced",
) -> int:
    """Fail closed for the removed raw-to-figure recomputation bypass."""
    del compute_kind, input_paths
    destination = ROOT / output_dir / f"{producer}.status.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": "OBJECT_ADAPTER_V1",
        "object_id": producer,
        "claim_ids": claim_ids,
        "required_cohorts": required_cohorts,
        "status": "HOLD_RAW_RECOMPUTE_RETIRED",
        "hold_reason": "historical raw-to-figure recomputation bypass is retired; invoke reproduce_all or an object adapter",
        "adapter_only": True,
    }
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"{producer}: {payload['status']}")
    return 2
