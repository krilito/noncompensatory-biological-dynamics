"""Compatibility entrypoint for one object adapter; no producer execution."""
from __future__ import annotations
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
from meld_icb.object_adapters import run_object_adapter

def main(object_id: str) -> int:
    return run_object_adapter(object_id, shared_dir=ROOT / "outputs" / "shared_chains", output_dir=ROOT / "figures" / "reproduced", repo_root=ROOT)
