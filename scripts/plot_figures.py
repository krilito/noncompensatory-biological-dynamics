"""Quantitative redraw entrypoint for Figures 2-5 from committed source_data.

This is not analysis recomputation. Editorial Adobe layouts are out of scope.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from reproduce_all import main as reproduce_main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(reproduce_main(["--mode", "figures", *sys.argv[1:]]))
