"""Analysis recomputation entrypoint (A--D producers + adapters; no figures)."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from reproduce_all import main as reproduce_main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(reproduce_main(["--mode", "analysis", *sys.argv[1:]]))
