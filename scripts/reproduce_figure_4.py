"""Adapter status plus quantitative Figure 4 render from source_data."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from object_adapter_entrypoint import main as adapter_main

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    code = adapter_main("Figure4")
    plot = ROOT / "scripts" / "plot_figure4.py"
    if plot.exists():
        completed = subprocess.run([sys.executable, str(plot)], cwd=ROOT)
        if completed.returncode != 0:
            return completed.returncode
    return code


if __name__ == "__main__":
    raise SystemExit(main())
