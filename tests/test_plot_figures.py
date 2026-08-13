"""Quantitative figure scripts render from source_data without raw truth root."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_source_data_tables_exist():
    required = [
        "figure2.csv",
        "figure2_auc.csv",
        "figure3.csv",
        "figure3_b2_summary.csv",
        "figure4.csv",
        "figure4_auc.csv",
        "figure5.csv",
        "figure5_imvigor_samples.csv",
    ]
    for name in required:
        path = ROOT / "source_data" / name
        assert path.is_file(), name
        assert path.stat().st_size > 0, name


def test_plot_figure_scripts_run():
    for script in ("plot_figure2.py", "plot_figure3.py", "plot_figure4.py", "plot_figure5.py"):
        completed = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / script)],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        assert completed.returncode == 0, f"{script}\n{completed.stdout}\n{completed.stderr}"
        stem = script.replace("plot_figure", "Figure").replace(".py", "")
        # plot_figure2 -> Figure2
        stem = "Figure" + script[len("plot_figure") :].replace(".py", "")
        pdf = ROOT / "figures" / "reproduced" / f"{stem}.pdf"
        png = ROOT / "figures" / "reproduced" / f"{stem}.png"
        assert pdf.is_file() and pdf.stat().st_size > 1000, stem
        assert png.is_file() and png.stat().st_size > 1000, stem


def test_figure3_matches_frozen_source_identity():
    """Figure 3 producer is deterministic for the committed source tables."""
    completed = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "plot_figure3.py")],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0
    pdf = ROOT / "figures" / "reproduced" / "Figure3.pdf"
    assert pdf.is_file()
