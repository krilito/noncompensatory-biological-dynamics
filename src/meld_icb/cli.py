"""Command-line entry point for the public reproduction preflight."""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import load_freeze


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", default=Path(__file__).resolve().parents[2])
    args = parser.parse_args(argv)
    freeze = load_freeze(args.repo_root)
    boundary = freeze["boundary"]
    print(f"loaded {boundary['schema_version']} freeze")
    print(f"threshold={boundary['threshold']}")
    print("EXTERNAL_ACQUISITION_REQUIRED: see docs/DATA_ACCESS.md and data/README.md")
    return 0
