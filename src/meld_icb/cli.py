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
    print(f"loaded {freeze['boundary']['representation']} freeze")
    print("EXTERNAL_ACQUISITION_REQUIRED: see data/external contracts")
    return 0
