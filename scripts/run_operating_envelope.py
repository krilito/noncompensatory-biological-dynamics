"""Run one or all declared D operating-envelope contracts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from qualified_dynamics.operating_envelope.pipeline import run_declared_cohorts, run_operating_envelope


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--data-root", type=Path, default=Path("data/external"))
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config-dir", type=Path, default=Path("configs/cohorts"))
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--representation", choices=("gene_expression", "canonical_axis_z"))
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    if args.all:
        payload = run_declared_cohorts(
            config_dir=repo_root / args.config_dir,
            data_root=args.data_root,
            output_dir=args.output,
            source_root=args.source_root,
            repo_root=repo_root,
            representation=args.representation,
        )
    else:
        if args.config is None:
            parser.error("--config is required unless --all is set")
        payload = run_operating_envelope(
            config_path=repo_root / args.config,
            input_path=args.input,
            output_path=args.output,
            source_root=args.source_root,
            repo_root=repo_root,
            data_root=args.data_root,
            representation=args.representation,
        )
    print(json.dumps(payload, indent=2, sort_keys=True, default=str))
    return 0 if payload.get("status") == "COMPUTED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
