"""Run the sole A scientific entry point with explicit local paths."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from meld_icb.frozen_state_transfer import run_declared_cohorts, run_frozen_state_transfer


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference-output", type=Path)
    parser.add_argument("--reference-input", type=Path)
    parser.add_argument("--representation", choices=("gene_expression", "canonical_axis_z"), required=False)
    parser.add_argument("--analysis-unit", choices=("sample", "patient"), required=False)
    parser.add_argument("--cohort-id", default="UNDECLARED")
    parser.add_argument("--cohort-dir", type=Path, help="run every cohort contract declared in the method config")
    args = parser.parse_args(argv)
    if args.cohort_dir is not None:
        if args.representation is None:
            parser.error("--representation is required with --cohort-dir")
        payload = run_declared_cohorts(config_path=args.config, data_dir=args.cohort_dir, output_dir=args.output.parent, input_representation=args.representation, analysis_unit=args.analysis_unit)
    else:
        payload = run_frozen_state_transfer(
            config_path=args.config,
            input_path=args.input,
            output_path=args.output,
            input_representation=args.representation,
            analysis_unit=args.analysis_unit,
            cohort_id=args.cohort_id,
            reference_output_path=args.reference_output,
            reference_input_path=args.reference_input,
        )
    print(payload["status"])
    return 0 if payload["status"] == "COMPUTED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
