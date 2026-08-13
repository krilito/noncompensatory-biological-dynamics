"""Run the audit-only falsification-first state/movement analysis."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from meld_icb.falsification_first import run_falsification_first


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--bootstrap-repeats", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260807)
    args = parser.parse_args(argv)
    result = run_falsification_first(
        truth_root=args.truth_root,
        repository_root=ROOT,
        output_dir=args.output_dir,
        bootstrap_repeats=args.bootstrap_repeats,
        seed=args.seed,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
