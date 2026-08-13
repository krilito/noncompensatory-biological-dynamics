"""Run the frozen IFNG6 paired information decomposition audit."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from meld_icb.ifng6_information_decomposition import (
    BOOTSTRAP_REPEATS,
    PRIMARY_SEED,
    run_ifng6_information_decomposition,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prior-audit-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--bootstrap-repeats", type=int, default=BOOTSTRAP_REPEATS)
    parser.add_argument("--seed", type=int, default=PRIMARY_SEED)
    args = parser.parse_args(argv)
    result = run_ifng6_information_decomposition(
        prior_audit_dir=args.prior_audit_dir,
        repository_root=ROOT,
        output_dir=args.output_dir,
        bootstrap_repeats=args.bootstrap_repeats,
        seed=args.seed,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
