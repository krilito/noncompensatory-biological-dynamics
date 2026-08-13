"""Run the public C carrier-grounding producer."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from meld_icb.carrier_grounding import run_carrier_grounding  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth-root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    result = run_carrier_grounding(truth_root=args.truth_root, config_path=args.config)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "C_PUBLIC_NUMERIC_OUTPUT.json").write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    print(result.get("verdict", result.get("status", "HOLD")))
    return 0 if result.get("verdict") == "GATES_PASSED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
