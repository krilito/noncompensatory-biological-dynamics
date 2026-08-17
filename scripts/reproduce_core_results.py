"""LEVEL 2: reconcile committed Source Data with the manuscript numeric contract.

This is not analysis recomputation from raw transcriptomes. It reads the
committed tables and compares them with the rounded preprint values.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def load_contract(root: Path) -> dict[str, Any]:
    path = root / "manifests" / "manuscript_numeric_contract.yaml"
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != "MANUSCRIPT_NUMERIC_CONTRACT_V1":
        raise ValueError(f"unsupported numeric contract: {path}")
    return payload


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _select(rows: list[dict[str, str]], filters: dict[str, Any] | None) -> dict[str, str]:
    selected = rows
    for key, value in (filters or {}).items():
        selected = [row for row in selected if str(row.get(key, "")) == str(value)]
    if len(selected) != 1:
        raise ValueError(f"expected one row for {filters}, got {len(selected)}")
    return selected[0]


def reconcile(root: Path) -> list[dict[str, Any]]:
    contract = load_contract(root)
    rows: list[dict[str, Any]] = []
    for item in contract["results"]:
        table = _read_csv(root / item["source_path"])
        source_row = _select(table, item.get("filter"))
        recomputed = float(source_row[item["column"]])
        manuscript = float(item["manuscript"])
        abs_diff = abs(recomputed - manuscript)
        tolerance = float(item["tolerance"])
        rows.append(
            {
                "RESULT": item["id"],
                "MANUSCRIPT_VALUE": manuscript,
                "RECOMPUTED_VALUE": recomputed,
                "ABS_DIFF": abs_diff,
                "TOLERANCE": tolerance,
                "STATUS": "PASS" if abs_diff <= tolerance else "RELEASE_BLOCKER",
                "CLAIM_STATUS": item["status"],
            }
        )
    return rows


def render(rows: list[dict[str, Any]]) -> str:
    header = ["RESULT", "MANUSCRIPT_VALUE", "RECOMPUTED_VALUE", "ABS_DIFF", "TOLERANCE", "STATUS"]
    lines = ["\t".join(header)]
    for row in rows:
        lines.append(
            "\t".join(
                [
                    str(row["RESULT"]),
                    f"{row['MANUSCRIPT_VALUE']}",
                    f"{row['RECOMPUTED_VALUE']}",
                    f"{row['ABS_DIFF']}",
                    f"{row['TOLERANCE']}",
                    str(row["STATUS"]),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    rows = reconcile(args.repo_root)
    sys.stdout.write(render(rows))
    blockers = [row["RESULT"] for row in rows if row["STATUS"] != "PASS"]
    if blockers:
        print(f"RELEASE_BLOCKER: {', '.join(blockers)}", file=sys.stderr)
        return 1
    print("MANUSCRIPT_SOURCE_RECONCILIATION_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
