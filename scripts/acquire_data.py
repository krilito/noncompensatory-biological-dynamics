"""Emit an acquisition checklist; never copy data without explicit terms."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", help="cohort_id from manifests/cohorts.tsv")
    parser.add_argument("--output", default="data/external/acquisition_checklist.tsv")
    args = parser.parse_args()
    rows = list(csv.DictReader((root / "manifests" / "cohorts.tsv").open(encoding="utf-8"), delimiter="\t"))
    if args.cohort:
        rows = [row for row in rows if row["cohort_id"] == args.cohort]
        if not rows:
            raise SystemExit(f"unknown cohort: {args.cohort}")
    out = root / args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["cohort_id", "accession", "source", "download_method", "expected_schema", "checksum", "terms", "status"]
    with out.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(fields) + "\n")
        for row in rows:
            handle.write("\t".join(row[field] if field in row else "PENDING" for field in fields) + "\n")
    print(f"wrote acquisition checklist: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
