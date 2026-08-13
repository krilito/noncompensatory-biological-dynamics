"""Strict tabular schemas used by public producers."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable


def read_tsv(path: str | Path, required_columns: Iterable[str]) -> list[dict[str, str]]:
    file_path = Path(path)
    with file_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"missing header: {file_path}")
        missing = [column for column in required_columns if column not in reader.fieldnames]
        if missing:
            raise ValueError(f"missing columns {missing}: {file_path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"empty data table: {file_path}")
    return rows
