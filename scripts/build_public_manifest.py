"""Build the deterministic public file hash manifest."""

from __future__ import annotations

import hashlib
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    manifest = root / "manifests" / "public_files.sha256"
    excluded = {manifest.resolve()}
    excluded_dirs = {".git", "__pycache__", ".pytest_cache"}
    excluded_output = root / "figures" / "reproduced"
    rows: list[tuple[str, str]] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.resolve() in excluded:
            continue
        if any(part in excluded_dirs for part in path.parts):
            continue
        if excluded_output in path.parents:
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append((digest, path.relative_to(root).as_posix()))
    manifest.write_text("\n".join(f"{digest}  {relative}" for digest, relative in rows) + "\n", encoding="utf-8")
    print(f"wrote {len(rows)} public file hashes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
