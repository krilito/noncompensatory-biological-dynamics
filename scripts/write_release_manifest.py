"""Write RELEASE_MANIFEST.sha256 from the Git candidate file set.

Uses `git ls-files --cached --others --exclude-standard`, so ignored
build metadata such as `*.egg-info/` is never hashed. Files deleted in
the working tree but still in the index are omitted.
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_NAME = "RELEASE_MANIFEST.sha256"


def candidate_files(root: Path = ROOT) -> list[str]:
    completed = subprocess.run(
        [
            "git",
            "--no-optional-locks",
            "-C",
            str(root),
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
            "-z",
        ],
        check=True,
        capture_output=True,
    )
    relatives: list[str] = []
    for raw in completed.stdout.split(b"\0"):
        if not raw:
            continue
        relative = raw.decode("utf-8")
        if relative == MANIFEST_NAME:
            continue
        path = root / relative
        if path.is_file():
            relatives.append(Path(relative).as_posix())
    return sorted(set(relatives))


def write_manifest(root: Path = ROOT) -> Path:
    rows: list[str] = []
    for relative in candidate_files(root):
        digest = hashlib.sha256((root / relative).read_bytes()).hexdigest()
        rows.append(f"{digest}  {relative}")
    destination = root / MANIFEST_NAME
    destination.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return destination


def main() -> int:
    path = write_manifest(ROOT)
    count = sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line.strip())
    print(f"wrote {count} hashes to {path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
