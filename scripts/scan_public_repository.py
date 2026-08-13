"""Fail closed on private, sensitive, or machine-bound public-repository content."""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKIP_DIRS = {
    ".git", "__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache",
    ".venv", "venv", "figures/reproduced", "outputs", "data/external",
    "data/permitted_derived",
}
FORBIDDEN_SUFFIXES = {
    ".ai", ".ait", ".psd", ".eps", ".pptx", ".key", ".tex", ".bib",
    ".docx", ".rtf", ".odt", ".h5ad", ".rds",
}
FORBIDDEN_COMPONENTS = {
    "private_assets", "governance", "agent_logs", "agent-traces", "prompts",
    "manuscript", "submission", "preprint", "reviewer", "recovery",
}
TEXT_SUFFIXES = {
    "", ".md", ".txt", ".py", ".json", ".yaml", ".yml", ".toml", ".tsv",
    ".csv", ".cff", ".ini", ".cfg",
}
SECRET_RE = re.compile(
    r"(?i)(?:api[_-]?key|access[_-]?token|password|secret|authorization\s*:\s*bearer|"
    r"ghp_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,})\s*[:=]\s*[^\s#]+"
)
PRIVATE_KEY_RE = re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----")
ABSOLUTE_PATH_RE = re.compile(r"(?:[A-Za-z]:\\|[A-Za-z]:/(?!/)|\\\\|/home/|/Users/|/mnt/)")
PRIVATE_IDENTITY_RE = re.compile(
    r"(?i)(?:qualified-dynamics-audit|private-workbench|D:\\Research|D:/Research|"
    r"governance/ai|agent_logs|agent-traces)"
)
REQUIRED = {
    "README.md", "LICENSE", "CITATION.cff", "THIRD_PARTY_NOTICES.md",
    "pyproject.toml", "environment.yml", "RELEASE_MANIFEST.sha256",
    "scripts/plot_figures.py", "scripts/scan_public_repository.py",
    "source_data/README.md", "docs/DATA_ACCESS.md", ".github/workflows/ci.yml",
}
CONTENT_EXEMPT = {"scripts/scan_public_repository.py"}


def _relative_files(root: Path) -> list[tuple[str, Path]]:
    files: list[tuple[str, Path]] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        if any(relative == item or relative.startswith(item + "/") for item in SKIP_DIRS):
            continue
        files.append((relative, path))
    return files


def scan_public_tree(root: Path = ROOT) -> list[str]:
    errors: list[str] = []
    files = _relative_files(root)
    actual = {relative for relative, _ in files}
    for missing in sorted(REQUIRED - actual):
        errors.append(f"MISSING_REQUIRED\t{missing}")
    for relative, path in files:
        parts = {part.lower() for part in Path(relative).parts[:-1]}
        if path.suffix.lower() in FORBIDDEN_SUFFIXES:
            errors.append(f"FORBIDDEN_EXTENSION\t{relative}")
        if parts & FORBIDDEN_COMPONENTS:
            errors.append(f"FORBIDDEN_COMPONENT\t{relative}")
        if path.suffix.lower() not in TEXT_SUFFIXES or path.stat().st_size > 5 * 1024 * 1024:
            continue
        if relative in CONTENT_EXEMPT:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            errors.append(f"TEXT_DECODE_ERROR\t{relative}")
            continue
        if SECRET_RE.search(text) or PRIVATE_KEY_RE.search(text):
            errors.append(f"SECRET_PATTERN\t{relative}")
        if ABSOLUTE_PATH_RE.search(text):
            errors.append(f"ABSOLUTE_PATH\t{relative}")
        if PRIVATE_IDENTITY_RE.search(text):
            errors.append(f"PRIVATE_IDENTITY\t{relative}")
    return errors


def main() -> int:
    errors = scan_public_tree(ROOT)
    if errors:
        print("\n".join(errors))
        return 1
    print("PUBLIC_REPOSITORY_SCAN_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
