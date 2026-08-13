from pathlib import Path

from scripts.scan_public_repository import scan_public_tree


def test_public_repository_boundary_passes() -> None:
    root = Path(__file__).resolve().parents[1]
    assert scan_public_tree(root) == []


def test_public_repository_boundary_rejects_private_and_machine_bound_content(tmp_path) -> None:
    required = {
        "README.md": "public",
        "LICENSE": "MIT License",
        "CITATION.cff": "cff-version: 1.2.0",
        "THIRD_PARTY_NOTICES.md": "notices",
        "pyproject.toml": "[project]",
        "environment.yml": "name: test",
        "RELEASE_MANIFEST.sha256": "manifest",
        "scripts/plot_figures.py": "pass",
        "scripts/scan_public_repository.py": "pass",
        "source_data/README.md": "data",
        "docs/DATA_ACCESS.md": "data",
        ".github/workflows/ci.yml": "name: test",
    }
    for relative, text in required.items():
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    private = tmp_path / "private_assets" / "Figure1.ai"
    private.parent.mkdir()
    private.write_bytes(b"private")
    bound = tmp_path / "notes.txt"
    bound.write_text("root = " + "C:" + chr(92) + "private", encoding="utf-8")
    errors = scan_public_tree(tmp_path)
    assert any(error.startswith("FORBIDDEN_EXTENSION") for error in errors)
    assert any(error.startswith("FORBIDDEN_COMPONENT") for error in errors)
    assert any(error.startswith("ABSOLUTE_PATH") for error in errors)
