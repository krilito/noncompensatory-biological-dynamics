"""Run analysis chains and/or quantitative figure redraws.

Modes
-----
analysis
    Execute A--D producers (when a truth root is available) and object
    adapters only. Never draws figures. Missing truth root exits 2 unless
    ``--allow-hold`` is set (still does not plot).

figures
    Redraw Figures 2--5 from committed ``source_data/`` tables only.
    Does not run producers. This is a quantitative redraw, not analysis
    recomputation.

all
    Run analysis, then figures only when explicitly requested with truth
    root present, or after ``--allow-hold`` status-only analysis.
    Prefer calling the modes separately for clarity.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from meld_icb.object_adapters import CHAIN_FILES, all_object_ids, run_object_adapter  # noqa: E402


PRODUCERS = [
    "reproduce_figure_1.py", "reproduce_figure_2.py", "reproduce_figure_3.py", "reproduce_figure_4.py", "reproduce_figure_5.py",
    *[f"reproduce_supplementary_s{i:02d}.py" for i in range(1, 13)],
    *[f"reproduce_table_st{i:02d}.py" for i in range(1, 9)],
]

PLOT_SCRIPTS = [
    ROOT / "scripts" / "plot_figure2.py",
    ROOT / "scripts" / "plot_figure3.py",
    ROOT / "scripts" / "plot_figure4.py",
    ROOT / "scripts" / "plot_figure5.py",
]


def _write(path: Path, payload: dict[str, object]) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return payload


def _hold(chain: str, reason: str, *, error_type: str | None = None) -> dict[str, object]:
    payload: dict[str, object] = {
        "producer": f"{chain}_SHARED_CHAIN",
        "status": "HOLD",
        "hold_reason": reason,
        "chain": chain,
    }
    if error_type is not None:
        payload["error_type"] = error_type
    return payload


def run_shared_chains(
    *,
    truth_root: Path | None = None,
    shared_dir: Path | None = None,
    strict: bool = False,
) -> tuple[dict[str, dict[str, object]], list[str]]:
    """Execute each scientific chain at most once and persist canonical output.

    Returns (outputs, hard_errors). hard_errors are producer exceptions (not missing data).
    """
    shared = shared_dir or ROOT / "outputs" / "shared_chains"
    shared.mkdir(parents=True, exist_ok=True)
    truth = truth_root or (Path(os.environ["MELD_TRUTH_ROOT"]) if os.environ.get("MELD_TRUTH_ROOT") else None)
    outputs: dict[str, dict[str, object]] = {}
    hard_errors: list[str] = []
    if truth is None or not truth.exists():
        for chain, filename in (
            ("A", "A_PUBLIC_NUMERIC_OUTPUT.json"),
            ("B", "B_PUBLIC_NUMERIC_OUTPUT.json"),
            ("C", "C_PUBLIC_NUMERIC_OUTPUT.json"),
            ("D", "D_PUBLIC_NUMERIC_OUTPUT.json"),
        ):
            outputs[chain] = _write(shared / filename, _hold(chain, "truth root is not declared or does not exist"))
        return outputs, hard_errors

    def _run_chain(chain: str, filename: str, fn, *, already_written: bool = False) -> None:
        try:
            result = fn()
            if already_written and isinstance(result, dict):
                outputs[chain] = result
                if not (shared / filename).exists():
                    _write(shared / filename, result)
            else:
                outputs[chain] = _write(shared / filename, result if isinstance(result, dict) else {"status": "COMPUTED", "payload": result})
        except Exception as exc:
            tb_path = shared / f"{chain}_EXCEPTION.txt"
            tb_path.write_text(traceback.format_exc(), encoding="utf-8")
            msg = f"{type(exc).__name__}: {exc} (traceback: {tb_path.as_posix()})"
            hard_errors.append(f"{chain}: {msg}")
            outputs[chain] = _write(
                shared / filename,
                _hold(chain, f"{chain} producer failed closed: {msg}", error_type=type(exc).__name__),
            )
            if strict:
                raise

    from meld_icb.frozen_state_inputs import run_truth_cohorts
    from meld_icb.carrier_grounding import run_carrier_grounding
    from scripts.run_paired_dynamics import run as run_b
    from qualified_dynamics.operating_envelope.pipeline import run_declared_cohorts

    _run_chain(
        "A",
        "A_PUBLIC_NUMERIC_OUTPUT.json",
        lambda: run_truth_cohorts(
            truth_root=truth,
            audit_dir=shared / "A_cohorts",
            config_path=ROOT / "configs" / "frozen_state_transfer.yaml",
            historical_comparison=None,
        ),
    )
    _run_chain(
        "B",
        "B_PUBLIC_NUMERIC_OUTPUT.json",
        lambda: run_b(truth_root=truth, config_path=ROOT / "configs" / "paired_dynamics.yaml", output_dir=shared),
        already_written=True,
    )
    _run_chain(
        "C",
        "C_PUBLIC_NUMERIC_OUTPUT.json",
        lambda: run_carrier_grounding(truth_root=truth, config_path=ROOT / "configs" / "carrier_grounding.yaml"),
    )
    _run_chain(
        "D",
        "D_PUBLIC_NUMERIC_OUTPUT.json",
        lambda: run_declared_cohorts(
            config_dir=ROOT / "configs" / "cohorts",
            data_root=truth,
            output_dir=shared / "D_cohorts",
            source_root=None,
            repo_root=ROOT,
        ),
    )

    return outputs, hard_errors


def run_adapters(*, shared: Path, output_dir: Path) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    object_ids = all_object_ids(ROOT)
    known_ids = set(object_ids)
    for stale in output_dir.glob("*.status.json"):
        if stale.stem not in known_ids:
            stale.unlink()
    hard_failures = 0
    nonfatal_holds = {"HOLD_COMPONENT_NOT_AVAILABLE", "HOLD_PARTIAL_COHORTS", "RETIRED", "EXCLUDED_BY_POLICY", "CONCEPTUAL_ASSET_CONTRACT"}
    for object_id in object_ids:
        return_code = run_object_adapter(object_id, shared_dir=shared, output_dir=output_dir, repo_root=ROOT)
        if return_code == 0:
            continue
        status_path = output_dir / f"{object_id}.status.json"
        try:
            status = json.loads(status_path.read_text(encoding="utf-8")).get("status")
        except (OSError, json.JSONDecodeError, AttributeError):
            status = None
        malformed_chain = False
        if status in nonfatal_holds:
            try:
                chain = json.loads(status_path.read_text(encoding="utf-8")).get("shared_chain")
                for name in CHAIN_FILES.get(str(chain), ()):
                    candidate = shared / name
                    if not candidate.exists():
                        continue
                    try:
                        json.loads(candidate.read_text(encoding="utf-8"))
                    except (OSError, json.JSONDecodeError):
                        malformed_chain = True
            except (OSError, json.JSONDecodeError, AttributeError):
                malformed_chain = True
        if status not in nonfatal_holds or malformed_chain:
            hard_failures += 1
    return hard_failures


def run_figures() -> int:
    """Quantitative redraw from committed source_data only."""
    failures = 0
    for plot in PLOT_SCRIPTS:
        if not plot.exists():
            failures += 1
            continue
        completed = subprocess.run([sys.executable, str(plot)], cwd=ROOT)
        if completed.returncode != 0:
            failures += 1
    return failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-root", type=Path)
    parser.add_argument("--fail-on-missing", action="store_true", help="Exit 2 if any adapter is not terminal COMPUTED/RETIRED/EXCLUDED/CONCEPTUAL")
    parser.add_argument(
        "--mode",
        choices=("analysis", "figures", "all"),
        default="analysis",
        help="analysis=recompute chains+adapters (no plots); figures=source_data redraw only; all=analysis then figures",
    )
    parser.add_argument(
        "--allow-hold",
        action="store_true",
        help="In analysis mode without truth root, write HOLD statuses and exit 0 if adapters are only expected holds",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Propagate producer exceptions as non-zero exit (still write HOLD payload before exit)",
    )
    args = parser.parse_args(argv)

    if args.mode == "figures":
        failures = run_figures()
        return 1 if failures else 0

    shared = ROOT / "outputs" / "shared_chains"
    output_dir = ROOT / "figures" / "reproduced"
    truth = args.truth_root
    if (truth is None or not Path(truth).exists()) and not args.allow_hold:
        # Explicit non-success: analysis cannot run without inputs.
        run_shared_chains(truth_root=truth, shared_dir=shared, strict=False)
        run_adapters(shared=shared, output_dir=output_dir)
        print("ANALYSIS_HOLD: truth root is not declared or does not exist; no figures drawn", file=sys.stderr)
        return 2

    try:
        _, hard_errors = run_shared_chains(truth_root=truth, shared_dir=shared, strict=args.strict)
    except Exception as exc:
        print(f"STRICT_PRODUCER_FAILURE: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1

    hard_failures = run_adapters(shared=shared, output_dir=output_dir)
    if hard_errors:
        hard_failures += len(hard_errors)
        for err in hard_errors:
            print(f"PRODUCER_ERROR: {err}", file=sys.stderr)

    if args.mode == "all":
        # Only after analysis path; figures remain a separate quantitative redraw.
        fig_failures = run_figures()
        if fig_failures:
            hard_failures += fig_failures

    if args.fail_on_missing:
        terminal = {
            "COMPUTED",
            "COMPUTED_PRIMARY_PLUS_COMPARISON_ONLY",
            "CONCEPTUAL_ASSET_CONTRACT",
            "EXCLUDED_BY_POLICY",
            "RETIRED",
        }
        object_ids = all_object_ids(ROOT)
        known_ids = set(object_ids)
        status_paths = {path.stem: path for path in output_dir.glob("*.status.json")}
        if set(status_paths) != known_ids or len(status_paths) != 25:
            return 2
        statuses = [json.loads(status_paths[object_id].read_text(encoding="utf-8")).get("status") for object_id in object_ids]
        if any(status not in terminal for status in statuses):
            return 2
    return 1 if hard_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
