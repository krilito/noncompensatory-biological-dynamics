"""Run the public B paired-movement/specificity/incremental producers."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from meld_icb.config import load_freeze  # noqa: E402
from meld_icb.incremental_value import fold_safe_incremental_metrics, influence_leave_one_patient_reruns  # noqa: E402
from meld_icb.paired_movement import paired_movement_metrics  # noqa: E402
from meld_icb.paired_movement_inputs import (  # noqa: E402
    PairedInputError,
    build_paired_records,
    load_paired_expression,
    load_reference_direction,
)
from meld_icb.response_specificity import response_specificity_metrics  # noqa: E402


def run(*, truth_root: str | Path, config_path: str | Path, output_dir: str | Path) -> dict[str, object]:
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    try:
        config_file = Path(config_path).resolve()
        freeze = load_freeze(config_file.parents[1])
        records, provenance = build_paired_records(truth_root, config_path=config_file)
        raw_axes, z_axes, normalization = load_paired_expression(truth_root, records, config_path=config_file, freeze=freeze)
        reference, reference_provenance = load_reference_direction(truth_root, config_path=config_file)
        b2 = paired_movement_metrics(records, z_axes, reference, freeze)
        movement_rows = b2.pop("_rows")
        b3 = response_specificity_metrics(records, movement_rows, seed=42)
        b4 = fold_safe_incremental_metrics(records, raw_axes, freeze, seed=42)
        influence = influence_leave_one_patient_reruns(records, raw_axes, freeze, seed=42)
        result = {
            "producer": "B_PUBLIC_PAIRED_DYNAMICS",
            "status": "COMPUTED",
            "b2_movement": b2,
            "b3_response_specificity": b3,
            "b4_incremental_foldsafe": {key: value for key, value in b4.items() if key != "patient_delta_L"},
            "b4_influence_leave_one_patient": {
                key: influence[key]
                for key in (
                    "analysis_id",
                    "rebuilds_folds",
                    "refits_normalization",
                    "refits_models",
                    "n_reruns",
                    "all_means_negative",
                    "mean_delta_L_min",
                    "mean_delta_L_max",
                )
            },
            "provenance": {**provenance, **reference_provenance, **normalization, "truth_root_not_serialized": True},
            "comparison_only": {"status": "EXCLUDED_FROM_PRODUCER", "objects": ["B2_CANONICAL_MOVEMENT_RESULT.json", "B3_CANONICAL_SPECIFICITY_RESULT.json", "B4_FOLDSAFE_RESULT.json", "PAIRED_RECORD_REGISTRY.csv"]},
        }
    except (PairedInputError, OSError, ValueError, KeyError) as exc:
        result = {"producer": "B_PUBLIC_PAIRED_DYNAMICS", "status": "HOLD", "hold_reason": str(exc)}
    (output / "B_PUBLIC_NUMERIC_OUTPUT.json").write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth-root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    result = run(truth_root=args.truth_root, config_path=args.config, output_dir=args.output_dir)
    print(result["status"])
    return 0 if result["status"] == "COMPUTED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
