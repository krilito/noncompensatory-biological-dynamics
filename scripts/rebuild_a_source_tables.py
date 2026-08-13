"""Rebuild active A-chain Source Data tables from the live producer.

Writes repository-relative figure2_auc and ST2 external-support rows from
run_truth_cohorts output. Historical AUC columns are retained only as labeled
provenance drift fields when present in the previous table; they are never used
as active estimates.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _ci(metrics: dict[str, Any]) -> tuple[float | None, float | None]:
    ci = metrics.get("ci_95")
    if isinstance(ci, (list, tuple)) and len(ci) == 2:
        return float(ci[0]), float(ci[1])
    return None, None


def _labels(metrics: dict[str, Any]) -> tuple[int | None, int | None]:
    labels = metrics.get("labels") or {}
    r = labels.get("responder")
    nr = labels.get("nonresponder")
    return (int(r) if r is not None else None, int(nr) if nr is not None else None)


def rebuild(*, truth_root: Path, output_json: Path | None = None) -> dict[str, Any]:
    from meld_icb.frozen_state_inputs import run_truth_cohorts

    audit = ROOT / "outputs" / "shared_chains" / "A_cohorts"
    audit.mkdir(parents=True, exist_ok=True)
    result = run_truth_cohorts(
        truth_root=truth_root,
        audit_dir=audit,
        config_path=ROOT / "configs" / "frozen_state_transfer.yaml",
        historical_comparison=None,
    )
    out_json = output_json or (ROOT / "outputs" / "shared_chains" / "A_PUBLIC_NUMERIC_OUTPUT.json")
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")

    cohorts = result.get("cohorts") or {}
    prjeb = (cohorts.get("PRJEB23709") or {}).get("metrics") or {}
    morrison = (cohorts.get("MORRISON-1-public") or {}).get("metrics") or {}
    mgh = (cohorts.get("MGH") or {}).get("metrics") or {}
    mgh_patient = ((cohorts.get("MGH") or {}).get("patient_metrics")
                   or (cohorts.get("MGH") or {}).get("patient_level_metrics")
                   or {})
    # Some producers nest patient sensitivity under populations or sibling keys.
    if not mgh_patient:
        for key, value in (cohorts.get("MGH") or {}).items():
            if isinstance(value, dict) and value.get("analysis_unit") == "patient":
                mgh_patient = value
                break
            if key.lower().startswith("patient") and isinstance(value, dict) and "auc" in value:
                mgh_patient = value
                break

    # Prefer nested structure used by current producer if present.
    mgh_blob = cohorts.get("MGH") or {}
    if isinstance(mgh_blob.get("patient"), dict) and "auc" in mgh_blob["patient"]:
        mgh_patient = mgh_blob["patient"]
    if isinstance(mgh_blob.get("metrics_patient"), dict):
        mgh_patient = mgh_blob["metrics_patient"]

    # Read previous historical columns for drift audit only.
    prev_path = ROOT / "source_data" / "figure2_auc.csv"
    hist: dict[str, dict[str, str]] = {}
    if prev_path.exists():
        with prev_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                hist[row["cohort"]] = row

    def row_for(name: str, metrics: dict[str, Any], *, analysis_unit: str, role: str, display: str) -> dict[str, Any]:
        lo, hi = _ci(metrics)
        n_r, n_nr = _labels(metrics)
        n = metrics.get("n_records") or metrics.get("n_patients") or metrics.get("n")
        prev = hist.get(display) or hist.get(name) or {}
        return {
            "cohort": display,
            "analysis_unit": analysis_unit,
            "n": n,
            "n_R": n_r if n_r is not None else "",
            "n_NR": n_nr if n_nr is not None else "",
            "auc_canonical": metrics.get("auc"),
            "p_canonical": metrics.get("p_two_sided"),
            "auc_historical": prev.get("auc_historical", ""),
            "p_historical": prev.get("p_historical", ""),
            "ci_low": lo if lo is not None else "",
            "ci_high": hi if hi is not None else "",
            "role": role,
        }

    # If patient metrics missing from this run, keep previous patient row AUC/p only if
    # we can recompute — fail closed instead.
    if not mgh_patient or mgh_patient.get("auc") is None:
        # Try A_MGH.json audit file
        mgh_audit = audit / "A_MGH.json"
        if mgh_audit.exists():
            blob = json.loads(mgh_audit.read_text(encoding="utf-8"))
            for key in ("patient_metrics", "patient", "metrics_patient", "patient_level"):
                if isinstance(blob.get(key), dict) and blob[key].get("auc") is not None:
                    mgh_patient = blob[key]
                    break
            if not mgh_patient:
                # search nested
                for key, value in blob.items():
                    if isinstance(value, dict) and value.get("analysis_unit") == "patient" and "auc" in value:
                        mgh_patient = value
                        break

    f2_rows = [
        row_for("PRJEB23709", prjeb, analysis_unit="biopsy sample", role="core", display="PRJEB23709"),
        row_for("MORRISON-1-public", morrison, analysis_unit="biopsy sample", role="core", display="MORRISON-1-public"),
        row_for("MGH", mgh, analysis_unit="biopsy sample", role="supportive", display="MGH sample-level"),
    ]
    if mgh_patient and mgh_patient.get("auc") is not None:
        f2_rows.append(
            row_for(
                "MGH patient-level",
                mgh_patient,
                analysis_unit="patient",
                role="supportive_different_unit",
                display="MGH patient-level",
            )
        )
    else:
        raise RuntimeError("MGH patient-level metrics missing from A producer output; cannot rebuild Source Data")

    # Write figure2_auc to active path only (source_data/figure2_auc.csv)
    f2_path = ROOT / "source_data" / "figure2_auc.csv"
    fieldnames = [
        "cohort", "analysis_unit", "n", "n_R", "n_NR", "auc_canonical", "p_canonical",
        "auc_historical", "p_historical", "ci_low", "ci_high", "role",
    ]
    with f2_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in f2_rows:
            writer.writerow(row)

    # Rebuild ST2 external-support block; keep non-external rows from previous file.
    st2_path = ROOT / "source_data" / "supplementary_tables" / "table2.csv"
    other_rows: list[dict[str, str]] = []
    if st2_path.exists():
        with st2_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if row.get("section") != "External response-state support":
                    other_rows.append(row)

    def st2_row(cohort: str, level: str, metrics: dict[str, Any], role: str) -> dict[str, Any]:
        lo, hi = _ci(metrics)
        n_r, n_nr = _labels(metrics)
        n = metrics.get("n_records") or metrics.get("n_patients") or metrics.get("n")
        ci = f"[{lo}, {hi}]" if lo is not None and hi is not None else ""
        return {
            "section": "External response-state support",
            "cohort": cohort,
            "analysis_level": level,
            "score": "MELD-ICB boundary",
            "n": n,
            "n_R": n_r if n_r is not None else "",
            "n_NR": n_nr if n_nr is not None else "",
            "AUC": metrics.get("auc"),
            "AUC_95_CI": ci,
            "permutation_P": metrics.get("p_two_sided"),
            "evidence_role": role,
        }

    pooled = result.get("pooled") or {}
    external = [
        st2_row("PRJEB23709 / Gide", "Sample", prjeb, "External validation"),
        st2_row("MORRISON-1-public", "Sample", morrison, "External validation"),
        st2_row("MGH", "Sample", mgh, "Directional external validation"),
        st2_row("MGH", "Patient", mgh_patient, "Sensitivity analysis"),
        st2_row("PRJEB23709 + MORRISON-1-public + MGH", "Biopsy", pooled, "Pooled summary"),
    ]
    st2_fields = [
        "section", "cohort", "analysis_level", "score", "n", "n_R", "n_NR",
        "AUC", "AUC_95_CI", "permutation_P", "evidence_role",
    ]
    with st2_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=st2_fields)
        writer.writeheader()
        for row in external:
            writer.writerow(row)
        for row in other_rows:
            writer.writerow(row)

    summary = {
        "PRJEB23709": {"auc": prjeb.get("auc"), "ci": _ci(prjeb), "p": prjeb.get("p_two_sided")},
        "MORRISON-1-public": {"auc": morrison.get("auc"), "ci": _ci(morrison), "p": morrison.get("p_two_sided")},
        "MGH_sample": {"auc": mgh.get("auc"), "ci": _ci(mgh), "p": mgh.get("p_two_sided")},
        "MGH_patient": {"auc": mgh_patient.get("auc"), "ci": _ci(mgh_patient), "p": mgh_patient.get("p_two_sided")},
        "pooled": {"auc": pooled.get("auc"), "ci": _ci(pooled), "p": pooled.get("p_two_sided")},
        "wrote": [str(f2_path.relative_to(ROOT)), str(st2_path.relative_to(ROOT)), str(out_json.relative_to(ROOT))],
    }
    print(json.dumps(summary, indent=2, default=str))
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-root", type=Path, required=True)
    args = parser.parse_args(argv)
    rebuild(truth_root=args.truth_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
