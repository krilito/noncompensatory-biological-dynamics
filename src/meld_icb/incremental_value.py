"""Held-out incremental-value estimands."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from .statistics import binary_log_loss


def delta_loss(baseline_probabilities: Sequence[float], augmented_probabilities: Sequence[float], labels: Sequence[int]) -> float:
    """Compute loss(baseline) - loss(augmented); negative means adding change hurt."""
    return binary_log_loss(baseline_probabilities, labels) - binary_log_loss(augmented_probabilities, labels)


def exact_sign_flip_distribution(differences: Sequence[float], *, tolerance: float = 1e-15) -> dict[str, Any]:
    """Enumerate the frozen 2^15 patient sign-flip null exactly."""
    import itertools
    import numpy as np

    values = np.asarray(differences, dtype=float)
    if values.size == 0:
        raise ValueError("sign flips require non-empty patient differences")
    null_means = np.asarray([float(np.mean(values * np.asarray(signs, dtype=float))) for signs in itertools.product((1.0, -1.0), repeat=len(values))])
    observed = float(np.mean(values))
    return {
        "n_sign_patterns": int(len(null_means)),
        "P_improvement": float(np.mean(null_means >= observed - tolerance)),
        "P_worsening": float(np.mean(null_means <= observed + tolerance)),
        "P_two_sided": float(np.mean(np.abs(null_means) >= abs(observed) - tolerance)),
    }


def fold_safe_incremental_metrics(
    records: Sequence[Mapping[str, Any]],
    raw_axes: Mapping[str, Mapping[str, float]],
    freeze: Mapping[str, Any],
    *,
    seed: int = 42,
    canonical_scope: bool = True,
) -> dict[str, Any]:
    """Recompute B4 with patient-LOPO normalization and downstream scaling."""
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler
    from .frozen_state_transfer import frozen_boundary_score, fit_axis_normalization, transform_axis_normalization

    patients: list[Any] = []
    for record in records:
        if record["patient_id"] not in patients:
            patients.append(record["patient_id"])
    if canonical_scope:
        if len(records) != 16 or len(patients) != 15:
            raise ValueError(f"B4 requires 16 records/15 patients, got {len(records)}/{len(patients)}")
    elif len(patients) < 2:
        raise ValueError(f"B4 influence reruns require at least two patients, got {len(patients)}")
    sample_ids = [str(record[key]) for record in records for key in ("pre_sample_id", "edt_sample_id")]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("B4 sample leakage: a PRE/EDT sample is bound to more than one record")
    oof: list[dict[str, Any]] = []
    losses: list[float] = []
    loss_x0_values: list[float] = []
    loss_x1_values: list[float] = []
    for held_patient in sorted(patients):
        training_records = [record for record in records if record["patient_id"] != held_patient]
        heldout_records = [record for record in records if record["patient_id"] == held_patient]
        training_ids = [sample_id for record in training_records for sample_id in (str(record["pre_sample_id"]), str(record["edt_sample_id"]))]
        held_ids = [sample_id for record in heldout_records for sample_id in (str(record["pre_sample_id"]), str(record["edt_sample_id"]))]
        parameters = fit_axis_normalization([raw_axes[sample_id] for sample_id in training_ids], training_ids, scope_id="B4_fold_training_patients_PRE_EDT", ddof=1)
        leaked = set(held_ids) & set(parameters.fitted_sample_ids)
        if leaked:
            raise ValueError(f"B4 fold leakage: held-out samples entered normalization fit: {sorted(leaked)}")
        fold_ids = list(dict.fromkeys(training_ids + [sample_id for record in heldout_records for sample_id in (str(record["pre_sample_id"]), str(record["edt_sample_id"]))]))
        fold_z = {sample_id: row for sample_id, row in zip(fold_ids, transform_axis_normalization([raw_axes[sample_id] for sample_id in fold_ids], parameters))}
        fold_rows: list[dict[str, Any]] = []
        for record in records:
            pre_score = frozen_boundary_score(fold_z[str(record["pre_sample_id"])], freeze)
            on_score = frozen_boundary_score(fold_z[str(record["edt_sample_id"])], freeze)
            fold_rows.append({**record, "pre_score": pre_score, "on_score": on_score, "delta_score": on_score - pre_score})
        train_rows = [row for row in fold_rows if row["patient_id"] != held_patient]
        test_rows = [row for row in fold_rows if row["patient_id"] == held_patient]
        y_train = np.asarray([int(row["y_true"]) for row in train_rows], dtype=int)
        y_test = np.asarray([int(row["y_true"]) for row in test_rows], dtype=int)

        def predict(columns: list[str]) -> np.ndarray:
            x_train = np.asarray([[float(row[column]) for column in columns] for row in train_rows], dtype=float)
            x_test = np.asarray([[float(row[column]) for column in columns] for row in test_rows], dtype=float)
            if len(np.unique(y_train)) < 2:
                return np.full(len(test_rows), 0.5)
            scaler = StandardScaler()
            classifier = LogisticRegression(C=1.0, solver="liblinear", random_state=seed, max_iter=2000)
            return classifier.fit(scaler.fit_transform(x_train), y_train).predict_proba(scaler.transform(x_test))[:, 1]

        probabilities_x0 = np.clip(predict(["pre_score"]), 1e-6, 1 - 1e-6)
        probabilities_x1 = np.clip(predict(["pre_score", "delta_score"]), 1e-6, 1 - 1e-6)
        loss_x0 = -y_test * np.log(probabilities_x0) - (1 - y_test) * np.log(1 - probabilities_x0)
        loss_x1 = -y_test * np.log(probabilities_x1) - (1 - y_test) * np.log(1 - probabilities_x1)
        loss_x0_patient = float(np.mean(loss_x0))
        loss_x1_patient = float(np.mean(loss_x1))
        loss_x0_values.append(loss_x0_patient)
        loss_x1_values.append(loss_x1_patient)
        losses.append(loss_x0_patient - loss_x1_patient)
        for index, row in enumerate(test_rows):
            oof.append({"y_true": int(row["y_true"]), "probability_X0": float(probabilities_x0[index]), "probability_X1": float(probabilities_x1[index])})
    loss_x0 = float(np.mean(loss_x0_values))
    loss_x1 = float(np.mean(loss_x1_values))
    delta_l = loss_x0 - loss_x1
    flips = exact_sign_flip_distribution(losses)
    auc_x0 = float(roc_auc_score([row["y_true"] for row in oof], [row["probability_X0"] for row in oof]))
    auc_x1 = float(roc_auc_score([row["y_true"] for row in oof], [row["probability_X1"] for row in oof]))
    held_order = sorted(patients)
    return {
        "analysis_id": "B4_CANONICAL_FOLDSAFE_PRIMARY" if canonical_scope else "B4_FOLDSAFE_RERUN",
        "n_patients": len(patients),
        "n_records": len(records),
        "equal_patient_weighting": True,
        "axis_normalization_train_fold_only": True,
        "L_baseline_X0_pre": loss_x0,
        "L_augmented_X1_pre_plus_delta": loss_x1,
        "delta_L": delta_l,
        **flips,
        "oof_auc_X0": auc_x0,
        "oof_auc_X1": auc_x1,
        "secondary_delta_auc": auc_x1 - auc_x0,
        "canonical_status": "SUPPORTED" if delta_l > 0 and flips["P_improvement"] <= 0.05 else "NOT_SUPPORTED",
        "loss_sign_convention": "positive delta_L favors augmented; negative delta_L means worsening",
        "patient_delta_L": [
            {
                "patient_id": patient_id,
                "L_baseline": float(loss_x0_values[index]),
                "L_augmented": float(loss_x1_values[index]),
                "delta_L": float(losses[index]),
            }
            for index, patient_id in enumerate(held_order)
        ],
    }


def influence_leave_one_patient_reruns(
    records: Sequence[Mapping[str, Any]],
    raw_axes: Mapping[str, Mapping[str, float]],
    freeze: Mapping[str, Any],
    *,
    seed: int = 42,
) -> dict[str, Any]:
    """Re-run the full fold-safe procedure after omitting each patient.

    This is not deletion of one entry from a finished result array. Each rerun
    rebuilds remaining folds, refits axis normalization on the remaining
    training patients, and refits both logistic models.
    """
    import numpy as np

    patients: list[Any] = []
    for record in records:
        if record["patient_id"] not in patients:
            patients.append(record["patient_id"])
    rows: list[dict[str, Any]] = []
    for omitted in patients:
        remaining = [record for record in records if record["patient_id"] != omitted]
        result = fold_safe_incremental_metrics(
            remaining,
            raw_axes,
            freeze,
            seed=seed,
            canonical_scope=False,
        )
        rows.append(
            {
                "omitted_patient": omitted,
                "n_remaining": result["n_patients"],
                "n_records_remaining": result["n_records"],
                "mean_delta_L": result["delta_L"],
                "median_delta_L": float(np.median([item["delta_L"] for item in result["patient_delta_L"]])),
                "L_baseline": result["L_baseline_X0_pre"],
                "L_augmented": result["L_augmented_X1_pre_plus_delta"],
                "P_two_sided": result["P_two_sided"],
                "P_worsening": result["P_worsening"],
                "P_improvement": result["P_improvement"],
                "n_sign_patterns": result["n_sign_patterns"],
                "oof_auc_X0": result["oof_auc_X0"],
                "oof_auc_X1": result["oof_auc_X1"],
                "mean_changes_sign": False,
            }
        )
    means = [float(row["mean_delta_L"]) for row in rows]
    return {
        "analysis_id": "B4_INFLUENCE_LEAVE_ONE_PATIENT_RERUN",
        "rebuilds_folds": True,
        "refits_normalization": True,
        "refits_models": True,
        "n_reruns": len(rows),
        "all_means_negative": bool(means) and all(value < 0 for value in means),
        "mean_delta_L_min": min(means) if means else None,
        "mean_delta_L_max": max(means) if means else None,
        "rows": rows,
    }
