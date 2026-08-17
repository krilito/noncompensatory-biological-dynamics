"""Label-free residualization and low-dimensional LOPO response models."""

from __future__ import annotations

import warnings
from typing import Any, Mapping, Sequence

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from .falsification_statistics import patient_cluster_bootstrap_indices
from .incremental_value import exact_sign_flip_distribution
from .statistics import auc, binary_log_loss


MODEL_COLUMNS = {
    "M_PRE": ("IFNG6_PRE",),
    "M_ON": ("IFNG6_ON",),
    "M_PRE_ON": ("IFNG6_PRE", "IFNG6_ON"),
    "M_ON_DELTA": ("IFNG6_ON", "IFNG6_RAW_DELTA"),
}


def residualized_change(
    pre: Sequence[float], on: Sequence[float], patient_ids: Sequence[str]
) -> dict[str, Any]:
    """Fit ON ~ PRE without response labels, both in-sample and patient-LOPO."""
    pre_values = np.asarray(pre, dtype=float)
    on_values = np.asarray(on, dtype=float)
    patients = np.asarray([str(value) for value in patient_ids], dtype=object)
    if not (len(pre_values) == len(on_values) == len(patients)):
        raise ValueError("residualization inputs must be record-aligned")

    design = np.column_stack([np.ones(len(pre_values)), pre_values])
    full_coefficients, *_ = np.linalg.lstsq(design, on_values, rcond=None)
    full_fitted = design @ full_coefficients
    full_residual = on_values - full_fitted

    lopo_residual = np.empty(len(pre_values), dtype=float)
    folds = []
    for held_patient in dict.fromkeys(patients.tolist()):
        train = patients != held_patient
        test = ~train
        fold_design = np.column_stack([np.ones(int(train.sum())), pre_values[train]])
        coefficients, *_ = np.linalg.lstsq(fold_design, on_values[train], rcond=None)
        predicted = coefficients[0] + coefficients[1] * pre_values[test]
        lopo_residual[test] = on_values[test] - predicted
        folds.append({
            "held_out_patient": held_patient,
            "n_training_records": int(train.sum()),
            "n_held_out_records": int(test.sum()),
            "alpha": float(coefficients[0]),
            "beta": float(coefficients[1]),
        })
    return {
        "full_sample": {
            "alpha": float(full_coefficients[0]),
            "beta": float(full_coefficients[1]),
            "r_squared": float(1.0 - np.sum(full_residual**2) / np.sum((on_values - np.mean(on_values)) ** 2)),
            "residuals": full_residual.tolist(),
        },
        "lopo": {"residuals": lopo_residual.tolist(), "folds": folds},
    }


def _lopo_predictions(
    rows: Sequence[Mapping[str, Any]], seed: int
) -> tuple[dict[str, np.ndarray], list[dict[str, Any]], list[str]]:
    patients = [str(row["patient_id"]) for row in rows]
    ordered_patients = list(dict.fromkeys(patients))
    labels = np.asarray([int(row["response_binary"]) for row in rows], dtype=int)
    predictions = {model: np.empty(len(rows), dtype=float) for model in MODEL_COLUMNS}
    warnings_seen: list[str] = []
    fold_rows: list[dict[str, Any]] = []

    for held_patient in ordered_patients:
        train_indices = [index for index, patient in enumerate(patients) if patient != held_patient]
        test_indices = [index for index, patient in enumerate(patients) if patient == held_patient]
        y_train = labels[train_indices]
        if set(y_train.tolist()) != {0, 1}:
            raise ValueError(f"LOPO training fold lacks both response classes: {held_patient}")
        for model, columns in MODEL_COLUMNS.items():
            x_train = np.asarray(
                [[float(rows[index][column]) for column in columns] for index in train_indices], dtype=float
            )
            x_test = np.asarray(
                [[float(rows[index][column]) for column in columns] for index in test_indices], dtype=float
            )
            scaler = StandardScaler()
            classifier = LogisticRegression(
                C=1.0, solver="liblinear", random_state=seed, max_iter=2000
            )
            with warnings.catch_warnings(record=True) as captured:
                warnings.simplefilter("always", ConvergenceWarning)
                fitted = classifier.fit(scaler.fit_transform(x_train), y_train)
            for warning in captured:
                if issubclass(warning.category, ConvergenceWarning):
                    warnings_seen.append(f"{held_patient}:{model}:{warning.message}")
            probability = np.clip(
                fitted.predict_proba(scaler.transform(x_test))[:, 1], 1e-6, 1.0 - 1e-6
            )
            predictions[model][test_indices] = probability
            fold_rows.append({
                "held_out_patient": held_patient,
                "model": model,
                "n_training_records": len(train_indices),
                "n_held_out_records": len(test_indices),
                "coefficient_l2_norm": float(np.linalg.norm(fitted.coef_)),
            })
    return predictions, fold_rows, warnings_seen


def _percentile(values: Sequence[float]) -> list[float]:
    array = np.asarray(values, dtype=float)
    return [float(np.percentile(array, 2.5)), float(np.percentile(array, 97.5))]


def lopo_incremental_information(
    rows: Sequence[Mapping[str, Any]], *, bootstrap_repeats: int, seed: int
) -> dict[str, Any]:
    """Compare fixed low-dimensional ridge-logistic models by patient-LOPO."""
    labels = np.asarray([int(row["response_binary"]) for row in rows], dtype=int)
    patient_ids = [str(row["patient_id"]) for row in rows]
    predictions, folds, convergence_warnings = _lopo_predictions(rows, seed)
    model_auc = {model: float(auc(values.tolist(), labels.tolist())) for model, values in predictions.items()}
    model_auc_draws = {model: [] for model in predictions}
    valid_draws = 0
    bootstrap_indices = patient_cluster_bootstrap_indices(patient_ids, bootstrap_repeats, seed + 1)
    for indices in bootstrap_indices:
        drawn_labels = labels[indices]
        if set(drawn_labels.tolist()) != {0, 1}:
            continue
        for model, values in predictions.items():
            model_auc_draws[model].append(float(auc(values[indices].tolist(), drawn_labels.tolist())))
        valid_draws += 1

    ordered_patients = list(dict.fromkeys(patient_ids))
    patient_losses: dict[str, dict[str, float]] = {}
    for patient in ordered_patients:
        indices = [index for index, value in enumerate(patient_ids) if value == patient]
        patient_losses[patient] = {
            model: binary_log_loss(values[indices].tolist(), labels[indices].tolist())
            for model, values in predictions.items()
        }
    model_summary = {
        model: {
            "oof_auc": model_auc[model],
            "oof_auc_patient_cluster_ci_95": _percentile(model_auc_draws[model]),
            "patient_equal_weight_log_loss": float(np.mean([patient_losses[p][model] for p in ordered_patients])),
        }
        for model in predictions
    }

    rng = np.random.default_rng(seed + 2)
    patient_draws = rng.choice(ordered_patients, size=(bootstrap_repeats, len(ordered_patients)), replace=True)

    def compare(baseline: str, augmented: str) -> dict[str, Any]:
        loss_improvements = np.asarray(
            [patient_losses[patient][baseline] - patient_losses[patient][augmented] for patient in ordered_patients],
            dtype=float,
        )
        boot_loss = [
            float(np.mean([patient_losses[str(patient)][baseline] - patient_losses[str(patient)][augmented] for patient in draw]))
            for draw in patient_draws
        ]
        auc_differences = [
            augmented_auc - baseline_auc
            for augmented_auc, baseline_auc in zip(model_auc_draws[augmented], model_auc_draws[baseline])
        ]
        loss_ci = _percentile(boot_loss)
        auc_ci = _percentile(auc_differences)
        adds = loss_ci[0] > 0 and auc_ci[0] > 0
        degrades = loss_ci[1] < 0 and auc_ci[1] < 0
        if adds:
            status = "ADDS_HELD_OUT_INFORMATION"
        elif degrades:
            status = "DEGRADES_HELD_OUT_INFORMATION"
        else:
            status = "INCREMENTAL_ANALYSIS_UNSTABLE"
        return {
            "baseline_model": baseline,
            "augmented_model": augmented,
            "patient_equal_weight_delta_log_loss": float(np.mean(loss_improvements)),
            "delta_log_loss_patient_bootstrap_ci_95": loss_ci,
            "auc_difference": model_auc[augmented] - model_auc[baseline],
            "auc_difference_patient_cluster_ci_95": auc_ci,
            "patient_sign_flip": exact_sign_flip_distribution(loss_improvements.tolist()),
            "status": status,
        }

    oof_rows = [
        {
            "patient_id": row["patient_id"],
            "therapy_record_id": row["therapy_record_id"],
            "response_binary": row["response_binary"],
            **{f"probability_{model}": float(predictions[model][index]) for model in predictions},
        }
        for index, row in enumerate(rows)
    ]
    return {
        "analysis_role": "small-sample decomposition sensitivity; not a prediction-model training claim",
        "model_contract": "patient-LOPO; train-fold StandardScaler; L2 logistic C=1.0; liblinear; fixed seed",
        "equal_patient_weighting_for_log_loss": True,
        "record_level_auc_with_patient_cluster_interval": True,
        "bootstrap_repeats_valid": valid_draws,
        "models": model_summary,
        "incremental_on_given_pre": compare("M_PRE", "M_PRE_ON"),
        "incremental_delta_given_on": compare("M_ON", "M_ON_DELTA"),
        "convergence_warnings": convergence_warnings,
        "fold_diagnostics": folds,
        "oof_predictions": oof_rows,
    }
