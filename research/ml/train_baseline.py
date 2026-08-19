from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
    precision_score,
    recall_score,
    confusion_matrix,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "research" / "ml" / "generated"


FEATURE_SETS = {
    "A_conventional": [
        "abs_logfc",
        "neg_log10_p",
    ],
    "B_robustness": [
        "source_direction_stability",
        "source_rank_stability",
    ],
    "C_combined": [
        "abs_logfc",
        "neg_log10_p",
        "source_direction_stability",
        "source_rank_stability",
    ],
}


def load_dataset(name: str) -> pd.DataFrame:
    df = pd.read_csv(
        DATA_DIR / f"{name}_deseq2_to_edger.csv"
    )

    df["abs_logfc"] = df[
        "source_log_fold_change"
    ].abs()

    clipped_p = df[
        "source_p_value"
    ].clip(lower=1e-300)

    df["neg_log10_p"] = -np.log10(clipped_p)

    return df


def evaluate(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
    feature_names: list[str],
):
    X_train = train_df[feature_names]
    y_train = train_df["target_method_significant"]

    X_test = test_df[feature_names]
    y_test = test_df["target_method_significant"]

    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            max_iter=2000,
            class_weight="balanced",
        ),
    )

    model.fit(X_train, y_train)

    probabilities = model.predict_proba(X_test)[:, 1]
    predictions = model.predict(X_test)

    return {
        "pr_auc": average_precision_score(
            y_test,
            probabilities,
        ),
        "roc_auc": roc_auc_score(
            y_test,
            probabilities,
        ),
        "precision": precision_score(
            y_test,
            predictions,
            zero_division=0,
        ),
        "recall": recall_score(
            y_test,
            predictions,
            zero_division=0,
        ),
        "confusion_matrix": confusion_matrix(
            y_test,
            predictions,
        ),
        "coefficients": model.named_steps[
            "logisticregression"
        ].coef_[0],
    }


def run_direction(
    train_name: str,
    test_name: str,
):
    train_df = load_dataset(train_name)
    test_df = load_dataset(test_name)

    print(
        f"\n===== {train_name.upper()} -> "
        f"{test_name.upper()} ====="
    )

    print(
        "Test positive prevalence:",
        round(
            test_df[
                "target_method_significant"
            ].mean(),
            4,
        ),
    )

    for model_name, features in FEATURE_SETS.items():
        result = evaluate(
            train_df,
            test_df,
            features,
        )

        print(f"\n--- {model_name} ---")
        print("features:", features)
        print(
            "PR-AUC:",
            round(result["pr_auc"], 4),
        )
        print(
            "ROC-AUC:",
            round(result["roc_auc"], 4),
        )
        print(
            "precision:",
            round(result["precision"], 4),
        )
        print(
            "recall:",
            round(result["recall"], 4),
        )
        print(
            "confusion matrix:"
        )
        print(result["confusion_matrix"])

        print("coefficients:")
        for feature, coefficient in zip(
            features,
            result["coefficients"],
        ):
            print(
                f"  {feature}: "
                f"{coefficient:.4f}"
            )


if __name__ == "__main__":
    run_direction(
        "airway",
        "pasilla",
    )

    run_direction(
        "pasilla",
        "airway",
    )