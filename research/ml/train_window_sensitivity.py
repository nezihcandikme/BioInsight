from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "research" / "ml" / "generated"

WINDOWS = [
    (0.01, 0.05),
    (0.02, 0.08),
    (0.01, 0.10),
    (0.025, 0.075),
    (0.05, 0.10),
]


def load_dataset(name: str) -> pd.DataFrame:
    df = pd.read_csv(DATA_DIR / f"{name}_deseq2_to_edger.csv")
    df["abs_logfc"] = df["source_log_fold_change"].abs()
    p = df["source_p_value"].clip(lower=1e-300)
    df["neg_log10_p"] = -np.log10(p)
    return df


def evaluate(train_df, test_df, features):
    X_train = train_df[features]
    y_train = train_df["target_method_significant"]

    X_test = test_df[features]
    y_test = test_df["target_method_significant"]

    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            max_iter=2000,
            class_weight="balanced",
        ),
    )

    model.fit(X_train, y_train)
    probs = model.predict_proba(X_test)[:, 1]

    return {
        "pr_auc": average_precision_score(y_test, probs),
        "roc_auc": roc_auc_score(y_test, probs),
    }


def run_direction(train_name: str, test_name: str):
    train_all = load_dataset(train_name)
    test_all = load_dataset(test_name)

    print(f"\n===== {train_name.upper()} -> {test_name.upper()} =====")

    for low, high in WINDOWS:
        train_df = train_all[
            (train_all["source_adjusted_p_value"] >= low)
            & (train_all["source_adjusted_p_value"] <= high)
        ].copy()

        test_df = test_all[
            (test_all["source_adjusted_p_value"] >= low)
            & (test_all["source_adjusted_p_value"] <= high)
        ].copy()

        print(f"\nWindow: {low:.3f} - {high:.3f}")
        print("Train genes:", len(train_df))
        print("Test genes:", len(test_df))
        print(
            "Train prevalence:",
            round(train_df["target_method_significant"].mean(), 4),
        )
        print(
            "Test prevalence:",
            round(test_df["target_method_significant"].mean(), 4),
        )

        # Avoid invalid evaluation if a window contains only one class.
        if (
            train_df["target_method_significant"].nunique() < 2
            or test_df["target_method_significant"].nunique() < 2
        ):
            print("Skipped: only one target class present.")
            continue

        conventional = evaluate(
            train_df,
            test_df,
            ["abs_logfc", "neg_log10_p"],
        )

        plus_rank = evaluate(
            train_df,
            test_df,
            [
                "abs_logfc",
                "neg_log10_p",
                "source_rank_stability",
            ],
        )

        delta_pr = plus_rank["pr_auc"] - conventional["pr_auc"]
        delta_roc = plus_rank["roc_auc"] - conventional["roc_auc"]

        print(
            "A conventional PR-AUC:",
            round(conventional["pr_auc"], 4),
        )
        print(
            "+ rank PR-AUC:",
            round(plus_rank["pr_auc"], 4),
        )
        print(
            "Δ PR-AUC:",
            round(delta_pr, 4),
        )

        print(
            "A conventional ROC-AUC:",
            round(conventional["roc_auc"], 4),
        )
        print(
            "+ rank ROC-AUC:",
            round(plus_rank["roc_auc"], 4),
        )
        print(
            "Δ ROC-AUC:",
            round(delta_roc, 4),
        )


if __name__ == "__main__":
    run_direction("airway", "pasilla")
    run_direction("pasilla", "airway")