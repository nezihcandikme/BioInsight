from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "research" / "ml" / "generated"


def load_dataset(name: str) -> pd.DataFrame:
    df = pd.read_csv(DATA_DIR / f"{name}_deseq2_to_edger.csv")

    df["abs_logfc"] = df["source_log_fold_change"].abs()

    p = df["source_p_value"].clip(lower=1e-300)
    df["neg_log10_p"] = -np.log10(p)

    return df[
        (df["source_adjusted_p_value"] >= 0.01)
        & (df["source_adjusted_p_value"] <= 0.10)
    ].copy()


def add_residual_rank(train_df, test_df):
    base_features = ["abs_logfc", "neg_log10_p"]

    rank_model = make_pipeline(
        StandardScaler(),
        LinearRegression(),
    )

    rank_model.fit(
        train_df[base_features],
        train_df["source_rank_stability"],
    )

    train_pred = rank_model.predict(train_df[base_features])
    test_pred = rank_model.predict(test_df[base_features])

    train_df = train_df.copy()
    test_df = test_df.copy()

    train_df["residual_rank_stability"] = (
        train_df["source_rank_stability"] - train_pred
    )

    test_df["residual_rank_stability"] = (
        test_df["source_rank_stability"] - test_pred
    )

    return train_df, test_df


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
        "coefficients": model.named_steps[
            "logisticregression"
        ].coef_[0],
    }


def run_direction(train_name, test_name):
    train_df = load_dataset(train_name)
    test_df = load_dataset(test_name)

    train_df, test_df = add_residual_rank(
        train_df,
        test_df,
    )

    feature_sets = {
        "A_conventional": [
            "abs_logfc",
            "neg_log10_p",
        ],
        "E_plus_residual_rank": [
            "abs_logfc",
            "neg_log10_p",
            "residual_rank_stability",
        ],
    }

    print(
        f"\n===== {train_name.upper()} -> "
        f"{test_name.upper()} ====="
    )

    print(
        "Residual rank correlation with target:",
        round(
            test_df[
                [
                    "residual_rank_stability",
                    "target_method_significant",
                ]
            ]
            .corr(method="spearman")
            .iloc[0, 1],
            4,
        ),
    )

    for model_name, features in feature_sets.items():
        result = evaluate(
            train_df,
            test_df,
            features,
        )

        print(f"\n--- {model_name} ---")
        print("PR-AUC:", round(result["pr_auc"], 4))
        print("ROC-AUC:", round(result["roc_auc"], 4))

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
    run_direction("airway", "pasilla")
    run_direction("pasilla", "airway")
    