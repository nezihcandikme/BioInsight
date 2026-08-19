from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "research" / "ml" / "generated"


FEATURE_SETS = {
    "A_conventional": [
        "abs_logfc",
        "neg_log10_p",
    ],
    "D_plus_rank": [
        "abs_logfc",
        "neg_log10_p",
        "source_rank_stability",
    ],
}


def load_dataset(name: str) -> pd.DataFrame:
    df = pd.read_csv(
        DATA_DIR / f"{name}_deseq2_to_edger.csv"
    )

    df["abs_logfc"] = df["source_log_fold_change"].abs()

    p = df["source_p_value"].clip(lower=1e-300)
    df["neg_log10_p"] = -np.log10(p)

    # Exploratory borderline subset
    return df[
        (df["source_adjusted_p_value"] >= 0.01)
        & (df["source_adjusted_p_value"] <= 0.10)
    ].copy()


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

    coefficients = model.named_steps[
        "logisticregression"
    ].coef_[0]

    return {
        "pr_auc": average_precision_score(y_test, probs),
        "roc_auc": roc_auc_score(y_test, probs),
        "coefficients": coefficients,
    }


def run_direction(train_name, test_name):
    train_df = load_dataset(train_name)
    test_df = load_dataset(test_name)

    print(
        f"\n===== {train_name.upper()} -> "
        f"{test_name.upper()} ====="
    )

    for model_name, features in FEATURE_SETS.items():
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