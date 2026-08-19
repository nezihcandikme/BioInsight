from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "research" / "ml" / "generated"


FEATURES = [
    "abs_logfc",
    "neg_log10_p",
    "source_rank_stability",
]


def load_dataset(name: str) -> pd.DataFrame:
    df = pd.read_csv(
        DATA_DIR / f"{name}_deseq2_to_edger.csv"
    )

    df["abs_logfc"] = df["source_log_fold_change"].abs()

    p = df["source_p_value"].clip(lower=1e-300)
    df["neg_log10_p"] = -np.log10(p)

    return df[
        (df["source_adjusted_p_value"] >= 0.01)
        & (df["source_adjusted_p_value"] <= 0.10)
    ].copy()


def run_direction(train_name: str, test_name: str) -> None:
    train_df = load_dataset(train_name)
    test_df = load_dataset(test_name)

    X_train = train_df[FEATURES]
    y_train = train_df["target_method_significant"]

    X_test = test_df[FEATURES]
    y_test = test_df["target_method_significant"]

    print(f"\n===== {train_name.upper()} -> {test_name.upper()} =====")

    for C in [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]:
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(
                penalty="l1",
                solver="liblinear",
                C=C,
                class_weight="balanced",
                max_iter=3000,
            ),
        )

        model.fit(X_train, y_train)

        probs = model.predict_proba(X_test)[:, 1]

        coefficients = model.named_steps[
            "logisticregression"
        ].coef_[0]

        print(f"\nC = {C}")
        print(
            "PR-AUC:",
            round(
                average_precision_score(y_test, probs),
                4,
            ),
        )
        print(
            "ROC-AUC:",
            round(
                roc_auc_score(y_test, probs),
                4,
            ),
        )

        print("coefficients:")
        for feature, coefficient in zip(FEATURES, coefficients):
            print(f"  {feature}: {coefficient:.4f}")


if __name__ == "__main__":
    run_direction("airway", "pasilla")
    run_direction("pasilla", "airway")