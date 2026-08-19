from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "research" / "ml" / "generated"

C_VALUES = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]


def load_dataset(name: str) -> pd.DataFrame:
    df = pd.read_csv(DATA_DIR / f"{name}_deseq2_to_edger.csv")

    df["abs_logfc"] = df["source_log_fold_change"].abs()
    p = df["source_p_value"].clip(lower=1e-300)
    df["neg_log10_p"] = -np.log10(p)

    return df[
        (df["source_adjusted_p_value"] >= 0.01)
        & (df["source_adjusted_p_value"] <= 0.10)
    ].copy()


def evaluate(train_df, test_df, features, C):
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            solver="liblinear",
            penalty="l1",
            C=C,
            class_weight="balanced",
            max_iter=3000,
        ),
    )

    model.fit(
        train_df[features],
        train_df["target_method_significant"],
    )

    probs = model.predict_proba(
        test_df[features]
    )[:, 1]

    return {
        "pr_auc": average_precision_score(
            test_df["target_method_significant"],
            probs,
        ),
        "roc_auc": roc_auc_score(
            test_df["target_method_significant"],
            probs,
        ),
        "coefficients": model.named_steps[
            "logisticregression"
        ].coef_[0],
    }


def run_direction(train_name, test_name):
    train_df = load_dataset(train_name)
    test_df = load_dataset(test_name)

    conventional = [
        "abs_logfc",
        "neg_log10_p",
    ]

    combined = [
        "abs_logfc",
        "neg_log10_p",
        "source_rank_stability",
    ]

    print(f"\n===== {train_name.upper()} -> {test_name.upper()} =====")

    for C in C_VALUES:
        a = evaluate(
            train_df,
            test_df,
            conventional,
            C,
        )

        d = evaluate(
            train_df,
            test_df,
            combined,
            C,
        )

        print(f"\nC = {C}")

        print(
            "Conventional PR-AUC:",
            round(a["pr_auc"], 4),
        )

        print(
            "+ rank PR-AUC:",
            round(d["pr_auc"], 4),
        )

        print(
            "Δ PR-AUC:",
            round(
                d["pr_auc"] - a["pr_auc"],
                4,
            ),
        )

        print(
            "Conventional ROC-AUC:",
            round(a["roc_auc"], 4),
        )

        print(
            "+ rank ROC-AUC:",
            round(d["roc_auc"], 4),
        )

        print(
            "Rank coefficient:",
            round(
                d["coefficients"][-1],
                4,
            ),
        )


if __name__ == "__main__":
    run_direction("airway", "pasilla")
    run_direction("pasilla", "airway")