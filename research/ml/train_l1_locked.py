from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "research" / "ml" / "generated"

C_VALUES = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]


def load_dataset(name):
    df = pd.read_csv(DATA_DIR / f"{name}_deseq2_to_edger.csv")

    df["abs_logfc"] = df["source_log_fold_change"].abs()
    df["neg_log10_p"] = -np.log10(
        df["source_p_value"].clip(lower=1e-300)
    )

    return df[
        (df["source_adjusted_p_value"] >= 0.01)
        & (df["source_adjusted_p_value"] <= 0.10)
    ].copy()


def make_model(C):
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(
            solver="liblinear",
            penalty="l1",
            C=C,
            class_weight="balanced",
            max_iter=3000,
        ),
    )


def choose_c(train_df, features):
    X = train_df[features]
    y = train_df["target_method_significant"]

    cv = StratifiedKFold(
        n_splits=5,
        shuffle=True,
        random_state=42,
    )

    scores = {}

    for C in C_VALUES:
        fold_scores = []

        for train_idx, val_idx in cv.split(X, y):
            model = make_model(C)

            model.fit(
                X.iloc[train_idx],
                y.iloc[train_idx],
            )

            probs = model.predict_proba(
                X.iloc[val_idx]
            )[:, 1]

            score = average_precision_score(
                y.iloc[val_idx],
                probs,
            )

            fold_scores.append(score)

        scores[C] = np.mean(fold_scores)

    best_c = max(scores, key=scores.get)

    return best_c, scores


def evaluate_locked(train_df, test_df, features, C):
    model = make_model(C)

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

    print(
        f"\n===== {train_name.upper()} -> "
        f"{test_name.upper()} ====="
    )

    for label, features in [
        ("Conventional", conventional),
        ("+ Rank", combined),
    ]:
        best_c, cv_scores = choose_c(
            train_df,
            features,
        )

        result = evaluate_locked(
            train_df,
            test_df,
            features,
            best_c,
        )

        print(f"\n{label}")
        print("CV scores:")

        for C, score in cv_scores.items():
            print(
                f"  C={C:<5} "
                f"mean CV PR-AUC={score:.4f}"
            )

        print("Selected C:", best_c)
        print(
            "Locked test PR-AUC:",
            round(result["pr_auc"], 4),
        )
        print(
            "Locked test ROC-AUC:",
            round(result["roc_auc"], 4),
        )

        print("Coefficients:")

        for feature, coef in zip(
            features,
            result["coefficients"],
        ):
            print(
                f"  {feature}: {coef:.4f}"
            )


def _available_datasets():
    # Auto-discovers whatever *_deseq2_to_edger.csv tables build_tables.py
    # has actually produced, so this script covers every dataset currently
    # available (2 today, 3 once zebrafish's real R output lands) without
    # needing another edit -- the frozen experiment design should apply to
    # every dataset pair, not just the two it happened to be written
    # against first.
    names = []
    for path in sorted(DATA_DIR.glob("*_deseq2_to_edger.csv")):
        names.append(path.name.removesuffix("_deseq2_to_edger.csv"))
    return names


if __name__ == "__main__":
    datasets = _available_datasets()

    for train_name in datasets:
        for test_name in datasets:
            if train_name != test_name:
                run_direction(train_name, test_name)