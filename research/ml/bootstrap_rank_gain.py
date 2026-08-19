from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
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

N_BOOTSTRAPS = 1000
RANDOM_SEED = 42


def load_dataset(name):
    df = pd.read_csv(DATA_DIR / f"{name}_deseq2_to_edger.csv")

    df["abs_logfc"] = df["source_log_fold_change"].abs()

    p = df["source_p_value"].clip(lower=1e-300)
    df["neg_log10_p"] = -np.log10(p)

    return df


def fit_predict(train_df, test_df, features):
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            max_iter=2000,
            class_weight="balanced",
        ),
    )

    model.fit(
        train_df[features],
        train_df["target_method_significant"],
    )

    return model.predict_proba(test_df[features])[:, 1]


def bootstrap_delta(y, probs_a, probs_rank):
    rng = np.random.default_rng(RANDOM_SEED)

    deltas = []

    n = len(y)

    for _ in range(N_BOOTSTRAPS):
        idx = rng.integers(0, n, size=n)

        y_boot = y[idx]

        # PR-AUC needs both classes present.
        if len(np.unique(y_boot)) < 2:
            continue

        a = average_precision_score(
            y_boot,
            probs_a[idx],
        )

        rank = average_precision_score(
            y_boot,
            probs_rank[idx],
        )

        deltas.append(rank - a)

    deltas = np.array(deltas)

    return {
        "mean": deltas.mean(),
        "low": np.percentile(deltas, 2.5),
        "high": np.percentile(deltas, 97.5),
        "positive_fraction": np.mean(deltas > 0),
    }


def run_direction(train_name, test_name):
    train_all = load_dataset(train_name)
    test_all = load_dataset(test_name)

    print(f"\n===== {train_name.upper()} -> {test_name.upper()} =====")

    for low, high in WINDOWS:
        train_df = train_all[
            train_all["source_adjusted_p_value"].between(low, high)
        ].copy()

        test_df = test_all[
            test_all["source_adjusted_p_value"].between(low, high)
        ].copy()

        if (
            train_df["target_method_significant"].nunique() < 2
            or test_df["target_method_significant"].nunique() < 2
        ):
            continue

        conventional_features = [
            "abs_logfc",
            "neg_log10_p",
        ]

        rank_features = [
            "abs_logfc",
            "neg_log10_p",
            "source_rank_stability",
        ]

        probs_a = fit_predict(
            train_df,
            test_df,
            conventional_features,
        )

        probs_rank = fit_predict(
            train_df,
            test_df,
            rank_features,
        )

        y = test_df["target_method_significant"].to_numpy()

        result = bootstrap_delta(
            y,
            probs_a,
            probs_rank,
        )

        print(f"\nWindow {low:.3f} - {high:.3f}")
        print(
            "Mean Δ PR-AUC:",
            round(result["mean"], 4),
        )
        print(
            "95% bootstrap CI:",
            (
                round(result["low"], 4),
                round(result["high"], 4),
            ),
        )
        print(
            "Bootstrap samples with Δ > 0:",
            round(result["positive_fraction"], 4),
        )


if __name__ == "__main__":
    run_direction("airway", "pasilla")
    run_direction("pasilla", "airway")