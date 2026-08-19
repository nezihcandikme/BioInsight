"""Build experimental per-gene ML tables from existing DE result files.

This module is research-only. It does not train models or split data.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class MethodSchema:
    gene_id: str
    log_fold_change: str
    p_value: str
    adjusted_p_value: str


METHOD_SCHEMAS = {
    "DESeq2": MethodSchema("gene_id", "log2FoldChange", "pvalue", "padj"),
    "edgeR": MethodSchema("gene_id", "logFC", "pvalue", "FDR"),
    "limma-voom": MethodSchema("gene_id", "logFC", "pvalue", "adj.P.Val"),
}


def _load_method_table(
    path: str | Path,
    schema: MethodSchema,
    method: str,
) -> pd.DataFrame:
    table = pd.read_csv(path, na_values=["NA"])

    required = [
        schema.gene_id,
        schema.log_fold_change,
        schema.p_value,
        schema.adjusted_p_value,
    ]

    missing = [
        column
        for column in required
        if column not in table.columns
    ]

    if missing:
        raise ValueError(
            f"{method} table is missing required columns: {missing}"
        )

    if table[schema.gene_id].duplicated().any():
        raise ValueError(
            f"{method} table contains duplicate gene IDs"
        )

    table = table[required].set_index(schema.gene_id)

    return table.rename(
        columns={
            schema.log_fold_change: "log_fold_change",
            schema.p_value: "p_value",
            schema.adjusted_p_value: "adjusted_p_value",
        }
    )


def _load_robustness_table(
    path: str | Path,
) -> pd.DataFrame:
    robustness = pd.read_csv(path)

    required = [
        "gene_id",
        "frac_significant",
        "direction_stability",
        "rank_stability",
    ]

    missing = [
        column
        for column in required
        if column not in robustness.columns
    ]

    if missing:
        raise ValueError(
            f"Robustness table is missing required columns: {missing}"
        )

    if robustness["gene_id"].duplicated().any():
        raise ValueError(
            "Robustness table contains duplicate gene IDs"
        )

    robustness = robustness[required].set_index("gene_id")

    return robustness.rename(
        columns={
            "frac_significant": "source_significance_stability",
            "direction_stability": "source_direction_stability",
            "rank_stability": "source_rank_stability",
        }
    )


def build_feature_table(
    source_path: str | Path,
    target_path: str | Path,
    *,
    dataset_id: str,
    source_method: str = "DESeq2",
    target_method: str = "edgeR",
    robustness_path: str | Path | None = None,
    alpha: float = 0.05,
    threshold_alphas: tuple[float, ...] = (0.01, 0.05, 0.1),
    schemas: dict[str, MethodSchema] | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Return an aligned source-feature/target table and alignment report."""

    if not 0 < alpha <= 1:
        raise ValueError("alpha must be in (0, 1]")

    if (
        not threshold_alphas
        or any(
            not 0 < value <= 1
            for value in threshold_alphas
        )
    ):
        raise ValueError(
            "threshold_alphas must contain values in (0, 1]"
        )

    schemas = METHOD_SCHEMAS if schemas is None else schemas

    if (
        source_method not in schemas
        or target_method not in schemas
    ):
        raise ValueError(
            "source_method and target_method must have configured schemas"
        )

    source_schema = schemas[source_method]
    target_schema = schemas[target_method]

    source = _load_method_table(
        source_path,
        source_schema,
        source_method,
    )

    target = _load_method_table(
        target_path,
        target_schema,
        target_method,
    )

    # Keep only genes present in both DE methods.
    overlap = source.index.intersection(
        target.index,
        sort=False,
    )

    aligned = source.loc[overlap].join(
        target.loc[overlap],
        how="inner",
        lsuffix="_source",
        rsuffix="_target",
    )

    source_required = [
        "log_fold_change_source",
        "p_value_source",
        "adjusted_p_value_source",
    ]

    target_adjusted = "adjusted_p_value_target"

    required_for_row = (
        source_required
        + [target_adjusted]
    )

    missing_required = aligned[
        required_for_row
    ].isna().any(axis=1)

    usable = aligned.loc[
        ~missing_required
    ].copy()

    result = pd.DataFrame(
        index=usable.index
    )

    result.index.name = "gene_id"

    # Metadata
    result["dataset_id"] = dataset_id
    result["source_method"] = source_method
    result["target_method"] = target_method

    # Conventional DE features
    result["source_log_fold_change"] = (
        usable["log_fold_change_source"]
    )

    result["source_p_value"] = (
        usable["p_value_source"]
    )

    result["source_adjusted_p_value"] = (
        usable["adjusted_p_value_source"]
    )

    # Existing threshold-derived feature.
    # This is retained for research comparison but should not be treated
    # as an independent robustness signal.
    result["source_threshold_significance_fraction"] = (
        sum(
            result["source_adjusted_p_value"]
            .lt(value)
            .astype(int)
            for value in threshold_alphas
        )
        / len(threshold_alphas)
    )

    # Target label
    result["target_method_significant"] = (
        usable[target_adjusted]
        .lt(alpha)
        .astype("int8")
    )

    # Optional DeConcord v0.18 robustness features
    if robustness_path is not None:
        robustness = _load_robustness_table(
            robustness_path
        )

        result = result.join(
            robustness,
            how="left",
        )

    result = result.reset_index()

    report = {
        "dataset_id": dataset_id,
        "source_method": source_method,
        "target_method": target_method,
        "source_gene_count": int(len(source)),
        "target_gene_count": int(len(target)),
        "overlap_gene_count": int(len(overlap)),
        "usable_gene_count": int(len(result)),
        "source_only_gene_count": int(
            len(
                source.index.difference(
                    target.index
                )
            )
        ),
        "target_only_gene_count": int(
            len(
                target.index.difference(
                    source.index
                )
            )
        ),
        "overlap_dropped_missing_required_count": int(
            missing_required.sum()
        ),
        "overlap_missing_by_required_field": {
            column: int(
                aligned[column]
                .isna()
                .sum()
            )
            for column in required_for_row
        },
        "positive_target_count": int(
            result[
                "target_method_significant"
            ].sum()
        ),
        "negative_target_count": int(
            (
                result[
                    "target_method_significant"
                ]
                == 0
            ).sum()
        ),
        "positive_prevalence": float(
            result[
                "target_method_significant"
            ].mean()
        ),
        "missingness_by_feature": {
            column: int(
                result[column]
                .isna()
                .sum()
            )
            for column in result.columns
        },
        "threshold_alphas": list(
            threshold_alphas
        ),
        "robustness_features_included": (
            robustness_path is not None
        ),
    }

    return result, report
