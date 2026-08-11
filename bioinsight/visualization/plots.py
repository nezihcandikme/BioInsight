import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

# Adjusted p-values of exactly 0 are possible (e.g. a constant-expression
# gene with a clear between-group difference — see compute_pvalues). -log10(0)
# is -inf, which would otherwise stretch the y-axis to infinity and print a
# RuntimeWarning. Flooring at a tiny epsilon keeps the plot readable without
# hiding how significant the point actually is.
_MIN_PLOTTABLE_PVALUE = 1e-300


def plot_volcano(results_df: pd.DataFrame) -> plt.Figure:
    """
    Volcano plot: log fold change (x) vs. -log10 adjusted p-value (y).

    Parameters
    ----------
    results_df : pd.DataFrame
        Output of ``run_differential_expression`` — needs ``log_fold_change``,
        ``adjusted_p_value``, and ``significant`` columns.

    Returns
    -------
    plt.Figure
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axhline(y=-np.log10(0.05), color="black", linestyle="--", label="Significance threshold (0.05)")
    ax.set_title("Volcano Plot")

    significant = results_df["significant"]
    plottable_pvalues = results_df["adjusted_p_value"].clip(lower=_MIN_PLOTTABLE_PVALUE)

    ax.scatter(
        results_df.loc[~significant, "log_fold_change"],
        -np.log10(plottable_pvalues.loc[~significant]),
        alpha=0.5, color="gray", label="Not significant",
    )
    ax.scatter(
        results_df.loc[significant, "log_fold_change"],
        -np.log10(plottable_pvalues.loc[significant]),
        alpha=0.5, color="red", label="Significant",
    )

    ax.axvline(1, color="black", linestyle="--")
    ax.axvline(-1, color="black", linestyle="--")
    ax.legend()
    ax.set_xlabel("Log Fold Change")
    ax.set_ylabel("-Log10 Adjusted P-value")
    return fig


def plot_pca(df: pd.DataFrame, labels: pd.Series) -> plt.Figure:
    """
    2-component PCA scatter plot over samples, colored by ``labels``.

    Parameters
    ----------
    df : pd.DataFrame
        Expression matrix, genes as rows, samples as columns. Typically
        the normalized, log-transformed matrix rather than raw counts,
        so a handful of very highly expressed genes don't dominate the
        variance PCA is built to summarize.
    labels : pd.Series
        Group label per sample, indexed to match ``df.columns``.

    Returns
    -------
    plt.Figure

    Raises
    ------
    ValueError
        If ``df`` has fewer than 2 samples or fewer than 2 genes — there's
        no such thing as 2 principal components of less than 2 dimensions,
        and scikit-learn's own error for this is not exactly self-explanatory.
    """
    n_genes, n_samples = df.shape
    if n_samples < 2 or n_genes < 2:
        raise ValueError(
            f"plot_pca needs at least 2 samples and 2 genes to compute 2 "
            f"principal components; got {n_samples} sample(s) and "
            f"{n_genes} gene(s)."
        )

    pca = PCA(n_components=2)
    coords = pca.fit_transform(df.T)
    variance_ratio = pca.explained_variance_ratio_

    fig, ax = plt.subplots(figsize=(8, 6))
    for label in labels.unique():
        mask = (labels == label).values
        ax.scatter(coords[mask, 0], coords[mask, 1], label=str(label), alpha=0.7)
    ax.legend()

    for i, sample_name in enumerate(df.columns):
        ax.annotate(sample_name, (coords[i, 0], coords[i, 1]))

    ax.set_title("PCA Plot")
    ax.set_xlabel(f"Principal Component 1 ({variance_ratio[0]*100:.1f}% variance)")
    ax.set_ylabel(f"Principal Component 2 ({variance_ratio[1]*100:.1f}% variance)")

    return fig
