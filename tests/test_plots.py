import numpy as np
import pandas as pd
import pytest

from omicforge.differential_expression.methods import run_differential_expression
from omicforge.visualization.plots import plot_volcano, plot_pca


def test_plot_volcano_returns_figure():
    df = pd.DataFrame({
        "sample1": [10, 5], "sample2": [20, 15],
        "sample3": [30, 25], "sample4": [28, 24],
    }, index=["gene1", "gene2"])
    results_df = run_differential_expression(df, ["sample1", "sample2"], ["sample3", "sample4"])
    fig = plot_volcano(results_df)
    assert fig is not None


def test_plot_volcano_handles_zero_adjusted_pvalue():
    # An adjusted p-value of exactly 0 is a real possibility now (see the
    # constant-expression special case in compute_pvalues) and used to
    # send -log10 to -inf. This should plot cleanly instead of warning or
    # producing an unbounded axis.
    results_df = pd.DataFrame({
        "log_fold_change": [3.0, 0.2],
        "adjusted_p_value": [0.0, 0.6],
        "significant": [True, False],
    }, index=["gene1", "gene2"])

    with np.errstate(divide="raise"):
        fig = plot_volcano(results_df)

    assert fig is not None
    assert np.isfinite(fig.axes[0].collections[1].get_offsets()[:, 1]).all()


def test_plot_pca_returns_figure():
    df = pd.DataFrame({
        "sample1": [10, 5, 2, 8], "sample2": [12, 6, 3, 9],
        "sample3": [30, 25, 20, 28], "sample4": [28, 24, 22, 26],
    }, index=["gene1", "gene2", "gene3", "gene4"])
    labels = pd.Series(["control", "control", "treatment", "treatment"], index=df.columns)

    fig = plot_pca(df, labels)
    assert fig is not None


def test_plot_pca_too_few_samples_raises():
    df = pd.DataFrame({"sample1": [10, 5, 2]}, index=["gene1", "gene2", "gene3"])
    labels = pd.Series(["control"], index=df.columns)

    with pytest.raises(ValueError, match=r"got 1 sample"):
        plot_pca(df, labels)


def test_plot_pca_too_few_genes_raises():
    df = pd.DataFrame({"sample1": [10], "sample2": [20]}, index=["gene1"])
    labels = pd.Series(["control", "treatment"], index=df.columns)

    with pytest.raises(ValueError, match=r"got 2 sample\(s\) and 1 gene"):
        plot_pca(df, labels)
