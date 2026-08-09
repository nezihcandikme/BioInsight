import pandas as pd

from bioinsight.differential_expression.methods import run_differential_expression
from bioinsight.visualization.plots import plot_volcano, plot_pca


def test_plot_volcano_returns_figure():
    df = pd.DataFrame({
        "sample1": [10, 5], "sample2": [20, 15],
        "sample3": [30, 25], "sample4": [28, 24],
    }, index=["gene1", "gene2"])
    results_df = run_differential_expression(df, ["sample1", "sample2"], ["sample3", "sample4"])
    fig = plot_volcano(results_df)
    assert fig is not None

def test_plot_pca_returns_figure():
    df = pd.DataFrame({
        "sample1": [10, 5, 2, 8], "sample2": [12, 6, 3, 9],
        "sample3": [30, 25, 20, 28], "sample4": [28, 24, 22, 26],
    }, index=["gene1", "gene2", "gene3", "gene4"])
    labels = pd.Series(["control", "control", "treatment", "treatment"], index=df.columns)

    fig = plot_pca(df, labels)
    assert fig is not None


