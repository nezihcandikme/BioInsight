import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
def plot_volcano(results_df: pd.DataFrame) -> plt.Figure:
    """
    Create a volcano plot from the results of differential expression analysis.

    Parameters
    ----------
    results_df : pd.DataFrame
        A DataFrame containing the results of differential expression analysis with columns:
        - 'log_fold_change': Log fold change values for each gene.
        - 'adjusted_p_value': Adjusted p-values for each gene.
        - 'significant': Boolean indicating whether the gene is significantly differentially expressed.

    Returns
    -------
    plt.Figure
        A matplotlib Figure object containing the volcano plot.
    """
    
    # Create a new figure
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.axhline(y=-np.log10(0.05), color="black", linestyle="--", label="Significance threshold (0.05)")
    ax.set_title("Volcano Plot")

    # Add horizontal line for significance threshold (e.g., adjusted p-value < 0.05)
    significant = results_df["significant"]

    ax.scatter(
        results_df.loc[~significant, "log_fold_change"],
        -np.log10(results_df.loc[~significant, "adjusted_p_value"]),
        alpha=0.5, color="gray", label="Not significant",
    )
    ax.scatter(
        results_df.loc[significant, "log_fold_change"],
        -np.log10(results_df.loc[significant, "adjusted_p_value"]),
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
    Create a PCA plot from the given DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing the data to be used for PCA.
    labels : pd.Series
        A Series containing the labels for each sample.

    Returns
    -------
    plt.Figure
        A matplotlib Figure object containing the PCA plot.
    """
    
    # Perform PCA
    pca = PCA(n_components=2)
    coords = pca.fit_transform(df.T)
    variance_ratio = pca.explained_variance_ratio_

    # Create a new figure
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