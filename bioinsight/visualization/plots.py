import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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