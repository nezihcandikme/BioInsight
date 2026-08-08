import pandas as pd
import numpy as np

def compute_log_fold_change(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> pd.Series:
    """
    Compute the log fold change between two groups of samples.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing gene expression data with genes as rows and samples as columns.
    group_1 : list[str]
        A list of sample names for the first group.
    group_2 : list[str]
        A list of sample names for the second group.

    Returns
    -------
    pd.Series
        A Series containing the log fold change values for each gene.
    """
    # Calculate the mean expression for each group
    mean_group_1 = df[group_1].mean(axis=1)
    mean_group_2 = df[group_2].mean(axis=1)

    # Compute log fold change
    log_fold_change = mean_group_1 - mean_group_2

    return log_fold_change
    