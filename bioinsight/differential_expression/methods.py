import pandas as pd
import numpy as np
from scipy import stats


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

def compute_pvalues(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> pd.Series:
    """
    Compute the p-values for each gene between two groups of samples using a t-test.

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
        A Series containing the p-values for each gene.
    """
    def test_one_gene(row):
        group_1_values = row[group_1]
        group_2_values = row[group_2]
        t_stat, p_value = stats.ttest_ind(group_1_values, group_2_values, equal_var=False)
        return p_value
    df.apply(test_one_gene, axis=1)

    return df.apply(test_one_gene, axis=1)
