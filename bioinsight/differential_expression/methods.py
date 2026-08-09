import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

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
    if len(group_1) < 2 or len(group_2) < 2:
        raise ValueError("Each group must have at least 2 samples to perform a t-test.")
    def test_one_gene(row):
        group_1_values = row[group_1]
        group_2_values = row[group_2]
        t_stat, p_value = stats.ttest_ind(group_1_values, group_2_values, equal_var=False)
        return p_value

    return df.apply(test_one_gene, axis=1)

def compute_adjusted_pvalues(pvalues: pd.Series, method: str = 'fdr_bh') -> pd.Series:
    """
    Adjust p-values for multiple testing using the specified method.

    Parameters
    ----------
    pvalues : pd.Series
        A Series containing the p-values to be adjusted.
    method : str, optional
        The method to use for adjusting p-values. Default is 'fdr_bh' (Benjamini/Hochberg).

    Returns
    -------
    pd.Series
        A Series containing the adjusted p-values.
    """
    adjusted_pvalues = multipletests(pvalues, method=method)[1]
    return pd.Series(adjusted_pvalues, index=pvalues.index)
def run_differential_expression(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> pd.DataFrame:
    """
    Run differential expression analysis on the given DataFrame.

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
    pd.DataFrame
        A DataFrame containing the log fold change, p-values, adjusted p-values, and significance for each gene.
    """
    log_fold_change = compute_log_fold_change(df, group_1, group_2)
    pvalues = compute_pvalues(df, group_1, group_2)
    adjusted_pvalues = compute_adjusted_pvalues(pvalues)

    results_df = pd.DataFrame({
        'log_fold_change': log_fold_change,
        'p_value': pvalues,
        'adjusted_p_value': adjusted_pvalues
    })
    results_df["significant"] = (results_df["adjusted_p_value"] < 0.05) & (results_df["log_fold_change"].abs() > 1)

    return results_df
