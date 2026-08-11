"""
Differential expression analysis (basic, exploratory method).

IMPORTANT: The functions in this module implement a simple mean-difference
t-test approach (a "poor man's DESeq2"), NOT a count-based negative-binomial
model like DESeq2 or edgeR. They do not model biological dispersion or the
mean-variance relationship of raw RNA-seq counts, and they do not perform
independent normalization internally.

Practical implications:
- ``df`` is expected to already be normalized and log-transformed (e.g. via
  ``bioinsight.normalization.methods.compute_cpm`` followed by
  ``log2_transform``) before being passed in here. If you pass in raw counts,
  ``log_fold_change`` will just be a difference of raw count means, not a
  true log fold change.
- This method is best suited for quick, exploratory looks at a dataset. For
  publication-quality or clinically meaningful differential expression
  results, use an established count-based tool such as DESeq2 or edgeR,
  which BioInsight does not currently replace.
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def compute_log_fold_change(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> pd.Series:
    """
    Compute the difference in mean expression between two groups of samples.

    This is only a true "log fold change" if ``df`` already contains
    log-transformed values (e.g. log2-CPM). This is a simple exploratory
    metric, not a model-based estimate like the shrunken log fold changes
    produced by DESeq2 or edgeR.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing gene expression data with genes as rows and
        samples as columns. Should already be normalized and log-transformed.
    group_1 : list[str]
        A list of sample names for the first group.
    group_2 : list[str]
        A list of sample names for the second group.

    Returns
    -------
    pd.Series
        A Series containing the mean-difference "log fold change" values for
        each gene.
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
        group_1_values = row[group_1].astype(float)
        group_2_values = row[group_2].astype(float)

        # Welch's t-test is undefined (returns NaN) when both groups have
        # zero variance, which happens for genes with constant expression
        # (e.g. all-zero rows). Handle that case explicitly instead of
        # silently propagating NaN p-values downstream.
        if group_1_values.var(ddof=1) == 0 and group_2_values.var(ddof=1) == 0:
            means_equal = np.isclose(group_1_values.mean(), group_2_values.mean())
            return 1.0 if means_equal else 0.0

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
    Run a basic, exploratory differential expression analysis on the given DataFrame.

    NOTE: This uses a mean-difference metric plus a Welch's t-test on
    per-gene values, not a count-based negative-binomial model. It is
    intended for quick exploratory analysis, not as a replacement for
    DESeq2/edgeR. ``df`` should already be normalized and log-transformed
    (e.g. log2-CPM) — see ``bioinsight.normalization.methods``.

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
