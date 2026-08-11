import numpy as np
import pandas as pd

from bioinsight.core import library_size


def compute_cpm(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize raw counts to Counts Per Million (CPM): each sample's counts
    divided by that sample's total count, times one million.

    CPM corrects for sequencing depth (library size) so samples become
    comparable to each other. It does *not* correct for gene length, RNA
    composition effects, or other technical biases — for those, TPM or a
    method like DESeq2's median-of-ratios normalization would be needed
    instead.

    Parameters
    ----------
    df : pd.DataFrame
        Raw count matrix, genes as rows, samples as columns.

    Returns
    -------
    pd.DataFrame
        Same shape as ``df``, with CPM values.

    Raises
    ------
    ValueError
        If any sample has a total count of zero. CPM divides by the
        sample's total, so a zero-count sample isn't just noisy data —
        it makes the whole column mathematically undefined (0/0). That
        sample needs to be dropped or investigated before normalizing,
        not silently turned into NaN or inf.
    """
    total_counts = library_size(df)

    zero_total_samples = total_counts[total_counts == 0].index.tolist()
    if zero_total_samples:
        raise ValueError(
            f"Sample(s) {zero_total_samples} have zero total counts; "
            "CPM is undefined for these samples. Remove them or check "
            "the input count matrix before normalizing."
        )

    return df.div(total_counts, axis=1) * 1e6


def log2_transform(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply a log2(x + 1) transformation.

    The +1 pseudocount avoids log2(0), which is undefined, at the cost of
    a small compression of low counts. This is the conventional tradeoff
    for RNA-seq-style data and matches what most exploratory tools do.

    Parameters
    ----------
    df : pd.DataFrame
        Non-negative expression values (typically CPM output from
        ``compute_cpm``, not raw counts).

    Returns
    -------
    pd.DataFrame
        Same shape as ``df``, log2(x + 1)-transformed.

    Raises
    ------
    ValueError
        If ``df`` contains negative values.
    """
    if (df < 0).any().any():
        raise ValueError("Input DataFrame contains negative values. Log2 transformation is not defined for negative numbers.")
    return np.log2(df + 1)
