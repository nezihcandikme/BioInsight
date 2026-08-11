"""
Pre-DE gene filtering.

Testing a gene that was barely detected in anyone is close to pure noise,
and it isn't free: every gene tested adds to the multiple-testing burden
in ``compute_adjusted_pvalues``, making the Benjamini-Hochberg correction
stricter for every *other* gene too. Dropping genes that never reached a
meaningful count anywhere, before differential expression runs, is
standard practice (edgeR's ``filterByExpr`` and DESeq2's own
recommendations both do a version of this) — not an optional nicety.

This module intentionally does a simpler version of that idea: a fixed
count threshold met in a minimum number of samples, rather than
edgeR's design-aware ``filterByExpr`` logic. Good enough to stop testing
obviously-dead genes; not a reimplementation of that algorithm.
"""

import pandas as pd


def filter_low_count_genes(df: pd.DataFrame, min_count: int = 10, min_samples: int = 2) -> pd.DataFrame:
    """
    Drop genes that never reach ``min_count`` raw reads in at least
    ``min_samples`` samples.

    Run this on raw counts, before normalization — CPM/log2-transformed
    values are on the wrong scale for a raw-count threshold to mean
    anything.

    Parameters
    ----------
    df : pd.DataFrame
        Raw count matrix, genes as rows, samples as columns.
    min_count : int, optional
        Minimum raw count for a sample to "count" for a gene. Default 10.
    min_samples : int, optional
        Minimum number of samples that need to clear ``min_count`` for the
        gene to be kept. Default 2 — a gene expressed in just one sample
        out of a whole experiment is usually noise, not biology, but
        picking the right number here depends on your design (e.g. use
        your smallest group size, not a fixed constant, if group sizes
        vary a lot).

    Returns
    -------
    pd.DataFrame
        Subset of ``df`` (same columns, fewer rows).

    Raises
    ------
    ValueError
        If ``min_count`` is negative, ``min_samples`` is less than 1,
        ``min_samples`` exceeds the number of samples in ``df`` (which
        would filter out every gene by construction), or the filter
        would remove every gene given the actual data — each of these
        is more useful caught here than as a confusing empty-matrix
        failure two steps later in normalization or DE.
    """
    if min_count < 0:
        raise ValueError(f"min_count must be >= 0; got {min_count}.")
    if min_samples < 1:
        raise ValueError(f"min_samples must be >= 1; got {min_samples}.")
    if min_samples > df.shape[1]:
        raise ValueError(
            f"min_samples ({min_samples}) is larger than the number of "
            f"samples in the matrix ({df.shape[1]}) — no gene could ever "
            "pass this filter."
        )

    keep = (df >= min_count).sum(axis=1) >= min_samples

    if not keep.any():
        raise ValueError(
            f"No genes pass min_count={min_count}, min_samples={min_samples} "
            f"on this data (started with {df.shape[0]} genes). Lower one or "
            "both thresholds, or double check the matrix isn't already "
            "normalized (this filter expects raw counts)."
        )

    return df.loc[keep]
