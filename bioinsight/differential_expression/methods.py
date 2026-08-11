"""
Differential expression analysis — basic, exploratory method.

The functions here implement a mean-difference metric plus a per-gene
Welch's t-test. This is NOT a count-aware negative-binomial model like
DESeq2 or edgeR: it doesn't model biological dispersion or the
mean-variance relationship that raw RNA-seq counts actually have, and it
performs no normalization of its own.

Practical implications:
- ``df`` is expected to already be normalized and log-transformed (e.g.
  ``compute_cpm`` followed by ``log2_transform`` from
  ``bioinsight.normalization.methods``) before being passed in here. Feed
  it raw counts and "log fold change" becomes a difference of raw count
  means — a real number, but not a fold change of anything.
- Treat this module's output as a fast, exploratory first look. For
  results meant to support a real biological conclusion, use DESeq2 or
  edgeR, which BioInsight does not replace.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def _validate_groups(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> None:
    """
    Shared input checks for group_1/group_2, used by every public function
    in this module. Exists so a typo'd sample name or an accidentally
    reused sample fails loudly and specifically, instead of quietly
    producing a nonsense comparison or a confusing pandas KeyError.
    """
    if not group_1:
        raise ValueError("group_1 is empty — it needs at least one sample name.")
    if not group_2:
        raise ValueError("group_2 is empty — it needs at least one sample name.")

    for group_name, group in (("group_1", group_1), ("group_2", group_2)):
        seen = set()
        for sample in group:
            if sample in seen:
                raise ValueError(
                    f"Sample '{sample}' appears more than once in {group_name}."
                )
            seen.add(sample)

    overlap = set(group_1) & set(group_2)
    if overlap:
        raise ValueError(
            f"Sample(s) {sorted(overlap)} appear in both group_1 and group_2 "
            "— a sample can't be compared against itself."
        )

    for group_name, group in (("group_1", group_1), ("group_2", group_2)):
        for sample in group:
            if sample not in df.columns:
                raise ValueError(
                    f"Sample '{sample}' was requested in {group_name} but is "
                    "not present in the count matrix."
                )


def compute_log_fold_change(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> pd.Series:
    """
    Compute the difference in mean expression between two groups, per gene.

    This is only a genuine "log fold change" if ``df`` is already on a
    log2 scale — a difference of logs is a log ratio. Pass in raw or
    linear-scale CPM values and this becomes a difference of means, which
    is a different (and less standard) quantity. This is also a simple
    plug-in estimate, not the shrunken effect size DESeq2/edgeR report.

    Parameters
    ----------
    df : pd.DataFrame
        Log-transformed, normalized expression matrix (genes as rows,
        samples as columns).
    group_1, group_2 : list[str]
        Sample names (column names in ``df``) for each comparison group.

    Returns
    -------
    pd.Series
        Per-gene ``mean(group_1) - mean(group_2)``, indexed by gene.
    """
    _validate_groups(df, group_1, group_2)

    mean_group_1 = df[group_1].mean(axis=1)
    mean_group_2 = df[group_2].mean(axis=1)

    return mean_group_1 - mean_group_2


def compute_pvalues(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> pd.Series:
    """
    Per-gene Welch's t-test between two groups.

    Welch's t-test (unequal-variance t-test) is used instead of the
    standard Student's t-test because there's no reason to assume the two
    groups have equal variance — and RNA-seq expression variance is
    routinely condition-dependent, so assuming equal variance would be the
    wrong default.

    Parameters
    ----------
    df : pd.DataFrame
        Log-transformed, normalized expression matrix (genes as rows,
        samples as columns).
    group_1, group_2 : list[str]
        Sample names (column names in ``df``) for each comparison group.
        Each group needs at least 2 samples — a t-test needs within-group
        variance to estimate, and one sample can't provide that.

    Returns
    -------
    pd.Series
        Raw (not multiple-testing-corrected) p-value per gene.

    Raises
    ------
    ValueError
        If either group has fewer than 2 samples, or the group inputs fail
        the shared validation in ``_validate_groups`` (missing/duplicate/
        overlapping sample names).
    """
    _validate_groups(df, group_1, group_2)

    if len(group_1) < 2 or len(group_2) < 2:
        raise ValueError("Each group must have at least 2 samples to perform a t-test.")

    def test_one_gene(row):
        group_1_values = row[group_1].astype(float)
        group_2_values = row[group_2].astype(float)

        # Welch's t-test is undefined (NaN) when both groups have zero
        # variance — e.g. an all-zero or otherwise constant-expression
        # gene. Rather than let NaN silently leak into the results table,
        # resolve it directly: identical constants means "not different"
        # (p=1), differing constants means "as different as it gets" (p=0).
        if group_1_values.var(ddof=1) == 0 and group_2_values.var(ddof=1) == 0:
            means_equal = np.isclose(group_1_values.mean(), group_2_values.mean())
            return 1.0 if means_equal else 0.0

        _, p_value = stats.ttest_ind(group_1_values, group_2_values, equal_var=False)
        return p_value

    return df.apply(test_one_gene, axis=1)


def compute_adjusted_pvalues(pvalues: pd.Series, method: str = 'fdr_bh') -> pd.Series:
    """
    Correct p-values for multiple testing.

    Testing thousands of genes at once means thousands of chances for a
    p < 0.05 result by pure chance — without correction, a typical RNA-seq
    experiment would report hundreds of "significant" genes even if
    nothing in the biology changed. Benjamini-Hochberg (the default here)
    controls the expected false discovery rate rather than the much
    stricter per-test error rate, which is the standard, less
    conservative choice for exploratory genomics.

    Parameters
    ----------
    pvalues : pd.Series
        Raw per-gene p-values.
    method : str, optional
        Passed through to ``statsmodels.stats.multitest.multipletests``.
        Default ``'fdr_bh'`` (Benjamini-Hochberg).

    Returns
    -------
    pd.Series
        Adjusted p-values, same index as ``pvalues``.

    Raises
    ------
    ValueError
        If ``pvalues`` contains any NaN. ``statsmodels`` doesn't just skip
        the offending gene here — a single NaN input silently turns every
        adjusted p-value in the whole Series into NaN, which is exactly
        the kind of quiet corruption this function should never allow
        through. ``compute_pvalues`` in this module already avoids
        producing NaN, so this mainly guards against p-values coming from
        somewhere else.
    """
    # statsmodels divides by len(pvals) internally and doesn't guard
    # against an empty array itself (ZeroDivisionError), so handle the
    # "nothing to correct" case here instead of letting that surface.
    if pvalues.empty:
        return pd.Series(dtype=float, index=pvalues.index)

    if pvalues.isna().any():
        bad_genes = pvalues[pvalues.isna()].index.tolist()
        raise ValueError(
            f"pvalues contains NaN for gene(s) {bad_genes}. Multiple-testing "
            "correction can't proceed — a single NaN p-value poisons every "
            "adjusted p-value in the batch, not just its own row."
        )

    adjusted_pvalues = multipletests(pvalues, method=method)[1]
    return pd.Series(adjusted_pvalues, index=pvalues.index)


def run_differential_expression(
    df: pd.DataFrame,
    group_1: list[str],
    group_2: list[str],
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
) -> pd.DataFrame:
    """
    Run the full exploratory DE workflow: log fold change, Welch's t-test,
    and Benjamini-Hochberg correction, combined into one results table.

    Not a substitute for DESeq2/edgeR — see the module docstring. ``df``
    should already be normalized and log-transformed (e.g. log2-CPM via
    ``bioinsight.normalization.methods``).

    Parameters
    ----------
    df : pd.DataFrame
        Log-transformed, normalized expression matrix (genes as rows,
        samples as columns).
    group_1, group_2 : list[str]
        Sample names (column names in ``df``) for each comparison group.
    alpha : float, optional
        Adjusted p-value cutoff for the ``significant`` column. Default
        0.05. There's nothing statistically special about 0.05 — it's a
        convention, not a law — so this is exposed rather than buried.
    lfc_threshold : float, optional
        Minimum ``abs(log_fold_change)`` for the ``significant`` column.
        Default 1.0 (i.e. a 2x difference on a log2 scale). Only meaningful
        as a fold-change cutoff if ``df`` is actually log2-scaled — see the
        module docstring.

    Returns
    -------
    pd.DataFrame
        Indexed by gene, with columns:
        - ``log_fold_change``: mean(group_1) - mean(group_2)
        - ``p_value``: raw Welch's t-test p-value
        - ``adjusted_p_value``: Benjamini-Hochberg corrected p-value
        - ``significant``: ``adjusted_p_value < alpha`` and
          ``abs(log_fold_change) > lfc_threshold`` — a starting point for
          further inspection, not a verdict, regardless of what you set
          the thresholds to.

    Raises
    ------
    ValueError
        If ``alpha`` is not in (0, 1], or ``lfc_threshold`` is negative.
    """
    if not (0 < alpha <= 1):
        raise ValueError(f"alpha must be in (0, 1]; got {alpha}.")
    if lfc_threshold < 0:
        raise ValueError(f"lfc_threshold must be >= 0; got {lfc_threshold}.")

    log_fold_change = compute_log_fold_change(df, group_1, group_2)
    pvalues = compute_pvalues(df, group_1, group_2)
    adjusted_pvalues = compute_adjusted_pvalues(pvalues)

    results_df = pd.DataFrame({
        'log_fold_change': log_fold_change,
        'p_value': pvalues,
        'adjusted_p_value': adjusted_pvalues
    })
    results_df["significant"] = (results_df["adjusted_p_value"] < alpha) & (results_df["log_fold_change"].abs() > lfc_threshold)

    return results_df
