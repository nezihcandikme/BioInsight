"""Moderated (empirical Bayes) t-test differential expression:
compute_moderated_pvalues and run_differential_expression's
method="moderated" path. See test_differential_expression_welch.py for the
plain Welch's t-test and test_differential_expression_covariates.py for the
covariate-adjusted linear model.
"""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from deconcord.differential_expression.methods import (
    compute_moderated_pvalues,
    run_differential_expression,
)


def _background_and_noisy_gene_dataset():
    # 24 "background" genes: tight, consistent within-group variance, no
    # real difference between groups -- these are what the prior gets fit
    # from. Values are deterministic (no randomness) so the test is stable.
    background = {}
    for i in range(24):
        base = 10 + i * 0.3
        background[f"background_{i}"] = [
            base + 0.05, base - 0.03, base + 0.02,   # group 1: tight
            base + 0.04, base - 0.02, base + 0.03,   # group 2: tight, same center
        ]
    # One gene with a real mean difference but one wildly noisy replicate
    # in group 2 -- inflates its own pooled variance enough that an
    # unmoderated pooled t-test misses the real effect.
    noisy_gene = [10.0, 10.0, 10.0, 12.0, 12.0, 50.0]

    data = {**background, "noisy_gene": noisy_gene}
    df = pd.DataFrame(data).T
    df.columns = ["s1", "s2", "s3", "s4", "s5", "s6"]
    return df, ["s1", "s2", "s3"], ["s4", "s5", "s6"]


def test_compute_moderated_pvalues_returns_series_indexed_like_input():
    df, group_1, group_2 = _background_and_noisy_gene_dataset()
    pvalues = compute_moderated_pvalues(df, group_1, group_2)

    assert list(pvalues.index) == list(df.index)
    assert not pvalues.isna().any()
    assert ((pvalues >= 0) & (pvalues <= 1)).all()


def test_compute_moderated_pvalues_clear_difference_is_significant():
    df = pd.DataFrame({
        "sample1": [1, 1, 1, 1, 1],
        "sample2": [1.1, 0.9, 1.0, 1.05, 0.95],
        "sample3": [1000, 1000, 1000, 1000, 1000],
        "sample4": [999, 1001, 998, 1002, 1000],
    }, index=["gene1", "gene2", "gene3", "gene4", "gene5"])

    pvalues = compute_moderated_pvalues(df, ["sample1", "sample2"], ["sample3", "sample4"])

    assert all(pvalues < 0.05)


def test_compute_moderated_pvalues_shrinkage_beats_unmoderated_for_noisy_gene():
    # The core claim of the moderated method: a gene with a real effect but
    # an inflated own-variance (from one noisy replicate) should score
    # better under shrinkage -- borrowing the tight variance of the many
    # background genes -- than it does under a plain pooled t-test using
    # only its own (noisy) variance estimate.
    df, group_1, group_2 = _background_and_noisy_gene_dataset()

    moderated = compute_moderated_pvalues(df, group_1, group_2)

    noisy = df.loc["noisy_gene"]
    n1, n2 = len(group_1), len(group_2)
    var_1 = noisy[group_1].var(ddof=1)
    var_2 = noisy[group_2].var(ddof=1)
    df_resid = n1 + n2 - 2
    pooled_var = ((n1 - 1) * var_1 + (n2 - 1) * var_2) / df_resid
    mean_diff = noisy[group_1].mean() - noisy[group_2].mean()
    t_stat = mean_diff / np.sqrt(pooled_var * (1 / n1 + 1 / n2))
    unmoderated_p = 2 * stats.t.sf(abs(t_stat), df=df_resid)

    assert moderated["noisy_gene"] < unmoderated_p


def test_compute_moderated_pvalues_constant_expression_no_nan():
    # Unlike compute_pvalues, this needs no ad hoc zero-variance special
    # case -- shrinkage toward a nonzero prior variance handles it for free.
    df, group_1, group_2 = _background_and_noisy_gene_dataset()
    df.loc["constant_gene"] = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]

    pvalues = compute_moderated_pvalues(df, group_1, group_2)

    assert not pvalues.isna().any()
    assert pvalues["constant_gene"] == pytest.approx(1.0)


def test_compute_moderated_pvalues_insufficient_samples_raises():
    df, group_1, _ = _background_and_noisy_gene_dataset()

    with pytest.raises(ValueError, match="at least 2 samples"):
        compute_moderated_pvalues(df, group_1, ["s4"])


def test_compute_moderated_pvalues_too_few_genes_raises():
    # Single gene: nothing to borrow variance information from.
    df = pd.DataFrame({
        "s1": [1.0], "s2": [1.1], "s3": [3.0], "s4": [3.1],
    }, index=["only_gene"])

    with pytest.raises(ValueError, match="at least 2 genes"):
        compute_moderated_pvalues(df, ["s1", "s2"], ["s3", "s4"])


def test_run_differential_expression_method_moderated():
    df, group_1, group_2 = _background_and_noisy_gene_dataset()

    results_df = run_differential_expression(df, group_1, group_2, method="moderated")

    assert "log_fold_change" in results_df.columns
    assert "p_value" in results_df.columns
    assert list(results_df.index) == list(df.index)


def test_run_differential_expression_invalid_method_raises():
    df = pd.DataFrame({
        "sample1": [1, 2], "sample2": [1, 2], "sample3": [3, 4], "sample4": [3, 4],
    }, index=["gene1", "gene2"])

    with pytest.raises(ValueError, match="method"):
        run_differential_expression(df, ["sample1", "sample2"], ["sample3", "sample4"], method="deseq2")

