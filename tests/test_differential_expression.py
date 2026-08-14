import numpy as np
import pandas as pd
import pytest
from scipy import stats

from omicforge.differential_expression.methods import (
    compute_log_fold_change,
    compute_pvalues,
    compute_moderated_pvalues,
    compute_adjusted_pvalues,
    run_differential_expression,
    run_differential_expression_with_covariates,
)


def test_compute_log_fold_change():
    df = pd.DataFrame({
        "sample1": [10, 5, 2],
        "sample2": [20, 15, 4],
        "sample3": [30, 25, 6],
    }, index=["gene1", "gene2", "gene3"])

    group_1 = ['sample1', 'sample2']
    group_2 = ['sample3']

    log_fc = compute_log_fold_change(df, group_1, group_2)

    expected_log_fc = pd.Series({
        'gene1': (10 + 20) / 2 - 30,
        'gene2': (5 + 15) / 2 - 25,
        'gene3': (2 + 4) / 2 - 6
    })

    pd.testing.assert_series_equal(log_fc, expected_log_fc)


def test_compute_pvalues_clear_difference():
    df = pd.DataFrame({
        "sample1": [1, 1, 1],
        "sample2": [2, 2, 2],
        "sample3": [1000, 1000, 1000],
        "sample4": [999, 999, 999],
    }, index=["gene1", "gene2", "gene3"])

    pvalues = compute_pvalues(df, ["sample1", "sample2"], ["sample3", "sample4"])

    assert all(pvalues < 0.05)


def test_compute_pvalues_no_difference():
    df = pd.DataFrame({
        "sample1": [9, 19, 4],
        "sample2": [11, 21, 6],
        "sample3": [8, 18, 3],
        "sample4": [12, 22, 7],
    }, index=["gene1", "gene2", "gene3"])

    pvalues = compute_pvalues(df, ["sample1", "sample2"], ["sample3", "sample4"])

    assert all(pvalues.apply(lambda p: p == pytest.approx(1.0, abs=0.01)))


def test_compute_adjusted_pvalues_single():
    pvalues = pd.Series({"gene1": 0.03})
    adjusted = compute_adjusted_pvalues(pvalues)
    assert adjusted["gene1"] == pytest.approx(0.03)


def test_compute_adjusted_pvalues_never_smaller():
    pvalues = pd.Series({"gene1": 0.01, "gene2": 0.2, "gene3": 0.04, "gene4": 0.5})
    adjusted = compute_adjusted_pvalues(pvalues)
    assert all(adjusted >= pvalues)


def test_compute_adjusted_pvalues_equal_spacing():
    pvalues = pd.Series({
        "gene1": 0.01, "gene2": 0.02, "gene3": 0.03, "gene4": 0.04, "gene5": 0.05
    })
    adjusted = compute_adjusted_pvalues(pvalues)
    assert all(adjusted.apply(lambda p: p == pytest.approx(0.05, abs=1e-6)))


def test_compute_adjusted_pvalues_nan_input_raises():
    # A single NaN silently poisons every adjusted p-value via statsmodels
    # — that should never happen quietly, so this must raise instead.
    import numpy as np

    pvalues = pd.Series({"gene1": 0.01, "gene2": np.nan, "gene3": 0.03})

    with pytest.raises(ValueError, match="gene2"):
        compute_adjusted_pvalues(pvalues)


def test_run_differential_expression():
    df = pd.DataFrame({
        "sample1": [10, 5, 2],
        "sample2": [20, 15, 4],
        "sample3": [30, 25, 6],
        "sample4": [28, 24, 5],
    }, index=["gene1", "gene2", "gene3"])

    results_df = run_differential_expression(df, ["sample1", "sample2"], ["sample3", "sample4"])

    assert "log_fold_change" in results_df.columns
    assert "p_value" in results_df.columns
    assert "adjusted_p_value" in results_df.columns
    assert "significant" in results_df.columns


def test_run_differential_expression_custom_thresholds_change_significance():
    df = pd.DataFrame({
        "sample1": [1, 1, 1],
        "sample2": [2, 2, 2],
        "sample3": [3, 3, 3],
        "sample4": [4, 4, 4],
    }, index=["gene1", "gene2", "gene3"])

    default_results = run_differential_expression(df, ["sample1", "sample2"], ["sample3", "sample4"])
    lenient_results = run_differential_expression(
        df, ["sample1", "sample2"], ["sample3", "sample4"], alpha=1.0, lfc_threshold=0.0
    )

    # Same underlying numbers either way — only which rows count as
    # "significant" should move when the thresholds move.
    pd.testing.assert_series_equal(default_results["p_value"], lenient_results["p_value"])
    assert lenient_results["significant"].sum() >= default_results["significant"].sum()


def test_run_differential_expression_invalid_alpha_raises():
    df = pd.DataFrame({
        "sample1": [1, 2], "sample2": [1, 2], "sample3": [3, 4], "sample4": [3, 4],
    }, index=["gene1", "gene2"])

    with pytest.raises(ValueError, match="alpha"):
        run_differential_expression(df, ["sample1", "sample2"], ["sample3", "sample4"], alpha=0)

    with pytest.raises(ValueError, match="alpha"):
        run_differential_expression(df, ["sample1", "sample2"], ["sample3", "sample4"], alpha=1.5)


def test_run_differential_expression_negative_lfc_threshold_raises():
    df = pd.DataFrame({
        "sample1": [1, 2], "sample2": [1, 2], "sample3": [3, 4], "sample4": [3, 4],
    }, index=["gene1", "gene2"])

    with pytest.raises(ValueError, match="lfc_threshold"):
        run_differential_expression(df, ["sample1", "sample2"], ["sample3", "sample4"], lfc_threshold=-1)


def test_compute_pvalues_constant_expression_equal_means():
    df = pd.DataFrame({
        "sample1": [5, 5],
        "sample2": [5, 5],
        "sample3": [5, 5],
        "sample4": [5, 5],
    }, index=["gene1", "gene2"])

    pvalues = compute_pvalues(df, ["sample1", "sample2"], ["sample3", "sample4"])

    assert not pvalues.isna().any()
    assert all(pvalues.apply(lambda p: p == pytest.approx(1.0)))


def test_compute_pvalues_constant_expression_different_means():
    df = pd.DataFrame({
        "sample1": [1, 1],
        "sample2": [1, 1],
        "sample3": [9, 9],
        "sample4": [9, 9],
    }, index=["gene1", "gene2"])

    pvalues = compute_pvalues(df, ["sample1", "sample2"], ["sample3", "sample4"])

    assert not pvalues.isna().any()
    assert all(pvalues.apply(lambda p: p == pytest.approx(0.0)))


def test_run_differential_expression_no_nan_leaks_through_for_constant_genes():
    # Regression test tying the constant-expression fix in compute_pvalues
    # to the NaN guard in compute_adjusted_pvalues: a dataset mixing a
    # constant gene with normal ones should produce a fully usable table.
    df = pd.DataFrame({
        "sample1": [5, 1, 3],
        "sample2": [5, 2, 4],
        "sample3": [5, 1000, 15],
        "sample4": [5, 999, 13],
    }, index=["constant_gene", "clear_gene", "ordinary_gene"])

    results_df = run_differential_expression(df, ["sample1", "sample2"], ["sample3", "sample4"])

    assert not results_df["p_value"].isna().any()
    assert not results_df["adjusted_p_value"].isna().any()


def test_compute_pvalues_insufficient_samples():
    df = pd.DataFrame({
        "sample1": [10, 5, 2],
        "sample2": [20, 15, 4],
        "sample3": [30, 25, 6],
    }, index=["gene1", "gene2", "gene3"])

    with pytest.raises(ValueError):
        compute_pvalues(df, ["sample1", "sample2"], ["sample3"])


def test_compute_pvalues_missing_sample_names_the_sample():
    df = pd.DataFrame({
        "sample1": [10, 5], "sample2": [20, 15],
    }, index=["gene1", "gene2"])

    with pytest.raises(ValueError, match="control_1"):
        compute_pvalues(df, ["sample1", "control_1"], ["sample2"])


def test_compute_log_fold_change_duplicate_sample_in_same_group_raises():
    df = pd.DataFrame({
        "sample1": [10, 5], "sample2": [20, 15],
    }, index=["gene1", "gene2"])

    with pytest.raises(ValueError, match="sample1"):
        compute_log_fold_change(df, ["sample1", "sample1"], ["sample2"])


def test_compute_log_fold_change_sample_in_both_groups_raises():
    df = pd.DataFrame({
        "sample1": [10, 5], "sample2": [20, 15],
    }, index=["gene1", "gene2"])

    with pytest.raises(ValueError, match="sample1"):
        compute_log_fold_change(df, ["sample1"], ["sample1", "sample2"])


def test_compute_log_fold_change_empty_group_raises():
    df = pd.DataFrame({
        "sample1": [10, 5], "sample2": [20, 15],
    }, index=["gene1", "gene2"])

    with pytest.raises(ValueError, match="group_2"):
        compute_log_fold_change(df, ["sample1"], [])


# ---------------------------------------------------------------------
# Moderated (empirical Bayes) t-test
# ---------------------------------------------------------------------

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


# ---------------------------------------------------------------------
# Covariate-adjusted linear model DE test
# ---------------------------------------------------------------------

def _confounded_batch_dataset():
    # group_1 = s1,s2,s3 / group_2 = s4,s5,s6, but batch membership is
    # unevenly split across the two condition groups (2 group_1 samples
    # and 1 group_2 sample are batch_B) -- a real, if imperfect,
    # confound. True condition effect is +2.0; true batch_B effect is
    # +5.0 on top of a baseline of 10. Small deterministic perturbations
    # are added so the model doesn't fit exactly (residual variance would
    # otherwise be zero).
    group_1 = ["s1", "s2", "s3"]
    group_2 = ["s4", "s5", "s6"]
    batch = {"s1": "batch_A", "s2": "batch_A", "s3": "batch_B",
             "s4": "batch_B", "s5": "batch_B", "s6": "batch_A"}
    noiseless = {"s1": 12.0, "s2": 12.0, "s3": 17.0, "s4": 15.0, "s5": 15.0, "s6": 10.0}
    perturbation = {"s1": 0.05, "s2": -0.05, "s3": 0.02, "s4": -0.02, "s5": 0.03, "s6": -0.03}
    values = {s: noiseless[s] + perturbation[s] for s in noiseless}

    df = pd.DataFrame({"gene1": values}).T
    metadata = pd.DataFrame({"batch": batch})
    return df, group_1, group_2, metadata


def test_run_differential_expression_with_covariates_matches_manual_ols():
    df, group_1, group_2, metadata = _confounded_batch_dataset()
    samples = group_1 + group_2

    results = run_differential_expression_with_covariates(
        df, group_1, group_2, metadata, ["batch"], moderated=False,
    )

    # Rebuild the same design matrix independently (intercept, condition,
    # batch_B dummy -- batch_A is the alphabetically-first, dropped level)
    # and solve it directly, to check the function's result against a
    # from-scratch computation rather than another layer of the same code.
    condition = np.array([1.0 if s in group_1 else 0.0 for s in samples])
    batch_b = np.array([1.0 if metadata.loc[s, "batch"] == "batch_B" else 0.0 for s in samples])
    X = np.column_stack([np.ones(6), condition, batch_b])
    y = df.loc["gene1", samples].to_numpy(dtype=float)

    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    expected_lfc = beta[1]

    assert results.loc["gene1", "log_fold_change"] == pytest.approx(expected_lfc)
    assert results.loc["gene1", "log_fold_change"] == pytest.approx(2.0, abs=0.2)


def test_run_differential_expression_with_covariates_recovers_effect_naive_test_confounds():
    # The point of the feature: the plain two-group comparison conflates
    # the batch effect into the condition estimate (batch is unevenly
    # split across the groups), while the covariate-adjusted model
    # isolates the true ~2.0 condition effect.
    df, group_1, group_2, metadata = _confounded_batch_dataset()

    naive_lfc = compute_log_fold_change(df, group_1, group_2)["gene1"]
    adjusted = run_differential_expression_with_covariates(
        df, group_1, group_2, metadata, ["batch"], moderated=False,
    )
    adjusted_lfc = adjusted.loc["gene1", "log_fold_change"]

    assert abs(adjusted_lfc - 2.0) < abs(naive_lfc - 2.0)


def test_run_differential_expression_with_covariates_numeric_covariate():
    df, group_1, group_2, metadata = _confounded_batch_dataset()
    metadata["numeric_batch"] = metadata["batch"].map({"batch_A": 0.0, "batch_B": 1.0})

    results = run_differential_expression_with_covariates(
        df, group_1, group_2, metadata, ["numeric_batch"], moderated=False,
    )

    assert results.loc["gene1", "log_fold_change"] == pytest.approx(2.0, abs=0.2)


def test_run_differential_expression_with_covariates_returns_expected_columns():
    df, group_1, group_2, metadata = _confounded_batch_dataset()

    results = run_differential_expression_with_covariates(
        df, group_1, group_2, metadata, ["batch"], moderated=False,
    )

    assert list(results.columns) == ["log_fold_change", "p_value", "adjusted_p_value", "significant"]
    assert list(results.index) == list(df.index)


def test_run_differential_expression_with_covariates_empty_covariate_cols_raises():
    df, group_1, group_2, metadata = _confounded_batch_dataset()

    with pytest.raises(ValueError, match="covariate_cols"):
        run_differential_expression_with_covariates(df, group_1, group_2, metadata, [])


def test_run_differential_expression_with_covariates_missing_sample_raises():
    df, group_1, group_2, metadata = _confounded_batch_dataset()
    metadata = metadata.drop("s3")

    with pytest.raises(ValueError, match="s3"):
        run_differential_expression_with_covariates(df, group_1, group_2, metadata, ["batch"])


def test_run_differential_expression_with_covariates_missing_column_raises():
    df, group_1, group_2, metadata = _confounded_batch_dataset()

    with pytest.raises(ValueError, match="not_a_real_column"):
        run_differential_expression_with_covariates(
            df, group_1, group_2, metadata, ["not_a_real_column"],
        )


def test_run_differential_expression_with_covariates_nan_covariate_raises():
    df, group_1, group_2, metadata = _confounded_batch_dataset()
    metadata.loc["s2", "batch"] = None

    with pytest.raises(ValueError, match="s2"):
        run_differential_expression_with_covariates(df, group_1, group_2, metadata, ["batch"])


def test_run_differential_expression_with_covariates_constant_covariate_raises():
    df, group_1, group_2, metadata = _confounded_batch_dataset()
    metadata["batch"] = "same_batch_for_everyone"

    with pytest.raises(ValueError, match="only one distinct value"):
        run_differential_expression_with_covariates(df, group_1, group_2, metadata, ["batch"], moderated=False)


def test_run_differential_expression_with_covariates_constant_numeric_covariate_raises():
    df, group_1, group_2, metadata = _confounded_batch_dataset()
    metadata["numeric_batch"] = 3.0

    with pytest.raises(ValueError, match="rank-deficient"):
        run_differential_expression_with_covariates(
            df, group_1, group_2, metadata, ["numeric_batch"], moderated=False,
        )


def test_run_differential_expression_with_covariates_collinear_with_condition_raises():
    df, group_1, group_2, metadata = _confounded_batch_dataset()
    # A covariate that's just a relabeling of group membership is perfectly
    # collinear with the condition term -- their individual effects can't
    # be separated.
    metadata["batch"] = ["batch_A" if s in group_1 else "batch_B" for s in metadata.index]

    with pytest.raises(ValueError, match="rank-deficient"):
        run_differential_expression_with_covariates(df, group_1, group_2, metadata, ["batch"])


def test_run_differential_expression_with_covariates_moderated_needs_two_genes():
    df, group_1, group_2, metadata = _confounded_batch_dataset()

    with pytest.raises(ValueError, match="at least 2 genes"):
        run_differential_expression_with_covariates(df, group_1, group_2, metadata, ["batch"], moderated=True)


def test_run_differential_expression_with_covariates_moderated_multi_gene():
    df, group_1, group_2, metadata = _confounded_batch_dataset()
    # Add background genes with tight, consistent variance so the
    # empirical Bayes prior has something to fit.
    for i in range(6):
        base = 10 + i * 0.4
        df.loc[f"background_{i}"] = [
            base + 0.03, base - 0.02, base + 0.04, base - 0.01, base + 0.02, base - 0.03,
        ]

    results = run_differential_expression_with_covariates(
        df, group_1, group_2, metadata, ["batch"], moderated=True,
    )

    assert not results["p_value"].isna().any()
    assert ((results["p_value"] >= 0) & (results["p_value"] <= 1)).all()


def test_run_differential_expression_with_covariates_invalid_alpha_raises():
    df, group_1, group_2, metadata = _confounded_batch_dataset()

    with pytest.raises(ValueError, match="alpha"):
        run_differential_expression_with_covariates(df, group_1, group_2, metadata, ["batch"], alpha=0)


def test_run_differential_expression_with_covariates_negative_lfc_threshold_raises():
    df, group_1, group_2, metadata = _confounded_batch_dataset()

    with pytest.raises(ValueError, match="lfc_threshold"):
        run_differential_expression_with_covariates(
            df, group_1, group_2, metadata, ["batch"], lfc_threshold=-1,
        )


def test_run_differential_expression_with_covariates_no_covariates_message():
    df, group_1, group_2, metadata = _confounded_batch_dataset()

    with pytest.raises(ValueError, match="run_differential_expression instead"):
        run_differential_expression_with_covariates(df, group_1, group_2, metadata, [])
