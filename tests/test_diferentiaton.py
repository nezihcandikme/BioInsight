import pandas as pd
import numpy as np
from scipy import stats
import pytest
from statsmodels.stats.multitest import multipletests

from bioinsight.differential_expression.methods import compute_log_fold_change, compute_pvalues, compute_adjusted_pvalues, run_differential_expression

def test_compute_log_fold_change():
    # Create a sample DataFrame
    df = pd.DataFrame({
    "sample1": [10, 5, 2],
    "sample2": [20, 15, 4],
    "sample3": [30, 25, 6],
    }, index=["gene1", "gene2", "gene3"])

    # Define groups
    group_1 = ['sample1', 'sample2']
    group_2 = ['sample3']

    # Compute log fold change
    log_fc = compute_log_fold_change(df, group_1, group_2)

    # Expected log fold change values
    expected_log_fc = pd.Series({
        'gene1': (10 + 20) / 2 - 30,
        'gene2': (5 + 15) / 2 - 25,
        'gene3': (2 + 4) / 2 - 6
    })

    # Assert that the computed log fold change matches the expected values
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


def test_compute_pvalues_insufficient_samples():
    df = pd.DataFrame({
        "sample1": [10, 5, 2],
        "sample2": [20, 15, 4],
        "sample3": [30, 25, 6],
    }, index=["gene1", "gene2", "gene3"])

    with pytest.raises(ValueError):
        compute_pvalues(df, ["sample1", "sample2"], ["sample3"])
