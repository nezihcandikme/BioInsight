import pandas as pd
import numpy as np
from scipy import stats
import pytest

from bioinsight.differential_expression.methods import compute_log_fold_change, compute_pvalues

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