import pandas as pd
import numpy as np

from bioinsight.differential_expression.methods import compute_log_fold_change

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