import pandas as pd
import pytest
from bioinsight.core import library_size, genes_detected
from bioinsight.normalization.methods import compute_cpm

def test_compute_cpm():
    # Create a sample DataFrame with raw counts
    data = {
        'Sample1': [100, 200, 300],
        'Sample2': [400, 500, 600],
        'Sample3': [700, 800, 900]
    }
    df = pd.DataFrame(data)

    # Compute CPM using the function
    cpm_df = compute_cpm(df)

    # Calculate expected CPM values manually
    total_counts = library_size(df)
    assert cpm_df.iloc[0, 0] == pytest.approx(166666.67, rel=1e-2)
    assert cpm_df.iloc[1, 0] == pytest.approx(333333.33, rel=1e-2)
    assert cpm_df.iloc[2, 0] == pytest.approx(500000.00, rel=1e-2)
    assert cpm_df.iloc[0, 1] == pytest.approx(266666.67, rel=1e-2)
    assert cpm_df.iloc[0, 2] == pytest.approx(291666.67, rel=1e-2)
 