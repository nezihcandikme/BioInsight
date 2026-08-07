import pandas as pd
import pytest

from bioinsight.qc._init_ import library_size, genes_detected

def test_library_size():
    df = pd.DataFrame({"sample_A": [1, 2, 3], "sample_B": [4, 5, 6]})
    expected = pd.Series({"sample_A": 6, "sample_B": 15})
    result = library_size(df)
    pd.testing.assert_series_equal(result, expected)

def test_genes_detected():
    df = pd.DataFrame({"sample_A": [1, 0, 3], "sample_B": [0, 5, 0]})
    expected = pd.Series({"sample_A": 2, "sample_B": 1})
    result = genes_detected(df)
    pd.testing.assert_series_equal(result, expected)
