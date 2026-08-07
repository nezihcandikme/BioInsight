import pandas as pd
import pytest

from bioinsight.qc.metrics import library_size, genes_detected, flag_outlier_samples, run_sample_qc

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
def test_flag_outlier_samples():
    sizes = pd.Series({"sample_A": 10, "sample_B": 20, "sample_C": 30, "sample_D": 100})
    expected = pd.Series({"sample_A": False, "sample_B": False, "sample_C": False, "sample_D": True})
    result = flag_outlier_samples(sizes)
    pd.testing.assert_series_equal(result, expected)
def test_flag_outlier_samples_mad_zero():
    sizes = pd.Series({"sample_A": 10, "sample_B": 10, "sample_C": 10})
    expected = pd.Series({"sample_A": False, "sample_B": False, "sample_C": False})
    result = flag_outlier_samples(sizes)
    pd.testing.assert_series_equal(result, expected)
def test_run_sample_qc():
    df = pd.DataFrame({"sample_A": [1, 2, 3], "sample_B": [4, 5, 6]})
    result = run_sample_qc(df)
    assert list(result["library_size"]) == [6, 15]
    assert list(result["genes_detected"]) == [3, 3]
    assert list(result["is_outlier"]) == [False, False]
    assert list(result["library_size_outlier"]) == [False, False]
    assert list(result["genes_detected_outlier"]) == [False, False]
