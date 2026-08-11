import pandas as pd
import pytest

from bioinsight.normalization.methods import compute_cpm, log2_transform


def test_compute_cpm():
    data = {
        'Sample1': [100, 200, 300],
        'Sample2': [400, 500, 600],
        'Sample3': [700, 800, 900]
    }
    df = pd.DataFrame(data)

    cpm_df = compute_cpm(df)

    assert cpm_df.iloc[0, 0] == pytest.approx(166666.67, rel=1e-2)
    assert cpm_df.iloc[1, 0] == pytest.approx(333333.33, rel=1e-2)
    assert cpm_df.iloc[2, 0] == pytest.approx(500000.00, rel=1e-2)
    assert cpm_df.iloc[0, 1] == pytest.approx(266666.67, rel=1e-2)
    assert cpm_df.iloc[0, 2] == pytest.approx(291666.67, rel=1e-2)


def test_compute_cpm_zero_total_raises():
    df = pd.DataFrame({
        'Sample1': [100, 200, 300],
        'Sample2': [0, 0, 0],
    })

    with pytest.raises(ValueError, match="Sample2"):
        compute_cpm(df)


def test_compute_cpm_zero_total_error_names_the_sample():
    df = pd.DataFrame({
        'good_sample': [10, 20, 30],
        'empty_sample': [0, 0, 0],
    })

    with pytest.raises(ValueError) as excinfo:
        compute_cpm(df)

    assert "empty_sample" in str(excinfo.value)
    assert "good_sample" not in str(excinfo.value)


def test_log2_transform():
    data = {
        'Sample1': [0, 1, 2],
        'Sample2': [3, 4, 5],
        'Sample3': [6, 7, 8]
    }
    df = pd.DataFrame(data)

    log2_df = log2_transform(df)

    expected_log2_df = pd.DataFrame({
        'Sample1': [0.0, 1.0, 1.584962500721156],
        'Sample2': [2.0, 2.321928094887362, 2.584962500721156],
        'Sample3': [2.807354922057604, 3.0, 3.169925001442312]
    })
    pd.testing.assert_frame_equal(log2_df, expected_log2_df, check_exact=False, rtol=1e-3)


def test_log2_transform_negative_values_raises():
    df = pd.DataFrame({'Sample1': [1, -2, 3]})

    with pytest.raises(ValueError):
        log2_transform(df)
