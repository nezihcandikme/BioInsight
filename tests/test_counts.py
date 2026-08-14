import pandas as pd
import pytest

from deconcord.io.counts import (
    CountMatrixError,
    check_all_integer,
    validate_counts,
    load_count_matrix,
)


def test_validate_counts_invalid_raises():
    df = pd.DataFrame({"sample_A": [1.5, 2, 3]})

    with pytest.raises(CountMatrixError):
        validate_counts(df)


def test_check_all_integer_valid():
    df = pd.DataFrame({"sample_A": [1, 2, 3]})

    assert check_all_integer(df) is True


def test_check_all_integer_invalid():
    df = pd.DataFrame({"sample_A": [1.5, 2, 3]})

    assert check_all_integer(df) is False


def test_validate_counts_valid_does_not_raise():
    df = pd.DataFrame({"sample_A": [1, 2, 3]})

    validate_counts(df)


def test_validate_counts_multiple_errors_raises():
    df = pd.DataFrame({"sample_A": [-1.5, 2, 3]})

    with pytest.raises(CountMatrixError):
        validate_counts(df)


def test_load_count_matrix_invalid_raises():
    with pytest.raises(CountMatrixError):
        load_count_matrix("tests/fixtures/bad_counts.csv")


def test_load_count_matrix_valid():
    df = load_count_matrix("tests/fixtures/toy_counts.csv")
    assert df.shape == (4, 3)


def test_validate_counts_warns_on_likely_transposed_matrix():
    # Wide and short: many more "genes" (columns) than "samples" (rows) —
    # the opposite of what a real count matrix looks like, and a strong
    # hint this DataFrame was loaded transposed.
    df = pd.DataFrame(
        [[1] * 50, [2] * 50],
        columns=[f"gene{i}" for i in range(50)],
    )

    with pytest.warns(UserWarning, match="transpose"):
        validate_counts(df)


def test_validate_counts_normal_shape_does_not_warn(recwarn):
    # More genes (rows) than samples (columns), the normal RNA-seq shape —
    # the orientation heuristic should stay quiet here.
    df = pd.DataFrame({"sample_A": [1, 2, 3, 4], "sample_B": [5, 6, 7, 8]})

    validate_counts(df)

    assert len(recwarn) == 0
