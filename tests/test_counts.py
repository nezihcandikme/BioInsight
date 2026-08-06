import pandas as pd
import pytest

from bioinsight.io.counts import (
    NonIntegerCountsError,
    check_all_integer,
    validate_counts,
)


def test_check_all_integer_valid():
    df = pd.DataFrame({"sample_A": [1, 2, 3]})

    assert check_all_integer(df) is True


def test_check_all_integer_invalid():
    df = pd.DataFrame({"sample_A": [1.5, 2, 3]})

    assert check_all_integer(df) is False


def test_validate_counts_valid_does_not_raise():
    df = pd.DataFrame({"sample_A": [1, 2, 3]})

    validate_counts(df)


def test_validate_counts_invalid_raises():
    df = pd.DataFrame({"sample_A": [1.5, 2, 3]})

    with pytest.raises(NonIntegerCountsError):
        validate_counts(df)