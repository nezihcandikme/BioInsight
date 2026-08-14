import numpy as np
import pandas as pd
import pytest

from deconcord.io.counts import (
    CountMatrixError,
    NonUniqueIndexError,
    check_all_integer,
    check_unique_index,
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


def test_check_unique_index_true_for_unique():
    df = pd.DataFrame({"sample_A": [1, 2, 3]}, index=["gene1", "gene2", "gene3"])
    assert check_unique_index(df) is True


def test_check_unique_index_false_for_duplicates():
    df = pd.DataFrame({"sample_A": [1, 2, 3]}, index=["gene1", "gene1", "gene3"])
    assert check_unique_index(df) is False


def test_validate_counts_duplicate_gene_ids_raises():
    # Two rows sharing a gene ID is a real, easy-to-hit input mistake
    # (e.g. concatenating two annotation releases) -- validate_counts must
    # catch it rather than silently letting duplicate-indexed rows through
    # to functions that assume a unique index (like .loc lookups downstream).
    df = pd.DataFrame(
        {"sample_A": [1, 2, 3], "sample_B": [4, 5, 6]},
        index=["gene1", "gene1", "gene2"],
    )

    with pytest.raises(NonUniqueIndexError):
        validate_counts(df)


def test_validate_counts_infinite_values_raises():
    # An infinite value can't be an integer count, so this is already
    # rejected today -- but via the generic "non-integer" path, since
    # float('inf') gives the column a float dtype. Pinned down explicitly
    # here so a future refactor of check_all_integer can't silently let
    # inf slip through as if it were a valid float count.
    df = pd.DataFrame({"sample_A": [1.0, np.inf, 3.0]})

    with pytest.raises(CountMatrixError):
        validate_counts(df)


def test_validate_counts_all_zero_column_is_valid():
    # A sample where every gene has zero counts is unusual but not
    # malformed on its own -- validate_counts only checks structural/type
    # validity. (Normalization downstream, e.g. compute_cpm, is where an
    # all-zero sample is rejected, since CPM is undefined for it.)
    df = pd.DataFrame({"sample_A": [0, 0, 0], "sample_B": [1, 2, 3]})

    validate_counts(df)
