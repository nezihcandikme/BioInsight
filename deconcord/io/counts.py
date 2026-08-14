import warnings

import numpy as np
import pandas as pd

# A real RNA-seq count matrix is lopsided: thousands of genes, a much
# smaller handful of samples. If a matrix has more columns than rows and
# very few rows, it's a decent bet someone loaded it transposed.
_ORIENTATION_WARNING_MAX_ROWS = 30


class CountMatrixError(Exception):
    """Base exception for problems found in a count matrix."""
    pass


class NonIntegerCountsError(CountMatrixError):
    """Raised when a count matrix contains non-integer columns."""
    pass


class NonUniqueIndexError(CountMatrixError):
    """Raised when a count matrix has a non-unique (duplicated) gene index."""
    pass


class NegativeCountsError(CountMatrixError):
    """Raised when a count matrix contains negative values."""
    pass


class MissingValuesError(CountMatrixError):
    """Raised when a count matrix contains missing (NaN) values."""
    pass


def check_no_missing_values(df: pd.DataFrame) -> bool:
    """True if the DataFrame has no missing values."""
    return not df.isna().any().any()


def check_all_integer(df: pd.DataFrame) -> bool:
    """True if every column has an integer dtype."""
    for dtype in df.dtypes:
        if not np.issubdtype(dtype, np.integer):
            return False

    return True


def check_unique_index(df: pd.DataFrame) -> bool:
    """True if the DataFrame's index (gene IDs) has no duplicates."""
    return df.index.is_unique


def check_nonnegative_counts(df: pd.DataFrame) -> bool:
    """True if every value in the DataFrame is >= 0."""
    return (df >= 0).all().all()


def check_likely_correct_orientation(df: pd.DataFrame) -> bool:
    """
    Heuristic check that ``df`` looks like genes-as-rows, samples-as-columns
    (the orientation every other function in DEConcord assumes).

    This can't be known for certain from the numbers alone — it's a sanity
    check, not a proof. A matrix with more columns than rows and only a
    handful of rows is flagged as "probably transposed" because that shape
    is far more consistent with samples-as-rows than with a real gene list.
    A legitimately tiny custom gene panel would also trip this, so treat a
    False result as "go double check your orientation," not as evidence of
    a bug on its own.
    """
    n_rows, n_cols = df.shape
    return not (n_cols > n_rows and n_rows <= _ORIENTATION_WARNING_MAX_ROWS)


def validate_counts(df: pd.DataFrame) -> None:
    """
    Validate the structure and values of a raw count matrix.

    Expects genes as rows and samples as columns, with non-negative integer
    values and no missing data. All problems found are collected and raised
    together (rather than stopping at the first one) so a single failed run
    tells you everything that's wrong instead of one error at a time.

    Also emits a ``UserWarning`` (not a hard failure) if the matrix's shape
    looks like it might be transposed — see ``check_likely_correct_orientation``.

    Parameters
    ----------
    df : pd.DataFrame
        The count matrix to validate.

    Raises
    ------
    CountMatrixError
        If the matrix contains missing values, non-integer columns,
        duplicate gene IDs, negative counts, or a combination of these.
    """
    errors = []

    if not check_no_missing_values(df):
        errors.append(MissingValuesError("Counts must not contain missing values."))

    if not check_all_integer(df):
        bad_cols = df.dtypes[
            ~df.dtypes.apply(
                lambda dtype: np.issubdtype(dtype, np.integer)
            )
        ].index.tolist()

        errors.append(NonIntegerCountsError(
            f"Columns {bad_cols} must have integer data types."
        ))

    if not check_unique_index(df):
        errors.append(NonUniqueIndexError("The index of the DataFrame must be unique."))

    if not check_nonnegative_counts(df):
        errors.append(NegativeCountsError("Counts must be non-negative."))

    if len(errors) == 1:
        raise errors[0]
    elif errors:
        raise CountMatrixError("\n".join(str(e) for e in errors))

    if not check_likely_correct_orientation(df):
        warnings.warn(
            f"This count matrix has {df.shape[1]} columns and only "
            f"{df.shape[0]} rows. DEConcord expects genes as rows and "
            "samples as columns — if that's backwards here, transpose "
            "your input (df.T) before passing it in.",
            stacklevel=2,
        )


def load_count_matrix(path: str) -> pd.DataFrame:
    """
    Load a count matrix from a CSV file and validate it via ``validate_counts``.

    The first column of the CSV is used as the gene ID index; every other
    column is treated as a sample.

    Parameters
    ----------
    path : str
        Path to the CSV file.

    Returns
    -------
    pd.DataFrame
        The validated, raw (not normalized) count matrix.

    Raises
    ------
    CountMatrixError
        If the loaded DataFrame fails validation.
    """
    df = pd.read_csv(path, index_col=0)
    validate_counts(df)
    return df
