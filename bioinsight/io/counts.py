import numpy as np
import pandas as pd


class CountMatrixError(Exception):
    """
    Base exception for errors found in a count matrix.

    Specific count-matrix validation errors may inherit from this class.
    """
    pass


class NonIntegerCountsError(CountMatrixError):
    """
    Exception raised when a count matrix contains non-integer columns.
    """
    pass


class NonUniqueIndexError(CountMatrixError):
    """
    Exception raised when a count matrix has a non-unique index.
    """
    pass


class NegativeCountsError(CountMatrixError):
    """
    Exception raised when a count matrix contains negative values.
    """
    pass
class MissingValuesError(CountMatrixError):
    """
    Exception raised when a count matrix contains missing values.
    """
    pass



def check_no_missing_values(df: pd.DataFrame) -> bool:
    """
    Check whether a DataFrame contains no missing values.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to validate.

    Returns
    -------
    bool
        True if there are no missing values, otherwise False.
    """
    return not df.isna().any().any()


def check_all_integer(df: pd.DataFrame) -> bool:
    """
    Check whether every column in a DataFrame has an integer data type.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to validate.

    Returns
    -------
    bool
        True if every column has an integer dtype, otherwise False.
    """
    for dtype in df.dtypes:
        if not np.issubdtype(dtype, np.integer):
            return False

    return True


def check_unique_index(df: pd.DataFrame) -> bool:
    """
    Check whether the index of a DataFrame is unique.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to validate.

    Returns
    -------
    bool
        True if the index is unique, otherwise False.
    """
    return df.index.is_unique


def check_nonnegative_counts(df: pd.DataFrame) -> bool:
    """
    Check whether all values in a DataFrame are non-negative.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to validate.

    Returns
    -------
    bool
        True if every value is greater than or equal to zero,
        otherwise False.
    """
    return (df >= 0).all().all()


def validate_counts(df: pd.DataFrame) -> None:
    """
    Validate the structure and values of a count DataFrame.

    All detected validation problems are collected and reported together.

    Parameters
    ----------
    df : pd.DataFrame
        The count DataFrame to validate.

    Raises
    ------
    CountMatrixError
        If the DataFrame contains missing values, non-integer columns,
        duplicate index values, negative counts, or multiple problems.
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
    if not errors:
        return
    elif len(errors) == 1:
        raise errors[0]
    else:
        raise CountMatrixError("\n".join(str(e) for e in errors))   
def load_count_matrix(path: str) -> pd.DataFrame:
    """
    Load a count matrix from a CSV file and validate its structure and values.

    Parameters
    ----------
    path : str
        The file path to the CSV file containing the count matrix.

    Returns
    -------
    pd.DataFrame
        The validated count matrix as a DataFrame.

    Raises
    ------
    CountMatrixError
        If the loaded DataFrame fails validation checks.
    """
    df = pd.read_csv(path, index_col=0)
    validate_counts(df)
    return df