import numpy as np
import pandas as pd


class NonIntegerCountsError(Exception):
    """
    Exception raised when a count matrix contains non-integer values.

    This error is used to indicate that one or more columns in a DataFrame
    contain values that are not of an integer data type, even though
    integer counts are required.
    """
    pass


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
        True if all columns have integer dtypes, otherwise False.
    """
    for dtype in df.dtypes:
        if not np.issubdtype(dtype, np.integer):
            return False
    return True


def validate_counts(df: pd.DataFrame) -> None:
    """
    Validate that a count DataFrame contains only integer-valued columns.

    Parameters
    ----------
    df : pd.DataFrame
        The count DataFrame to validate.

    Raises
    ------
    NonIntegerCountsError
        If any column does not have an integer data type.
    """
    if not check_all_integer(df):
        bad_cols = df.dtypes[
            ~df.dtypes.apply(lambda x: np.issubdtype(x, np.integer))
        ].index.tolist()

        raise NonIntegerCountsError(
            f"Counts must be integers. Non-integer columns: {bad_cols}"
        )