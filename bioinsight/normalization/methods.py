import pandas as pd
import numpy as np
from bioinsight.core import library_size

def compute_cpm(df: pd.DataFrame) -> pd.DataFrame:
    """
Compute Counts Per Million (CPM) for a given DataFrame.

Parameters
----------
df : pd.DataFrame
    Input DataFrame with raw counts.

Returns
-------
pd.DataFrame
    DataFrame with CPM values.
"""
    total_counts = library_size(df)
    
    cpm = df.div(total_counts, axis=1) * 1e6
    
    return cpm

def log2_transform(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply log2 transformation to a given DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with raw counts.

    Returns
    -------
    pd.DataFrame
        DataFrame with log2-transformed values.
    """
    if (df < 0).any().any():
        raise ValueError("Input DataFrame contains negative values. Log2 transformation is not defined for negative numbers.")
    return np.log2(df + 1)
