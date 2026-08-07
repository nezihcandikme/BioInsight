import pandas as pd
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
