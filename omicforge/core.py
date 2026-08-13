import pandas as pd


def library_size(df: pd.DataFrame) -> pd.Series:
    """Total counts per sample (column sums)."""
    return df.sum(axis=0)


def genes_detected(df: pd.DataFrame) -> pd.Series:
    """Number of genes with a nonzero count in each sample."""
    return (df > 0).sum(axis=0)
