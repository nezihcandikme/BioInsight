import pandas as pd

def library_size(df: pd.DataFrame) -> pd.Series:
    """
    Calculate the library size for each sample in the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing gene expression data with samples as columns.

    Returns
    -------
    pd.Series
        A Series containing the library size for each sample.
    """
    return df.sum(axis=0)
def genes_detected(df: pd.DataFrame) -> pd.Series:
    """
    Calculate the number of genes detected for each sample in the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing gene expression data with samples as columns.

    Returns
    -------
    pd.Series
        A Series containing the number of genes detected for each sample.
    """
    return (df > 0).sum(axis=0)
