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
def flag_outlier_samples(sizes: pd.Series, threshold: float = 3.0) -> pd.Series:
    """
    Flag outlier samples based on library sizes using the modified Z-score method.

    Parameters
    ----------
    sizes : pd.Series
        A Series containing the library sizes for each sample.
    threshold : float, optional
        The threshold for flagging outliers. Samples with a modified Z-score greater than this value will be flagged as outliers. Default is 3.0.

    Returns
    -------
    pd.Series
        A Series of boolean values indicating whether each sample is an outlier (True) or not (False).
    """
    deviations = (sizes - sizes.median()).abs()
    mad = deviations.median()
    if mad == 0: 
        return pd.Series(False, index=sizes.index)
    modified_z_scores = 0.6745 * deviations / mad
    return modified_z_scores > threshold
def run_sample_qc(df: pd.DataFrame) -> pd.DataFrame:
    """
    Run quality control metrics on the given DataFrame of gene expression data.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing gene expression data with samples as columns.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the QC metrics for each sample, including library size, genes detected, and outlier flags.
    """
    sizes = library_size(df)
    genes = genes_detected(df)

    lib_outliers = flag_outlier_samples(sizes)
    genes_outliers = flag_outlier_samples(genes)

    qc_metrics = pd.DataFrame({
    "library_size": sizes,
    "genes_detected": genes,
    "library_size_outlier": lib_outliers,
    "genes_detected_outlier": genes_outliers,
    "is_outlier": lib_outliers | genes_outliers,
    })

    return qc_metrics

   