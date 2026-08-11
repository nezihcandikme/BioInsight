import pandas as pd

from bioinsight.core import library_size, genes_detected


def flag_outlier_samples(sizes: pd.Series, threshold: float = 3.0) -> pd.Series:
    """
    Flag samples whose value is an outlier relative to the rest, using the
    median absolute deviation (MAD)-based modified Z-score.

    MAD is used instead of standard deviation because a couple of extreme
    samples (the exact thing we're trying to catch) would otherwise drag
    the mean and standard deviation along with them and hide themselves.

    Parameters
    ----------
    sizes : pd.Series
        Per-sample values to check (e.g. library size or genes detected).
    threshold : float, optional
        Modified Z-score above which a sample is flagged. Default 3.0,
        a common rule-of-thumb cutoff (Iglewicz & Hoaglin).

    Returns
    -------
    pd.Series
        Boolean flags, indexed the same as ``sizes``.
    """
    deviations = (sizes - sizes.median()).abs()
    mad = deviations.median()

    # If every sample has (near) the same value, MAD is 0 and dividing by
    # it would produce inf/NaN "outliers" for meaningless noise. There's
    # nothing to flag when nothing actually varies.
    if mad == 0:
        return pd.Series(False, index=sizes.index)

    modified_z_scores = 0.6745 * deviations / mad
    return modified_z_scores > threshold


def run_sample_qc(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-sample QC metrics: library size, genes detected, and
    MAD-based outlier flags for each.

    Parameters
    ----------
    df : pd.DataFrame
        Raw count matrix, genes as rows, samples as columns.

    Returns
    -------
    pd.DataFrame
        Indexed by sample, with columns ``library_size``, ``genes_detected``,
        ``library_size_outlier``, ``genes_detected_outlier``, and
        ``is_outlier`` (the OR of the two outlier flags).
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
