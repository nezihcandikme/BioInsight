"""
Differential expression analysis — basic, exploratory method.

The functions here implement a mean-difference metric plus a per-gene
Welch's t-test. This is NOT a count-aware negative-binomial model like
DESeq2 or edgeR: it doesn't model biological dispersion or the
mean-variance relationship that raw RNA-seq counts actually have, and it
performs no normalization of its own.

Practical implications:
- ``df`` is expected to already be normalized and log-transformed (e.g.
  ``compute_cpm`` followed by ``log2_transform`` from
  ``omicforge.normalization.methods``) before being passed in here. Feed
  it raw counts and "log fold change" becomes a difference of raw count
  means — a real number, but not a fold change of anything.
- Treat this module's output as a fast, exploratory first look. For
  results meant to support a real biological conclusion, use DESeq2 or
  edgeR, which OmicForge does not replace.
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq
from scipy.special import digamma, polygamma
from statsmodels.stats.multitest import multipletests


def _validate_groups(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> None:
    """
    Shared input checks for group_1/group_2, used by every public function
    in this module. Exists so a typo'd sample name or an accidentally
    reused sample fails loudly and specifically, instead of quietly
    producing a nonsense comparison or a confusing pandas KeyError.
    """
    if not group_1:
        raise ValueError("group_1 is empty — it needs at least one sample name.")
    if not group_2:
        raise ValueError("group_2 is empty — it needs at least one sample name.")

    for group_name, group in (("group_1", group_1), ("group_2", group_2)):
        seen = set()
        for sample in group:
            if sample in seen:
                raise ValueError(
                    f"Sample '{sample}' appears more than once in {group_name}."
                )
            seen.add(sample)

    overlap = set(group_1) & set(group_2)
    if overlap:
        raise ValueError(
            f"Sample(s) {sorted(overlap)} appear in both group_1 and group_2 "
            "— a sample can't be compared against itself."
        )

    for group_name, group in (("group_1", group_1), ("group_2", group_2)):
        for sample in group:
            if sample not in df.columns:
                raise ValueError(
                    f"Sample '{sample}' was requested in {group_name} but is "
                    "not present in the count matrix."
                )


def compute_log_fold_change(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> pd.Series:
    """
    Compute the difference in mean expression between two groups, per gene.

    This is only a genuine "log fold change" if ``df`` is already on a
    log2 scale — a difference of logs is a log ratio. Pass in raw or
    linear-scale CPM values and this becomes a difference of means, which
    is a different (and less standard) quantity. This is also a simple
    plug-in estimate, not the shrunken effect size DESeq2/edgeR report.

    Parameters
    ----------
    df : pd.DataFrame
        Log-transformed, normalized expression matrix (genes as rows,
        samples as columns).
    group_1, group_2 : list[str]
        Sample names (column names in ``df``) for each comparison group.

    Returns
    -------
    pd.Series
        Per-gene ``mean(group_1) - mean(group_2)``, indexed by gene.
    """
    _validate_groups(df, group_1, group_2)

    mean_group_1 = df[group_1].mean(axis=1)
    mean_group_2 = df[group_2].mean(axis=1)

    return mean_group_1 - mean_group_2


def compute_pvalues(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> pd.Series:
    """
    Per-gene Welch's t-test between two groups.

    Welch's t-test (unequal-variance t-test) is used instead of the
    standard Student's t-test because there's no reason to assume the two
    groups have equal variance — and RNA-seq expression variance is
    routinely condition-dependent, so assuming equal variance would be the
    wrong default.

    Parameters
    ----------
    df : pd.DataFrame
        Log-transformed, normalized expression matrix (genes as rows,
        samples as columns).
    group_1, group_2 : list[str]
        Sample names (column names in ``df``) for each comparison group.
        Each group needs at least 2 samples — a t-test needs within-group
        variance to estimate, and one sample can't provide that.

    Returns
    -------
    pd.Series
        Raw (not multiple-testing-corrected) p-value per gene.

    Raises
    ------
    ValueError
        If either group has fewer than 2 samples, or the group inputs fail
        the shared validation in ``_validate_groups`` (missing/duplicate/
        overlapping sample names).
    """
    _validate_groups(df, group_1, group_2)

    if len(group_1) < 2 or len(group_2) < 2:
        raise ValueError("Each group must have at least 2 samples to perform a t-test.")

    def test_one_gene(row):
        group_1_values = row[group_1].astype(float)
        group_2_values = row[group_2].astype(float)

        # Welch's t-test is undefined (NaN) when both groups have zero
        # variance — e.g. an all-zero or otherwise constant-expression
        # gene. Rather than let NaN silently leak into the results table,
        # resolve it directly: identical constants means "not different"
        # (p=1), differing constants means "as different as it gets" (p=0).
        if group_1_values.var(ddof=1) == 0 and group_2_values.var(ddof=1) == 0:
            means_equal = np.isclose(group_1_values.mean(), group_2_values.mean())
            return 1.0 if means_equal else 0.0

        _, p_value = stats.ttest_ind(group_1_values, group_2_values, equal_var=False)
        return p_value

    return df.apply(test_one_gene, axis=1)


def _pooled_variance_and_df(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> tuple[pd.Series, int]:
    """
    Per-gene pooled (equal-variance-assumed) sample variance and its
    residual degrees of freedom — the input the moderated t-test needs.

    Unlike ``compute_pvalues``'s Welch's t-test, this assumes the two
    groups share one within-gene variance. That's a real assumption this
    method makes and Welch's doesn't; it's what makes the variance a
    single per-gene number that can be shrunk toward a value shared
    across genes in the first place.
    """
    n1, n2 = len(group_1), len(group_2)
    df_resid = n1 + n2 - 2

    var_1 = df[group_1].var(axis=1, ddof=1)
    var_2 = df[group_2].var(axis=1, ddof=1)
    pooled_var = ((n1 - 1) * var_1 + (n2 - 1) * var_2) / df_resid

    return pooled_var, df_resid


def _fit_ebayes_prior(sample_variances: pd.Series, df_resid: int) -> tuple[float, float]:
    """
    Method-of-moments fit of the empirical Bayes prior (d0, s0^2) that the
    moderated t-test shrinks every gene's variance toward.

    This is the same idea behind limma's ``squeezeVar``/``eBayes`` (Smyth,
    2004, *Statistical Applications in Genetics and Molecular Biology*):
    model each gene's true variance as drawn from a shared scaled-inverse-
    chi-squared prior, and estimate that prior's degrees of freedom (d0)
    and scale (s0^2) from how much the observed per-gene variances vary
    across genes versus how much variation is expected from sampling noise
    alone. Simplified relative to limma's actual ``fitFDist``: this uses a
    direct method-of-moments solve rather than an iterative MLE, and
    assumes every gene shares the same residual degrees of freedom (true
    for OmicForge's fixed two-group design).

    Parameters
    ----------
    sample_variances : pd.Series
        Per-gene pooled variance, e.g. from ``_pooled_variance_and_df``.
    df_resid : int
        Residual degrees of freedom shared by every gene's variance
        estimate.

    Returns
    -------
    tuple[float, float]
        ``(d0, s0_sq)`` — prior degrees of freedom and prior variance.

    Raises
    ------
    ValueError
        If fewer than 2 genes have nonzero variance. The whole point of
        this method is borrowing information *across* genes, so it has
        nothing to borrow from with only one usable gene.
    """
    positive = sample_variances[sample_variances > 0]
    if len(positive) < 2:
        raise ValueError(
            "Need at least 2 genes with nonzero within-group variance to fit "
            f"an empirical Bayes prior for the moderated t-test; got {len(positive)}."
        )

    log_variances = np.log(positive.to_numpy())
    mean_log_variance = log_variances.mean()
    var_of_log_variance = log_variances.var(ddof=1)

    # trigamma(df_resid / 2) is the sampling variance of log(s^2) expected
    # from residual degrees of freedom alone, with no real gene-to-gene
    # variability at all. Anything above that is attributed to genes
    # genuinely having different true variances.
    trigamma_resid = polygamma(1, df_resid / 2)
    excess_variance = var_of_log_variance - trigamma_resid

    d0_cap = 1e6  # stands in for "no detectable cross-gene variability, shrink fully"
    if excess_variance <= 1e-12:
        d0 = d0_cap
    else:
        def trigamma_gap(candidate_d0):
            return polygamma(1, candidate_d0 / 2) - excess_variance

        d0 = brentq(trigamma_gap, 1e-6, d0_cap)

    log_s0_sq = (
        mean_log_variance
        - digamma(df_resid / 2) + np.log(df_resid / 2)
        - digamma(d0 / 2) + np.log(d0 / 2)
    )
    s0_sq = np.exp(log_s0_sq)

    return d0, s0_sq


def compute_moderated_pvalues(df: pd.DataFrame, group_1: list[str], group_2: list[str]) -> pd.Series:
    """
    Per-gene moderated (empirical-Bayes-shrunk) t-test between two groups.

    Where ``compute_pvalues`` tests every gene in isolation, this borrows
    statistical power across genes: each gene's variance estimate is
    shrunk toward a value fit from every other gene's variance, exactly
    the mechanism DESeq2 and edgeR use to detect real effects that a
    naive per-gene test misses in small experiments (see
    ``benchmarks/README.md`` for the measured size of that gap). The
    tradeoff for the extra power is the assumption ``compute_pvalues``
    deliberately avoids: that a gene's two groups share one variance,
    not each their own. Genes with few, noisy replicates benefit most;
    genes with many consistent replicates are barely shrunk at all
    (their own data already dominates the prior).

    A side effect worth knowing about: because every gene's variance is
    shrunk toward a prior that's never exactly zero, a constant-expression
    gene no longer needs the ad hoc zero-variance special case
    ``compute_pvalues`` requires — it gets a well-defined, usually large
    but finite, p-value automatically.

    Parameters
    ----------
    df : pd.DataFrame
        Log-transformed, normalized expression matrix (genes as rows,
        samples as columns).
    group_1, group_2 : list[str]
        Sample names (column names in ``df``) for each comparison group.
        Each group needs at least 2 samples, and at least 2 genes overall
        need nonzero variance for the prior fit.

    Returns
    -------
    pd.Series
        Raw (not multiple-testing-corrected) p-value per gene.

    Raises
    ------
    ValueError
        If either group has fewer than 2 samples, the group inputs fail
        ``_validate_groups``, or fewer than 2 genes have nonzero variance
        (see ``_fit_ebayes_prior``).
    """
    _validate_groups(df, group_1, group_2)

    if len(group_1) < 2 or len(group_2) < 2:
        raise ValueError("Each group must have at least 2 samples to perform a moderated t-test.")

    pooled_var, df_resid = _pooled_variance_and_df(df, group_1, group_2)
    d0, s0_sq = _fit_ebayes_prior(pooled_var, df_resid)

    moderated_var = (d0 * s0_sq + df_resid * pooled_var) / (d0 + df_resid)

    n1, n2 = len(group_1), len(group_2)
    mean_diff = df[group_1].mean(axis=1) - df[group_2].mean(axis=1)
    standard_error = np.sqrt(moderated_var * (1 / n1 + 1 / n2))

    t_statistic = mean_diff / standard_error
    moderated_df = d0 + df_resid
    p_values = 2 * stats.t.sf(np.abs(t_statistic), df=moderated_df)

    return pd.Series(p_values, index=df.index)


def compute_adjusted_pvalues(pvalues: pd.Series, method: str = 'fdr_bh') -> pd.Series:
    """
    Correct p-values for multiple testing.

    Testing thousands of genes at once means thousands of chances for a
    p < 0.05 result by pure chance — without correction, a typical RNA-seq
    experiment would report hundreds of "significant" genes even if
    nothing in the biology changed. Benjamini-Hochberg (the default here)
    controls the expected false discovery rate rather than the much
    stricter per-test error rate, which is the standard, less
    conservative choice for exploratory genomics.

    Parameters
    ----------
    pvalues : pd.Series
        Raw per-gene p-values.
    method : str, optional
        Passed through to ``statsmodels.stats.multitest.multipletests``.
        Default ``'fdr_bh'`` (Benjamini-Hochberg).

    Returns
    -------
    pd.Series
        Adjusted p-values, same index as ``pvalues``.

    Raises
    ------
    ValueError
        If ``pvalues`` contains any NaN. ``statsmodels`` doesn't just skip
        the offending gene here — a single NaN input silently turns every
        adjusted p-value in the whole Series into NaN, which is exactly
        the kind of quiet corruption this function should never allow
        through. ``compute_pvalues`` in this module already avoids
        producing NaN, so this mainly guards against p-values coming from
        somewhere else.
    """
    # statsmodels divides by len(pvals) internally and doesn't guard
    # against an empty array itself (ZeroDivisionError), so handle the
    # "nothing to correct" case here instead of letting that surface.
    if pvalues.empty:
        return pd.Series(dtype=float, index=pvalues.index)

    if pvalues.isna().any():
        bad_genes = pvalues[pvalues.isna()].index.tolist()
        raise ValueError(
            f"pvalues contains NaN for gene(s) {bad_genes}. Multiple-testing "
            "correction can't proceed — a single NaN p-value poisons every "
            "adjusted p-value in the batch, not just its own row."
        )

    adjusted_pvalues = multipletests(pvalues, method=method)[1]
    return pd.Series(adjusted_pvalues, index=pvalues.index)


_PVALUE_METHODS = {
    "welch": compute_pvalues,
    "moderated": compute_moderated_pvalues,
}


def run_differential_expression_with_covariates(
    df: pd.DataFrame,
    group_1: list[str],
    group_2: list[str],
    metadata: pd.DataFrame,
    covariate_cols: list[str],
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
    moderated: bool = True,
) -> pd.DataFrame:
    """
    Per-gene multiple linear regression DE test that adjusts for one or
    more covariates, instead of comparing group_1 vs group_2 in isolation.

    ``compute_pvalues``/``compute_moderated_pvalues`` treat every
    within-group difference as noise and have no way to account for a
    known source of variation that isn't the condition of interest -- a
    batch, a subject's sex, an RNA extraction date. If group membership
    happens to correlate with one of those, the simple two-group test
    can't tell the real effect from the confound. This function instead
    fits, per gene::

        expression ~ intercept + condition + covariate_1 + covariate_2 + ...

    by ordinary least squares, and tests the *condition* coefficient --
    the group_1-vs-group_2 effect *holding the covariates fixed*.
    Categorical covariates (non-numeric columns) are one-hot encoded with
    the alphabetically first level dropped as the reference; numeric
    covariates are used as-is (a linear covariate effect is assumed, not
    tested).

    There's no direct equivalent of ``compute_pvalues``'s Welch's t-test
    here: Welch's test is specifically a two-independent-samples test and
    doesn't generalize to "compare two levels while holding other
    variables fixed" the way a linear model does. What this function
    offers instead is the same equal-variance-per-gene linear model
    ``compute_moderated_pvalues`` uses, generalized from a plain two-group
    design to an arbitrary design matrix, with the same optional
    empirical-Bayes variance shrinkage (Smyth, 2004) across genes.

    Parameters
    ----------
    df : pd.DataFrame
        Log-transformed, normalized expression matrix (genes as rows,
        samples as columns).
    group_1, group_2 : list[str]
        Sample names for the two condition levels being compared. The
        reported log fold change is the group_1-vs-group_2 coefficient
        after adjusting for the covariates -- same sign convention as
        ``compute_log_fold_change`` (group_1 relative to group_2), but not
        the same quantity: this is a partial regression coefficient, not
        a raw mean difference, unless every covariate happens to be
        constant.
    metadata : pd.DataFrame
        Indexed by sample name, covering every sample in
        ``group_1 + group_2``. Must contain every column named in
        ``covariate_cols``.
    covariate_cols : list[str]
        Column names in ``metadata`` to include as covariates. At least
        one is required -- with none, use ``run_differential_expression``
        instead, which is simpler and doesn't carry a linear model's
        extra assumptions.
    alpha : float, optional
        Adjusted p-value cutoff for the ``significant`` column. Default 0.05.
    lfc_threshold : float, optional
        Minimum ``abs(log_fold_change)`` for the ``significant`` column.
        Default 1.0.
    moderated : bool, optional
        Whether to shrink each gene's residual variance toward a prior fit
        across every other gene's residual variance, the same mechanism
        ``compute_moderated_pvalues`` uses. Default True. Set False for a
        plain per-gene OLS t-test with no cross-gene borrowing -- useful
        with very few genes, where a prior can't be fit (see
        ``_fit_ebayes_prior``).

    Returns
    -------
    pd.DataFrame
        Indexed by gene, with columns ``log_fold_change``, ``p_value``,
        ``adjusted_p_value``, ``significant`` -- same meaning as
        ``run_differential_expression``, except ``log_fold_change`` is the
        covariate-adjusted coefficient described above.

    Raises
    ------
    ValueError
        If ``covariate_cols`` is empty, a sample from ``group_1``/
        ``group_2`` is missing from ``metadata``, a covariate column is
        missing or contains NaN for a used sample, or the resulting
        design matrix is rank-deficient (e.g. a covariate that's constant,
        or collinear with group membership -- a linear model can't
        separate two effects that never vary independently in the data).
    """
    if not (0 < alpha <= 1):
        raise ValueError(f"alpha must be in (0, 1]; got {alpha}.")
    if lfc_threshold < 0:
        raise ValueError(f"lfc_threshold must be >= 0; got {lfc_threshold}.")

    _validate_groups(df, group_1, group_2)

    if not covariate_cols:
        raise ValueError(
            "covariate_cols is empty -- with no covariates, use "
            "run_differential_expression instead."
        )

    samples = group_1 + group_2

    missing_samples = [s for s in samples if s not in metadata.index]
    if missing_samples:
        raise ValueError(
            f"Sample(s) {missing_samples} appear in group_1/group_2 but are "
            "not present in metadata's index."
        )

    missing_cols = [c for c in covariate_cols if c not in metadata.columns]
    if missing_cols:
        raise ValueError(f"covariate_cols {missing_cols} not found in metadata's columns.")

    meta = metadata.loc[samples, covariate_cols]
    if meta.isna().any().any():
        bad_samples = meta.index[meta.isna().any(axis=1)].tolist()
        raise ValueError(f"metadata has missing covariate values for sample(s) {bad_samples}.")

    condition = pd.Series(
        {**{s: 1.0 for s in group_1}, **{s: 0.0 for s in group_2}}
    )[samples]

    design = pd.DataFrame({"condition": condition}, index=samples)
    for col in covariate_cols:
        column = meta[col]
        if pd.api.types.is_numeric_dtype(column):
            design[col] = column.astype(float)
        else:
            dummies = pd.get_dummies(column.astype(str), prefix=col, drop_first=True, dtype=float)
            if dummies.shape[1] == 0:
                # A single-level categorical produces zero dummy columns
                # (drop_first drops the only level there is) -- it would
                # otherwise silently vanish from the design matrix instead
                # of being caught as "nothing to adjust for."
                raise ValueError(
                    f"Covariate '{col}' has only one distinct value across "
                    "these samples -- there's nothing to adjust for."
                )
            dummies.index = samples
            design = design.join(dummies)
    design.insert(0, "intercept", 1.0)

    X = design.to_numpy(dtype=float)
    n_samples, n_params = X.shape
    rank = int(np.linalg.matrix_rank(X))
    if rank < n_params:
        raise ValueError(
            f"The design matrix is rank-deficient (rank {rank} < {n_params} "
            "parameters) -- a covariate is likely constant or collinear "
            "with group_1/group_2 membership, so their individual effects "
            "can't be separated."
        )

    df_resid = n_samples - rank
    if df_resid < 1:
        raise ValueError(
            f"Not enough samples ({n_samples}) for {n_params} design "
            "parameters -- need at least one residual degree of freedom."
        )

    condition_idx = design.columns.get_loc("condition")

    Y = df[samples].to_numpy(dtype=float).T  # samples x genes
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ Y  # n_params x genes

    residuals = Y - X @ beta
    sigma_sq = np.sum(residuals ** 2, axis=0) / df_resid  # per gene

    if moderated:
        d0, s0_sq = _fit_ebayes_prior(pd.Series(sigma_sq, index=df.index), df_resid)
        variance_for_test = (d0 * s0_sq + df_resid * sigma_sq) / (d0 + df_resid)
        test_df = d0 + df_resid
    else:
        variance_for_test = sigma_sq
        test_df = df_resid

    log_fold_change = beta[condition_idx]
    standard_error = np.sqrt(variance_for_test * XtX_inv[condition_idx, condition_idx])
    t_statistic = log_fold_change / standard_error
    p_values = 2 * stats.t.sf(np.abs(t_statistic), df=test_df)

    results_df = pd.DataFrame({
        "log_fold_change": pd.Series(log_fold_change, index=df.index),
        "p_value": pd.Series(p_values, index=df.index),
    })
    results_df["adjusted_p_value"] = compute_adjusted_pvalues(results_df["p_value"])
    results_df["significant"] = (
        (results_df["adjusted_p_value"] < alpha)
        & (results_df["log_fold_change"].abs() > lfc_threshold)
    )

    return results_df


def run_differential_expression(
    df: pd.DataFrame,
    group_1: list[str],
    group_2: list[str],
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
    method: str = "welch",
) -> pd.DataFrame:
    """
    Run the full exploratory DE workflow: log fold change, a per-gene
    significance test, and Benjamini-Hochberg correction, combined into
    one results table.

    Not a substitute for DESeq2/edgeR — see the module docstring. ``df``
    should already be normalized and log-transformed (e.g. log2-CPM via
    ``omicforge.normalization.methods``).

    Parameters
    ----------
    df : pd.DataFrame
        Log-transformed, normalized expression matrix (genes as rows,
        samples as columns).
    group_1, group_2 : list[str]
        Sample names (column names in ``df``) for each comparison group.
    alpha : float, optional
        Adjusted p-value cutoff for the ``significant`` column. Default
        0.05. There's nothing statistically special about 0.05 — it's a
        convention, not a law — so this is exposed rather than buried.
    lfc_threshold : float, optional
        Minimum ``abs(log_fold_change)`` for the ``significant`` column.
        Default 1.0 (i.e. a 2x difference on a log2 scale). Only meaningful
        as a fold-change cutoff if ``df`` is actually log2-scaled — see the
        module docstring.
    method : str, optional
        ``"welch"`` (default) tests every gene in isolation with Welch's
        t-test — no assumptions borrowed across genes, but limited
        statistical power in small experiments. ``"moderated"`` shrinks
        each gene's variance toward a value estimated from every other
        gene (see ``compute_moderated_pvalues``) — more power, at the cost
        of assuming the two groups share one variance per gene rather than
        each their own. Neither is "the correct one"; they trade different
        assumptions for different power, and ``benchmarks/`` has real
        numbers on the difference.

    Returns
    -------
    pd.DataFrame
        Indexed by gene, with columns:
        - ``log_fold_change``: mean(group_1) - mean(group_2)
        - ``p_value``: raw per-gene p-value (Welch's or moderated,
          depending on ``method``)
        - ``adjusted_p_value``: Benjamini-Hochberg corrected p-value
        - ``significant``: ``adjusted_p_value < alpha`` and
          ``abs(log_fold_change) > lfc_threshold`` — a starting point for
          further inspection, not a verdict, regardless of what you set
          the thresholds to.

    Raises
    ------
    ValueError
        If ``alpha`` is not in (0, 1], ``lfc_threshold`` is negative, or
        ``method`` isn't ``"welch"`` or ``"moderated"``.
    """
    if not (0 < alpha <= 1):
        raise ValueError(f"alpha must be in (0, 1]; got {alpha}.")
    if lfc_threshold < 0:
        raise ValueError(f"lfc_threshold must be >= 0; got {lfc_threshold}.")
    if method not in _PVALUE_METHODS:
        raise ValueError(f"method must be one of {sorted(_PVALUE_METHODS)}; got {method!r}.")

    log_fold_change = compute_log_fold_change(df, group_1, group_2)
    pvalues = _PVALUE_METHODS[method](df, group_1, group_2)
    adjusted_pvalues = compute_adjusted_pvalues(pvalues)

    results_df = pd.DataFrame({
        'log_fold_change': log_fold_change,
        'p_value': pvalues,
        'adjusted_p_value': adjusted_pvalues
    })
    results_df["significant"] = (results_df["adjusted_p_value"] < alpha) & (results_df["log_fold_change"].abs() > lfc_threshold)

    return results_df
