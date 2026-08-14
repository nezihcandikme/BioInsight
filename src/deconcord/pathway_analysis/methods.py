import pandas as pd
from scipy.stats import hypergeom

from deconcord.differential_expression.methods import compute_adjusted_pvalues


def compute_enrichment_pvalue(significant_genes: set[str], gene_set: set[str], background_genes: set[str]) -> float:
    """
    Hypergeometric test: is a gene set enriched among the significant genes,
    more than expected by chance given the tested background?

    The result is only meaningful if ``background_genes`` is the actual
    universe of genes that could have been called significant (e.g. every
    gene that passed expression filtering in the DE step) — not every gene
    in the genome. An inflated or mismatched background silently changes
    what "enriched" means.

    Parameters
    ----------
    significant_genes : set[str]
        Genes called significant by the DE step.
    gene_set : set[str]
        Genes belonging to the pathway/gene set being tested.
    background_genes : set[str]
        The full tested gene universe.

    Returns
    -------
    float
        One-sided hypergeometric p-value (probability of seeing at least
        this much overlap by chance).

    Raises
    ------
    ValueError
        If ``background_genes`` is empty — a hypergeometric test against
        an empty universe has no defined answer (it silently evaluates to
        NaN otherwise, which is worse than failing loudly).
    """
    if not background_genes:
        raise ValueError(
            "background_genes is empty — enrichment needs a non-empty "
            "gene universe to test against."
        )

    # Genes outside the tested background were never eligible to be
    # "significant" in the first place, so they shouldn't be able to
    # inflate the overlap count or the significant-gene total.
    significant_genes = significant_genes & background_genes
    gene_set = gene_set & background_genes

    k = len(significant_genes & gene_set)
    n = len(significant_genes)
    K = len(gene_set)
    N = len(background_genes)

    return hypergeom.sf(k - 1, N, K, n)


def run_pathway_enrichment_analysis(significant_genes: set[str], gene_sets: dict[str, set[str]], background_genes: set[str]) -> pd.DataFrame:
    """
    Run ``compute_enrichment_pvalue`` across many gene sets at once and
    Benjamini-Hochberg-correct the results.

    Parameters
    ----------
    significant_genes : set[str]
        Genes called significant by the DE step.
    gene_sets : dict[str, set[str]]
        Gene set name -> genes in that set (e.g. pathway definitions).
    background_genes : set[str]
        The full tested gene universe.

    Returns
    -------
    pd.DataFrame
        One row per gene set, columns ``pathway``, ``p_value``,
        ``adjusted_p_value``, ``overlap_count``.

    Raises
    ------
    ValueError
        If ``background_genes`` is empty (see ``compute_enrichment_pvalue``).
    """
    if not background_genes:
        raise ValueError(
            "background_genes is empty — enrichment needs a non-empty "
            "gene universe to test against."
        )

    background_significant = significant_genes & background_genes

    p_values = {}
    overlap_counts = {}
    for pathway_name, gene_set in gene_sets.items():
        p_values[pathway_name] = compute_enrichment_pvalue(significant_genes, gene_set, background_genes)
        overlap_counts[pathway_name] = len(background_significant & (gene_set & background_genes))

    pvalues_series = pd.Series(p_values)
    adjusted = compute_adjusted_pvalues(pvalues_series)

    return pd.DataFrame({
        'pathway': list(gene_sets.keys()),
        'p_value': list(p_values.values()),
        'adjusted_p_value': list(adjusted),
        'overlap_count': list(overlap_counts.values()),
    })
