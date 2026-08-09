from scipy.stats import hypergeom
import pandas as pd

from bioinsight.differential_expression.methods import compute_adjusted_pvalues

def compute_enrichment_pvalue(significant_genes: set[str], gene_set: set[str], background_genes: set[str]) -> float:
    """
Compute the enrichment p-value using the hypergeometric test.

Parameters
----------
significant_genes : set[str]
    A set of significant genes.
gene_set : set[str]
    A set of genes in the pathway or gene set of interest.
background_genes : set[str]
    A set of all background genes (the tested universe).

Returns
-------
float
    The enrichment p-value.
"""
    # Number of significant genes in the gene set
    k = len(significant_genes.intersection(gene_set))
    
    # Total number of significant genes
    n = len(significant_genes)
    
    # Total number of genes in the gene set
    K = len(gene_set.intersection(background_genes))
    
    # Total number of background genes
    N = len(background_genes)
    
    # Calculate the p-value using the hypergeometric distribution
    p_value = hypergeom.sf(k - 1, N, K, n)
    
    return p_value
def run_pathway_enrichment_analysis(significant_genes: set[str], gene_sets: dict[str, set[str]], background_genes: set[str]) -> pd.DataFrame:
    """
Run pathway enrichment analysis for a set of significant genes against multiple gene sets.

Parameters
----------
significant_genes : set[str]
    A set of significant genes.
gene_sets : dict[str, set[str]]
    A dictionary where keys are gene set names and values are sets of genes in those gene sets.
background_genes : set[str]
    A set of all background genes (the tested universe).

Returns
-------
pd.DataFrame
    A DataFrame containing the gene set names and their corresponding enrichment p-values.
"""
    p_values = {}
    overlap_counts = {}
    for pathway_name, gene_set in gene_sets.items():
        p_values[pathway_name] = compute_enrichment_pvalue(significant_genes, gene_set, background_genes)
        overlap_counts[pathway_name] = len(significant_genes.intersection(gene_set))
    pvalues_series = pd.Series(p_values)
    adjusted = compute_adjusted_pvalues(pvalues_series)
    results_df = pd.DataFrame({
        'Pathway': list(gene_sets.keys()),
        'P-value': list(p_values.values()),
        'Adjusted P-value': list(adjusted),
        'Overlap Count': list(overlap_counts.values())
    })
    return results_df


