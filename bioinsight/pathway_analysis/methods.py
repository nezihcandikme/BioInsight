from scipy.stats import hypergeom
import pandas as pd

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