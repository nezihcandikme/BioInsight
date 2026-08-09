import pandas as pd
from scipy.stats import hypergeom
import pytest
from bioinsight.pathway_analysis.methods import compute_enrichment_pvalue



def test_compute_enrichment_pvalue_no_overlap():
    background = {f"gene{i}" for i in range(1, 101)}
    gene_set = {f"gene{i}" for i in range(1, 11)}
    significant = {f"gene{i}" for i in range(91, 101)}

    p_value = compute_enrichment_pvalue(significant, gene_set, background)
    assert p_value == pytest.approx(1.0)




def test_compute_enrichment_pvalue_perfect_overlap():
    background = {f"gene{i}" for i in range(1, 101)}
    gene_set = {f"gene{i}" for i in range(1, 11)}
    significant = {f"gene{i}" for i in range(1, 11)}

    p_value = compute_enrichment_pvalue(significant, gene_set, background)
    assert p_value < 0.01