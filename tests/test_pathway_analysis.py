import pandas as pd
from scipy.stats import hypergeom
import pytest
from bioinsight.pathway_analysis.methods import compute_enrichment_pvalue, run_pathway_enrichment_analysis



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

def test_run_pathway_enrichment_analysis():
    background = {f"gene{i}" for i in range(1, 101)}
    significant = {f"gene{i}" for i in range(1, 11)}

    gene_sets = {
        "pathway_A": {f"gene{i}" for i in range(1, 11)},   # tam örtüşme
        "pathway_B": {f"gene{i}" for i in range(91, 101)},  # hiç örtüşme yok
    }

    results_df = run_pathway_enrichment_analysis(significant, gene_sets, background)

    assert results_df.loc["pathway_A" if "Pathway" not in results_df.columns else 0]

def test_run_pathway_enrichment_analysis():
    background = {f"gene{i}" for i in range(1, 101)}
    significant = {f"gene{i}" for i in range(1, 11)}

    gene_sets = {
        "pathway_A": {f"gene{i}" for i in range(1, 11)},
        "pathway_B": {f"gene{i}" for i in range(91, 101)},
    }

    results_df = run_pathway_enrichment_analysis(significant, gene_sets, background)
    results_df = results_df.set_index("Pathway")

    assert results_df.loc["pathway_A", "P-value"] < 0.01
    assert results_df.loc["pathway_B", "P-value"] == pytest.approx(1.0)