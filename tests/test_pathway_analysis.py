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

def test_compute_enrichment_pvalue_ignores_genes_outside_background():
    background = {f"gene{i}" for i in range(1, 21)}
    gene_set = {f"gene{i}" for i in range(1, 11)} | {"not_in_background_1"}
    significant = {f"gene{i}" for i in range(1, 11)} | {"not_in_background_2"}

    with_extra = compute_enrichment_pvalue(significant, gene_set, background)
    without_extra = compute_enrichment_pvalue(
        {f"gene{i}" for i in range(1, 11)},
        {f"gene{i}" for i in range(1, 11)},
        background,
    )

    assert with_extra == pytest.approx(without_extra)


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
