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


def test_compute_enrichment_pvalue_empty_background_raises():
    with pytest.raises(ValueError, match="background_genes"):
        compute_enrichment_pvalue({"gene1"}, {"gene1"}, set())


def test_compute_enrichment_pvalue_empty_gene_set_is_never_enriched():
    background = {f"gene{i}" for i in range(1, 101)}
    significant = {f"gene{i}" for i in range(1, 11)}

    p_value = compute_enrichment_pvalue(significant, set(), background)

    assert p_value == pytest.approx(1.0)


def test_compute_enrichment_pvalue_empty_significant_genes_is_never_enriched():
    background = {f"gene{i}" for i in range(1, 101)}
    gene_set = {f"gene{i}" for i in range(1, 11)}

    p_value = compute_enrichment_pvalue(set(), gene_set, background)

    assert p_value == pytest.approx(1.0)


def test_run_pathway_enrichment_analysis():
    background = {f"gene{i}" for i in range(1, 101)}
    significant = {f"gene{i}" for i in range(1, 11)}

    gene_sets = {
        "pathway_A": {f"gene{i}" for i in range(1, 11)},
        "pathway_B": {f"gene{i}" for i in range(91, 101)},
    }

    results_df = run_pathway_enrichment_analysis(significant, gene_sets, background)
    results_df = results_df.set_index("pathway")

    assert results_df.loc["pathway_A", "p_value"] < 0.01
    assert results_df.loc["pathway_B", "p_value"] == pytest.approx(1.0)


def test_run_pathway_enrichment_analysis_overlap_count_ignores_non_background_genes():
    background = {f"gene{i}" for i in range(1, 21)}
    significant = {f"gene{i}" for i in range(1, 6)} | {"outside_background"}

    gene_sets = {
        "pathway_A": {f"gene{i}" for i in range(1, 6)} | {"also_outside_background"},
    }

    results_df = run_pathway_enrichment_analysis(significant, gene_sets, background)

    assert results_df.loc[results_df["pathway"] == "pathway_A", "overlap_count"].iloc[0] == 5


def test_run_pathway_enrichment_analysis_empty_background_raises():
    with pytest.raises(ValueError, match="background_genes"):
        run_pathway_enrichment_analysis({"gene1"}, {"pathway_A": {"gene1"}}, set())


def test_run_pathway_enrichment_analysis_empty_gene_sets_returns_empty_table():
    background = {f"gene{i}" for i in range(1, 11)}
    significant = {"gene1", "gene2"}

    results_df = run_pathway_enrichment_analysis(significant, {}, background)

    assert results_df.empty
    assert list(results_df.columns) == ["pathway", "p_value", "adjusted_p_value", "overlap_count"]
