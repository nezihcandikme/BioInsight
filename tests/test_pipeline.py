import pandas as pd
import pytest

from bioinsight.pipeline import run_analysis


def _toy_counts():
    return pd.DataFrame({
        "sample1": [10, 5, 2, 8],
        "sample2": [13, 7, 3, 11],
        "sample3": [200, 150, 220, 4],
        "sample4": [175, 142, 208, 7],
    }, index=["gene1", "gene2", "gene3", "gene4"])


def test_run_analysis_end_to_end():
    df = _toy_counts()

    results = run_analysis(
        df,
        group_1=["sample1", "sample2"],
        group_2=["sample3", "sample4"],
    )

    assert "qc" in results
    assert "normalized" in results
    assert "differential_expression" in results
    assert "volcano_fig" in results
    assert "pca_fig" in results

    de = results["differential_expression"]
    assert list(de.index) == list(df.index)
    assert "significant" in de.columns


def test_run_analysis_with_pathway_enrichment():
    df = _toy_counts()

    gene_sets = {
        "pathway_A": {"gene1", "gene2"},
        "pathway_B": {"gene3", "gene4"},
    }
    background_genes = set(df.index)

    results = run_analysis(
        df,
        group_1=["sample1", "sample2"],
        group_2=["sample3", "sample4"],
        gene_sets=gene_sets,
        background_genes=background_genes,
        generate_plots=False,
    )

    assert "pathway_enrichment" in results
    assert "volcano_fig" not in results


def test_run_analysis_missing_sample_raises():
    df = _toy_counts()

    with pytest.raises(ValueError, match="does_not_exist"):
        run_analysis(
            df,
            group_1=["sample1", "does_not_exist"],
            group_2=["sample3", "sample4"],
        )


def test_run_analysis_passes_through_de_thresholds():
    df = _toy_counts()

    strict = run_analysis(
        df,
        group_1=["sample1", "sample2"],
        group_2=["sample3", "sample4"],
        alpha=0.001,
        generate_plots=False,
    )
    lenient = run_analysis(
        df,
        group_1=["sample1", "sample2"],
        group_2=["sample3", "sample4"],
        alpha=1.0,
        lfc_threshold=0.0,
        generate_plots=False,
    )

    strict_sig = strict["differential_expression"]["significant"].sum()
    lenient_sig = lenient["differential_expression"]["significant"].sum()
    assert lenient_sig >= strict_sig


def test_run_analysis_min_count_filters_low_count_genes():
    df = _toy_counts()
    # gene4 = [8, 11, 4, 7] only clears a count of 10 in one sample
    # (sample2) -- below min_samples=2, so it should be the only gene
    # filtered out here.

    unfiltered = run_analysis(
        df,
        group_1=["sample1", "sample2"],
        group_2=["sample3", "sample4"],
        generate_plots=False,
    )
    filtered = run_analysis(
        df,
        group_1=["sample1", "sample2"],
        group_2=["sample3", "sample4"],
        min_count=10,
        min_samples=2,
        generate_plots=False,
    )

    assert "filtered_counts" not in unfiltered
    assert list(unfiltered["differential_expression"].index) == list(df.index)

    assert "filtered_counts" in filtered
    assert "gene4" not in filtered["filtered_counts"].index
    assert list(filtered["differential_expression"].index) == ["gene1", "gene2", "gene3"]


def test_run_analysis_duplicate_sample_across_groups_raises():
    df = _toy_counts()

    with pytest.raises(ValueError, match="sample1"):
        run_analysis(
            df,
            group_1=["sample1", "sample2"],
            group_2=["sample1", "sample4"],
        )


def test_run_analysis_gene_sets_without_background_raises():
    df = _toy_counts()

    with pytest.raises(ValueError):
        run_analysis(
            df,
            group_1=["sample1", "sample2"],
            group_2=["sample3", "sample4"],
            gene_sets={"pathway_A": {"gene1"}},
            generate_plots=False,
        )
