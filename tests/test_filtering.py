import pandas as pd
import pytest

from deconcord.filtering.methods import filter_low_count_genes


def test_filter_low_count_genes_keeps_well_expressed_genes():
    df = pd.DataFrame({
        "sample1": [50, 2, 0],
        "sample2": [40, 1, 0],
        "sample3": [60, 3, 1],
    }, index=["expressed_gene", "barely_expressed_gene", "dead_gene"])

    filtered = filter_low_count_genes(df, min_count=10, min_samples=2)

    assert list(filtered.index) == ["expressed_gene"]


def test_filter_low_count_genes_boundary_is_inclusive():
    # A gene hitting exactly min_count in exactly min_samples samples
    # should survive, not get treated as "not quite enough."
    df = pd.DataFrame({
        "sample1": [10, 9],
        "sample2": [10, 9],
        "sample3": [0, 9],
    }, index=["boundary_gene", "just_under_gene"])

    filtered = filter_low_count_genes(df, min_count=10, min_samples=2)

    assert list(filtered.index) == ["boundary_gene"]


def test_filter_low_count_genes_negative_min_count_raises():
    df = pd.DataFrame({"sample1": [10, 20]})

    with pytest.raises(ValueError, match="min_count"):
        filter_low_count_genes(df, min_count=-1)


def test_filter_low_count_genes_zero_min_samples_raises():
    df = pd.DataFrame({"sample1": [10, 20]})

    with pytest.raises(ValueError, match="min_samples"):
        filter_low_count_genes(df, min_samples=0)


def test_filter_low_count_genes_min_samples_exceeds_matrix_raises():
    df = pd.DataFrame({"sample1": [10, 20], "sample2": [10, 20]})

    with pytest.raises(ValueError, match="min_samples"):
        filter_low_count_genes(df, min_samples=3)


def test_filter_low_count_genes_all_genes_dropped_raises():
    df = pd.DataFrame({
        "sample1": [1, 2],
        "sample2": [1, 2],
    }, index=["gene1", "gene2"])

    with pytest.raises(ValueError, match="No genes pass"):
        filter_low_count_genes(df, min_count=1000, min_samples=1)


def test_filter_low_count_genes_preserves_columns():
    df = pd.DataFrame({
        "sample1": [50, 0],
        "sample2": [40, 0],
    }, index=["gene1", "gene2"])

    filtered = filter_low_count_genes(df, min_count=10, min_samples=1)

    assert list(filtered.columns) == ["sample1", "sample2"]
