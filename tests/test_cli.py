import pandas as pd

from bioinsight.cli import main


def test_cli_runs_end_to_end_and_writes_outputs(tmp_path, capsys):
    out_dir = tmp_path / "out"

    exit_code = main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "sample1,sample2",
        "--group2", "sample3,sample4",
        "--out", str(out_dir),
    ])

    assert exit_code == 0
    assert (out_dir / "differential_expression.csv").exists()
    assert (out_dir / "qc.csv").exists()
    assert (out_dir / "volcano.png").exists()
    assert (out_dir / "pca.png").exists()

    de = pd.read_csv(out_dir / "differential_expression.csv", index_col=0)
    assert list(de.index) == ["gene1", "gene2", "gene3", "gene4"]

    captured = capsys.readouterr()
    assert "genes tested" in captured.out


def test_cli_no_plots_skips_figures(tmp_path):
    out_dir = tmp_path / "out"

    main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "sample1,sample2",
        "--group2", "sample3,sample4",
        "--no-plots",
        "--out", str(out_dir),
    ])

    assert not (out_dir / "volcano.png").exists()
    assert not (out_dir / "pca.png").exists()


def test_cli_with_gmt_writes_pathway_enrichment(tmp_path):
    out_dir = tmp_path / "out"

    main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "sample1,sample2",
        "--group2", "sample3,sample4",
        "--gmt", "tests/fixtures/sample_gene_sets.gmt",
        "--no-plots",
        "--out", str(out_dir),
    ])

    # sample_gene_sets.gmt doesn't share gene names with cli_counts.csv,
    # so this exercises the "runs without crashing on zero overlap" path --
    # the real assertion is that it still produced a table.
    assert (out_dir / "pathway_enrichment.csv").exists()


def test_cli_min_count_filters_low_count_genes(tmp_path):
    out_dir = tmp_path / "out"

    main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "sample1,sample2",
        "--group2", "sample3,sample4",
        "--min-count", "10",
        "--min-samples", "2",
        "--no-plots",
        "--out", str(out_dir),
    ])

    de = pd.read_csv(out_dir / "differential_expression.csv", index_col=0)
    # gene4 = [8, 11, 4, 7] only clears 10 reads in one sample -- filtered out.
    assert "gene4" not in list(de.index)


def test_cli_missing_sample_exits_nonzero(tmp_path, capsys):
    exit_code = main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "sample1,does_not_exist",
        "--group2", "sample3,sample4",
        "--out", str(tmp_path / "out"),
    ])

    assert exit_code == 1
    captured = capsys.readouterr()
    assert "does_not_exist" in captured.err


def test_cli_missing_file_exits_nonzero():
    exit_code = main([
        "tests/fixtures/does_not_exist.csv",
        "--group1", "sample1",
        "--group2", "sample2",
    ])

    assert exit_code == 1


def test_cli_empty_group_raises_via_exit_code():
    exit_code = main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "",
        "--group2", "sample3,sample4",
    ])

    assert exit_code == 1
