from unittest.mock import patch

import pandas as pd
import pytest

from omicforge.cli import main


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


def test_cli_moderated_method(tmp_path):
    out_dir = tmp_path / "out"

    exit_code = main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "sample1,sample2",
        "--group2", "sample3,sample4",
        "--method", "moderated",
        "--no-plots",
        "--out", str(out_dir),
    ])

    assert exit_code == 0
    de = pd.read_csv(out_dir / "differential_expression.csv", index_col=0)
    assert list(de.index) == ["gene1", "gene2", "gene3", "gene4"]


def test_cli_invalid_method_exits_nonzero(tmp_path):
    # argparse itself rejects an unknown --method choice via SystemExit,
    # before main()'s own try/except ever runs.
    with pytest.raises(SystemExit) as exc_info:
        main([
            "tests/fixtures/cli_counts.csv",
            "--group1", "sample1,sample2",
            "--group2", "sample3,sample4",
            "--method", "not-a-real-method",
            "--out", str(tmp_path / "out"),
        ])

    assert exc_info.value.code != 0


def test_cli_annotation_adds_gene_symbol_column(tmp_path):
    out_dir = tmp_path / "out"

    exit_code = main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "sample1,sample2",
        "--group2", "sample3,sample4",
        "--annotation", "tests/fixtures/cli_gene_annotation.csv",
        "--no-plots",
        "--out", str(out_dir),
    ])

    assert exit_code == 0
    de = pd.read_csv(out_dir / "differential_expression.csv", index_col=0)
    assert de.loc["gene1", "gene_symbol"] == "GENE_ONE"
    assert de.loc["gene2", "gene_symbol"] == "GENE_TWO"
    assert pd.isna(de.loc["gene3", "gene_symbol"])


def test_cli_without_annotation_has_no_symbol_column(tmp_path):
    out_dir = tmp_path / "out"

    main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "sample1,sample2",
        "--group2", "sample3,sample4",
        "--no-plots",
        "--out", str(out_dir),
    ])

    de = pd.read_csv(out_dir / "differential_expression.csv", index_col=0)
    assert "gene_symbol" not in de.columns


def test_cli_bad_annotation_file_exits_nonzero(tmp_path):
    bad_annotation = tmp_path / "bad_annotation.csv"
    bad_annotation.write_text("gene_id,gene_symbol\ngene1,GENE_ONE\ngene1,DIFFERENT\n")

    exit_code = main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "sample1,sample2",
        "--group2", "sample3,sample4",
        "--annotation", str(bad_annotation),
        "--out", str(tmp_path / "out"),
    ])

    assert exit_code == 1


def test_cli_missing_file_exits_nonzero():
    exit_code = main([
        "tests/fixtures/does_not_exist.csv",
        "--group1", "sample1",
        "--group2", "sample2",
    ])

    assert exit_code == 1


def test_cli_live_enrichment_organism_writes_output(tmp_path):
    out_dir = tmp_path / "out"
    fake_result = pd.DataFrame({"source": ["GO:BP"], "name": ["fake term"], "p_value": [0.01]})

    with patch(
        "omicforge.pathway_analysis.live_enrichment.run_gprofiler_enrichment",
        return_value=fake_result,
    ) as mock_query:
        exit_code = main([
            "tests/fixtures/cli_counts.csv",
            "--group1", "sample1,sample2",
            "--group2", "sample3,sample4",
            "--live-enrichment-organism", "hsapiens",
            "--no-plots",
            "--out", str(out_dir),
        ])

    assert exit_code == 0
    mock_query.assert_called_once()
    assert (out_dir / "live_pathway_enrichment.csv").exists()


def test_cli_without_live_enrichment_organism_skips_it(tmp_path):
    out_dir = tmp_path / "out"

    with patch("omicforge.pathway_analysis.live_enrichment.run_gprofiler_enrichment") as mock_query:
        main([
            "tests/fixtures/cli_counts.csv",
            "--group1", "sample1,sample2",
            "--group2", "sample3,sample4",
            "--no-plots",
            "--out", str(out_dir),
        ])

    mock_query.assert_not_called()
    assert not (out_dir / "live_pathway_enrichment.csv").exists()


def test_cli_live_enrichment_network_failure_exits_nonzero(tmp_path):
    with patch(
        "omicforge.pathway_analysis.live_enrichment.run_gprofiler_enrichment",
        side_effect=ConnectionError("Could not reach g:Profiler"),
    ):
        exit_code = main([
            "tests/fixtures/cli_counts.csv",
            "--group1", "sample1,sample2",
            "--group2", "sample3,sample4",
            "--live-enrichment-organism", "hsapiens",
            "--no-plots",
            "--out", str(tmp_path / "out"),
        ])

    assert exit_code == 1


def test_cli_empty_group_raises_via_exit_code():
    exit_code = main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "",
        "--group2", "sample3,sample4",
    ])

    assert exit_code == 1
