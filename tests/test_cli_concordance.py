"""`deconcord concordance` subcommand: compares two existing DE result table
CSVs with compute_de_concordance, without running the pipeline. See
tests/test_cli.py for the original `deconcord counts.csv ...` pipeline
invocation, unaffected by this subcommand's addition.
"""
import json

import pandas as pd
import pytest

from deconcord.cli import main

# Fixtures: concordance_a.csv / concordance_b.csv, both DEConcord's own
# default column names (log_fold_change, adjusted_p_value), 4 genes.
# At alpha=0.05: sig_a={g1,g2,g4}, sig_b={g1,g4} (g2's padj=0.5 in b).
# sig_in_both={g1,g4}, both same-signed in a and b -> concordant, no
# discordant genes. only_in_method_a={g2}, only_in_method_b={}.
# jaccard = |{g1,g4}| / |{g1,g2,g4}| = 2/3.


def test_concordance_runs_end_to_end_and_writes_outputs(tmp_path, capsys):
    out_dir = tmp_path / "out"

    exit_code = main([
        "concordance",
        "tests/fixtures/concordance_a.csv",
        "tests/fixtures/concordance_b.csv",
        "--out", str(out_dir),
    ])

    assert exit_code == 0
    assert (out_dir / "summary.json").exists()
    assert (out_dir / "merged.csv").exists()
    assert (out_dir / "concordant_genes.csv").exists()
    assert (out_dir / "discordant_genes.csv").exists()
    assert (out_dir / "only_in_method_a.csv").exists()
    assert (out_dir / "only_in_method_b.csv").exists()
    assert (out_dir / "run_metadata.json").exists()

    captured = capsys.readouterr()
    assert "genes compared" in captured.out


def test_concordance_summary_values(tmp_path):
    out_dir = tmp_path / "out"

    main([
        "concordance",
        "tests/fixtures/concordance_a.csv",
        "tests/fixtures/concordance_b.csv",
        "--out", str(out_dir),
    ])

    summary = json.loads((out_dir / "summary.json").read_text())
    assert summary["genes_compared"] == 4
    assert summary["significant_method_a"] == 3
    assert summary["significant_method_b"] == 2
    assert summary["significant_in_both"] == 2
    assert summary["jaccard_index"] == pytest.approx(2 / 3)


def test_concordance_gene_lists(tmp_path):
    out_dir = tmp_path / "out"

    main([
        "concordance",
        "tests/fixtures/concordance_a.csv",
        "tests/fixtures/concordance_b.csv",
        "--out", str(out_dir),
    ])

    concordant = pd.read_csv(out_dir / "concordant_genes.csv")
    discordant = pd.read_csv(out_dir / "discordant_genes.csv")
    only_a = pd.read_csv(out_dir / "only_in_method_a.csv")
    only_b = pd.read_csv(out_dir / "only_in_method_b.csv")

    assert list(concordant["gene_id"]) == ["g1", "g4"]
    assert len(discordant) == 0
    assert list(only_a["gene_id"]) == ["g2"]
    assert len(only_b) == 0


def test_concordance_merged_table_has_expected_shape(tmp_path):
    out_dir = tmp_path / "out"

    main([
        "concordance",
        "tests/fixtures/concordance_a.csv",
        "tests/fixtures/concordance_b.csv",
        "--out", str(out_dir),
    ])

    merged = pd.read_csv(out_dir / "merged.csv", index_col="gene_id")
    assert list(merged.index) == ["g1", "g2", "g3", "g4"]
    assert "log_fold_change_method_a" in merged.columns
    assert "significant_method_a" in merged.columns


def test_concordance_custom_names_change_output_filenames(tmp_path):
    out_dir = tmp_path / "out"

    exit_code = main([
        "concordance",
        "tests/fixtures/concordance_a.csv",
        "tests/fixtures/concordance_b.csv",
        "--name-a", "deseq2",
        "--name-b", "edger",
        "--out", str(out_dir),
    ])

    assert exit_code == 0
    assert (out_dir / "only_in_deseq2.csv").exists()
    assert (out_dir / "only_in_edger.csv").exists()

    summary = json.loads((out_dir / "summary.json").read_text())
    assert "significant_deseq2" in summary
    assert "significant_edger" in summary


def test_concordance_custom_alpha_changes_significance(tmp_path):
    out_dir = tmp_path / "out"

    # A stricter alpha than the default 0.05: g2's padj=0.02 in a no longer
    # clears it, so significant_method_a drops from 3 to 2.
    main([
        "concordance",
        "tests/fixtures/concordance_a.csv",
        "tests/fixtures/concordance_b.csv",
        "--alpha", "0.01",
        "--out", str(out_dir),
    ])

    summary = json.loads((out_dir / "summary.json").read_text())
    assert summary["significant_method_a"] == 2


def test_concordance_configurable_column_names(tmp_path):
    out_dir = tmp_path / "out"

    exit_code = main([
        "concordance",
        "tests/fixtures/concordance_a_deseq2_style.csv",
        "tests/fixtures/concordance_b.csv",
        "--lfc-col-a", "log2FoldChange",
        "--pvalue-col-a", "padj",
        "--out", str(out_dir),
    ])

    assert exit_code == 0
    summary = json.loads((out_dir / "summary.json").read_text())
    # Same underlying numbers as the default-column-name fixture, just
    # under DESeq2-style headers -- results should match exactly.
    assert summary["genes_compared"] == 4
    assert summary["jaccard_index"] == pytest.approx(2 / 3)


def test_concordance_explicit_gene_id_col(tmp_path):
    out_dir = tmp_path / "out"

    exit_code = main([
        "concordance",
        "tests/fixtures/concordance_a.csv",
        "tests/fixtures/concordance_b.csv",
        "--gene-id-col-a", "gene_id",
        "--gene-id-col-b", "gene_id",
        "--out", str(out_dir),
    ])

    assert exit_code == 0
    merged = pd.read_csv(out_dir / "merged.csv", index_col="gene_id")
    assert list(merged.index) == ["g1", "g2", "g3", "g4"]


def test_concordance_missing_file_exits_nonzero(tmp_path):
    exit_code = main([
        "concordance",
        "tests/fixtures/does_not_exist.csv",
        "tests/fixtures/concordance_b.csv",
        "--out", str(tmp_path / "out"),
    ])

    assert exit_code == 1


def test_concordance_missing_column_exits_nonzero(tmp_path, capsys):
    exit_code = main([
        "concordance",
        "tests/fixtures/concordance_a.csv",
        "tests/fixtures/concordance_b.csv",
        "--lfc-col-a", "not_a_real_column",
        "--out", str(tmp_path / "out"),
    ])

    assert exit_code == 1
    captured = capsys.readouterr()
    assert "not_a_real_column" in captured.err


def test_concordance_bad_gene_id_col_exits_nonzero(tmp_path):
    exit_code = main([
        "concordance",
        "tests/fixtures/concordance_a.csv",
        "tests/fixtures/concordance_b.csv",
        "--gene-id-col-a", "not_a_real_column",
        "--out", str(tmp_path / "out"),
    ])

    assert exit_code == 1


def test_concordance_invalid_alpha_exits_nonzero(tmp_path):
    exit_code = main([
        "concordance",
        "tests/fixtures/concordance_a.csv",
        "tests/fixtures/concordance_b.csv",
        "--alpha", "1.5",
        "--out", str(tmp_path / "out"),
    ])

    assert exit_code == 1


def test_concordance_run_metadata_fields(tmp_path):
    out_dir = tmp_path / "out"

    main([
        "concordance",
        "tests/fixtures/concordance_a.csv",
        "tests/fixtures/concordance_b.csv",
        "--out", str(out_dir),
    ])

    metadata = json.loads((out_dir / "run_metadata.json").read_text())
    assert "deconcord_version" in metadata
    assert metadata["input_summary"]["results_a_file"] == "tests/fixtures/concordance_a.csv"
    assert metadata["parameters"]["alpha"] == 0.05
    assert metadata["command"].startswith("deconcord concordance tests/fixtures/concordance_a.csv")


def test_bare_pipeline_invocation_still_works_unaffected(tmp_path):
    # Regression guard: adding the "concordance" subcommand must not touch
    # the original bare-invocation pipeline path (deconcord counts.csv
    # --group1 ... --group2 ...), which has no subcommand prefix at all.
    out_dir = tmp_path / "out"

    exit_code = main([
        "tests/fixtures/cli_counts.csv",
        "--group1", "sample1,sample2",
        "--group2", "sample3,sample4",
        "--no-plots",
        "--out", str(out_dir),
    ])

    assert exit_code == 0
    assert (out_dir / "tables" / "differential_expression.csv").exists()
