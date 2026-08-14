import numpy as np
import pandas as pd
import pytest

from deconcord.annotation.methods import annotate_de_results, load_gene_annotation


def test_load_gene_annotation_parses_mapping():
    annotation = load_gene_annotation("tests/fixtures/sample_gene_annotation.csv")

    assert annotation == {
        "ENSG00000141510": "TP53",
        "ENSG00000146648": "EGFR",
        "ENSG00000012048": "BRCA1",
    }


def test_load_gene_annotation_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_gene_annotation("tests/fixtures/does_not_exist.csv")


def test_load_gene_annotation_missing_column_raises(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("gene_id,symbol\nENSG1,TP53\n")

    with pytest.raises(ValueError, match="gene_symbol"):
        load_gene_annotation(str(path))


def test_load_gene_annotation_empty_file_raises(tmp_path):
    path = tmp_path / "empty.csv"
    path.write_text("gene_id,gene_symbol\n")

    with pytest.raises(ValueError, match="no data rows"):
        load_gene_annotation(str(path))


def test_load_gene_annotation_missing_value_raises(tmp_path):
    path = tmp_path / "incomplete.csv"
    path.write_text("gene_id,gene_symbol\nENSG1,TP53\nENSG2,\n")

    with pytest.raises(ValueError, match="missing"):
        load_gene_annotation(str(path))


def test_load_gene_annotation_duplicate_same_symbol_is_fine(tmp_path):
    path = tmp_path / "dup_ok.csv"
    path.write_text("gene_id,gene_symbol\nENSG1,TP53\nENSG1,TP53\nENSG2,EGFR\n")

    annotation = load_gene_annotation(str(path))

    assert annotation == {"ENSG1": "TP53", "ENSG2": "EGFR"}


def test_load_gene_annotation_duplicate_conflicting_symbol_raises(tmp_path):
    path = tmp_path / "dup_conflict.csv"
    path.write_text("gene_id,gene_symbol\nENSG1,TP53\nENSG1,EGFR\n")

    with pytest.raises(ValueError, match="ENSG1"):
        load_gene_annotation(str(path))


def test_load_gene_annotation_custom_column_names(tmp_path):
    path = tmp_path / "custom.csv"
    path.write_text("ensembl_id,hgnc_symbol\nENSG1,TP53\n")

    annotation = load_gene_annotation(str(path), id_col="ensembl_id", symbol_col="hgnc_symbol")

    assert annotation == {"ENSG1": "TP53"}


def _de_results():
    return pd.DataFrame({
        "log_fold_change": [2.1, -1.4, 0.3],
        "p_value": [0.001, 0.02, 0.5],
        "adjusted_p_value": [0.003, 0.03, 0.6],
        "significant": [True, True, False],
    }, index=["ENSG00000141510", "ENSG00000146648", "ENSG_UNKNOWN"])


def test_annotate_de_results_adds_symbol_column():
    annotation = load_gene_annotation("tests/fixtures/sample_gene_annotation.csv")
    de = _de_results()

    annotated = annotate_de_results(de, annotation)

    assert list(annotated.columns)[0] == "gene_symbol"
    assert annotated.loc["ENSG00000141510", "gene_symbol"] == "TP53"
    assert annotated.loc["ENSG00000146648", "gene_symbol"] == "EGFR"


def test_annotate_de_results_unmapped_gene_is_nan():
    annotation = load_gene_annotation("tests/fixtures/sample_gene_annotation.csv")
    de = _de_results()

    annotated = annotate_de_results(de, annotation)

    assert pd.isna(annotated.loc["ENSG_UNKNOWN", "gene_symbol"])


def test_annotate_de_results_empty_annotation_no_crash():
    de = _de_results()

    annotated = annotate_de_results(de, {})

    assert annotated["gene_symbol"].isna().all()
    assert len(annotated) == len(de)


def test_annotate_de_results_does_not_mutate_input():
    annotation = load_gene_annotation("tests/fixtures/sample_gene_annotation.csv")
    de = _de_results()
    original_columns = list(de.columns)

    annotate_de_results(de, annotation)

    assert list(de.columns) == original_columns


def test_annotate_de_results_preserves_original_columns_and_index():
    annotation = load_gene_annotation("tests/fixtures/sample_gene_annotation.csv")
    de = _de_results()

    annotated = annotate_de_results(de, annotation)

    assert list(annotated.index) == list(de.index)
    for col in de.columns:
        assert col in annotated.columns
        pd.testing.assert_series_equal(annotated[col], de[col])


def test_annotate_de_results_on_empty_de_results_no_crash():
    empty_de = pd.DataFrame(columns=["log_fold_change", "p_value"], dtype=float)
    annotated = annotate_de_results(empty_de, {"ENSG1": "TP53"})

    assert "gene_symbol" in annotated.columns
    assert len(annotated) == 0


def test_load_gene_annotation_returns_np_nan_type_consistently():
    # Sanity check that unmapped values are real NaN (isna-detectable),
    # not the string "nan" or None slipping through inconsistently.
    de = pd.DataFrame({"log_fold_change": [1.0]}, index=["NOT_IN_ANNOTATION"])
    annotated = annotate_de_results(de, {"OTHER_GENE": "X"})

    assert np.isnan(annotated["gene_symbol"].iloc[0]) or pd.isna(annotated["gene_symbol"].iloc[0])
