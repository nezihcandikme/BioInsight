from pathlib import Path

import pandas as pd

from feature_table import MethodSchema, build_feature_table


SCHEMAS = {
    "source": MethodSchema("id", "effect", "raw_p", "source_q"),
    "target": MethodSchema("id", "effect", "raw_p", "target_q"),
}


def _write(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def test_alignment_target_and_identity(tmp_path):
    source = tmp_path / "source.csv"
    target = tmp_path / "target.csv"
    _write(source, [
        {"id": "g1", "effect": 1.0, "raw_p": 0.001, "source_q": 0.02},
        {"id": "g2", "effect": -1.0, "raw_p": 0.2, "source_q": 0.3},
        {"id": "source_only", "effect": 0.0, "raw_p": 0.8, "source_q": 0.9},
    ])
    _write(target, [
        {"id": "g2", "effect": -0.8, "raw_p": 0.01, "target_q": 0.049},
        {"id": "g1", "effect": 0.7, "raw_p": 0.02, "target_q": 0.05},
        {"id": "target_only", "effect": 0.0, "raw_p": 0.9, "target_q": 0.9},
    ])

    table, report = build_feature_table(
        source, target, dataset_id="toy", source_method="source", target_method="target", schemas=SCHEMAS
    )

    assert table["gene_id"].tolist() == ["g1", "g2"]
    assert table["target_method_significant"].tolist() == [0, 1]
    assert set(table["dataset_id"]) == {"toy"}
    assert set(table["source_method"]) == {"source"}
    assert set(table["target_method"]) == {"target"}
    assert report["source_only_gene_count"] == 1
    assert report["target_only_gene_count"] == 1


def test_method_specific_adjusted_p_value_and_missing_drop(tmp_path):
    source = tmp_path / "source.csv"
    target = tmp_path / "target.csv"
    _write(source, [
        {"id": "kept", "effect": 1.0, "raw_p": 0.9, "source_q": 0.01},
        {"id": "missing_source", "effect": 1.0, "raw_p": 0.1, "source_q": None},
    ])
    _write(target, [
        {"id": "kept", "effect": 1.0, "raw_p": 0.0001, "target_q": 0.2},
        {"id": "missing_source", "effect": 1.0, "raw_p": 0.1, "target_q": 0.01},
    ])

    table, report = build_feature_table(
        source, target, dataset_id="toy", source_method="source", target_method="target", schemas=SCHEMAS
    )

    assert table.loc[0, "source_adjusted_p_value"] == 0.01
    assert table.loc[0, "target_method_significant"] == 0
    assert report["overlap_dropped_missing_required_count"] == 1

