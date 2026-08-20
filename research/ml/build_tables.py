from __future__ import annotations

import json
from pathlib import Path

from feature_table import build_feature_table


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent / "generated"

DATASETS = {
    "airway": (
        "deseq2_results.csv",
        "edger_results.csv",
        "airway_robustness.csv",
    ),
    "pasilla": (
        "pasilla_deseq2_results.csv",
        "pasilla_edger_results.csv",
        "pasilla_robustness.csv",
    ),
    "zebrafish": (
        "zebrafish_deseq2_results.csv",
        "zebrafish_edger_results.csv",
        "zebrafish_robustness.csv",
    ),
    "bottomly": (
        "bottomly_deseq2_results.csv",
        "bottomly_edger_results.csv",
        "bottomly_robustness.csv",
    ),
}


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)

    reports = {}

    for dataset_id, (source_name, target_name, robustness_name) in DATASETS.items():
        table, report = build_feature_table(
            ROOT / "benchmarks" / "results" / source_name,
            ROOT / "benchmarks" / "results" / target_name,
            dataset_id=dataset_id,
            robustness_path=OUTPUT / robustness_name,
        )

        table.to_csv(
            OUTPUT / f"{dataset_id}_deseq2_to_edger.csv",
            index=False,
        )

        reports[dataset_id] = report

    (OUTPUT / "build_report.json").write_text(
        json.dumps(reports, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()