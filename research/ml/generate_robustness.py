from pathlib import Path

import pandas as pd
import deconcord as dc


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "benchmarks" / "data"
OUTPUT_DIR = ROOT / "research" / "ml" / "generated"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def run_airway():
    counts = pd.read_csv(DATA_DIR / "airway_counts.csv", index_col=0)

    group_1 = [
        "SRR1039509",
        "SRR1039513",
        "SRR1039517",
        "SRR1039521",
    ]

    group_2 = [
        "SRR1039508",
        "SRR1039512",
        "SRR1039516",
        "SRR1039520",
    ]

    cpm = dc.compute_cpm(counts)
    log2_cpm = dc.log2_transform(cpm)

    result = dc.compute_resampling_stability(
        log2_cpm,
        group_1,
        group_2,
        resample_method="leave_one_out",
    )

    gene_stability = result["gene_stability"].reset_index()
    gene_stability.to_csv(
        OUTPUT_DIR / "airway_robustness.csv",
        index=False,
    )


def run_pasilla():
    counts = pd.read_csv(DATA_DIR / "pasilla_counts.csv", index_col=0)

    group_1 = [
        "treated1",
        "treated2",
        "treated3",
    ]

    group_2 = [
        "untreated1",
        "untreated2",
        "untreated3",
        "untreated4",
    ]

    cpm = dc.compute_cpm(counts)
    log2_cpm = dc.log2_transform(cpm)

    result = dc.compute_resampling_stability(
        log2_cpm,
        group_1,
        group_2,
        resample_method="leave_one_out",
    )

    gene_stability = result["gene_stability"].reset_index()
    gene_stability.to_csv(
        OUTPUT_DIR / "pasilla_robustness.csv",
        index=False,
    )


if __name__ == "__main__":
    run_airway()
    run_pasilla()

    print("Robustness tables written successfully.")