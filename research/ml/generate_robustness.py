"""Generate per-gene resampling-stability robustness tables (frac_significant,
direction_stability, rank_stability) for the ML research track's feature
tables. One dataset at a time or all of them; see parse_args() for CLI usage.
"""

import argparse
from pathlib import Path

import pandas as pd
import deconcord as dc


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "benchmarks" / "data"
OUTPUT_DIR = ROOT / "research" / "ml" / "generated"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def _run_resampling(dataset_id, counts_filename, group_1, group_2, output_filename, force):
    output_path = OUTPUT_DIR / output_filename

    if output_path.exists() and not force:
        print(f"[{dataset_id}] {output_path.name} already exists -- skipping (pass --force to recompute).")
        return

    print(f"[{dataset_id}] loading {counts_filename} ...")
    counts = pd.read_csv(DATA_DIR / counts_filename, index_col=0)

    print(f"[{dataset_id}] group_1 (n={len(group_1)}): {group_1}")
    print(f"[{dataset_id}] group_2 (n={len(group_2)}): {group_2}")

    cpm = dc.compute_cpm(counts)
    log2_cpm = dc.log2_transform(cpm)

    # Exact same call every dataset gets: leave-one-out resampling, every
    # other compute_resampling_stability parameter left at its package
    # default (alpha=0.05, lfc_threshold=1.0, de_method="welch",
    # stability_threshold=0.9). Not something to tune per dataset -- the
    # whole point of the frozen experiment design is that every dataset
    # gets measured the same way.
    print(f"[{dataset_id}] running compute_resampling_stability (leave_one_out) ...")
    result = dc.compute_resampling_stability(
        log2_cpm,
        group_1,
        group_2,
        resample_method="leave_one_out",
    )

    gene_stability = result["gene_stability"].reset_index()
    gene_stability.to_csv(output_path, index=False)
    print(f"[{dataset_id}] wrote {output_path} ({len(gene_stability)} genes).")


def run_airway(force: bool = False):
    _run_resampling(
        "airway",
        "airway_counts.csv",
        group_1=["SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521"],
        group_2=["SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520"],
        output_filename="airway_robustness.csv",
        force=force,
    )


def run_pasilla(force: bool = False):
    _run_resampling(
        "pasilla",
        "pasilla_counts.csv",
        group_1=["treated1", "treated2", "treated3"],
        group_2=["untreated1", "untreated2", "untreated3", "untreated4"],
        output_filename="pasilla_robustness.csv",
        force=force,
    )


def run_zebrafish(force: bool = False):
    # Corrected sample names -- the real zfGenes knockdown replicates are
    # Trt9/Trt11/Trt13, not Trt1/Trt3/Trt5 (a wrong guess from an earlier
    # pass, since this sandbox has no R to check the real column names
    # against). Controls (Ctl1/Ctl3/Ctl5) were already correct.
    _run_resampling(
        "zebrafish",
        "zebrafish_counts.csv",
        group_1=["Trt9", "Trt11", "Trt13"],
        group_2=["Ctl1", "Ctl3", "Ctl5"],
        output_filename="zebrafish_robustness.csv",
        force=force,
    )


DATASET_RUNNERS = {
    "airway": run_airway,
    "pasilla": run_pasilla,
    "zebrafish": run_zebrafish,
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset", choices=[*sorted(DATASET_RUNNERS), "all"], default="all",
        help="Which dataset to generate robustness features for. Default 'all' runs every dataset in "
             "DATASET_RUNNERS.",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Recompute even if the dataset's robustness CSV already exists. Default: skip datasets that "
             "already have one, since a leave-one-out resampling run is not free and this script is meant to "
             "be safe to rerun after adding a new dataset without redoing the existing ones.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    dataset_ids = sorted(DATASET_RUNNERS) if args.dataset == "all" else [args.dataset]

    print(f"Generating robustness tables for: {', '.join(dataset_ids)}")
    for dataset_id in dataset_ids:
        DATASET_RUNNERS[dataset_id](force=args.force)

    print("Done.")


if __name__ == "__main__":
    main()
