"""
Compare DEConcord's differential expression output against DESeq2 and
edgeR on the same input.

Supports two datasets (see DATASETS below): `airway` (human, the original
benchmark) and `pasilla` (fly, added later specifically to check whether
the pattern found on `airway` -- precision 1.0, recall ~0.08-0.10 against
each tool's own defaults -- is a property of the method or just a property
of that one dataset). Run the matching R script first
(`run_deseq2_edger.R` or `run_deseq2_edger_pasilla.R`) -- this script reads
what that produces (benchmarks/data/<dataset>_counts.csv,
benchmarks/data/<dataset>_metadata.csv, and each tool's results CSV) and
runs DEConcord on the identical count matrix and groups, then reports
where the three agree and disagree.

This is not trying to make DEConcord look good -- it's trying to find out
where a simple Welch's-t-test-on-log2-CPM approach actually diverges from
count-based negative-binomial modeling, and by how much. Both kinds of
result are worth knowing.

Usage (from the repo root, after the matching R script has run):
    python benchmarks/compare_results.py
    python benchmarks/compare_results.py --dataset pasilla
    python benchmarks/compare_results.py --dataset pasilla --method moderated
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from deconcord.pipeline import run_analysis

DATA_DIR = Path(__file__).resolve().parent / "data"
RESULTS_DIR = Path(__file__).resolve().parent / "results"

# Per-dataset config. `result_prefix` controls output filenames -- "airway"
# keeps the original, unprefixed names (deseq2_results.csv,
# lfc_vs_deseq2.png, ...) so the files already committed from the first
# benchmark run stay untouched; any dataset added after it gets its name
# prefixed instead, to avoid the two runs' outputs colliding.
DATASETS = {
    "airway": {
        "r_script": "run_deseq2_edger.R",
        "group_col": "dex",
        "group_1_value": "trt",
        "group_2_value": "untrt",
        "result_prefix": "",
    },
    "pasilla": {
        "r_script": "run_deseq2_edger_pasilla.R",
        "group_col": "condition",
        "group_1_value": "treated",
        "group_2_value": "untreated",
        "result_prefix": "pasilla_",
    },
}


def _require(path: Path, hint: str) -> Path:
    if not path.exists():
        raise SystemExit(f"Missing {path}. {hint}")
    return path


def load_inputs(dataset: str):
    cfg = DATASETS[dataset]
    hint = f"Run `Rscript benchmarks/{cfg['r_script']}` first."
    counts_path = _require(DATA_DIR / f"{dataset}_counts.csv", hint)
    meta_path = _require(DATA_DIR / f"{dataset}_metadata.csv", hint)
    deseq2_path = _require(RESULTS_DIR / f"{cfg['result_prefix']}deseq2_results.csv", hint)
    edger_path = _require(RESULTS_DIR / f"{cfg['result_prefix']}edger_results.csv", hint)

    counts = pd.read_csv(counts_path, index_col=0)
    meta = pd.read_csv(meta_path)
    deseq2 = pd.read_csv(deseq2_path, index_col="gene_id")
    edger = pd.read_csv(edger_path, index_col="gene_id")
    return counts, meta, deseq2, edger


def run_deconcord(counts: pd.DataFrame, meta: pd.DataFrame, dataset: str, method: str = "welch") -> pd.DataFrame:
    cfg = DATASETS[dataset]
    # DEConcord's log_fold_change is mean(group_1) - mean(group_2), so
    # group_1 has to be the treated samples to match each tool's own
    # "treated vs untreated"-direction contrast. Get this backwards and
    # every correlation below comes out negative for no statistical reason
    # at all.
    group_1 = meta.loc[meta[cfg["group_col"]] == cfg["group_1_value"], "sample"].tolist()
    group_2 = meta.loc[meta[cfg["group_col"]] == cfg["group_2_value"], "sample"].tolist()

    results = run_analysis(
        counts,
        group_1=group_1,
        group_2=group_2,
        min_count=10,
        min_samples=4,  # at least one full group's worth of samples
        method=method,
        generate_plots=False,
    )
    de = results["differential_expression"]
    print(f"DEConcord ({method}): {len(de)} genes tested (post min-count filtering), "
          f"{int(de['significant'].sum())} significant (adjusted p < 0.05, |log2FC| > 1)")
    return de


def compare(deconcord: pd.DataFrame, other: pd.DataFrame, other_name: str,
            other_lfc_col: str, other_padj_col: str, file_suffix: str = ""):
    merged = deconcord.join(other, how="inner", rsuffix=f"_{other_name}")
    merged = merged.dropna(subset=["log_fold_change", other_lfc_col])

    pearson_r, _ = pearsonr(merged["log_fold_change"], merged[other_lfc_col])
    spearman_r, _ = spearmanr(merged["log_fold_change"], merged[other_lfc_col])

    bi_sig = set(merged.index[merged["significant"]])
    other_sig = set(merged.index[merged[other_padj_col] < 0.05])
    intersection = bi_sig & other_sig
    union = bi_sig | other_sig
    jaccard = len(intersection) / len(union) if union else float("nan")
    precision = len(intersection) / len(bi_sig) if bi_sig else float("nan")
    recall = len(intersection) / len(other_sig) if other_sig else float("nan")

    summary = {
        "comparison": f"DEConcord vs {other_name}",
        "genes_in_common": len(merged),
        "pearson_r_lfc": round(pearson_r, 3),
        "spearman_r_lfc": round(spearman_r, 3),
        "deconcord_significant": len(bi_sig),
        "other_significant": len(other_sig),
        "significant_in_both": len(intersection),
        "jaccard_index": round(jaccard, 3),
        "precision_vs_other": round(precision, 3),
        "recall_vs_other": round(recall, 3),
    }

    fig, ax = plt.subplots(figsize=(7, 6))
    colors = merged.apply(
        lambda row: "#3654ff" if row.name in intersection
        else ("#ff5c5c" if row.name in bi_sig else ("#86868b" if row.name in other_sig else "#d2d2d7")),
        axis=1,
    )
    ax.scatter(merged[other_lfc_col], merged["log_fold_change"], c=colors, alpha=0.5, s=14)
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, color="black", linestyle="--", linewidth=1, label="y = x")
    ax.set_xlabel(f"{other_name} log2 fold change")
    ax.set_ylabel("DEConcord log fold change")
    ax.set_title(f"DEConcord vs {other_name} (Pearson r = {pearson_r:.3f})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"lfc_vs_{other_name.lower()}{file_suffix}.png", dpi=150)
    plt.close(fig)

    return summary


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset", choices=sorted(DATASETS), default="airway",
        help="Which benchmark dataset to compare on. Default 'airway' -- matches the "
             "original run and its committed, unprefixed output filenames. Other "
             "datasets write to suffixed filenames (e.g. lfc_vs_deseq2_pasilla.png) "
             "so they don't collide with the airway run's committed results.",
    )
    parser.add_argument(
        "--method", choices=["welch", "moderated"], default="welch",
        help="Which DEConcord DE method to benchmark. Default 'welch' -- matches the "
             "original run and its committed output filenames. 'moderated' writes to "
             "suffixed filenames (e.g. lfc_vs_deseq2_moderated.png) so it doesn't "
             "overwrite the welch results.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    dataset_suffix = "" if args.dataset == "airway" else f"_{args.dataset}"
    method_suffix = "" if args.method == "welch" else f"_{args.method}"
    file_suffix = dataset_suffix + method_suffix

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    counts, meta, deseq2, edger = load_inputs(args.dataset)
    deconcord = run_deconcord(counts, meta, args.dataset, method=args.method)

    deseq2_summary = compare(deconcord, deseq2, "DESeq2", "log2FoldChange", "padj", file_suffix)
    edger_summary = compare(deconcord, edger, "edgeR", "logFC", "FDR", file_suffix)

    summary_df = pd.DataFrame([deseq2_summary, edger_summary])
    print()
    print(summary_df.to_string(index=False))

    md_path = RESULTS_DIR / f"comparison_summary{file_suffix}.md"
    with open(md_path, "w") as f:
        f.write(f"# DEConcord ({args.method}) vs DESeq2 / edgeR — {args.dataset} dataset\n\n")
        f.write("Real numbers from a real run. See `benchmarks/README.md` for methodology "
                "and honest caveats before reading anything into these on their own.\n\n")
        # Plain code block instead of a real markdown table -- avoids pulling
        # in `tabulate` (pandas.to_markdown's dependency) just to format text
        # that's only ever read as a monospace block anyway.
        f.write("```\n")
        f.write(summary_df.to_string(index=False))
        f.write("\n```\n")

    print(f"\nWrote {md_path}")
    print(f"Wrote {RESULTS_DIR / f'lfc_vs_deseq2{file_suffix}.png'}")
    print(f"Wrote {RESULTS_DIR / f'lfc_vs_edger{file_suffix}.png'}")


if __name__ == "__main__":
    main()
