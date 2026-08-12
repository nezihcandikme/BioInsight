"""
Compare BioInsight's differential expression output against DESeq2 and
edgeR on the same input (the `airway` dataset).

Run `Rscript benchmarks/run_deseq2_edger.R` first -- this script reads
what that produces (benchmarks/data/airway_counts.csv,
benchmarks/data/airway_metadata.csv, benchmarks/results/deseq2_results.csv,
benchmarks/results/edger_results.csv) and runs BioInsight on the identical
count matrix and groups, then reports where the three agree and disagree.

This is not trying to make BioInsight look good -- it's trying to find out
where a simple Welch's-t-test-on-log2-CPM approach actually diverges from
count-based negative-binomial modeling, and by how much. Both kinds of
result are worth knowing.

Usage (from the repo root, after the R script has run):
    python benchmarks/compare_results.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bioinsight.pipeline import run_analysis

DATA_DIR = Path(__file__).resolve().parent / "data"
RESULTS_DIR = Path(__file__).resolve().parent / "results"


def _require(path: Path, hint: str) -> Path:
    if not path.exists():
        raise SystemExit(f"Missing {path}. {hint}")
    return path


def load_inputs():
    counts_path = _require(DATA_DIR / "airway_counts.csv", "Run `Rscript benchmarks/run_deseq2_edger.R` first.")
    meta_path = _require(DATA_DIR / "airway_metadata.csv", "Run `Rscript benchmarks/run_deseq2_edger.R` first.")
    deseq2_path = _require(RESULTS_DIR / "deseq2_results.csv", "Run `Rscript benchmarks/run_deseq2_edger.R` first.")
    edger_path = _require(RESULTS_DIR / "edger_results.csv", "Run `Rscript benchmarks/run_deseq2_edger.R` first.")

    counts = pd.read_csv(counts_path, index_col=0)
    meta = pd.read_csv(meta_path)
    deseq2 = pd.read_csv(deseq2_path, index_col="gene_id")
    edger = pd.read_csv(edger_path, index_col="gene_id")
    return counts, meta, deseq2, edger


def run_bioinsight(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    # BioInsight's log_fold_change is mean(group_1) - mean(group_2), so
    # group_1 has to be the treated samples to match DESeq2's
    # contrast=c("dex","trt","untrt") sign convention (positive = higher
    # in treatment). Get this backwards and every correlation below comes
    # out negative for no statistical reason at all.
    group_1 = meta.loc[meta["dex"] == "trt", "sample"].tolist()
    group_2 = meta.loc[meta["dex"] == "untrt", "sample"].tolist()

    results = run_analysis(
        counts,
        group_1=group_1,
        group_2=group_2,
        min_count=10,
        min_samples=4,  # at least one full group's worth of samples
        generate_plots=False,
    )
    de = results["differential_expression"]
    print(f"BioInsight: {len(de)} genes tested (post min-count filtering), "
          f"{int(de['significant'].sum())} significant (adjusted p < 0.05, |log2FC| > 1)")
    return de


def compare(bioinsight: pd.DataFrame, other: pd.DataFrame, other_name: str,
            other_lfc_col: str, other_padj_col: str):
    merged = bioinsight.join(other, how="inner", rsuffix=f"_{other_name}")
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
        "comparison": f"BioInsight vs {other_name}",
        "genes_in_common": len(merged),
        "pearson_r_lfc": round(pearson_r, 3),
        "spearman_r_lfc": round(spearman_r, 3),
        "bioinsight_significant": len(bi_sig),
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
    ax.set_ylabel("BioInsight log fold change")
    ax.set_title(f"BioInsight vs {other_name} (Pearson r = {pearson_r:.3f})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"lfc_vs_{other_name.lower()}.png", dpi=150)
    plt.close(fig)

    return summary


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    counts, meta, deseq2, edger = load_inputs()
    bioinsight = run_bioinsight(counts, meta)

    deseq2_summary = compare(bioinsight, deseq2, "DESeq2", "log2FoldChange", "padj")
    edger_summary = compare(bioinsight, edger, "edgeR", "logFC", "FDR")

    summary_df = pd.DataFrame([deseq2_summary, edger_summary])
    print()
    print(summary_df.to_string(index=False))

    md_path = RESULTS_DIR / "comparison_summary.md"
    with open(md_path, "w") as f:
        f.write("# BioInsight vs DESeq2 / edgeR — airway dataset\n\n")
        f.write("Real numbers from a real run. See `benchmarks/README.md` for methodology "
                "and honest caveats before reading anything into these on their own.\n\n")
        f.write(summary_df.to_markdown(index=False))
        f.write("\n")

    print(f"\nWrote {md_path}")
    print(f"Wrote {RESULTS_DIR / 'lfc_vs_deseq2.png'}")
    print(f"Wrote {RESULTS_DIR / 'lfc_vs_edger.png'}")


if __name__ == "__main__":
    main()
