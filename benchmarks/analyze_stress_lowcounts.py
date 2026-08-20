"""
Stress test, regime 3: low counts. Does method concordance (DESeq2 vs
edgeR vs limma-voom) hold up as real sequencing depth drops?

Runs no DE method itself. Loads whatever `run_stress_lowcounts.R` already
produced per depth level (`full` = the real, unmodified Bottomly counts;
`moderate`/`low` = the same real counts thinned to 25%/5% of their real
depth by binomial thinning -- see that script's header for why this is
real read removal, not fabricated signal) and hands each of the three
tool pairs, at each depth, to
`deconcord.concordance.methods.compute_de_concordance` -- the same
function `method_concordance.py`, `analyze_stress_unbalanced.py`, and
`analyze_stress_batch.py` already use.

Usage (from the repo root, after `Rscript benchmarks/run_stress_lowcounts.R`):
    python benchmarks/analyze_stress_lowcounts.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from deconcord.concordance.methods import compute_de_concordance

RESULTS_DIR = Path(__file__).resolve().parent / "results" / "stress_lowcounts"

METHODS = {
    "DESeq2": {"file_suffix": "deseq2_results.csv", "lfc_col": "log2FoldChange", "pvalue_col": "padj"},
    "edgeR": {"file_suffix": "edger_results.csv", "lfc_col": "logFC", "pvalue_col": "FDR"},
    "limma-voom": {"file_suffix": "limma_voom_results.csv", "lfc_col": "logFC", "pvalue_col": "adj.P.Val"},
}
PAIRS = [("DESeq2", "edgeR"), ("DESeq2", "limma-voom"), ("edgeR", "limma-voom")]

# Ordered so the summary table reads as an actual depth gradient (full
# real depth down to the shallowest), not alphabetically.
DEPTHS = ["full", "moderate", "low"]


def _require(path: Path) -> Path:
    if not path.exists():
        raise SystemExit(
            f"Missing {path}. Run `Rscript benchmarks/run_stress_lowcounts.R` first "
            "(which itself needs benchmarks/data/bottomly_counts.csv / "
            "bottomly_metadata.csv from run_deseq2_edger_bottomly.R)."
        )
    return path


def _total_reads(depth: str) -> int:
    counts_path = _require(RESULTS_DIR / f"{depth}_counts.csv")
    counts = pd.read_csv(counts_path, index_col="gene_id")
    return int(counts.to_numpy().sum())


def _run_pair(depth: str, name_a: str, name_b: str) -> dict:
    path_a = _require(RESULTS_DIR / f"{depth}_{METHODS[name_a]['file_suffix']}")
    path_b = _require(RESULTS_DIR / f"{depth}_{METHODS[name_b]['file_suffix']}")
    df_a = pd.read_csv(path_a, index_col="gene_id")
    df_b = pd.read_csv(path_b, index_col="gene_id")

    result = compute_de_concordance(
        df_a, df_b, name_a=name_a, name_b=name_b,
        lfc_col_a=METHODS[name_a]["lfc_col"], pvalue_col_a=METHODS[name_a]["pvalue_col"],
        lfc_col_b=METHODS[name_b]["lfc_col"], pvalue_col_b=METHODS[name_b]["pvalue_col"],
    )

    merged = result["merged"]
    sig_a = merged[f"significant_{name_a}"]
    sig_b = merged[f"significant_{name_b}"]
    sig_either = sig_a | sig_b
    lfc_a = merged[f"log_fold_change_{name_a}"]
    lfc_b = merged[f"log_fold_change_{name_b}"]
    same_direction_either = np.sign(lfc_a[sig_either]) == np.sign(lfc_b[sig_either])

    result["summary"]["significant_union"] = int(sig_either.sum())
    result["summary"]["directional_agreement_any_significant"] = (
        float(same_direction_either.mean()) if sig_either.any() else float("nan")
    )
    result["summary"]["opposite_direction_any_significant"] = int((~same_direction_either).sum())

    return result


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    md_path = RESULTS_DIR.parent / "stress_lowcounts_summary.md"
    md_parts = [
        "# Stress test: low counts\n\n"
        "Real numbers from a real run: DESeq2, edgeR, and limma-voom, each with "
        "its own default settings, compared pairwise with "
        "`deconcord.concordance.methods.compute_de_concordance` at three real "
        "sequencing-depth levels of the Bottomly dataset (`full` = real, "
        "unmodified depth; `moderate`/`low` = the same real reads binomially "
        "thinned to 25%/5%). See `benchmarks/run_stress_lowcounts.R` for exactly "
        "how thinning was done and why it doesn't fabricate any signal.\n\n"
    ]

    rows = []
    for depth in DEPTHS:
        total_reads = _total_reads(depth)
        print(f"\n===== Depth: {depth} ({total_reads:,} total reads) =====")
        md_parts.append(f"## Depth: {depth} ({total_reads:,} total reads)\n\n")

        for name_a, name_b in PAIRS:
            result = _run_pair(depth, name_a, name_b)
            summary = result["summary"]

            print(f"{name_a} vs {name_b}")
            for key, value in summary.items():
                print(f"  {key}: {value:.4f}" if isinstance(value, float) else f"  {key}: {value}")
            print(f"  concordant genes: {len(result['concordant_genes'])}")
            print(f"  discordant genes: {len(result['discordant_genes'])}")

            md_parts.append(f"**{name_a} vs {name_b}**\n\n```\n")
            for key, value in summary.items():
                md_parts.append(f"{key}: {value:.4f}\n" if isinstance(value, float) else f"{key}: {value}\n")
            md_parts.append(f"concordant_genes: {len(result['concordant_genes'])}\n")
            md_parts.append(f"discordant_genes: {len(result['discordant_genes'])}\n")
            md_parts.append("```\n\n")

            rows.append({
                "depth": depth, "total_reads": total_reads,
                "pair": f"{name_a} vs {name_b}",
                "genes_compared": summary["genes_compared"],
                "jaccard_index": summary["jaccard_index"],
                "pearson_r_lfc": summary["pearson_r_lfc"],
                "directional_agreement": summary["directional_agreement"],
                "discordant_genes": len(result["discordant_genes"]),
            })

    summary_df = pd.DataFrame(rows)
    print("\n===== Gradient across depth levels =====")
    print(summary_df.to_string(index=False))

    md_parts.append("## Gradient across depth levels\n\n```\n")
    md_parts.append(summary_df.to_string(index=False))
    md_parts.append("\n```\n")

    with open(md_path, "w") as f:
        f.write("".join(md_parts))

    print(f"\nWrote {md_path}")


if __name__ == "__main__":
    main()
