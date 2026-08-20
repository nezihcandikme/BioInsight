"""
Stress test, regime 2: batch effects. Does method concordance (DESeq2 vs
edgeR vs limma-voom) change once there's a real technical covariate that
either does or doesn't get modeled?

Runs no DE method itself. Loads whatever `run_stress_batch.R` already
produced for the "naive" (~condition, ignoring the real per-sample
library-size batch split) and "adjusted" (~depth_batch + condition,
modeling it) design and hands each of the three tool pairs, under each
model, to `deconcord.concordance.methods.compute_de_concordance` -- the
same function `method_concordance.py` and `analyze_stress_unbalanced.py`
already use.

Usage (from the repo root, after `Rscript benchmarks/run_stress_batch.R`):
    python benchmarks/analyze_stress_batch.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from deconcord.concordance.methods import compute_de_concordance

RESULTS_DIR = Path(__file__).resolve().parent / "results" / "stress_batch"

METHODS = {
    "DESeq2": {"file_suffix": "deseq2_results.csv", "lfc_col": "log2FoldChange", "pvalue_col": "padj"},
    "edgeR": {"file_suffix": "edger_results.csv", "lfc_col": "logFC", "pvalue_col": "FDR"},
    "limma-voom": {"file_suffix": "limma_voom_results.csv", "lfc_col": "logFC", "pvalue_col": "adj.P.Val"},
}
PAIRS = [("DESeq2", "edgeR"), ("DESeq2", "limma-voom"), ("edgeR", "limma-voom")]
MODELS = ["naive", "adjusted"]


def _require(path: Path) -> Path:
    if not path.exists():
        raise SystemExit(
            f"Missing {path}. Run `Rscript benchmarks/run_stress_batch.R` first "
            "(which itself needs benchmarks/data/bottomly_counts.csv / "
            "bottomly_metadata.csv from run_deseq2_edger_bottomly.R)."
        )
    return path


def _batch_condition_table() -> pd.DataFrame:
    meta_path = _require(RESULTS_DIR / "batch_metadata.csv")
    return pd.read_csv(meta_path)


def _run_pair(model: str, name_a: str, name_b: str) -> dict:
    path_a = _require(RESULTS_DIR / f"{model}_{METHODS[name_a]['file_suffix']}")
    path_b = _require(RESULTS_DIR / f"{model}_{METHODS[name_b]['file_suffix']}")
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
    md_path = RESULTS_DIR.parent / "stress_batch_summary.md"
    md_parts = [
        "# Stress test: batch effects\n\n"
        "Real numbers from a real run: DESeq2, edgeR, and limma-voom, each with "
        "its own default settings, compared pairwise with "
        "`deconcord.concordance.methods.compute_de_concordance` under two design "
        "formulas -- `naive` (~condition, ignoring the real per-sample "
        "library-size batch split) and `adjusted` (~depth_batch + condition, "
        "modeling it). See `benchmarks/run_stress_batch.R` for exactly how the "
        "real, data-derived batch split was constructed and why.\n\n"
    ]

    batch_meta = _batch_condition_table()
    contingency = pd.crosstab(batch_meta["depth_batch"], batch_meta["condition"])
    print("Real batch x condition contingency table (from run_stress_batch.R's actual split):")
    print(contingency.to_string())
    md_parts.append("## Real batch x condition split\n\n```\n" + contingency.to_string() + "\n```\n\n")

    rows = []
    for model in MODELS:
        print(f"\n===== Model: {model} =====")
        md_parts.append(f"## Model: {model}\n\n")

        for name_a, name_b in PAIRS:
            result = _run_pair(model, name_a, name_b)
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
                "model": model, "pair": f"{name_a} vs {name_b}",
                "genes_compared": summary["genes_compared"],
                "jaccard_index": summary["jaccard_index"],
                "pearson_r_lfc": summary["pearson_r_lfc"],
                "directional_agreement": summary["directional_agreement"],
                "discordant_genes": len(result["discordant_genes"]),
            })

    summary_df = pd.DataFrame(rows)
    print("\n===== naive vs adjusted, side by side =====")
    print(summary_df.to_string(index=False))

    md_parts.append("## naive vs adjusted, side by side\n\n```\n")
    md_parts.append(summary_df.to_string(index=False))
    md_parts.append("\n```\n")

    with open(md_path, "w") as f:
        f.write("".join(md_parts))

    print(f"\nWrote {md_path}")


if __name__ == "__main__":
    main()
