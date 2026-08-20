"""
Stress test, regime 1: unbalanced groups. Does method concordance (DESeq2 vs
edgeR vs limma-voom, the same three tools and the same
`compute_de_concordance` question `method_concordance.py` already checks on
airway/pasilla) hold up once group sizes stop being clean and balanced?

Runs no DE method itself. Loads whatever `run_stress_unbalanced.R` already
produced per imbalance level (`balanced` = 10 B6 vs 11 D2, the real full
Bottomly design; `moderate` = 5 vs 11; `severe` = 2 vs 11 -- see that
script's header for why these three and why D2 is held fixed) and hands
each of the three tool pairs, at each level, to
`deconcord.concordance.methods.compute_de_concordance` -- the same library
function `method_concordance.py` uses, applied across a harder axis than
the two clean, balanced datasets checked so far.

Usage (from the repo root, after `Rscript benchmarks/run_stress_unbalanced.R`):
    python benchmarks/analyze_stress_unbalanced.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from deconcord.concordance.methods import compute_de_concordance

RESULTS_DIR = Path(__file__).resolve().parent / "results" / "stress_unbalanced"

# Same three tools, same file-suffix/column convention as
# method_concordance.py -- deliberately not redefined differently here, a
# reader who already knows that script's METHODS dict shouldn't have to
# learn a second one.
METHODS = {
    "DESeq2": {"file_suffix": "deseq2_results.csv", "lfc_col": "log2FoldChange", "pvalue_col": "padj"},
    "edgeR": {"file_suffix": "edger_results.csv", "lfc_col": "logFC", "pvalue_col": "FDR"},
    "limma-voom": {"file_suffix": "limma_voom_results.csv", "lfc_col": "logFC", "pvalue_col": "adj.P.Val"},
}
PAIRS = [("DESeq2", "edgeR"), ("DESeq2", "limma-voom"), ("edgeR", "limma-voom")]

# Ordered so the summary table reads as an actual imbalance gradient, not
# alphabetically.
LEVELS = ["balanced", "moderate", "severe"]


def _require(path: Path) -> Path:
    if not path.exists():
        raise SystemExit(
            f"Missing {path}. Run `Rscript benchmarks/run_stress_unbalanced.R` first "
            "(which itself needs benchmarks/data/bottomly_counts.csv / "
            "bottomly_metadata.csv from run_deseq2_edger_bottomly.R)."
        )
    return path


def _level_group_sizes(level: str) -> tuple[int, int]:
    meta_path = _require(RESULTS_DIR / f"{level}_metadata.csv")
    meta = pd.read_csv(meta_path)
    n_b6 = int((meta["condition"] == "B6").sum())
    n_d2 = int((meta["condition"] == "D2").sum())
    return n_b6, n_d2


def _run_pair(level: str, name_a: str, name_b: str) -> dict:
    path_a = _require(RESULTS_DIR / f"{level}_{METHODS[name_a]['file_suffix']}")
    path_b = _require(RESULTS_DIR / f"{level}_{METHODS[name_b]['file_suffix']}")
    df_a = pd.read_csv(path_a, index_col="gene_id")
    df_b = pd.read_csv(path_b, index_col="gene_id")

    result = compute_de_concordance(
        df_a, df_b, name_a=name_a, name_b=name_b,
        lfc_col_a=METHODS[name_a]["lfc_col"], pvalue_col_a=METHODS[name_a]["pvalue_col"],
        lfc_col_b=METHODS[name_b]["lfc_col"], pvalue_col_b=METHODS[name_b]["pvalue_col"],
    )

    # Same widened-to-"either significant" directional check
    # method_concordance.py adds, kept consistent so the two scripts'
    # numbers are directly comparable rather than subtly different metrics
    # that happen to share a name.
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
    md_path = RESULTS_DIR.parent / "stress_unbalanced_summary.md"
    md_parts = [
        "# Stress test: unbalanced groups\n\n"
        "Real numbers from a real run: DESeq2, edgeR, and limma-voom, each with "
        "its own default settings, compared pairwise with "
        "`deconcord.concordance.methods.compute_de_concordance` at three real "
        "sample-subset levels of the Bottomly dataset. See "
        "`benchmarks/run_stress_unbalanced.R` for exactly which samples make up "
        "each level and why. Not tuned toward any particular outcome; this is "
        "the plain output of three established tools disagreeing or agreeing "
        "with each other as one group's sample size shrinks.\n\n"
    ]

    rows = []
    for level in LEVELS:
        n_b6, n_d2 = _level_group_sizes(level)
        print(f"\n===== Level: {level} ({n_b6} B6 vs {n_d2} D2) =====")
        md_parts.append(f"## Level: {level} ({n_b6} B6 vs {n_d2} D2)\n\n")

        for name_a, name_b in PAIRS:
            result = _run_pair(level, name_a, name_b)
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
                "level": level, "n_b6": n_b6, "n_d2": n_d2,
                "pair": f"{name_a} vs {name_b}",
                "genes_compared": summary["genes_compared"],
                "jaccard_index": summary["jaccard_index"],
                "pearson_r_lfc": summary["pearson_r_lfc"],
                "directional_agreement": summary["directional_agreement"],
                "discordant_genes": len(result["discordant_genes"]),
            })

    summary_df = pd.DataFrame(rows)
    print("\n===== Gradient across imbalance levels =====")
    print(summary_df.to_string(index=False))

    md_parts.append("## Gradient across imbalance levels\n\n```\n")
    md_parts.append(summary_df.to_string(index=False))
    md_parts.append("\n```\n")

    with open(md_path, "w") as f:
        f.write("".join(md_parts))

    print(f"\nWrote {md_path}")


if __name__ == "__main__":
    main()
