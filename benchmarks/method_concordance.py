"""
DEConcord's actual core question, run for real: how much do independently
developed DE tools agree with *each other* on the same count matrix?
Covers all three pairwise comparisons among DESeq2, edgeR, and limma-voom.

`compare_results.py` asks a different, older question -- "how close is
DEConcord's own simple t-test to DESeq2/edgeR" -- left in place because
it's still a real, useful check (see benchmarks/README.md). This script
is the new one: it runs no DE method itself, just loads whatever
DESeq2/edgeR/limma-voom results `run_deseq2_edger.R` (or
`run_deseq2_edger_pasilla.R`) and `run_limma_voom.R` already produced and
hands each pair to `deconcord.concordance.methods.compute_de_concordance`
-- the same library function anyone using DEConcord would call on their
own two result tables.

Usage (from the repo root, after the matching R scripts have run):
    python benchmarks/method_concordance.py
    python benchmarks/method_concordance.py --dataset pasilla
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from deconcord.concordance.methods import compute_de_concordance

RESULTS_DIR = Path(__file__).resolve().parent / "results"

DATASETS = {
    "airway": {"result_prefix": ""},
    "pasilla": {"result_prefix": "pasilla_"},
}

# Each tool's result file suffix and its own raw column names -- same
# names benchmarks/run_deseq2_edger*.R and benchmarks/run_limma_voom.R
# write, unchanged, so nothing here renames a tool's native output.
METHODS = {
    "DESeq2": {"file_suffix": "deseq2_results.csv", "lfc_col": "log2FoldChange", "pvalue_col": "padj"},
    "edgeR": {"file_suffix": "edger_results.csv", "lfc_col": "logFC", "pvalue_col": "FDR"},
    "limma-voom": {"file_suffix": "limma_voom_results.csv", "lfc_col": "logFC", "pvalue_col": "adj.P.Val"},
}

# DESeq2 vs edgeR first and in this position, matching the original,
# single-pair version of this script.
PAIRS = [("DESeq2", "edgeR"), ("DESeq2", "limma-voom"), ("edgeR", "limma-voom")]


def _require(path: Path, hint: str) -> Path:
    if not path.exists():
        raise SystemExit(f"Missing {path}. {hint}")
    return path


def _run_pair(name_a: str, name_b: str, prefix: str, hint: str) -> dict:
    path_a = _require(RESULTS_DIR / f"{prefix}{METHODS[name_a]['file_suffix']}", hint)
    path_b = _require(RESULTS_DIR / f"{prefix}{METHODS[name_b]['file_suffix']}", hint)
    df_a = pd.read_csv(path_a, index_col="gene_id")
    df_b = pd.read_csv(path_b, index_col="gene_id")

    result = compute_de_concordance(
        df_a, df_b, name_a=name_a, name_b=name_b,
        lfc_col_a=METHODS[name_a]["lfc_col"], pvalue_col_a=METHODS[name_a]["pvalue_col"],
        lfc_col_b=METHODS[name_b]["lfc_col"], pvalue_col_b=METHODS[name_b]["pvalue_col"],
    )

    # compute_de_concordance's own "directional_agreement" is restricted to
    # genes significant in *both* tables. Add the significant-set union and
    # the same directional check widened to genes significant in *either*
    # table, since a gene only one tool calls significant is exactly the
    # marginal, threshold-adjacent case this benchmark is meant to surface.
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=sorted(DATASETS), default="airway")
    args = parser.parse_args()

    prefix = DATASETS[args.dataset]["result_prefix"]
    r_script = "run_deseq2_edger.R" if args.dataset == "airway" else "run_deseq2_edger_pasilla.R"
    hint = f"Run `Rscript benchmarks/{r_script}` and `Rscript benchmarks/run_limma_voom.R` first."

    suffix = "" if args.dataset == "airway" else f"_{args.dataset}"
    md_path = RESULTS_DIR / f"method_concordance{suffix}.md"
    md_parts = []

    for name_a, name_b in PAIRS:
        result = _run_pair(name_a, name_b, prefix, hint)
        summary = result["summary"]

        print(f"{name_a} vs {name_b} -- {args.dataset} dataset")
        for key, value in summary.items():
            print(f"  {key}: {value:.4f}" if isinstance(value, float) else f"  {key}: {value}")
        print(f"  concordant genes: {len(result['concordant_genes'])}")
        print(f"  discordant genes: {len(result['discordant_genes'])}")
        print(f"  only in {name_a}: {len(result[f'only_in_{name_a}'])}")
        print(f"  only in {name_b}: {len(result[f'only_in_{name_b}'])}")
        print()

        section = [f"# {name_a} vs {name_b} concordance — {args.dataset} dataset\n\n"]
        section.append(
            "Real numbers from a real run, computed by "
            "`deconcord.concordance.methods.compute_de_concordance` on the "
            "committed DE result tables. Both tools ran with their own "
            "default settings on the identical count matrix and design "
            "-- neither is being validated against the other here, they're "
            "being checked for whether they agree.\n\n"
        )
        section.append("```\n")
        for key, value in summary.items():
            section.append(f"{key}: {value:.4f}\n" if isinstance(value, float) else f"{key}: {value}\n")
        section.append(f"concordant_genes: {len(result['concordant_genes'])}\n")
        section.append(f"discordant_genes: {len(result['discordant_genes'])}\n")
        section.append(f"only_in_{name_a}: {len(result[f'only_in_{name_a}'])}\n")
        section.append(f"only_in_{name_b}: {len(result[f'only_in_{name_b}'])}\n")
        section.append("```\n")
        if len(result["discordant_genes"]) > 0:
            section.append(f"\nDiscordant genes (significant in both, opposite-signed log fold change): "
                            f"{list(result['discordant_genes'])}\n")
        md_parts.append("".join(section))

    with open(md_path, "w") as f:
        f.write("\n".join(md_parts))

    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
