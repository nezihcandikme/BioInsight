"""
DEConcord's actual core question, run for real: how much do DESeq2 and
edgeR agree with *each other* on the same count matrix?

`compare_results.py` asks a different, older question -- "how close is
DEConcord's own simple t-test to DESeq2/edgeR" -- left in place because
it's still a real, useful check (see benchmarks/README.md). This script
is the new one: it runs no DE method itself, just loads whatever
DESeq2/edgeR results `run_deseq2_edger.R` (or `run_deseq2_edger_pasilla.R`)
already produced and hands them to
`deconcord.concordance.methods.compute_de_concordance` -- the same
library function anyone using DEConcord would call on their own two
result tables.

Usage (from the repo root, after the matching R script has run):
    python benchmarks/method_concordance.py
    python benchmarks/method_concordance.py --dataset pasilla
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from deconcord.concordance.methods import compute_de_concordance

RESULTS_DIR = Path(__file__).resolve().parent / "results"

DATASETS = {
    "airway": {"result_prefix": ""},
    "pasilla": {"result_prefix": "pasilla_"},
}


def _require(path: Path, hint: str) -> Path:
    if not path.exists():
        raise SystemExit(f"Missing {path}. {hint}")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=sorted(DATASETS), default="airway")
    args = parser.parse_args()

    cfg = DATASETS[args.dataset]
    r_script = "run_deseq2_edger.R" if args.dataset == "airway" else "run_deseq2_edger_pasilla.R"
    hint = f"Run `Rscript benchmarks/{r_script}` first."

    deseq2_path = _require(RESULTS_DIR / f"{cfg['result_prefix']}deseq2_results.csv", hint)
    edger_path = _require(RESULTS_DIR / f"{cfg['result_prefix']}edger_results.csv", hint)

    deseq2 = pd.read_csv(deseq2_path, index_col="gene_id")
    edger = pd.read_csv(edger_path, index_col="gene_id")

    result = compute_de_concordance(
        deseq2, edger, name_a="DESeq2", name_b="edgeR",
        lfc_col_a="log2FoldChange", pvalue_col_a="padj",
        lfc_col_b="logFC", pvalue_col_b="FDR",
    )
    summary = result["summary"]

    print(f"DESeq2 vs edgeR -- {args.dataset} dataset")
    for key, value in summary.items():
        print(f"  {key}: {value:.4f}" if isinstance(value, float) else f"  {key}: {value}")
    print(f"  concordant genes: {len(result['concordant_genes'])}")
    print(f"  discordant genes: {len(result['discordant_genes'])}")
    print(f"  only in DESeq2: {len(result['only_in_DESeq2'])}")
    print(f"  only in edgeR: {len(result['only_in_edgeR'])}")

    suffix = "" if args.dataset == "airway" else f"_{args.dataset}"
    md_path = RESULTS_DIR / f"method_concordance{suffix}.md"
    with open(md_path, "w") as f:
        f.write(f"# DESeq2 vs edgeR concordance — {args.dataset} dataset\n\n")
        f.write(
            "Real numbers from a real run, computed by "
            "`deconcord.concordance.methods.compute_de_concordance` on the "
            "committed DESeq2/edgeR result tables. Both tools ran with their "
            "own default settings on the identical count matrix and design "
            "-- neither is being validated against the other here, they're "
            "being checked for whether they agree.\n\n"
        )
        f.write("```\n")
        for key, value in summary.items():
            f.write(f"{key}: {value:.4f}\n" if isinstance(value, float) else f"{key}: {value}\n")
        f.write(f"concordant_genes: {len(result['concordant_genes'])}\n")
        f.write(f"discordant_genes: {len(result['discordant_genes'])}\n")
        f.write(f"only_in_DESeq2: {len(result['only_in_DESeq2'])}\n")
        f.write(f"only_in_edgeR: {len(result['only_in_edgeR'])}\n")
        f.write("```\n")
        if len(result["discordant_genes"]) > 0:
            f.write(f"\nDiscordant genes (significant in both, opposite-signed log fold change): "
                     f"{list(result['discordant_genes'])}\n")

    print(f"\nWrote {md_path}")


if __name__ == "__main__":
    main()
