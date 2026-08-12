"""
Command-line entry point for BioInsight.

This wraps ``bioinsight.pipeline.run_analysis`` so the pipeline can be run
without writing a Python script: point it at a count matrix CSV, name your
two groups, and it writes the differential expression table (and, if asked,
plots and a pathway enrichment table) to an output directory.

It is deliberately a thin wrapper. All the actual logic -- validation, QC,
filtering, normalization, DE, enrichment -- lives in the library modules;
this module's only job is argument parsing, wiring, and turning exceptions
that are expected to happen on bad input (a missing file, an unknown sample
name, a malformed GMT file) into a clean error message instead of a
traceback.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from bioinsight.io.counts import CountMatrixError
from bioinsight.pathway_analysis.gmt import load_gmt
from bioinsight.pipeline import run_analysis


def _parse_sample_list(raw: str) -> list[str]:
    samples = [s.strip() for s in raw.split(",") if s.strip()]
    if not samples:
        raise ValueError("Expected a comma-separated list of sample names, got an empty one.")
    return samples


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bioinsight",
        description="Run the BioInsight RNA-seq pipeline on a count matrix CSV.",
    )
    parser.add_argument("counts", help="Path to a count matrix CSV (genes as rows, samples as columns).")
    parser.add_argument("--group1", required=True, help="Comma-separated sample names for the first comparison group.")
    parser.add_argument("--group2", required=True, help="Comma-separated sample names for the second comparison group.")
    parser.add_argument("--gmt", help="Path to a .gmt gene-set file to run pathway enrichment against.")
    parser.add_argument(
        "--background",
        help="Comma-separated background gene list for enrichment. Defaults to every gene in the "
             "(filtered, if --min-count is set) count matrix if --gmt is given without this.",
    )
    parser.add_argument("--min-count", type=int, default=None, help="Minimum raw read count for a gene to pass pre-DE filtering. Off by default.")
    parser.add_argument("--min-samples", type=int, default=2, help="Number of samples that must clear --min-count. Default 2. Ignored without --min-count.")
    parser.add_argument("--alpha", type=float, default=0.05, help="Adjusted p-value cutoff for significance. Default 0.05.")
    parser.add_argument("--lfc-threshold", type=float, default=1.0, help="Absolute log fold change cutoff for significance. Default 1.0.")
    parser.add_argument("--no-plots", action="store_true", help="Skip generating the volcano and PCA plots.")
    parser.add_argument("--explain", action="store_true", help="Generate a plain-language explanation of the results (requires ANTHROPIC_API_KEY).")
    parser.add_argument("--out", default="bioinsight_output", help="Output directory. Default 'bioinsight_output'.")
    return parser


def _run(argv: list[str]) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    group_1 = _parse_sample_list(args.group1)
    group_2 = _parse_sample_list(args.group2)

    gene_sets = None
    background_genes = None
    if args.gmt is not None:
        gene_sets = load_gmt(args.gmt)
        if args.background is not None:
            background_genes = set(_parse_sample_list(args.background))
        else:
            # No explicit background given: fall back to every gene in the
            # count matrix, which is the honest "everything that could have
            # been tested" universe for this run.
            background_genes = set(pd.read_csv(args.counts, index_col=0).index)

    results = run_analysis(
        args.counts,
        group_1=group_1,
        group_2=group_2,
        gene_sets=gene_sets,
        background_genes=background_genes,
        min_count=args.min_count,
        min_samples=args.min_samples,
        alpha=args.alpha,
        lfc_threshold=args.lfc_threshold,
        generate_plots=not args.no_plots,
        explain_results=args.explain,
    )

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    de = results["differential_expression"]
    de.to_csv(out_dir / "differential_expression.csv")
    results["qc"].to_csv(out_dir / "qc.csv")

    if "pathway_enrichment" in results:
        results["pathway_enrichment"].to_csv(out_dir / "pathway_enrichment.csv", index=False)

    if "volcano_fig" in results:
        results["volcano_fig"].savefig(out_dir / "volcano.png", dpi=150, bbox_inches="tight")
    if "pca_fig" in results:
        results["pca_fig"].savefig(out_dir / "pca.png", dpi=150, bbox_inches="tight")

    if "explanation" in results:
        (out_dir / "explanation.txt").write_text(results["explanation"])

    n_sig = int(de["significant"].sum())
    print(f"{len(de)} genes tested, {n_sig} significant at alpha={args.alpha}, |log2FC|>{args.lfc_threshold}.")
    print(f"Results written to {out_dir}/")

    return 0


def main(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    try:
        return _run(argv)
    except (CountMatrixError, ValueError, FileNotFoundError) as exc:
        print(f"bioinsight: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
