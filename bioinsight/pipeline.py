"""
End-to-end orchestration for BioInsight.

``run_analysis`` wires together the individual modules (I/O validation, QC,
normalization, differential expression, visualization, and optional pathway
enrichment / AI explanation) into a single call, so users don't have to
manually chain every step together.

NOTE: the differential expression step is a basic, exploratory
mean-difference t-test (see ``bioinsight.differential_expression.methods``),
not a replacement for count-based tools like DESeq2 or edgeR.
"""

from __future__ import annotations

import pandas as pd

from bioinsight.io.counts import load_count_matrix, validate_counts
from bioinsight.qc.metrics import run_sample_qc
from bioinsight.normalization.methods import compute_cpm, log2_transform
from bioinsight.differential_expression.methods import run_differential_expression
from bioinsight.visualization.plots import plot_volcano, plot_pca
from bioinsight.pathway_analysis.methods import run_pathway_enrichment_analysis


def run_analysis(
    counts: str | pd.DataFrame,
    group_1: list[str],
    group_2: list[str],
    gene_sets: dict[str, set[str]] | None = None,
    background_genes: set[str] | None = None,
    generate_plots: bool = True,
    explain_results: bool = False,
) -> dict:
    """
    Run the full BioInsight workflow: load/validate, QC, normalize,
    differential expression, and (optionally) plots, pathway enrichment,
    and an AI-generated explanation.

    Parameters
    ----------
    counts : str or pd.DataFrame
        Either a path to a CSV count matrix (genes as rows, samples as
        columns) or an already-loaded raw count DataFrame. If a DataFrame
        is passed, it is still validated via ``validate_counts``.
    group_1 : list[str]
        Sample names (column names in ``counts``) belonging to the first
        comparison group.
    group_2 : list[str]
        Sample names (column names in ``counts``) belonging to the second
        comparison group.
    gene_sets : dict[str, set[str]], optional
        Pathway/gene-set definitions for enrichment analysis. If provided
        together with ``background_genes``, pathway enrichment is run on
        the significant genes from the differential expression step.
    background_genes : set[str], optional
        The background gene universe for enrichment analysis. Required if
        ``gene_sets`` is provided.
    generate_plots : bool, optional
        Whether to generate a volcano plot and a PCA plot. Default True.
    explain_results : bool, optional
        Whether to generate a natural-language explanation of the
        differential expression results using
        ``bioinsight.ai_explanation.methods.explain_de_results``. Requires
        an ``ANTHROPIC_API_KEY`` to be set. Default False.

    Returns
    -------
    dict
        A dictionary that may contain the following keys:
        - "raw_counts": the validated raw count matrix
        - "qc": per-sample QC metrics (library size, genes detected, outlier flags)
        - "normalized": CPM-normalized, log2-transformed expression matrix
        - "differential_expression": DE results (log fold change, p-values, significance)
        - "volcano_fig", "pca_fig": matplotlib Figures (if generate_plots=True)
        - "pathway_enrichment": enrichment results (if gene_sets/background_genes provided)
        - "explanation": AI-generated summary text (if explain_results=True)
    """
    if isinstance(counts, str):
        raw_counts = load_count_matrix(counts)
    else:
        validate_counts(counts)
        raw_counts = counts

    missing_samples = [s for s in group_1 + group_2 if s not in raw_counts.columns]
    if missing_samples:
        raise ValueError(
            f"The following samples in group_1/group_2 are not columns of "
            f"the count matrix: {missing_samples}"
        )

    results: dict = {"raw_counts": raw_counts}

    results["qc"] = run_sample_qc(raw_counts)

    cpm = compute_cpm(raw_counts)
    normalized = log2_transform(cpm)
    results["normalized"] = normalized

    de_results = run_differential_expression(normalized, group_1, group_2)
    results["differential_expression"] = de_results

    if generate_plots:
        involved_samples = group_1 + group_2
        labels = pd.Series(
            {**{s: "group_1" for s in group_1}, **{s: "group_2" for s in group_2}}
        )[involved_samples]

        results["volcano_fig"] = plot_volcano(de_results)
        results["pca_fig"] = plot_pca(normalized[involved_samples], labels)

    if gene_sets is not None:
        if background_genes is None:
            raise ValueError(
                "background_genes must be provided when gene_sets is given."
            )
        significant_genes = set(de_results.index[de_results["significant"]])
        results["pathway_enrichment"] = run_pathway_enrichment_analysis(
            significant_genes, gene_sets, background_genes
        )

    if explain_results:
        # Imported lazily: this pulls in the anthropic SDK / python-dotenv,
        # which shouldn't be required unless AI explanations are requested.
        from bioinsight.ai_explanation.methods import explain_de_results

        results["explanation"] = explain_de_results(de_results)

    return results
