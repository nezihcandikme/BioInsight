"""
End-to-end orchestration for DEConcord.

``run_analysis`` wires together the individual modules (I/O validation, QC,
normalization, differential expression, visualization, and optional pathway
enrichment / AI explanation) into a single call, so users don't have to
manually chain every step together.

NOTE: the differential expression step is a basic, exploratory
mean-difference t-test (see ``deconcord.differential_expression.methods``),
not a replacement for count-based tools like DESeq2 or edgeR.
"""

from __future__ import annotations

import pandas as pd

from deconcord.io.counts import load_count_matrix, validate_counts
from deconcord.qc.metrics import run_sample_qc
from deconcord.filtering.methods import filter_low_count_genes
from deconcord.normalization.methods import compute_cpm, log2_transform
from deconcord.differential_expression.methods import (
    run_differential_expression,
    run_differential_expression_with_covariates,
)
from deconcord.visualization.plots import plot_volcano, plot_pca, plot_pathway_enrichment
from deconcord.pathway_analysis.methods import run_pathway_enrichment_analysis
from deconcord.annotation.methods import annotate_de_results, load_gene_annotation


def run_analysis(
    counts: str | pd.DataFrame,
    group_1: list[str],
    group_2: list[str],
    gene_sets: dict[str, set[str]] | None = None,
    background_genes: set[str] | None = None,
    min_count: int | None = None,
    min_samples: int = 2,
    alpha: float = 0.05,
    lfc_threshold: float = 1.0,
    method: str = "welch",
    metadata: pd.DataFrame | None = None,
    covariate_cols: list[str] | None = None,
    moderated_covariates: bool = True,
    gene_annotation: dict[str, str] | str | None = None,
    generate_plots: bool = True,
    explain_results: bool = False,
    live_enrichment_organism: str | None = None,
) -> dict:
    """
    Run the full DEConcord workflow: load/validate, QC, normalize,
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
    min_count : int, optional
        If given, genes are filtered via ``filter_low_count_genes`` (kept
        only if at least ``min_samples`` samples have >= ``min_count`` raw
        reads) before normalization and DE. Default ``None`` — no
        filtering, to keep existing callers' behavior unchanged. Turning
        this on is recommended for real data: untested, barely-expressed
        genes cost statistical power for every other gene through the
        multiple-testing correction, for no benefit.
    min_samples : int, optional
        Passed through to ``filter_low_count_genes`` when ``min_count`` is
        set. Default 2. Ignored if ``min_count`` is ``None``.
    alpha : float, optional
        Adjusted p-value cutoff passed through to ``run_differential_expression``.
        Default 0.05.
    lfc_threshold : float, optional
        Log fold change cutoff passed through to ``run_differential_expression``.
        Default 1.0.
    method : str, optional
        ``"welch"`` (default, per-gene, no cross-gene assumptions) or
        ``"moderated"`` (empirical-Bayes variance shrinkage across genes —
        more statistical power, assumes each gene's two groups share one
        variance). Passed through to ``run_differential_expression``; see
        its docstring and ``benchmarks/`` for the measured tradeoff. Ignored
        if ``covariate_cols`` is given — see below.
    metadata : pd.DataFrame, optional
        Sample metadata indexed by sample name, covering every sample in
        ``group_1 + group_2``. Only used when ``covariate_cols`` is also
        given (e.g. as returned by
        ``deconcord.io.geo.fetch_geo_sample_metadata``, or your own CSV
        loaded with the sample name as the index).
    covariate_cols : list[str], optional
        Column names in ``metadata`` to adjust for — a batch, sex, or any
        other measured source of variation that isn't the condition of
        interest. If given, differential expression switches from
        ``run_differential_expression`` (a plain two-group test, ``method``
        selects Welch's or moderated) to
        ``run_differential_expression_with_covariates`` (a per-gene linear
        model with the condition and these covariates as terms — see its
        docstring for what "log fold change" means in that case, and why
        there's no Welch's-test option once covariates are involved).
        Default ``None`` — no covariate adjustment, existing callers see no
        change.
    moderated_covariates : bool, optional
        Only used when ``covariate_cols`` is given. Whether to shrink each
        gene's residual variance toward a cross-gene prior, same idea as
        ``method="moderated"`` for the plain two-group case. Default True.
    gene_annotation : dict[str, str] or str, optional
        Either a pre-loaded gene ID -> gene symbol dict, or a path to a
        two-column CSV to load one from (see
        ``deconcord.annotation.methods.load_gene_annotation``). If given,
        the differential expression results gain a ``gene_symbol`` column
        (``NaN`` for any gene not in the mapping). Default ``None`` — no
        annotation, no new column, existing callers see no change.
    generate_plots : bool, optional
        Whether to generate a volcano plot, a PCA plot, and (if any
        enrichment was run) a pathway enrichment dot plot per enrichment
        source. Default True.
    explain_results : bool, optional
        Whether to generate a natural-language explanation of the
        differential expression results using
        ``deconcord.ai_explanation.methods.explain_de_results``. Requires
        an ``ANTHROPIC_API_KEY`` to be set. Default False.
    live_enrichment_organism : str, optional
        If given, also runs pathway enrichment against the live g:Profiler
        API (``deconcord.pathway_analysis.live_enrichment.run_gprofiler_enrichment``)
        on the significant genes, using this as the g:Profiler organism
        code (e.g. ``"hsapiens"``). This is independent of ``gene_sets``/
        ``background_genes`` — both can be used together, or either alone.
        Requires outbound internet access; raises ``ConnectionError`` if
        g:Profiler can't be reached. Default ``None`` — no live call, no
        network dependency, existing callers see no change.

    Returns
    -------
    dict
        A dictionary that may contain the following keys:
        - "raw_counts": the validated raw count matrix
        - "qc": per-sample QC metrics (library size, genes detected, outlier flags)
        - "filtered_counts": raw counts after low-count filtering (only if
          ``min_count`` is set)
        - "normalized": CPM-normalized, log2-transformed expression matrix
          (computed from filtered counts if filtering was applied)
        - "differential_expression": DE results (log fold change, p-values,
          significance, plus a "gene_symbol" column if gene_annotation was given)
        - "volcano_fig", "pca_fig": matplotlib Figures (if generate_plots=True)
        - "pathway_enrichment": local GMT-based enrichment results (if
          gene_sets/background_genes provided)
        - "pathway_enrichment_fig": dot plot of the above (if
          generate_plots=True and the result isn't empty)
        - "live_pathway_enrichment": g:Profiler enrichment results (if
          live_enrichment_organism was given)
        - "live_pathway_enrichment_fig": dot plot of the above (if
          generate_plots=True and the result isn't empty)
        - "explanation": AI-generated summary text (if explain_results=True)
    """
    if isinstance(counts, str):
        raw_counts = load_count_matrix(counts)
    else:
        validate_counts(counts)
        raw_counts = counts

    # group_1/group_2 are re-checked (missing/duplicate/overlapping sample
    # names) inside run_differential_expression via the shared
    # _validate_groups helper, so there's no separate check needed here.

    results: dict = {"raw_counts": raw_counts}

    # QC runs on the raw, unfiltered matrix — dropping low-count genes
    # first would quietly change library size and genes-detected numbers
    # out from under the sample-quality check.
    results["qc"] = run_sample_qc(raw_counts)

    counts_for_analysis = raw_counts
    if min_count is not None:
        counts_for_analysis = filter_low_count_genes(raw_counts, min_count=min_count, min_samples=min_samples)
        results["filtered_counts"] = counts_for_analysis

    cpm = compute_cpm(counts_for_analysis)
    normalized = log2_transform(cpm)
    results["normalized"] = normalized

    if covariate_cols is not None:
        if metadata is None:
            raise ValueError("metadata must be provided when covariate_cols is given.")
        de_results = run_differential_expression_with_covariates(
            normalized, group_1, group_2, metadata, covariate_cols,
            alpha=alpha, lfc_threshold=lfc_threshold, moderated=moderated_covariates,
        )
    else:
        de_results = run_differential_expression(normalized, group_1, group_2, alpha=alpha, lfc_threshold=lfc_threshold, method=method)
    results["differential_expression"] = de_results

    if gene_annotation is not None:
        annotation = load_gene_annotation(gene_annotation) if isinstance(gene_annotation, str) else gene_annotation
        results["differential_expression"] = annotate_de_results(de_results, annotation)

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
        if generate_plots and not results["pathway_enrichment"].empty:
            results["pathway_enrichment_fig"] = plot_pathway_enrichment(results["pathway_enrichment"])

    if live_enrichment_organism is not None:
        # Imported lazily, same reasoning as explain_de_results below: this
        # is the one function in the package that makes a live network
        # call, so nothing pays for that possibility unless it's actually
        # requested.
        from deconcord.pathway_analysis.live_enrichment import run_gprofiler_enrichment

        significant_genes = set(de_results.index[de_results["significant"]])
        results["live_pathway_enrichment"] = run_gprofiler_enrichment(
            significant_genes, organism=live_enrichment_organism
        )
        if generate_plots and not results["live_pathway_enrichment"].empty:
            results["live_pathway_enrichment_fig"] = plot_pathway_enrichment(
                results["live_pathway_enrichment"],
                name_col="name", pvalue_col="p_value", size_col="intersection_size",
            )

    if explain_results:
        # Imported lazily: this pulls in the anthropic SDK / python-dotenv,
        # which shouldn't be required unless AI explanations are requested --
        # they're an optional extra (`pip install deconcord[ai]`), not a
        # core dependency.
        try:
            from deconcord.ai_explanation.methods import explain_de_results
        except ImportError as exc:
            raise ImportError(
                "explain_results=True requires the optional 'ai' extra "
                "(anthropic, python-dotenv). Install it with: "
                "pip install deconcord[ai]"
            ) from exc

        results["explanation"] = explain_de_results(de_results)

    return results
