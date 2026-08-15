"""
DEConcord: robustness and concordance analysis for RNA-seq differential
expression.

This top-level module exposes the functions most users need without
requiring you to know which internal submodule they live in. Every name
below is also importable from its actual submodule (e.g.
``deconcord.concordance.methods.compute_de_concordance``) if you'd rather
be explicit -- nothing here is exclusive to the top level.

Typical usage::

    import deconcord as dc

    # Compare two existing DE result tables (DEConcord's core question):
    result = dc.compute_de_concordance(deseq2_results, edger_results, ...)

    # Or run the underlying RNA-seq pipeline that produces a DE table:
    results = dc.run_analysis("counts.csv", group_1=[...], group_2=[...])
"""

from importlib.metadata import PackageNotFoundError, version

from deconcord.concordance.methods import compute_de_concordance
from deconcord.concordance.pathway_stability import compute_pathway_stability
from deconcord.concordance.threshold_sensitivity import compute_threshold_sensitivity
from deconcord.differential_expression.methods import (
    run_differential_expression,
    run_differential_expression_with_covariates,
)
from deconcord.io.counts import load_count_matrix, validate_counts
from deconcord.normalization.methods import compute_cpm, log2_transform
from deconcord.pipeline import run_analysis
from deconcord.qc.metrics import run_sample_qc
from deconcord.visualization.plots import plot_pathway_enrichment, plot_pca, plot_volcano

try:
    __version__ = version("deconcord")
except PackageNotFoundError:  # pragma: no cover - only hit if run from an uninstalled checkout
    __version__ = "0.0.0+unknown"

__all__ = [
    "__version__",
    "compute_cpm",
    "compute_de_concordance",
    "compute_pathway_stability",
    "compute_threshold_sensitivity",
    "load_count_matrix",
    "log2_transform",
    "plot_pathway_enrichment",
    "plot_pca",
    "plot_volcano",
    "run_analysis",
    "run_differential_expression",
    "run_differential_expression_with_covariates",
    "run_sample_qc",
    "validate_counts",
]
