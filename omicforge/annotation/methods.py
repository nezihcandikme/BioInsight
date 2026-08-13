"""
Gene ID -> gene symbol annotation.

Differential expression results are indexed by whatever gene ID the input
count matrix used -- commonly an Ensembl ID like ``ENSG00000141510``,
correct but not exactly something you can eyeball. This module doesn't
fetch symbols from Ensembl, BioMart, or any other live service: it reads
a plain two-column mapping file you already have (or exported once,
outside OmicForge), the same "point it at a local file" approach
``omicforge.pathway_analysis.gmt`` takes for gene sets. No network
dependency, no API to go down, and no surprise about which annotation
build/release the symbols came from -- that's on whatever file you hand
it.
"""

import os

import pandas as pd


def load_gene_annotation(path: str, id_col: str = "gene_id", symbol_col: str = "gene_symbol") -> dict[str, str]:
    """
    Load a two-column gene ID -> gene symbol mapping from a CSV file.

    Expected format (header required)::

        gene_id,gene_symbol
        ENSG00000141510,TP53
        ENSG00000146648,EGFR

    Parameters
    ----------
    path : str
        Path to the mapping CSV.
    id_col, symbol_col : str, optional
        Column names to read the gene ID and gene symbol from. Defaults
        match the format above; override if your export uses different
        header names.

    Returns
    -------
    dict[str, str]
        Gene ID -> gene symbol, ready to pass to ``annotate_de_results``
        or ``omicforge.pipeline.run_analysis``.

    Raises
    ------
    FileNotFoundError
        If ``path`` doesn't exist.
    ValueError
        If ``id_col`` or ``symbol_col`` isn't a column in the file, if
        the file has a header but no data rows, if any row is missing an
        ID or a symbol, or if the same gene ID maps to two different
        symbols in different rows -- each of these means the resulting
        dict would be incomplete or ambiguous in a way that's better
        caught here than discovered later as a confusing NaN or a wrong
        symbol downstream.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"No such gene annotation file: {path}")

    df = pd.read_csv(path, dtype=str)

    missing_cols = [c for c in (id_col, symbol_col) if c not in df.columns]
    if missing_cols:
        raise ValueError(
            f"{path}: missing column(s) {missing_cols}. Found columns: {list(df.columns)}."
        )

    if df.empty:
        raise ValueError(f"{path} has a header but no data rows.")

    incomplete = df[df[id_col].isna() | df[symbol_col].isna()]
    if not incomplete.empty:
        bad_rows = [i + 2 for i in incomplete.index]  # +2: 1-indexed, plus the header line
        raise ValueError(
            f"{path}: row(s) {bad_rows} have a missing {id_col} or {symbol_col} value."
        )

    annotation: dict[str, str] = {}
    for gene_id, group in df.groupby(id_col)[symbol_col]:
        symbols = set(group)
        if len(symbols) > 1:
            raise ValueError(
                f"{path}: gene ID '{gene_id}' maps to more than one symbol: "
                f"{sorted(symbols)}. Each gene ID must map to exactly one symbol."
            )
        annotation[gene_id] = group.iloc[0]

    return annotation


def annotate_de_results(de_results: pd.DataFrame, annotation: dict[str, str]) -> pd.DataFrame:
    """
    Attach a ``gene_symbol`` column to a differential expression results
    table, without touching the original gene IDs.

    Parameters
    ----------
    de_results : pd.DataFrame
        Output of ``run_differential_expression`` (or anything else
        indexed by gene ID) -- not modified in place.
    annotation : dict[str, str]
        Gene ID -> gene symbol, e.g. from ``load_gene_annotation``.

    Returns
    -------
    pd.DataFrame
        A copy of ``de_results`` with ``gene_symbol`` inserted as the
        first column. Genes not present in ``annotation`` get ``NaN`` --
        deliberately not the original ID, so "unmapped" stays
        distinguishable from "symbol happens to equal the ID." The index
        (the original gene ID) is unchanged.
    """
    annotated = de_results.copy()
    gene_symbols = annotated.index.to_series().map(annotation)
    annotated.insert(0, "gene_symbol", gene_symbols)
    return annotated
