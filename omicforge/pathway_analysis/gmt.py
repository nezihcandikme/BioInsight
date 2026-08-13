"""
GMT gene-set file loading.

GMT (Gene Matrix Transposed) is the plain tab-separated format MSigDB and
most pathway databases ship gene sets in: one gene set per line, the set's
name in the first column, a description or source URL in the second, and
every gene in the set after that. This is a thin parser, not a client for
any particular database — point it at any .gmt file (MSigDB, an Enrichr
export, an in-house list) and it becomes the same ``dict[str, set[str]]``
``run_pathway_enrichment_analysis`` already expects for ``gene_sets``.
"""

import os


def load_gmt(path: str) -> dict[str, set[str]]:
    """
    Parse a .gmt file into a gene-set dictionary.

    Parameters
    ----------
    path : str
        Path to a GMT file. Each line: ``name<TAB>description<TAB>gene1<TAB>gene2...``.
        Blank lines are skipped.

    Returns
    -------
    dict[str, set[str]]
        Gene set name -> set of gene IDs, ready to pass as ``gene_sets`` to
        ``run_pathway_enrichment_analysis`` or ``omicforge.pipeline.run_analysis``.

    Raises
    ------
    FileNotFoundError
        If ``path`` doesn't exist.
    ValueError
        If a non-blank line has fewer than 3 tab-separated fields (a gene
        set needs a name, a description, and at least one gene), if a gene
        set name appears more than once in the file, or if the file has no
        gene sets at all — each of these means whatever downstream
        enrichment result you get would be silently built on incomplete or
        ambiguous input.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"No such GMT file: {path}")

    gene_sets: dict[str, set[str]] = {}

    with open(path) as f:
        for line_num, raw_line in enumerate(f, start=1):
            line = raw_line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue

            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"{path}, line {line_num}: found {len(fields)} tab-separated "
                    "field(s), but a GMT line needs at least a name, a "
                    "description, and one gene."
                )

            name, _description, *gene_fields = fields
            genes = {g for g in gene_fields if g}

            if not genes:
                raise ValueError(f"{path}, line {line_num}: gene set '{name}' has no genes.")

            if name in gene_sets:
                raise ValueError(
                    f"{path}, line {line_num}: gene set name '{name}' already "
                    "appeared earlier in this file — gene set names must be unique."
                )

            gene_sets[name] = genes

    if not gene_sets:
        raise ValueError(f"{path} contains no gene sets.")

    return gene_sets
