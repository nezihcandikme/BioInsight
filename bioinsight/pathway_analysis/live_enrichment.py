"""
Live pathway enrichment via the g:Profiler g:GOSt API.

This is deliberately separate from ``gmt.py``/``methods.py``: those work
entirely from a local gene-set file and a local hypergeometric test, no
network required — the same "no live external service" design
``bioinsight.annotation`` uses, for the same reason (this project is built
and tested in a sandbox with no general internet access). This module is
the one place in pathway analysis that breaks that pattern on purpose,
because g:Profiler's curated, versioned pathway databases (GO, KEGG,
Reactome, WikiPathways, and more) are broader and more current than any
single GMT file someone happens to have on disk, and reimplementing that
database locally isn't a reasonable thing for a project this size to take
on.

It is opt-in and self-contained: nothing else in the package imports this
module automatically, and ``run_analysis`` only calls it if
``live_enrichment_organism`` is explicitly passed. It requires outbound
internet access to ``biit.cs.ut.ee``, which this project's own development
sandbox does not have (confirmed: every one of g:Profiler, Enrichr, and
NCBI's endpoints time out from here) — run it from an environment that can
actually reach the public internet, and expect it to fail loudly, not
silently, if it can't.
"""

import json
import urllib.error
import urllib.request

import pandas as pd

_GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

_RESULT_COLUMNS = [
    "source", "native", "name", "p_value", "term_size", "query_size", "intersection_size",
]


def run_gprofiler_enrichment(
    gene_list: list[str] | set[str],
    organism: str = "hsapiens",
    sources: list[str] | None = None,
    significance_threshold_method: str = "g_SCS",
    user_threshold: float = 0.05,
    timeout: float = 30.0,
) -> pd.DataFrame:
    """
    Query g:Profiler's g:GOSt enrichment endpoint for a gene list.

    Parameters
    ----------
    gene_list : list[str] or set[str]
        Genes to test for enrichment (e.g. the significant genes from a
        differential expression result). Gene ID format must match what
        ``organism`` expects — g:Profiler accepts Ensembl IDs, gene
        symbols, and several other identifier types, and resolves them
        itself; it does not require IDs to be pre-converted.
    organism : str, optional
        g:Profiler organism code. Default ``"hsapiens"`` (human). See
        https://biit.cs.ut.ee/gprofiler/page/organism-list for the full list.
    sources : list[str], optional
        Restrict results to specific databases (e.g. ``["GO:BP", "KEGG"]``).
        Default ``None`` — g:Profiler's own default set of sources.
    significance_threshold_method : str, optional
        g:Profiler's multiple-testing correction method. Default
        ``"g_SCS"`` (g:Profiler's own default, a set-counting correction
        tuned for this exact kind of enrichment test). Other options
        include ``"fdr"`` and ``"bonferroni"``.
    user_threshold : float, optional
        Significance cutoff applied by g:Profiler itself. Default 0.05.
    timeout : float, optional
        Request timeout in seconds. Default 30.

    Returns
    -------
    pd.DataFrame
        One row per significant term, sorted by p-value, with columns
        ``source``, ``native`` (the term's database-native ID), ``name``,
        ``p_value``, ``term_size``, ``query_size``, ``intersection_size``.
        Empty (but correctly-columned) if g:Profiler found nothing
        significant — that's a real result, not treated as an error.

    Raises
    ------
    ValueError
        If ``gene_list`` is empty.
    ConnectionError
        If g:Profiler can't be reached at all (DNS failure, timeout,
        connection refused) — this is the expected failure mode in a
        network-restricted environment, so it's raised as a specific,
        identifiable error rather than letting a generic ``URLError``
        propagate.
    RuntimeError
        If g:Profiler is reachable but returns a non-2xx HTTP response
        (e.g. a malformed request or an unrecognized organism code).
    """
    genes = sorted(set(gene_list))
    if not genes:
        raise ValueError("gene_list is empty — nothing to query g:Profiler with.")

    payload = {
        "organism": organism,
        "query": genes,
        "significance_threshold_method": significance_threshold_method,
        "user_threshold": user_threshold,
    }
    if sources:
        payload["sources"] = sources

    request = urllib.request.Request(
        _GPROFILER_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )

    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            body = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        raise RuntimeError(
            f"g:Profiler returned HTTP {e.code} for this request "
            f"(organism='{organism}', {len(genes)} genes): {e.reason}. "
            "This usually means an unrecognized organism code or a "
            "malformed gene list, not a network problem."
        ) from e
    except urllib.error.URLError as e:
        raise ConnectionError(
            f"Could not reach g:Profiler ({_GPROFILER_URL}): {e.reason}. "
            "This function needs outbound internet access to biit.cs.ut.ee; "
            "it will not work in a network-restricted environment."
        ) from e

    terms = body.get("result", [])
    if not terms:
        return pd.DataFrame(columns=_RESULT_COLUMNS)

    df = pd.DataFrame(terms)
    present_columns = [c for c in _RESULT_COLUMNS if c in df.columns]
    return df[present_columns].sort_values("p_value").reset_index(drop=True)
