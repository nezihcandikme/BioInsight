import json
import urllib.error
from unittest.mock import patch

import pytest

from bioinsight.pathway_analysis.live_enrichment import run_gprofiler_enrichment


class _FakeResponse:
    """Minimal stand-in for the object urllib.request.urlopen returns,
    usable as a context manager the way the real one is."""

    def __init__(self, body: dict):
        self._body = json.dumps(body).encode("utf-8")

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


_SAMPLE_RESULT = {
    "result": [
        {
            "source": "GO:BP", "native": "GO:0006915", "name": "apoptotic process",
            "p_value": 0.02, "term_size": 500, "query_size": 20, "intersection_size": 5,
            "extra_field_gprofiler_might_add": "ignored",
        },
        {
            "source": "KEGG", "native": "KEGG:04210", "name": "Apoptosis",
            "p_value": 0.001, "term_size": 100, "query_size": 20, "intersection_size": 8,
        },
    ]
}


def test_empty_gene_list_raises():
    with pytest.raises(ValueError, match="empty"):
        run_gprofiler_enrichment([])


@patch("urllib.request.urlopen")
def test_successful_query_returns_sorted_dataframe(mock_urlopen):
    mock_urlopen.return_value = _FakeResponse(_SAMPLE_RESULT)

    result = run_gprofiler_enrichment(["TP53", "EGFR"], organism="hsapiens")

    assert list(result.columns) == [
        "source", "native", "name", "p_value", "term_size", "query_size", "intersection_size",
    ]
    assert len(result) == 2
    # Sorted by p_value ascending -- KEGG (0.001) before GO:BP (0.02).
    assert result.iloc[0]["native"] == "KEGG:04210"
    assert result.iloc[1]["native"] == "GO:0006915"


@patch("urllib.request.urlopen")
def test_empty_result_returns_empty_dataframe_with_columns(mock_urlopen):
    mock_urlopen.return_value = _FakeResponse({"result": []})

    result = run_gprofiler_enrichment(["TP53"])

    assert result.empty
    assert list(result.columns) == [
        "source", "native", "name", "p_value", "term_size", "query_size", "intersection_size",
    ]


@patch("urllib.request.urlopen")
def test_request_payload_includes_genes_and_organism(mock_urlopen):
    mock_urlopen.return_value = _FakeResponse({"result": []})

    run_gprofiler_enrichment(["TP53", "EGFR", "TP53"], organism="mmusculus", sources=["GO:BP"])

    sent_request = mock_urlopen.call_args[0][0]
    payload = json.loads(sent_request.data.decode("utf-8"))
    assert payload["organism"] == "mmusculus"
    assert payload["query"] == ["EGFR", "TP53"]  # deduplicated and sorted
    assert payload["sources"] == ["GO:BP"]


@patch("urllib.request.urlopen")
def test_network_failure_raises_connection_error(mock_urlopen):
    mock_urlopen.side_effect = urllib.error.URLError("timed out")

    with pytest.raises(ConnectionError, match="g:Profiler"):
        run_gprofiler_enrichment(["TP53"])


@patch("urllib.request.urlopen")
def test_http_error_raises_runtime_error(mock_urlopen):
    mock_urlopen.side_effect = urllib.error.HTTPError(
        url="https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
        code=400, msg="Bad Request", hdrs=None, fp=None,
    )

    with pytest.raises(RuntimeError, match="HTTP 400"):
        run_gprofiler_enrichment(["TP53"], organism="not_a_real_organism")
