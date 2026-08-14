"""
Fetching data from NCBI GEO by series accession (e.g. "GSE52778").

Scope, deliberately narrow: this finds and downloads the supplementary
files a GEO series' submitter already uploaded (frequently a processed
count matrix, sometimes per-sample files), and parses the series matrix
file's own sample metadata into a DataFrame. It does **not** do anything
with raw sequencing reads (FASTQ/SRA) -- turning those into counts means
alignment/quantification, a different and much heavier job that belongs to
a dedicated pipeline tool (this project's own environment has an
nf-core/rnaseq-based skill for exactly that), not this small module.
Whether a given series even has a usable pre-computed count matrix in its
supplementary files varies by submitter -- this module does the tedious,
uniform part (finding the right NCBI FTP path, listing what's there,
downloading it, parsing the metadata) and leaves the "is this file actually
a count matrix, and does it need reshaping" judgment call to the caller,
because that part genuinely isn't uniform across GEO submissions.

Every function here makes a live network call to
``ftp.ncbi.nlm.nih.gov`` and needs outbound internet access -- confirmed
unavailable in this project's own development sandbox (see
``deconcord.pathway_analysis.live_enrichment`` for the same caveat, and
the reason both are opt-in, self-contained modules rather than something
``run_analysis`` calls automatically).
"""

import gzip
import os
import re
import shutil
import urllib.error
import urllib.request

import pandas as pd

_GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"

_GSE_PATTERN = re.compile(r"GSE\d+")


def _series_dir_url(gse_id: str) -> str:
    """
    Build the NCBI FTP directory URL for a GEO series, e.g. "GSE52778" ->
    ".../series/GSE52nnn/GSE52778" -- GEO buckets series into directories
    of 1000 by truncating the last 3 digits of the accession number and
    replacing them with "nnn".
    """
    if not _GSE_PATTERN.fullmatch(gse_id):
        raise ValueError(
            f"'{gse_id}' doesn't look like a GEO series accession -- expected "
            "something like 'GSE52778' (a GEO sample accession like 'GSM...' "
            "or a GEO dataset accession like 'GDS...' won't work here)."
        )
    digits = gse_id[3:]
    truncated = digits[:-3] if len(digits) > 3 else ""
    bucket = f"GSE{truncated}nnn"
    return f"{_GEO_FTP_BASE}/{bucket}/{gse_id}"


def _wrap_network_error(exc: urllib.error.URLError, url: str) -> Exception:
    if isinstance(exc, urllib.error.HTTPError):
        return RuntimeError(
            f"GEO returned HTTP {exc.code} for {url}: {exc.reason}. "
            "This usually means the accession doesn't exist or has no data "
            "at this path, not a network problem."
        )
    return ConnectionError(
        f"Could not reach GEO ({url}): {exc.reason}. This function needs "
        "outbound internet access to ftp.ncbi.nlm.nih.gov; it will not "
        "work in a network-restricted environment."
    )


def list_geo_supplementary_files(gse_id: str, timeout: float = 30.0) -> list[str]:
    """
    List the supplementary files a GEO series' submitter uploaded.

    Parameters
    ----------
    gse_id : str
        A GEO series accession, e.g. ``"GSE52778"``.
    timeout : float, optional
        Request timeout in seconds. Default 30.

    Returns
    -------
    list[str]
        Filenames available in the series' ``suppl/`` directory (not full
        URLs -- pass one of these to ``download_geo_supplementary_file``).
        Empty if the series has no supplementary files at all.

    Raises
    ------
    ValueError
        If ``gse_id`` isn't a well-formed GEO series accession.
    ConnectionError
        If GEO can't be reached at all.
    RuntimeError
        If GEO is reachable but returns a non-2xx response (e.g. the
        accession doesn't exist).
    """
    url = f"{_series_dir_url(gse_id)}/suppl/"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            html = response.read().decode("utf-8", errors="replace")
    except urllib.error.URLError as e:
        raise _wrap_network_error(e, url) from e

    hrefs = re.findall(r'href="([^"]+)"', html)
    return sorted({h for h in hrefs if not h.endswith("/") and not h.startswith("?")})


def download_geo_supplementary_file(gse_id: str, filename: str, dest_dir: str = ".", timeout: float = 120.0) -> str:
    """
    Download one of a GEO series' supplementary files to a local directory.

    Parameters
    ----------
    gse_id : str
        A GEO series accession, e.g. ``"GSE52778"``.
    filename : str
        A filename as returned by ``list_geo_supplementary_files``.
    dest_dir : str, optional
        Local directory to write the file into (created if it doesn't
        exist). Default the current directory.
    timeout : float, optional
        Request timeout in seconds. Default 120 (supplementary files can
        be large).

    Returns
    -------
    str
        The local path the file was written to.

    Raises
    ------
    ConnectionError
        If GEO can't be reached at all.
    RuntimeError
        If GEO is reachable but returns a non-2xx response (e.g. the
        filename doesn't exist for this series).
    """
    url = f"{_series_dir_url(gse_id)}/suppl/{filename}"
    os.makedirs(dest_dir, exist_ok=True)
    dest_path = os.path.join(dest_dir, filename)

    try:
        with urllib.request.urlopen(url, timeout=timeout) as response, open(dest_path, "wb") as out_file:
            shutil.copyfileobj(response, out_file)
    except urllib.error.URLError as e:
        raise _wrap_network_error(e, url) from e

    return dest_path


def fetch_geo_sample_metadata(gse_id: str, timeout: float = 30.0) -> pd.DataFrame:
    """
    Fetch and parse a GEO series' sample metadata (the ``!Sample_*`` lines
    of its series matrix file) into a DataFrame -- useful for figuring out
    which samples belong to which experimental group before calling
    ``deconcord.pipeline.run_analysis``.

    Parameters
    ----------
    gse_id : str
        A GEO series accession, e.g. ``"GSE52778"``.
    timeout : float, optional
        Request timeout in seconds. Default 30.

    Returns
    -------
    pd.DataFrame
        One row per sample, indexed by GEO sample accession (GSM...).
        Always includes a ``title`` column if present; one column per
        distinct ``!Sample_characteristics_ch1`` field (named from the
        "key: value" text GEO submitters write, e.g. a line reading
        "tissue: lung" becomes a ``tissue`` column). Column sets vary by
        series, since submitters choose their own characteristic fields.

    Raises
    ------
    ValueError
        If ``gse_id`` isn't a well-formed GEO series accession, or if the
        downloaded series matrix has no recognizable sample accession
        line at all.
    ConnectionError
        If GEO can't be reached at all.
    RuntimeError
        If GEO is reachable but returns a non-2xx response.
    """
    url = f"{_series_dir_url(gse_id)}/matrix/{gse_id}_series_matrix.txt.gz"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            raw = response.read()
    except urllib.error.URLError as e:
        raise _wrap_network_error(e, url) from e

    text = gzip.decompress(raw).decode("utf-8", errors="replace")

    accessions: list[str] = []
    fields: dict[str, list[str]] = {}

    for line in text.splitlines():
        if not line.startswith("!Sample_"):
            continue

        parts = line.split("\t")
        label = parts[0].lstrip("!")
        values = [v.strip().strip('"') for v in parts[1:]]

        if label == "Sample_geo_accession":
            accessions = values
        elif label == "Sample_title":
            fields["title"] = values
        elif label == "Sample_characteristics_ch1":
            if values and ":" in values[0]:
                column = values[0].split(":", 1)[0].strip()
            else:
                column = f"characteristic_{len(fields) + 1}"
            # A series can repeat this line once per characteristic;
            # guard against two characteristics happening to share a name.
            unique_column = column
            suffix = 2
            while unique_column in fields:
                unique_column = f"{column}_{suffix}"
                suffix += 1
            fields[unique_column] = [v.split(":", 1)[1].strip() if ":" in v else v for v in values]

    if not accessions:
        raise ValueError(
            f"Could not find a '!Sample_geo_accession' line in {gse_id}'s "
            "series matrix file -- unexpected format for this series."
        )

    metadata = pd.DataFrame(fields, index=pd.Index(accessions, name="sample"))
    return metadata
