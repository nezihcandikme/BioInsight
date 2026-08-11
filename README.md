# BioInsight

BioInsight is a small RNA-seq analysis pipeline, built in Python, that takes a raw gene count matrix and walks it through validation, quality control, normalization, exploratory differential expression, visualization, and pathway enrichment — with an optional layer that asks an LLM to summarize the results in plain language.

It is a personal, actively developed project (currently v0.3.2), not a maintained scientific tool. If you're evaluating whether to depend on it: don't, yet. If you want to read clean-ish code that implements the standard early steps of an RNA-seq workflow, or you're learning the same things I was learning while writing it, it might be genuinely useful.

## What BioInsight is

BioInsight is what happens when "let's just validate a CSV of gene counts" doesn't stop being interesting. It started as a single function that checked whether a count matrix was well-formed, and it kept growing because every step of a real RNA-seq analysis raised the next question: now that the matrix is valid, is it any good? Now that I know it's decent, how do I make samples with different sequencing depths comparable? Now that they're comparable, which genes actually differ between conditions — and how do I stop myself from lying to myself about which differences are real?

Each answer became a module. The result is a working, testable, incrementally-built pipeline that mirrors the actual order a bioinformatician thinks in: load → check → normalize → test → visualize → interpret.

## Why I started building it

I'm a highschool student, and RNA-seq analysis was one of those things I understood in outline — "you count reads per gene, compare conditions, get a list of genes that changed" — without understanding any of the machinery underneath. The fastest way I know how to actually learn something is to implement it badly, watch it break, and fix it until it stops lying to me.

So BioInsight began as an exercise in learning how a count matrix moves through a pipeline: what "valid" even means for one, why library size differs sample to sample and what that does to your comparisons, why a plain difference of means isn't a fold change unless you're on a log scale, why testing thousands of genes at once makes p-values misleading unless you correct for it, and why "the code ran without an error" and "the biology this represents makes sense" are two completely different bars to clear. This project is the record of working through those questions in order, one module at a time.

## Current pipeline

```
raw counts (CSV)
      │
      ▼
  validate  ──────────  bioinsight.io.counts
      │
      ▼
  sample QC  ─────────  bioinsight.qc.metrics
      │
      ▼
  CPM + log2  ─────────  bioinsight.normalization.methods
      │
      ▼
  differential expression ── bioinsight.differential_expression.methods
      │
      ├──▶ volcano / PCA plots ── bioinsight.visualization.plots
      │
      ├──▶ pathway enrichment ── bioinsight.pathway_analysis.methods
      │
      └──▶ plain-language summary (optional) ── bioinsight.ai_explanation.methods
```

`bioinsight.pipeline.run_analysis()` runs the whole thing in one call. Every step is also usable on its own if you want more control (or you're studying how one piece works).

## Installation

Requires Python >= 3.9.

```bash
git clone https://github.com/nezihcandikme/BioInsight.git
cd BioInsight
pip install -e .
```

This pulls in pandas, numpy, scipy, statsmodels, matplotlib, and scikit-learn for the core pipeline, plus `anthropic` and `python-dotenv` for the optional AI summary layer.

### Optional: AI-generated summaries

`explain_de_results()` / `run_analysis(..., explain_results=True)` calls the Anthropic API. Set a key in a `.env` file at the project root:

```
ANTHROPIC_API_KEY=your_key_here
```

Nothing else in the pipeline needs this — it's purely a "narrate my results" convenience on top of numbers that are already computed.

## Quick-start example

```python
from bioinsight.pipeline import run_analysis

results = run_analysis(
    counts="counts.csv",                 # path to a CSV, or a pandas DataFrame
    group_1=["control_1", "control_2"],
    group_2=["treated_1", "treated_2"],
    gene_sets={                          # optional: pathway enrichment
        "pathway_1": {"ENSG001", "ENSG003"},
        "pathway_2": {"ENSG002", "ENSG004"},
    },
    background_genes={"ENSG001", "ENSG002", "ENSG003", "ENSG004"},
    generate_plots=True,
    explain_results=False,               # set True to also get an AI summary
)

results["qc"]                       # per-sample QC metrics
results["normalized"]               # CPM-normalized, log2-transformed matrix
results["differential_expression"]  # log_fold_change, p_value, adjusted_p_value, significant
results["volcano_fig"].savefig("volcano.png")
results["pca_fig"].savefig("pca.png")
results["pathway_enrichment"]       # pathway, p_value, adjusted_p_value, overlap_count
```

Or call the modules individually if you want to inspect each intermediate step:

```python
from bioinsight.io.counts import load_count_matrix
from bioinsight.qc.metrics import run_sample_qc
from bioinsight.normalization.methods import compute_cpm, log2_transform
from bioinsight.differential_expression.methods import run_differential_expression
from bioinsight.visualization.plots import plot_volcano

counts = load_count_matrix("counts.csv")
qc = run_sample_qc(counts)

normalized = log2_transform(compute_cpm(counts))
de_results = run_differential_expression(
    normalized, group_1=["control_1", "control_2"], group_2=["treated_1", "treated_2"]
)

fig = plot_volcano(de_results)
```

## Input format

A raw count matrix as a CSV (or `pandas.DataFrame`):

- Genes as rows, identified by a unique gene ID in the first column (used as the index).
- Samples as columns, with non-negative integer counts and no missing values.

```csv
gene_id,control_1,control_2,treated_1,treated_2
ENSG001,10,15,8,120
ENSG002,0,3,1,98
ENSG003,120,98,140,15
ENSG004,5,4,0,25
```

`validate_counts` / `load_count_matrix` will reject non-integer values, negative values, missing values, and duplicate gene IDs, and will warn (not fail) if the matrix's shape looks like it might be transposed.

## Available modules

| Module | What it does |
|---|---|
| `bioinsight.io.counts` | Load and validate a raw count matrix; exception hierarchy for specific validation failures. |
| `bioinsight.qc.metrics` | Library size, genes detected, and MAD-based outlier flagging per sample. |
| `bioinsight.normalization.methods` | CPM normalization and log2(x + 1) transformation. |
| `bioinsight.differential_expression.methods` | Mean-difference "log fold change," per-gene Welch's t-test, Benjamini-Hochberg correction. |
| `bioinsight.visualization.plots` | Volcano plot and PCA plot. |
| `bioinsight.pathway_analysis.methods` | Hypergeometric pathway/gene-set enrichment against a defined background. |
| `bioinsight.ai_explanation.methods` | Optional plain-language summary of DE results via the Anthropic API. |
| `bioinsight.pipeline` | `run_analysis()` — chains all of the above into one call. |

## Statistical assumptions

Worth reading before trusting any output:

- **CPM corrects for sequencing depth only.** It does not correct for gene length, RNA composition effects, or other technical/compositional biases. It is not TPM, and it is not DESeq2's median-of-ratios normalization.
- **`log_fold_change` is `mean(group_1) - mean(group_2)`.** This is a genuine log2 fold change *only* if the input is already on a log2 scale (which is why `run_differential_expression` expects log2-CPM input, not raw counts). Feed it raw or linear-scale data and you get a difference of means — a real, computable number, just not the quantity its name implies.
- **Differential expression uses Welch's t-test per gene**, not a count-based negative-binomial model. It does not model biological dispersion or the count-specific mean-variance relationship that tools like DESeq2 and edgeR are built around. It's a legitimate exploratory method; it is not a substitute for those tools when the result needs to support a real biological or clinical claim.
- **Multiple-testing correction (Benjamini-Hochberg) is not optional.** Testing thousands of genes without it means a meaningful fraction of "significant" genes are noise, regardless of what the biology actually did.
- **Pathway enrichment is only as correct as the background you give it.** The hypergeometric test's answer depends entirely on `background_genes` being the actual tested universe (e.g. every gene that passed expression filtering) — not "every gene in the genome." Get the background wrong and the p-values are answering a different question than the one you think you're asking.
- **A p-value is not a biological conclusion.** BioInsight computes numbers; it doesn't know what your experiment means.

## Current limitations

- The DE method has not been benchmarked against DESeq2 or edgeR on a real dataset. Until that happens, treat its output as directionally interesting, not confirmed.
- No independent filtering step before DE (e.g. dropping very-low-count genes), which real pipelines typically do before testing.
- No batch-effect handling.
- `run_differential_expression`'s significance thresholds (`adjusted_p_value < 0.05`, `|log_fold_change| > 1`) are hardcoded, not configurable.
- Pathway enrichment takes gene sets as plain Python `dict`/`set` input — there's no built-in loader for standard formats (GMT, MSigDB, Enrichr) yet.
- The AI explanation layer narrates existing numbers; it doesn't verify them, and it shouldn't be trusted more than the analysis it's describing.

## Development history

Built in a short, dense stretch — the git history reads almost one module a day:

- **Aug 6** — Project scaffolding, then `validate_counts`: the exception hierarchy, integer/negativity/uniqueness/missing-value checks. This is where "what does a *valid* count matrix even mean" got answered in code for the first time.
- **Aug 7** — `load_count_matrix`, then sample QC: library size, genes detected, MAD-based outlier flagging (v0.0.2). CPM normalization followed the same day (v0.0.3 start) — the first time "sample A has more reads than sample B" had to become "and here's how I correct for that."
- **Aug 8** — log2 transformation, then the first differential expression code: log fold change and Welch's t-test. This is roughly where "the code runs" stopped being good enough and "does this number mean what I think it means" took over.
- **Aug 9** — Benjamini-Hochberg correction and a minimum-sample-size check for the t-test (v0.1); volcano and PCA plots (v0.2); hypergeometric pathway enrichment (v0.3).
- **Aug 11** — A full correctness pass (v0.3.1): fixed broken package files, restricted enrichment to the tested background, guarded CPM and the t-test against zero-count and zero-variance edge cases, added `run_analysis()` to tie every module together, and wrote the tests that should have existed for all of the above from the start. Same day, a consistency/documentation pass (v0.3.2): more edge-case guards (empty gene sets, empty backgrounds, PCA with too few samples, NaN leaking through multiple-testing correction), per-sample error messages instead of generic ones, and this README rewrite.

## Roadmap

1. Benchmark DE output against DESeq2/edgeR on a real public dataset — this is the one that actually matters before calling the DE module trustworthy.
2. A pre-DE low-count filtering step.
3. Configurable significance thresholds.
4. GMT/MSigDB gene set loading for pathway enrichment.
5. Optional integrations: Enrichr/g:Profiler for enrichment, gene annotation lookups, GEO/SRA dataset retrieval.

Current mission, in short: make the statistical layer harder to fool before making it do more things.

## Running tests

```bash
pip install -e ".[dev]"
pytest
```

56 tests as of v0.3.2, covering both the happy path and the edge cases that actually break real analyses: zero-count samples, constant-expression genes, NaN propagation through multiple-testing correction, missing/duplicate/overlapping sample names, empty gene sets and empty backgrounds, and a couple of plotting edge cases (too few samples for PCA, a p-value of exactly zero on a volcano plot). CI runs the suite on Python 3.10–3.12 on every push.
