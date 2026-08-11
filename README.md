# BioInsight

BioInsight is a long-term project that aims to make bioinformatics analysis more accessible through the use of artificial intelligence.

The idea behind the project is not to replace existing bioinformatics tools, but to build an intelligent layer on top of them. Instead of requiring researchers to manually perform complex analyses and interpret dozens of outputs, BioInsight aims to gradually become capable of understanding biological datasets, assisting with data analysis, and helping users interpret their results.

## Why?

Modern biology produces enormous amounts of data, yet many analyses remain difficult to perform without programming experience.

I started BioInsight because I want to help researchers focus on understanding biology instead of struggling with computational tools. At the same time, this project is my own journey into computational biology, bioinformatics, software engineering, and artificial intelligence.

## Status

BioInsight is an early-stage, actively developed project (v0.3). It currently supports:

- Count-matrix loading and validation
- Sample QC and MAD-based outlier detection
- CPM normalization and log2 transformation
- A basic, exploratory differential expression analysis (mean difference + Welch's t-test, Benjamini-Hochberg correction)
- Volcano and PCA plots
- Hypergeometric pathway enrichment
- An optional AI explanation layer (Anthropic Claude) that summarizes DE results in plain language

**A note on the differential expression method:** BioInsight's built-in DE analysis is a simple mean-difference t-test on normalized, log-transformed expression values. It does **not** model count dispersion or the mean-variance relationship the way DESeq2 or edgeR do, and it is not a substitute for those tools. Treat its output as a quick, exploratory look at your data, not a publication-ready differential expression result.

## Installation

Requires Python >= 3.9.

```bash
git clone https://github.com/nezihcandikme/BioInsight.git
cd BioInsight
pip install -e .
```

This installs the core dependencies (pandas, numpy, scipy, statsmodels, matplotlib, scikit-learn) along with `anthropic` and `python-dotenv`, which are only needed if you use the optional AI explanation feature.

### Optional: AI explanations

To use `explain_de_results` / `run_analysis(..., explain_results=True)`, set an `ANTHROPIC_API_KEY` in a `.env` file at the project root:

```
ANTHROPIC_API_KEY=your_key_here
```

## Input format

BioInsight expects a **raw count matrix** as a CSV file (or `pandas.DataFrame`) with:

- Genes as rows, identified by a unique gene ID in the first column (used as the index)
- Samples as columns, with non-negative integer counts, no missing values

Example (`counts.csv`):

```csv
gene_id,sample_A,sample_B,sample_C,sample_D
ENSG001,10,15,8,120
ENSG002,0,3,1,98
ENSG003,120,98,140,15
ENSG004,5,4,0,25
```

## Usage

### Full pipeline

```python
from bioinsight.pipeline import run_analysis

results = run_analysis(
    counts="counts.csv",                 # path to a CSV, or a pandas DataFrame
    group_1=["sample_A", "sample_B"],
    group_2=["sample_C", "sample_D"],
    gene_sets={                          # optional: pathway enrichment
        "pathway_1": {"ENSG001", "ENSG003"},
        "pathway_2": {"ENSG002", "ENSG004"},
    },
    background_genes={"ENSG001", "ENSG002", "ENSG003", "ENSG004"},
    generate_plots=True,
    explain_results=False,               # set True to also get an AI summary
)

results["qc"]                     # per-sample QC metrics
results["normalized"]             # CPM-normalized, log2-transformed matrix
results["differential_expression"]  # log fold change, p-value, adjusted p-value, significant
results["volcano_fig"].savefig("volcano.png")
results["pca_fig"].savefig("pca.png")
results["pathway_enrichment"]     # enrichment p-values per gene set
```

### Using individual modules

You can also call each step directly if you want more control:

```python
from bioinsight.io.counts import load_count_matrix
from bioinsight.qc.metrics import run_sample_qc
from bioinsight.normalization.methods import compute_cpm, log2_transform
from bioinsight.differential_expression.methods import run_differential_expression
from bioinsight.visualization.plots import plot_volcano, plot_pca

counts = load_count_matrix("counts.csv")
qc = run_sample_qc(counts)

normalized = log2_transform(compute_cpm(counts))
de_results = run_differential_expression(
    normalized, group_1=["sample_A", "sample_B"], group_2=["sample_C", "sample_D"]
)

fig = plot_volcano(de_results)
```

## Testing

```bash
pip install -e ".[dev]"  # or just: pip install pytest
pytest
```

## Roadmap

1. Validate DE output against DESeq2/edgeR on a real public dataset
2. Continuous integration (automated tests on every push)
3. Optional integrations for pathway databases (Enrichr, g:Profiler), gene annotation, and public dataset retrieval (GEO/SRA)
