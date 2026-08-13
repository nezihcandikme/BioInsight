# OmicForge

A small Python pipeline for bulk RNA-seq: validation, sample QC, normalization, differential expression, gene annotation, pathway enrichment, and plots — as a library or a CLI. v0.10.0.

Not a replacement for DESeq2 or edgeR — checked against both (see [Validation](#validation)), explicit about where it falls short.

## Quick start

```bash
git clone https://github.com/nezihcandikme/BioInsight.git
cd BioInsight
pip install -e .
```

CLI, on a raw count matrix (genes as rows, samples as columns):

```bash
omicforge counts.csv \
  --group1 control_1,control_2,control_3 \
  --group2 treated_1,treated_2,treated_3 \
  --min-count 10
```

Writes `differential_expression.csv`, `qc.csv`, a volcano plot, and a PCA plot to `omicforge_output/`.

As a library, for the intermediate objects instead of files:

```python
from omicforge.pipeline import run_analysis

results = run_analysis(
    "counts.csv",
    group_1=["control_1", "control_2", "control_3"],
    group_2=["treated_1", "treated_2", "treated_3"],
    min_count=10,
)

de = results["differential_expression"]  # log fold change, p-values, significance
```

Every stage (`validate_counts`, `run_sample_qc`, `compute_cpm`, `run_differential_expression`, ...) is also importable on its own.

## Why

I'm a high-school student learning statistics and computational biology by building something slightly beyond what I already know. I wanted to understand the actual mechanics between "here's a count matrix" and "these genes changed" — what makes a count matrix valid, why normalization matters, why testing thousands of genes at once needs correction, why significance isn't the same as biological importance.

The code running is step one. The numbers meaning what I think they mean is the actual objective. The dated account of what changed, why, and what didn't work is in [`DEVLOG.md`](DEVLOG.md).

## What it does

```
load → validate → QC → normalize → (filter) → differential expression → (annotate) → (enrich) → plot
```

- **Validation**: non-negative integer counts, no duplicate IDs, no missing values, a transposed-matrix heuristic — specific error messages, checked first.
- **Sample QC**: library size, genes detected, MAD-based outlier flags.
- **Normalization**: CPM + log2 scale, with optional low-count gene filtering first.
- **Differential expression**: Welch's t-test per gene (default), or an opt-in empirical-Bayes moderated t-test that borrows variance across genes (`method="moderated"`). Benjamini-Hochberg correction either way.
- **Gene annotation** (optional): attach symbols from a local `gene_id,gene_symbol` CSV. No live API; unmapped genes get `NaN`, not a guess.
- **Pathway enrichment** (optional): local hypergeometric test against a GMT file and an explicit background, or a live g:Profiler query (GO/KEGG/Reactome/WikiPathways) if you have internet access.
- **GEO acquisition** (optional): list/download a GEO series' supplementary files and parse its sample metadata from a `GSE` accession.
- **Plots**: volcano, PCA, and an enrichment dot plot (significance vs. gene overlap) when enrichment was run — all from the actual result objects.
- **Plain-language summary** (optional, off by default): an LLM narrates a result table — doesn't verify anything, isn't the point of this project. Requires `ANTHROPIC_API_KEY`.

## Validation

Checked against real DESeq2 and real edgeR (own default settings) on two independent public datasets: [`airway`](https://bioconductor.org/packages/release/data/experiment/html/airway.html) (human) and [`pasilla`](https://bioconductor.org/packages/release/data/experiment/html/pasilla.html) (*Drosophila*) — different organisms, labs, and designs. Methodology, caveats, and reproduction scripts: [`benchmarks/`](benchmarks/).

| dataset | method | precision (DESeq2 / edgeR) | recall (DESeq2 / edgeR) | log2FC Pearson r (DESeq2 / edgeR) |
|---|---|---:|---:|---:|
| airway | welch | 1.0 / 1.0 | 0.075 / 0.098 | 0.938 / 0.947 |
| airway | moderated | 1.0 / 1.0 | 0.165 / 0.217 | 0.938 / 0.947 |
| pasilla | welch | 1.0 / 1.0 | 0.085 / 0.104 | 0.974 / 0.975 |
| pasilla | moderated | 1.0 / 1.0 | 0.174 / 0.215 | 0.974 / 0.975 |

Precision: of the genes OmicForge calls significant, the fraction DESeq2/edgeR also call significant. Recall: the reverse. No false-positive calls relative to the reference DESeq2/edgeR calls were observed in these two benchmarks — not the same claim as a universal zero false-positive rate, and two datasets don't establish one. Welch's t-test tests each gene in isolation and finds a fraction of what DESeq2/edgeR find, since those tools borrow statistical power across genes with similar expression and a plain per-gene test doesn't. The moderated method closes some of that gap the same way, roughly doubling recall at the same precision — a different variance assumption, not a strictly better default.

## Limitations

- Bulk RNA-seq count matrices only. No single-cell, no other omics types, despite the name's broader ambition.
- Default DE method (Welch's t-test) has measurably lower recall than DESeq2/edgeR on every dataset checked so far — see the table above. Expected and explained, not a bug, but real effects get missed.
- Two benchmark datasets is a real check, not a large one. Other tissues, organisms, or designs (unbalanced groups, batch effects, more than two conditions) haven't been measured.
- Two-group comparisons only — no covariates, paired samples, or multi-factor models.
- Local pathway enrichment has no gene-set-size or overlap correction beyond Benjamini-Hochberg; only as good as the background you give it.
- The LLM layer narrates numbers, doesn't verify them — a summary, not an additional check.

## Roadmap

**Near term**: strengthen bulk RNA-seq validation, more benchmark datasets, better experimental-design/metadata handling, better enrichment reporting.

**Later**: additional transcriptomic workflows; broader omics support where scientifically justified, not by default.

## Development

`DEVLOG.md` has the day-by-day account: what changed, what broke, what got rejected, and why. Git commit messages say *what* changed; that file is for *why*.
