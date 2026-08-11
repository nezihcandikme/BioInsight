BioInsight

BioInsight is a Python pipeline for taking an RNA-seq count matrix from “okay, I have a CSV” to actual exploratory results: validation, sample QC, normalization, differential expression, PCA and volcano plots, pathway enrichment, and—if you explicitly ask for it—an LLM-generated explanation of the output.

I’m building it because I wanted to understand what happens between getting biological data and claiming that it means something. Turns out there are approximately seventeen ways to produce convincing nonsense before breakfast, so BioInsight is currently focused on making those mistakes loud, testable, and difficult to ignore.

Current version: v0.3.2. It works, it has tests, and it is actively improving. It is also an educational project—not a production replacement for DESeq2, edgeR, or an actual bioinformatician who has seen your experimental design.

So what is this, exactly?

BioInsight started with one extremely glamorous problem: validating a CSV.

Then the CSV raised questions.

Are these really non-negative integer counts?

Are genes rows and samples columns, or did someone rotate reality by 90 degrees?

Does one sample have dramatically fewer reads than the others?

How do I compare samples with different sequencing depths?

Which genes differ between conditions?

If I test thousands of genes, how many “discoveries” are just statistical noise wearing a lab coat?

Do those genes converge on any interesting pathways?

Each question became a module. The result is a small, testable RNA-seq pipeline that follows the same broad reasoning order as an analysis: load → validate → inspect → normalize → compare → visualize → interpret.

The code running is step one. The numbers meaning what I think they mean is the actual objective.

Why I’m building it

I’m a high-school student learning programming, statistics, and computational biology by building things that are slightly beyond what I currently know how to build.

Before BioInsight, I understood RNA-seq in the dangerously comfortable summary-version: count reads per gene, compare two conditions, find changed genes. That explanation is technically related to reality, but it skips nearly everything capable of ruining an analysis.

Building the pipeline forced me to confront those missing layers directly: what makes a count matrix valid, why library size matters, why normalization is not one universal operation, why a difference of means only becomes a log fold change on the right scale, why thousands of simultaneous tests require correction, and why statistical significance is not the same thing as biological importance.

BioInsight is the record of learning those lessons in code—one function, one broken test, and one increasingly specific error message at a time.

The pipeline

raw counts (CSV or DataFrame)
          │
          ▼
      validation  ─────────────  bioinsight.io.counts
          │
          ▼
      sample QC  ──────────────  bioinsight.qc.metrics
          │
          ▼
      CPM + log2  ─────────────  bioinsight.normalization.methods
          │
          ▼
      exploratory DE  ─────────  bioinsight.differential_expression.methods
          │
          ├──▶ volcano + PCA  ─  bioinsight.visualization.plots
          │
          ├──▶ enrichment  ─────  bioinsight.pathway_analysis.methods
          │
          └──▶ AI summary  ─────  bioinsight.ai_explanation.methods

bioinsight.pipeline.run_analysis() connects the full workflow. Every module can also be called independently when you want to inspect an intermediate result, replace a step, or figure out exactly where the numbers went off the rails.

Installation

BioInsight requires Python 3.9 or newer.

git clone https://github.com/nezihcandikme/BioInsight.git
cd BioInsight
pip install -e .

The core analysis uses pandas, NumPy, SciPy, statsmodels, Matplotlib, and scikit-learn.

Optional AI explanations

explain_de_results() and run_analysis(..., explain_results=True) can send a compact summary of the differential-expression table to the Anthropic API.

Add your key to a .env file in the project root:

ANTHROPIC_API_KEY=your_key_here

This layer does not perform the statistics, verify the experiment, or magically discover biology. It narrates results that BioInsight has already computed. A confident paragraph generated from a weak analysis is still a weak analysis—just with better punctuation.

Quick start

from bioinsight.pipeline import run_analysis

results = run_analysis(
    counts="counts.csv",                 # CSV path or pandas DataFrame
    group_1=["control_1", "control_2"],
    group_2=["treated_1", "treated_2"],
    gene_sets={                          # optional pathway enrichment
        "pathway_1": {"ENSG001", "ENSG003"},
        "pathway_2": {"ENSG002", "ENSG004"},
    },
    background_genes={"ENSG001", "ENSG002", "ENSG003", "ENSG004"},
    generate_plots=True,
    explain_results=False,
)

results["qc"]
results["normalized"]
results["differential_expression"]
results["pathway_enrichment"]

results["volcano_fig"].savefig("volcano.png")
results["pca_fig"].savefig("pca.png")

Want to see the gears instead of pressing the large pipeline button?

from bioinsight.io.counts import load_count_matrix
from bioinsight.qc.metrics import run_sample_qc
from bioinsight.normalization.methods import compute_cpm, log2_transform
from bioinsight.differential_expression.methods import run_differential_expression
from bioinsight.visualization.plots import plot_volcano

counts = load_count_matrix("counts.csv")
qc = run_sample_qc(counts)

normalized = log2_transform(compute_cpm(counts))
de_results = run_differential_expression(
    normalized,
    group_1=["control_1", "control_2"],
    group_2=["treated_1", "treated_2"],
)

fig = plot_volcano(de_results)

Input format

BioInsight expects a raw count matrix supplied as a CSV file or pandas.DataFrame:

Genes are rows.

Samples are columns.

The first CSV column contains unique gene identifiers and becomes the index.

Counts are non-negative integers.

Missing values are not allowed.

gene_id,control_1,control_2,treated_1,treated_2
ENSG001,10,15,8,120
ENSG002,0,3,1,98
ENSG003,120,98,140,15
ENSG004,5,4,0,25

validate_counts() rejects non-integer counts, negative values, missing values, and duplicated gene IDs. It also warns when the matrix shape looks suspiciously transposed. That check is heuristic, because a DataFrame cannot explain its own experimental design no matter how intensely we stare at it.

Modules

Module

Responsibility

bioinsight.io.counts

Load count matrices, validate their structure, and fail with specific errors.

bioinsight.qc.metrics

Calculate library size and detected genes; flag sample-level outliers using MAD.

bioinsight.normalization.methods

Compute CPM and apply log2(x + 1).

bioinsight.differential_expression.methods

Compute mean-difference log fold changes, Welch’s t-tests, and Benjamini–Hochberg correction.

bioinsight.visualization.plots

Generate volcano and PCA plots.

bioinsight.pathway_analysis.methods

Run hypergeometric pathway enrichment against a defined background universe.

bioinsight.ai_explanation.methods

Optionally summarize an existing DE table through the Anthropic API.

bioinsight.pipeline

Connect the modules through run_analysis().

Statistical reality check

This is the section to read before trusting a beautiful volcano plot.

CPM corrects for sequencing depth—only sequencing depth. It does not correct gene length, RNA-composition effects, batch effects, or every other unpleasant thing hiding inside an experiment. It is not TPM, and it is not DESeq2’s median-of-ratios normalization.

log_fold_change is calculated as mean(group_1) - mean(group_2). That represents a log2 fold change only when the input is already on a log2 scale. On raw or linear-scale data it is simply a difference of means with an inaccurately exciting name.

Differential expression currently uses a per-gene Welch’s t-test. This is an exploratory comparison, not a count-aware negative-binomial model. It does not model the RNA-seq mean–variance relationship or biological dispersion the way DESeq2 and edgeR do. Use it to investigate data, not to support clinical claims.

Benjamini–Hochberg correction is essential. Test enough genes and random chance will happily manufacture “significant” results for you. Multiple-testing correction controls that problem; it does not eliminate every possible false discovery.

Enrichment depends on the background universe. background_genes should contain the genes that could actually have been selected by the analysis—usually the tested genes that passed filtering. Using every known gene in existence answers a different statistical question.

A p-value is evidence under a model, not a biological verdict. BioInsight can compute and organize evidence. It cannot rescue a weak experimental design or decide what a biological result means in context.

Current limitations

BioInsight is useful enough to run and early enough to distrust responsibly.

The DE implementation has not yet been benchmarked against DESeq2 or edgeR on a real public dataset.

There is no independent low-count filtering step before differential expression.

Batch effects are not modeled.

The significance thresholds—adjusted p < 0.05 and |log_fold_change| > 1—are currently hardcoded.

Pathway enrichment expects Python dict/set inputs; there is no GMT, MSigDB, or Enrichr loader yet.

The AI explanation layer summarizes the existing table but does not independently validate it.

BioInsight currently understands one scientific-data modality: bulk RNA-seq count matrices. The long-term goal is broader, but pretending it can analyze “any dataset” today would be an excellent way to make the README more advanced than the software.

Development log

The repository grew quickly because every completed step exposed the next missing one.

August 6 — the CSV era. Created the package structure and count-matrix validation: integers, non-negativity, unique gene IDs, missing values, and an exception hierarchy. Apparently “read the file” was already several functions.

August 7 — the data has entered the building. Added CSV loading, library size, detected-gene counts, MAD-based sample outliers, and CPM normalization. This was the point where different sequencing depths stopped being an abstract caveat and became something the code had to handle.

August 8 — statistics have consequences. Added log2(x + 1), mean-difference effect sizes, and Welch’s t-tests. Also discovered that naming a number correctly requires knowing what scale produced it. Rude but fair.

August 9 — thousands of p-values appeared. Added Benjamini–Hochberg correction, minimum sample-size validation, volcano plots, PCA, and hypergeometric enrichment. BioInsight officially became too large to describe as “the CSV thing.”

August 11 — make it harder to fool. Repaired package initialization, restricted enrichment calculations to the tested background, guarded zero-count and zero-variance edge cases, connected everything through run_analysis(), expanded the test suite, added more specific error messages, and rewrote the documentation to state what the pipeline does—and what it absolutely does not do.

Roadmap

The long-term idea is larger than RNA-seq: a modality-aware system that can inspect scientific data, determine what kind of analysis is appropriate, ask for missing metadata, run validated methods, and produce a reproducible report.

That goal is extremely far from “upload arbitrary CSV, receive truth,” so the next steps stay concrete:

Benchmark BioInsight’s DE output against DESeq2 and edgeR on a real public RNA-seq dataset.

Add independent low-count filtering before differential expression.

Make significance and effect-size thresholds configurable.

Add GMT/MSigDB gene-set loading.

Add optional Enrichr or g enrichment, gene-annotation lookup, and GEO/SRA retrieval.

Separate dataset understanding, analysis planning, deterministic computation, and explanation into explicit layers.

Add new scientific-data modalities one validated adapter at a time.

Current mission: make the statistical layer harder to fool before teaching it new tricks.

Tests

pip install -e ".[dev]"
pytest

BioInsight currently has 56 tests covering both the normal workflow and the cases that tend to become 2 a.m. debugging sessions: zero-count samples, constant-expression genes, NaN propagation, missing or duplicated sample names, overlapping groups, empty gene sets, empty backgrounds, invalid PCA inputs, and plotting edge cases.

GitHub Actions runs the suite on Python 3.10–3.12 after every push. Green tests do not prove the biology is correct, but red tests are at least considerate enough to tell us something is definitely wrong.
