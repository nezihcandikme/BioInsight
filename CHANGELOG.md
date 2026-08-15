# Changelog

All notable changes to DEConcord are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/). Versioning follows
[Semantic Versioning](https://semver.org/). Pre-1.0, breaking changes can
happen in a minor version bump.

This file is user-facing: what changed and what it means for you. For the
reasoning behind each decision, including rejected alternatives, see
[`DEVLOG.md`](DEVLOG.md).

## [0.13.2] — 2026-08-14

### Changed
- Rewrote the prose in `README.md`, `docs/index.html`, `CHANGELOG.md`,
  and `DEVLOG.md` into shorter, plainer sentences and removed em-dashes
  throughout. No factual content changed: same numbers, same decisions,
  same test-count checkpoints.

## [0.13.1] — 2026-08-14

### Fixed
- The GitHub repository was renamed `BioInsight` to `DeConcord`, but every
  URL in `README.md`, `CONTRIBUTING.md`, `SECURITY.md`, `CITATION.cff`,
  `pyproject.toml`, and `docs/index.html` still pointed at the old
  `BioInsight` path, including the README's own `git clone` command.
  Updated all of them. GitHub's redirect from the old name masked this
  for anyone clicking through, but not for anyone reading the raw text.
  `DEVLOG.md`'s historical entries are unaffected; they describe what
  was true when written. `pyproject.toml`'s `Homepage` now also points
  at the live GitHub Pages site instead of the repository itself.

## [0.13.0] — 2026-08-14

Professionalization pass. No new analytical capability. The focus was
making the repository's structure, docs, tests, and packaging match a
project that expects outside users. See `DEVLOG.md` for the full
reasoning.

### Changed
- **Breaking: CLI output layout.** `deconcord_output/` now organizes
  results into `tables/` (CSVs) and `figures/` (PNGs) instead of a flat
  directory, plus a new `run_metadata.json`. Existing scripts that read
  `deconcord_output/differential_expression.csv` directly need to read
  `deconcord_output/tables/differential_expression.csv` instead. Same
  rename for `qc.csv`, `pathway_enrichment.csv`,
  `live_pathway_enrichment.csv`, and the `.png` figures.
- **Package layout**: source moved to `src/deconcord/` (was `deconcord/`
  at the repo root). A `pip install -e .` from a clean checkout is
  unaffected. This only matters if you had tooling pointed at the old
  path.
- `anthropic`/`python-dotenv` are now an optional `[ai]` extra
  (`pip install deconcord[ai]`) instead of core dependencies. Only
  needed for the opt-in `explain_results=True` AI summary. Calling it
  without the extra now raises a clear `ImportError` naming the fix,
  instead of a bare `ModuleNotFoundError`.
- `tests/test_differential_expression.py` (554 lines) split into
  `test_differential_expression_welch.py`, `_moderated.py`, and
  `_covariates.py`. No test content changed, only organization.

### Added
- `deconcord/__init__.py` now exposes a curated public API
  (`import deconcord as dc`; `dc.run_analysis`, `dc.compute_de_concordance`,
  etc.) and `deconcord.__version__`.
- `deconcord --version` CLI flag.
- CLI writes `run_metadata.json` alongside its normal output: DEConcord
  version, Python version, core dependency versions, the exact command,
  an input summary, and the parameters used, for reproducing a run.
- `LICENSE` (MIT), `CITATION.cff`, `CONTRIBUTING.md`, `SECURITY.md`, and
  GitHub issue/PR templates.
- `METHODOLOGY.md`: what each concordance metric computes, its source
  if it's a standard statistic, and its limitations.
- `examples/quickstart.py`: a single runnable script covering the full
  current workflow. Synthetic data, load and validate, two DE methods,
  concordance, plot. No network or R required.
- Test coverage for previously untested edge cases: duplicate gene IDs
  and infinite values in `validate_counts`.

## [0.12.0] — 2026-08-14

### Changed
- **Project renamed from OmicForge to DEConcord**, and scope narrowed
  from "general bulk RNA-seq pipeline" to "robustness/concordance
  analysis for RNA-seq differential expression." Package import path is
  now `deconcord` (was `omicforge`); CLI command is now `deconcord`
  (was `omicforge`). See the README's "Why" section and the Aug 14
  DEVLOG entry for the full reasoning.

### Added
- `deconcord.concordance.methods.compute_de_concordance`: the project's
  new core function. Given two DE result tables (any tool, configurable
  column names), it computes significant-gene overlap, Jaccard index,
  Pearson/Spearman log fold change correlation, directional agreement,
  and explicit concordant, discordant, and method-specific gene lists.
- `benchmarks/method_concordance.py`: runs the new function against
  real, committed DESeq2/edgeR result tables (airway, pasilla).

## [0.11.0] — 2026-08-13

### Added
- Covariate-adjusted differential expression:
  `run_differential_expression_with_covariates`, a per-gene linear model
  that tests group_1 vs group_2 while holding one or more covariates
  (batch, sex, etc.) fixed. Wired into `run_analysis`
  (`metadata`/`covariate_cols`/`moderated_covariates`) and the CLI
  (`--metadata`/`--covariates`/`--no-moderated-covariates`).

## [0.10.0] — 2026-08-13

### Added
- `plot_pathway_enrichment`: a dot plot (term vs. -log10 p-value, dot
  size and color tied to gene overlap) for both local and live pathway
  enrichment results. Generated automatically by `run_analysis` and the
  CLI when enrichment was run.

## [0.9.1] — 2026-08-13

### Changed
- Finished the OmicForge rename at the code level: package import path,
  CLI command, output directory name. The previous release only renamed
  user-facing text.
- README rewritten. Real, verified quick-start examples moved near the
  top, a benchmark results table added, an explicit Limitations section
  added, narrative tone reduced (about 30% shorter).

### Removed
- `requirements.txt`, unmaintained since v0.2 and duplicated by
  `pyproject.toml`.

## [0.9.0] — 2026-08-13

### Added
- Second independent benchmark dataset (`pasilla`, *Drosophila*)
  alongside the original `airway` (human) benchmark, to check whether
  the precision/recall pattern found on one dataset held on another. It
  did.
- Live pathway enrichment against the g:Profiler API
  (`run_gprofiler_enrichment`): GO, KEGG, Reactome, WikiPathways.
  Opt-in, requires outbound internet access.
- GEO series data acquisition (`deconcord.io.geo`): list and download a
  GEO series' supplementary files, and parse its sample metadata from a
  `GSE` accession.

## [0.8.0] — 2026-08-12

### Added
- Gene ID to symbol annotation (`load_gene_annotation`,
  `annotate_de_results`) from a local two-column CSV. Wired into
  `run_analysis` (`gene_annotation`) and the CLI (`--annotation`).

## [0.7.0] — 2026-08-12

### Added
- Empirical-Bayes moderated t-test (`compute_moderated_pvalues`),
  selectable via `method="moderated"`. Shrinks each gene's variance
  estimate toward a value fit from every other gene's variance, the
  same idea behind limma's `eBayes`. More statistical power than the
  default per-gene Welch's t-test at the same precision.

## [0.6.0] — 2026-08-12

### Added
- Command-line interface (`deconcord counts.csv --group1 ... --group2
  ...`) as a real console script, installed via `pip install -e .`.

## [0.5.0] — 2026-08-11

### Added
- `load_gmt`: a parser for the standard tab-separated GMT gene-set file
  format, so pathway enrichment can be run against a real MSigDB or
  Enrichr download instead of a hand-built Python dict.

## [0.4.0] — 2026-08-11

### Added
- Pre-DE low-count gene filtering (`filter_low_count_genes`, opt-in via
  `min_count`). Drops genes that never had a realistic chance of
  reaching significance, improving multiple-testing correction for
  every other gene.

### Changed
- `alpha` and `lfc_threshold` (previously hardcoded) are now real,
  validated parameters on `run_differential_expression` and
  `run_analysis`.

## [0.3.x] — 2026-08-09 to 2026-08-11

### Fixed
- A correctness pass across the whole pipeline (v0.3.1 to v0.3.3):
  several bugs found by re-reading the code closely while rewriting
  docstrings, including NaN handling in constant-expression genes and
  multiple-testing correction.

## [0.3.0] — 2026-08-09

### Added
- Pathway enrichment: hypergeometric test against an explicit
  background gene universe (`run_pathway_enrichment_analysis`).

## [0.2.0] — 2026-08-08

### Added
- Visualization: volcano plot and PCA plot.

## [0.1.0] — 2026-08-07

### Added
- Differential expression: log fold change, per-gene Welch's t-test,
  Benjamini-Hochberg correction (`run_differential_expression`).

## [0.0.1] to [0.0.3] — 2026-08-06 to 2026-08-07

### Added
- Count matrix validation (non-negative integers, unique gene IDs, no
  missing values).
- Sample QC (library size, genes detected, MAD-based outlier
  detection).
- Normalization (CPM, log2 transform).
