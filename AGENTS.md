# Repository Guidelines

## Project Structure & Module Organization

DEConcord is a Python 3.9+ package using a `src` layout. Production code lives in `src/deconcord/`; major areas include `concordance/`, `differential_expression/`, `pathway_analysis/`, `io/`, and `visualization/`. The CLI entry point is `src/deconcord/cli.py`, while `pipeline.py` coordinates the full workflow. Tests mirror features in `tests/test_*.py`, with small sample inputs under `tests/fixtures/`. Use `examples/quickstart.py` for an end-to-end example. Benchmark scripts and generated reports live in `benchmarks/`; static documentation is in `docs/`.

## Build, Test, and Development Commands

Create an isolated environment, then install the package and development tools:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

- `pytest`: run the complete test suite.
- `pytest tests/test_concordance.py`: run one focused test module.
- `ruff check src/ tests/ benchmarks/ examples/`: check syntax, correctness, and Python 3.9-compatible annotations.
- `python examples/quickstart.py`: exercise the library locally without network or R dependencies.
- `deconcord --help`: inspect the installed CLI commands and options.

Install `.[ai]` only when working on the optional AI explanation feature.

## Coding Style & Naming Conventions

Use four-space indentation and keep lines within Ruff's 120-character limit. Follow standard Python naming: `snake_case` for functions, variables, and modules; `PascalCase` for classes; `UPPER_CASE` for constants. Add type hints and clear docstrings to public functions. Prefer early validation and actionable errors that identify the offending sample, column, or value. Run Ruff before submitting; its configured rules are `E9`, `F`, and `FA`.

## Testing Guidelines

Tests use pytest. Name files `test_<feature>.py` and test functions `test_<behavior>`. Every bug fix should include a regression test that fails before the fix. Keep fixtures small and deterministic, placing reusable data in `tests/fixtures/`. All tests must pass; no explicit coverage threshold is configured.

## Research / ML Work

Keep experimental ML work under `research/ml/`; do not import it from or expose it through the public `deconcord` package unless explicitly approved later. Exploratory work must not change production APIs, README claims, package versions, or public documentation unless requested. Use only robustness metrics DEConcord genuinely computes, report unavailable features clearly, and never invent features because they sound useful. Preserve `dataset_id`, `gene_id`, and method identity in ML tables for dataset-aware validation. Random gene-level train/test splits are not evidence of cross-dataset generalization. Flag features derived from the prediction target or target method as potential leakage before use. Establish interpretable baselines before adding expressive models, and do not tune toward a desired positive result: null results and failure to improve prediction are valid. Label generated outputs as experimental, not validated biological replication predictions.

## Commit & Pull Request Guidelines

Recent commits use concise, imperative subjects, often with a release tag, for example: `Add resampling stability (v0.15.0)`. Keep each commit and PR focused on one concern. PRs should explain what changed and why, link relevant issues, list tests run, and include screenshots for visual output changes. Update `CHANGELOG.md` under `[Unreleased]` for user-facing changes. For real releases, keep `CHANGELOG.md` and `DEVLOG.md` aligned, and synchronize version bumps in `pyproject.toml` and `CITATION.cff`. Preserve the project's narrow focus on robustness and concordance analysis for RNA-seq differential expression; open an issue before proposing broader scope.
