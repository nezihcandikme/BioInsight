# BioInsight — Devlog

The README tells you what BioInsight is. This tells you why it looks the way it does — the actual decisions, the bugs that forced them, and the things I considered and didn't do. Git commit messages say *what* changed. This is for the *why*, in more than one line.

Organized by day, not by session — a day can have several working sessions in it, but the entry below is one consolidated record of what changed and why, not a session-by-session transcript. Early entries are shorter because they're reconstructed from the commits and the code itself — I didn't start keeping a real devlog until Aug 11, so I'm not going to pretend I remember internal debates from Aug 6 that I didn't write down.

---

## Aug 6 — the CSV era

Started with the single least glamorous problem in the whole pipeline: is this CSV even a valid count matrix. Built the exception hierarchy (`CountMatrixError` and specific subclasses for non-integer, negative, duplicate-index, and missing-value cases) instead of one generic error, because "your data is bad" is useless and "column X has negative values" is actionable.

Decision worth naming: `validate_counts` collects *every* problem it finds and raises them together, instead of stopping at the first one. Nobody wants to fix one validation error, rerun, hit the next one, rerun again — better to get the whole list at once.

## Aug 7 — the data has entered the building

Added `load_count_matrix` (read CSV, validate, done), then the first QC metrics: library size and genes detected per sample. Outlier detection for those used MAD (median absolute deviation) instead of standard deviation — the whole point is catching a sample that's wildly different from the rest, and a couple of extreme samples will drag a mean and standard deviation right along with them, hiding themselves in the process. MAD doesn't have that problem.

CPM normalization also went in this day. This is the first place "sample A got sequenced deeper than sample B" stopped being a caveat and became something the code actually had to correct for.

## Aug 8 — statistics have consequences

log2(x + 1) transformation, then the first differential expression code: a mean-difference "log fold change" and a per-gene Welch's t-test. Welch's over the standard Student's t-test because there's no good reason to assume both groups have equal variance, and RNA-seq expression variance is routinely condition-dependent — assuming equal variance would've been the wrong default, not just a simplifying one.

This is roughly where "the code runs" stopped being the bar and "does this number mean what I think it means" took over. It would keep not being enough for a while.

## Aug 9 — thousands of p-values appeared

Benjamini-Hochberg correction, a minimum-sample-size check for the t-test, volcano and PCA plots, and hypergeometric pathway enrichment. This is the day BioInsight stopped being "the CSV thing" and became an actual small pipeline with an end product.

---

## Aug 11 — make it harder to fool, then make it honest about what it is

The longest day so far, by a wide margin. Grouped by topic below, not by sitting, since the topics are what matter a year from now — not how many times I got up for coffee.

**Correctness pass (v0.3.1).** A full walk through the whole pipeline looking for the gap between "the code doesn't crash" and "the code is right." Found more than expected:

- Every `__init__.py` in the project except one was misspelled (`_init_.py`, `__init_.py`, even `__init.py`), and the top-level `bioinsight/__init__.py` didn't exist at all. Python's namespace-package fallback was quietly papering over this — exactly the kind of thing that works until it doesn't. Fixed all of them.
- Pathway enrichment wasn't restricting to the tested background before computing overlaps, so a gene that was never even eligible to be "significant" could still inflate the enrichment p-value. Fixed by intersecting `significant_genes` and each `gene_set` with `background_genes` before touching the hypergeometric math.
- `compute_cpm` divided by zero silently for any sample with zero total counts. Now it raises, and names the sample.
- Constant-expression genes produced NaN p-values (Welch's t-test is undefined when both groups have zero variance). Special-cased it: identical constants means "not different" (p=1), differing constants means "as different as it gets" (p=0).
- A duplicate test function name in `test_pathway_analysis.py` had Python silently letting the second definition shadow the first, so half the intended coverage was never actually running.
- Wrote `bioinsight.pipeline.run_analysis()` so the modules could be run as one call instead of manually chained every time.

**Humanizing the docs surfaced more bugs (v0.3.2).** Rewriting docstrings to actually sound like someone who understands the function wrote them forced a close enough re-read to find real problems, not just prose ones:

- `compute_adjusted_pvalues` silently turned *every* adjusted p-value in the whole batch to NaN if even one input p-value was NaN (a `statsmodels` quirk). Now it checks and raises, naming the offending gene.
- Same function crashed outright (`ZeroDivisionError`) on an empty p-value Series — found while writing a test for "what if `run_pathway_enrichment_analysis` gets an empty `gene_sets` dict," a real edge case, not a hypothetical one.
- `plot_pca` handed you scikit-learn's fairly cryptic error below 2 samples or genes. Added a guard that says what's actually wrong.
- `plot_volcano` could plot `-log10(0) = -inf` for a gene with an exact-zero adjusted p-value — newly possible once the constant-expression fix above started returning literal 0.0. Floored it before plotting.
- `validate_counts` now warns (doesn't fail — it's a heuristic) if a matrix looks like it might be transposed.
- Pathway enrichment's result columns were `'Pathway'`, `'P-value'`, etc. (Title Case), while DE results use snake_case. Renamed the pathway columns to match — a deliberate breaking change, made while nothing outside this repo depended on the old names yet.
- DE functions now reject missing sample names, a sample listed twice in one group, and a sample appearing in both groups — each error names the actual sample, instead of a bare pandas `KeyError` or a silently nonsensical comparison.

**Configurable significance thresholds (v0.3.3).** `adjusted_p_value < 0.05` and `abs(log_fold_change) > 1` were hardcoded into `run_differential_expression`. Neither number is statistically special — they're conventions — so baking them in meant anyone who disagreed had to edit the source. Added `alpha` and `lfc_threshold` as real parameters (same defaults, validated on the way in), threaded through `run_analysis` too.

**Pre-DE low-count filtering (v0.4.0).** Added `bioinsight.filtering.filter_low_count_genes(df, min_count=10, min_samples=2)` — every gene handed to `compute_pvalues` gets tested, and every test tightens the Benjamini-Hochberg correction for every *other* gene, so a gene that was one read away from zero everywhere was only ever going to be dead weight. Wired into `run_analysis` as opt-in (`min_count=None` by default) — the function itself always filters when called directly, but defaulting the *pipeline* to filter would have silently changed which genes get a p-value for every existing caller. QC still runs on the raw, unfiltered matrix first, so filtering can't skew a sample's own quality numbers.

**Splitting the docs: README as the pitch, this file as the record.** The README had grown into something that tried to be a landing page, an install guide, an API reference, and a statistics disclaimer at once — readable, but not *explanatory* the way a first-time visitor needs. Stripped it down to what BioInsight is, why it exists, and where it's headed; moved every code block, module table, and exact threshold out. This file exists so none of that reasoning actually gets lost — it just moved.

**The website (`docs/`, GitHub Pages).** Built a static site — plain HTML/CSS and a little vanilla JS, deliberately, because GitHub Pages needs zero build step and the site shouldn't need a toolchain the rest of the project doesn't have. Design direction was Apple × scientific software × developer tool: black/white/gray, one accent color, system fonts, no gradients or icon clutter. The volcano and PCA images on the site are genuine output — generated by calling BioInsight's own `plot_volcano`/`plot_pca` on a synthetic 4,000-gene dataset, with only the color palette adjusted for the marketing image (the library's actual plotting code wasn't touched). Every claim and number on the page traces back to the README or this file — nothing on the site exists that isn't true of the code.

**GMT gene-set loading (v0.5.0).** Pathway enrichment took `gene_sets` as a plain Python `dict[str, set[str]]`, which meant translating an actual MSigDB or Enrichr download into Python by hand before you could use it. Added `bioinsight.pathway_analysis.gmt.load_gmt(path)` — a thin parser for the standard tab-separated GMT format (name, description, then genes), returning exactly the dict shape `run_pathway_enrichment_analysis` already expects. Kept it in its own file instead of `methods.py`, since parsing a file format and running a statistical test are different jobs that don't need to share a module just because they're both "pathway analysis." Rejects malformed lines, duplicate gene set names, and empty files by name and line number instead of producing a dict that's quietly missing or merging gene sets.

Test count across the day: 27 → 75.

---

## Aug 12 — a CLI, and cleaning up after myself

**Command-line entry point (v0.6.0).** Every single use of BioInsight so far meant writing a Python script that imports `run_analysis` and calls it — fine for me, annoying for literally anyone else who just wants to point the pipeline at a CSV and get results back. Added `bioinsight/cli.py`: an argparse wrapper around `run_analysis`, wired up as a real console script (`[project.scripts]` in `pyproject.toml`), so `pip install`-ing the package gives you a `bioinsight` command, not just an importable package.

It's deliberately thin — no logic of its own beyond argument parsing and writing outputs to disk (DE table and QC table as CSV, plots as PNG, pathway enrichment as CSV if a `--gmt` file was given, an explanation as a text file if `--explain` was passed). All the actual pipeline logic still lives in `run_analysis`; the CLI's only real design decision was what to do about a background gene universe for enrichment when the user doesn't hand-specify one with `--background` — it defaults to every gene in the input matrix, which is the same "everything that could have been tested" reasoning the library docstring already commits to elsewhere, just applied automatically instead of requiring the caller to know to do it.

Errors a user is likely to actually hit — a missing file, an unknown sample name, a malformed GMT line — get caught and printed as a one-line `bioinsight: ...` message with exit code 1, instead of a Python traceback. Anything outside that expected set of exceptions is left to crash loudly, on purpose: a CLI that swallows every exception into a generic "something went wrong" is worse than one that occasionally shows you a traceback for a real bug.

Tested it two ways: a normal pytest suite (`tests/test_cli.py`, 7 tests — happy path, `--no-plots`, `--gmt`, `--min-count` filtering, and three failure modes) calling `main()` directly, and then actually installing the package into a clean venv and running the real `bioinsight` command against a real CSV, because a CLI's argument parser and entry-point wiring are exactly the kind of thing that can pass every unit test and still be broken the moment a real shell invokes it.

**Also found and removed two committed accidents.** `test_pca.png` and `test_volcano.png` had been sitting at the repo root since the very first plotting commit (Aug 9, `v0.2`) — leftover output from manually eyeballing a plot locally, never meant to be tracked. Removed them and added `/test_*.png` and `/bioinsight_output/` (the CLI's default output folder) to `.gitignore` so it doesn't happen again, either by hand or by running the new CLI inside the repo.

**Site refresh.** Version bump aside, regenerated the volcano and PCA demo images on `docs/` by actually running the new `bioinsight` CLI against a fresh synthetic dataset, in-process, and copying its real output PNGs straight onto the site — more honest than the previous approach of calling the plotting functions directly, since it's now the exact same code path a real user's terminal would hit. The old demo happened to land on 207 significant genes; this run landed on 515 out of 4,000 with a different random seed's worth of synthetic effect sizes — both numbers are real, neither is picked to look good, and the hero terminal panel and capabilities grid (now 10 items, laid out 5x2) got updated to match: a CLI card, a GMT-loading card, and the hero code block now leads with the actual `bioinsight ...` command and its real stdout instead of a Python-only example.

Test count: 75 → 82.

**Simplifying the site.** Went back through `docs/` and cut the one piece of the design that didn't earn its place: a scroll-triggered fade-up animation on every section (an `IntersectionObserver` in `script.js` plus matching CSS). It looked fine, but it's exactly the kind of thing that reads as "trying too hard" on a page whose whole pitch is restraint — content should just be there when you land on it, not perform an entrance. Removed it outright rather than toning it down. Kept the animations that are actually functional (button press feedback, the mobile nav's slide-open) since those communicate something; cut the one that was purely decorative. Also swept up some small things while in there: two pieces of dead CSS nothing referenced, two inline `style=""` attributes moved into real classes, a softer drop shadow on the terminal panel, and a one-line Python-version note next to the install command, since that's a real detail worth having now that the README doesn't carry install specifics anymore.

**Started the DESeq2/edgeR benchmark harness.** The README has said since the docs split that the near-term plan is "check the current results against DESeq2 and edgeR on a real public dataset" — everything else in this project has been building toward being able to say that honestly, so it was time to actually do it. Sandbox reality check first: this environment has no root, no `apt`/`conda`, and network access restricted to `pypi.org` and `github.com` — R, CRAN, and Bioconductor are all unreachable from here. Rather than fake it with a lighter substitute, running the actual R-based comparison locally instead, so built the harness for that: `benchmarks/run_deseq2_edger.R` pulls the `airway` dataset (the same one DESeq2's own vignette uses — 8 human airway smooth muscle samples, 4 cell lines, untreated vs. dexamethasone-treated), runs DESeq2 and edgeR with each package's own standard, untuned workflow, and writes out both their results and the raw input matrix. `benchmarks/compare_results.py` then runs BioInsight on the identical matrix and groups and reports log fold change correlation and significant-gene overlap against each tool, plus scatter plots.

One sign-convention bug worth naming so it doesn't bite again: `compute_log_fold_change` returns `mean(group_1) - mean(group_2)`, so getting `group_1`/`group_2` backwards relative to DESeq2's `contrast=c("dex","trt","untrt")` would silently flip every correlation negative — not a code bug, just an easy way to misread real results as a disagreement that isn't there. Documented directly in the comparison script's own comments instead of leaving it as a trap.

Smoke-tested `compare_results.py` end to end against fabricated mock CSVs shaped like real DESeq2/edgeR output (never committed, deleted immediately after) to catch parsing and column-name bugs before handing this off — the real run, with real R installed, is next. Whatever those numbers turn out to be, favorable or not, they get written up here as-run.

**The real run.** Results are in `benchmarks/results/` (`comparison_summary.md`, `lfc_vs_deseq2.png`, `lfc_vs_edger.png`). The actual numbers, on the real `airway` dataset, real DESeq2 and real edgeR, both run with their own default settings:

- **Fold change direction and magnitude track well**: Pearson r = 0.938 against DESeq2, 0.947 against edgeR; Spearman ρ = 0.982 / 0.983. The two methods are looking at the same underlying signal and mostly agree on how big it is.
- **Precision is 1.0 against both.** Every single gene BioInsight calls significant (196 of them, out of 16,139 tested after `min_count=10` filtering), DESeq2 and edgeR also call significant. Zero disagreement in that direction — nothing BioInsight flags turns out to be a false alarm relative to either reference tool.
- **Recall is low: 0.075 against DESeq2 (2,618 significant), 0.098 against edgeR (1,990 significant).** BioInsight finds under a tenth of what the count-based tools find. This is the real, quantified cost of testing each gene in isolation with a t-test instead of a negative-binomial GLM that borrows statistical power across genes with similar expression levels (via empirical Bayes dispersion shrinkage) — exactly the mechanism `benchmarks/README.md` predicted would cause this gap, before the run happened.
- **One more thing the scatter plots show that wasn't predicted going in**: for the genes with the largest observed differences (DESeq2/edgeR log2FC of 5-10), BioInsight's estimate systematically undershoots theirs — points bend below the y=x line on the positive side, above it on the negative side. Most likely explanation: `results()` without `lfcShrink()` returns DESeq2's unshrunken MLE, which can produce very large coefficients for genes with near-zero counts in one group (this is precisely why DESeq2's own docs recommend `lfcShrink()` for reporting effect sizes rather than raw `results()` output). BioInsight's `log2(CPM + 1)` transform bounds this by construction — a pseudocount limits how large a ratio involving a near-zero count can get. Not claiming BioInsight is "more correct" here, just noting where and plausibly why the two disagree most.

So: the thing the README said BioInsight hadn't done yet — check it against the field's actual tools instead of describing itself as "not validated" and leaving it there — is done, on one real dataset. The honest one-line summary is precision 1.0, recall ~0.08-0.10 against DESeq2/edgeR defaults: when BioInsight finds something, it's very likely real; it just doesn't find most of what's there. That's a specific, falsifiable claim now, not a vague disclaimer.

---

## Decisions I looked at and didn't make

Worth writing down what got *rejected*, not just what shipped:

- **Silently supporting both old and new pathway column names**, instead of a clean breaking rename. Would've avoided the breaking change, but permanently — every future reader would have to know two names meant the same thing for no active reason. Rejected.
- **Validating matrix orientation as a hard error** instead of a warning. There's no way to be *certain* a wide, short matrix is transposed rather than just a small custom gene panel — it's a heuristic, and heuristics that hard-fail on a guess are worse than ones that speak up and let you decide.
- **Assigning an arbitrary small p-value (like 1e-10) instead of exactly 0.0** for constant-expression genes with a clear between-group difference. Went with exactly 0.0 because it's honest about what the test is actually saying ("as significant as this method can express"), and handled the downstream consequence (the volcano plot's `-log10(0)`) directly instead of avoiding it by fudging the number.
- **Reimplementing edgeR's `filterByExpr`** (which accounts for group sizes and normalized CPM, not just a flat raw-count cutoff) instead of a simple fixed threshold. `filterByExpr` is the more defensible method, but it's also a meaningfully bigger piece of statistical machinery to get right and test. Shipped the simple version now; the module docstring says outright that it isn't `filterByExpr`, so nobody mistakes "good enough to stop testing dead genes" for "the field-standard algorithm."
- **Duplicating install/usage commands in both the README and this file**, for convenience. Rejected — two copies of the same command drift the moment one changes and nobody updates the other. The README has the one copy that matters for actually running the thing.
- **A separate repo for the website** (e.g. a dedicated site repo, or the classic `username.github.io` pattern). Rejected in favor of `docs/` inside this repo: one source of truth, and a design/content update to the site can't quietly fall out of sync with the code it's describing.
- **`click` or `typer` instead of `argparse`** for the CLI. Both are nicer to write. Neither is worth a new dependency for a handful of flags that the standard library already handles fine — `argparse` stays until the CLI's surface area actually outgrows it.
