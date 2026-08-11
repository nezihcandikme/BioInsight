# BioInsight — Devlog

The README tells you what BioInsight is. This tells you why it looks the way it does — the actual decisions, the bugs that forced them, and the things I considered and didn't do. Git commit messages say *what* changed. This is for the *why*, in more than one line.

Organized by day. Early entries are shorter because they're reconstructed from the commits and the code itself — I didn't start keeping a real devlog until the Aug 11 session, so I'm not going to pretend I remember internal debates from Aug 6 that I didn't write down.

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

## Aug 11 — session 1: make it harder to fool (v0.3.1)

A full pass through the whole pipeline looking for the gap between "the code doesn't crash" and "the code is right." Found more than expected:

- **Every `__init__.py` in the project except one was misspelled** (`_init_.py`, `__init_.py`, even `__init.py`) and the top-level `bioinsight/__init__.py` didn't exist at all. Python's namespace-package fallback was quietly papering over this, which is exactly the kind of thing that works until it doesn't. Fixed all of them, added the missing one.
- **Pathway enrichment wasn't restricting to the tested background before computing overlaps.** A gene that was never even eligible to be "significant" (outside the background universe) could still count toward the overlap and inflate the enrichment p-value. Fixed by intersecting `significant_genes` and each `gene_set` with `background_genes` before touching the hypergeometric math.
- **`compute_cpm` divided by zero silently** for any sample with zero total counts, producing NaN/inf instead of an error. Now it raises, and names the sample.
- **Constant-expression genes produced NaN p-values.** Welch's t-test is mathematically undefined when both groups have zero variance (e.g. an all-zero gene) — division by zero in the t-statistic. Rather than let NaN leak into the results table, special-cased it directly: identical constants means "not different" (p=1), differing constants means "as different as it gets" (p=0).
- Wrote `bioinsight.pipeline.run_analysis()` so the five separate modules could be run as one call instead of manually chaining them every time.
- Test count: 27 → 35.

Also fixed a duplicate test function name in `test_pathway_analysis.py` — Python silently let the second definition shadow the first, so half the intended coverage was never actually running. That kind of bug is exactly why "the tests pass" isn't proof of anything on its own.

## Aug 11 — session 2: humanize the docs, harden the remaining edges (v0.3.2)

Two things happened together here, and they ended up reinforcing each other: rewriting the docs to actually sound like a person wrote them, and rewriting the docs *forced* a re-read of every function closely enough to find more bugs.

**Bugs found while writing better docstrings:**
- `compute_adjusted_pvalues` would silently turn *every* adjusted p-value in the whole batch to NaN if even one input p-value was NaN (a `statsmodels` quirk, not obvious from the outside). Now it checks and raises, naming the offending gene.
- Same function crashed outright (`ZeroDivisionError`) on an empty p-value Series — found this while writing a test for "what if `run_pathway_enrichment_analysis` gets an empty `gene_sets` dict," which is a real edge case, not a hypothetical one.
- `plot_pca` handed you scikit-learn's fairly cryptic error if you gave it fewer than 2 samples or genes. Added a guard with an error that says what's actually wrong.
- `plot_volcano` could plot `-log10(0) = -inf` for a gene with an exact-zero adjusted p-value — which became newly *possible* once the constant-expression fix above started returning literal 0.0. Floored it before plotting.
- `validate_counts` now warns (doesn't fail — it's a heuristic, not a proof) if a matrix has way more columns than rows and very few rows, since that shape usually means someone loaded a transposed matrix.

**A deliberate breaking change:** pathway enrichment's result columns were `'Pathway'`, `'P-value'`, `'Adjusted P-value'`, `'Overlap Count'` — Title Case with spaces, while the differential expression results use snake_case (`p_value`, `adjusted_p_value`). Renamed the pathway columns to match. This breaks anything already depending on the old names, but there was nothing outside this repo depending on them yet, and shipping two naming conventions in the same codebase forever felt worse than fixing it once while it's still cheap.

**Group validation:** DE functions now reject missing sample names, a sample listed twice in the same group, and a sample appearing in both groups — each with a message naming the actual sample, instead of a bare pandas `KeyError` or a silently nonsensical comparison (a sample compared against itself).

Test count: 35 → 56 (constant-expression edge cases, empty gene sets, empty backgrounds, zero-library-size samples, PCA with too few samples, volcano with a zero p-value, the NaN-in-correction guard, and more).

## Aug 11 — session 3: configurable significance thresholds (v0.3.3)

`run_differential_expression` had `adjusted_p_value < 0.05` and `abs(log_fold_change) > 1` hardcoded into the significance call. Neither number is statistically special — they're conventions — so baking them into the function body meant anyone who disagreed with the convention had to edit the source. Added `alpha` and `lfc_threshold` as real parameters (same defaults, so nothing changes if you don't touch them), validated on the way in (`alpha` has to be in `(0, 1]`, `lfc_threshold` can't be negative), and threaded through `run_analysis` too.

Test count: 56 → 60.

## Aug 11 — session 4: pre-DE low-count filtering (v0.4.0)

Added the second item off the roadmap: a `filtering` module with `filter_low_count_genes(df, min_count=10, min_samples=2)` — drop genes that never hit `min_count` raw reads in at least `min_samples` samples. Reasoning, not just "real pipelines do this": every gene handed to `compute_pvalues` gets tested, and every test tightens the Benjamini-Hochberg correction for every *other* gene. A gene that was one read away from zero in every sample was never going anywhere except into that correction as dead weight.

Wired it into `run_analysis` as `min_count` (default `None`) and `min_samples` (default 2). QC still runs on the full, unfiltered matrix first — filtering before QC would let a dropped gene quietly change a sample's library size or genes-detected count without it showing up as what it is.

**Default-off was deliberate, not an oversight.** The function itself always filters when you call it directly. But wiring it into `run_analysis` *on* by default would have silently changed which genes get a p-value for every existing caller, without them asking for it — and it would have broken every existing pipeline test's assumption that the output has the same genes as the input. A statistically-better default that surprises the people already using the tool isn't obviously better. Opt-in for now; worth revisiting once there's a real default threshold worth defending instead of an arbitrary 10/2.

Test count: 60 → 68 (the filtering module directly: boundary conditions, invalid thresholds, min_samples bigger than the matrix, a filter that would remove every gene; plus one pipeline-level test that filtering actually removes the gene it's supposed to and leaves everything else alone).

---

## Decisions I looked at and didn't make

Worth writing down what got *rejected*, not just what shipped:

- **Silently supporting both old and new pathway column names**, instead of a clean breaking rename. Would've avoided the breaking change, but permanently — every future reader would have to know two names meant the same thing for no active reason. Rejected.
- **Validating matrix orientation as a hard error** instead of a warning. There's no way to be *certain* a wide, short matrix is transposed rather than just a small custom gene panel — it's a heuristic, and heuristics that hard-fail on a guess are worse than ones that speak up and let you decide.
- **Assigning an arbitrary small p-value (like 1e-10) instead of exactly 0.0** for constant-expression genes with a clear between-group difference. Went with exactly 0.0 because it's honest about what the test is actually saying ("as significant as this method can express"), and handled the downstream consequence (the volcano plot's `-log10(0)`) directly instead of avoiding it by fudging the number.
- **Reimplementing edgeR's `filterByExpr`** (which accounts for group sizes and normalized CPM, not just a flat raw-count cutoff) instead of a simple fixed threshold. `filterByExpr` is the more defensible method, but it's also a meaningfully bigger piece of statistical machinery to get right and test. Shipped the simple version now; the module docstring says outright that it isn't `filterByExpr`, so nobody mistakes "good enough to stop testing dead genes" for "the field-standard algorithm."
