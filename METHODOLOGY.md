# Methodology

What DEConcord actually computes, why those specific metrics and defaults,
what they assume, and — most importantly — what they don't mean. This
file exists because a metric without a documented interpretation is worse
than no metric at all: it invites reading a number as stronger evidence
than it is.

If a metric described here is a standard, established statistic, that's
noted explicitly with a source. DEConcord has not introduced any novel
statistical method (see [Limitations](#limitations-of-this-methodology)
below and the README's own Roadmap) — everything currently implemented is
a well-understood quantity, applied to a comparison of DE result tables.

## Method concordance

Implemented in `deconcord.concordance.methods.compute_de_concordance`.
Takes two differential expression result tables for the *same* comparison
(same samples, same two groups) — produced by any two tools, methods, or
analytical choices — and reports where they agree and disagree.

### Inputs

Each table needs, at minimum, a log fold change column and an adjusted
p-value column, indexed by gene ID. Column names are configurable because
DESeq2 (`log2FoldChange`, `padj`), edgeR (`logFC`, `FDR`), and DEConcord's
own DE output (`log_fold_change`, `adjusted_p_value`) all use different
ones.

Genes are restricted to the intersection of both tables with a
non-missing log fold change and p-value in both (`genes_compared` in the
output). A gene present in only one table, or with a `NA` p-value in
either (e.g. DESeq2's independent filtering), is excluded from every
metric below rather than treated as "not significant" — an excluded gene
and a tested-but-not-significant gene are different things, and
conflating them would understate genuine method-specific missingness.

### Significance

A gene is called significant in a table if `adjusted_p_value < alpha`
(default `alpha = 0.05`). This is the standard, if arbitrary,
convention — there is nothing statistically special about 0.05, and a
different threshold is a legitimate, different analytical choice (this is
exactly what "threshold sensitivity," on the roadmap, is meant to
evaluate directly instead of picking one number and moving on).

### Metrics

**Significant-gene overlap and Jaccard index.** Given the significant-gene
sets from table A and table B, `significant_in_both = |A ∩ B|`, and the
Jaccard index is `|A ∩ B| / |A ∪ B|` — a standard set-similarity measure
(Jaccard, 1912; see e.g. [Real & Vargas, 1996](https://doi.org/10.1093/sysbio/45.3.380)
for a review of its use as a similarity coefficient). Ranges 0 (no
overlap) to 1 (identical significant-gene sets). `NaN` if neither table
calls anything significant — there is no meaningful ratio when the union
is empty.

**Log fold change correlation.** Pearson and Spearman correlation
(standard statistics — see any introductory statistics reference, e.g.
[Pearson, 1895](https://doi.org/10.1098/rspl.1895.0041) and
[Spearman, 1904](https://doi.org/10.2307/1412159)) of the log fold change
values across *every* gene compared, not just the significant ones. This
measures effect-size agreement independent of whichever significance
threshold was chosen — two tools can have near-identical fold change
estimates while still disagreeing about which specific genes clear a
p-value cutoff, and that distinction matters for interpretation.

**Directional agreement.** Among genes significant in *both* tables, the
fraction where `sign(log_fold_change_A) == sign(log_fold_change_B)`. This
is not a standard named statistic — it's a direct, literal count, kept
deliberately simple rather than wrapped in a more elaborate agreement
coefficient (e.g. Cohen's kappa), because the interesting cases (0 or a
small number of genes) are better read directly as "these N genes
disagree on direction" than compressed into a corrected-for-chance score
that needs its own explanation.

**Concordant / discordant / method-specific gene sets.** Explicit gene
lists, not just counts: `concordant_genes` (significant in both,
same-signed effect), `discordant_genes` (significant in both,
opposite-signed effect — the case actually worth inspecting gene-by-gene,
since it means the analyses disagree about the *direction* of a real
effect, not just its statistical significance), and `only_in_<name>`
(significant in one table only).

### Interpretation

**Disagreement is information, not evidence that one method is wrong.**
DESeq2 and edgeR use different statistical models (different dispersion
estimation, different default filtering); a gene one calls significant
and the other doesn't is frequently a marginal case sitting near both
tools' own significance boundary, not a mistake by either tool. High
concordance is reassuring; low concordance on a specific gene is a signal
to look at that gene's data more closely, not automatically a bug report
against either method.

**A high overall concordance score does not certify every individual
finding.** Two datasets with 0.75+ Jaccard overlap (see the README's
Validation section) still have method-specific genes on each side — the
aggregate number describes the comparison as a whole, not any one gene's
result.

**Concordance is not correctness.** Two methods can agree with each other
and both be wrong in the same way (e.g. both missing a batch effect
neither was told about). Concordance measures agreement between the
methods actually checked, not agreement with ground truth — there usually
isn't ground truth available for a real RNA-seq experiment.

## Threshold sensitivity

Implemented in `deconcord.concordance.threshold_sensitivity.compute_threshold_sensitivity`.
Sweeps `compute_de_concordance` above across a range of `alpha` values
between the same two DE result tables, to check whether the concordance
conclusion depends on which threshold in that range was picked.

### What it computes

For each `alpha` in a supplied list (default `(0.01, 0.05, 0.1)`), runs
`compute_de_concordance` at that threshold and collects the summary
metrics (`jaccard_index`, `significant_in_both`, `directional_agreement`,
etc.) into a table indexed by alpha (`by_alpha`). Log fold change
correlation isn't included there since it doesn't depend on alpha at
all — it's computed across every compared gene, not just the significant
ones, so it would just repeat the same number on every row.

Separately, for each gene, computes what fraction of the swept alphas
call it significant in each table (`gene_stability`). A fraction of 1.0
means the gene is significant at every threshold tried; 0.0 means it's
significant at none of them. Either extreme means the gene's own
significance call doesn't depend on which threshold in the swept range
was used — it's threshold-stable. Anything strictly between 0 and 1
means the call flips somewhere in the range: whether that gene counts as
a "finding" in that table depends on an arbitrary choice, not just on
the data.

A gene is reported as `stable` only if *both* fractions (in table A and
in table B) are 0.0 or 1.0. This is a statement about each table's own
sensitivity to threshold choice, not about whether the two tables agree
with each other — a gene that's significant at every alpha in table A
and non-significant at every alpha in table B is threshold-stable (its
status doesn't waver with the threshold) even though the two tables
disagree about it. Threshold stability and method concordance are
different questions; conflating them would hide genes whose disagreement
is a real, consistent disagreement rather than a borderline call.

### Interpretation

A gene significant across the whole swept range in both tables is a
stronger finding than one that's only significant at the loosest
threshold in one table — the latter is one threshold change away from
disappearing from the "significant" list entirely. A high `jaccard_index`
at every swept alpha means the two tables' agreement isn't an artifact of
having picked exactly 0.05; a `jaccard_index` that swings a lot across
the range means the agreement number itself shouldn't be quoted without
saying which threshold it's at.

The range swept is still a choice. `(0.01, 0.05, 0.1)` is a reasonable
default spread around the conventional cutoff, not a claim that
everything outside it is unreasonable — a caller checking a genuinely
different convention (e.g. a stricter Bonferroni-style threshold) should
pass their own `alphas`.

## Pathway stability

Implemented in `deconcord.concordance.pathway_stability.compute_pathway_stability`.
The same overlap-and-Jaccard idea as method concordance's significant-gene
comparison, applied to two pathway enrichment result tables instead of
two gene-level DE result tables.

### What it computes

Given two enrichment tables (e.g. one from a DE run using Welch's t-test,
another from the same data using the moderated t-test — or one from
DESeq2's significant genes, one from edgeR's), restricts to pathways
tested in both, calls a pathway significant in a table if its p-value is
below `alpha` (default 0.05), and reports the significant-pathway
overlap and Jaccard index, plus explicit `stable_pathways` (significant
in both) and `only_in_<name>` lists.

This is deliberately narrower than method concordance. There's no
directional-agreement or fold-change-correlation analogue here, because
pathway enrichment tables don't have a shared, comparable effect-size
column the way DE tables share log fold change — local enrichment
reports an overlap count against whatever background was supplied, live
g:Profiler reports an intersection size against its own background, and
the two aren't on the same scale. Comparing presence/absence of
significance is the part that's actually well-defined across both
sources.

### Interpretation

A pathway enriched under every DE configuration checked is a more
robust biological conclusion than one that only shows up under one
specific choice of DE method or threshold — the latter could be a real,
subtle effect that only one method has the power to detect, or it could
be noise that happened to clear one method's particular threshold.
Pathway stability doesn't distinguish between those two explanations by
itself; it only tells you which case you're in, so you know whether
further investigation is warranted.

Pathway-level (in)stability inherits every uncertainty already present
in the gene-level DE calls and the enrichment step's own background
choice — see [Limitations](#limitations-of-this-methodology) below.

## Resampling stability

Implemented in `deconcord.concordance.resampling_stability.compute_resampling_stability`.
Method concordance, threshold sensitivity, and pathway stability all ask
whether a conclusion survives a different reasonable *analytical* choice.
Resampling stability asks whether it survives a different *sample set*:
does a gene's significance call hold up if the analysis had happened to
include (or exclude) one particular sample?

### What it computes

Runs `run_differential_expression` once on the full sample set (the
baseline), then reruns it many times on perturbed sample sets, and
reports what fraction of reruns still call each gene significant. Two
resampling schemes are offered:

- `"subsample"`: repeatedly draws a random subset of each group (without
  replacement) and reruns DE. Randomized (seeded via `random_state`);
  needs enough iterations for the per-gene fractions to be meaningful.
- `"leave_one_out"`: drops each sample exactly once and reruns DE.
  Deterministic and exhaustive, but a gentler, more limited perturbation
  than subsampling, and the iteration count is fixed at the sample count.

Bootstrap resampling (drawing samples *with* replacement) is deliberately
not offered. `run_differential_expression` refuses a group with a
repeated sample name, on purpose, to catch a real class of caller bugs.
A bootstrap draw routinely produces repeats, and supporting it would mean
reworking how the DE functions index into the expression matrix to
tolerate duplicates — a bigger change than this function's scope
justifies. Subsampling and leave-one-out give the same "how much does
the exact sample set matter" signal without that rework.

A gene's significance call is `stable` if it reappears in at least
`stability_threshold` (default 0.9) of reruns when the baseline called it
significant, or in at most `1 - stability_threshold` of reruns when the
baseline didn't. Everything else is `sensitive`: baseline-significant
genes that don't reliably replicate, which is the case worth the closest
look.

### Interpretation

A gene significant in the baseline and in nearly every resampled rerun
is a broadly reproducible signal, not an artifact of which samples
happened to be collected. A gene significant in the baseline but only
sometimes in reruns is a much weaker finding — it may be real but
underpowered, or it may be driven by one or two samples. Resampling
stability doesn't tell you which; it tells you the call is fragile and
worth a closer look at the individual samples before treating it as
settled.

This is a per-gene diagnostic, not a whole-dataset one. A dataset can
have many stable genes and a handful of sensitive ones; the fragile genes
are the ones this function exists to surface, not a verdict on the
dataset as a whole.

## Limitations of this methodology

- Every metric here assumes both DE result tables were run on the *same*
  samples and the *same* two-group comparison. Comparing tables from
  different experiments isn't meaningful and isn't guarded against beyond
  requiring overlapping gene IDs — the caller is responsible for making
  sure the comparison itself makes sense.
- `compute_de_concordance` itself still evaluates significance at one
  threshold at a time; `concordant_genes`/`discordant_genes` from that
  function are specific to whatever `alpha` was passed. Threshold
  sensitivity (above) is a separate function you call in addition, not
  something `compute_de_concordance` does automatically.
- Small sample sizes reduce the reliability of everything downstream of
  them (the p-values and fold changes being compared, not the concordance
  computation itself) — a low-power comparison can look either falsely
  concordant (both methods underpowered the same way) or falsely
  discordant (both near their own significance boundary in opposite
  directions by chance). Concordance numbers should be read alongside the
  sample size and effect sizes involved, not in isolation.
- Pathway enrichment depends entirely on the gene-set database and
  background gene universe supplied — a pathway's enrichment status can
  change with either, independent of anything about the DE methods being
  compared. Pathway stability inherits this: two enrichment tables built
  on different backgrounds or gene-set collections aren't measuring the
  same thing, and `compute_pathway_stability` has no way to detect that
  from the tables alone.
- Technical confounding (batch effects, subtle covariates) cannot always
  be detected automatically. DEConcord's covariate-adjustment feature
  (see the README) can *correct for* a known, measured confound, but
  concordance between two methods that share the same *unmeasured*
  confound will not reveal it.
- `compute_resampling_stability`'s `"subsample"` mode and default
  `stability_threshold` of 0.9 are both choices, not laws — a smaller
  `subsample_fraction` is a harsher perturbation and will flag more genes
  as sensitive, and a lower `stability_threshold` will flag fewer. A
  dataset with very few samples per group has little room to subsample
  from in the first place, which limits how informative the result can
  be regardless of settings.
- None of this establishes biological causality. RNA-seq differential
  expression, however robust or concordant across methods, is
  observational — a stable, method-independent finding is a stronger
  empirical claim than a fragile one, not proof of a causal mechanism.
