# Model Findings: Does `rank_stability` Predict Cross-Method Agreement?

This is an exploratory research finding, not a production claim and not a
biological result. It lives under `research/ml/` on purpose, see the top of
`EXPERIMENTAL_REPORT.md` and `INVENTORY.md` for why this track stays
separate from the public `deconcord` package. Nothing here is merged into
`src/deconcord/`, marketed, or asserted as validated.

**Status: the airway/pasilla discovery phase for this hypothesis is
complete.** The locked evaluation below is the strongest and most rigorous
result produced, and the one to cite if you cite one. Everything after it
in this document is the supporting work that built up to it. Next priority
is external dataset expansion under the same frozen experiment design, not
a more complex model, see "Next" at the bottom.

## Headline result: locked L1-logistic evaluation

The most rigorous experiment run on this hypothesis, and the one this
finding should be summarized by. `train_l1_locked.py`, L1-penalized
logistic regression (`liblinear`, `class_weight="balanced"`), regularization
strength `C` chosen by 5-fold stratified cross-validation *on the training
dataset only* (`C` in `{0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0}`), then locked
and evaluated once on the held-out dataset. The held-out dataset is never
touched during model selection, this is the closest this track has come to
a real train/select/test split rather than a PR-AUC comparison at a fixed,
unselected hyperparameter.

Borderline DESeq2 genes only (adjusted p-value 0.01-0.10). Conventional
features: `abs_logfc`, `neg_log10_p`. Combined: same, plus
`source_rank_stability`.

**Locked `pasilla -> airway` evaluation** (trained on `pasilla`, `C`
selected by CV on `pasilla` alone, scored once on `airway`):

| Model | Selected C | PR-AUC | ROC-AUC |
|---|---:|---:|---:|
| Conventional | 3 | 0.5017 | 0.8520 |
| + `source_rank_stability` | 10 | 0.6504 | 0.8843 |

Delta PR-AUC: **+0.1487**. `source_rank_stability` coefficient in the
locked model: **+2.4305** (positive, and by far the largest-magnitude
coefficient in that model, on standardized features).

This is a stronger, more carefully validated version of the borderline
result reported further down this document, same direction, same
"`rank_stability` specifically, not the other robustness columns" pattern,
but with the model-selection step itself kept honest about what it was and
wasn't allowed to see.

## Setup

Task: predict `target_method_significant` (edgeR `FDR < 0.05`) for a gene
from DESeq2's own per-gene output on that gene, on `airway` and `pasilla`.
Model: logistic regression (`StandardScaler` + `LogisticRegression`,
`class_weight="balanced"`), chosen over a more expressive model for the same
reason stated when this track was scoped, understanding what the model
learned mattered more than maximizing a score.

Three feature groups, defined in `train_borderline.py` and reused across the
other scripts:

- **Conventional DE**: `abs_logfc` (`|source_log_fold_change|`),
  `neg_log10_p` (`-log10(source_p_value)`).
- **Robustness**: `source_direction_stability`, `source_rank_stability`
  (see `INVENTORY.md` Section C; `source_significance_stability` was tried
  and dropped, see below).
- **Combined**: both groups together.

All cross-dataset numbers below are trained on one dataset and evaluated on
the other (`airway -> pasilla` means fit on `airway`, score on `pasilla`),
never a within-dataset split, per the leakage concern flagged when this
track was designed: adjusted p-values are computed jointly across every
gene in a single DE run, so genes within one dataset aren't independent
rows, and a random within-dataset split would overstate how well this
generalizes.

## Why `source_significance_stability` was dropped

It's the weakest of the three robustness columns in practice. DEConcord's
own baseline significance calls (the Welch's-t-test workflow the robustness
features are computed on, see the methodological note in `INVENTORY.md`)
are sparse, especially on `airway`, which leaves too few
baseline-significant genes for a fraction-of-reruns-significant feature to
discriminate on. It's kept in the generated tables for completeness but
excluded from the feature sets below.

## All-gene results

Across the full gene set (no adjusted-p windowing), conventional DE features
alone already predict edgeR agreement well: PR-AUC around 0.97 in both
train/test directions. A robustness-only model (direction + rank stability,
no logFC or p-value at all) still performs strongly but is clearly weaker
than conventional features. Adding robustness to conventional features
moves PR-AUC only slightly. Reading this plainly: for the typical gene, the
DESeq2 p-value and effect size already say almost everything there is to
say about whether edgeR will agree. There isn't much room left for a
robustness feature to add.

## Earlier supporting result: unlocked borderline-window comparison

The experiment that led to the headline result above. Same borderline
regime, same feature groups, but an ordinary (unregularized) logistic
regression at default settings rather than the locked, CV-selected L1 model.
Reported here because it's what motivated everything that followed, not as
the number to cite on its own.

Genes DESeq2 itself is unsure about, adjusted p-value between 0.01 and 0.10
(`train_borderline.py`). Genes DESeq2 is confident about, deep in the tails,
need less help; genes near its own decision boundary are exactly where an
independent robustness signal could matter most, which is also the premise
this whole research track was built to test.

| Direction | Conventional PR-AUC | + `rank_stability` PR-AUC |
|---|---:|---:|
| airway -> pasilla | 0.6676 | 0.7437 |
| pasilla -> airway | 0.5019 | 0.6373 |

`direction_stability` added nothing in this regime; it's exactly 1.0 for
every gene in these borderline subsets (a borderline-significance gene
apparently still keeps a consistent sign across DEConcord's leave-one-out
reruns, even when its significance call itself wobbles), so it carries zero
variance to split on. The entire gain above comes from `rank_stability`.

## Residualization check

The obvious objection: maybe `rank_stability` is just repackaging
`abs_logfc`/`neg_log10_p` under a different name, in which case the gain
above would be spurious. `train_residual_rank.py` tests this directly: fit
`rank_stability ~ abs_logfc + neg_log10_p` (linear regression, fit on train,
applied to both train and test) and use the residual, the part of
`rank_stability` conventional DE features can't already explain, as the
extra feature instead of the raw value.

| Direction | Conventional PR-AUC | + residual `rank_stability` PR-AUC |
|---|---:|---:|
| airway -> pasilla | 0.6676 | 0.7420 |
| pasilla -> airway | 0.5019 | 0.6482 |

The gain survives almost unchanged. Whatever `rank_stability` is capturing
in the borderline regime, it isn't simply a restatement of effect size and
p-value.

## L1 regularization checks

Two more checks, both supporting the headline result rather than replacing
it. `train_l1_ablation.py` sweeps the same `C` grid used by the locked
evaluation on the combined feature set alone and confirms `rank_stability`
keeps a nonzero, non-trivial coefficient across the full range, an L1
penalty is specifically willing to zero out a feature it finds redundant,
and doesn't here. `train_l1_matched.py` compares conventional-vs-combined
at each fixed `C` (rather than separately tuning and locking one `C` per
feature set) and finds the combined model ahead of the conventional one at
matched regularization strength across the grid, consistent with the locked
result rather than an artifact of `C` selection favoring the combined model.

## Window sensitivity

Five overlapping adjusted-p windows were checked (`train_window_sensitivity.py`,
same conventional-vs-plus-rank comparison as above): 0.01-0.05, 0.02-0.08,
0.01-0.10, 0.025-0.075, 0.05-0.10. Adding `rank_stability` improved PR-AUC
in 9 of the 10 window x direction combinations. This wasn't one lucky
threshold choice.

## Bootstrap uncertainty

`bootstrap_rank_gain.py` resamples the test genes (1,000 bootstrap draws,
`random_state=42`) and reports the 95% interval of delta-PR-AUC
(rank-augmented minus conventional-only) for each of the same 10 window x
direction combinations:

- 7/10 settings: 95% interval clearly positive (excludes zero).
- 2/10 settings: interval crosses zero.
- 1/10 setting: clearly negative, `pasilla -> airway`, window 0.05-0.10.

That negative window is also the sparsest one tested, roughly a dozen
positive-class genes in the test set at that point, small enough that a
single bootstrap draw can swing the estimate substantially. It's reported
here rather than dropped, a result that doesn't survive its own uncertainty
check in every setting is still real information about where this
signal is (and isn't) reliable.

Important scope note on these intervals: they quantify sampling
uncertainty *over test genes within one dataset*, given the model already
fit. They are not, and shouldn't be read as, evidence about generalization
to a third dataset. With `n=2` datasets total, no bootstrap procedure can
substitute for that.

## Current interpretation

- Robustness features add essentially nothing when predicting agreement
  across the whole gene set, conventional DE statistics already carry
  almost all the signal there.
- `rank_stability` specifically, not `direction_stability`, not
  `significance_stability`, adds real signal in the regime where DESeq2's
  own call is already borderline. It survives residualization against
  conventional DE strength, survives L1 regularization across a range of
  penalty strengths, and produces the largest gain of any experiment run
  in the locked, CV-selected evaluation above.
- This is the shape of result the project's stated long-term aim expects:
  robustness information should matter most exactly when the ordinary
  DE conclusion is ambiguous, not when it's already obvious.
- This predicts **cross-method significance agreement** (does edgeR agree
  with a borderline DESeq2 call), **not biological replication**. Those are
  different questions; see Caveats.

## Caveats (not a checklist to relax later, an honest accounting now)

- **Two datasets.** `airway` and `pasilla` are the only data this has been
  checked on. Every cross-dataset number above is real, but two datasets
  is not a claim of generality.
- **Post-hoc window selection.** The 0.01-0.10 borderline window was chosen
  after seeing that the all-gene result was uninteresting, then five
  neighboring windows were checked to see if the finding was an artifact of
  that specific choice. It wasn't, per the sensitivity sweep, but the
  original window itself was picked by looking at the data first, worth
  saying plainly rather than presenting it as pre-registered.
- **Correlated genes.** Genes within a dataset are not independent
  (co-regulation, shared pathways), so effective sample size is smaller
  than raw gene count for both the model fit and the bootstrap CIs.
- **Bootstrap CIs are within-dataset uncertainty, not cross-dataset proof.**
  See above.
- **This predicts cross-method agreement, not biological replication.**
  DESeq2-agrees-with-edgeR is a real, useful, but different question from
  the project's actual long-term target, whether a discovery-dataset
  finding replicates in an independent dataset. Reaching that target needs
  paired discovery/replication datasets DEConcord doesn't currently have
  access to, this result is evidence the underlying idea (robustness
  signals carry independent information near a decision boundary) is worth
  pursuing toward that target, not a stand-in for it.
- **Robustness features are computed on DEConcord's own DE workflow, not
  DESeq2's.** See `INVENTORY.md`'s methodological note. `rank_stability` is
  a proxy signal, not literally "how stable is DESeq2's own call."
- **Gene-level CV/bootstrap is not study-level generalization.** The locked
  evaluation's cross-validation, and the bootstrap intervals further up,
  both quantify uncertainty over genes within the datasets already used.
  Neither substitutes for testing on an independent third dataset, which is
  the actual next step, see below.

## Next: external dataset expansion before more complex models

The airway/pasilla discovery phase for this hypothesis is done. The next
scientific priority is testing the same, now-frozen experiment design
(borderline window, conventional + `rank_stability` features, locked L1
logistic regression with train-only CV for `C`) against a third,
independent dataset, not building a more complex model on the same two
datasets. Two datasets can show a pattern is real and survives several
robustness checks, which this one has; they can't show it generalizes past
those two datasets, and a more expressive model fit to the same two
datasets would risk mistaking dataset-specific quirks for a real effect
before that question is even asked. Dataset selection and acquisition for
this step haven't started yet.

## Reproduce

From the repository root, with the research dependencies available:

```bash
.venv/bin/python research/ml/build_tables.py
.venv/bin/python research/ml/generate_robustness.py
.venv/bin/python research/ml/train_baseline.py
.venv/bin/python research/ml/train_borderline.py
.venv/bin/python research/ml/train_residual_rank.py
.venv/bin/python research/ml/train_window_sensitivity.py
.venv/bin/python research/ml/bootstrap_rank_gain.py
.venv/bin/python research/ml/train_l1_ablation.py
.venv/bin/python research/ml/train_l1_matched.py
.venv/bin/python research/ml/train_l1_locked.py
```

None of these scripts currently write their results to a file, the numbers
in this document were transcribed from their stdout output. A natural next
step, if this track continues, is having them write structured results
(JSON or CSV) the way `build_tables.py` already does with
`build_report.json`, so future runs can be diffed instead of re-transcribed.
