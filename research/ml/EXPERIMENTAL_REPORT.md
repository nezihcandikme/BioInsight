# Experimental Feature-table Build Report

These research outputs are experimental data-engineering artifacts. They are not validated biological replication predictions.

## Feature table schema

Each row is one aligned gene for one dataset/source/target combination. Identity columns are `gene_id`, `dataset_id`, `source_method`, and `target_method`. Candidate predictor columns are `source_log_fold_change`, `source_p_value`, `source_adjusted_p_value`, `source_threshold_significance_fraction`, and, as of v0.18.0, three per-gene robustness columns joined in from `compute_resampling_stability`: `source_significance_stability`, `source_direction_stability`, `source_rank_stability`. The sole target is `target_method_significant` (edgeR `FDR < 0.05`). No edgeR measurement is exposed as a predictor. See `INVENTORY.md` Section C for what the robustness columns are computed on (DEConcord's own DE workflow, not DESeq2 itself) and for which of the three actually carries model-usable signal.

## Build results

| Dataset | Usable | Positive | Negative | Prevalence | Source-only | Target-only | Missing-required drops |
|---|---:|---:|---:|---:|---:|---:|---:|
| airway | 15,896 | 1,990 | 13,906 | 12.519% | 47,751 | 0 | 30 |
| pasilla | 7,919 | 684 | 7,235 | 8.637% | 6,680 | 0 | 0 |

All generated feature columns have zero missing values after alignment. Source-only genes are absent from the edgeR result table (edgeR's benchmark workflow applies expression filtering), so no target can be constructed for them. The 30 airway overlap drops have missing DESeq2 source values; the machine-readable report gives field-level counts. No pasilla overlap genes were dropped.

Four target-independent per-gene robustness features are now derivable: the DESeq2 threshold-significance fraction over adjusted-p thresholds 0.01, 0.05, and 0.1, plus the three `compute_resampling_stability` columns added in v0.18.0. All four depend on dataset-wide statistics (adjusted p-values or gene-relative rank within a resampling run), not just the single gene in isolation. See `INVENTORY.md` Section C for which of the four actually contributes model signal, and `FINDINGS.md` for the model results.

## Reproduce

From the repository root, run:

```bash
.venv/bin/python research/ml/build_tables.py
.venv/bin/pytest research/ml/test_feature_table.py -q
```

