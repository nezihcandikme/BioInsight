# Per-gene Data Inventory

This inventory describes committed artifacts inspected in `benchmarks/results/` and the current DeConcord output schemas. It makes no biological or predictive claim.

## A. Identifiers / metadata

- Input tables: `gene_id`; dataset and method identity are encoded by filenames.
- Generated ML tables: `gene_id`, `dataset_id`, `source_method`, `target_method`.

## B. Conventional DE features

- DESeq2: `log2FoldChange`, `pvalue`, `padj`.
- edgeR: `logFC`, `pvalue`, `FDR`.
- limma-voom: `logFC`, `pvalue`, `adj.P.Val`.
- The initial table exposes only renamed DESeq2 values: `source_log_fold_change`, `source_p_value`, and `source_adjusted_p_value`.

## C. DeConcord robustness features

- Available per gene from current inputs: `source_threshold_significance_fraction`, calculated over adjusted-p thresholds 0.01, 0.05, and 0.1.
- Available per gene, as of v0.18.0's `compute_resampling_stability`, run via `generate_robustness.py` with `resample_method="leave_one_out"` and joined into the ML tables by `gene_id`: `source_significance_stability` (from `frac_significant`), `source_direction_stability`, and `source_rank_stability`. See `src/deconcord/concordance/resampling_stability.py` for exact definitions (sign-match rate and one-minus-mean-percentile-deviation, both ranked/signed on log fold change, not p-value).
- Important methodological detail: these robustness features are computed on DEConcord's own from-scratch Welch's-t-test DE workflow (CPM to log2 transform, `compute_resampling_stability` with default `de_method="welch"`), run on the same raw counts and sample groups, not on DESeq2 itself. They're joined onto the DESeq2-vs-edgeR table by `gene_id` as a proxy signal, DEConcord's own DE robustness standing in for a robustness measurement DESeq2 doesn't expose. `FINDINGS.md` treats this as a real caveat on interpretation, not a data error.
- In practice, `source_significance_stability` has turned out to carry little usable signal in the current setup, DEConcord's own baseline significance calls are sparse (especially on `airway`), which limits how much a fraction-of-reruns-significant feature can discriminate. `source_direction_stability` is informative across the full gene set but collapses to exactly 1.0 within the borderline DESeq2 subsets used in `train_borderline.py` and related scripts, so it adds no discrimination there. `source_rank_stability` is the feature carrying the model gains reported in `FINDINGS.md`.
- Computable but excluded as target leakage: edgeR threshold fraction; pairwise significant-in-both, source-only/target-only, concordant/discordant status, and direction agreement derived from DESeq2 plus edgeR.
- Dataset-level only: Jaccard index, LFC Pearson/Spearman correlations, and aggregate directional agreement.

Bootstrap resampling is not implemented in the package (see `resampling_stability.py`'s module docstring for why), so no bootstrap-based robustness feature exists or is planned under the current package API. Cross-dataset replication labels and pathway-stability features mapped to genes remain absent, same as before.

## D. Candidate target variables

- Initial target: `target_method_significant`, defined as edgeR `FDR < 0.05`.
- Future method pairs may use the configured method-specific adjusted-p-value column.

## E. Leakage / target-derived features

- Every edgeR value is target-derived, including `logFC`, `pvalue`, `FDR`, ranks, signs, and threshold fractions.
- Every DESeq2/edgeR concordance, overlap, or direction feature uses the target method.
- Significance booleans and threshold fractions directly encode significance decisions; the source-only threshold fraction is retained but flagged.
- Adjusted p-values and features derived from them depend on the full set of genes tested within a dataset. Their rows are not statistically independent gene-level measurements.
- Random gene-level splits cannot demonstrate transfer to unseen datasets; validation must preserve `dataset_id` groups.

## Missing Robustness Features

Superseded as of v0.18.0 for significance/direction/rank stability, see Section C above. Still absent: any per-gene resampling stability computed on DESeq2 or edgeR's own model (only DEConcord's own DE workflow is resampled, see the methodological note in Section C), bootstrap stability values (not implemented in the package, see `resampling_stability.py`), cross-dataset replication labels, pathway-stability features mapped to genes, and target-independent pairwise concordance features. These must continue to be reported unavailable, not synthesized.

