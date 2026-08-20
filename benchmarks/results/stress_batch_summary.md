# Stress test: batch effects

Real numbers from a real run: DESeq2, edgeR, and limma-voom, each with its own default settings, compared pairwise with `deconcord.concordance.methods.compute_de_concordance` under two design formulas -- `naive` (~condition, ignoring the real per-sample library-size batch split) and `adjusted` (~depth_batch + condition, modeling it). See `benchmarks/run_stress_batch.R` for exactly how the real, data-derived batch split was constructed and why.

## Real batch x condition split

```
condition    B6  D2
depth_batch        
high_depth    4   6
low_depth     6   5
```

## Model: naive

**DESeq2 vs edgeR**

```
genes_compared: 21662
significant_DESeq2: 3432
significant_edgeR: 3264
significant_in_both: 3226
jaccard_index: 0.9297
pearson_r_lfc: 0.9998
spearman_r_lfc: 0.9999
directional_agreement: 1.0000
significant_union: 3470
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 3226
discordant_genes: 0
```

**DESeq2 vs limma-voom**

```
genes_compared: 21662
significant_DESeq2: 3432
significant_limma-voom: 3136
significant_in_both: 3075
jaccard_index: 0.8803
pearson_r_lfc: 0.9847
spearman_r_lfc: 0.9749
directional_agreement: 1.0000
significant_union: 3493
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 3075
discordant_genes: 0
```

**edgeR vs limma-voom**

```
genes_compared: 21662
significant_edgeR: 3264
significant_limma-voom: 3136
significant_in_both: 3069
jaccard_index: 0.9213
pearson_r_lfc: 0.9842
spearman_r_lfc: 0.9745
directional_agreement: 1.0000
significant_union: 3331
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 3069
discordant_genes: 0
```

## Model: adjusted

**DESeq2 vs edgeR**

```
genes_compared: 21662
significant_DESeq2: 3349
significant_edgeR: 3162
significant_in_both: 3112
jaccard_index: 0.9156
pearson_r_lfc: 0.9986
spearman_r_lfc: 0.9973
directional_agreement: 1.0000
significant_union: 3399
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 3112
discordant_genes: 0
```

**DESeq2 vs limma-voom**

```
genes_compared: 21662
significant_DESeq2: 3349
significant_limma-voom: 3015
significant_in_both: 2969
jaccard_index: 0.8745
pearson_r_lfc: 0.9848
spearman_r_lfc: 0.9740
directional_agreement: 1.0000
significant_union: 3395
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 2969
discordant_genes: 0
```

**edgeR vs limma-voom**

```
genes_compared: 21662
significant_edgeR: 3162
significant_limma-voom: 3015
significant_in_both: 2963
jaccard_index: 0.9219
pearson_r_lfc: 0.9863
spearman_r_lfc: 0.9776
directional_agreement: 1.0000
significant_union: 3214
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 2963
discordant_genes: 0
```

## naive vs adjusted, side by side

```
   model                 pair  genes_compared  jaccard_index  pearson_r_lfc  directional_agreement  discordant_genes
   naive      DESeq2 vs edgeR           21662       0.929683       0.999790                    1.0                 0
   naive DESeq2 vs limma-voom           21662       0.880332       0.984676                    1.0                 0
   naive  edgeR vs limma-voom           21662       0.921345       0.984166                    1.0                 0
adjusted      DESeq2 vs edgeR           21662       0.915563       0.998617                    1.0                 0
adjusted DESeq2 vs limma-voom           21662       0.874521       0.984814                    1.0                 0
adjusted  edgeR vs limma-voom           21662       0.921904       0.986308                    1.0                 0
```
