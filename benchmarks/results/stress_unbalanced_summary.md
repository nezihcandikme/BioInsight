# Stress test: unbalanced groups

Real numbers from a real run: DESeq2, edgeR, and limma-voom, each with its own default settings, compared pairwise with `deconcord.concordance.methods.compute_de_concordance` at three real sample-subset levels of the Bottomly dataset. See `benchmarks/run_stress_unbalanced.R` for exactly which samples make up each level and why. Not tuned toward any particular outcome; this is the plain output of three established tools disagreeing or agreeing with each other as one group's sample size shrinks.

## Level: balanced (10 B6 vs 11 D2)

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

## Level: moderate (5 B6 vs 11 D2)

**DESeq2 vs edgeR**

```
genes_compared: 22650
significant_DESeq2: 2401
significant_edgeR: 2086
significant_in_both: 2053
jaccard_index: 0.8435
pearson_r_lfc: 0.9992
spearman_r_lfc: 0.9998
directional_agreement: 1.0000
significant_union: 2434
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 2053
discordant_genes: 0
```

**DESeq2 vs limma-voom**

```
genes_compared: 22650
significant_DESeq2: 2401
significant_limma-voom: 1930
significant_in_both: 1898
jaccard_index: 0.7801
pearson_r_lfc: 0.9714
spearman_r_lfc: 0.9620
directional_agreement: 1.0000
significant_union: 2433
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 1898
discordant_genes: 0
```

**edgeR vs limma-voom**

```
genes_compared: 22652
significant_edgeR: 2087
significant_limma-voom: 1931
significant_in_both: 1878
jaccard_index: 0.8776
pearson_r_lfc: 0.9677
spearman_r_lfc: 0.9616
directional_agreement: 1.0000
significant_union: 2140
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 1878
discordant_genes: 0
```

## Level: severe (2 B6 vs 11 D2)

**DESeq2 vs edgeR**

```
genes_compared: 24126
significant_DESeq2: 1797
significant_edgeR: 1441
significant_in_both: 1338
jaccard_index: 0.7042
pearson_r_lfc: 0.9981
spearman_r_lfc: 0.9997
directional_agreement: 1.0000
significant_union: 1900
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 1338
discordant_genes: 0
```

**DESeq2 vs limma-voom**

```
genes_compared: 24126
significant_DESeq2: 1797
significant_limma-voom: 1319
significant_in_both: 1171
jaccard_index: 0.6021
pearson_r_lfc: 0.9593
spearman_r_lfc: 0.9668
directional_agreement: 1.0000
significant_union: 1945
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 1171
discordant_genes: 0
```

**edgeR vs limma-voom**

```
genes_compared: 24133
significant_edgeR: 1445
significant_limma-voom: 1323
significant_in_both: 1249
jaccard_index: 0.8223
pearson_r_lfc: 0.9513
spearman_r_lfc: 0.9671
directional_agreement: 1.0000
significant_union: 1519
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 1249
discordant_genes: 0
```

## Gradient across imbalance levels

```
   level  n_b6  n_d2                 pair  genes_compared  jaccard_index  pearson_r_lfc  directional_agreement  discordant_genes
balanced    10    11      DESeq2 vs edgeR           21662       0.929683       0.999790                    1.0                 0
balanced    10    11 DESeq2 vs limma-voom           21662       0.880332       0.984676                    1.0                 0
balanced    10    11  edgeR vs limma-voom           21662       0.921345       0.984166                    1.0                 0
moderate     5    11      DESeq2 vs edgeR           22650       0.843468       0.999186                    1.0                 0
moderate     5    11 DESeq2 vs limma-voom           22650       0.780107       0.971402                    1.0                 0
moderate     5    11  edgeR vs limma-voom           22652       0.877570       0.967709                    1.0                 0
  severe     2    11      DESeq2 vs edgeR           24126       0.704211       0.998146                    1.0                 0
  severe     2    11 DESeq2 vs limma-voom           24126       0.602057       0.959287                    1.0                 0
  severe     2    11  edgeR vs limma-voom           24133       0.822251       0.951299                    1.0                 0
```
