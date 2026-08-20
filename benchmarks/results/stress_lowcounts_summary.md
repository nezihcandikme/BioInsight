# Stress test: low counts

Real numbers from a real run: DESeq2, edgeR, and limma-voom, each with its own default settings, compared pairwise with `deconcord.concordance.methods.compute_de_concordance` at three real sequencing-depth levels of the Bottomly dataset (`full` = real, unmodified depth; `moderate`/`low` = the same real reads binomially thinned to 25%/5%). See `benchmarks/run_stress_lowcounts.R` for exactly how thinning was done and why it doesn't fabricate any signal.

## Depth: full (749,860,180 total reads)

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

## Depth: moderate (187,473,037 total reads)

**DESeq2 vs edgeR**

```
genes_compared: 17585
significant_DESeq2: 2372
significant_edgeR: 2206
significant_in_both: 2179
jaccard_index: 0.9083
pearson_r_lfc: 0.9994
spearman_r_lfc: 0.9999
directional_agreement: 1.0000
significant_union: 2399
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 2179
discordant_genes: 0
```

**DESeq2 vs limma-voom**

```
genes_compared: 17585
significant_DESeq2: 2372
significant_limma-voom: 2134
significant_in_both: 2096
jaccard_index: 0.8697
pearson_r_lfc: 0.9910
spearman_r_lfc: 0.9874
directional_agreement: 1.0000
significant_union: 2410
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 2096
discordant_genes: 0
```

**edgeR vs limma-voom**

```
genes_compared: 17585
significant_edgeR: 2206
significant_limma-voom: 2134
significant_in_both: 2084
jaccard_index: 0.9238
pearson_r_lfc: 0.9895
spearman_r_lfc: 0.9874
directional_agreement: 1.0000
significant_union: 2256
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 2084
discordant_genes: 0
```

## Depth: low (37,493,306 total reads)

**DESeq2 vs edgeR**

```
genes_compared: 12933
significant_DESeq2: 852
significant_edgeR: 832
significant_in_both: 793
jaccard_index: 0.8900
pearson_r_lfc: 0.9996
spearman_r_lfc: 1.0000
directional_agreement: 1.0000
significant_union: 891
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 793
discordant_genes: 0
```

**DESeq2 vs limma-voom**

```
genes_compared: 12933
significant_DESeq2: 852
significant_limma-voom: 824
significant_in_both: 774
jaccard_index: 0.8581
pearson_r_lfc: 0.9910
spearman_r_lfc: 0.9900
directional_agreement: 1.0000
significant_union: 902
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 774
discordant_genes: 0
```

**edgeR vs limma-voom**

```
genes_compared: 12933
significant_edgeR: 832
significant_limma-voom: 824
significant_in_both: 788
jaccard_index: 0.9078
pearson_r_lfc: 0.9900
spearman_r_lfc: 0.9900
directional_agreement: 1.0000
significant_union: 868
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 788
discordant_genes: 0
```

## Gradient across depth levels

```
   depth  total_reads                 pair  genes_compared  jaccard_index  pearson_r_lfc  directional_agreement  discordant_genes
    full    749860180      DESeq2 vs edgeR           21662       0.929683       0.999790                    1.0                 0
    full    749860180 DESeq2 vs limma-voom           21662       0.880332       0.984676                    1.0                 0
    full    749860180  edgeR vs limma-voom           21662       0.921345       0.984166                    1.0                 0
moderate    187473037      DESeq2 vs edgeR           17585       0.908295       0.999394                    1.0                 0
moderate    187473037 DESeq2 vs limma-voom           17585       0.869710       0.990952                    1.0                 0
moderate    187473037  edgeR vs limma-voom           17585       0.923759       0.989524                    1.0                 0
     low     37493306      DESeq2 vs edgeR           12933       0.890011       0.999593                    1.0                 0
     low     37493306 DESeq2 vs limma-voom           12933       0.858093       0.991027                    1.0                 0
     low     37493306  edgeR vs limma-voom           12933       0.907834       0.989952                    1.0                 0
```
