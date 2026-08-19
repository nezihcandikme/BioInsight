# DESeq2 vs edgeR concordance — airway dataset

Real numbers from a real run, computed by `deconcord.concordance.methods.compute_de_concordance` on the committed DE result tables. Both tools ran with their own default settings on the identical count matrix and design -- neither is being validated against the other here, they're being checked for whether they agree.

```
genes_compared: 15896
significant_DESeq2: 2596
significant_edgeR: 1990
significant_in_both: 1977
jaccard_index: 0.7578
pearson_r_lfc: 0.9995
spearman_r_lfc: 0.9992
directional_agreement: 1.0000
significant_union: 2609
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 1977
discordant_genes: 0
only_in_DESeq2: 619
only_in_edgeR: 13
```

# DESeq2 vs limma-voom concordance — airway dataset

Real numbers from a real run, computed by `deconcord.concordance.methods.compute_de_concordance` on the committed DE result tables. Both tools ran with their own default settings on the identical count matrix and design -- neither is being validated against the other here, they're being checked for whether they agree.

```
genes_compared: 15896
significant_DESeq2: 2596
significant_limma-voom: 1882
significant_in_both: 1871
jaccard_index: 0.7177
pearson_r_lfc: 0.9901
spearman_r_lfc: 0.9884
directional_agreement: 1.0000
significant_union: 2607
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 1871
discordant_genes: 0
only_in_DESeq2: 725
only_in_limma-voom: 11
```

# edgeR vs limma-voom concordance — airway dataset

Real numbers from a real run, computed by `deconcord.concordance.methods.compute_de_concordance` on the committed DE result tables. Both tools ran with their own default settings on the identical count matrix and design -- neither is being validated against the other here, they're being checked for whether they agree.

```
genes_compared: 15926
significant_edgeR: 1991
significant_limma-voom: 1882
significant_in_both: 1848
jaccard_index: 0.9126
pearson_r_lfc: 0.9897
spearman_r_lfc: 0.9889
directional_agreement: 1.0000
significant_union: 2025
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 1848
discordant_genes: 0
only_in_edgeR: 143
only_in_limma-voom: 34
```
