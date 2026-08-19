# DESeq2 vs edgeR concordance — pasilla dataset

Real numbers from a real run, computed by `deconcord.concordance.methods.compute_de_concordance` on the committed DE result tables. Both tools ran with their own default settings on the identical count matrix and design -- neither is being validated against the other here, they're being checked for whether they agree.

```
genes_compared: 7919
significant_DESeq2: 840
significant_edgeR: 684
significant_in_both: 670
jaccard_index: 0.7845
pearson_r_lfc: 0.9988
spearman_r_lfc: 0.9982
directional_agreement: 1.0000
significant_union: 854
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 670
discordant_genes: 0
only_in_DESeq2: 170
only_in_edgeR: 14
```

# DESeq2 vs limma-voom concordance — pasilla dataset

Real numbers from a real run, computed by `deconcord.concordance.methods.compute_de_concordance` on the committed DE result tables. Both tools ran with their own default settings on the identical count matrix and design -- neither is being validated against the other here, they're being checked for whether they agree.

```
genes_compared: 7919
significant_DESeq2: 840
significant_limma-voom: 670
significant_in_both: 655
jaccard_index: 0.7661
pearson_r_lfc: 0.9876
spearman_r_lfc: 0.9883
directional_agreement: 1.0000
significant_union: 855
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 655
discordant_genes: 0
only_in_DESeq2: 185
only_in_limma-voom: 15
```

# edgeR vs limma-voom concordance — pasilla dataset

Real numbers from a real run, computed by `deconcord.concordance.methods.compute_de_concordance` on the committed DE result tables. Both tools ran with their own default settings on the identical count matrix and design -- neither is being validated against the other here, they're being checked for whether they agree.

```
genes_compared: 7919
significant_edgeR: 684
significant_limma-voom: 670
significant_in_both: 660
jaccard_index: 0.9510
pearson_r_lfc: 0.9891
spearman_r_lfc: 0.9892
directional_agreement: 1.0000
significant_union: 694
directional_agreement_any_significant: 1.0000
opposite_direction_any_significant: 0
concordant_genes: 660
discordant_genes: 0
only_in_edgeR: 24
only_in_limma-voom: 10
```
