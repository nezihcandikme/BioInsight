# DESeq2 vs edgeR concordance — pasilla dataset

Real numbers from a real run, computed by `deconcord.concordance.methods.compute_de_concordance` on the committed DESeq2/edgeR result tables. Both tools ran with their own default settings on the identical count matrix and design -- neither is being validated against the other here, they're being checked for whether they agree.

```
genes_compared: 7919
significant_DESeq2: 840
significant_edgeR: 684
significant_in_both: 670
jaccard_index: 0.7845
pearson_r_lfc: 0.9988
spearman_r_lfc: 0.9982
directional_agreement: 1.0000
concordant_genes: 670
discordant_genes: 0
only_in_DESeq2: 170
only_in_edgeR: 14
```
