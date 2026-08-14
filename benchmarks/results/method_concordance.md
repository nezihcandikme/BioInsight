# DESeq2 vs edgeR concordance — airway dataset

Real numbers from a real run, computed by `deconcord.concordance.methods.compute_de_concordance` on the committed DESeq2/edgeR result tables. Both tools ran with their own default settings on the identical count matrix and design -- neither is being validated against the other here, they're being checked for whether they agree.

```
genes_compared: 15896
significant_DESeq2: 2596
significant_edgeR: 1990
significant_in_both: 1977
jaccard_index: 0.7578
pearson_r_lfc: 0.9995
spearman_r_lfc: 0.9992
directional_agreement: 1.0000
concordant_genes: 1977
discordant_genes: 0
only_in_DESeq2: 619
only_in_edgeR: 13
```
