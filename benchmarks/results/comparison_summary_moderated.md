# BioInsight (moderated) vs DESeq2 / edgeR — airway dataset

Real numbers from a real run. See `benchmarks/README.md` for methodology and honest caveats before reading anything into these on their own.

```
          comparison  genes_in_common  pearson_r_lfc  spearman_r_lfc  bioinsight_significant  other_significant  significant_in_both  jaccard_index  precision_vs_other  recall_vs_other
BioInsight vs DESeq2            16139          0.938           0.982                     431               2618                  431          0.165                 1.0            0.165
 BioInsight vs edgeR            15890          0.947           0.983                     431               1990                  431          0.217                 1.0            0.217
```
