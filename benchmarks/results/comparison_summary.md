# DEConcord (welch) vs DESeq2 / edgeR — airway dataset

Real numbers from a real run. See `benchmarks/README.md` for methodology and honest caveats before reading anything into these on their own.

```
         comparison  genes_in_common  pearson_r_lfc  spearman_r_lfc  deconcord_significant  other_significant  significant_in_both  jaccard_index  precision_vs_other  recall_vs_other
DEConcord vs DESeq2            16139          0.938           0.982                    196               2618                  196          0.075                 1.0            0.075
 DEConcord vs edgeR            15890          0.947           0.983                    196               1990                  196          0.098                 1.0            0.098
```
