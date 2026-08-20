# DEConcord (welch) vs DESeq2 / edgeR — bottomly dataset

Real numbers from a real run. See `benchmarks/README.md` for methodology and honest caveats before reading anything into these on their own.

```
         comparison  genes_in_common  pearson_r_lfc  spearman_r_lfc  deconcord_significant  other_significant  significant_in_both  jaccard_index  precision_vs_other  recall_vs_other
DEConcord vs DESeq2            23828          0.780           0.956                    273               3718                  273          0.073                 1.0            0.073
 DEConcord vs edgeR            21662          0.859           0.975                    272               3264                  272          0.083                 1.0            0.083
```
