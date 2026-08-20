# DEConcord (welch) vs DESeq2 / edgeR — zebrafish dataset

Real numbers from a real run. See `benchmarks/README.md` for methodology and honest caveats before reading anything into these on their own.

```
         comparison  genes_in_common  pearson_r_lfc  spearman_r_lfc  deconcord_significant  other_significant  significant_in_both  jaccard_index  precision_vs_other  recall_vs_other
DEConcord vs DESeq2            16837          0.908           0.923                      0                 57                    0            0.0                 NaN              0.0
 DEConcord vs edgeR            16833          0.909           0.924                      0                  0                    0            NaN                 NaN              NaN
```
