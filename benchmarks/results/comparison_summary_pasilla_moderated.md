# DEConcord (moderated) vs DESeq2 / edgeR — pasilla dataset

Real numbers from a real run. See `benchmarks/README.md` for methodology and honest caveats before reading anything into these on their own.

```
         comparison  genes_in_common  pearson_r_lfc  spearman_r_lfc  deconcord_significant  other_significant  significant_in_both  jaccard_index  precision_vs_other  recall_vs_other
DEConcord vs DESeq2             7920          0.974           0.991                    146                837                  146          0.174                 1.0            0.174
 DEConcord vs edgeR             7872          0.975           0.988                    146                680                  146          0.215                 1.0            0.215
```
