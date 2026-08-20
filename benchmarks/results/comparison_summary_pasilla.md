# DEConcord (welch) vs DESeq2 / edgeR — pasilla dataset

Real numbers from a real run. See `benchmarks/README.md` for methodology and honest caveats before reading anything into these on their own.

```
         comparison  genes_in_common  pearson_r_lfc  spearman_r_lfc  deconcord_significant  other_significant  significant_in_both  jaccard_index  precision_vs_other  recall_vs_other
DEConcord vs DESeq2             7920          0.974           0.991                     71                837                   71          0.085                 1.0            0.085
 DEConcord vs edgeR             7872          0.975           0.988                     71                680                   71          0.104                 1.0            0.104
```
