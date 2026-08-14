#!/usr/bin/env Rscript
#
# Runs DESeq2 and edgeR on the `airway` dataset (Himes et al. 2014, PLoS ONE;
# GEO GSE52778) using each package's own standard, idiomatic workflow --
# not tuned or nudged toward agreeing with DEConcord in any way. Also
# writes out the raw count matrix and sample metadata as CSVs so DEConcord
# can be run on the exact same input.
#
# `airway`: 8 human airway smooth muscle samples, 4 cell lines, each with an
# untreated and a dexamethasone-treated sample. It's the dataset DESeq2's
# own vignette uses as its worked example, which makes it a fair, standard
# choice here rather than something cherry-picked to make DEConcord look
# good.
#
# Usage:
#   Rscript benchmarks/run_deseq2_edger.R
#
# Requires: BiocManager, DESeq2, edgeR, airway (installed automatically
# below if missing -- comment that block out if you'd rather manage your
# own R library).

options(warn = 1)

required_bioc <- c("DESeq2", "edgeR", "airway")
missing_pkgs <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(airway)
})

data_dir <- file.path("benchmarks", "data")
results_dir <- file.path("benchmarks", "results")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

data(airway)
counts <- assay(airway)
meta <- as.data.frame(colData(airway))[, c("SampleName", "cell", "dex")]
meta$dex <- factor(meta$dex, levels = c("untrt", "trt"))

# --- write raw inputs so DEConcord can be run on the identical matrix ---
counts_out <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
write.csv(counts_out, file.path(data_dir, "airway_counts.csv"), row.names = FALSE)
write.csv(data.frame(sample = colnames(counts), dex = meta$dex, cell = meta$cell),
          file.path(data_dir, "airway_metadata.csv"), row.names = FALSE)

cat(sprintf(
  "airway: %d genes x %d samples. Groups -- untrt: %s | trt: %s\n",
  nrow(counts), ncol(counts),
  paste(colnames(counts)[meta$dex == "untrt"], collapse = ", "),
  paste(colnames(counts)[meta$dex == "trt"], collapse = ", ")
))

# =====================================================================
# DESeq2 -- standard workflow, default settings, contrast trt vs untrt
# =====================================================================
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~dex)
dds <- DESeq(dds)
deseq2_res <- results(dds, contrast = c("dex", "trt", "untrt"))
deseq2_df <- as.data.frame(deseq2_res)
deseq2_df$gene_id <- rownames(deseq2_df)
deseq2_df <- deseq2_df[, c("gene_id", "log2FoldChange", "pvalue", "padj")]
write.csv(deseq2_df, file.path(results_dir, "deseq2_results.csv"), row.names = FALSE)
cat(sprintf("DESeq2: %d genes tested, %d significant (padj < 0.05)\n",
            sum(!is.na(deseq2_df$padj)), sum(deseq2_df$padj < 0.05, na.rm = TRUE)))

# =====================================================================
# edgeR -- standard workflow: filterByExpr, TMM, quasi-likelihood F-test
# =====================================================================
y <- DGEList(counts = counts, group = meta$dex)
keep <- filterByExpr(y, group = meta$dex)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~dex, data = meta)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "dextrt")
edger_df <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none"))
edger_df$gene_id <- rownames(edger_df)
edger_df <- edger_df[, c("gene_id", "logFC", "PValue", "FDR")]
colnames(edger_df) <- c("gene_id", "logFC", "pvalue", "FDR")
write.csv(edger_df, file.path(results_dir, "edger_results.csv"), row.names = FALSE)
cat(sprintf("edgeR: %d genes tested (post filterByExpr), %d significant (FDR < 0.05)\n",
            nrow(edger_df), sum(edger_df$FDR < 0.05, na.rm = TRUE)))

cat("\nDone. Raw inputs in benchmarks/data/, results in benchmarks/results/.\n")
cat("Next: python benchmarks/compare_results.py\n")
