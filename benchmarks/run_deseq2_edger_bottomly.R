#!/usr/bin/env Rscript
#
# Fourth external benchmark dataset -- same purpose as run_deseq2_edger.R,
# run_deseq2_edger_pasilla.R, and run_deseq2_edger_zebrafish.R, picked
# specifically to give the research/ml/ rank-stability finding (see
# research/ml/FINDINGS.md) a real chance at an evaluable external test.
# zebrafish (3 vs 3) turned out non-evaluable for the frozen ML endpoint:
# edgeR found zero significant genes on that dataset at all, not because
# of anything about rank_stability, just thin replication against a real
# significance bar. Bottomly has roughly triple the replication.
#
# `countsBottomly` (dexus package): Bottomly et al. 2011, mouse striatum
# RNA-seq, two inbred strains -- C57BL/6J ("B6") vs DBA/2J ("D2"), 21
# samples total (10 B6, 11 D2). Bioconductor's own documentation traces
# this object back to the original ReCount Bottomly data; it's a
# real, published, widely-used DE-methods benchmark dataset (used in
# Soneson & Delorenzi 2013 and Rapaport et al. 2013, among others),
# picked here specifically because it has more replication than any
# dataset checked so far, not because it's a new organism or design.
#
# Group membership is derived from countsBottomly's own column names
# below, NOT hardcoded from a remembered sample-ID list -- exactly the
# thing that went wrong for zebrafish (Trt1/Trt3/Trt5 guessed instead of
# the real Trt9/Trt11/Trt13, without R access to check). The script prints
# the loaded column names, classifies each by a B6/D2 pattern match, and
# refuses to proceed (stop()) unless that classification lands on exactly
# 10 B6 and 11 D2 -- the verified design -- with a diagnostic dump of the
# actual column names either way, since that's the only thing that lets
# whoever runs this catch a wrong guess before it silently produces a
# wrong benchmark instead of after.
#
# Usage:
#   Rscript benchmarks/run_deseq2_edger_bottomly.R
#
# Requires: BiocManager, DESeq2, edgeR, dexus (installed automatically
# below if missing).

options(warn = 1)

required_bioc <- c("DESeq2", "edgeR", "dexus")
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
  library(dexus)
})

data_dir <- file.path("benchmarks", "data")
results_dir <- file.path("benchmarks", "results")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- load counts, print exactly what was loaded before doing anything with it ---
data(countsBottomly)
counts <- as.matrix(countsBottomly)
sample_names <- colnames(counts)

cat(sprintf("Loaded countsBottomly: %d genes x %d samples.\n", nrow(counts), ncol(counts)))
cat("Column names:\n")
print(sample_names)

# --- derive B6/D2 from the column names, validate strictly, refuse to guess ---
is_b6 <- grepl("b6|c57", sample_names, ignore.case = TRUE)
is_d2 <- grepl("d2|dba", sample_names, ignore.case = TRUE)

if (any(is_b6 & is_d2)) {
  stop(
    "Ambiguous strain classification: at least one column name matched both B6 ",
    "and D2 patterns. Column names: ", paste(sample_names, collapse = ", ")
  )
}
if (any(!is_b6 & !is_d2)) {
  unmatched <- sample_names[!is_b6 & !is_d2]
  stop(
    "Could not classify ", length(unmatched), " column name(s) as B6 or D2 from ",
    "a name-pattern match: ", paste(unmatched, collapse = ", "),
    ". Full column names: ", paste(sample_names, collapse = ", "),
    ". Inspect colnames(countsBottomly) (and any attached phenotype data, e.g. via ",
    "Biobase::pData if it turns out to be an ExpressionSet rather than a plain ",
    "matrix) and adjust the classification pattern above before rerunning."
  )
}

condition <- ifelse(is_b6, "B6", "D2")
n_b6 <- sum(condition == "B6")
n_d2 <- sum(condition == "D2")
cat(sprintf("Classified from column names: %d B6, %d D2 (verified design: 10 B6, 11 D2).\n", n_b6, n_d2))

if (n_b6 != 10 || n_d2 != 11) {
  stop(sprintf(
    paste(
      "Sample counts after classification (%d B6, %d D2) do not match the",
      "verified Bottomly design (10 B6, 11 D2). Refusing to proceed with a",
      "possibly-wrong grouping. Column names were: %s"
    ),
    n_b6, n_d2, paste(sample_names, collapse = ", ")
  ))
}

coldata <- data.frame(row.names = sample_names, condition = factor(condition, levels = c("B6", "D2")))

# --- write raw inputs so DEConcord can be run on the identical matrix ---
counts_out <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
write.csv(counts_out, file.path(data_dir, "bottomly_counts.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(coldata), condition = coldata$condition),
          file.path(data_dir, "bottomly_metadata.csv"), row.names = FALSE)

cat(sprintf(
  "bottomly: %d genes x %d samples. Groups -- B6: %s | D2: %s\n",
  nrow(counts), ncol(counts),
  paste(rownames(coldata)[coldata$condition == "B6"], collapse = ", "),
  paste(rownames(coldata)[coldata$condition == "D2"], collapse = ", ")
))

# =====================================================================
# DESeq2 -- standard workflow, default settings, contrast D2 vs B6.
# B6 (C57BL/6J) is the field's usual reference strain, so it's the
# denominator/reference level here, matching that convention rather than
# an arbitrary pick.
# =====================================================================
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
dds <- DESeq(dds)
deseq2_res <- results(dds, contrast = c("condition", "D2", "B6"))
deseq2_df <- as.data.frame(deseq2_res)
deseq2_df$gene_id <- rownames(deseq2_df)
deseq2_df <- deseq2_df[, c("gene_id", "log2FoldChange", "pvalue", "padj")]
write.csv(deseq2_df, file.path(results_dir, "bottomly_deseq2_results.csv"), row.names = FALSE)
cat(sprintf("DESeq2: %d genes tested, %d significant (padj < 0.05)\n",
            sum(!is.na(deseq2_df$padj)), sum(deseq2_df$padj < 0.05, na.rm = TRUE)))

# =====================================================================
# edgeR -- standard workflow: filterByExpr, TMM, quasi-likelihood F-test
# =====================================================================
y <- DGEList(counts = counts, group = coldata$condition)
keep <- filterByExpr(y, group = coldata$condition)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~condition, data = coldata)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "conditionD2")
edger_df <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none"))
edger_df$gene_id <- rownames(edger_df)
edger_df <- edger_df[, c("gene_id", "logFC", "PValue", "FDR")]
colnames(edger_df) <- c("gene_id", "logFC", "pvalue", "FDR")
write.csv(edger_df, file.path(results_dir, "bottomly_edger_results.csv"), row.names = FALSE)
cat(sprintf("edgeR: %d genes tested (post filterByExpr), %d significant (FDR < 0.05)\n",
            nrow(edger_df), sum(edger_df$FDR < 0.05, na.rm = TRUE)))

cat("\nDone. Raw inputs in benchmarks/data/, results in benchmarks/results/.\n")
cat("Next: python benchmarks/compare_results.py --dataset bottomly\n")
