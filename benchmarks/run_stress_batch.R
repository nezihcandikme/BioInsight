#!/usr/bin/env Rscript
#
# Second "harder data regime" in the stress-test suite: batch effects.
# See run_stress_unbalanced.R's header for the same framing this shares
# -- this is DEConcord's core method-concordance question
# (method_concordance.py) applied outside the two clean, balanced
# two-group datasets checked so far.
#
# A batch effect needs a real confound to be a real test -- inventing one
# by scaling counts for an arbitrary sample subset would fabricate the
# exact signal being tested for, not stress-test anything. Instead, this
# script uses a real, data-derived technical grouping that's already
# present in the real Bottomly counts: per-sample library size (total
# read count), computed directly from the real count matrix, no
# fabrication anywhere. Samples are split into "low_depth"/"high_depth"
# batches by a median split on their real library size -- a genuine
# technical property of a real sequencing run, exactly the kind of thing
# that produces real batch effects in practice (different sequencing
# runs, different lanes, different library prep batches routinely differ
# in depth).
#
# The actual regime: does correctly modeling this real technical
# covariate change method concordance, relative to ignoring it? Two
# models per tool, same real counts, same real groups, only the design
# formula differs:
#   - naive:    ~condition               (batch structure ignored)
#   - adjusted: ~depth_batch + condition (batch structure modeled)
# If depth_batch happens to correlate with condition (B6 vs D2) in the
# real data, that's reported directly below, not hidden or corrected for
# -- a real confound between a technical and biological variable is
# itself part of what this regime is checking.
#
# Runs on Bottomly's real, full "balanced" design (10 B6, 11 D2) --
# same real counts benchmarks/run_deseq2_edger_bottomly.R and
# run_stress_unbalanced.R's "balanced" level already use, not a new
# acquisition.
#
# Usage:
#   Rscript benchmarks/run_stress_batch.R
#
# Requires: BiocManager, DESeq2, edgeR, limma (installed automatically
# below if missing). Requires benchmarks/data/bottomly_counts.csv and
# bottomly_metadata.csv to already exist -- run
# `Rscript benchmarks/run_deseq2_edger_bottomly.R` first if they don't.

options(warn = 1)

required_bioc <- c("DESeq2", "edgeR", "limma")
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
  library(limma)
})

data_dir <- file.path("benchmarks", "data")
results_dir <- file.path("benchmarks", "results", "stress_batch")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(data_dir, "bottomly_counts.csv")
metadata_path <- file.path(data_dir, "bottomly_metadata.csv")
if (!file.exists(counts_path) || !file.exists(metadata_path)) {
  stop(sprintf(
    paste(
      "Missing %s and/or %s. Run `Rscript benchmarks/run_deseq2_edger_bottomly.R`",
      "first -- this script builds the batch split from Bottomly's real counts",
      "and doesn't reacquire them itself."
    ),
    counts_path, metadata_path
  ))
}

counts <- read.csv(counts_path, row.names = "gene_id", check.names = FALSE)
counts <- as.matrix(counts)
meta <- read.csv(metadata_path)

cat(sprintf("Loaded bottomly_counts.csv: %d genes x %d samples.\n", nrow(counts), ncol(counts)))
cat(sprintf("Loaded bottomly_metadata.csv: %d samples.\n", nrow(meta)))

if (sum(meta$condition == "B6") != 10 || sum(meta$condition == "D2") != 11) {
  stop(sprintf(
    "Expected 10 B6 / 11 D2 (the verified Bottomly design), found %d B6 / %d D2. Refusing to proceed.",
    sum(meta$condition == "B6"), sum(meta$condition == "D2")
  ))
}

# Order metadata to match the count matrix's own column order -- the
# real, data-derived thing this script does is compute library size per
# column, so the two need to line up before anything else happens.
meta <- meta[match(colnames(counts), meta$sample), ]
if (any(is.na(meta$sample)) || !identical(meta$sample, colnames(counts))) {
  stop("bottomly_metadata.csv's samples don't match bottomly_counts.csv's columns after reordering -- inspect both files.")
}

# --- real, data-derived batch: median split on real per-sample library size ---
library_size <- colSums(counts)
cat("Per-sample library size (real, computed from the real count matrix):\n")
print(data.frame(sample = names(library_size), condition = meta$condition, library_size = library_size))

median_size <- median(library_size)
depth_batch <- ifelse(library_size > median_size, "high_depth", "low_depth")
# Ties (a sample landing exactly on the median) go to low_depth by this
# rule -- deterministic, not adjustable after seeing the split, and
# printed below either way so it's checkable.
n_high <- sum(depth_batch == "high_depth")
n_low <- sum(depth_batch == "low_depth")
cat(sprintf("\nMedian library size: %.0f. Split: %d high_depth, %d low_depth.\n", median_size, n_high, n_low))

# Report, don't correct for, any real correlation between the technical
# batch and the biological condition -- this is information about the
# real data, not something to rebalance.
cat("\nBatch x condition contingency table (real data, reported as-is):\n")
print(table(depth_batch, meta$condition))

coldata <- data.frame(
  row.names = colnames(counts),
  condition = factor(meta$condition, levels = c("B6", "D2")),
  depth_batch = factor(depth_batch, levels = c("low_depth", "high_depth"))
)

write.csv(
  data.frame(sample = rownames(coldata), condition = coldata$condition, depth_batch = coldata$depth_batch,
             library_size = library_size),
  file.path(results_dir, "batch_metadata.csv"),
  row.names = FALSE
)

# =====================================================================
# Per-tool, per-model run. "naive" ignores the real batch structure
# (~condition only); "adjusted" models it (~depth_batch + condition).
# Same real counts and groups in both -- only the design formula differs.
# =====================================================================
run_deseq2 <- function(counts, coldata, design, out_path) {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = design)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "D2", "B6"))
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df <- df[, c("gene_id", "log2FoldChange", "pvalue", "padj")]
  write.csv(df, out_path, row.names = FALSE)
  cat(sprintf("  DESeq2: %d genes tested, %d significant (padj < 0.05)\n",
              sum(!is.na(df$padj)), sum(df$padj < 0.05, na.rm = TRUE)))
}

run_edger <- function(counts, coldata, design, out_path) {
  y <- DGEList(counts = counts, group = coldata$condition)
  keep <- filterByExpr(y, group = coldata$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "conditionD2")
  df <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none"))
  df$gene_id <- rownames(df)
  df <- df[, c("gene_id", "logFC", "PValue", "FDR")]
  colnames(df) <- c("gene_id", "logFC", "pvalue", "FDR")
  write.csv(df, out_path, row.names = FALSE)
  cat(sprintf("  edgeR: %d genes tested (post filterByExpr), %d significant (FDR < 0.05)\n",
              nrow(df), sum(df$FDR < 0.05, na.rm = TRUE)))
}

run_limma <- function(counts, coldata, design, out_path) {
  y <- DGEList(counts = counts, group = coldata$condition)
  keep <- filterByExpr(y, group = coldata$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  v <- voom(y, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  df <- topTable(fit, coef = "conditionD2", n = Inf, sort.by = "none")
  df$gene_id <- rownames(df)
  df <- df[, c("gene_id", "logFC", "P.Value", "adj.P.Val")]
  colnames(df) <- c("gene_id", "logFC", "pvalue", "adj.P.Val")
  write.csv(df, out_path, row.names = FALSE)
  cat(sprintf("  limma-voom: %d genes tested (post filterByExpr), %d significant (adj.P.Val < 0.05)\n",
              nrow(df), sum(df$adj.P.Val < 0.05, na.rm = TRUE)))
}

models <- list(
  naive = ~condition,
  adjusted = ~depth_batch + condition
)

for (model_name in names(models)) {
  design_formula <- models[[model_name]]
  cat(sprintf("\n===== Model: %s (%s) =====\n", model_name, deparse(design_formula)))

  run_deseq2(counts, coldata, design_formula, file.path(results_dir, sprintf("%s_deseq2_results.csv", model_name)))

  design_matrix <- model.matrix(design_formula, data = coldata)
  run_edger(counts, coldata, design_matrix, file.path(results_dir, sprintf("%s_edger_results.csv", model_name)))
  run_limma(counts, coldata, design_matrix, file.path(results_dir, sprintf("%s_limma_voom_results.csv", model_name)))
}

cat("\nDone. batch_metadata.csv and per-model DESeq2/edgeR/limma-voom results in benchmarks/results/stress_batch/.\n")
cat("Next: python benchmarks/analyze_stress_batch.py\n")
