#!/usr/bin/env Rscript
#
# First "harder data regime" in the stress-test suite the README's roadmap
# named: unbalanced groups. The question isn't "does DEConcord's own DE
# match DESeq2/edgeR/limma-voom" (that's compare_results.py, already
# answered on airway/pasilla) -- it's "how much do independently developed
# tools still agree with *each other* once group sizes stop being clean
# and balanced," which is DEConcord's actual core question
# (method_concordance.py) applied outside the two clean, balanced
# two-group datasets checked so far.
#
# Data: real Bottomly et al. 2011 counts (benchmarks/data/bottomly_counts.csv
# / bottomly_metadata.csv, produced by run_deseq2_edger_bottomly.R -- run
# that first). Bottomly was picked as the base dataset specifically
# because, with 21 real samples (10 B6, 11 D2), it's the only dataset in
# hand with enough replicates to construct a genuine imbalance sweep
# without going below what any of the three tools can estimate a
# within-group variance from at all.
#
# No expression values are fabricated anywhere in this script. Every
# "regime" below is a real subset of Bottomly's real samples -- what
# changes is which and how many of them go into each group, not the
# counts themselves.
#
# Three imbalance levels, D2 held at its full 11 samples throughout and
# B6 shrunk, not both sides shrunk together -- shrinking both at once
# would confound "the groups are unbalanced" with "there's less total
# data," two different stresses. Holding one side fixed isolates the
# imbalance ratio as the actual variable being tested:
#   - balanced:  10 B6 vs 11 D2 (the real, full Bottomly design -- ratio ~0.91)
#   - moderate:   5 B6 vs 11 D2 (ratio ~0.45)
#   - severe:     2 B6 vs 11 D2 (ratio ~0.18)
# 2 is the floor, not 1: every tool here (DESeq2's dispersion estimate,
# edgeR's, limma-voom's within-group variance in lmFit) needs at least 2
# replicates in a group to estimate any within-group variance at all. n=1
# wouldn't be a harder stress point, it would be a different, degenerate
# question (can these tools even run at all), which isn't what this
# regime is checking.
#
# B6 samples for the "moderate" and "severe" levels are chosen as the
# first N of the full B6 sample set in sorted (lexicographic) order --
# computed from the metadata actually loaded below, not a remembered or
# guessed list, and printed in full for every level before anything is
# run on it. Sorted order is a neutral, reproducible selection rule, not
# one tuned toward any particular outcome.
#
# Usage:
#   Rscript benchmarks/run_stress_unbalanced.R
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
results_dir <- file.path("benchmarks", "results", "stress_unbalanced")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(data_dir, "bottomly_counts.csv")
metadata_path <- file.path(data_dir, "bottomly_metadata.csv")
if (!file.exists(counts_path) || !file.exists(metadata_path)) {
  stop(sprintf(
    paste(
      "Missing %s and/or %s. Run `Rscript benchmarks/run_deseq2_edger_bottomly.R`",
      "first -- this script builds imbalance levels from Bottomly's real counts",
      "and doesn't reacquire them itself."
    ),
    counts_path, metadata_path
  ))
}

counts_full <- read.csv(counts_path, row.names = "gene_id", check.names = FALSE)
counts_full <- as.matrix(counts_full)
meta_full <- read.csv(metadata_path)

cat(sprintf("Loaded bottomly_counts.csv: %d genes x %d samples.\n", nrow(counts_full), ncol(counts_full)))
cat(sprintf("Loaded bottomly_metadata.csv: %d samples.\n", nrow(meta_full)))

b6_all <- sort(meta_full$sample[meta_full$condition == "B6"])
d2_all <- sort(meta_full$sample[meta_full$condition == "D2"])

if (length(b6_all) != 10 || length(d2_all) != 11) {
  stop(sprintf(
    paste(
      "Expected 10 B6 / 11 D2 in bottomly_metadata.csv (the verified Bottomly",
      "design), found %d B6 / %d D2. Refusing to build imbalance levels from",
      "an unexpected group size -- inspect bottomly_metadata.csv before rerunning."
    ),
    length(b6_all), length(d2_all)
  ))
}

cat(sprintf("Full B6 (sorted, n=%d): %s\n", length(b6_all), paste(b6_all, collapse = ", ")))
cat(sprintf("Full D2 (sorted, n=%d): %s\n", length(d2_all), paste(d2_all, collapse = ", ")))

levels_spec <- list(
  balanced = length(b6_all),
  moderate = 5,
  severe = 2
)

# =====================================================================
# Per-level, per-tool run. Each tool's own standard, default workflow --
# identical to run_deseq2_edger_bottomly.R and run_limma_voom.R, just
# parameterized over which sample subset is fed in. Nothing about any
# tool's own settings changes across levels; only the input group sizes do.
# =====================================================================
run_deseq2 <- function(counts, coldata, out_path) {
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "D2", "B6"))
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df <- df[, c("gene_id", "log2FoldChange", "pvalue", "padj")]
  write.csv(df, out_path, row.names = FALSE)
  cat(sprintf("  DESeq2: %d genes tested, %d significant (padj < 0.05)\n",
              sum(!is.na(df$padj)), sum(df$padj < 0.05, na.rm = TRUE)))
}

run_edger <- function(counts, coldata, out_path) {
  y <- DGEList(counts = counts, group = coldata$condition)
  keep <- filterByExpr(y, group = coldata$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~condition, data = coldata)
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

run_limma <- function(counts, coldata, out_path) {
  y <- DGEList(counts = counts, group = coldata$condition)
  keep <- filterByExpr(y, group = coldata$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~condition, data = coldata)
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

for (level_name in names(levels_spec)) {
  n_b6 <- levels_spec[[level_name]]
  b6_subset <- b6_all[seq_len(n_b6)]
  d2_subset <- d2_all  # held fixed at full 11 across every level, by design

  cat(sprintf("\n===== Level: %s (%d B6 vs %d D2) =====\n", level_name, length(b6_subset), length(d2_subset)))
  cat(sprintf("B6 samples used: %s\n", paste(b6_subset, collapse = ", ")))
  cat(sprintf("D2 samples used: %s\n", paste(d2_subset, collapse = ", ")))

  sample_subset <- c(b6_subset, d2_subset)
  counts_subset <- counts_full[, sample_subset, drop = FALSE]
  coldata_subset <- data.frame(
    row.names = sample_subset,
    condition = factor(
      ifelse(sample_subset %in% b6_subset, "B6", "D2"),
      levels = c("B6", "D2")
    )
  )

  # Write the exact subset used, for anyone checking a level's result
  # against the samples that actually produced it.
  write.csv(
    data.frame(sample = rownames(coldata_subset), condition = coldata_subset$condition),
    file.path(results_dir, sprintf("%s_metadata.csv", level_name)),
    row.names = FALSE
  )

  run_deseq2(counts_subset, coldata_subset, file.path(results_dir, sprintf("%s_deseq2_results.csv", level_name)))
  run_edger(counts_subset, coldata_subset, file.path(results_dir, sprintf("%s_edger_results.csv", level_name)))
  run_limma(counts_subset, coldata_subset, file.path(results_dir, sprintf("%s_limma_voom_results.csv", level_name)))
}

cat("\nDone. Per-level metadata and DESeq2/edgeR/limma-voom results in benchmarks/results/stress_unbalanced/.\n")
cat("Next: python benchmarks/analyze_stress_unbalanced.py\n")
