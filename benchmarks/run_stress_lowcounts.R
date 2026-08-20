#!/usr/bin/env Rscript
#
# Third "harder data regime" in the stress-test suite: low counts. Same
# framing as run_stress_unbalanced.R and run_stress_batch.R -- DEConcord's
# core method-concordance question, applied outside the two clean,
# balanced, well-sequenced two-group datasets checked so far.
#
# No expression values are fabricated. Depth is reduced by binomial
# thinning: each real read in the real Bottomly count matrix is kept
# with a fixed probability p and dropped otherwise
# (rbinom(1, real_count, p) per gene per sample), the standard,
# well-established way to simulate a shallower sequencing run from real
# data (the same logic behind rarefaction and tools like
# DropletUtils::downsampleMatrix) -- it only ever removes real reads,
# never invents a signal that wasn't there. Three depths: full (p = 1,
# the real, unmodified data -- the reference point), moderate (p = 0.25),
# low (p = 0.05). A fixed RNG seed makes the thinning reproducible.
#
# Runs on Bottomly's real, full "balanced" design (10 B6, 11 D2) -- same
# real counts benchmarks/run_deseq2_edger_bottomly.R and
# run_stress_unbalanced.R's "balanced" level already use, not a new
# acquisition. Sample sizes and groups are held constant across depths
# here, deliberately -- this regime isolates sequencing depth as the one
# variable that changes, the same way run_stress_unbalanced.R isolated
# group-size ratio by holding D2 fixed.
#
# Usage:
#   Rscript benchmarks/run_stress_lowcounts.R
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
results_dir <- file.path("benchmarks", "results", "stress_lowcounts")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(data_dir, "bottomly_counts.csv")
metadata_path <- file.path(data_dir, "bottomly_metadata.csv")
if (!file.exists(counts_path) || !file.exists(metadata_path)) {
  stop(sprintf(
    paste(
      "Missing %s and/or %s. Run `Rscript benchmarks/run_deseq2_edger_bottomly.R`",
      "first -- this script downsamples Bottomly's real counts and doesn't",
      "reacquire them itself."
    ),
    counts_path, metadata_path
  ))
}

counts_full <- read.csv(counts_path, row.names = "gene_id", check.names = FALSE)
counts_full <- as.matrix(counts_full)
meta <- read.csv(metadata_path)
meta <- meta[match(colnames(counts_full), meta$sample), ]
if (any(is.na(meta$sample)) || !identical(meta$sample, colnames(counts_full))) {
  stop("bottomly_metadata.csv's samples don't match bottomly_counts.csv's columns after reordering -- inspect both files.")
}
if (sum(meta$condition == "B6") != 10 || sum(meta$condition == "D2") != 11) {
  stop(sprintf(
    "Expected 10 B6 / 11 D2 (the verified Bottomly design), found %d B6 / %d D2. Refusing to proceed.",
    sum(meta$condition == "B6"), sum(meta$condition == "D2")
  ))
}

cat(sprintf("Loaded bottomly_counts.csv: %d genes x %d samples.\n", nrow(counts_full), ncol(counts_full)))
cat(sprintf("Total real reads (sum of all counts): %.0f\n", sum(as.numeric(counts_full))))

coldata <- data.frame(row.names = colnames(counts_full), condition = factor(meta$condition, levels = c("B6", "D2")))

set.seed(20260820)  # fixed seed -- thinning is reproducible, not re-rollable toward a preferred outcome

depths <- list(full = 1.0, moderate = 0.25, low = 0.05)

thin_counts <- function(counts, p) {
  if (p >= 1.0) return(counts)  # p = 1 is the real, unmodified matrix -- no binomial draw needed or wanted
  thinned <- matrix(
    rbinom(n = length(counts), size = as.integer(counts), prob = p),
    nrow = nrow(counts), ncol = ncol(counts),
    dimnames = dimnames(counts)
  )
  storage.mode(thinned) <- "integer"
  thinned
}

# =====================================================================
# Per-depth, per-tool run. Same standard, default workflow each tool
# already uses elsewhere in this repo -- only the input depth changes.
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

for (depth_name in names(depths)) {
  p <- depths[[depth_name]]
  counts_thinned <- thin_counts(counts_full, p)
  total_reads <- sum(as.numeric(counts_thinned))

  cat(sprintf("\n===== Depth: %s (p = %.2f) =====\n", depth_name, p))
  cat(sprintf("Total reads after thinning: %.0f (%.1f%% of real total)\n",
              total_reads, 100 * total_reads / sum(as.numeric(counts_full))))

  write.csv(
    data.frame(gene_id = rownames(counts_thinned), counts_thinned, check.names = FALSE),
    file.path(results_dir, sprintf("%s_counts.csv", depth_name)),
    row.names = FALSE
  )

  run_deseq2(counts_thinned, coldata, file.path(results_dir, sprintf("%s_deseq2_results.csv", depth_name)))
  run_edger(counts_thinned, coldata, file.path(results_dir, sprintf("%s_edger_results.csv", depth_name)))
  run_limma(counts_thinned, coldata, file.path(results_dir, sprintf("%s_limma_voom_results.csv", depth_name)))
}

cat("\nDone. Per-depth counts and DESeq2/edgeR/limma-voom results in benchmarks/results/stress_lowcounts/.\n")
cat("Next: python benchmarks/analyze_stress_lowcounts.py\n")
