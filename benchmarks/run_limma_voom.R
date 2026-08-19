#!/usr/bin/env Rscript
#
# Runs limma-voom on the `airway` and `pasilla` datasets -- the same two
# datasets and the same group comparisons `run_deseq2_edger.R` and
# `run_deseq2_edger_pasilla.R` already run DESeq2 and edgeR on -- using
# limma's own standard, idiomatic voom workflow. Not tuned or nudged
# toward agreeing with DEConcord, DESeq2, or edgeR in any way. A third,
# independently developed tool to check concordance against, the same
# role DESeq2 and edgeR already play for each other (see
# `benchmarks/method_concordance.py`).
#
# Data loading mirrors the two existing scripts exactly (same count
# matrix, same sample metadata, same group factor levels and contrast
# direction), duplicated here rather than shared so this script stays
# runnable on its own without depending on the other two having run
# first. It does not rewrite `benchmarks/data/` -- that's already
# `run_deseq2_edger.R`/`run_deseq2_edger_pasilla.R`'s job, and the raw
# counts would be byte-for-byte identical.
#
# Usage:
#   Rscript benchmarks/run_limma_voom.R
#
# Requires: BiocManager, limma, edgeR, airway, pasilla (installed
# automatically below if missing -- comment that block out if you'd
# rather manage your own R library).

options(warn = 1)

required_bioc <- c("limma", "edgeR", "airway", "pasilla")
missing_pkgs <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
}

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(airway)
  library(pasilla)
})

results_dir <- file.path("benchmarks", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# =====================================================================
# limma-voom -- standard workflow: filterByExpr, TMM, voom, lmFit, eBayes
# =====================================================================
run_limma_voom <- function(counts, group, design, coef, out_path, label) {
  y <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(y, group = group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  v <- voom(y, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  voom_df <- topTable(fit, coef = coef, n = Inf, sort.by = "none")
  voom_df$gene_id <- rownames(voom_df)
  voom_df <- voom_df[, c("gene_id", "logFC", "P.Value", "adj.P.Val")]
  colnames(voom_df) <- c("gene_id", "logFC", "pvalue", "adj.P.Val")
  write.csv(voom_df, out_path, row.names = FALSE)
  cat(sprintf("limma-voom (%s): %d genes tested (post filterByExpr), %d significant (adj.P.Val < 0.05)\n",
              label, nrow(voom_df), sum(voom_df$adj.P.Val < 0.05, na.rm = TRUE)))
}

# --- airway: same loading as run_deseq2_edger.R ---
data(airway)
airway_counts <- assay(airway)
airway_meta <- as.data.frame(colData(airway))[, c("SampleName", "cell", "dex")]
airway_meta$dex <- factor(airway_meta$dex, levels = c("untrt", "trt"))

cat(sprintf(
  "airway: %d genes x %d samples. Groups -- untrt: %s | trt: %s\n",
  nrow(airway_counts), ncol(airway_counts),
  paste(colnames(airway_counts)[airway_meta$dex == "untrt"], collapse = ", "),
  paste(colnames(airway_counts)[airway_meta$dex == "trt"], collapse = ", ")
))

airway_design <- model.matrix(~dex, data = airway_meta)
run_limma_voom(
  airway_counts, airway_meta$dex, airway_design, coef = "dextrt",
  out_path = file.path(results_dir, "limma_voom_results.csv"), label = "airway"
)

# --- pasilla: same loading as run_deseq2_edger_pasilla.R ---
pas_cts_path <- system.file("extdata", "pasilla_gene_counts.tsv", package = "pasilla", mustWork = TRUE)
pas_anno_path <- system.file("extdata", "pasilla_sample_annotation.csv", package = "pasilla", mustWork = TRUE)

pasilla_counts <- as.matrix(read.csv(pas_cts_path, sep = "\t", row.names = "gene_id"))
pasilla_coldata <- read.csv(pas_anno_path, row.names = 1)
pasilla_coldata <- pasilla_coldata[, c("condition", "type")]

rownames(pasilla_coldata) <- sub("fb$", "", rownames(pasilla_coldata))
colnames(pasilla_counts) <- sub("^X", "", colnames(pasilla_counts))
pasilla_counts <- pasilla_counts[, rownames(pasilla_coldata)]

pasilla_coldata$condition <- factor(pasilla_coldata$condition, levels = c("untreated", "treated"))

cat(sprintf(
  "pasilla: %d genes x %d samples. Groups -- untreated: %s | treated: %s\n",
  nrow(pasilla_counts), ncol(pasilla_counts),
  paste(rownames(pasilla_coldata)[pasilla_coldata$condition == "untreated"], collapse = ", "),
  paste(rownames(pasilla_coldata)[pasilla_coldata$condition == "treated"], collapse = ", ")
))

pasilla_design <- model.matrix(~condition, data = pasilla_coldata)
run_limma_voom(
  pasilla_counts, pasilla_coldata$condition, pasilla_design, coef = "conditiontreated",
  out_path = file.path(results_dir, "pasilla_limma_voom_results.csv"), label = "pasilla"
)

cat("\nDone. Results in benchmarks/results/.\n")
cat("Next: hand limma_voom_results.csv / pasilla_limma_voom_results.csv to compute_de_concordance\n")
cat("alongside the matching deseq2_results.csv / edger_results.csv, e.g.:\n")
cat('  compute_de_concordance(deseq2, limma_voom, lfc_col_a="log2FoldChange", pvalue_col_a="padj",\n')
cat('                          lfc_col_b="logFC", pvalue_col_b="adj.P.Val")\n')
