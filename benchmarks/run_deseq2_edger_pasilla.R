#!/usr/bin/env Rscript
#
# Second, independent benchmark dataset -- same purpose as
# run_deseq2_edger.R, different organism and design, to check whether the
# precision/recall pattern found on `airway` is a property of DEConcord's
# method or just a property of that one dataset.
#
# `pasilla`: Brooks et al. 2011 (GEO GSE18508) -- Drosophila melanogaster
# S2-DRSC cells, RNAi knockdown of the *pasilla* gene (the fly homolog of
# human NOVA1/2) vs. untreated control, 7 samples total (3 treated, 4
# untreated). This is the dataset the DESeq2 vignette itself uses for its
# "count matrix input" walkthrough, so the loading code below follows that
# vignette's own pattern deliberately -- group labels come from the
# package's own sample annotation file, not from guessing column-name
# suffixes.
#
# Usage:
#   Rscript benchmarks/run_deseq2_edger_pasilla.R
#
# Requires: BiocManager, DESeq2, edgeR, pasilla (installed automatically
# below if missing).

options(warn = 1)

required_bioc <- c("DESeq2", "edgeR", "pasilla")
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
  library(pasilla)
})

data_dir <- file.path("benchmarks", "data")
results_dir <- file.path("benchmarks", "results")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- load counts + sample annotation the same way the DESeq2 vignette does ---
pas_cts_path <- system.file("extdata", "pasilla_gene_counts.tsv", package = "pasilla", mustWork = TRUE)
pas_anno_path <- system.file("extdata", "pasilla_sample_annotation.csv", package = "pasilla", mustWork = TRUE)

counts <- as.matrix(read.csv(pas_cts_path, sep = "\t", row.names = "gene_id"))
coldata <- read.csv(pas_anno_path, row.names = 1)
coldata <- coldata[, c("condition", "type")]

# Sample names differ slightly between the two files (counts columns have an
# "X" prefix from read.csv on numeric-leading names; annotation rownames have
# an "fb" suffix) -- strip both so they line up, then reorder to match.
rownames(coldata) <- sub("fb$", "", rownames(coldata))
colnames(counts) <- sub("^X", "", colnames(counts))
counts <- counts[, rownames(coldata)]

coldata$condition <- factor(coldata$condition, levels = c("untreated", "treated"))

# --- write raw inputs so DEConcord can be run on the identical matrix ---
counts_out <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
write.csv(counts_out, file.path(data_dir, "pasilla_counts.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(coldata), condition = coldata$condition, type = coldata$type),
          file.path(data_dir, "pasilla_metadata.csv"), row.names = FALSE)

cat(sprintf(
  "pasilla: %d genes x %d samples. Groups -- untreated: %s | treated: %s\n",
  nrow(counts), ncol(counts),
  paste(rownames(coldata)[coldata$condition == "untreated"], collapse = ", "),
  paste(rownames(coldata)[coldata$condition == "treated"], collapse = ", ")
))

# =====================================================================
# DESeq2 -- standard workflow, default settings, contrast treated vs untreated
# =====================================================================
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
dds <- DESeq(dds)
deseq2_res <- results(dds, contrast = c("condition", "treated", "untreated"))
deseq2_df <- as.data.frame(deseq2_res)
deseq2_df$gene_id <- rownames(deseq2_df)
deseq2_df <- deseq2_df[, c("gene_id", "log2FoldChange", "pvalue", "padj")]
write.csv(deseq2_df, file.path(results_dir, "pasilla_deseq2_results.csv"), row.names = FALSE)
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
qlf <- glmQLFTest(fit, coef = "conditiontreated")
edger_df <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none"))
edger_df$gene_id <- rownames(edger_df)
edger_df <- edger_df[, c("gene_id", "logFC", "PValue", "FDR")]
colnames(edger_df) <- c("gene_id", "logFC", "pvalue", "FDR")
write.csv(edger_df, file.path(results_dir, "pasilla_edger_results.csv"), row.names = FALSE)
cat(sprintf("edgeR: %d genes tested (post filterByExpr), %d significant (FDR < 0.05)\n",
            nrow(edger_df), sum(edger_df$FDR < 0.05, na.rm = TRUE)))

cat("\nDone. Raw inputs in benchmarks/data/, results in benchmarks/results/.\n")
cat("Next: python benchmarks/compare_results.py --dataset pasilla\n")
