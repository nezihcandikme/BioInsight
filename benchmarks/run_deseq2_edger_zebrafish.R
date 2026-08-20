#!/usr/bin/env Rscript
#
# Third, independent benchmark dataset -- same purpose as run_deseq2_edger.R
# and run_deseq2_edger_pasilla.R, a different organism and experimental
# design again, specifically to give the research/ml/ rank-stability
# finding (see research/ml/FINDINGS.md) a genuinely external test: airway
# and pasilla are both "explant/cell line plus a drug or RNAi knockdown"
# designs. This one is neither.
#
# `zfGenes` (zebrafishRNASeq package): Ferreira et al., Danio rerio
# (zebrafish), gsk3 knockdown vs. wild-type control, 6 samples (3 control:
# Ctl1/Ctl3/Ctl5, 3 knockdown: Trt9/Trt11/Trt13). This is the exact dataset
# the edgeR User's Guide itself uses for its "count matrix input with ERCC
# spike-ins" walkthrough, so the loading and filtering steps below follow
# that guide's own pattern deliberately, same reasoning run_deseq2_edger_
# pasilla.R gives for following the DESeq2 vignette's pattern on pasilla.
#
# One real difference from airway/pasilla worth naming: this count matrix
# includes 92 ERCC spike-in control rows (synthetic RNA added at a known
# concentration, not real zebrafish genes) alongside the ~32,000 real
# genes. These are dropped before DE analysis, same as the edgeR guide
# does, since a spike-in isn't a biological gene and including it would
# make "genes tested" and "genes significant" counts mean something
# different from every other dataset this project has benchmarked. Row
# count and dropped-row count are both printed below so this is checkable,
# not just asserted.
#
# Usage:
#   Rscript benchmarks/run_deseq2_edger_zebrafish.R
#
# Requires: BiocManager, DESeq2, edgeR, zebrafishRNASeq (installed
# automatically below if missing).

options(warn = 1)

required_bioc <- c("DESeq2", "edgeR", "zebrafishRNASeq")
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
  library(zebrafishRNASeq)
})

data_dir <- file.path("benchmarks", "data")
results_dir <- file.path("benchmarks", "results")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- load counts, drop ERCC spike-ins, build sample metadata ---
data(zfGenes)
counts_raw <- as.matrix(zfGenes)

is_ercc <- grepl("^ERCC", rownames(counts_raw))
cat(sprintf(
  "zfGenes: %d rows total, %d ERCC spike-ins dropped, %d real genes kept.\n",
  nrow(counts_raw), sum(is_ercc), sum(!is_ercc)
))
counts <- counts_raw[!is_ercc, , drop = FALSE]

sample_names <- colnames(counts)
condition <- ifelse(grepl("^Ctl", sample_names), "control", "knockdown")
coldata <- data.frame(row.names = sample_names, condition = factor(condition, levels = c("control", "knockdown")))

# --- write raw inputs so DEConcord can be run on the identical matrix ---
counts_out <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
write.csv(counts_out, file.path(data_dir, "zebrafish_counts.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(coldata), condition = coldata$condition),
          file.path(data_dir, "zebrafish_metadata.csv"), row.names = FALSE)

cat(sprintf(
  "zebrafish: %d genes x %d samples. Groups -- control: %s | knockdown: %s\n",
  nrow(counts), ncol(counts),
  paste(rownames(coldata)[coldata$condition == "control"], collapse = ", "),
  paste(rownames(coldata)[coldata$condition == "knockdown"], collapse = ", ")
))

# =====================================================================
# DESeq2 -- standard workflow, default settings, contrast knockdown vs control
# =====================================================================
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
dds <- DESeq(dds)
deseq2_res <- results(dds, contrast = c("condition", "knockdown", "control"))
deseq2_df <- as.data.frame(deseq2_res)
deseq2_df$gene_id <- rownames(deseq2_df)
deseq2_df <- deseq2_df[, c("gene_id", "log2FoldChange", "pvalue", "padj")]
write.csv(deseq2_df, file.path(results_dir, "zebrafish_deseq2_results.csv"), row.names = FALSE)
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
qlf <- glmQLFTest(fit, coef = "conditionknockdown")
edger_df <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none"))
edger_df$gene_id <- rownames(edger_df)
edger_df <- edger_df[, c("gene_id", "logFC", "PValue", "FDR")]
colnames(edger_df) <- c("gene_id", "logFC", "pvalue", "FDR")
write.csv(edger_df, file.path(results_dir, "zebrafish_edger_results.csv"), row.names = FALSE)
cat(sprintf("edgeR: %d genes tested (post filterByExpr), %d significant (FDR < 0.05)\n",
            nrow(edger_df), sum(edger_df$FDR < 0.05, na.rm = TRUE)))

cat("\nDone. Raw inputs in benchmarks/data/, results in benchmarks/results/.\n")
cat("Next: python benchmarks/compare_results.py --dataset zebrafish\n")
