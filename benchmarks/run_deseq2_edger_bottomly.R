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
# Bottomly et al. 2011: mouse striatum RNA-seq, two inbred strains --
# C57BL/6J ("B6") vs DBA/2J ("D2"), 21 samples total (10 B6, 11 D2). A
# real, published, widely-used DE-methods benchmark dataset (used in
# Soneson & Delorenzi 2013 and Rapaport et al. 2013, among others).
#
# ACQUISITION ROUTE CHANGED FROM dexus TO recount3 (Aug 2026):
# The original version of this script loaded `dexus::countsBottomly`.
# `dexus` left Bioconductor after release 3.12 and is not installable on
# current Bioconductor (3.23 at time of writing) without downgrading R/
# Bioconductor, which was explicitly ruled out. This version instead
# pulls the same underlying study straight from SRA via `recount3`
# (organism = "mouse", data_source = "sra", project = "SRP004777" --
# verified directly against recount3's own hosted files and against
# NCBI's SRA listing for SRP004777 before writing this script: 21 runs,
# sample titles following a `{strain}_{replicate}_{lane}` pattern, e.g.
# B6_1_7 / D2_2_4, consistent with the 10 B6 / 11 D2 design).
#
# IMPORTANT COUNT-UNIT DIFFERENCE FROM THE OLD ROUTE: recount3's
# gene-level `raw_counts` assay is base-pair coverage (AUC-based), NOT
# literal read counts -- unlike `countsBottomly`, which was already a
# literal read-count matrix. Feeding `raw_counts` directly into DESeq2/
# edgeR would be statistically wrong. recount3's own documentation
# prescribes `transform_counts(rse, by = "auc", targetSize = 4e7, L =
# 100, round = TRUE)` to produce scaled, rounded (integer-like) counts
# suitable for count-based DE tools, and that is applied below before
# anything touches DESeq2 or edgeR. This is a genuine methodological
# difference from the dexus-based run, not a cosmetic one, and is worth
# remembering if Bottomly's numbers are ever compared against a
# hypothetical dexus-based rerun.
#
# Group membership is derived from the RSE's own colData() at runtime,
# NOT hardcoded from a remembered sample-ID list -- exactly the thing
# that went wrong for zebrafish (Trt1/Trt3/Trt5 guessed instead of the
# real Trt9/Trt11/Trt13). Since the exact colData column carrying the
# B6/D2 strain label isn't something this script's author could verify
# without R access, the script does not assume a specific column name:
# it scans every colData() column for one that classifies all 21 samples
# unambiguously into a B6/D2 pattern, prints every column name plus the
# resolved sample identifiers either way, and refuses to proceed (stop())
# unless exactly one such column is found and it yields exactly 10 B6 and
# 11 D2 -- with a full metadata dump on failure, since that's the only
# thing that lets whoever runs this catch a wrong guess before it
# silently produces a wrong benchmark instead of after.
#
# Usage:
#   Rscript benchmarks/run_deseq2_edger_bottomly.R
#
# Requires: BiocManager, DESeq2, edgeR, recount3 (installed automatically
# below if missing).

options(warn = 1)

required_bioc <- c("DESeq2", "edgeR", "recount3")
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
  library(recount3)
  library(SummarizedExperiment)
})

data_dir <- file.path("benchmarks", "data")
results_dir <- file.path("benchmarks", "results")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# --- locate and load the project ---
cat("Querying recount3 for available mouse SRA projects ...\n")
mouse_projects <- available_projects(organism = "mouse")

proj_info <- subset(
  mouse_projects,
  project == "SRP004777" & project_type == "data_sources" & organism == "mouse"
)

if (nrow(proj_info) != 1) {
  stop(sprintf(
    paste(
      "Expected exactly 1 matching recount3 project row for SRP004777",
      "(mouse, data_sources), found %d. This means either the project",
      "moved/changed in recount3's index, or the filter above is wrong.",
      "Inspect available_projects(organism = \"mouse\") and its 'project'/",
      "'project_type'/'organism' columns before proceeding."
    ),
    nrow(proj_info)
  ))
}

cat("Found project row:\n")
print(proj_info)

cat("Creating RSE (this downloads gene-level coverage counts + metadata) ...\n")
rse_gene <- create_rse(proj_info)

sample_names <- colnames(rse_gene)
cat(sprintf("Loaded rse_gene: %d genes x %d samples.\n", nrow(rse_gene), ncol(rse_gene)))
cat("Resolved sample identifiers (colnames(rse_gene)):\n")
print(sample_names)

# --- print every colData() column name before doing anything with it ---
cd <- as.data.frame(colData(rse_gene))
cat(sprintf("colData(rse_gene) has %d columns:\n", ncol(cd)))
print(colnames(cd))

# --- find whichever column(s) classify all 21 samples unambiguously into B6/D2 ---
# Deliberately not assuming a specific recount3 column name (e.g.
# "sra.sample_title") -- scan every column as character text and keep
# only the ones where every single sample matches exactly one of the two
# strain patterns. This is the runtime-discovery equivalent of the
# dexus script's column-name regex match, just applied across all
# available metadata columns instead of one assumed one.
classify_column <- function(values) {
  values <- as.character(values)
  is_b6 <- grepl("b6|c57", values, ignore.case = TRUE)
  is_d2 <- grepl("d2|dba", values, ignore.case = TRUE)
  list(is_b6 = is_b6, is_d2 = is_d2)
}

candidate_cols <- character(0)
candidate_conditions <- list()

for (col in colnames(cd)) {
  cls <- classify_column(cd[[col]])
  unambiguous <- !(cls$is_b6 & cls$is_d2) & (cls$is_b6 | cls$is_d2)
  if (all(unambiguous)) {
    candidate_cols <- c(candidate_cols, col)
    candidate_conditions[[col]] <- ifelse(cls$is_b6, "B6", "D2")
  }
}

if (length(candidate_cols) == 0) {
  cat("No colData() column unambiguously classified all samples as B6 or D2.\n")
  cat("Full colData() dump for diagnosis:\n")
  print(cd)
  stop(
    "Could not recover B6/D2 group labels from any recount3 colData() column ",
    "for SRP004777. Inspect the colData() dump above (and consider ",
    "sra.sample_attributes / sra.experiment_title / sra.sample_title-style ",
    "columns specifically) and adjust the classification logic above -- do ",
    "not hardcode a remembered B6/D2 sample list instead."
  )
}

cat(sprintf(
  "Column(s) that unambiguously classify all samples as B6/D2: %s\n",
  paste(candidate_cols, collapse = ", ")
))

# If multiple columns qualify, they must all agree with each other --
# disagreement means the pattern match is picking up something other than
# genuine strain identity (e.g. a free-text field that happens to contain
# "d2" for an unrelated reason), and that has to stop the script rather
# than pick one arbitrarily.
reference <- candidate_conditions[[candidate_cols[1]]]
for (col in candidate_cols[-1]) {
  if (!identical(candidate_conditions[[col]], reference)) {
    cat("Full colData() dump for diagnosis:\n")
    print(cd)
    stop(sprintf(
      paste(
        "Columns '%s' and '%s' both look like B6/D2 classifiers but disagree",
        "on at least one sample's strain. Refusing to guess which is correct --",
        "inspect the colData() dump above."
      ),
      candidate_cols[1], col
    ))
  }
}

condition <- reference
n_b6 <- sum(condition == "B6")
n_d2 <- sum(condition == "D2")
cat(sprintf(
  "Classified from column '%s': %d B6, %d D2 (verified design: 10 B6, 11 D2).\n",
  candidate_cols[1], n_b6, n_d2
))
cat("Per-sample classification:\n")
print(data.frame(sample = sample_names, condition = condition))

if (n_b6 != 10 || n_d2 != 11) {
  cat("Full colData() dump for diagnosis:\n")
  print(cd)
  stop(sprintf(
    paste(
      "Sample counts after classification (%d B6, %d D2) do not match the",
      "verified Bottomly design (10 B6, 11 D2). Refusing to proceed with a",
      "possibly-wrong grouping. Sample identifiers were: %s"
    ),
    n_b6, n_d2, paste(sample_names, collapse = ", ")
  ))
}

coldata <- data.frame(row.names = sample_names, condition = factor(condition, levels = c("B6", "D2")))

# --- AUC coverage counts -> DESeq2/edgeR-appropriate counts ---
# See the header comment: raw_counts here is base-pair coverage, not read
# counts. This is the recount3-documented conversion, applied before the
# transformed assay is ever written to disk or handed to DESeq2/edgeR.
cat("Converting AUC-based raw_counts to scaled read counts via transform_counts() ...\n")
assay(rse_gene, "counts") <- transform_counts(
  rse_gene,
  by = "auc",
  targetSize = 4e7,
  L = 100,
  round = TRUE
)
counts <- assay(rse_gene, "counts")
storage.mode(counts) <- "integer"

# --- write raw inputs so DEConcord can be run on the identical matrix ---
counts_out <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
write.csv(counts_out, file.path(data_dir, "bottomly_counts.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(coldata), condition = coldata$condition),
          file.path(data_dir, "bottomly_metadata.csv"), row.names = FALSE)

cat(sprintf(
  "bottomly: %d genes x %d samples (post AUC->count transform). Groups -- B6: %s | D2: %s\n",
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
