# =============================================================================
# OSD-168 | Spaceflight vs Ground Control | Mouse liver RNA-seq (DESeq2)
# -----------------------------------------------------------------------------
# Dataset : NASA Open Science Data Repository (OSDR), OSD-168 (formerly GLDS-168)
#           C57BL/6J mouse liver, ISS mission SpaceX CRS-4 (RR1)
# Question: How does spaceflight change liver gene expression?
# Contrast: Space Flight (RR1, with ERCC spike-in) vs Ground Control (RR1, with spike-in)
#
# HOW TO RUN
#   1. Download these files from OSDR (OSD-168) into a folder called "data/"
#      next to this script:
#        - GLDS-168_rna_seq_Unnormalized_Counts.csv   (required)
#        - GLDS-168_rna_seq_SampleTable.csv           (required)
#        - GLDS-168_rna_seq_differential_expression.csv (optional, only used
#          to get gene symbols for plot labels)
#   2. In RStudio: open this script, then
#      Session > Set Working Directory > To Source File Location
#   3. Click "Source" (or run it line by line).
#   All outputs are written to the "results/" folder.
# =============================================================================


# ---- 0. Packages -------------------------------------------------------------
# Installs anything missing (first run can take a few minutes).
cran_pkgs <- c("ggplot2", "pheatmap")
bioc_pkgs <- c("DESeq2", "EnhancedVolcano")

missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(EnhancedVolcano)
  library(pheatmap)
})


# ---- 1. Settings -------------------------------------------------------------
data_dir    <- "data"
out_dir     <- "results"

counts_file <- file.path(data_dir, "GLDS-168_rna_seq_Unnormalized_Counts.csv")
meta_file   <- file.path(data_dir, "GLDS-168_rna_seq_SampleTable.csv")
de_file     <- file.path(data_dir, "GLDS-168_rna_seq_differential_expression.csv")  # optional

padj_cutoff <- 0.05   # adjusted p-value threshold
lfc_cutoff  <- 1      # |log2 fold change| > 1  (i.e. at least 2-fold)
min_count   <- 10     # keep genes with total count >= 10 across all samples
n_heatmap   <- 40     # number of top DEGs shown in the heatmap

dir.create(out_dir, showWarnings = FALSE)

if (!file.exists(counts_file)) stop("Counts file not found: ", counts_file)
if (!file.exists(meta_file))   stop("Sample table not found: ", meta_file)


# ---- 2. Load data ------------------------------------------------------------
raw_counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
metadata   <- read.csv(meta_file,   row.names = 1, check.names = FALSE)

cat("Count matrix:", nrow(raw_counts), "genes x", ncol(raw_counts), "samples\n")

if (!"condition" %in% colnames(metadata)) {
  stop("No 'condition' column in the sample table. Columns found: ",
       paste(colnames(metadata), collapse = ", "))
}
if (!all(colnames(raw_counts) %in% rownames(metadata))) {
  stop("Some count-matrix samples are missing from the sample table.")
}

# Put metadata rows in the same order as count columns (DESeq2 requires this)
metadata <- metadata[colnames(raw_counts), , drop = FALSE]
stopifnot(identical(colnames(raw_counts), rownames(metadata)))

counts_matrix      <- round(as.matrix(raw_counts))   # DESeq2 needs integer counts
metadata$condition <- factor(make.names(metadata$condition))

cat("\nExperimental groups in this dataset:\n")
print(table(metadata$condition))


# ---- 3. Choose the two groups to compare -------------------------------------
# Picks the Space Flight and Ground Control groups (RR1, with spike-in) by name.
flt_level <- grep("^Space\\.Flight.*RR1\\.+With\\.+Spike", levels(metadata$condition), value = TRUE)
gc_level  <- grep("^Ground\\.Control.*RR1\\.+With\\.+Spike", levels(metadata$condition), value = TRUE)

if (length(flt_level) != 1 || length(gc_level) != 1) {
  stop("Could not identify exactly one Space Flight and one Ground Control group.\n",
       "Check the group names printed above and edit flt_level / gc_level manually.")
}
cat("\nComparing:\n  Space Flight   =", flt_level, "\n  Ground Control =", gc_level, "\n")


# ---- 4. DESeq2: normalisation + differential expression ----------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = metadata,
  design    = ~ condition
)

# Remove very low-count genes
dds <- dds[rowSums(counts(dds)) >= min_count, ]
cat("Genes kept after low-count filter:", nrow(dds), "\n")

# Median-of-ratios normalisation + model fit + Wald test
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(out_dir, "normalized_counts.csv"))

res <- results(dds, contrast = c("condition", flt_level, gc_level))
res <- res[order(res$padj), ]
summary(res)

res_df      <- as.data.frame(res)
res_df$gene <- rownames(res_df)


# ---- 5. Gene symbols (optional, for readable labels) -------------------------
# Falls back to the original gene IDs if the optional file is missing.
res_df$symbol <- res_df$gene

if (file.exists(de_file)) {
  hdr <- read.csv(de_file, nrows = 1, check.names = FALSE)
  if ("SYMBOL" %in% colnames(hdr)) {
    cls <- rep("NULL", ncol(hdr))                       # skip every column...
    cls[c(1, which(colnames(hdr) == "SYMBOL"))] <- NA   # ...except ID and SYMBOL
    map <- read.csv(de_file, colClasses = cls, check.names = FALSE)

    ids <- map[[1]]
    sym <- map[["SYMBOL"]]
    ok  <- !is.na(sym) & sym != ""
    lookup <- setNames(sym[ok], ids[ok])

    hit <- res_df$gene %in% names(lookup)
    res_df$symbol[hit] <- lookup[res_df$gene[hit]]
    cat("Gene symbols found for", sum(hit), "of", nrow(res_df), "genes\n")
  }
}


# ---- 6. Significant genes ----------------------------------------------------
sig <- subset(res_df, !is.na(padj) & padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff)
up   <- subset(sig, log2FoldChange > 0)
down <- subset(sig, log2FoldChange < 0)

cat("\n===== RESULTS: Space Flight vs Ground Control =====\n")
cat("Genes tested (non-NA padj):", sum(!is.na(res_df$padj)), "\n")
cat("Significant DEGs          :", nrow(sig),  "\n")
cat("  Upregulated in flight   :", nrow(up),   "\n")
cat("  Downregulated in flight :", nrow(down), "\n\n")

cols <- c("symbol", "log2FoldChange", "padj")
cat("Top 10 upregulated:\n");   print(head(up[order(-up$log2FoldChange), cols], 10))
cat("\nTop 10 downregulated:\n"); print(head(down[order(down$log2FoldChange), cols], 10))

write.csv(res_df, file.path(out_dir, "DGE_FLT_vs_GC_all_genes.csv"), row.names = FALSE)
write.csv(sig,    file.path(out_dir, "DGE_FLT_vs_GC_significant.csv"), row.names = FALSE)


# ---- 7. PCA plot (all samples) -----------------------------------------------
vsd      <- vst(dds, blind = FALSE)   # variance-stabilising transform
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3.5, alpha = 0.85) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  ggtitle("PCA: OSD-168 mouse liver RNA-seq") +
  guides(color = guide_legend(ncol = 2)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 7))

ggsave(file.path(out_dir, "pca_plot.png"), p_pca, width = 9, height = 7, dpi = 300)


# ---- 8. Volcano plot ---------------------------------------------------------
volc_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]  # already sorted by padj

p_volcano <- EnhancedVolcano(
  volc_df,
  lab            = volc_df$symbol,
  selectLab      = head(volc_df$symbol, 15),      # label the 15 most significant genes
  x              = "log2FoldChange",
  y              = "padj",
  pCutoff        = padj_cutoff,
  FCcutoff       = lfc_cutoff,
  title          = "Spaceflight vs Ground Control (OSD-168)",
  subtitle       = "Mouse liver RNA-seq",
  pointSize      = 2,
  labSize        = 3.5,
  legendPosition = "top"
)

ggsave(file.path(out_dir, "volcano_plot.png"), p_volcano, width = 9, height = 8, dpi = 300)


# ---- 9. Heatmap: top DEGs, Space Flight vs Ground Control only ---------------
top_genes  <- head(rownames(sig)[order(sig$padj)], n_heatmap)
fg_samples <- rownames(metadata)[metadata$condition %in% c(flt_level, gc_level)]

mat <- assay(vsd)[top_genes, fg_samples]
rownames(mat) <- make.unique(res_df[top_genes, "symbol"])

anno_col <- data.frame(
  Group = factor(ifelse(metadata[fg_samples, "condition"] == flt_level,
                        "Spaceflight", "Ground control")),
  row.names = fg_samples
)

pheatmap(
  mat,
  scale          = "row",                         # z-score per gene
  annotation_col = anno_col,
  show_colnames  = FALSE,
  fontsize_row   = 7,
  color          = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
  main           = paste0("Top ", length(top_genes), " DEGs: Spaceflight vs Ground Control"),
  filename       = file.path(out_dir, "heatmap_top_DEGs.png"),
  width          = 8,
  height         = 9
)


# ---- 10. Save session info (for reproducibility) -----------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nDone. All results saved in the '", out_dir, "' folder.\n", sep = "")
