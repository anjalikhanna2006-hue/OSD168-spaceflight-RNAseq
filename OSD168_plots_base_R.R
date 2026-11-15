# =============================================================================
# OSD-168 | PCA, volcano plot and heatmap using BASE R only
# -----------------------------------------------------------------------------
# Use this if OSD168_analysis.R finishes the DESeq2 part but crashes/blocks at
# the plotting step (e.g. Windows blocking cli.dll used by ggplot2).
#
# It needs NO extra packages. It reads the files that OSD168_analysis.R already
# saved in "results/", so run OSD168_analysis.R first.
#
# Needs:
#   data/GLDS-168_rna_seq_SampleTable.csv
#   results/normalized_counts.csv
#   results/DGE_FLT_vs_GC_all_genes.csv
# Writes: results/pca_plot.png, results/volcano_plot.png, results/heatmap_top_DEGs.png
#
# Note: PCA and heatmap here use log2(normalised counts + 1), not DESeq2's VST,
# so they look very similar but are not numerically identical to the main script.
# =============================================================================

data_dir <- "data"
out_dir  <- "results"
padj_cutoff <- 0.05
lfc_cutoff  <- 1
n_heatmap   <- 40

meta_file <- file.path(data_dir, "GLDS-168_rna_seq_SampleTable.csv")
norm_file <- file.path(out_dir,  "normalized_counts.csv")
dge_file  <- file.path(out_dir,  "DGE_FLT_vs_GC_all_genes.csv")
for (f in c(meta_file, norm_file, dge_file)) if (!file.exists(f)) stop("File not found: ", f)


# ---- Load --------------------------------------------------------------------
meta <- read.csv(meta_file, row.names = 1, check.names = FALSE)
norm <- read.csv(norm_file, row.names = 1, check.names = FALSE)
res  <- read.csv(dge_file, check.names = FALSE)

meta <- meta[colnames(norm), , drop = FALSE]
grp  <- factor(make.names(meta$condition))
names(grp) <- rownames(meta)

flt_level <- grep("^Space\\.Flight.*RR1\\.+With\\.+Spike", levels(grp), value = TRUE)
gc_level  <- grep("^Ground\\.Control.*RR1\\.+With\\.+Spike", levels(grp), value = TRUE)
if (length(flt_level) != 1 || length(gc_level) != 1) stop("Could not identify the Space Flight / Ground Control groups.")

logn <- log2(as.matrix(norm) + 1)


# ---- 1. PCA (all samples, top 500 most variable genes) ----------------------
top_var <- order(apply(logn, 1, var), decreasing = TRUE)[1:500]
pca     <- prcomp(t(logn[top_var, ]), scale. = FALSE)
pct     <- round(100 * summary(pca)$importance[2, 1:2])

pal <- hcl.colors(nlevels(grp), "Dark 3")

png(file.path(out_dir, "pca_plot.png"), width = 2400, height = 2400, res = 300)
layout(matrix(1:2, ncol = 1), heights = c(3, 1.3))
par(mar = c(4.5, 4.5, 3, 1))
plot(pca$x[, 1], pca$x[, 2], col = pal[as.integer(grp)], pch = 19, cex = 1.6,
     xlab = paste0("PC1: ", pct[1], "% variance"),
     ylab = paste0("PC2: ", pct[2], "% variance"),
     main = "PCA: OSD-168 mouse liver RNA-seq")
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = gsub("\\.+", " ", levels(grp)), col = pal, pch = 19,
       cex = 0.6, ncol = 2, bty = "n")
dev.off()


# ---- 2. Volcano plot ---------------------------------------------------------
v   <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
sig <- v$padj < padj_cutoff & abs(v$log2FoldChange) > lfc_cutoff
up  <- sig & v$log2FoldChange > 0
dn  <- sig & v$log2FoldChange < 0
ypl <- pmin(-log10(v$padj), 300)          # cap extreme values for plotting

cols <- rep("grey70", nrow(v))
cols[up] <- "#e74c3c"
cols[dn] <- "#3498db"

png(file.path(out_dir, "volcano_plot.png"), width = 2400, height = 2100, res = 300)
par(mar = c(4.5, 4.5, 3, 1))
plot(v$log2FoldChange, ypl, col = cols, pch = 19, cex = 0.5,
     xlab = "log2 fold change (Spaceflight vs Ground Control)",
     ylab = "-log10 adjusted p-value",
     main = "Spaceflight vs Ground Control (OSD-168)")
abline(v = c(-lfc_cutoff, lfc_cutoff), h = -log10(padj_cutoff), lty = 2, col = "grey40")
lab <- order(v$padj)[1:15]
text(v$log2FoldChange[lab], ypl[lab], labels = v$symbol[lab], cex = 0.5, pos = 3)
legend("topright", bty = "n", pch = 19, cex = 0.8,
       col = c("#e74c3c", "#3498db", "grey70"),
       legend = c(paste0("Up (", sum(up), ")"), paste0("Down (", sum(dn), ")"), "Not significant"))
dev.off()


# ---- 3. Heatmap: top DEGs, Space Flight vs Ground Control only --------------
top_genes <- head(v$gene[sig][order(v$padj[sig])], n_heatmap)
fg <- names(grp)[grp %in% c(flt_level, gc_level)]

mat <- logn[top_genes, fg, drop = FALSE]
mat <- t(scale(t(mat)))                    # z-score per gene
mat[is.nan(mat)] <- 0
rownames(mat) <- make.unique(v$symbol[match(top_genes, v$gene)])

side <- ifelse(grp[fg] == flt_level, "#e74c3c", "#3498db")

png(file.path(out_dir, "heatmap_top_DEGs.png"), width = 2400, height = 2800, res = 300)
heatmap(mat, scale = "none", ColSideColors = side, labCol = rep("", ncol(mat)),
        col = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
        margins = c(2, 9), cexRow = 0.7,
        main = paste0("Top ", length(top_genes), " DEGs (red bar = spaceflight, blue = ground)"))
dev.off()

cat("Done. Plots saved in the '", out_dir, "' folder.\n", sep = "")
