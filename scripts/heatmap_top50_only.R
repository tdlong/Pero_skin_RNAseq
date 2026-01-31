#!/usr/bin/env Rscript
# Regenerate only the top-50 heatmap from saved DESeq2 objects.
# Run locally after downloading deseq2_results/ from the server.
# Usage: Rscript scripts/heatmap_top50_only.R [path_to_deseq2_results]
#   Default: deseq2_results (relative to current dir). Output: <that_path>/figures/heatmap_top50_infected.pdf
# Requires: pheatmap (install.packages("pheatmap"))

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else "deseq2_results"
saved_path <- file.path(results_dir, "saved", "deseq2_objects.rds")
fig_dir    <- file.path(results_dir, "figures")

if (!file.exists(saved_path)) {
  stop("Saved objects not found: ", saved_path, "\nRun from project root or pass path, e.g. Rscript scripts/heatmap_top50_only.R deseq2_results")
}

dir.create(fig_dir, showWarnings = FALSE)

library(DESeq2)
obj <- readRDS(saved_path)
vsd          <- obj$vsd
meta         <- obj$meta
res_infected <- obj$res_infected
loc_add_to_gene2 <- obj$loc_add_to_gene2

# Top 50 genes by padj (infected contrast)
top_n   <- min(50, nrow(res_infected))
top_genes <- rownames(res_infected)[order(res_infected$padj)][seq_len(top_n)]

mat <- assay(vsd)[top_genes, ]
mat <- t(scale(t(mat)))

# Order columns by treatment: infected males, infected females, uninfected males, uninfected females (no clustering by expression)
treatment_order <- factor(
  paste(meta[colnames(mat), "infected"], meta[colnames(mat), "sex"]),
  levels = c("infected Male", "infected Female", "uninfected Male", "uninfected Female")
)
ord <- order(treatment_order)
mat <- mat[, ord, drop = FALSE]

# Column annotation: two rows of squares (Infected, Sex) — no sample IDs on x-axis
ann_col <- as.data.frame(meta[colnames(mat), c("infected", "sex"), drop = FALSE])
colnames(ann_col) <- c("Infected", "Sex")
ann_colors <- list(
  Infected = c(uninfected = "grey95", infected = "darkred"),
  Sex      = c(Female = "red", Male = "blue")
)

# Row labels: display names when available
rownames(mat) <- if (length(loc_add_to_gene2) > 0) {
  ifelse(rownames(mat) %in% names(loc_add_to_gene2), loc_add_to_gene2[rownames(mat)], rownames(mat))
} else {
  rownames(mat)
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("Install pheatmap: install.packages(\"pheatmap\")")
}

# Shrink row font so all gene labels fit; x-axis = two annotation rows (Infected, Sex), no sample IDs
n_genes <- nrow(mat)
pdf(file.path(fig_dir, "heatmap_top50_infected.pdf"), width = 8, height = max(8, n_genes * 0.2))
pheatmap::pheatmap(mat,
                   scale = "none",
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   annotation_col = ann_col,
                   annotation_colors = ann_colors,
                   show_colnames = FALSE,
                   show_rownames = TRUE,
                   fontsize_row = 3,
                   col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                   main = "Top DE genes (infected vs uninfected)")
dev.off()

message("Saved: ", file.path(fig_dir, "heatmap_top50_infected.pdf"))
