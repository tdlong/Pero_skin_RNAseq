# DESeq2 differential expression analysis
# Uses count matrix from featureCounts and exp_design.txt.
# Run in R (interactive or Rscript) after 03_featurecounts.slurm.
#
# R/4.4.2: If DESeq2 is not installed:
#   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("DESeq2")
# Optional: install.packages("pheatmap") for heatmap column annotations (infected, sex)
#
# Experimental design: full model ~ sex + infected; sex is covariate, infection is effect of interest.
# Infection p-values from LRT: full (~ sex + infected) vs reduced (~ sex), so the test is "effect of infection" adjusted for sex.
# Writes tables (incl. results_infected_by_pvalue.csv), saves R objects and default plots in deseq2_results/saved and deseq2_results/figures.
# To make custom figures later: obj <- readRDS("deseq2_results/saved/deseq2_objects.rds"); then use obj$dds, obj$vsd, obj$res_infected, obj$normalized_counts, obj$meta, obj$loc_add_to_gene2.

library(DESeq2)

# Paths (adjust if not running from project root)
project_root <- "."
count_file  <- file.path(project_root, "counts/counts.txt")
meta_file   <- file.path(project_root, "exp_design.txt")
name_file   <- file.path(project_root, "ref/cds_Pleuc_names_v6.txt")  # loc_add (NCBI/GTF id) -> gene2 (bespoke name)
out_dir     <- file.path(project_root, "deseq2_results")
saved_dir   <- file.path(out_dir, "saved")   # R objects to load for custom figures
fig_dir     <- file.path(out_dir, "figures") # default plots
dir.create(out_dir, showWarnings = FALSE)
dir.create(saved_dir, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)

# Build loc_add -> gene2: GTF gene ids (loc_add) get display name gene2 when present; else keep NCBI name
if (file.exists(name_file)) {
  L <- read.delim(name_file, check.names = FALSE)
  L <- L[L$loc_add != "", c("gene", "loc_add")]
  L2 <- aggregate(gene ~ loc_add, data = L, FUN = function(x) paste(x, collapse = "_"))
  names(L2)[names(L2) == "gene"] <- "gene2"
  loc_add_to_gene2 <- setNames(L2$gene2, L2$loc_add)
} else {
  loc_add_to_gene2 <- character(0)
  message("Name file not found: ", name_file, " (skipping display names)")
}

# Design: sex + infected (main effects only). Infection p-values from LRT vs reduced model ~ sex.
DESIGN <- ~ sex + infected
DESIGN_REDUCED <- ~ sex   # reduced model for LRT; LRT tests effect of infection

# Read count matrix (featureCounts: first column is Geneid, then Length, then sample columns)
cnt <- read.table(count_file, header = TRUE, comment.char = "#", check.names = FALSE)
gene_id <- cnt$Geneid
# Sample columns: skip Geneid and Length (and Chr/Start/End/Strand if present)
sample_cols <- setdiff(colnames(cnt), c("Geneid", "Length", "Chr", "Start", "End", "Strand"))
count_matrix <- as.matrix(cnt[, sample_cols])
rownames(count_matrix) <- gene_id

# Clean sample names if they are full paths (e.g. bams/S-25666.sorted.bam -> S-25666)
sample_names <- basename(sample_cols)
sample_names <- sub("\\.sorted\\.bam$", "", sample_names)
colnames(count_matrix) <- sample_names

# Metadata: must have "animal" column matching count matrix column names
meta <- read.delim(meta_file, check.names = FALSE)
rownames(meta) <- meta$animal
meta <- meta[colnames(count_matrix), , drop = FALSE]
if (any(is.na(meta$animal))) stop("Some samples in count matrix are not in exp_design")

# Create infected from 16sCt: infected if 16sCt > 10 (column name needs backticks)
meta$infected <- ifelse(meta$`16sCt` > 10, "infected", "uninfected")
meta$infected <- factor(meta$infected, levels = c("uninfected", "infected"))

# Build DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = meta,
  design    = DESIGN
)

# Filter low-count genes (optional but often recommended)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq2: LRT full (~ sex + infected) vs reduced (~ sex) to get p-values for effect of infection
dds <- DESeq(dds, test = "LRT", reduced = DESIGN_REDUCED)

# ---- Main result: effect of infection (adjusted for sex); p-values from LRT (reduced vs full) ----
res_infected <- results(dds)
res_infected <- res_infected[order(res_infected$padj, res_infected$pvalue), ]
summary(res_infected)

# Primary results table: infected effect, sorted by p-value
res_inf_df <- as.data.frame(res_infected)
res_inf_df$gene_id <- rownames(res_inf_df)
res_inf_df$display_name <- if (length(loc_add_to_gene2) > 0) {
  ifelse(res_inf_df$gene_id %in% names(loc_add_to_gene2), loc_add_to_gene2[res_inf_df$gene_id], res_inf_df$gene_id)
} else {
  res_inf_df$gene_id
}
res_inf_df <- res_inf_df[, c("gene_id", "display_name", setdiff(colnames(res_inf_df), c("gene_id", "display_name")))]
write.csv(res_inf_df, file.path(out_dir, "results_infected_by_pvalue.csv"), row.names = FALSE)
message("Primary table: results_infected_by_pvalue.csv (infected vs uninfected, sorted by padj)")

# Normalized counts (with display_name)
nc <- counts(dds, normalized = TRUE)
nc_df <- as.data.frame(nc)
nc_df$gene_id <- rownames(nc_df)
nc_df$display_name <- if (length(loc_add_to_gene2) > 0) {
  ifelse(nc_df$gene_id %in% names(loc_add_to_gene2), loc_add_to_gene2[nc_df$gene_id], nc_df$gene_id)
} else {
  nc_df$gene_id
}
nc_df <- nc_df[, c("gene_id", "display_name", setdiff(colnames(nc_df), c("gene_id", "display_name")))]
write.csv(nc_df, file.path(out_dir, "normalized_counts.csv"), row.names = FALSE)

# ---- Save R objects for later (load and make custom figures) ----
vsd <- vst(dds, blind = FALSE)  # variance-stabilized counts for PCA/heatmaps
saveRDS(list(
  dds = dds,
  vsd = vsd,
  res_infected = res_infected,
  normalized_counts = nc,
  meta = meta,
  loc_add_to_gene2 = loc_add_to_gene2
), file.path(saved_dir, "deseq2_objects.rds"))
message("Saved: ", file.path(saved_dir, "deseq2_objects.rds"), " (load with readRDS(...) for custom figures)")

# ---- Default plots (saved in figures/) ----
# PCA: infected (main) and sex (covariate)
pdf(file.path(fig_dir, "PCA_infected_sex.pdf"), width = 6, height = 5)
print(plotPCA(vsd, intgroup = c("infected", "sex")) + ggplot2::theme_bw())
dev.off()
message("Saved: figures/PCA_infected_sex.pdf")

# MA plot: infected vs uninfected
pdf(file.path(fig_dir, "MA_infected.pdf"), width = 5, height = 5)
DESeq2::plotMA(res_infected, main = "Infected vs uninfected")
dev.off()
message("Saved: figures/MA_infected.pdf")

# Volcano: infected vs uninfected, label top 20 genes by padj
n_label <- 20
vd <- as.data.frame(res_infected)
vd$gene_id <- rownames(vd)
vd$display_name <- if (length(loc_add_to_gene2) > 0) {
  ifelse(vd$gene_id %in% names(loc_add_to_gene2), loc_add_to_gene2[vd$gene_id], vd$gene_id)
} else {
  vd$gene_id
}
vd$neglog10padj <- -log10(vd$padj + 1e-300)
vd$neglog10padj <- pmin(vd$neglog10padj, 50)  # cap for plot
vd$sig <- !is.na(vd$padj) & vd$padj < 0.05
top20 <- rownames(res_infected)[order(res_infected$padj)][seq_len(min(n_label, nrow(res_infected)))]
vd$label <- ifelse(vd$gene_id %in% top20, vd$display_name, "")

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- ggplot2::ggplot(vd, ggplot2::aes(x = log2FoldChange, y = neglog10padj, color = sig)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = Inf, box.padding = 0.5) +
    ggplot2::scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick3"), name = "padj < 0.05") +
    ggplot2::labs(x = "log2 fold change", y = "-log10(adj p)", title = "Infected vs uninfected (top 20 by padj)") +
    ggplot2::theme_bw()
} else {
  p <- ggplot2::ggplot(vd, ggplot2::aes(x = log2FoldChange, y = neglog10padj, color = sig)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3, hjust = -0.1, check_overlap = TRUE) +
    ggplot2::scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick3"), name = "padj < 0.05") +
    ggplot2::labs(x = "log2 fold change", y = "-log10(adj p)", title = "Infected vs uninfected (top 20 by padj)") +
    ggplot2::theme_bw()
}
pdf(file.path(fig_dir, "volcano_infected.pdf"), width = 7, height = 6)
print(p)
dev.off()
message("Saved: figures/volcano_infected.pdf (install ggrepel for nicer labels)")

# Heatmap: top 50 genes by padj (infected contrast), VST-transformed, row-scaled
# X-axis: two rows (Infected, Sex); columns ordered by treatment: infected M, infected F, uninfected M, uninfected F (no sample IDs).
top_n <- min(50, nrow(res_infected))
top_genes <- rownames(res_infected)[order(res_infected$padj)][seq_len(top_n)]
if (length(top_genes) >= 5) {
  mat <- assay(vsd)[top_genes, ]
  mat <- t(scale(t(mat)))
  # Order columns by treatment: infected males, infected females, uninfected males, uninfected females
  treatment_order <- factor(
    paste(meta[colnames(mat), "infected"], meta[colnames(mat), "sex"]),
    levels = c("infected Male", "infected Female", "uninfected Male", "uninfected Female")
  )
  ord <- order(treatment_order)
  mat <- mat[, ord, drop = FALSE]
  ann_col <- as.data.frame(meta[colnames(mat), c("infected", "sex"), drop = FALSE])
  colnames(ann_col) <- c("Infected", "Sex")
  ann_colors <- list(
    Infected = c(uninfected = "grey95", infected = "darkred"),
    Sex      = c(Female = "red", Male = "blue")
  )

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pdf(file.path(fig_dir, "heatmap_top50_infected.pdf"), width = 8, height = max(4, length(top_genes) * 0.12))
    pheatmap::pheatmap(mat,
                       scale = "none",
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       annotation_col = ann_col,
                       annotation_colors = ann_colors,
                       show_colnames = FALSE,
                       col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                       main = "Top DE genes (infected vs uninfected)",
                       fontsize_row = 6)
    dev.off()
    message("Saved: figures/heatmap_top50_infected.pdf (annotation: Infected, Sex; no sample IDs)")
  } else {
    pdf(file.path(fig_dir, "heatmap_top50_infected.pdf"), width = 8, height = max(4, length(top_genes) * 0.12))
    heatmap(mat, scale = "none", Rowv = NA, Colv = NA,
            col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
            main = "Top DE genes (infected vs uninfected)")
    dev.off()
    message("Saved: figures/heatmap_top50_infected.pdf (install pheatmap for column annotations: infected, sex)")
  }
}

message("Results, tables, and figures written to ", out_dir)
