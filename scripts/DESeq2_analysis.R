# DESeq2 differential expression analysis
# Uses count matrix from featureCounts and exp_design.txt.
# Run in R (interactive or Rscript) after 03_featurecounts.slurm.
#
# R/4.4.2: If DESeq2 is not installed:
#   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("DESeq2")
#
# Experimental design: infected (from 16sCt > 10) + sex + infected:sex.
# Primary interest: infected vs uninfected; sex is covariate.
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

# Design: infected (16sCt > 10) + sex + infected:sex
DESIGN <- ~ infected + sex + infected:sex

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

# Run DESeq2
dds <- DESeq(dds)

# List results available (one per coefficient)
resultsNames(dds)

# Extract and write results for each coefficient (add display_name from bespoke names)
for (coef in resultsNames(dds)[-1]) {  # skip intercept
  res <- results(dds, name = coef)
  res <- res[order(res$padj), ]
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$display_name <- if (length(loc_add_to_gene2) > 0) {
    ifelse(res_df$gene_id %in% names(loc_add_to_gene2), loc_add_to_gene2[res_df$gene_id], res_df$gene_id)
  } else {
    res_df$gene_id
  }
  res_df <- res_df[, c("gene_id", "display_name", setdiff(colnames(res_df), c("gene_id", "display_name")))]
  fname <- paste0("results_", gsub("[^a-zA-Z0-9]", "_", coef), ".csv")
  write.csv(res_df, file.path(out_dir, fname), row.names = FALSE)
  message("Wrote ", fname)
}

# ---- Main contrast: infected vs uninfected (primary interest; sex is covariate) ----
coef_infected <- resultsNames(dds)[grepl("^infected_", resultsNames(dds))][1]
res_infected <- results(dds, name = coef_infected)
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
top_n <- min(50, nrow(res_infected))
top_genes <- rownames(res_infected)[order(res_infected$padj)][seq_len(top_n)]
if (length(top_genes) >= 5) {
  mat <- assay(vsd)[top_genes, ]
  mat <- t(scale(t(mat)))
  pdf(file.path(fig_dir, "heatmap_top50_infected.pdf"), width = 8, height = max(4, length(top_genes) * 0.12))
  heatmap(mat, scale = "none", Rowv = NA, Colv = NA,
          col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
          main = "Top DE genes (infected vs uninfected)")
  dev.off()
  message("Saved: figures/heatmap_top50_infected.pdf")
}

message("Results, tables, and figures written to ", out_dir)
