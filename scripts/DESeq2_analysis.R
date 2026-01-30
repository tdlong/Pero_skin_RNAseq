# DESeq2 differential expression analysis
# Uses count matrix from featureCounts and exp_design.txt.
# Run in R (interactive or Rscript) after 03_featurecounts.slurm.
#
# R/4.4.2: If DESeq2 is not installed:
#   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("DESeq2")
#
# Before running:
#   - Edit DESIGN below to your comparison of interest (e.g. ~ sex, ~ mating, ~ exp).
#   - Ensure column names in the count matrix match the "animal" column in exp_design.

library(DESeq2)

# Paths (adjust if not running from project root)
project_root <- "."
count_file  <- file.path(project_root, "counts/counts.txt")
meta_file   <- file.path(project_root, "exp_design.txt")
out_dir     <- file.path(project_root, "deseq2_results")
dir.create(out_dir, showWarnings = FALSE)

# Design formula: change to your comparison of interest
# Examples: ~ sex  |  ~ mating  |  ~ exp  |  ~ sex + mating
DESIGN <- ~ sex

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

# List results available (one per coefficient when design has factors)
resultsNames(dds)

# Extract results for the first (or main) comparison
# Example: for design ~ sex, first coefficient is typically "sex_Male_vs_Female"
res <- results(dds, name = resultsNames(dds)[2])  # or specify: name = "sex_Male_vs_Female"
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), file.path(out_dir, "results.csv"), row.names = TRUE)

# Summary
summary(res)

# Optional: normalized counts for plotting
nc <- counts(dds, normalized = TRUE)
write.csv(nc, file.path(out_dir, "normalized_counts.csv"), row.names = TRUE)

message("Results written to ", out_dir)
