#!/bin/bash
#SBATCH --job-name=deseq2
#SBATCH -A tdlong_lab
#SBATCH -p highmem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2:00:00
#SBATCH --output=logs/deseq2/%x_%j.out
#SBATCH --error=logs/deseq2/%x_%j.err

# Run DESeq2 analysis (infected + sex + infected:sex). Run after 03_featurecounts.slurm.
# Requires: R/4.4.2 with DESeq2 installed.
# Output: deseq2_results/ (tables, saved R objects, figures).

set -e
mkdir -p logs/deseq2

PROJECT_ROOT="${PROJECT_ROOT:-.}"
cd "${PROJECT_ROOT}"

module load R/4.4.2

Rscript scripts/DESeq2_analysis.R

echo "DESeq2 done. Check deseq2_results/"
