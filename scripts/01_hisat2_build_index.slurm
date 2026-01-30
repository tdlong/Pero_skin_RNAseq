#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5G
#SBATCH --time=2:00:00
#SBATCH --output=logs/index/%x_%j.out
#SBATCH --error=logs/index/%x_%j.err

# Build HISAT2 index (no splice/exon hints: ~8G RAM, fast; mammalian with --ss/--exon needs ~200G and days).
# Requires: hisat2

set -e
mkdir -p logs/index ref

REF_DIR="ref"
REF_FASTA="${REF_DIR}/GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna"
INDEX_BASE="${REF_DIR}/PerLeu_trans"

module load hisat2/2.2.1

hisat2-build -p 4 "${REF_FASTA}" "${INDEX_BASE}"

echo "Index build complete: ${INDEX_BASE}"
