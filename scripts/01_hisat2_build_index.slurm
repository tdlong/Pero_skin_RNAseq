#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH -A tdlong_lab
#SBATCH -p highmem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=4:00:00
#SBATCH --output=logs/hisat2_build_%j.out
#SBATCH --error=logs/hisat2_build_%j.err

# Build HISAT2 index with splice sites and exons from GTF (one-time job).
# Requires: hisat2, python (for hisat2_extract_splice_sites.py and hisat2_extract_exons.py)
# module load hisat2/2.2.1  # or whatever your cluster provides

set -e
mkdir -p logs ref

REF_DIR="ref"
REF_FASTA="${REF_DIR}/GCF_004664715.2_UCI_PerLeu_2.1_genomic.fna"
GTF="${REF_DIR}/GCF_004664715.2_UCI_PerLeu_2.1_genomic.gtf"
INDEX_BASE="${REF_DIR}/PerLeu_trans"

module load hisat2/2.2.1

# Extract splice sites and exons (HISAT2 scripts are in PATH with the module)
hisat2_extract_splice_sites.py "${GTF}" > "${REF_DIR}/PerLeu.ss"
hisat2_extract_exons.py "${GTF}"     > "${REF_DIR}/PerLeu.exon"

# Build index (--ss and --exon improve RNA-seq alignment)
hisat2-build -p 4 --exon "${REF_DIR}/PerLeu.exon" --ss "${REF_DIR}/PerLeu.ss" \
  "${REF_FASTA}" "${INDEX_BASE}"

echo "Index build complete: ${INDEX_BASE}"
