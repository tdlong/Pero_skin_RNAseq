#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2:00:00
#SBATCH --output=logs/counts/%x_%j.out
#SBATCH --error=logs/counts/%x_%j.err

# Generate gene-level count matrix from sorted BAMs using featureCounts (Subread).
# Run after 02_hisat2_align.slurm. Output: counts/counts.txt (and .summary).
# Requires: subread (featureCounts)
# module load subread/2.0.x  # or your cluster's module name

set -e
mkdir -p logs/counts counts

PROJECT_ROOT="${PROJECT_ROOT:-.}"
BAM_DIR="${PROJECT_ROOT}/bams"
GTF="${PROJECT_ROOT}/ref/GCF_004664715.2_UCI_PerLeu_2.1_genomic.gtf"
OUT_PREFIX="${PROJECT_ROOT}/counts/counts"

module load subread/2.0.3

# All sorted BAMs (one per animal)
BAMS=("${BAM_DIR}"/*.sorted.bam)
if [ ! -f "${BAMS[0]}" ]; then
  echo "No BAMs found in ${BAM_DIR}" >&2
  exit 1
fi

featureCounts -T 4 -t exon -g gene_id \
  -a "${GTF}" \
  -o "${OUT_PREFIX}.txt" \
  --countReadPairs \
  "${BAMS[@]}"

echo "Count matrix: ${OUT_PREFIX}.txt"
echo "Summary: ${OUT_PREFIX}.txt.summary"
