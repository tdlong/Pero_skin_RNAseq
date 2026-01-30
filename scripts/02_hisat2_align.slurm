#!/bin/bash
#SBATCH --job-name=hisat2_align
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=24:00:00
#SBATCH --array=1-22
#SBATCH --output=logs/align/%x_%A_%a.out
#SBATCH --error=logs/align/%x_%A_%a.err

# Align paired-end RNA-seq with HISAT2; output sorted BAMs named by animal ID.
# Prereqs: run 01_hisat2_build_index.slurm once, then run make_sample_list.sh to create samples_to_align.txt.
# Set NSAMPLES in the script or leave as 22 and ensure samples_to_align.txt has that many lines (excluding header).
# Requires: hisat2, samtools
# module load hisat2/2.2.1
# module load samtools/1.x

set -e
mkdir -p logs/align bams

# Paths (run from project root)
PROJECT_ROOT="${PROJECT_ROOT:-.}"
SAMPLE_LIST="${PROJECT_ROOT}/samples_to_align.txt"
INDEX_BASE="${PROJECT_ROOT}/ref/PerLeu_trans"
BAM_DIR="${PROJECT_ROOT}/bams"

module load hisat2/2.2.1
module load samtools/1.15.1

# Array index: line number in sample list (skip header)
line=$((SLURM_ARRAY_TASK_ID + 1))
animal=$(awk -v n="$line" 'NR==n {print $1}' "$SAMPLE_LIST")
R1=$(awk -v n="$line" 'NR==n {print $2}' "$SAMPLE_LIST")
R2=$(awk -v n="$line" 'NR==n {print $3}' "$SAMPLE_LIST")

if [ -z "$animal" ] || [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
  echo "Skipping task $SLURM_ARRAY_TASK_ID (no sample or missing files)" >&2
  exit 0
fi

out_bam="${BAM_DIR}/${animal}.sorted.bam"

hisat2 -p 2 -x "${INDEX_BASE}" -1 "$R1" -2 "$R2" | \
  samtools view -bS - | \
  samtools sort -@ 1 -o "${out_bam}" -
samtools index "${out_bam}"

echo "Done: ${out_bam}"
