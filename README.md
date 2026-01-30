# Pero skin RNA-seq

HISAT2 → featureCounts → DESeq2 pipeline for skin RNA-seq. Scripts are in `scripts/` and are meant to run on a SLURM cluster.

**On the cluster:** Clone this repo, then put your `data/` (fastq) and `ref/` (genome + GTF) in the project directory. See `scripts/README.md` for run order and module names.

## Quick start on cluster

```bash
git clone https://github.com/YOUR_USERNAME/Pero_skin_RNAseq.git
cd Pero_skin_RNAseq
# Add data/, ref/, and exp_design.txt (or symlink)
bash scripts/make_sample_list.sh
sbatch scripts/01_hisat2_build_index.slurm
# after index finishes:
sbatch scripts/02_hisat2_align.slurm
sbatch scripts/03_featurecounts.slurm
# then run DESeq2 in R (see scripts/DESeq2_analysis.R)
```
