# Pero skin RNA-seq

HISAT2 → featureCounts → DESeq2 pipeline for skin RNA-seq. Scripts are in `scripts/` and are meant to run on a SLURM cluster.

**Sync workflow (Mac ↔ cluster):** Edit scripts on your Mac in Cursor, then commit and push so the cluster has the latest code. On the cluster, run `git pull` before running jobs.

```bash
# Mac (after editing)
git add -A && git commit -m "your message" && git push

# Cluster (before running pipeline)
git pull
```

**On the cluster:** If you already have a project folder with `data/` and `ref/`, clone the repo *into a different name* so you can move the repo contents into that folder (see below). See `scripts/README.md` for run order and module names.

## Quick start on cluster

**If you already have a folder with `data/`, `ref/`, and `exp_design.txt`:** clone into a temp dir, then move the repo contents up so `scripts/` sits next to `data/` and `ref/`:

```bash
cd /path/to/your/project   # folder that contains data/, ref/, exp_design.txt
git clone https://github.com/tdlong/Pero_skin_RNAseq.git repo
mv repo/.git repo/.gitignore repo/* .
rmdir repo
```

Then run the pipeline from that same folder:

```bash
bash scripts/make_sample_list.sh
sbatch scripts/01_hisat2_build_index.slurm
# after index finishes:
sbatch scripts/02_hisat2_align.slurm
sbatch scripts/03_featurecounts.slurm
# then run DESeq2 in R (see scripts/DESeq2_analysis.R)
```
