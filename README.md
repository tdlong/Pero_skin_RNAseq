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

**Logs:** SLURM output goes to `logs/index/`, `logs/align/`, `logs/counts/`, and `logs/deseq2/` by step; filenames include the job name (e.g. `hisat2_index_47947.out`) so you can tell what each log is. To thin the pile: keep only recent, or archive by date, e.g. `mkdir -p archive/logs && mv logs/* archive/logs/$(date +%Y-%m-%d)/` then rerun; or delete old ones: `find logs -name "*.out" -mtime +30 -delete`.

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
sbatch scripts/04_deseq2.slurm
# or run DESeq2 in R: Rscript scripts/DESeq2_analysis.R
```

**featureCounts:** The step 3 logs report “Successfully assigned alignments” as a fraction of total alignments. Only about 20–25% of read pairs are assigned to features (exons/genes) in this pipeline; that is what featureCounts reports and is expected (many alignments fall in introns, intergenic, or multimapping and are not counted toward gene expression).

---

## Analysis walkthrough (conceptual)

1. **Sample list** – `make_sample_list.sh` matches FASTQ filenames (i7/i5) to `exp_design.txt` and writes `samples_to_align.txt` (animal, R1, R2).

2. **Index** – HISAT2 builds a genome index from the reference (no splice/exon hints; plain index for speed and low RAM). Output: `ref/PerLeu_trans*`.

3. **Align** – HISAT2 aligns paired-end reads; BAMs are sorted and indexed and named by animal ID. Output: `bams/<animal>.sorted.bam`.

4. **Counts** – featureCounts assigns reads to genes using the GTF. Output: `counts/counts.txt` (gene × sample matrix).

5. **Differential expression** – DESeq2 full model: **sex + infected** (sex as covariate). We want the p-value for the **effect of infection**; that comes from a likelihood-ratio test: full model (~ sex + infected) vs reduced model (~ sex). So we fit sex + infected and test “does adding infection improve the model?” — that’s the infection effect adjusted for sex. Primary table: infected vs uninfected (LRT padj, log2FC). Outputs live under `deseq2_results/`.

**Model and testing (documented):**
- **Full model:** `~ sex + infected`. Sex and infection are the two variables of interest; sex is fitted as a covariate so the infection effect is adjusted for it.
- **Infection p-values:** From LRT: compare full (~ sex + infected) to reduced (~ sex). The LRT p-value is "effect of infection" (one per gene). No Wald on a single coefficient; we use the reduced-model comparison so the test is exactly "does infection matter, given sex?"
- **Heatmap:** Top 50 DE genes (infected contrast). Columns = samples **grouped by infection, then by sex within infection** (uninfected first, then infected; within each, ordered by sex). Column annotation bars show infected and sex. Requires R package `pheatmap`.

### Key outputs (so you know what each file is)

**Tables (in `deseq2_results/`)**

| File | What it is |
|------|------------|
| `results_infected_by_pvalue.csv` | **Primary table.** Infected vs uninfected, one row per gene, sorted by adjusted p-value. Columns include `gene_id`, `display_name` (bespoke name when in `ref/cds_Pleuc_names_v6.txt`), log2FC, padj, etc. |
| `normalized_counts.csv` | Normalized counts (gene × sample) plus `gene_id` and `display_name`, for plotting. |
| (other CSVs) | Primary output is `results_infected_by_pvalue.csv` (infection effect from LRT). |

**Saved R objects (in `deseq2_results/saved/`)**

| File | What it is |
|------|------------|
| `deseq2_objects.rds` | Single RDS containing `dds`, `vsd`, `res_infected`, `normalized_counts`, `meta`, `loc_add_to_gene2`. Load with `readRDS(...)` to make custom figures without re-running DESeq2. |

**Figures (in `deseq2_results/figures/`)**

| File | What it is |
|------|------------|
| `PCA_infected_sex.pdf` | PCA of samples (VST-transformed counts), colored by **infected**, shaped by **sex**. Quick check that infection and sex separate samples. |
| `MA_infected.pdf` | MA plot for infected vs uninfected: log2 fold change vs mean expression; DE genes stand out. |
| `volcano_infected.pdf` | Volcano plot: log2FC vs -log10(adj p). Points colored by padj &lt; 0.05; **top 20 genes by padj labeled** (display_name). Install `ggrepel` for non-overlapping labels. |
| `heatmap_top50_infected.pdf` | Heatmap of the top 50 genes by padj (infected contrast). Rows = genes (row-scaled), columns = samples **ordered by infection then sex** (uninfected then infected; within each group by sex). Column annotation bars show **infected** and **sex**. VST-transformed counts. Requires `pheatmap`. |
