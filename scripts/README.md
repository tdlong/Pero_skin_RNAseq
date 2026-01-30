# RNA-seq pipeline scripts (SLURM)

Pipeline: HISAT2 index → align (BAMs by animal ID) → featureCounts → DESeq2 in R.

## Order of operations

1. **Build HISAT2 index** (once):  
   `sbatch scripts/01_hisat2_build_index.slurm`

2. **Build sample list** (once, from project root):  
   `bash scripts/make_sample_list.sh`  
   Creates `samples_to_align.txt` (animal, R1, R2). Set `--array=1-N` in step 3 to match the number of data lines in that file (e.g. 21).

3. **Align** (array job):  
   `sbatch scripts/02_hisat2_align.slurm`  
   Output: `bams/<animal>.sorted.bam`

4. **Count** (single job):  
   `sbatch scripts/03_featurecounts.slurm`  
   Output: `counts/counts.txt`

5. **DESeq2** (in R):  
   Edit `scripts/DESeq2_analysis.R` (set `DESIGN` and paths), then run in R or `Rscript scripts/DESeq2_analysis.R`.  
   Results: `deseq2_results/results.csv`, `normalized_counts.csv`.

## Programs / modules

| Step | Module(s) |
|------|-----------|
| 01   | `hisat2/2.2.1` |
| 02   | `hisat2/2.2.1`, `samtools/1.15.1` |
| 03   | `subread/2.0.3` |
| 05   | `R/4.4.2` + DESeq2 (see below) |

**DESeq2 in R/4.4.2:** If DESeq2 is not already installed, load R and run once:
```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
```

## Paths

Scripts assume you submit from the **project root** (directory that contains `data/`, `ref/`, `exp_design.txt`). They use `PROJECT_ROOT=.` by default; set `PROJECT_ROOT=/path/to/project` when submitting from elsewhere.

## Array size (02)

After `make_sample_list.sh`, check:  
`wc -l samples_to_align.txt`  
Use that number minus 1 (header) as the array size, or leave `--array=1-22` and let extra tasks no-op.
