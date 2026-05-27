# DiRT v2 — Dicistronic Transcript Detection Pipeline

A standardized pipeline for detecting **dicistronic transcripts** (tRNA-mRNA co-transcripts)
from paired-end Illumina RNA-seq data of any plant species. Inputs are a reference
genome FASTA, a GFF3 annotation (e.g. from EnsemblPlants), and either raw FASTQ
files or pre-aligned BAMs.

The pipeline supports **two entry points**:

| Use this when… | Run this | Section |
|---|---|---|
| You have raw FASTQ files (or just SRA accessions) | `nextflow run main.nf` | [Path A](#path-a-from-fastqsra--nextflow) |
| You already have sorted+indexed BAMs (one per replicate) | `Rscript -e "rmarkdown::render(...)"` | [Path B](#path-b-from-pre-aligned-bams--rscript) |

Both paths produce the same final result tables.

## Repository layout

```
DiRTv2_pipeline/
├── main.nf                       Nextflow DSL2 entry point
├── nextflow.config               Profiles + per-process resources
├── modules/                      Per-step DSL2 modules
│   ├── adapter_removal.nf
│   ├── fastqc.nf
│   ├── hisat2_build.nf
│   ├── hisat2.nf
│   ├── samtools_merge.nf
│   ├── faidx.nf
│   ├── trnascan.nf
│   └── dirtv2_analysis.nf
├── bin/
│   ├── DiRTv2_manual_optimized.Rmd      Parameterized R Markdown — the actual analysis
│   ├── 01_RNA_mapping.sh         Optional standalone bash mapper (no Nextflow)
│   ├── 02_generate_tRNA_bed.sh   Helper: tRNAscan-SE + BED conversion (for Path B)
│   └── optional_sra_download.sh  Optional SRA → FASTQ helper
├── conf/
│   ├── environment.yml           Full conda env (all tools — for Path A)
│   ├── environment_minimal.yml   Smaller conda env (BAM-input only — for Path B)
│   ├── environment_relaxed.yml   Same as full but with unpinned versions
│   └── mapping_config.sh         Config for bin/01_RNA_mapping.sh
├── assets/
│   └── samples.csv               Template sample sheet
└── README.md                     You are here
```

---

## 1. Prerequisites (one-time setup)

Requirements common to both paths:

- **Linux or WSL2 on Windows** (tested on Ubuntu 22.04, WSL2)
- **Miniconda or Anaconda** ([install instructions](https://docs.conda.io/projects/miniconda/en/latest/))
- **~15 GB free disk space** for the full conda env
- **For Path A only:** Nextflow ≥ 23.04 ([install instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation))

### 1a. Speed up conda's solver (highly recommended)

The classic conda solver is notoriously slow on bioconda + R/Bioconductor environments
(can take 30+ minutes or fail). Switch to the fast `libmamba` solver one time:

```bash
conda install -n base conda-libmamba-solver -y
conda config --set solver libmamba
```

Optionally install `mamba` (an even faster drop-in replacement for `conda`):

```bash
conda install -n base -c conda-forge mamba -y
```

### 1b. Create the conda environment

Pick one of two YAMLs depending on your input data:

**Option 1 — Full env (for Path A: FASTQ input via Nextflow):**

```bash
cd path/to/DiRTv2_pipeline
mamba env create -f conf/environment.yml      # or: conda env create -f conf/environment.yml
conda activate dirtv2
```

This installs: SRA-tools, AdapterRemoval, FastQC, HISAT2, samtools, bedtools,
tRNAscan-SE, R 4.3 + tidyverse + Bioconductor (GenomicAlignments, GenomicFeatures,
rtracklayer, txdbmaker), openxlsx, rmarkdown. **~3 GB.**

**Option 2 — Minimal env (for Path B: BAM input via Rscript):**

```bash
cd path/to/DiRTv2_pipeline
mamba env create -f conf/environment_minimal.yml
conda activate dirtv2
```

This drops the trimming/alignment tools you don't need when inputs are already BAMs.
Keeps: samtools (for `faidx`), bedtools, tRNAscan-SE, R + Bioconductor + tidyverse.
**~1.5 GB.**

### 1c. Verify the env is correctly set up

```bash
echo "$CONDA_DEFAULT_ENV"                # must print: dirtv2
which samtools                           # must contain /envs/dirtv2/
which bedtools tRNAscan-SE Rscript
Rscript -e 'library(tidyverse); library(GenomicAlignments); library(rtracklayer); library(GenomicFeatures); library(txdbmaker); library(openxlsx); cat("ALL R PACKAGES OK\n")'
```

If anything is missing, stop and fix the install before proceeding.

---

## Path A: from FASTQ/SRA → Nextflow

Use this when you have paired-end Illumina FASTQ files (or only SRA accessions).
Nextflow handles trimming, alignment, BAM merging by biological replicate, tRNAscan-SE,
bedtools closest, and the DiRT v2 statistical analysis.

### A.1. (Optional) Download FASTQ from SRA accessions

Skip this step if you already have FASTQ files.

```bash
# Create a file with one SRA accession per line
cat > accessions.txt <<EOF
ERR3338762
ERR3338763
ERR3338764
ERR3338765
EOF

# Download + extract + gzip them, and generate a samples.csv stub
conda activate dirtv2
./bin/optional_sra_download.sh accessions.txt /path/to/fastq_dir samples.csv
```

After this finishes, **edit `samples.csv`** to assign each row to a biological
replicate (the `sample` column). Multiple rows with the same `sample` value get merged.

### A.2. Prepare your `samples.csv`

Format: one row per FASTQ pair. Multiple pairs that share the same `sample` value
will be merged into one BAM (biological replicate):

```csv
sample,fastq_1,fastq_2
RSL1,/abs/path/run1_1.fastq.gz,/abs/path/run1_2.fastq.gz
RSL1,/abs/path/run2_1.fastq.gz,/abs/path/run2_2.fastq.gz
RSL2,/abs/path/run3_1.fastq.gz,/abs/path/run3_2.fastq.gz
RSL3,/abs/path/run4_1.fastq.gz,/abs/path/run4_2.fastq.gz
```

**Use absolute paths** — Nextflow runs each process in a temporary working
directory, so relative paths break.

### A.3. Run the pipeline

```bash
conda activate dirtv2

nextflow run main.nf -profile standard,conda \
    --samplesheet   samples.csv \
    --genome        /abs/path/to/genome.fasta \
    --gff           /abs/path/to/annotation.gff3 \
    --outdir        results

# With a stricter FDR threshold (default is 0.05):
nextflow run main.nf -profile standard,conda \
    --samplesheet   samples.csv \
    --genome        /abs/path/to/genome.fasta \
    --gff           /abs/path/to/annotation.gff3 \
    --fdr_threshold 0.01 \
    --outdir        results
```

Resume after an interruption — Nextflow caches each step:

```bash
nextflow run main.nf -profile standard,conda ...same args... -resume
```

### A.4. (Optional) Standalone bash mapping without Nextflow

If you don't want Nextflow, the standalone `bin/01_RNA_mapping.sh` does the same
trim+align+merge steps using your `samples.csv`:

```bash
cp conf/mapping_config.sh my_config.sh
# Edit my_config.sh to set WORK_DIR, SAMPLESHEET, GENOME_FASTA
conda activate dirtv2
./bin/01_RNA_mapping.sh my_config.sh
```

This produces merged BAMs in `<WORK_DIR>/bam_merged/`. Then continue with **Path B**
below to run the R analysis on those BAMs.

---

## Path B: from pre-aligned BAMs → Rscript

Use this when you already have one sorted+indexed BAM per biological replicate
and just want to run the DiRT v2 statistical analysis. Faster, no Nextflow required.

### B.1. Build the genome FASTA index (~10 sec)

```bash
conda activate dirtv2
cd /path/to/your_data_folder
samtools faidx genome.fasta                     # produces genome.fasta.fai
```

### B.2. Generate `tRNA.bed` from the genome FASTA (~2–4 hours on a plant genome)

The DiRTv2 Rmd needs a BED file of tRNA loci produced by tRNAscan-SE. Two equivalent
ways to make it — both produce an identical 6-column BED file:

**Option 1 — One-command helper (recommended):**

```bash
conda activate dirtv2
./bin/02_generate_tRNA_bed.sh /path/to/your_data_folder/genome.fasta
# Output: /path/to/your_data_folder/{tRNA.txt, tRNA_stats.txt, tRNA.bed, trnascan.log}
```

The helper is idempotent — if `tRNA.bed` already exists in the output directory, it
skips the slow tRNAscan-SE step.

**Option 2 — Manual two-step:**

```bash
cd /path/to/your_data_folder
nohup tRNAscan-SE -E -o tRNA.txt -m tRNA_stats.txt genome.fasta > trnascan.log 2>&1 &

# Monitor:  tail -f trnascan.log
# When finished, convert tabular output → BED (skip the 3-line header):
awk 'NR > 3 {print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $9 "\t" $6}' tRNA.txt > tRNA.bed
wc -l tRNA.bed                                  # sanity: hundreds-to-thousands of tRNAs
```

**Note for both options:** about half of all detected tRNAs will have `start > end`
in the BED file — that's normal tRNAscan-SE output for minus-strand tRNAs, and
the DiRTv2 Rmd normalizes the coordinates internally before use.

> Path A users do **not** need this step — `modules/trnascan.nf` runs the same
> tRNAscan-SE + awk pipeline automatically inside Nextflow.

### B.3. Sanity-check that chromosome names agree across FASTA, GFF3, BAM, tRNA.bed

This is the #1 silent failure mode. **Run all three commands and confirm matching naming:**

```bash
cd /path/to/your_data_folder
cut -f1 genome.fasta.fai                                            | sort -u
awk '$1 !~ /^#/ {print $1}' annotation.gff3                         | sort -u
samtools view -H BAM_FILE_1.sorted.bam | grep '^@SQ' | awk '{print $2}' | sed 's/^SN://'
cut -f1 tRNA.bed                                                    | sort -u
```

All four lists must use the same naming convention (e.g. all `chr01`, `chr02`, …
or all `1`, `2`, …). If they differ, fix the source files before running.

### B.4. Set your working directory

```bash
mkdir -p /path/to/dirtv2_run
cd       /path/to/dirtv2_run
cp /path/to/DiRTv2_pipeline/bin/DiRTv2_manual_optimized.Rmd .
```

No need to copy or symlink the BAMs — `bam_dir` below will point at them directly.

### B.5. Run the optimized R Markdown analysis

⚠ **Use absolute paths for every input.** Relative paths inside `rmarkdown::render`
can resolve differently than your shell's `cd`, leading to "BAM files not found" errors.

```bash
cd /path/to/dirtv2_run

Rscript -e "rmarkdown::render('DiRTv2_manual_optimized.Rmd', \
    params = list( \
        gff_file       = '/abs/path/to/annotation.gff3', \
        genome_fai     = '/abs/path/to/genome.fasta.fai', \
        tRNA_bed       = '/abs/path/to/tRNA.bed', \
        bam_dir        = '/abs/path/to/folder_containing_bams', \
        sample_names   = c('RSL1','RSL2','RSL3','RSL4'), \
        min_count      = 1, \
        fdr_threshold  = 0.05, \
        out_dir        = 'results', \
    ))" 2>&1 | tee dirtv2_render.log
```

#### Critical naming convention

`sample_names` MUST exactly match the BAM filename stem (everything before `.sorted.bam`):

| sample_names entry | Expected BAM filename in `bam_dir/` |
|---|---|
| `RSL1` | `RSL1.sorted.bam` (+ `RSL1.sorted.bam.bai`) |
| `LeafRep_A` | `LeafRep_A.sorted.bam` (+ `LeafRep_A.sorted.bam.bai`) |

### B.6. (Worked example with the test data referenced in this repo)

```bash
conda activate dirtv2
mkdir -p /mnt/f/PhD_Research/DiRTv2_optimized_run
cd       /mnt/f/PhD_Research/DiRTv2_optimized_run
cp /mnt/f/PhD_Research/.../DiRTv2_pipeline/bin/DiRTv2_manual_optimized.Rmd .

Rscript -e "rmarkdown::render('DiRTv2_manual_optimized.Rmd', \
    params = list( \
        gff_file       = '/mnt/f/PhD_Research/Thesis_writing/Manuscripts/Bio-Protocol/Data/PN40024_5.1_on_T2T_ref_with_names.gff3', \
        genome_fai     = '/mnt/f/PhD_Research/Thesis_writing/Manuscripts/Bio-Protocol/Data/T2T_ref.fasta.fai', \
        tRNA_bed       = '/mnt/f/PhD_Research/Thesis_writing/Manuscripts/Bio-Protocol/Data/tRNA.bed', \
        bam_dir        = '/mnt/f/PhD_Research/Thesis_writing/Manuscripts/Bio-Protocol/Data', \
        sample_names   = c('RSL1','RSL2','RSL3','RSL4'), \
        min_count      = 1, \
        fdr_threshold  = 0.05, \
        out_dir        = 'results', \
    ))" 2>&1 | tee dirtv2_render.log
```

---

## 2. Output files

After a successful run you'll find, under `results/` (Path B) or `<outdir>/dirtv2/results/` (Path A):

| File | What it is |
|---|---|
| **`Final_Result_Manual_optimized.xlsx`** | **The final result table.** One row per detected dicistronic-transcript candidate. Columns include: `id, chr, DT.start, DT.end, gene.id, tRNA.id, tRNA.start, tRNA.end, tRNA.strand, tRNA.upstream, tRNAFirst, interval.length, num.Intron, combined.continu.cov, <sample>_nZero, gene.chr, gene.start, gene.end, gene.strand, CDS.start, CDS.end, tRNA_Overlaps_Gene, tRNA_In_CDS, tRNA_Overlaps_CDS, Location_Type` |
| `DiRTv2_manual_optimized.html` | Full rendered analysis report — includes every plot and printed sanity check |

The working directory also retains the bedtools intermediates (`Vitis_T2T_*.txt`, `*_closest_*.bed`, `tRNA_CDS_Combination.txt`, etc.) — useful for debugging but not required as final outputs.

For Path A only, you also get, alongside the analysis outputs above:

- `fastq_trimmed/` — AdapterRemoval output
- `fastqc/` — FastQC HTML reports
- `bam_per_run/` — One BAM per FASTQ pair
- `bam_merged/` — One BAM per biological replicate
- `trnascan/tRNA.bed` — tRNAscan-SE BED (auto-generated from the genome FASTA)
- `pipeline_info/{execution_report.html,timeline.html,trace.txt,pipeline_dag.svg}`

### Location_Type values

The final `Location_Type` column classifies where each tRNA sits relative to its partner gene:

| Value | Meaning |
|---|---|
| `Intergenic` | tRNA does not overlap the partner gene at all (cleanest case) |
| `UTR` | tRNA overlaps the gene but not its CDS — likely 5'/3' UTR |
| `CDS-Internal` | tRNA is fully inside the CDS (rare; flag for inspection) |
| `CDS-Edge_Conflict` | tRNA partially overlaps the CDS edge (rare; flag for inspection) |

---

## 3. Tunable parameters

All parameters are settable from the CLI (Path A) or the `params = list(...)` block (Path B):

| Parameter | Default | Effect |
|---|---|---|
| `samplesheet` / sample sheet rows | — | CSV listing FASTQ pairs + replicate grouping (Path A only) |
| `genome` / `--genome` | — | Reference FASTA (Path A only — Path B uses `genome_fai`) |
| `gff` / `gff_file` | — | GFF3 annotation (any plant species, e.g. EnsemblPlants) |
| `genome_fai` | — | FASTA index (Path B; auto-built in Path A) |
| `tRNA_bed` | — | tRNAscan-SE output in BED format (Path B; auto-built in Path A) |
| `bam_dir` | — | Absolute path to folder containing the BAMs (Path B) |
| `sample_names` | — | Vector of replicate names matching BAM filename stems (Path B) |
| `hisat2_index` | (built from genome) | Optional prebuilt HISAT2 index prefix (Path A) |
| `min_count` | 1 | Minimum read count in **every** sample to call a tRNA/gene "expressed" |
| `fdr_threshold` | 0.05 | FDR cutoff for paired t-tests (intergenic vs intron, intron vs intron) |
| `min_length` | 100 | AdapterRemoval `--minlength` (Path A only) |
| `cpus` | 8 | CPUs per process |
| `memory` | 32 GB | RAM per process |
| `outdir` | `results` | Output directory root |

---

## 4. Execution profiles (Path A only)

Combine with commas in `-profile`:

| Profile | Use case |
|---|---|
| `standard` | Local execution (default) |
| `conda` | Activate the `dirtv2` conda env for every process |
| `docker` | Run each process in a Docker container |
| `singularity` | Run each process in Singularity (HPC-friendly) |
| `slurm` | Submit each process as a SLURM job |
| `test` | Pre-populated paths from `assets/samples.csv` for a smoke test |

Example: `-profile slurm,singularity` runs on SLURM with Singularity containers.

---

## 5. Pinned tool versions (`conf/environment.yml`)

- sra-tools 3.0.10 (optional — only for SRA download)
- AdapterRemoval 2.3.3
- FastQC 0.12.1
- HISAT2 2.2.1
- samtools 1.19
- bedtools 2.31.1
- tRNAscan-SE 2.0.12
- R 4.3 + tidyverse, magrittr, data.table, openxlsx, reshape2, rmarkdown, knitr
- Bioconductor: GenomicAlignments, GenomicFeatures, rtracklayer, txdbmaker

---

## 6. Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `EnvironmentNameNotFound: dirtv2` | env wasn't created (the install errored silently) | Re-run `mamba env create -f conf/environment.yml` and watch for solver errors |
| `Solving environment: working...` for > 30 min | Classic conda solver | Switch to libmamba: `conda config --set solver libmamba` |
| `Could not match all params$sample_names to BAM files` (Path B) | `bam_dir` is wrong or BAM filenames don't match `sample_names` | Use an **absolute** path; rename BAMs to match `<sample>.sorted.bam` exactly |
| Error from `intergenic-info` chunk: "names of metadata columns cannot be one of strand, …" | Old version of the Rmd before the `make_gr()` patch | Re-copy `bin/DiRTv2_manual_optimized.Rmd` from this repo |
| `unable to find an inherited method for function 'first' for signature 'x = "factor"'` | Old version of the Rmd before the namespace-qualification patch | Re-copy `bin/DiRTv2_manual_optimized.Rmd` |
| `Error in readGAlignments`: zero reads counted | BAM chromosome names don't match FASTA/GFF3 (e.g. `chr01` vs `1`) | Re-do Step B.3 chromosome-name check and reconcile |
| Out-of-RAM at `count-tRNA` / `count-genes` / `class-multi` | BAMs are larger than available RAM | Reduce concurrent samples; increase `--memory` |
| Render finishes with no error but `results/` is empty | Working dir mismatch — check `pwd` matches where you `cp`'d the .Rmd | Run Rscript from the same dir as the .Rmd |

---

## 7. Citation

This pipeline accompanies the Bio-Protocol manuscript:

> Zheng, F. et al. (2026). *DiRT v2: A Standardized Pipeline for Dicistronic Transcript
> Detection from Plant RNA-Seq Data.* Bio-Protocol.

Tool citations are included in the rendered HTML report's session-info section. If
you use the Nextflow path, please also cite Nextflow (Di Tommaso et al., 2017).

## 8. License

[Add your preferred license here, e.g. MIT or BSD-3-Clause.]
