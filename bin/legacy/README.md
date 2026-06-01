# Legacy Rmd

This folder contains the **original, unoptimized DiRT v2 R Markdown** as it was
written for interactive use in RStudio. It is kept here for reference only —
the supported, parameterized pipeline is `bin/DiRTv2_manual_optimized.Rmd`
(invoked by `main.nf` for the Nextflow workflow).

## Contents

- `DiRTv2_vitis_fromGff3_annotation.Rmd` — the original, hand-written DiRT v2
  Rmd. Hard-coded file paths, hard-coded sample names, designed for manual,
  chunk-by-chunk execution in RStudio. Useful if you want to inspect a single
  step of the pipeline (for example, the tRNA-gene combination construction
  from the GFF) or to compare the modernised pipeline against its origin.

## Not for production use

This file is **not** parameterised, **not** invoked by the Nextflow workflow,
and **not** maintained. Bug reports and PRs against this file will not be
addressed — the active codebase is the parameterised Rmd one folder up.

If you want a runnable workflow, use:

- **Path A** (FASTQ → result): `nextflow run main.nf …` (see README §Path A)
- **Path B** (BAM → result):   `Rscript - <<'RSCRIPT' … RSCRIPT` (see README §Path B)
