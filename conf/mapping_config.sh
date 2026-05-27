###############################################################################
# Config for bin/01_RNA_mapping.sh
# Generic mapping helper — works with any plant genome and any Illumina
# paired-end RNA-seq data. Edit values to match your environment, then run:
#   ./bin/01_RNA_mapping.sh conf/mapping_config.sh
###############################################################################

# Where everything lands (logs, trimmed FASTQ, BAMs, FastQC reports)
WORK_DIR="/mnt/f/PhD_Research/DiRTv2_run"

# CSV: sample,fastq_1,fastq_2  (see assets/samples.csv)
SAMPLESHEET="${WORK_DIR}/samples.csv"

# Reference genome FASTA (e.g. Vitis_vinifera.PN40024.v4.dna.toplevel.fa from
# EnsemblPlants). The matching GFF3 is consumed by the R analysis, not here.
GENOME_FASTA="/path/to/genome.fasta"

# Optional: prebuilt HISAT2 index prefix (no .ht2 suffix).
# Leave empty to build from GENOME_FASTA automatically.
HISAT2_INDEX=""

# Thread / memory tuning
THREADS=8
MIN_LENGTH=100
MEM_PER_THREAD="1G"
