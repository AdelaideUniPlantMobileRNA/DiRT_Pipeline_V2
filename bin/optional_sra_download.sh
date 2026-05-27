#!/usr/bin/env bash
###############################################################################
# OPTIONAL helper — only needed when starting from SRA accessions instead of
# local FASTQ files. Downloads + extracts paired-end FASTQ files for a list
# of accessions, then writes a samples.csv stub you can plug into the
# Nextflow pipeline.
#
# Usage:
#   ./optional_sra_download.sh <accessions.txt> <output_dir> [<sample_csv_out>]
#
# accessions.txt: one SRA accession per line (e.g. ERR3338762)
#
# After this completes, edit the generated samples.csv to assign each FASTQ
# pair to a biological replicate (`sample` column), then run:
#   nextflow run main.nf --samplesheet samples.csv --genome ... --gff ...
###############################################################################
set -euo pipefail

ACC_FILE="${1:?accessions file required}"
OUT_DIR="${2:?output directory required}"
SAMPLE_CSV="${3:-${OUT_DIR}/samples.csv}"
THREADS="${THREADS:-8}"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

echo "[$(date '+%F %T')] prefetch"
prefetch --option-file "$ACC_FILE"

echo "[$(date '+%F %T')] fasterq-dump"
while read -r acc; do
    [[ -z "$acc" || "$acc" =~ ^# ]] && continue
    [[ -f "${acc}_1.fastq" && -f "${acc}_2.fastq" ]] && continue
    fasterq-dump --split-3 --threads "$THREADS" "$acc"
done < "$ACC_FILE"

# Optionally gzip to save space
echo "[$(date '+%F %T')] gzipping FASTQ"
for f in *_1.fastq *_2.fastq; do
    [[ -f "$f" && ! -f "${f}.gz" ]] && gzip "$f"
done

# Generate samples.csv stub — defaults sample == accession; user should edit
# this file to group accessions into biological replicates.
{
    echo "sample,fastq_1,fastq_2"
    while read -r acc; do
        [[ -z "$acc" || "$acc" =~ ^# ]] && continue
        echo "${acc},$(pwd)/${acc}_1.fastq.gz,$(pwd)/${acc}_2.fastq.gz"
    done < "$ACC_FILE"
} > "$SAMPLE_CSV"

echo "Done. Wrote samples.csv stub to: $SAMPLE_CSV"
echo "Edit the 'sample' column to group runs into biological replicates."
