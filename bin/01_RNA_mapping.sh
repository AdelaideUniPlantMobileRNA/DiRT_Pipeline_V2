#!/usr/bin/env bash
###############################################################################
# DiRT v2 — Standalone mapping helper (no workflow manager required)
#
# Generic version: works with any paired-end Illumina FASTQ files and any
# reference genome / GFF3 (e.g. from EnsemblPlants). Use this if you want to
# run the upstream trim+map+merge steps without Nextflow.
#
# Inputs are described entirely in samples.csv:
#   sample,fastq_1,fastq_2
#   ERR1,/abs/path/to/RUN1_1.fq.gz,/abs/path/to/RUN1_2.fq.gz
#   ERR1,/abs/path/to/RUN2_1.fq.gz,/abs/path/to/RUN2_2.fq.gz
#   ERR2,/abs/path/to/RUN3_1.fq.gz,/abs/path/to/RUN3_2.fq.gz
#
# Multiple rows that share the same `sample` value are merged into one BAM
# (biological replicate). The set of unique `sample` values becomes the
# column names consumed by the DiRT v2 R analysis.
#
# Usage:
#   ./01_RNA_mapping.sh <config.sh>
###############################################################################
set -euo pipefail
IFS=$'\n\t'

CONFIG="${1:-conf/mapping_config.sh}"
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONFIG"

: "${WORK_DIR:?WORK_DIR not set}"
: "${SAMPLESHEET:?SAMPLESHEET not set}"
: "${GENOME_FASTA:?GENOME_FASTA not set}"
: "${THREADS:=8}"
: "${MIN_LENGTH:=100}"
: "${MEM_PER_THREAD:=1G}"
: "${HISAT2_INDEX:=}"   # optional; built from GENOME_FASTA if empty

mkdir -p "$WORK_DIR"/{logs,fastqc_out,fastq_trimmed,bam_per_run,bam_merged,hisat2_index}
cd "$WORK_DIR"

LOG="logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1
log()  { echo "[$(date '+%F %T')] $*"; }
fail() { echo "[$(date '+%F %T')] ERROR: $*" >&2; exit 1; }

# ----------------------------------------------------------------------------
# 1. Optionally build HISAT2 index from GENOME_FASTA
# ----------------------------------------------------------------------------
if [[ -z "$HISAT2_INDEX" ]]; then
    HISAT2_INDEX="${WORK_DIR}/hisat2_index/$(basename "${GENOME_FASTA%.*}")_index"
    if ! ls "${HISAT2_INDEX}".*.ht2 >/dev/null 2>&1; then
        log "Building HISAT2 index → ${HISAT2_INDEX}"
        hisat2-build -p "$THREADS" "$GENOME_FASTA" "$HISAT2_INDEX"
    else
        log "Re-using HISAT2 index at ${HISAT2_INDEX}"
    fi
fi

# ----------------------------------------------------------------------------
# 2. Trim + align every (sample, run) pair listed in SAMPLESHEET
# ----------------------------------------------------------------------------
log "Reading samplesheet: $SAMPLESHEET"
# CSV has header sample,fastq_1,fastq_2
declare -A SAMPLE_BAMS=()    # sample -> space-separated BAM paths

tail -n +2 "$SAMPLESHEET" | while IFS=, read -r sample r1 r2; do
    [[ -z "$sample" ]] && continue
    run_id="${sample}_$(basename "${r1%%.f*q*}")"
    trimmed_r1="fastq_trimmed/${run_id}_R1.trimmed.fq.gz"
    trimmed_r2="fastq_trimmed/${run_id}_R2.trimmed.fq.gz"
    bam_out="bam_per_run/${run_id}.sorted.bam"

    # AdapterRemoval
    if [[ ! -f "$trimmed_r1" || ! -f "$trimmed_r2" ]]; then
        log "Trim:  $run_id"
        AdapterRemoval \
            --file1 "$r1" --file2 "$r2" \
            --minlength "$MIN_LENGTH" --trimns --trimqualities --gzip \
            --threads "$THREADS" \
            --output1   "$trimmed_r1" \
            --output2   "$trimmed_r2" \
            --singleton "fastq_trimmed/${run_id}_singleton.trimmed.fq.gz"
    fi

    # HISAT2 → sorted, indexed BAM
    if [[ ! -f "${bam_out}.bai" ]]; then
        log "Align: $run_id"
        hisat2 -p "$THREADS" -x "$HISAT2_INDEX" \
            -1 "$trimmed_r1" -2 "$trimmed_r2" \
            --rna-strandness RF --dta --mp 4,2 --score-min L,0,-0.4 \
            --summary-file "logs/${run_id}_hisat2_summary.txt" \
          | samtools view -@ "$THREADS" -b -F 4 - \
          | samtools sort -@ "$THREADS" -m "$MEM_PER_THREAD" -o "$bam_out" -
        samtools index -@ "$THREADS" "$bam_out"
    fi
    # Persist the (sample → bam) mapping via temp file (subshell limitation)
    echo "${sample}|${bam_out}" >> "logs/.sample_bam_map.tsv"
done

# ----------------------------------------------------------------------------
# 3. FastQC on all trimmed reads
# ----------------------------------------------------------------------------
log "FastQC"
fastqc -t "$THREADS" -o fastqc_out fastq_trimmed/*trimmed.fq.gz || true

# ----------------------------------------------------------------------------
# 4. Merge per-sample BAMs (one merged BAM per unique sample id)
# ----------------------------------------------------------------------------
log "Merging per-sample BAMs"
awk -F'|' '{print $1}' logs/.sample_bam_map.tsv | sort -u | while read -r sample; do
    out="bam_merged/${sample}.sorted.bam"
    if [[ -f "${out}.bai" ]]; then
        log "Skip merge for $sample (exists)"
        continue
    fi
    mapfile -t bams < <(awk -F'|' -v s="$sample" '$1==s{print $2}' logs/.sample_bam_map.tsv)
    log "Merge $sample <- ${#bams[@]} BAMs"
    if [[ ${#bams[@]} -eq 1 ]]; then
        cp "${bams[0]}"      "$out"
        cp "${bams[0]}.bai"  "${out}.bai"
    else
        samtools merge -@ "$THREADS" -f -o "$out" "${bams[@]}"
        samtools index -@ "$THREADS" "$out"
    fi
done

rm -f logs/.sample_bam_map.tsv
log "Pipeline complete. Merged BAMs in: $WORK_DIR/bam_merged/"
log "Next step: run the DiRTv2 R analysis pointing --bam_dir at bam_merged/"
