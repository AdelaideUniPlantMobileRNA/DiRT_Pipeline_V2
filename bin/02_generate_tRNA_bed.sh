#!/usr/bin/env bash
###############################################################################
# DiRT v2 — Generate tRNA.bed from a genome FASTA (helper for Path B)
#
# Runs tRNAscan-SE on the genome, then reformats the tabular output into the
# 6-column BED file (chr, start, end, amino_acid, score, anticodon) that the
# DiRTv2 R Markdown expects.
#
# For Path A (Nextflow), this is done automatically by modules/trnascan.nf —
# you do NOT need to run this script.
#
# Usage:
#   ./02_generate_tRNA_bed.sh <genome.fasta> [output_dir]
#
# Output (in <output_dir>, defaulting to the genome's directory):
#   tRNA.txt          tRNAscan-SE tabular output (raw)
#   tRNA_stats.txt    tRNAscan-SE search statistics
#   tRNA.bed          6-column BED — pass this as `tRNA_bed` to the Rmd
#   trnascan.log      Run log
#
# Requirements:
#   tRNAscan-SE (from the dirtv2 conda env)
#   awk (any POSIX implementation)
#
# Wall time: ~2-4 hours on a typical plant genome (~500 Mb).
###############################################################################
set -euo pipefail

FASTA="${1:?Usage: $0 <genome.fasta> [output_dir]}"
OUT_DIR="${2:-$(dirname "$(readlink -f "$FASTA")")}"

[[ -f "$FASTA" ]] || { echo "ERROR: genome FASTA not found: $FASTA" >&2; exit 1; }
command -v tRNAscan-SE >/dev/null 2>&1 || {
    echo "ERROR: tRNAscan-SE not on PATH. Activate the dirtv2 conda env first:" >&2
    echo "       conda activate dirtv2" >&2
    exit 1
}

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

log() { echo "[$(date '+%F %T')] $*"; }

# Skip if already generated (idempotent)
if [[ -f tRNA.bed && -s tRNA.bed ]]; then
    log "tRNA.bed already exists ($(wc -l < tRNA.bed) entries). Skipping."
    log "Delete it if you want to regenerate."
    exit 0
fi

# 1. Run tRNAscan-SE
log "Running tRNAscan-SE on $FASTA (this takes ~2-4 hours on a plant genome)"
tRNAscan-SE -E \
    -o tRNA.txt \
    -m tRNA_stats.txt \
    "$FASTA" 2>&1 | tee trnascan.log

# 2. Convert tabular output → BED (skip the 3-line header)
log "Converting tRNA.txt → tRNA.bed"
awk 'NR > 3 {print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $9 "\t" $6}' tRNA.txt > tRNA.bed

# 3. Sanity check
n_trnas=$(wc -l < tRNA.bed)
n_minus=$(awk '$2 > $3' tRNA.bed | wc -l)
log "Done. Detected ${n_trnas} tRNAs (${n_minus} on the minus strand)."
log "Note: minus-strand tRNAs are reported with start > end — this is normal"
log "      tRNAscan-SE output. The DiRTv2 Rmd normalizes coordinates internally."
log "Output: $(pwd)/tRNA.bed"
