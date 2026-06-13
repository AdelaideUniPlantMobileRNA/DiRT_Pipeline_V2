#!/usr/bin/env bash
###############################################################################
# DiRT v2 — Convert GTF (or non-hierarchical GFF3) → properly hierarchical
#            GFF3 with gene→mRNA→exon parents (helper for Path B)
#
# The DiRTv2 R Markdown parses GFF3-style `ID=` / `Parent=` attributes and
# requires explicit `gene` rows that parent the `mRNA` rows. Plain GTF input
# (gene_id "..."; transcript_id "..." attribute style) and AUGUSTUS-style
# GFF3 output (transcript rows with no `gene` parent) both fail in the Rmd's
# `build-canonical-transcripts` chunk. This script normalises either input
# into the GFF3 format the Rmd expects using AGAT
# (Dainat et al., https://github.com/NBISweden/AGAT).
#
# For Path A (Nextflow), this is currently NOT done automatically — if your
# annotation needs conversion, run this script first, then pass the produced
# GFF3 as the `--gff` argument to `nextflow run main.nf`.
#
# Usage:
#   ./03_convert_gtf_to_gff3.sh <input.gtf_or_gff3> [output.gff3]
#
# Example:
#   ./03_convert_gtf_to_gff3.sh annotation.gtf annotation.AGAT.gff3
#
# Output:
#   <output.gff3>     A GFF3 with proper gene → mRNA → exon / CDS hierarchy
#   <output.agat.log> AGAT run log
#
# Requirements:
#   AGAT v1.0+ (from the dirtv2 conda env, or `mamba install -c bioconda agat`)
#
# Wall time: seconds to a few minutes on a typical plant annotation.
###############################################################################
set -euo pipefail

INPUT="${1:?Usage: $0 <input.gtf_or_gff3> [output.gff3]}"

# Default output name: same basename, force .AGAT.gff3 suffix
DEFAULT_OUT="$(dirname "$(readlink -f "$INPUT")")/$(basename "${INPUT%.*}").AGAT.gff3"
OUTPUT="${2:-$DEFAULT_OUT}"
LOGFILE="${OUTPUT%.gff3}.agat.log"

[[ -f "$INPUT" ]] || { echo "ERROR: input annotation not found: $INPUT" >&2; exit 1; }
command -v agat_convert_sp_gxf2gxf.pl >/dev/null 2>&1 || {
    echo "ERROR: AGAT not on PATH. Activate the dirtv2 conda env first:" >&2
    echo "       conda activate dirtv2" >&2
    echo "Or install AGAT into the active env:" >&2
    echo "       mamba install -c bioconda agat -y" >&2
    exit 1
}

log() { echo "[$(date '+%F %T')] $*"; }

# Skip if output already exists (idempotent)
if [[ -f "$OUTPUT" && -s "$OUTPUT" ]]; then
    log "Output already exists: $OUTPUT"
    log "Delete it if you want to regenerate."
    exit 0
fi

# 1. Run AGAT conversion
log "Converting $INPUT -> $OUTPUT"
log "(AGAT will add missing gene rows by grouping transcripts; warnings about"
log " 'No Parent attribute found' on AUGUSTUS-style input are expected and safe.)"
agat_convert_sp_gxf2gxf.pl \
    --gff "$INPUT" \
    -o   "$OUTPUT" 2>&1 | tee "$LOGFILE"

# 2. Sanity checks
n_gene=$(awk -F"\t" '$3 == "gene"' "$OUTPUT" | wc -l)
n_mrna=$(awk -F"\t" '$3 == "mRNA"' "$OUTPUT" | wc -l)
n_exon=$(awk -F"\t" '$3 == "exon"' "$OUTPUT" | wc -l)
n_cds=$(awk  -F"\t" '$3 == "CDS"'  "$OUTPUT" | wc -l)
log "Done."
log "  gene rows: $n_gene"
log "  mRNA rows: $n_mrna"
log "  exon rows: $n_exon"
log "  CDS  rows: $n_cds"

if [[ "$n_gene" -eq 0 ]]; then
    echo "WARNING: 0 gene rows in output — AGAT did not create the gene hierarchy." >&2
    echo "         The DiRTv2 Rmd will likely fail. Check $LOGFILE for AGAT errors." >&2
    exit 2
fi

log "Output GFF3: $OUTPUT"
log "Pass this file to the DiRTv2 Rmd as the \`gff_file\` parameter."
