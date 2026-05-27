process TRNASCAN_SE {
    tag "$fasta.simpleName"
    publishDir "${params.outdir}/trnascan", mode: 'copy'

    input:
    path fasta

    output:
    path "tRNA.bed"

    script:
    """
    tRNAscan-SE -E -o tRNA.txt ${fasta}
    awk 'NR > 3 {print \$1 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$9 "\\t" \$6}' tRNA.txt > tRNA.bed
    """
}
