process GENOME_FAIDX {
    tag "$fasta.simpleName"
    publishDir "${params.outdir}/genome", mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta}.fai"

    script:
    """
    samtools faidx ${fasta}
    """
}
