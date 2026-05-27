process HISAT2_BUILD {
    tag "$fasta.simpleName"
    publishDir "${params.outdir}/hisat2_index", mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta.simpleName}_index.*.ht2"

    script:
    """
    hisat2-build -p ${task.cpus} ${fasta} ${fasta.simpleName}_index
    """
}
