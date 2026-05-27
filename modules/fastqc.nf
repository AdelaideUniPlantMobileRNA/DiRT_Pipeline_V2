process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    path reads

    output:
    path "fastqc_out/*"

    script:
    """
    mkdir -p fastqc_out
    fastqc -t ${task.cpus} -o fastqc_out ${reads}
    """
}
