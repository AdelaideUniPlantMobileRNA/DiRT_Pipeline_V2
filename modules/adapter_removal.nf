process ADAPTER_REMOVAL {
    tag "$run_id"
    publishDir "${params.outdir}/fastq_trimmed", mode: 'copy'

    input:
    tuple val(sample), val(run_id), path(r1), path(r2)

    output:
    tuple val(sample), val(run_id),
          path("${run_id}_R1.trimmed.fq.gz"),
          path("${run_id}_R2.trimmed.fq.gz"),
          path("${run_id}_singleton.trimmed.fq.gz")

    script:
    """
    AdapterRemoval \\
        --file1 ${r1} \\
        --file2 ${r2} \\
        --minlength ${params.min_length} \\
        --trimns \\
        --trimqualities \\
        --gzip \\
        --threads ${task.cpus} \\
        --output1   ${run_id}_R1.trimmed.fq.gz \\
        --output2   ${run_id}_R2.trimmed.fq.gz \\
        --singleton ${run_id}_singleton.trimmed.fq.gz
    """
}
