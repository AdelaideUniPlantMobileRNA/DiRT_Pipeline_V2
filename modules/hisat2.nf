process HISAT2_ALIGN {
    tag "$run_id"
    publishDir "${params.outdir}/bam_per_run", mode: 'copy', pattern: '*.bam*'
    publishDir "${params.outdir}/logs",         mode: 'copy', pattern: '*_hisat2_summary.txt'

    input:
    tuple val(sample), val(run_id), path(r1), path(r2), path(singleton)
    path  index_files

    output:
    tuple val(sample), val(run_id),
          path("${run_id}.sorted.bam"),
          path("${run_id}.sorted.bam.bai")
    path  "${run_id}_hisat2_summary.txt"

    script:
    // Recover HISAT2 index prefix from staged files
    def idx_prefix = index_files[0].toString().replaceAll(/\.\d+\.ht2$/, '')
    """
    hisat2 -p ${task.cpus} \\
        -x ${idx_prefix} \\
        -1 ${r1} \\
        -2 ${r2} \\
        --rna-strandness RF \\
        --dta \\
        --mp 4,2 \\
        --score-min L,0,-0.4 \\
        --summary-file ${run_id}_hisat2_summary.txt \\
      | samtools view -@ ${task.cpus} -b -F 4 - \\
      | samtools sort -@ ${task.cpus} -m ${params.bam_mem_per_th} \\
                      -o ${run_id}.sorted.bam -

    samtools index -@ ${task.cpus} ${run_id}.sorted.bam
    """
}
