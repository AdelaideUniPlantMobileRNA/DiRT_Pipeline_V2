process SAMTOOLS_MERGE {
    tag "$sample"
    publishDir "${params.outdir}/bam_merged", mode: 'copy'

    input:
    tuple val(sample), path(bams), path(bais)

    output:
    tuple val(sample),
          path("${sample}.sorted.bam"),
          path("${sample}.sorted.bam.bai")

    script:
    """
    samtools merge -@ ${task.cpus} -f -o ${sample}.sorted.bam ${bams}
    samtools index -@ ${task.cpus} ${sample}.sorted.bam
    """
}
