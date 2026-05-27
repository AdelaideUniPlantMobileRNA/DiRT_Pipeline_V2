process DIRTV2_ANALYSIS {
    tag "DiRTv2 R analysis (manual optimized)"
    publishDir "${params.outdir}/dirtv2", mode: 'copy'

    input:
    path bam_files               // all merged BAMs + .bai (flat list)
    val  sample_names            // list of sample names (same order as bam_files pairs)
    path gff
    path fai
    path trna_bed
    path rmd_template

    output:
    path "results/*"
    path "DiRTv2_manual_optimized.html"

    script:
    def sample_str = sample_names.collect{ "\"${it}\"" }.join(', ')
    """
    # Stage all merged BAMs into a single directory so the .Rmd can scan it.
    # The .Rmd expects '<sample>.sorted.bam' as the basename — that's what
    # SAMTOOLS_MERGE produces upstream, so symlinking by basename is enough.
    mkdir -p bam_in
    for f in ${bam_files}; do
        ln -sf \$(readlink -f \$f) bam_in/\$(basename \$f)
    done

    cp ${rmd_template} DiRTv2_manual_optimized.Rmd

    Rscript -e "rmarkdown::render('DiRTv2_manual_optimized.Rmd', \\
        params = list( \\
            gff_file       = '${gff}', \\
            genome_fai     = '${fai}', \\
            tRNA_bed       = '${trna_bed}', \\
            bam_dir        = 'bam_in', \\
            bam_pattern    = '\\\\.sorted\\\\.bam\$', \\
            sample_names   = c(${sample_str}), \\
            min_count      = ${params.min_count}, \\
            fdr_threshold  = ${params.fdr_threshold}, \\
            out_dir        = 'results' \\
        ), \\
        output_file = 'DiRTv2_manual_optimized.html')"
    """
}
