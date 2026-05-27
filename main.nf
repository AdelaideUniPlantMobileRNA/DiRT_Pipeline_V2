#!/usr/bin/env nextflow
/*
 * DiRT v2 — Dicistronic Transcript Detection pipeline
 *
 * Generic version for any plant species. Inputs:
 *   --samplesheet  CSV with: sample,fastq_1,fastq_2     (paired-end Illumina RNA-seq)
 *   --genome       reference genome FASTA               (any plant species)
 *   --gff          gene annotation GFF3                 (e.g. from EnsemblPlants)
 *   --outdir       results directory                    (default: ./results)
 *
 * Optional (provide one to skip rebuild):
 *   --hisat2_index prebuilt HISAT2 index prefix (no .ht2 suffix)
 *
 * Multiple FASTQ pairs that share the same `sample` value in samples.csv are
 * automatically merged into one BAM (biological replicate). The set of unique
 * `sample` values becomes the column names used in the DiRTv2 R analysis.
 *
 * Usage:
 *   nextflow run main.nf -profile standard \
 *     --samplesheet samples.csv \
 *     --genome      /path/to/genome.fasta \
 *     --gff         /path/to/annotation.gff3
 */
nextflow.enable.dsl = 2

include { ADAPTER_REMOVAL  } from './modules/adapter_removal.nf'
include { FASTQC           } from './modules/fastqc.nf'
include { HISAT2_BUILD     } from './modules/hisat2_build.nf'
include { HISAT2_ALIGN     } from './modules/hisat2.nf'
include { SAMTOOLS_MERGE   } from './modules/samtools_merge.nf'
include { TRNASCAN_SE      } from './modules/trnascan.nf'
include { GENOME_FAIDX     } from './modules/faidx.nf'
include { DIRTV2_ANALYSIS  } from './modules/dirtv2_analysis.nf'

workflow {

    // --------------------------------------------------------------------
    // 0. Validate inputs
    // --------------------------------------------------------------------
    if (!params.samplesheet) error "Missing --samplesheet <CSV: sample,fastq_1,fastq_2>"
    if (!params.genome)      error "Missing --genome      <reference FASTA>"
    if (!params.gff)         error "Missing --gff         <annotation GFF3, e.g. EnsemblPlants>"

    // --------------------------------------------------------------------
    // 1. Samples channel: (sample_id, run_id, R1, R2)
    //    `run_id` is auto-generated so multiple lanes/runs of the same
    //    biological sample can be tracked independently before merging.
    // --------------------------------------------------------------------
    reads_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map{ row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2)
                error "samplesheet must have columns: sample,fastq_1,fastq_2"
            def r1 = file(row.fastq_1, checkIfExists: true)
            def r2 = file(row.fastq_2, checkIfExists: true)
            // run_id ensures uniqueness when one sample has multiple FASTQ pairs
            def run_id = "${row.sample}_${r1.simpleName}"
            tuple(row.sample, run_id, r1, r2)
        }

    // --------------------------------------------------------------------
    // 2. AdapterRemoval per run (paired-end)
    // --------------------------------------------------------------------
    trimmed_ch = ADAPTER_REMOVAL(reads_ch)

    // --------------------------------------------------------------------
    // 3. FastQC on trimmed reads (collected)
    // --------------------------------------------------------------------
    FASTQC(trimmed_ch.flatMap{ s, r, r1, r2, sg -> [r1, r2] }.collect())

    // --------------------------------------------------------------------
    // 4. HISAT2 index — build if not supplied
    // --------------------------------------------------------------------
    if (params.hisat2_index) {
        hisat2_idx = Channel.fromPath("${params.hisat2_index}*").collect()
    } else {
        hisat2_idx = HISAT2_BUILD(Channel.fromPath(params.genome)).collect()
    }

    // --------------------------------------------------------------------
    // 5. HISAT2 align → sorted, indexed BAM per run
    // --------------------------------------------------------------------
    aligned_ch = HISAT2_ALIGN(trimmed_ch, hisat2_idx)

    // --------------------------------------------------------------------
    // 6. Merge per-replicate BAMs (groups by sample_id)
    // --------------------------------------------------------------------
    bam_grouped = aligned_ch
        .map{ s, r, bam, bai -> tuple(s, bam, bai) }
        .groupTuple()                   // [sample, [bams...], [bais...]]

    merged_ch = SAMTOOLS_MERGE(bam_grouped)

    // --------------------------------------------------------------------
    // 7. Index genome FASTA + run tRNAscan-SE (once)
    // --------------------------------------------------------------------
    fasta_ch = Channel.fromPath(params.genome)
    fai_ch   = GENOME_FAIDX(fasta_ch)
    trna_ch  = TRNASCAN_SE(fasta_ch)

    // --------------------------------------------------------------------
    // 8. Collect all merged BAMs + sample names → DiRTv2 R analysis
    // --------------------------------------------------------------------
    // Sort merged BAMs by sample name so column order is deterministic
    merged_sorted = merged_ch
        .toSortedList{ a, b -> a[0] <=> b[0] }
        .flatMap()

    sample_names_ch = merged_sorted.map{ s, bam, bai -> s }.collect()
    bam_paths_ch    = merged_sorted.map{ s, bam, bai -> [bam, bai] }.flatten().collect()

    DIRTV2_ANALYSIS(
        bam_paths_ch,
        sample_names_ch,
        Channel.fromPath(params.gff).first(),
        fai_ch.first(),
        trna_ch.first(),
        Channel.fromPath("${projectDir}/bin/DiRTv2_manual_optimized.Rmd").first()
    )
}

workflow.onComplete {
    log.info "=============================================="
    log.info "DiRT v2 pipeline finished."
    log.info "Status:  ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Results: ${params.outdir}"
    log.info "Tables:  ${params.outdir}/dirtv2/results/"
    log.info "  - Final_Result_Manual_optimized.xlsx"
    log.info "Report:  ${params.outdir}/dirtv2/DiRTv2_manual_optimized.html"
    log.info "=============================================="
}
