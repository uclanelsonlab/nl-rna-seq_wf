nextflow.enable.dsl = 2

log.info """\
    R N A - S E Q _ W F   P I P E L I N E
    ===================================
    prefix         : ${params.prefix}
    family_id      : ${params.family_id}
    bucket_dir     : ${params.bucket_dir}
    output_bucket  : ${params.output_bucket}
    fastq_r1       : ${params.fastq_r1}
    fastq_r2       : ${params.fastq_r2}
    """
    .stripIndent(true)

include { DOWNLOAD_FASTQS; download_rna_ref as download_rrna; download_rna_ref as download_globinrna; download_human_ref } from './modules/download_files.nf'
include { run_fastp } from './modules/fastp.nf'
include { filter_fastq } from './modules/filters.nf'
include { bwa_mem as bwa_mem_rrna; bwa_mem as bwa_mem_globinrna } from './modules/bwa.nf'
include { samtools_view as samtools_view_rrna; samtools_flagstat as samtools_flagstat_rrna; samtools_index; samtools_cram } from './modules/samtools.nf'
include { samtools_view as samtools_view_globinrna; samtools_flagstat as samtools_flagstat_globinrna; samtools_view_sj } from './modules/samtools.nf'
include { check_star_reference; star_alignreads } from './modules/star.nf'
include { bam2sj } from './modules/bam2sj/main.nf'
include { prioritize_splice_junctions } from './modules/prioritize_splice_junctions/main.nf'
include { download_gencode; subread_featurecounts } from './modules/subreads.nf'
include { upload_files } from './modules/upload_outputs.nf'

workflow {
    download_fastqs_ch = DOWNLOAD_FASTQS(params.prefix, params.fastq_r1, params.fastq_r2)
    download_rrna_ch = download_rrna(params.rib_reference_path, "rrna")
    download_globinrna_ch = download_globinrna(params.rib_reference_path, "globinrna")

    // contamination check
    fastp_ch = run_fastp(download_fastqs_ch)
    filtered_fastq_ch = filter_fastq(fastp_ch)
    rrna_bwa_ch = bwa_mem_rrna(filtered_fastq_ch, download_rrna_ch, "human_rRNA_strict.fasta", "rrna")
    globinrna_bwa_ch = bwa_mem_globinrna(filtered_fastq_ch, download_globinrna_ch, "human_globinRNA.fa", "globinrna")
    rrna_samtools_view_ch = samtools_view_rrna(params.prefix, rrna_bwa_ch, "rrna")
    globinrna_samtools_view_ch = samtools_view_globinrna(params.prefix, globinrna_bwa_ch, "globinrna")
    rrna_samtools_flagstat_ch = samtools_flagstat_rrna(params.prefix, rrna_samtools_view_ch, "rrna")
    globinrna_samtools_flagstat_ch = samtools_flagstat_globinrna(params.prefix, globinrna_samtools_view_ch, "globinrna")

    // STAR alignment
    star_index_ref_ch = check_star_reference(download_fastqs_ch)
    star_alignreads_ch = star_alignreads(star_index_ref_ch, fastp_ch)
    samtools_index(star_alignreads_ch)

    // Create SJ tab file
    sam_ch = samtools_view_sj(params.prefix, star_alignreads_ch) 
    sj_tab_ch = bam2sj(params.prefix, sam_ch)

    // Create spreadsheet with prioritize splice junctions
    splice_junctions_ch = prioritize_splice_junctions(params.prefix, sj_tab_ch, params.latest_directory_date, params.latest_udn_id_key, params.sjdblist, params.genemap, params.constraint, params.gencode_cds)

    // Create counts by gene
    gencode_pc_ch = download_gencode(params.gencode_gtf_path)
    feature_counts_ch = subread_featurecounts(params.prefix, gencode_pc_ch, star_alignreads_ch)

    // Create CRAM files
    download_human_ref_ch = download_human_ref(params.human_fasta, params.human_fai, params.human_dict)
    cram_ch = samtools_cram(params.prefix, download_human_ref_ch, star_alignreads_ch)

    // Upload selected output files
    upload_files(
        params.family_id, 
        params.bucket_dir, 
        params.output_bucket,
        rrna_samtools_flagstat_ch, 
        globinrna_samtools_flagstat_ch, 
        star_alignreads_ch, 
        sj_tab_ch, 
        splice_junctions_ch, 
        feature_counts_ch, 
        cram_ch, 
        params.output_bucket)
}