nextflow.enable.dsl = 2

log.info """\
    R N A - S E Q _ W F   P I P E L I N E
    ===================================
    sample_name         : ${params.sample_name}
    library             : ${params.library}
    fastq_bucket        : ${params.fastq_bucket}
    rib_reference_path  : ${params.rib_reference_path}
    outdir              : ${params.outdir}
    bucket output       : ${params.output_bucket}
    """
    .stripIndent(true)

include { 
    DOWNLOAD_FASTQS; 
    DOWNLOAD_RNA_REF as DOWNLOAD_RRNA; 
    DOWNLOAD_RNA_REF as DOWNLOAD_GLOBINRNA; 
    DOWNLOAD_HUMAN_REF; DOWNLOAD_IR_REF; DOWNLOAD_BED } from './modules/download_files.nf'
include { RUN_FASTP } from './modules/fastp.nf'
include { FILTER_FASTQ } from './modules/filters.nf'
include { BWA_MEM as BWA_MEM_RRNA; BWA_MEM as BWA_MEM_GLOBINRNA } from './modules/bwa.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_RRNA; SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_RRNA; SAMTOOLS_INDEX; SAMTOOLS_CRAM } from './modules/samtools.nf'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_GLOBINRNA; SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_GLOBINRNA; SAMTOOLS_BAM2SAM } from './modules/samtools.nf'
include { check_star_reference; star_alignreads } from './modules/star.nf'
include { run_markdup } from './modules/picard.nf'
include { SAMBAMBA_MARKDUP } from './modules/sambamba.nf'
include { download_gencode as DOWNLOAD_GENCODE_NORMAL; download_gencode as DOWNLOAD_GENCODE_COLLAPSE; subread_featurecounts } from './modules/subreads.nf'
include { download_master_featureCounts; add_sample_counts_master; run_outrider } from './modules/outrider/main.nf'
include { RNASEQC } from './modules/rnaseqc.nf'
include { upload_files; UP_SJ } from './modules/upload_outputs.nf'
include { IRFINDER } from './modules/irfinder.nf'
include { BAM2SJ } from './modules/bam2sj/main.nf'
include { MOSDEPTH_BED } from './modules/mosdepth/main.nf'

workflow {
    DOWNLOAD_RRNA(params.rib_reference_path, "rrna") //DOWNLOAD_RRNA_ch
    DOWNLOAD_GLOBINRNA(params.rib_reference_path, "globinrna") //DOWNLOAD_GLOBINRNA_ch
    DOWNLOAD_BED(params.bed)
    
    if (params.use_cram) {
        // download cram 
        // transform it to BAM to be mark_dup_ch
    } else {
        DOWNLOAD_FASTQS(params.meta, params.library, params.fastq_bucket) //DOWNLOAD_FASTQS_ch
        // contamination check
        RUN_FASTP(DOWNLOAD_FASTQS.out.reads) //fastp_ch
        FILTER_FASTQ(RUN_FASTP.out.reads) //filtered_fastq_ch
        BWA_MEM_RRNA(FILTER_FASTQ.out.reads, DOWNLOAD_RRNA.out.reference_dir, "human_rRNA_strict.fasta", "rrna") //rrna_bwa_ch
        BWA_MEM_GLOBINRNA(FILTER_FASTQ.out.reads, DOWNLOAD_GLOBINRNA.out.reference_dir, "human_globinRNA.fa", "globinrna") //globinrna_bwa_ch
        SAMTOOLS_VIEW_RRNA(BWA_MEM_RRNA.out.bwa_bam, "rrna") //rrna_SAMTOOLS_VIEW_ch
        SAMTOOLS_VIEW_GLOBINRNA(BWA_MEM_GLOBINRNA.out.bwa_bam, "globinrna") //globinrna_SAMTOOLS_VIEW_ch
        // rrna_SAMTOOLS_FLAGSTAT_ch = SAMTOOLS_FLAGSTAT_RRNA(rrna_SAMTOOLS_VIEW_ch, "rrna")
        // globinrna_SAMTOOLS_FLAGSTAT_ch = SAMTOOLS_FLAGSTAT_GLOBINRNA(globinrna_SAMTOOLS_VIEW_ch, "globinrna")

        // STAR alignment
        // star_index_ref_ch = check_star_reference(DOWNLOAD_FASTQS_ch)
        // star_alignreads(star_index_ref_ch, fastp_ch) //star_alignreads_ch
        // SAMTOOLS_INDEX(star_alignreads.out.star_bam)
        // SAMBAMBA_MARKDUP(star_alignreads.out.star_bam) //mark_dup_ch
    }

    

    // Create counts by gene
    // gencode_pc_ch = DOWNLOAD_GENCODE_NORMAL(params.gencode_gtf_path)
    // feature_counts_ch = subread_featurecounts(gencode_pc_ch, mark_dup_ch)
    // featurecounts_master_ch = download_master_featureCounts(params.features_master_file)
    // featurecounts_updated_ch = add_sample_counts_master(featurecounts_master_ch, feature_counts_ch)
    // outrider_table_ch = run_outrider(featurecounts_updated_ch, params.tissue)

    // Run QC
    // gencode_collapse_ch = DOWNLOAD_GENCODE_COLLAPSE(params.gencode_gtf_collapse)
    // rnaseqc_ch = RNASEQC(gencode_collapse_ch, mark_dup_ch)

    // Create CRAM files
    // DOWNLOAD_HUMAN_REF(params.human_fasta, params.human_fai, params.human_dict)
    // cram_ch = SAMTOOLS_CRAM(DOWNLOAD_HUMAN_REF.out.human_fasta, DOWNLOAD_HUMAN_REF.out.human_fai, DOWNLOAD_HUMAN_REF.out.human_dict, mark_dup_ch)

    // Run IRFinder
    // ir_ref_ch = DOWNLOAD_IR_REF(params.ir_ref)
    // irfinder_ch = IRFINDER(ir_ref_ch, mark_dup_ch)

    // Calculate XBP1 coverage
    // MOSDEPTH_BED(DOWNLOAD_HUMAN_REF.out.human_fasta, DOWNLOAD_HUMAN_REF.out.human_fai, DOWNLOAD_HUMAN_REF.out.human_dict, DOWNLOAD_BED.out.bed, cram_ch)
    
    // CRAM to SAM 
    // SAMTOOLS_BAM2SAM(DOWNLOAD_HUMAN_REF.out.human_fasta, DOWNLOAD_HUMAN_REF.out.human_fai, DOWNLOAD_HUMAN_REF.out.human_dict, mark_dup_ch)
    
    // SAM to SJ
    // BAM2SJ(SAMTOOLS_BAM2SAM.out.rna_sam)
    
    // Upload selected output files
    // upload_files(
    //     params.sample_name, 
    //     params.proband, 
    //     params.tissue, 
    //     params.output_bucket, 
    //     rrna_SAMTOOLS_FLAGSTAT_ch, 
    //     globinrna_SAMTOOLS_FLAGSTAT_ch, 
    //     star_alignreads_ch, 
    //     feature_counts_ch, 
    //     outrider_table_ch, 
    //     rnaseqc_ch, 
    //     cram_ch, 
    //     irfinder_ch,
    //     fastp_ch
    // )

    // Uplaod SJ
    // UP_SJ(
    //     BAM2SJ.out.sj_tab_gz,
    //     MOSDEPTH_BED.out.global_dist, 
    //     MOSDEPTH_BED.out.region_dist, 
    //     MOSDEPTH_BED.out.summary, 
    //     MOSDEPTH_BED.out.perbase, 
    //     MOSDEPTH_BED.out.perbase_index, 
    //     MOSDEPTH_BED.out.regions_bed, 
    //     MOSDEPTH_BED.out.regions_bed_index,
    //     params.output_bucket
    // )
}
