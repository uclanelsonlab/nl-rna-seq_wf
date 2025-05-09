nextflow.enable.dsl = 2

log.info """\
    R N A - S E Q _ W F   P I P E L I N E
    ===================================
    prefix              : ${params.prefix}
    family_id           : ${params.family_id}
    bucket_dir          : ${params.bucket_dir}
    rib_reference_path  : ${params.rib_reference_path}
    outdir              : ${params.outdir}
    bucket output       : ${params.output_bucket}
    """
    .stripIndent(true)

include { RNASEQC } from './modules/rnaseqc.nf'
include { UPLOAD_FILES } from './modules/upload_outputs.nf'
include { IRFINDER } from './modules/irfinder.nf'
include { BAM2SJ } from './modules/bam2sj/main.nf'
include { MOSDEPTH_BED } from './modules/mosdepth/main.nf'
include { RUN_FASTP } from './modules/fastp.nf'
include { FILTER_FASTQ } from './modules/filters.nf'
include { CHECK_STAR_REF; STAR_ALIGNREADS } from './modules/star.nf'
include { run_markdup } from './modules/picard.nf'
include { SAMBAMBA_MARKDUP } from './modules/sambamba.nf'
include { BWA_MEM as BWA_MEM_RRNA; BWA_MEM as BWA_MEM_GLOBINRNA } from './modules/bwa.nf'
include { 
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_GLOBINRNA; 
    SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_GLOBINRNA } from './modules/samtools.nf'
include { 
    DOWNLOAD_GENCODE as DOWNLOAD_GENCODE_NORMAL; 
    DOWNLOAD_GENCODE as DOWNLOAD_GENCODE_COLLAPSE; 
    SUBREAD_FEATURECOUNTS } from './modules/subreads.nf'
include { 
    DOWNLOAD_MASTER_FEATURECOUNTS; 
    ADD_SAMPLE_COUNTS_MASTER; 
    RUN_OUTRIDER } from './modules/outrider/main.nf'
include { 
    GATK4_SPLITNCIGARREADS; 
    GATK4_BASERECALIBRATOR; 
    GATK4_APPLYBQSR; 
    GATK4_HAPLOTYPECALLER } from './modules/gatk/main.nf'
include { 
    DOWNLOAD_FASTQS; 
    DOWNLOAD_RNA_REF as DOWNLOAD_RRNA; 
    DOWNLOAD_RNA_REF as DOWNLOAD_GLOBINRNA; 
    DOWNLOAD_HUMAN_REF; 
    DOWNLOAD_IR_REF; 
    DOWNLOAD_BED;
    DOWNLOAD_CRAM } from './modules/download_files.nf'
include { 
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_RRNA; 
    SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_RRNA; 
    SAMTOOLS_INDEX; 
    SAMTOOLS_CRAM;
    SAMTOOLS_BAM2SAM;
    SAMTOOLS_CRAM2BAM } from './modules/samtools.nf'

workflow {
    DOWNLOAD_RRNA(params.rib_reference_path, "rrna") //DOWNLOAD_RRNA_ch
    DOWNLOAD_GLOBINRNA(params.rib_reference_path, "globinrna") //DOWNLOAD_GLOBINRNA_ch
    DOWNLOAD_HUMAN_REF(params.human_fasta, params.human_fai, params.human_dict)
    DOWNLOAD_BED(params.bed)

    ch_dbsnp = Channel.value([[id:"dbsnp138"], params.dbsnp138, params.dbsnp138_index])
    ch_known_indels = Channel.value([[id:"known_indels"], params.known_indels, params.known_indels_index])
    ch_indels_1000G = Channel.value([[id:"indels_1000G"], params.indels_1000G, params.indels_1000G_index])
    ch_af_only_gnomad = Channel.value([[id:"af_only_gnomad"], params.af_only_gnomad, params.af_only_gnomad_index])
    ch_small_exac_common_3 = Channel.value([[id:"small_exac_common_3"], params.small_exac_common_3, params.small_exac_common_3_index])

    
    if (params.use_cram) {
        // download cram 
        DOWNLOAD_CRAM(params.prefix, params.cram)
        // CRAM to BAM
        SAMTOOLS_CRAM2BAM(
            DOWNLOAD_HUMAN_REF.out.human_fasta,
            DOWNLOAD_HUMAN_REF.out.human_fai,
            DOWNLOAD_HUMAN_REF.out.human_dict,
            DOWNLOAD_CRAM.out.cram)
        // transform it to BAM to be mark_dup_ch
        SAMBAMBA_MARKDUP(SAMTOOLS_CRAM2BAM.out.bam)
    } else {
        DOWNLOAD_FASTQS(params.prefix, params.fastq_r1, params.fastq_r2) //DOWNLOAD_FASTQS_ch
        // contamination check
        RUN_FASTP(DOWNLOAD_FASTQS.out.reads) //fastp_ch
        FILTER_FASTQ(RUN_FASTP.out.reads) //filtered_fastq_ch
        BWA_MEM_RRNA(FILTER_FASTQ.out.reads, DOWNLOAD_RRNA.out.reference_dir, "human_rRNA_strict.fasta", "rrna") //rrna_bwa_ch
        BWA_MEM_GLOBINRNA(FILTER_FASTQ.out.reads, DOWNLOAD_GLOBINRNA.out.reference_dir, "human_globinRNA.fa", "globinrna") //globinrna_bwa_ch
        SAMTOOLS_VIEW_RRNA(BWA_MEM_RRNA.out.bwa_bam, "rrna") //rrna_SAMTOOLS_VIEW_ch
        SAMTOOLS_VIEW_GLOBINRNA(BWA_MEM_GLOBINRNA.out.bwa_bam, "globinrna") //globinrna_SAMTOOLS_VIEW_ch
        SAMTOOLS_FLAGSTAT_RRNA(SAMTOOLS_VIEW_RRNA.out.view_bam, "rrna") //rrna_SAMTOOLS_FLAGSTAT_ch
        SAMTOOLS_FLAGSTAT_GLOBINRNA(SAMTOOLS_VIEW_GLOBINRNA.out.view_bam, "globinrna") //globinrna_SAMTOOLS_FLAGSTAT_ch

        // STAR alignment
        CHECK_STAR_REF(DOWNLOAD_FASTQS.out.reads) //star_index_ref_ch
        STAR_ALIGNREADS(
            CHECK_STAR_REF.out.star_index,
            CHECK_STAR_REF.out.sjdb_overhang,
            RUN_FASTP.out.reads) //STAR_ALIGNREADS_ch
        SAMTOOLS_INDEX(STAR_ALIGNREADS.out.star_bam)
        SAMBAMBA_MARKDUP(STAR_ALIGNREADS.out.star_bam) //mark_dup_ch
    }

    // Create counts by gene
    DOWNLOAD_GENCODE_NORMAL(params.gencode_gtf_path) //gencode_pc_ch
    SUBREAD_FEATURECOUNTS(
        DOWNLOAD_GENCODE_NORMAL.out.gencode_gtf, 
        SAMBAMBA_MARKDUP.out.marked_bam) //feature_counts_ch
    
    // Run OUTRIDER
    if ( params.outrider ) {
        DOWNLOAD_MASTER_FEATURECOUNTS(params.features_master_file) //featurecounts_master_ch
        ADD_SAMPLE_COUNTS_MASTER(
            DOWNLOAD_MASTER_FEATURECOUNTS.out.featurecounts_master, 
            SUBREAD_FEATURECOUNTS.out.gene_counts_short) //featurecounts_updated_ch
        RUN_OUTRIDER(ADD_SAMPLE_COUNTS_MASTER.out.featurecounts_updated, params.tissue) //outrider_table_ch
    }
    // Run QC
    DOWNLOAD_GENCODE_COLLAPSE(params.gencode_gtf_collapse) //gencode_collapse_ch
    RNASEQC(
        DOWNLOAD_GENCODE_COLLAPSE.out.gencode_gtf, 
        SAMBAMBA_MARKDUP.out.marked_bam) //rnaseqc_ch

    // Create CRAM files
    SAMTOOLS_CRAM(
        DOWNLOAD_HUMAN_REF.out.human_fasta, 
        DOWNLOAD_HUMAN_REF.out.human_fai, 
        DOWNLOAD_HUMAN_REF.out.human_dict, 
        SAMBAMBA_MARKDUP.out.marked_bam) //cram_ch

    // Run IRFinder
    DOWNLOAD_IR_REF(params.ir_ref) //ir_ref_ch
    IRFINDER(
        DOWNLOAD_IR_REF.out.ir_reference_dir, 
        SAMBAMBA_MARKDUP.out.marked_bam) //irfinder_ch

    // Calculate XBP1 coverage
    MOSDEPTH_BED(
        DOWNLOAD_HUMAN_REF.out.human_fasta, 
        DOWNLOAD_HUMAN_REF.out.human_fai, 
        DOWNLOAD_HUMAN_REF.out.human_dict, 
        DOWNLOAD_BED.out.bed, 
        SAMTOOLS_CRAM.out.rna_cram,
        SAMTOOLS_CRAM.out.rna_crai)
    
    // CRAM to SAM 
    SAMTOOLS_BAM2SAM(
        DOWNLOAD_HUMAN_REF.out.human_fasta, 
        DOWNLOAD_HUMAN_REF.out.human_fai, 
        DOWNLOAD_HUMAN_REF.out.human_dict, 
        SAMBAMBA_MARKDUP.out.marked_bam)
    
    // SAM to SJ
    BAM2SJ(SAMTOOLS_BAM2SAM.out.rna_sam)

    // GATK variant calling
    GATK4_SPLITNCIGARREADS(
        DOWNLOAD_HUMAN_REF.out.human_fasta, 
        DOWNLOAD_HUMAN_REF.out.human_fai, 
        DOWNLOAD_HUMAN_REF.out.human_dict, 
        SAMBAMBA_MARKDUP.out.marked_bam)
    GATK4_BASERECALIBRATOR(
        GATK4_SPLITNCIGARREADS.out.bam, 
        DOWNLOAD_HUMAN_REF.out.human_fasta, 
        DOWNLOAD_HUMAN_REF.out.human_fai, 
        DOWNLOAD_HUMAN_REF.out.human_dict,  
        ch_dbsnp,
        ch_known_indels,
        ch_indels_1000G,
        ch_af_only_gnomad,
        ch_small_exac_common_3)
    GATK4_APPLYBQSR(
        GATK4_SPLITNCIGARREADS.out.bam, 
        GATK4_BASERECALIBRATOR.out.table, 
        DOWNLOAD_HUMAN_REF.out.human_fasta,
        DOWNLOAD_HUMAN_REF.out.human_fai, 
        DOWNLOAD_HUMAN_REF.out.human_dict)
    GATK4_HAPLOTYPECALLER(GATK4_APPLYBQSR.out.bam, 
        DOWNLOAD_HUMAN_REF.out.human_fasta, 
        DOWNLOAD_HUMAN_REF.out.human_fai, 
        DOWNLOAD_HUMAN_REF.out.human_dict,
        ch_dbsnp)
    
    // Upload selected output files
    UPLOAD_FILES(
        params.family_id, 
        params.bucket_dir, 
        params.output_bucket,
        RUN_FASTP.out.json, RUN_FASTP.out.html, RUN_FASTP.out.log, RUN_FASTP.out.versions, //fastp
        SAMTOOLS_FLAGSTAT_RRNA.out.flagstat_file, SAMTOOLS_FLAGSTAT_RRNA.out.versions, //flagstat_rrna
        SAMTOOLS_FLAGSTAT_GLOBINRNA.out.flagstat_file, SAMTOOLS_FLAGSTAT_GLOBINRNA.out.versions, //flagstat_globinrna
        STAR_ALIGNREADS.out.reads_gene, STAR_ALIGNREADS.out.reads_gene_log, STAR_ALIGNREADS.out.final_log, STAR_ALIGNREADS.out.sj_tab, STAR_ALIGNREADS.out.log, STAR_ALIGNREADS.out.versions, //star
        SAMBAMBA_MARKDUP.out.log, SAMBAMBA_MARKDUP.out.versions, //sambamba
        SUBREAD_FEATURECOUNTS.out.gene_counts, SUBREAD_FEATURECOUNTS.out.gene_counts_short, SUBREAD_FEATURECOUNTS.out.gene_counts_summary, SUBREAD_FEATURECOUNTS.out.log, SUBREAD_FEATURECOUNTS.out.versions, //featurecounts
        RNASEQC.out.coverage, RNASEQC.out.exon_cv, RNASEQC.out.exon_reads, RNASEQC.out.gene_fragments, RNASEQC.out.gene_reads, RNASEQC.out.gene_tpm, RNASEQC.out.metrics, RNASEQC.out.log, RNASEQC.out.versions, //rnaseqc
        SAMTOOLS_CRAM.out.rna_cram, SAMTOOLS_CRAM.out.rna_crai, SAMTOOLS_CRAM.out.log, SAMTOOLS_CRAM.out.versions, //cram
        IRFINDER.out.irfinder_chr_coverage, IRFINDER.out.irfinder_dir_val, IRFINDER.out.irfinder_dir, IRFINDER.out.irfinder_nondir_val, IRFINDER.out.irfinder_nondir, IRFINDER.out.irfinder_junc_count, IRFINDER.out.irfinder_roi, IRFINDER.out.irfinder_spans_point, IRFINDER.out.log, IRFINDER.out.versions, //irfinder
        SAMTOOLS_BAM2SAM.out.log, SAMTOOLS_BAM2SAM.out.versions, //sam
        BAM2SJ.out.sj_tab_gz,
        MOSDEPTH_BED.out.global_dist, MOSDEPTH_BED.out.region_dist, MOSDEPTH_BED.out.summary, MOSDEPTH_BED.out.perbase, MOSDEPTH_BED.out.perbase_index, MOSDEPTH_BED.out.regions_bed, MOSDEPTH_BED.out.regions_bed_index,
        GATK4_HAPLOTYPECALLER.out.vcf, GATK4_HAPLOTYPECALLER.out.tbi, GATK4_HAPLOTYPECALLER.out.versions
    )
}
