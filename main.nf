nextflow.enable.dsl = 2

params.sample_name = 'NB-8204-M'
params.library = "SN_7RNA_S-24-0479_XA044"
params.proband = "SAMPLE"
params.tissue = "fibroblast"
params.meta = params.sample_name + "-" + params.tissue
params.fastq_bucket = "s3://ucla-rare-diseases/UCLA-UDN/rnaseq/fastq"
params.rib_reference_path = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference"
params.gencode_gtf_path = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/gencode.v43.primary_assembly.annotation.gtf"
params.gencode_gtf_collapse = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/gencode.v43.primary_assembly.annotation.collapse.gtf"
params.outdir = "results"
params.human_fai = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.fa.fai"
params.human_dict = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.dict"
params.human_fasta = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.fa"
params.output_bucket = "s3://ucla-rare-diseases/UCLA-UDN/Analysis/UDN_cases"
params.features_master_file = "s3://ucla-rare-diseases/UCLA-UDN/gcarvalho_test/drop/test/fibroblast/featureCounts_fibroblast_24-07-22.tsv"
params.ir_ref = "s3://ucla-rare-diseases/UCLA-UDN/assets/IRFinder-1.3.1/REF/GRCh38.p13/"
params.bed = "s3://ucla-rare-diseases/UCLA-UDN/gcarvalho_test/rnaseq/hg38/XBP1-NM_001079539-hg38.bed"

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

include { download_fastqs; download_rna_ref as download_rrna; download_rna_ref as download_globinrna; download_human_ref; download_ir_ref; DOWNLOAD_BED } from './modules/download_files.nf'
include { run_fastp } from './modules/fastp.nf'
include { filter_fastq } from './modules/filters.nf'
include { bwa_mem as bwa_mem_rrna; bwa_mem as bwa_mem_globinrna } from './modules/bwa.nf'
include { samtools_view as samtools_view_rrna; samtools_flagstat as samtools_flagstat_rrna; samtools_index; samtools_cram } from './modules/samtools.nf'
include { samtools_view as samtools_view_globinrna; samtools_flagstat as samtools_flagstat_globinrna; SAMTOOLS_CRAM2SAM } from './modules/samtools.nf'
include { check_star_reference; star_alignreads } from './modules/star.nf'
include { run_markdup } from './modules/picard.nf'
include { SAMBAMBA_MARKDUP } from './modules/sambamba.nf'
include { download_gencode as download_gencode_normal; download_gencode as download_gencode_collapse; subread_featurecounts } from './modules/subreads.nf'
include { download_master_featureCounts; add_sample_counts_master; run_outrider } from './modules/outrider/main.nf'
include { RNASEQC } from './modules/rnaseqc.nf'
include { upload_files; UP_SJ } from './modules/upload_outputs.nf'
include {IRFINDER } from './modules/irfinder.nf'
include { BAM2SJ } from './modules/bam2sj/main.nf'
include { MOSDEPTH_BED } from './modules/mosdepth/main.nf'

workflow {
    download_fastqs_ch = download_fastqs(params.meta, params.library, params.fastq_bucket)
    download_rrna_ch = download_rrna(params.rib_reference_path, "rrna")
    download_globinrna_ch = download_globinrna(params.rib_reference_path, "globinrna")
    DOWNLOAD_BED(params.bed)

    // contamination check
    fastp_ch = run_fastp(download_fastqs_ch)
    filtered_fastq_ch = filter_fastq(fastp_ch)
    rrna_bwa_ch = bwa_mem_rrna(filtered_fastq_ch, download_rrna_ch, "human_rRNA_strict.fasta", "rrna")
    globinrna_bwa_ch = bwa_mem_globinrna(filtered_fastq_ch, download_globinrna_ch, "human_globinRNA.fa", "globinrna")
    rrna_samtools_view_ch = samtools_view_rrna(rrna_bwa_ch, "rrna")
    globinrna_samtools_view_ch = samtools_view_globinrna(globinrna_bwa_ch, "globinrna")
    rrna_samtools_flagstat_ch = samtools_flagstat_rrna(rrna_samtools_view_ch, "rrna")
    globinrna_samtools_flagstat_ch = samtools_flagstat_globinrna(globinrna_samtools_view_ch, "globinrna")

    // STAR alignment
    star_index_ref_ch = check_star_reference(download_fastqs_ch)
    star_alignreads_ch = star_alignreads(star_index_ref_ch, fastp_ch)
    samtools_index(star_alignreads_ch)
    mark_dup_ch = SAMBAMBA_MARKDUP(star_alignreads_ch)

    // Create counts by gene
    gencode_pc_ch = download_gencode_normal(params.gencode_gtf_path)
    feature_counts_ch = subread_featurecounts(gencode_pc_ch, mark_dup_ch)
    featurecounts_master_ch = download_master_featureCounts(params.features_master_file)
    featurecounts_updated_ch = add_sample_counts_master(featurecounts_master_ch, feature_counts_ch)
    outrider_table_ch = run_outrider(featurecounts_updated_ch, params.tissue)

    // Run QC
    gencode_collapse_ch = download_gencode_collapse(params.gencode_gtf_collapse)
    rnaseqc_ch = RNASEQC(gencode_collapse_ch, mark_dup_ch)

    // Create CRAM files
    download_human_ref_ch = download_human_ref(params.human_fasta, params.human_fai, params.human_dict)
    cram_ch = samtools_cram(download_human_ref_ch, mark_dup_ch)

    // Run IRFinder
    ir_ref_ch = download_ir_ref(params.ir_ref)
    irfinder_ch = IRFINDER(ir_ref_ch, mark_dup_ch)

    // Calculate XBP1 coverage
    MOSDEPTH_BED(download_human_ref_ch, DOWNLOAD_BED.out.bed, cram_ch)
    
    // CRAM to SAM 
    SAMTOOLS_CRAM2SAM(download_human_ref_ch, cram_ch)
    
    // SAM to SJ
    BAM2SJ(SAMTOOLS_CRAM2SAM.out.rna_sam)
    
    // Upload selected output files
    upload_files(
        params.sample_name, 
        params.proband, 
        params.tissue, 
        params.output_bucket, 
        rrna_samtools_flagstat_ch, 
        globinrna_samtools_flagstat_ch, 
        star_alignreads_ch, 
        feature_counts_ch, 
        outrider_table_ch, 
        rnaseqc_ch, 
        cram_ch, 
        irfinder_ch,
        fastp_ch
    )

    // Uplaod SJ
    UP_SJ(
        BAM2SJ.out.sj_tab_gz,
        MOSDEPTH_BED.out.global_dist, 
        MOSDEPTH_BED.out.region_dist, 
        MOSDEPTH_BED.out.summary, 
        MOSDEPTH_BED.out.perbase, 
        MOSDEPTH_BED.out.perbase_index, 
        MOSDEPTH_BED.out.regions_bed, 
        MOSDEPTH_BED.out.regions_bed_index,
        params.output_bucket
        )
}
