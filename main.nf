nextflow.enable.dsl = 2

log.info """\
    R N A - S E Q _ W F   P I P E L I N E
    ===================================
    sample_name         : ${params.sample_name}
    fastq_r1            : ${params.fastq_r1}
    fastq_r2            : ${params.fastq_r2}
    rrna_reference      : ${params.rrna_reference}
    globinrna_reference : ${params.globinrna_reference}
    """
    .stripIndent(true)

include { FASTP } from './modules/fastp/main.nf'
include { BAM2SJ } from './modules/bam2sj/main.nf'
include { RNASEQC } from './modules/rnaseqc/main.nf'
include { MULTIQC } from './modules/multiqc/main.nf'
include { IRFINDER } from './modules/irfinder/main.nf'
include { ARCASHLA } from './modules/arcashla/main.nf'
include { STAR_ALIGNREADS } from './modules/star/main.nf'
include { MOSDEPTH_BED; MOSDEPTH } from './modules/mosdepth/main.nf'
include { KALLISTO_QUANT } from './modules/kallisto/main.nf'
include { QUALIMAP_RNASEQ } from './modules/qualimap/main.nf'
include { SAMBAMBA_MARKDUP } from './modules/sambamba/main.nf'
include { SUBREAD_FEATURECOUNTS } from './modules/subreads/main.nf'
include { BEDTOOLS_MERGE_INTERSECT } from './modules/bedtools/main.nf'
include { DEEPVARIANT_RUNDEEPVARIANT } from './modules/deepvariant/main.nf'
include { 
    BWA_MEM as BWA_MEM_RRNA; 
    BWA_MEM as BWA_MEM_GLOBINRNA 
    } from './modules/bwa/main.nf'
include { 
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_RRNA; 
    SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_RRNA; 
    SAMTOOLS_VIEW as SAMTOOLS_VIEW_GLOBINRNA;
    SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_GLOBINRNA;
    SAMTOOLS_INDEX; 
    SAMTOOLS_CRAM;
    SAMTOOLS_BAM2SAM;
    SAMTOOLS_CRAM2BAM
    } from './modules/samtools/main.nf'

workflow {
    // Validate required parameters
    if (!params.fastq_r1) {
        error "Missing required parameter: --fastq_r1"
    }
    if (!params.fastq_r2) {
        error "Missing required parameter: --fastq_r2"
    }
    if (!params.rrna_reference) {
        error "Missing required parameter: --rrna_reference"
    }
    if (!params.globinrna_reference) {
        error "Missing required parameter: --globinrna_reference"
    }
    Channel
    .fromPath(params.fastq_r1)
    .map { fastq_r1 ->
        def meta = [:]
        meta.id = params.sample_name ?: fastq_r1.baseName
        def fastq_r2 = file(params.fastq_r2)
        [ meta, [fastq_r1, fastq_r2] ]
    }
    .set { ch_reads }

    Channel.fromPath(params.rrna_reference)
        .map { rrna_zip ->
            ["rrna", rrna_zip]
        }
        .set { ch_rrna_reference }

    Channel.fromPath(params.globinrna_reference)
        .map { globinrna_zip ->
            ["globinrna", globinrna_zip]
        }
        .set { ch_globinrna_reference }
        
    Channel.value([
            [id:"model"],
            file(params.model_data),
            file(params.model_index),
            file(params.model_meta),
            file(params.model_info),
        ])
        .set { ch_model }

    Channel.value([
            [id:"reference"], 
            file(params.human_fasta), 
            file(params.human_fai), 
            file(params.human_dict)
        ])
        .set { ch_reference }
    
    Channel.value(file(params.arcashla_reference))
        .set { ch_arcashla_ref }
    
    // contamination check
    FASTP(ch_reads)
    BWA_MEM_RRNA(FASTP.out.reads, ch_rrna_reference)
    BWA_MEM_GLOBINRNA(FASTP.out.reads, ch_globinrna_reference)
    SAMTOOLS_VIEW_RRNA(BWA_MEM_RRNA.out.bam, "rrna")
    SAMTOOLS_VIEW_GLOBINRNA(BWA_MEM_GLOBINRNA.out.bam, "globinrna")
    SAMTOOLS_FLAGSTAT_RRNA(SAMTOOLS_VIEW_RRNA.out.view_bam, "rrna")
    SAMTOOLS_FLAGSTAT_GLOBINRNA(SAMTOOLS_VIEW_GLOBINRNA.out.view_bam, "globinrna")
    
    // Kallisto quant
    KALLISTO_QUANT(
        ch_reads, 
        params.kallisto_index)

    // STAR alignment
    STAR_ALIGNREADS(
        ch_reads, 
        params.star_index_151, 
        params.star_index_150, 
        params.star_index_120, 
        params.star_index_100, 
        params.star_index_75, 
        params.star_index_69)
    SAMTOOLS_INDEX(STAR_ALIGNREADS.out.star_bam)
    SAMBAMBA_MARKDUP(STAR_ALIGNREADS.out.star_bam)
    
    // HLA typing with arcasHLA
    ARCASHLA(
        SAMBAMBA_MARKDUP.out.marked_bam,
        ch_arcashla_ref
    )

    // Create counts by gene
    SUBREAD_FEATURECOUNTS(
        params.gencode_gtf_path, 
        SAMBAMBA_MARKDUP.out.marked_bam)

    // QC
    QUALIMAP_RNASEQ(
        SAMBAMBA_MARKDUP.out.marked_bam,
        params.gencode_gtf_path)

    RNASEQC(
        params.gencode_gtf_collapse, 
        SAMBAMBA_MARKDUP.out.marked_bam)

    // Create CRAM files
    SAMTOOLS_CRAM(
        ch_reference, 
        SAMBAMBA_MARKDUP.out.marked_bam)

    // Run IRFinder
    IRFINDER(
        params.ir_ref, 
        SAMBAMBA_MARKDUP.out.marked_bam) 

    // Calculate XBP1 coverage
    MOSDEPTH_BED(
        ch_reference, 
        params.xbp1_bed, 
        params.mt_bed, 
        SAMTOOLS_CRAM.out.rna_cram,
        SAMTOOLS_CRAM.out.rna_crai)

    // Create CDS bed file and run variant calling
    MOSDEPTH(SAMBAMBA_MARKDUP.out.marked_bam, ch_reference)
    BEDTOOLS_MERGE_INTERSECT(MOSDEPTH.out.per_base_bed, params.gencode_bed, params.min_coverage)

    DEEPVARIANT_RUNDEEPVARIANT(
        SAMBAMBA_MARKDUP.out.marked_bam,
        ch_reference, 
        BEDTOOLS_MERGE_INTERSECT.out.cds_bed, 
        ch_model)

    // CRAM to SAM 
    SAMTOOLS_BAM2SAM(
        ch_reference,
        SAMBAMBA_MARKDUP.out.marked_bam)
    
    // SAM to SJ
    BAM2SJ(SAMTOOLS_BAM2SAM.out.rna_sam)

    // Run MultiQC
    MULTIQC(
        FASTP.out.json,
        SAMTOOLS_FLAGSTAT_RRNA.out.flagstat_file,
        SAMTOOLS_FLAGSTAT_GLOBINRNA.out.flagstat_file,
        STAR_ALIGNREADS.out.final_log,
        SUBREAD_FEATURECOUNTS.out.gene_counts_summary,
        QUALIMAP_RNASEQ.out.results,
        RNASEQC.out.gene_tpm,
        KALLISTO_QUANT.out.run_info_json
    )
}
