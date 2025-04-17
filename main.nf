nextflow.enable.dsl = 2

params.sample_name = 'UDN343710-S_fibroblast_rna'
params.cram = "s3://ucla-rare-diseases/UCLA-UDN/Analysis/UDN_cases/UDN018129/UDN343710-S_fibroblast_rna/hg38/alignment/UDN343710-S-fibroblast.hg38_rna.normal.cram"
params.cram_crai = "s3://ucla-rare-diseases/UCLA-UDN/Analysis/UDN_cases/UDN018129/UDN343710-S_fibroblast_rna/hg38/alignment/UDN343710-S-fibroblast.hg38_rna.normal.cram.crai"
params.human_fai = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.fa.fai"
params.human_dict = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.dict"
params.human_fasta = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.fa"
params.output_bucket = "s3://ucla-rare-diseases/UCLA-UDN/Analysis/UDN_cases/UDN018129/UDN343710-S_fibroblast_rna/hg38"
params.outdir = "results"
params.bed = "s3://ucla-rare-diseases/UCLA-UDN/gcarvalho_test/rnaseq/hg38/XBP1-NM_001079539-hg38.bed"

log.info """\
    R N A - S E Q _ W F   P I P E L I N E
    ===================================
    sample_name         : ${params.sample_name}
    cram                : ${params.cram}
    cram_crai           : ${params.cram_crai}
    bucket output       : ${params.output_bucket}
    """
    .stripIndent(true)

include { download_human_ref; DOWNLOAD_CRAM; DOWNLOAD_BED } from './modules/download_files.nf'
include { SAMTOOLS_CRAM2SAM } from './modules/samtools.nf'
include { BAM2SJ } from './modules/bam2sj/main.nf'
include { UP_SJ } from './modules/upload_outputs.nf'
include { MOSDEPTH_BED } from './modules/mosdepth/main.nf'


workflow {
    // Download CRAM and reference files
    download_human_ref(params.human_fasta, params.human_fai, params.human_dict)
    DOWNLOAD_CRAM(params.sample_name, params.cram, params.cram_crai)
    DOWNLOAD_BED(params.bed)
    // Calculate XBP1 coverage
    MOSDEPTH_BED(download_human_ref.out.human_ref, DOWNLOAD_CRAM.out.data, DOWNLOAD_BED.out.bed)
    // CRAM to SAM 
    SAMTOOLS_CRAM2SAM(download_human_ref.out.human_ref, DOWNLOAD_CRAM.out.data)
    // SAM to SJ
    BAM2SJ(SAMTOOLS_CRAM2SAM.out.rna_sam)
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
