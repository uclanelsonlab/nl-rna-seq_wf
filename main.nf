nextflow.enable.dsl = 2

params.sample_name = 'UDN343710-S_fibroblast_rna'
params.cram = "s3://ucla-rare-diseases/UCLA-UDN/Analysis/UDN_cases/UDN018129/UDN343710-S_fibroblast_rna/hg38/alignment/UDN343710-S-fibroblast.hg38_rna.normal.cram"
params.cram_crai = "s3://ucla-rare-diseases/UCLA-UDN/Analysis/UDN_cases/UDN018129/UDN343710-S_fibroblast_rna/hg38/alignment/UDN343710-S-fibroblast.hg38_rna.normal.cram.crai"
params.human_fai = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.fa.fai"
params.human_dict = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.dict"
params.human_fasta = "s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/GRCh38.primary_assembly.genome.fa"
params.output_bucket = "s3://ucla-rare-diseases/UCLA-UDN/Analysis/UDN_cases/UDN018129/UDN343710-S_fibroblast_rna/hg38/bam2sj/"
params.outdir = "results"

log.info """\
    R N A - S E Q _ W F   P I P E L I N E
    ===================================
    sample_name         : ${params.sample_name}
    cram                : ${params.cram}
    bucket output       : ${params.output_bucket}
    """
    .stripIndent(true)

include { download_human_ref; DOWNLOAD_CRAM } from './modules/download_files.nf'
include { SAMTOOLS_CRAM2SAM } from './modules/samtools.nf'


workflow {
    // Download CRAM and reference files
    download_human_ref(params.human_fasta, params.human_fai, params.human_dict)
    DOWNLOAD_CRAM(params.sample_name, params.cram, params.cram_crai)
    // CRAM to BAM 
    SAMTOOLS_CRAM2SAM(download_human_ref.out.human_ref, DOWNLOAD_CRAM.out.data)
    // BAM to SAM
    // SAM to SJ
    // Uplaod SJ
    // cram_ch = samtools_cram(download_human_ref_ch, mark_dup_ch)

    }
