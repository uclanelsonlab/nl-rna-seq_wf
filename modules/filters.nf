process FILTER_FASTQ {
    tag "Filter $meta FASTQ files"
    // publishDir params.outdir, mode:'symlink'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.4kreads.fastq.gz') , emit: reads

    script:
    """
    zcat ${reads[0]} | head -n 4000000 | gzip > ${meta}_R1.4kreads.fastq.gz
    zcat ${reads[1]} | head -n 4000000 | gzip > ${meta}_R2.4kreads.fastq.gz
    """
}