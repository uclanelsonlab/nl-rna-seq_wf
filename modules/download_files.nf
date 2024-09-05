process download_fastqs {
    tag "Download ${meta} FASTQ files"

    input:
    val meta
    val tissue
    val library
    val fastq_bucket

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads

    script:
    """
    aws s3 cp ${fastq_bucket}/${library}/ . --exclude "*" --recursive --include "${meta}*${tissue}*"
    """
}

process download_rna_ref {
    tag "Download rna reference files"

    input:
    val rna_reference_path
    val type

    output:
    path "${type}_reference", emit: rrna_reference_dir

    script:
    """
    mkdir ${type}_reference
    aws s3 cp ${rna_reference_path}/${type}/ ${type}_reference/ --recursive
    """
}

process download_human_ref {
    tag "Download rna reference files"

    input:
    val fasta
    val fai
    val dict

    output:
    path "human_g1k_v37.fasta", emit: human_fasta
    path "human_g1k_v37.fasta.fai", emit: human_fai
    path "human_g1k_v37.dict", emit: human_dict

    script:
    """
    aws s3 cp ${fasta} .
    aws s3 cp ${fai} .
    aws s3 cp ${dict} .
    """
}