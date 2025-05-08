process DOWNLOAD_FASTQS {
    tag "Download ${meta} FASTQ files"

    input:
    val meta
    val library
    val fastq_bucket

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads

    script:
    """
    aws s3 cp ${fastq_bucket}/${library}/ . --exclude "*" --recursive --include "${meta}*"
    """
}

process DOWNLOAD_RNA_REF {
    tag "Download rna reference files"

    input:
    val rna_reference_path
    val type

    output:
    path "${type}_reference", emit: reference_dir

    script:
    """
    mkdir ${type}_reference
    aws s3 cp ${rna_reference_path}/${type}/ ${type}_reference/ --recursive
    """
}

process DOWNLOAD_IR_REF {
    tag "Download IRFinder reference files"

    input:
    val ir_ref

    output:
    path "ir_reference", emit: ir_reference_dir

    script:
    """
    mkdir ir_reference
    aws s3 cp ${ir_ref} ir_reference/ --recursive
    """
}

process DOWNLOAD_HUMAN_REF {
    tag "Download rna reference files"

    input:
    val fasta
    val fai
    val dict

    output:
    path "GRCh38.primary_assembly.genome.fa",       emit: human_fasta
    path "GRCh38.primary_assembly.genome.fa.fai",   emit: human_fai
    path "GRCh38.primary_assembly.genome.dict",     emit: human_dict

    script:
    """
    aws s3 cp ${fasta} .
    aws s3 cp ${fai} .
    aws s3 cp ${dict} .
    """
}

process DOWNLOAD_BED {
    label "download_bed"
    tag "Download ${bed}"

    input:
    val bed

    output:
    path "*.bed", emit: bed

    script:
    """
    aws s3 cp ${bed} .
    """
}

process DOWNLOAD_CRAM {
    label "download_cram"
    tag "Download ${cram}"


    input:
    val meta
    val cram

    output:
    tuple val(meta), path("*.cram"), emit: cram

    script:
    """
    aws s3 cp ${cram} .
    """
}