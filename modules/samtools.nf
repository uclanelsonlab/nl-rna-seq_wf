process samtools_view {
    container "quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0"
    cpus 20
    tag "Samtools view on $meta"

    input:
    tuple val(meta), path(bwa_bam)
    tuple val(meta2), path(log)
    path versions
    val reference

    output:
    tuple val(meta), path("*.view.bam"), emit: view_bam
    tuple val(meta), path('*.log'),      emit: log
    path "*versions.yml",                emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    samtools view -@ $task.cpus -bS -F 2304 -o ${prefix}_${reference}.view.bam ${bwa_bam} 2> >(tee ${prefix}.view.log >&2)

    cat <<-END_VERSIONS > view_versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}

process samtools_flagstat {
    container "quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0"
    tag "Samtools flagstat on $meta"
    publishDir params.outdir, mode:'symlink'

    input:
    tuple val(meta), path(view_bam)
    tuple val(meta2), path(log)
    path versions
    val reference

    output:
    path "*.flagstat.txt",  emit: flagstat_file
    path "*versions.yml",   emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    samtools flagstat ${view_bam} -O json > ${prefix}_${reference}.flagstat.txt

    cat <<-END_VERSIONS > flagstat_${reference}_versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}

process samtools_index {
    container "quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0"
    cpus 20
    tag "Samtools index on $bam"

    input:
    tuple val(meta), path(reads_gene)
    tuple val(meta2), path(reads_gene_log)
    tuple val(meta3), path(final_log)
    tuple val(meta4), path(sj_tab)
    tuple val(meta5), path(bam)
    tuple val(meta6), path(log)
    path versions

    output:
    tuple val(meta), path("*.bai"), emit: bam_index
    path "*versions.yml",           emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools index -@ $task.cpus ${bam}

    cat <<-END_VERSIONS > index_versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}

process samtools_cram {
    container "quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0"
    cpus 40
    tag "Samtools view on $meta BAM to CRAM"
    publishDir params.outdir, mode:'symlink'

    input:
    path fasta
    path fai
    path dict
    tuple val(meta), path(bam)
    tuple val(meta2), path(log)
    path versions

    output:
    tuple val(meta), path("*.hg38_rna.normal.cram"),        emit: rna_cram
    tuple val(meta), path("*.hg38_rna.normal.cram.crai"),   emit: rna_crai
    path '*.log',                                           emit: log
    path "*versions.yml",                                   emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    samtools view -@ $task.cpus -T ${fasta} -C --output-fmt-option normal -o ${prefix}.hg38_rna.normal.cram ${bam} 2> >(tee ${prefix}.cram.log >&2)
    samtools index -@ $task.cpus ${prefix}.hg38_rna.normal.cram

    cat <<-END_VERSIONS > cram_versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}

process SAMTOOLS_CRAM2SAM {
    container "quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0"
    cpus 40
    tag "Samtools view on $meta CRAM to SAM"
    publishDir params.outdir, mode:'symlink'

    input:
    path fasta
    path fai
    path dict    
    tuple val(meta), path(cram)
    tuple val(meta2), path(crai)
    tuple val(meta3), path(cram_log)
    tuple val(meta4), path(cram_version)

    output:
    tuple val(meta), path("*.hg38_rna.normal.sam"),       emit: rna_sam
    path '*.log',                                                emit: log
    path "*versions.yml",                                        emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools view -@ $task.cpus -T ${fasta} --input-fmt-option normal -o ${meta}.hg38_rna.normal.sam ${cram} 2> >(tee ${meta}.cram.log >&2)

    cat <<-END_VERSIONS > cram_versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}