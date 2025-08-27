process SAMTOOLS_VIEW {
    label "samtools_view"

    input:
        tuple val(meta), path(bwa_bam)
        val reference

    output:
        tuple val(meta), path("*.view.bam"), emit: view_bam
        tuple val(meta), path('*.log'),      emit: log
        path "*versions.yml",                emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"

        """
        samtools view \\
            -@ $task.cpus \\
            -bS \\
            -F 2304 \\
            -o ${prefix}_${reference}.view.bam \\
            ${bwa_bam} \\
            2> >(tee ${prefix}.view.log >&2)

        cat <<-END_VERSIONS > view_versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}

process SAMTOOLS_FLAGSTAT {
    label "samtools_flagstat"

    input:
        tuple val(meta), path(view_bam)
        val reference

    output:
        path "*.flagstat.txt",  emit: flagstat_file
        path "*versions.yml",   emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"

        """
        samtools flagstat \\
            ${view_bam} \\
            -O json > ${prefix}_${reference}.flagstat.txt

        cat <<-END_VERSIONS > flagstat_${reference}_versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}

process SAMTOOLS_INDEX {
    label "samtools_index"

    input:
        tuple val(meta), path(bam), path(bai)

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

process SAMTOOLS_CRAM {
    label "samtools_cram"

    input:
        path fasta
        path fai
        path dict
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*.hg38_rna.normal.cram"),        emit: rna_cram
        tuple val(meta), path("*.hg38_rna.normal.cram.crai"),   emit: rna_crai
        path '*.log',                                           emit: log
        path "*versions.yml",                                   emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"

        """
        samtools view -@ $task.cpus -T ${fasta} -C --output-fmt-option normal -o ${prefix}.hg38_rna.normal.cram ${bam} 2> >(tee ${prefix}.cram.log >&2)
        samtools index -@ $task.cpus ${prefix}.hg38_rna.normal.cram

        cat <<-END_VERSIONS > cram_versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}

process SAMTOOLS_BAM2SAM {
    label "samtools_bam2sam"
    
    input:
        tuple val(meta2), path(fasta), path(fai), path (dict)    
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("*.hg38_rna.normal.sam"), emit: rna_sam
        path '*.sam.log',                                emit: log
        path "*versions.yml",                            emit: versions
        tuple val(meta2), path("*.hg38_rna.normal.sam"), emit: rna_sam
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        samtools view -@ $task.cpus -h -o ${prefix}.hg38_rna.normal.sam ${bam} 2> >(tee ${prefix}.sam.log >&2)

        cat <<-END_VERSIONS > bam2sam_versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}

process SAMTOOLS_CRAM2BAM {
    label "samtools_cram2bam"

    input:
        tuple val(meta2), path(fasta), path(fai), path (dict)    
        tuple val(meta), path(cram)

    output:
        tuple val(meta), path("*.hg38_rna.normal.bam"),  emit: bam
        path '*.bam.log',                                emit: log
        path "*versions.yml",                            emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        samtools view -@ $task.cpus -T ${fasta} -o ${prefix}.hg38_rna.normal.bam ${cram} 2> >(tee ${prefix}.bam.log >&2)

        cat <<-END_VERSIONS > cram_versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}