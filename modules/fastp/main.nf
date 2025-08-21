process FASTP {
    tag "$meta.id"
    label "fastp"

    input:
        tuple val(meta), path(reads)
    
    output:
        path "*.html"                               , emit: html
        path "*.json"                               , emit: json
        path "versions.yml"                         , emit: versions
        tuple val(meta), path('*.4kreads.fastq.gz') , emit: reads
    
    when:
        task.ext.when == null || task.ext.when
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fastq_r1 = reads[0]
        def fastq_r2 = reads[1]
        """
        fastp \\
            -i ${fastq_r1} \\
            -I ${fastq_r2} \\
            --thread ${task.cpus} \\
            -h ${prefix}_fastp.html \\
            -j ${prefix}_fastp.json \\
            --detect_adapter_for_pe
        
        zcat ${fastq_r1} | head -n 4000000 | gzip > ${prefix}_R1.4kreads.fastq.gz
        zcat ${fastq_r2} | head -n 4000000 | gzip > ${prefix}_R2.4kreads.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/^fastp //')
        END_VERSIONS
        """
    
    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}_fastp.html
        touch ${prefix}_fastp.json
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/^fastp //')
        END_VERSIONS
        """
}