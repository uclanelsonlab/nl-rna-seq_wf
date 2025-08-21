process BWA_MEM {
    tag "$meta.id"
    label 'bwa_mem'

    input:
        tuple val(meta), path(reads)
        tuple val(meta2), path(reference)
        
    output:
        tuple val(meta), path("*.bam")  , emit: bam
        path  "*versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def fastq_r1 = reads[0]
        def fastq_r2 = reads[1]

        """
        tar -xvzf ${reference}
        fasta_file=\$(find . -maxdepth 2 -name "*fa")
        bwa mem -t $task.cpus \$fasta_file ${fastq_r1} ${fastq_r2} > ${prefix}.${meta2}.bam

        cat <<-END_VERSIONS > bwa_versions.yml
        "${task.process}":
            bwa: \$(echo \$(bwa --version 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}