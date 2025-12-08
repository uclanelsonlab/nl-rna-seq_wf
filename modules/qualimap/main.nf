process QUALIMAP_RNASEQ {
    label 'qualimap_rnaseq'

    input:
        tuple val(meta), path(bam), path(bai)
        path gtf

    output:
        tuple val(meta), path("${prefix}"), emit: results
        path  "qualimap_versions.yml"     , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        prefix = task.ext.prefix ?: "${meta.id}"
        def memory = (task.memory.mega*0.8).intValue() + 'M'

        """
        unset DISPLAY
        mkdir -p tmp
        export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
        qualimap \\
            --java-mem-size=${memory} \\
            rnaseq \\
            -bam ${bam} \\
            -gtf ${gtf} \\
            -outdir ${prefix}

        cat <<-END_VERSIONS > qualimap_versions.yml
        "${task.process}":
            qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
        END_VERSIONS
        """

    stub:
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        mkdir ${prefix}

        cat <<-END_VERSIONS > qualimap_versions.yml
        "${task.process}":
            qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
        END_VERSIONS
        """
}