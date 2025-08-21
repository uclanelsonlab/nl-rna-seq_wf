process SAMBAMBA_MARKDUP {
    label "sambamba_markdup"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*.markdup.bam")  , emit: marked_bam
        tuple val(meta), path('*.log')          , emit: log
        path("*versions.yml")                   , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        sambamba markdup -t $task.cpus --tmpdir ./ $bam ${prefix}.markdup.bam 2> >(tee ${prefix}.markdup.log >&2)

        cat <<-END_VERSIONS > sambamba_versions.yml
        "${task.process}":
            sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}