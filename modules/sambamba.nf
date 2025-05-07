process SAMBAMBA_MARKDUP {
    container "quay.io/biocontainers/sambamba:1.0.1--h6f6fda4_1"
    cpus 40
    tag "Running sambamba markdup on $meta"
    publishDir params.outdir, mode:'symlink'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.markdup.bam"), emit: marked_bam
    tuple val(meta), path('*.log'),         emit: log
    path("*versions.yml"),                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta}"

    """
    sambamba markdup -t $task.cpus --tmpdir ./ $bam ${prefix}.markdup.bam 2> >(tee ${prefix}.markdup.log >&2)

    cat <<-END_VERSIONS > sambamba_versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}