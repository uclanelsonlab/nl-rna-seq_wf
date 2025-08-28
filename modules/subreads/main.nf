process SUBREAD_FEATURECOUNTS {
    label "subread_featurecounts"

    input:
        path gencode_pc
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("*.gene_id.exon.ct"),             emit: gene_counts
        tuple val(meta), path("*.gene_id.exon.ct.short.txt"),   emit: gene_counts_short
        tuple val(meta), path("*.gene_id.exon.ct.summary"),     emit: gene_counts_summary
        tuple val(meta), path("*.log"),                         emit: log
        path("*versions.yml"),                                  emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        featureCounts -T $task.cpus -t exon -g gene_id -a ${gencode_pc} -o ${prefix}.gene_id.exon.ct -p -C --primary ${bam} 2> >(tee ${prefix}.featureCounts.log >&2)
        awk -F \$' ' 'BEGIN {OFS=FS} { print \$1, \$7 }' ${prefix}.gene_id.exon.ct > ${prefix}.gene_id.exon.ct.short.txt
        
        cat <<-END_VERSIONS > featureCounts_versions.yml
        "${task.process}":
            featureCounts: \$(echo \$(featureCounts -v 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}