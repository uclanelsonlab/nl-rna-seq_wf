process DOWNLOAD_MASTER_FEATURECOUNTS {
    tag "Download featureCounts master file"

    input:
    val featurecounts_master_path

    output:
    path("*.tsv"), emit: featurecounts_master

    script:
    """
    aws s3 cp ${featurecounts_master_path} .
    """
}

process ADD_SAMPLE_COUNTS_MASTER {
    tag "Add sample's counts to featureCounts master file"
    conda 'pandas'
    publishDir params.outdir, mode:'symlink'

    input:
    path featurecounts_master_path
    tuple val(meta), path(gene_counts_short)

    output:
    tuple val(meta), path("*_updated.tsv"), emit: featurecounts_updated

    script:
    """
    add1_count.py -m ${featurecounts_master_path} -n ${gene_counts_short} -o .
    """
}

process RUN_OUTRIDER {
    tag "Rscript to run outrider"
    container 'gvcn/outrider:1.20.1'
    publishDir params.outdir, mode:'symlink'

    input:
    tuple val(meta), path(featurecounts_updated)
    val tissue

    output:
    tuple val(meta), path("*_results.tsv"),    emit: outrider_table
    tuple val(meta), path("*.rds"),             emit: outrider_project
    tuple val(meta), path("*.log"),             emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta}"
    """
    run_outrider.R -f ${featurecounts_updated} -t ${tissue} -o ${prefix}_${tissue}_results.tsv 2> >(tee ${prefix}.outrider.log >&2)
    """
}
