process download_master_featureCounts {
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

process add_sample_counts_master {
    tag "Add sample's counts to featureCounts master file"
    conda 'conda-forge::python=3.8 conda-forge::pandas=2.2.2 conda-forge::openpyxl'

    input:
    val featurecounts_master_path
    val sample_featurecounts

    output:
    path("*_updated.tsv"), emit: featurecounts_updated

    script:
    """
    python add1_count.py -m ${featurecounts_master_path} -n ${sample_featurecounts} -o .
    """
}