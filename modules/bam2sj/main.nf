process BAM2SJ {
    tag "Create SJ from BAM for $sample_name"
    publishDir params.outdir, mode:'symlink'

    input:
    tuple val(sample_name), path(sam_view)

    output:
    path "*.bam2SJ.out.tab.gz", emit: sj_tab_gz

    script:
    """
    bam2sj.pl ${sam_view} > ${sample_name}.bam2SJ.out.tab
    gzip ${sample_name}.bam2SJ.out.tab
    """
}