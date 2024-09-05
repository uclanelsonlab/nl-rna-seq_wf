process bam2sj {
    tag "Create SJ from BAM for $meta"
    publishDir params.outdir, mode:'symlink'

    input:
    val meta
    val tissue
    path sam_view

    output:
    path "*.bam2SJ.out.tab.gz", emit: sj_tab_gz

    script:
    """
    bam2sj.pl ${sam_view} > ${meta}-${tissue}.bam2SJ.out.tab
    gzip ${meta}-${tissue}.bam2SJ.out.tab
    """
}