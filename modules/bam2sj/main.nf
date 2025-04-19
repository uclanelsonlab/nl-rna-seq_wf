process BAM2SJ {
    tag "Create SJ from BAM for $sample_name"
    publishDir params.outdir, mode:'symlink'

    input:
    tuple val(sample_name), path(sam_view)

    output:
    path "*.bam2SJ.out.tab.gz", emit: sj_tab_gz

    script:
    """
    SCRIPT_PATH=\$(find ${projectDir} -name "bam2sj.pl" | head -n 1)
        if [ -z "\$SCRIPT_PATH" ]; then
            echo "Erro: Script 'bam2sj.pl' não encontrado"
            exit 1
        fi
    perl \$SCRIPT_PATH ${sam_view} > ${sample_name}.bam2SJ.out.tab
    gzip ${sample_name}.bam2SJ.out.tab
    """
}