process BAM2SJ {
    label "bam2sj"

    input:
        tuple val(meta), path(sam_view)

    output:
        path "*.bam2SJ.out.tab.gz", emit: sj_tab_gz

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        SCRIPT_PATH=\$(find ${projectDir} -name "bam2sj.pl" | head -n 1)
            if [ -z "\$SCRIPT_PATH" ]; then
                echo "Erro: Script 'bam2sj.pl' não encontrado"
                exit 1
            fi
        perl \$SCRIPT_PATH ${sam_view} > ${prefix}.bam2SJ.out.tab
        gzip ${prefix}.bam2SJ.out.tab
        """
}