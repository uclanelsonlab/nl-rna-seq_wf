process MULTIQC {  
    label 'multiqc'

    input:
        path fastp_json
        path samtools_flagstat_rrna
        path samtools_flagstat_globinrna
        tuple val(meta), path(star_final_log)
        tuple val(meta2), path(featurecounts_gene_counts_summary)
        tuple val(meta3), path(qualimap_results)
        path rnaseqc_gene_tpm
        tuple val(meta4), path(kallisto_run_info_json)
    
    output:
        path "*multiqc_report.html" , emit: report
        path "multiqc_versions.yml" , emit: versions
    
    when:
        task.ext.when == null || task.ext.when
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def args = task.ext.args ?: ''
        """
        multiqc \\
            $args \\
            . \\
            --filename ${prefix}_multiqc_report.html \\
            --outdir . \\
            --verbose
        
        cat <<-END_VERSIONS > multiqc_versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """
    
    stub:
        """
        touch ${prefix}_multiqc_report.html
        mkdir ${prefix}_multiqc_data
        touch ${prefix}_multiqc_data/${prefix}_multiqc_general_stats.txt
        
        cat <<-END_VERSIONS > multiqc_versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """
}