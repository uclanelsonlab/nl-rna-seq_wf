process DEEPVARIANT_RUNDEEPVARIANT {
    label 'deepvariant_rundeepvariant'

    input:
        tuple val(meta), path(input), path(index)
        tuple val(meta2), path(fasta), path(fai), path (dict)
        tuple val(meta3), path(bed)
        tuple val(meta4), path(model_data), path(model_index), path(model_meta), path(model_info)

    output:
        tuple val(meta), path("*.vcf.gz")               , emit: vcf
        tuple val(meta), path("*.vcf.gz.tbi")           , emit: vcf_tbi
        tuple val(meta), path("*.g.vcf.gz")             , emit: gvcf
        tuple val(meta), path("*.g.vcf.gz.tbi")         , emit: gvcf_tbi
        tuple val(meta), path("*.visual_report.html")   , emit: report
        path "deepvariant_versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WES \\
        --customized_model=model.ckpt \\
        --vcf_stats_report=true \\
        --ref=${fasta} \\
        --reads=${input} \\
        --output_vcf=${prefix}.vcf.gz \\
        --output_gvcf=${prefix}.g.vcf.gz \\
        --regions=${bed} \\
        --intermediate_results_dir=tmp \\
        --make_examples_extra_args="split_skip_reads=true,channels=''" \\
        --num_shards=${task.cpus}
    
    cat <<-END_VERSIONS > deepvariant_versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        echo "" | gzip > ${prefix}.vcf.gz
        touch ${prefix}.vcf.gz.tbi
        echo "" | gzip > ${prefix}.g.vcf.gz
        touch ${prefix}.g.vcf.gz.tbi

        cat <<-END_VERSIONS > deepvariant_versions.yml
        "${task.process}":
            deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
        END_VERSIONS
        """
}