process ARCASHLA {
    tag "$meta.id"
    label 'arcashla'
    
    input:
        tuple val(meta), path(bam), path(bai)
        path(arcashla_reference)
    
    output:
        tuple val(meta), path("*genotype.json"), emit: genotype_json
        tuple val(meta), path("*genotype.log"), emit: genotype_log
        path "versions.yml", emit: versions
    
    script:
        def prefix = meta.id
        def ref_dir = "/usr/local/share/arcas-hla-0.6.0-2/dat"
        """
        # untar reference
        mkdir -p ${ref_dir}
        tar -xzf ${arcashla_reference} -C ${ref_dir}

        mkdir -p ${prefix}

        arcasHLA extract \\
            ${bam} \\
            -o ${prefix} \\
            -t ${task.cpus} \\
            -v

        arcasHLA genotype \\
            ${prefix}/*.extracted.1.fq.gz \\
            ${prefix}/*.extracted.2.fq.gz \\
            -g A,B,C,DPB1,DQB1,DQA1,DRB1 \\
            -o ${prefix} \\
            -t ${task.cpus} \\
            -v
        mv ${prefix}/* .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            arcasHLA: \$(arcasHLA version 2>&1 | grep -oE 'v[0-9.]+' | sed 's/^v//' || echo "unknown")
        END_VERSIONS
        """
    
    stub:
        def prefix = meta.id
        """
        mkdir -p ${prefix}
        touch ${prefix}/${prefix}.extracted.1.fq.gz
        touch ${prefix}/${prefix}.extracted.2.fq.gz
        touch ${prefix}/${prefix}.genotype.json
        touch ${prefix}/${prefix}.genotype.log
        touch versions.yml
        """
}
