process KALLISTO_QUANT {
    label "kallisto_quant"

    input:
        tuple val(meta), path(reads)
        path index

    output:
        tuple val(meta), path("*abundance.tsv") , emit: abundance_tsv
        tuple val(meta), path("*run_info.json") , emit: run_info_json
        tuple val(meta), path("*abundance.h5")  , emit: abundance_h5
        path "*versions.yml"                    , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        kallisto quant \\
            --threads ${task.cpus} \\
            --index ${index} \\
            --output-dir ${prefix} \\
            ${reads}
        mv ${prefix}/abundance.tsv ${prefix}_abundance.tsv
        mv ${prefix}/abundance.h5 ${prefix}_abundance.h5
        mv ${prefix}/run_info.json ${prefix}_run_info.json

        cat <<-END_VERSIONS > kallisto_versions.yml
        "${task.process}":
            kallisto: \$(echo \$(kallisto version) | sed "s/kallisto, version //g" )
        END_VERSIONS
        """

    stub:
        prefix = task.ext.prefix ?: "${meta.id}"

        """
        mkdir -p $prefix
        touch ${prefix}.log
        touch ${prefix}.run_info.json

        cat <<-END_VERSIONS > kallisto_versions.yml
        "${task.process}":
            kallisto: \$(echo \$(kallisto version) | sed "s/kallisto, version //g" )
        END_VERSIONS
        """
}