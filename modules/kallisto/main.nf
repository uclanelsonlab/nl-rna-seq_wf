process KALLISTO_QUANT {
    container "quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1"
    cpus 20
    tag "Kallisto quantification for ${meta}"
    publishDir params.outdir, mode:'symlink'

    input:
        tuple val(meta), path(reads)
        path index

    output:
        path "${meta}/abundance.h5",    emit: abundance_h5
        path "${meta}/abundance.tsv",   emit: abundance_h5
        path "${meta}/run_info.json",   emit: abundance_h5
        path "${meta}.kallisto.log",    emit: log
        path "kallisto_versions.yml",   emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
    """
    kallisto quant -t $task.cpus -i ${index} -o ${meta}  ${reads[0]} ${reads[1]} 2> >(tee ${meta}.kallisto.log >&2)

    cat <<-END_VERSIONS > kallisto_versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}