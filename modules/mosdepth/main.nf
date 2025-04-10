process MOSDEPTH_BED {
    container "quay.io/biocontainers/mosdepth:0.3.3"
    cpus 40
    tag "Mosdepth BED coverage"
    publishDir params.outdir, mode:'symlink'

    input:
    tuple path(fasta), path(fai), path(dict)
    tuple val(sample_name), path(cram), path(crai)
    path bed

    output:
    path "*.mosdepth.global.dist.txt", emit: global_dist
    path "*.mosdepth.region.dist.txt", emit: region_dist
    path "*.mosdepth.summary.txt",     emit: summary
    path "*.per-base.bed.gz",          emit: per-base
    path "*.per-base.bed.gz.csi",      emit: per-base_index
    path "*.regions.bed.gz",           emit: regions_bed
    path "*.regions.bed.gz.csi",       emit: regions_bed_index
    path '*.log',                      emit: log
    path "*versions.yml",              emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mosdepth -t $task.cpus --by ${bed} -f ${fasta} ${sample_name}_mosdepth ${cram} 2> >(tee ${sample_name}.mosdepth.log >&2)

    cat <<-END_VERSIONS > cram_versions.yml
    "${task.process}":
        mosdepth: \$(echo \$(mosdepth --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}