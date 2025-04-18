process MOSDEPTH_BED {
    container "quay.io/biocontainers/mosdepth:0.3.10--h4e814b3_1"
    cpus 20
    tag "Mosdepth BED coverage"
    publishDir params.outdir, mode:'symlink'

    input:
    path fasta
    path fai
    path dict    
    path bed
    tuple val(meta), path(cram)
    tuple val(meta2), path(crai)
    tuple val(meta3), path(cram_log)
    tuple val(meta4), path(cram_version)


    output:
    path "*.mosdepth.global.dist.txt", emit: global_dist
    path "*.mosdepth.region.dist.txt", emit: region_dist
    path "*.mosdepth.summary.txt",     emit: summary
    path "*.per-base.bed.gz",          emit: perbase
    path "*.per-base.bed.gz.csi",      emit: perbase_index
    path "*.regions.bed.gz",           emit: regions_bed
    path "*.regions.bed.gz.csi",       emit: regions_bed_index
    path '*.log',                      emit: log
    path "*versions.yml",              emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ls ${crai}
    mosdepth -t $task.cpus --by ${bed} -f ${fasta} ${meta}_mosdepth ${cram} 2> >(tee ${sample_name}.mosdepth.log >&2)

    cat <<-END_VERSIONS > cram_versions.yml
    "${task.process}":
        mosdepth: \$(echo \$(mosdepth --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}