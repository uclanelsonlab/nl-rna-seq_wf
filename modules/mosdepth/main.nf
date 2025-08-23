process MOSDEPTH_BED {
    label "mosdepth"

    input:
        path fasta
        path fai
        path dict    
        path xbp1_bed
        path mt_bed
        tuple val(meta), path(cram)
        tuple val(meta2), path(crai)


    output:
        path "*.mosdepth.global.dist.txt", emit: global_dist
        path "*.mosdepth.region.dist.txt", emit: region_dist
        path "*.mosdepth.summary.txt",     emit: summary
        path "*.per-base.bed.gz",          emit: perbase
        path "*.per-base.bed.gz.csi",      emit: perbase_index
        path "*.regions.bed.gz",           emit: regions_bed
        path "*.regions.bed.gz.csi",       emit: regions_bed_index
        path "mosdepth_versions.yml",      emit: versions
    
    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        ls ${crai}
        # XBP1
        mosdepth -t ${task.cpus} --by ${xbp1_bed} -f ${fasta} ${prefix}_mosdepth ${cram}
        mosdepth -t ${task.cpus} -n --by ${mt_bed} -f ${fasta} ${prefix}_MT_genes_hg38 ${cram}
        mosdepth -t ${task.cpus} -n --fast-mode --by 500 -f ${fasta} ${prefix}_genome_hg38 ${cram}

        cat <<-END_VERSIONS > mosdepth_versions.yml
        "${task.process}":
            mosdepth: \$(echo \$(mosdepth --version 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}