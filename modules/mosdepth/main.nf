process MOSDEPTH_BED {
    label "mosdepth_bed"

    input:
        tuple val(meta), path(fasta), path(fai), path (dict)  
        path xbp1_bed
        path mt_bed
        tuple val(meta), path(cram)
        tuple val(meta2), path(crai)


    output:
        tuple val(meta), path("*.mosdepth.global.dist.txt"), emit: global_dist
        tuple val(meta), path("*.mosdepth.region.dist.txt"), emit: region_dist
        tuple val(meta), path("*.mosdepth.summary.txt"),     emit: summary
        tuple val(meta), path("*.per-base.bed.gz"),          emit: perbase
        tuple val(meta), path("*.per-base.bed.gz.csi"),      emit: perbase_index
        tuple val(meta), path("*.regions.bed.gz"),           emit: regions_bed
        tuple val(meta), path("*.regions.bed.gz.csi"),       emit: regions_bed_index
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

process MOSDEPTH {
    tag "$meta.id"
    label 'mosdepth'

    input:
        tuple val(meta),  path(aln), path(index)
        tuple val(meta2), path(fasta), path(fai), path (dict)

    output:
        tuple val(meta), path('*.global.dist.txt')      , emit: global_dist_txt
        tuple val(meta), path('*.summary.txt')          , emit: summary_txt
        tuple val(meta), path('*.per-base.bed.gz')      , emit: per_base_bed
        tuple val(meta), path('*.per-base.bed.gz.csi')  , emit: per_base_csi
        path  "mosdepth_versions.yml"                   , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        def reference = fasta ? "--fasta ${fasta}" : ""

        """
        mosdepth \\
            --threads $task.cpus \\
            $reference \\
            $prefix \\
            $aln

        cat <<-END_VERSIONS > mosdepth_versions.yml
        "${task.process}":
            mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.global.dist.txt
        touch ${prefix}.region.dist.txt
        touch ${prefix}.summary.txt
        touch ${prefix}.per-base.d4
        echo "" | gzip > ${prefix}.per-base.bed.gz
        touch ${prefix}.per-base.bed.gz.csi
        echo "" | gzip > ${prefix}.regions.bed.gz
        touch ${prefix}.regions.bed.gz.csi
        echo "" | gzip > ${prefix}.quantized.bed.gz
        touch ${prefix}.quantized.bed.gz.csi
        echo "" | gzip > ${prefix}.thresholds.bed.gz
        touch ${prefix}.thresholds.bed.gz.csi

        cat <<-END_VERSIONS > mosdepth_versions.yml
        "${task.process}":
            mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
        END_VERSIONS
        """
}