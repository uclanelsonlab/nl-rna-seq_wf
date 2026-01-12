process STAR_ALIGNREADS {
    label "star_alignreads"

    input:
        tuple val(meta), path(reads)
        path star_index_151
        path star_index_150
        path star_index_120
        path star_index_100
        path star_index_75
        path star_index_69
        path star_index
    output:
        tuple val(meta), path("*.Log.final.out")                , emit: final_log
        tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), emit: star_bam
        tuple val(meta), path('*.log')                          , emit: log
        path "*versions.yml"                                    , emit: versions
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mkdir -p star_index/ && tar -xvf ${star_index} -C star_index/

        STAR \\
            --runThreadN $task.cpus \\
            --genomeDir star_index/ \\
            --twopassMode Basic \\
            --readFilesIn ${reads[0]} ${reads[1]} \\
            --readFilesCommand zcat \\
            --outFileNamePrefix ${prefix}. \\
            --outSAMtype BAM SortedByCoordinate \\
            --outBAMcompression -1 \\
            --outSAMunmapped Within \\
            --genomeLoad NoSharedMemory \\
            --outBAMsortingThreadN $task.cpus \\
            --outSAMattrRGline ID:rg1 SM:${prefix} PL:Illumina LB:${prefix} 2> >(tee ${prefix}.star.log >&2)

        cat <<-END_VERSIONS > star_versions.yml
        "${task.process}":
            STAR: \$(echo \$(STAR --version 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}