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

    output:
        tuple val(meta), path("*.ReadsPerGene.out.tab.gz")      , emit: reads_gene
        tuple val(meta), path("*.ReadsPerGene.log.out")         , emit: reads_gene_log
        tuple val(meta), path("*.Log.final.out")                , emit: final_log
        tuple val(meta), path("*.SJ.out.tab.gz")                , emit: sj_tab
        tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), emit: star_bam
        tuple val(meta), path('*.log')                          , emit: log
        path "*versions.yml"                                    , emit: versions
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        read_length=\$(zcat ${reads[0]} | head -n 1000 | wc -L)
        case \$read_length in
            151)
                star_index=${star_index_151}
                sjdb_overhang=150
                ;;
            150)
                star_index=${star_index_150}
                sjdb_overhang=149
                ;;
            120)
                star_index=${star_index_120}
                sjdb_overhang=119
                ;;
            100)
                star_index=${star_index_100}
                sjdb_overhang=99
                ;;
            75)
                star_index=${star_index_75}
                sjdb_overhang=74
                ;;
            69)
                star_index=${star_index_69}
                sjdb_overhang=68
                ;;
            *)
                star_index=${star_index_100}
                sjdb_overhang=99
                ;;
        esac

        mkdir -p star_index/ && tar -xvf \$star_index -C star_index/

        STAR \\
            --runMode alignReads \\
            --runThreadN $task.cpus \\
            --genomeDir star_index/ \\
            --twopassMode Basic \\
            --sjdbOverhang \$sjdb_overhang \\
            --readFilesIn ${reads[0]} ${reads[1]} \\
            --readFilesCommand zcat \\
            --outFileNamePrefix ${prefix}. \\
            --alignSoftClipAtReferenceEnds Yes \\
            --quantMode GeneCounts \\
            --outSAMtype BAM SortedByCoordinate \\
            --outBAMcompression -1 \\
            --outSAMunmapped Within \\
            --genomeLoad NoSharedMemory \\
            --outBAMsortingThreadN $task.cpus \\
            --outSAMattrRGline ID:rg1 SM:${prefix} PL:Illumina LB:${prefix} 2> >(tee ${prefix}.star.log >&2)
        
        echo -e "Gene\t${prefix}.Unstranded\t${prefix}.Antisense\t${prefix}.Sense" > tempgene_counts
        tail -n 5 ${prefix}.ReadsPerGene.out.tab >> tempgene_counts
        echo -e "Gene\t${prefix}.Unstranded\t${prefix}.Antisense\t${prefix}.Sense" > tempgene_stats
        head -n 4 ${prefix}.ReadsPerGene.out.tab >> tempgene_stats
        mv tempgene_counts ${prefix}.ReadsPerGene.out.tab
        mv tempgene_stats ${prefix}.ReadsPerGene.log.out
        gzip ${prefix}.SJ.out.tab ${prefix}.ReadsPerGene.out.tab

        cat <<-END_VERSIONS > star_versions.yml
        "${task.process}":
            STAR: \$(echo \$(STAR --version 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}