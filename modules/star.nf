process check_star_reference {
    tag "Check the STAR index to download"

    input:
    tuple val(meta), path(reads)

    output:
    path "star_index", emit: star_index
    env sjdb_overhang
    
    script:
    """
    read_length=\$(zcat ${reads[0]} | head -n 1000 | wc --max-line-length)
    case \$read_length in
        151)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/star_index_151bp_gencode43.tar.gz"
            sjdb_overhang=150
            ;;
        150)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/star_index_150bp_gencode43.tar.gz"
            sjdb_overhang=149
            ;;
        120)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/star_index_120bp_gencode43.tar.gz"
            sjdb_overhang=119
            ;;
        100)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/star_index_100bp_gencode43.tar.gz"
            sjdb_overhang=99
            ;;
        75)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/star_index_75bp_gencode43.tar.gz"
            sjdb_overhang=74
            ;;
        69)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/star_index_69bp_gencode43.tar.gz"
            sjdb_overhang=68
            ;;
        *)
            star_index="s3://ucla-rare-diseases/UCLA-UDN/assets/reference/gencode43/GRCh38.p13/star_index_100bp_gencode43.tar.gz"
            sjdb_overhang=99
            ;;
    esac

    aws s3 cp \${star_index} .
    mkdir -p star_index/ && tar -xvf \$(basename \${star_index}) -C star_index/
    """
}

process star_alignreads {
    container 'quay.io/biocontainers/star:2.7.11b--h5ca1c30_4'
    cpus 32
    publishDir params.outdir, mode:'symlink'
    tag "STAR alignReads on $meta"   
    
    input:
    val reference
    val sjdb_overhang
    tuple val(meta), path(reads)
    tuple val(meta), path(json)
    tuple val(meta), path(html)
    tuple val(meta), path(log)
    path versions

    output:
    tuple val(meta), path("*.ReadsPerGene.out.tab.gz"),       emit: reads_gene
    tuple val(meta), path("*.ReadsPerGene.log.out"),          emit: reads_gene_log
    tuple val(meta), path("*.Log.final.out"),                 emit: final_log
    tuple val(meta), path("*.SJ.out.tab.gz"),                 emit: sj_tab
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"), emit: star_bam
    tuple val(meta), path('*.log'),                           emit: log
    path "*versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    STAR --runMode alignReads --runThreadN $task.cpus --genomeDir ${reference} --twopassMode Basic --sjdbOverhang ${sjdb_overhang} --readFilesIn ${reads[0]} ${reads[1]} --readFilesCommand zcat --outFileNamePrefix ${meta}. --alignSoftClipAtReferenceEnds Yes --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outBAMcompression -1 --outSAMunmapped Within --genomeLoad NoSharedMemory --outBAMsortingThreadN $task.cpus --outSAMattrRGline ID:rg1 SM:${prefix} PL:Illumina LB:${prefix} 2> >(tee ${prefix}.star.log >&2)
    
    echo -e "Gene\t${prefix}.Unstranded\t${prefix}.Antisense\t${prefix}.Sense" > tempgene_counts
    tail -n +5 ${prefix}.ReadsPerGene.out.tab >> tempgene_counts
    echo -e "Gene\t${prefix}.Unstranded\t${prefix}.Antisense\t${prefix}.Sense" > tempgene_stats
    head -n +4 ${prefix}.ReadsPerGene.out.tab >> tempgene_stats
    mv tempgene_counts ${prefix}.ReadsPerGene.out.tab
    mv tempgene_stats ${prefix}.ReadsPerGene.log.out
    gzip ${prefix}.SJ.out.tab ${prefix}.ReadsPerGene.out.tab

    cat <<-END_VERSIONS > star_versions.yml
    "${task.process}":
        STAR: \$(echo \$(STAR --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}
