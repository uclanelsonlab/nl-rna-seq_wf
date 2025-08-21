process IRFINDER {
    label "irfinder"

    input:
        path ir_reference_gz
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*IRFinder-ChrCoverage.txt"),   emit: irfinder_chr_coverage  
        tuple val(meta), path("*IRFinder-IR-dir-val.txt"),    emit: irfinder_dir_val  
        tuple val(meta), path("*IRFinder-IR-dir.txt"),        emit: irfinder_dir  
        tuple val(meta), path("*IRFinder-IR-nondir-val.txt"), emit: irfinder_nondir_val 
        tuple val(meta), path("*IRFinder-IR-nondir.txt"),     emit: irfinder_nondir  
        tuple val(meta), path("*IRFinder-JuncCount.txt"),     emit: irfinder_junc_count  
        tuple val(meta), path("*IRFinder-ROI.txt"),           emit: irfinder_roi
        tuple val(meta), path("*RFinder-SpansPoint.txt"),     emit: irfinder_spans_point
        path("*.log"),                                        emit: log
        path("*versions.yml"),                                emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        tar -xzf ${ir_reference_gz}
        genome=\$(find . -name "genome.fa")
        genome_dir=\$(dirname \$genome)
        IRFinder -m BAM -t $task.cpus -r \$genome_dir -d irfinder_output/ ${bam} 2> >(tee ${prefix}.IRFinder.log >&2)
        
        for output in \$(ls irfinder_output/*txt); do file=\$(basename \$output); mv \$output ${prefix}_\$file; done
        
        cat <<-END_VERSIONS > IRFinder_versions.yml
        "${task.process}":
            IRFinder: \$(echo \$(IRFinder -v 2>&1) | awk '{print \$2}' )
        END_VERSIONS
        """
}