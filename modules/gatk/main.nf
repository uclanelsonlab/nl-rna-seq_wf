process GATK4_SPLITNCIGARREADS {
    tag "${meta}"
    label "SplitNCigarReads"
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    memory 32.GB
    cpus 16

    input:
        path fasta
        path fai
        path dict    
        tuple val(meta), path(bam)
        
    output:
        tuple val(meta), path('*.bam'), path ('*.bai'), emit: bam
        path  "versions.yml",                           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta}"
        def avail_mem = 20000
        if (!task.memory) {
            log.info '[GATK SplitNCigarReads] Available memory not known - defaulting to 20GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
            SplitNCigarReads \\
            -R ${fasta} \\
            -I ${bam} \\
            -O ${prefix}.cigar.bam \\
            --tmp-dir .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta}"

        """
        touch ${prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process GATK4_BASERECALIBRATOR {
    tag "$meta"
    label 'BaseRecalibrator'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    memory 32.GB
    cpus 16

    input:
        tuple val(meta), path(input), path(input_index)
        path fasta
        path fai
        path dict
        tuple val(meta5), path(dbsnp138), path(dbsnp138_index)
        tuple val(meta6), path(known_indels), path(known_indels_index)
        tuple val(meta7), path(indels_1000G), path(indels_1000G_index)
        tuple val(meta8), path(af_only_gnomad), path(af_only_gnomad_index)
        tuple val(meta9), path(small_exac_common_3), path(small_exac_common_3_index)

    output:
        tuple val(meta), path("*.table"), emit: table
        path "versions.yml"             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta}"

        def avail_mem = 20000
        if (!task.memory) {
            log.info '[GATK BaseRecalibrator] Available memory not known - defaulting to 20GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
            BaseRecalibrator  \\
            --input $input \\
            --output ${prefix}.bqsr.table \\
            --reference $fasta \\
            --known-sites ${dbsnp138} \\
            --known-sites ${known_indels} \\
            --known-sites ${indels_1000G} \\
            --known-sites ${af_only_gnomad} \\
            --known-sites ${small_exac_common_3} \\
            --tmp-dir .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta}"
        """
        touch ${prefix}.table
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process GATK4_APPLYBQSR {
    tag "$meta"
    label 'ApplyBQSR'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    memory 32.GB
    cpus 16

    input:
        tuple val(meta), path(input), path(input_index)
        tuple val(meta2), path(bqsr_table)
        path fasta
        path fai
        path dict

    output:
        tuple val(meta), path("*.bam"), path("*.bai") , emit: bam
        path "versions.yml"            , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta}"
        def avail_mem = 20000
        if (!task.memory) {
            log.info '[GATK ApplyBQSR] Available memory not known - defaulting to 20GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
            ApplyBQSR \\
            --input $input \\
            --output ${prefix}.applybqsr.${input.getExtension()} \\
            --reference $fasta \\
            --bqsr-recal-file $bqsr_table \\
            --tmp-dir . --create-output-bam-index true

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def prefix = task.ext.prefix ?: "${meta}"
        def input_extension = input.getExtension()
        def output_extension = input_extension == 'bam' ? 'bam' : 'cram'
        """
        touch ${prefix}.${output_extension}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}

process GATK4_HAPLOTYPECALLER {
    tag "$meta"
    label 'process_low'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir params.outdir, mode:'symlink'
    memory 32.GB
    cpus 16

    input:
        tuple val(meta),  path(input), path(input_index)
        path fasta
        path fai
        path dict
        tuple val(meta5), path(dbsnp), path(dbsnp_tbi)

    output:
        tuple val(meta), path("*.vcf.gz")       , emit: vcf
        tuple val(meta), path("*.tbi")          , optional:true, emit: tbi
        path "versions.yml"                     , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta}"
        def dbsnp_command = dbsnp ? "--dbsnp $dbsnp" : ""

        def avail_mem = 20000
        if (!task.memory) {
            log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 20GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
            HaplotypeCaller \\
            --input $input \\
            --output ${prefix}.hc.vcf.gz \\
            --reference $fasta \\
            --native-pair-hmm-threads ${task.cpus} \\
            $dbsnp_command \\
            --tmp-dir .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """

    stub:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta}"
        def bamout_command = args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""

        def stub_realigned_bam = bamout_command ? "touch ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""
        """
        touch ${prefix}.vcf.gz
        touch ${prefix}.vcf.gz.tbi
        ${stub_realigned_bam}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}