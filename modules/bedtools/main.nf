process BEDTOOLS_MERGE_INTERSECT {
    tag "$meta.id"
    label 'bedtools_merge_intersect'

    input:
        tuple val(meta), path(per_base_bed)
        path gencode_bed
        val min_coverage

    output:
        tuple val(meta), path("*_cds_3x.bed"), emit: cds_bed
        path "bedtools_versions.yml"         , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        zcat ${per_base_bed} | egrep -v 'HLA|decoy|random|alt|chrUn|chrEBV' | awk -v OFS="\t" -v min_coverage=${min_coverage} '\$4 >= min_coverage { print }' > ${prefix}_filtered.bed
        bedtools merge -d 1 -c 4 -o mean -i ${prefix}_filtered.bed > ${prefix}_filtered_3x.bed
        bedtools intersect -a ${prefix}_filtered_3x.bed -b ${gencode_bed} > ${prefix}_cds_3x.bed

        cat <<-END_VERSIONS > bedtools_versions.yml
        "${task.process}":
            bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' ))
        END_VERSIONS
        """
}