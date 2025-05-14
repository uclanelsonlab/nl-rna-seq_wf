process upload_files {
    tag "Upload all the necessary output files"

    input:
    val family_id
    val bucket_dir
    val output_bucket
    // samtools
    path flagstat_rrna
    path flagstat_globinrna
    // star 
    path reads_gene
    path reads_gene_log
    path final_log
    path sj_tab
    path star_bam
    // bam2sj
    path sj_tab_gz
    // splice junctions
    path rare_junctions
    path all_rare_junctions
    // subreads
    path gene_counts
    path gene_counts_short
    path gene_counts_summary
    // cram
    path rna_cram
    path rna_crai

    script:
    """
    aws s3 cp ${flagstat_rrna} ${output_bucket}/${family_id}/${bucket_dir}/hg19/star_output/${flagstat_rrna}
    aws s3 cp ${flagstat_globinrna} ${output_bucket}/${family_id}/${bucket_dir}/hg19/star_output/${flagstat_globinrna}
    aws s3 cp ${reads_gene} ${output_bucket}/${family_id}/${bucket_dir}/hg19/star_output/${reads_gene}
    aws s3 cp ${reads_gene_log} ${output_bucket}/${family_id}/${bucket_dir}/hg19/star_output/${reads_gene_log}
    aws s3 cp ${final_log} ${output_bucket}/${family_id}/${bucket_dir}/hg19/star_output/${final_log}
    aws s3 cp ${sj_tab} ${output_bucket}/${family_id}/${bucket_dir}/hg19/star_output/${sj_tab}
    aws s3 cp ${star_bam} ${output_bucket}/${family_id}/${bucket_dir}/hg19/star_output/${star_bam}
    aws s3 cp ${sj_tab_gz} ${output_bucket}/${family_id}/${bucket_dir}/hg19/bam2sj/${sj_tab_gz}
    aws s3 cp ${all_rare_junctions} ${output_bucket}/${family_id}/${bucket_dir}/hg19/rare_junctions/${all_rare_junctions}
    aws s3 cp ${rare_junctions} ${output_bucket}/${family_id}/${bucket_dir}/hg19/rare_junctions/${rare_junctions}
    aws s3 cp ${gene_counts} ${output_bucket}/${family_id}/${bucket_dir}/hg19/featurecounts_pc_unstranded/${gene_counts}
    aws s3 cp ${gene_counts_short} ${output_bucket}/${family_id}/${bucket_dir}/hg19/featurecounts_pc_unstranded/${gene_counts_short}
    aws s3 cp ${gene_counts_summary} ${output_bucket}/${family_id}/${bucket_dir}/hg19/featurecounts_pc_unstranded/${gene_counts_summary}
    aws s3 cp ${rna_cram} ${output_bucket}/${family_id}/${bucket_dir}/hg19/star_output/${rna_cram}
    aws s3 cp ${rna_crai} ${output_bucket}/${family_id}/${bucket_dir}/hg19/star_output/${rna_crai}
    """
}