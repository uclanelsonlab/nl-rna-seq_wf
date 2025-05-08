process UPLOAD_FILES {
    tag "Upload all the necessary output files"

    input:
    val family_id
    val bucket_dir
    val output_bucket
    tuple val(meta), path(fastp_json)
    tuple val(meta2), path(fastp_html)
    tuple val(meta3), path(fastp_log)
    path fastp_versions
    path rrna_flagstat_file
    path rrna_versions
    path globinrna_flagstat_file
    path globinrna_versions
    tuple val(meta4), path(star_reads_gene)
    tuple val(meta5), path(star_reads_gene_log)
    tuple val(meta6), path(star_final_log)
    tuple val(meta7), path(star_sj_tab)
    tuple val(meta8), path(star_log)
    path star_versions
    tuple val(meta9), path(sambamba_log)
    path sambamba_versions
    tuple val(meta10), path(featurecounts_gene_counts)
    tuple val(meta11), path(featurecounts_gene_counts_short)
    tuple val(meta12), path(featurecounts_gene_counts_summary)
    tuple val(meta13), path(featurecounts_log)
    path featurecounts_versions
    path rnaseqc_coverage
    path rnaseqc_exon_cv
    path rnaseqc_exon_reads
    path rnaseqc_gene_fragments
    path rnaseqc_gene_reads
    path rnaseqc_gene_tpm
    path rnaseqc_metrics
    path rnaseqc_log
    path rnaseqc_versions
    tuple val(meta14), path(bam2cram_rna_cram)
    tuple val(meta15), path(bam2cram_rna_crai)
    path bam2cram_log
    path bam2cram_versions
    tuple val(meta16), path(irfinder_chr_coverage)
    tuple val(meta17), path(irfinder_dir_val)
    tuple val(meta18), path(irfinder_dir)
    tuple val(meta19), path(irfinder_nondir_val)
    tuple val(meta20), path(irfinder_nondir)
    tuple val(meta21), path(irfinder_junc_count)
    tuple val(meta22), path(irfinder_roi)
    tuple val(meta23), path(irfinder_spans_point)
    path irfinder_log
    path irfinder_versions
    path bam2sam_log
    path bam2sam_versions

    script:
    """
    # qc
    aws s3 cp ${fastp_json} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${fastp_json}
    aws s3 cp ${fastp_html} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${fastp_html}
    aws s3 cp ${fastp_log} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${fastp_log}
    aws s3 cp ${fastp_versions} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${fastp_versions}
    aws s3 cp ${rrna_flagstat_file} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${rrna_flagstat_file}
    aws s3 cp ${rrna_versions} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${rrna_versions}
    aws s3 cp ${globinrna_flagstat_file} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${globinrna_flagstat_file}
    aws s3 cp ${globinrna_versions} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${globinrna_versions}
    aws s3 cp ${rnaseqc_coverage} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${rnaseqc_coverage}
    aws s3 cp ${rnaseqc_exon_cv} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${rnaseqc_exon_cv}
    aws s3 cp ${rnaseqc_exon_reads} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${rnaseqc_exon_reads}
    aws s3 cp ${rnaseqc_gene_fragments} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${rnaseqc_gene_fragments}
    aws s3 cp ${rnaseqc_gene_reads} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${rnaseqc_gene_reads}
    aws s3 cp ${rnaseqc_gene_tpm} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${rnaseqc_gene_tpm}
    aws s3 cp ${rnaseqc_metrics} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${rnaseqc_metrics}
    aws s3 cp ${rnaseqc_log} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${rnaseqc_log}
    aws s3 cp ${rnaseqc_versions} ${output_bucket}/${family_id}/${bucket_dir}/hg38/qc/${rnaseqc_versions}
    # alignment
    aws s3 cp ${bam2sam_log} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${bam2sam_log}
    aws s3 cp ${bam2sam_versions} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${bam2sam_versions}
    aws s3 cp ${star_reads_gene} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${star_reads_gene}
    aws s3 cp ${star_reads_gene_log} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${star_reads_gene_log}
    aws s3 cp ${star_final_log} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${star_final_log}
    aws s3 cp ${star_sj_tab} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${star_sj_tab}
    aws s3 cp ${star_log} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${star_log}
    aws s3 cp ${star_versions} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${star_versions}
    aws s3 cp ${sambamba_log} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${sambamba_log}
    aws s3 cp ${sambamba_versions} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${sambamba_versions}
    aws s3 cp ${bam2cram_rna_cram} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${bam2cram_rna_cram}
    aws s3 cp ${bam2cram_rna_crai} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${bam2cram_rna_crai}
    aws s3 cp ${bam2cram_log} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${bam2cram_log}
    aws s3 cp ${bam2cram_versions} ${output_bucket}/${family_id}/${bucket_dir}/hg38/alignment/${bam2cram_versions}
    # counts
    aws s3 cp ${featurecounts_gene_counts} ${output_bucket}/${family_id}/${bucket_dir}/hg38/counts/${featurecounts_gene_counts}
    aws s3 cp ${featurecounts_gene_counts_short} ${output_bucket}/${family_id}/${bucket_dir}/hg38/counts/${featurecounts_gene_counts_short}
    aws s3 cp ${featurecounts_gene_counts_summary} ${output_bucket}/${family_id}/${bucket_dir}/hg38/counts/${featurecounts_gene_counts_summary}
    aws s3 cp ${featurecounts_log} ${output_bucket}/${family_id}/${bucket_dir}/hg38/counts/${featurecounts_log}
    aws s3 cp ${featurecounts_versions} ${output_bucket}/${family_id}/${bucket_dir}/hg38/counts/${featurecounts_versions}
    # irfinder
    aws s3 cp ${irfinder_chr_coverage} ${output_bucket}/${family_id}/${bucket_dir}/hg38/irfinder/${irfinder_chr_coverage}
    aws s3 cp ${irfinder_dir_val} ${output_bucket}/${family_id}/${bucket_dir}/hg38/irfinder/${irfinder_dir_val}
    aws s3 cp ${irfinder_dir} ${output_bucket}/${family_id}/${bucket_dir}/hg38/irfinder/${irfinder_dir}
    aws s3 cp ${irfinder_nondir_val} ${output_bucket}/${family_id}/${bucket_dir}/hg38/irfinder/${irfinder_nondir_val}
    aws s3 cp ${irfinder_nondir} ${output_bucket}/${family_id}/${bucket_dir}/hg38/irfinder/${irfinder_nondir}
    aws s3 cp ${irfinder_junc_count} ${output_bucket}/${family_id}/${bucket_dir}/hg38/irfinder/${irfinder_junc_count}
    aws s3 cp ${irfinder_roi} ${output_bucket}/${family_id}/${bucket_dir}/hg38/irfinder/${irfinder_roi}
    aws s3 cp ${irfinder_spans_point} ${output_bucket}/${family_id}/${bucket_dir}/hg38/irfinder/${irfinder_spans_point}
    aws s3 cp ${irfinder_log} ${output_bucket}/${family_id}/${bucket_dir}/hg38/irfinder/${irfinder_log}
    aws s3 cp ${irfinder_versions} ${output_bucket}/${family_id}/${bucket_dir}/hg38/irfinder/${irfinder_versions}
    """
}

process UP_SJ {
    tag "Upload SJ file"

    input:
    path sj_tab_gz
    path global_dist
    path region_dist
    path summary
    path perbase
    path perbase_index
    path regions_bed
    path regions_bed_index
    val output_bucket

    script:

    """
    aws s3 cp ${sj_tab_gz} ${output_bucket}/bam2sj/${sj_tab_gz}
    aws s3 cp ${global_dist} ${output_bucket}/mosdepth/${global_dist}
    aws s3 cp ${region_dist} ${output_bucket}/mosdepth/${region_dist}
    aws s3 cp ${summary} ${output_bucket}/mosdepth/${summary}
    aws s3 cp ${perbase} ${output_bucket}/mosdepth/${perbase}
    aws s3 cp ${perbase_index} ${output_bucket}/mosdepth/${perbase_index}
    aws s3 cp ${regions_bed} ${output_bucket}/mosdepth/${regions_bed}
    aws s3 cp ${regions_bed_index} ${output_bucket}/mosdepth/${regions_bed_index}
    """
}