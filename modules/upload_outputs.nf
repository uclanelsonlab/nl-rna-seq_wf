process upload_files {
    tag "Upload all the necessary output files"

    input:
    val sample_name
    val proband
    val tissue
    val output_bucket
    path flagstat_rrna
    path flagstat_rrna_v
    path flagstat_globinrna 
    path flagstat_globinrna_v
    tuple val(meta), path(reads_gene) 
    tuple val(meta), path(reads_gene_log)
    tuple val(meta), path(final_log)
    tuple val(meta), path(sj_tab)
    tuple val(meta), path(star_bam)
    tuple val(meta), path(star_log)
    path star_versions
    tuple val(meta), path(gene_counts)
    tuple val(meta), path(gene_counts_short)
    tuple val(meta), path(gene_counts_summary)
    tuple val(meta), path(subread_log)
    path subread_versions
    tuple val(meta), path(outrider_table)
    tuple val(meta), path(outrider_project)
    tuple val(meta), path(outrider_log)
    path qc_coverage
    path qc_exon_cv
    path qc_exon_reads
    path qc_gene_fragments
    path qc_gene_reads
    path qc_gene_tpm
    path qc_metrics
    path qc_log
    path qc_versions
    tuple val(meta), path(cram)
    tuple val(meta), path(crai)
    path cram_log
    path cram_versions
    tuple val(meta), path(irfinder_chr_coverage)
    tuple val(meta), path(irfinder_dir_val)
    tuple val(meta), path(irfinder_dir)
    tuple val(meta), path(irfinder_nondir_val)
    tuple val(meta), path(irfinder_nondir)
    tuple val(meta), path(irfinder_junc_count)
    tuple val(meta), path(irfinder_roi)
    tuple val(meta), path(irfinder_spans_point)
    path irfinder_log
    path irfinder_versions
    tuple val(meta), path(fastp_reads)
    tuple val(meta), path(fastp_json)
    tuple val(meta), path(fastp_html)
    tuple val(meta), path(fastp_log)
    path fastp_versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample_name}"

    """
    aws s3 cp ${flagstat_rrna} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${flagstat_rrna}
    aws s3 cp ${flagstat_rrna_v} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${flagstat_rrna_v}
    aws s3 cp ${flagstat_globinrna} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${flagstat_globinrna}
    aws s3 cp ${flagstat_globinrna_v} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${flagstat_globinrna_v}
    aws s3 cp ${reads_gene} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/alignment/${reads_gene}
    aws s3 cp ${reads_gene_log} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/alignment/${reads_gene_log}
    aws s3 cp ${final_log} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/alignment/${final_log}
    aws s3 cp ${sj_tab} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/alignment/${sj_tab}
    aws s3 cp ${star_log} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/alignment/${star_log}
    aws s3 cp ${star_bam} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/alignment/${star_bam}
    aws s3 cp ${star_versions} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/alignment/${star_versions}
    aws s3 cp ${gene_counts} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/counts/${gene_counts}
    aws s3 cp ${gene_counts_short} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/counts/${gene_counts_short}
    aws s3 cp ${gene_counts_summary} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/counts/${gene_counts_summary}
    aws s3 cp ${subread_log} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/counts/${subread_log}
    aws s3 cp ${outrider_table} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/outrider/${outrider_table}
    aws s3 cp ${outrider_project} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/outrider/${outrider_project}
    aws s3 cp ${outrider_log} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/outrider/${outrider_log}
    aws s3 cp ${qc_coverage} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${qc_coverage}
    aws s3 cp ${qc_exon_cv} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${qc_exon_cv}
    aws s3 cp ${qc_exon_reads} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${qc_exon_reads}
    aws s3 cp ${qc_gene_fragments} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${qc_gene_fragments}
    aws s3 cp ${qc_gene_reads} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${qc_gene_reads}
    aws s3 cp ${qc_gene_tpm} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${qc_gene_tpm}
    aws s3 cp ${qc_metrics} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${qc_metrics}
    aws s3 cp ${qc_log} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${qc_log}
    aws s3 cp ${qc_versions} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${qc_versions}
    aws s3 cp ${cram} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/alignment/${cram}
    aws s3 cp ${crai} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/alignment/${crai}
    aws s3 cp ${cram_log} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/alignment/${cram_log}
    aws s3 cp ${cram_versions} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/alignment/${cram_versions}
    aws s3 cp ${irfinder_chr_coverage} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/irfinder/${irfinder_chr_coverage}
    aws s3 cp ${irfinder_dir_val} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/irfinder/${irfinder_dir_val}
    aws s3 cp ${irfinder_dir} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/irfinder/${irfinder_dir}
    aws s3 cp ${irfinder_nondir_val} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/irfinder/${irfinder_nondir_val}
    aws s3 cp ${irfinder_nondir} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/irfinder/${irfinder_nondir}
    aws s3 cp ${irfinder_junc_count} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/irfinder/${irfinder_junc_count}
    aws s3 cp ${irfinder_roi} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/irfinder/${irfinder_roi}
    aws s3 cp ${irfinder_spans_point} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/irfinder/${irfinder_spans_point}
    aws s3 cp ${irfinder_log} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/irfinder/${irfinder_log}
    aws s3 cp ${irfinder_versions} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/irfinder/${irfinder_versions}
    
    aws s3 cp ${fastp_reads} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${fastp_reads}
    aws s3 cp ${fastp_json} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${fastp_json}
    aws s3 cp ${fastp_html} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${fastp_html}
    aws s3 cp ${fastp_log} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${fastp_log}
    aws s3 cp ${fastp_versions} ${output_bucket}/${proband}/${prefix}_${tissue}_rna/hg38/qc/${fastp_versions}
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