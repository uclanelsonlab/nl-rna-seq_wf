# rna-seq_wf
## RNA-seq pipeline documentation
- Make sure to have installed:
    - [Docker](https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository)
    ```bash
    apt-get install ca-certificates curl
    install -m 0755 -d /etc/apt/keyrings
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
    chmod a+r /etc/apt/keyrings/docker.asc
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null
    apt-get update
    apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
    docker ps
    ```
    - Conda
    ```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    ...
    exec bash
    ```
    - Add conda-forge as priority
    ```bash
    conda config --add channels conda-forge
    ```
- Requirements
    - Install nextflow
    - Install docker
    - Make sure your system has at least 36 nodes/CPUs

- Clone the repo
```bash
git clone https://github.com/uclanelsonlab/nl-rna-seq_wf.git
```

- Run the pipeline for your sample, it expects the FASTQ files to be at `s3://ucla-rare-diseases/UCLA-UDN/rnaseq/fastq` to download
```bash
cd nl-rna-seq_wf/
chmod u+x -R modules/
nextflow run main.nf --fastq_r1 s3://ucla-rare-diseases/UCLA-UDN/rnaseq/fastq/BG-2024-10-15/UDN748413-2931652-MGML0088-FBR1-R1_001.fast
q.gz --fastq_r2 s3://ucla-rare-diseases/UCLA-UDN/rnaseq/fastq/BG-2024-10-15/UDN748413-2931652-MGML0088-FBR1-R2_001.fastq.gz --prefix UDN748413-2931652-MGML0088-FBR1 --family_id UDN748413 --bucket_dir UDN748413-P_fibroblast_rnaseq
```
> If something happens and you need to rerun the pipeline from where it stopped, remember you can use `-resume` for those cases.

- If you want to add variant calling to your outputs just add the parameter: `--varcall true`

- If you want to add OUTRIDER to your outputs just add the parameters: `--outrider true --tissue ${SAMPLE_TISSUE}`

- Example of for loop to run samples over a samplesheet:
```bash
while IFS=, read fastq1 fastq2 prefix output_bucket family_id output_directory; do
	nextflow run main.nf --fastq_r1 ${fastq1} --fastq_r2 ${fastq2} --prefix ${prefix} --family_id ${family_id} --bucket_dir ${output_directory} --output_bucket ${output_bucket}
    rm -r work/ results/
done < "2025-05-09_rna-seq_processing-request.csv"
```

- Check if you have your outputs on S3:
```bash
aws s3 ls s3://ucla-rare-diseases/UCLA-UDN/Analysis/UDN_cases/ --recursive | grep UDN748413
```

- The expected outputs:
    - flagstat_rrna: Ribosomal contamination stats for human_rRNA_strict.fasta (`*.rrna.flagstat.txt`)
    - flagstat_globinrna: Ribosomal contamination stats for human_globinRNA.fa (`*.globinrna.flagstat.txt`)
    - reads_gene: STAR's tab output (`*.ReadsPerGene.out.tab.gz`)
    - reads_gene_log: STAR's log output (`*.ReadsPerGene.log.out`)
    - final_log: STAR's final log (`*.Log.final.out`)
    - sj_tab: STAR's junctions output (`*.SJ.out.tab.gz`)
    - sj_tab_gz: Home made script output of junctions recreated from STAR's output (`*.bam2SJ.out.tab.gz`)
    - all_rare_junctions: Home made script output of all junctions in TSV (`*_rare_junctions_all.tsv`) 
    - rare_junctions: Home made script output of junctions in spreadsheet (`*_rare_junctions_filtered.xlsx`)
    - gene_counts: Subreads featureCounts output (`*.gene_id.exon.ct`)
    - gene_counts_short: Subreads featureCounts output (`*.gene_id.exon.ct.short.txt`)
    - gene_counts_summary: Subreads featureCounts output (`*.gene_id.exon.ct.summary`)
    - rna_cram: CRAM aligned file output from STAR (`*.hg19_rna.normal.cram`)
    - rna_crai: CRAM's index aligned file output from STAR (`*.hg19_rna.normal.cram.crai`)
