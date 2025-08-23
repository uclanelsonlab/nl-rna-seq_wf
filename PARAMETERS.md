# Pipeline Parameters Reference

Complete reference for all parameters used in the RNA-seq workflow pipeline.

## 📋 Required Parameters

These parameters must be provided for the pipeline to run successfully.

### Input Files

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `--fastq_r1` | String | Forward FASTQ file path (R1) | `/path/to/sample_R1.fastq.gz` |
| `--fastq_r2` | String | Reverse FASTQ file path (R2) | `/path/to/sample_R2.fastq.gz` |
| `--sample_name` | String | Sample identifier/name | `SAMPLE001` |

### Reference Files

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `--rrna_reference` | String | rRNA reference genome path | `/refs/rrna.fa` |
| `--globinrna_reference` | String | Globin RNA reference path | `/refs/globin.fa` |
| `--gencode_gtf_path` | String | Gencode GTF annotation path | `/refs/gencode.v43.gtf` |
| `--gencode_gtf_collapse` | String | Collapsed GTF annotation path | `/refs/gencode.v43.collapse.gtf` |
| `--human_fasta` | String | Human reference genome FASTA | `/refs/GRCh38.fa` |
| `--human_fai` | String | Human genome FASTA index | `/refs/GRCh38.fa.fai` |
| `--human_dict` | String | Human genome dictionary | `/refs/GRCh38.dict` |

### Analysis Tools

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `--ir_ref` | String | IRFinder reference directory | `/refs/irfinder_ref/` |
| `--xbp1_bed` | String | XBP1 BED file for coverage | `/refs/xbp1.bed` |
| `--mt_bed` | String | Mitochondrial BED file | `/refs/mt.bed` |

### STAR Indices

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `--star_index_151` | String | STAR index for 151bp reads | `/refs/star_151/` |
| `--star_index_150` | String | STAR index for 150bp reads | `/refs/star_150/` |
| `--star_index_120` | String | STAR index for 120bp reads | `/refs/star_120/` |
| `--star_index_100` | String | STAR index for 100bp reads | `/refs/star_100/` |
| `--star_index_75` | String | STAR index for 75bp reads | `/refs/star_75/` |
| `--star_index_69` | String | STAR index for 69bp reads | `/refs/star_69/` |

### Quantification

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `--kallisto_index` | String | Kallisto transcript index | `/refs/kallisto.idx` |

## 🔧 Optional Parameters

These parameters have default values but can be customized.

### Docker Images

| Parameter | Default | Description | Example |
|-----------|---------|-------------|---------|
| `--bwa_docker` | `bwa:latest` | BWA Docker image | `your-registry/bwa:v0.7.17` |
| `--star_docker` | `star:latest` | STAR Docker image | `your-registry/star:2.7.10a` |
| `--samtools_docker` | `samtools:latest` | Samtools Docker image | `your-registry/samtools:1.17` |
| `--fastp_docker` | `fastp:latest` | Fastp Docker image | `your-registry/fastp:0.23.4` |
| `--sambamba_docker` | `sambamba:latest` | Sambamba Docker image | `your-registry/sambamba:0.8.2` |
| `--subread_docker` | `subread:latest` | Subread Docker image | `your-registry/subread:2.0.6` |
| `--irfinder_docker` | `irfinder:latest` | IRFinder Docker image | `your-registry/irfinder:1.3.1` |
| `--mosdepth_docker` | `mosdepth:latest` | Mosdepth Docker image | `your-registry/mosdepth:0.3.3` |
| `--rnaseqc_docker` | `rnaseqc:latest` | RNA-SeQC Docker image | `your-registry/rnaseqc:2.4.2` |
| `--kallisto_docker` | `kallisto:latest` | Kallisto Docker image | `your-registry/kallisto:0.48.0` |

### Output Configuration

| Parameter | Default | Description | Example |
|-----------|---------|-------------|---------|
| `--outdir` | `/mnt/workflow/pubdir` | Output directory | `/path/to/results` |
| `--publish_dir_mode` | `copy` | Output file mode | `copy`, `symlink`, `rellink` |

## 📝 Parameter Examples

### Basic Run
```bash
nextflow run main.nf \
    --fastq_r1 /data/sample_R1.fastq.gz \
    --fastq_r2 /data/sample_R2.fastq.gz \
    --sample_name UDN748413 \
    --rrna_reference /refs/rrna.fa \
    --globinrna_reference /refs/globin.fa \
    --gencode_gtf_path /refs/gencode.v43.gtf \
    --human_fasta /refs/GRCh38.fa \
    --star_index_151 /refs/star_151 \
    --kallisto_index /refs/kallisto.idx
```

### Complete Run with All References
```bash
nextflow run main.nf \
    --fastq_r1 /data/sample_R1.fastq.gz \
    --fastq_r2 /data/sample_R2.fastq.gz \
    --sample_name UDN748413 \
    --rrna_reference /refs/rrna.fa \
    --globinrna_reference /refs/globin.fa \
    --gencode_gtf_path /refs/gencode.v43.gtf \
    --gencode_gtf_collapse /refs/gencode.v43.collapse.gtf \
    --human_fasta /refs/GRCh38.fa \
    --human_fai /refs/GRCh38.fa.fai \
    --human_dict /refs/GRCh38.dict \
    --ir_ref /refs/irfinder_ref \
    --xbp1_bed /refs/xbp1.bed \
    --mt_bed /refs/mt.bed \
    --star_index_151 /refs/star_151 \
    --star_index_150 /refs/star_150 \
    --star_index_120 /refs/star_120 \
    --star_index_100 /refs/star_100 \
    --star_index_75 /refs/star_75 \
    --star_index_69 /refs/star_69 \
    --kallisto_index /refs/kallisto.idx
```

### Custom Docker Images
```bash
nextflow run main.nf \
    --fastq_r1 /data/sample_R1.fastq.gz \
    --fastq_r2 /data/sample_R2.fastq.gz \
    --sample_name SAMPLE001 \
    --rrna_reference /refs/rrna.fa \
    --globinrna_reference /refs/globin.fa \
    --gencode_gtf_path /refs/gencode.v43.gtf \
    --human_fasta /refs/GRCh38.fa \
    --star_index_151 /refs/star_151 \
    --kallisto_index /refs/kallisto.idx \
    --bwa_docker "my-registry/bwa:v0.7.17" \
    --star_docker "my-registry/star:2.7.10a" \
    --fastp_docker "my-registry/fastp:0.23.4"
```

### Custom Output Directory
```bash
nextflow run main.nf \
    --fastq_r1 /data/sample_R1.fastq.gz \
    --fastq_r2 /data/sample_R2.fastq.gz \
    --sample_name SAMPLE001 \
    --rrna_reference /refs/rrna.fa \
    --globinrna_reference /refs/globin.fa \
    --gencode_gtf_path /refs/gencode.v43.gtf \
    --human_fasta /refs/GRCh38.fa \
    --star_index_151 /refs/star_151 \
    --kallisto_index /refs/kallisto.idx \
    --outdir /path/to/custom/results \
    --publish_dir_mode symlink
```

## 🔍 Parameter Validation

The pipeline automatically validates required parameters:

```groovy
// Parameter validation in main.nf
if (!params.fastq_r1) {
    error "Missing required parameter: --fastq_r1"
}
if (!params.fastq_r2) {
    error "Missing required parameter: --fastq_r2"
}
if (!params.rrna_reference) {
    error "Missing required parameter: --rrna_reference"
}
if (!params.globinrna_reference) {
    error "Missing required parameter: --globinrna_reference"
}
```

## 📁 File Format Requirements

### Input Files
- **FASTQ**: Compressed (.gz) or uncompressed
- **Reference FASTA**: Standard FASTA format with .fai index
- **GTF**: Gencode annotation format
- **BED**: Standard BED format
- **STAR Index**: Directory containing STAR index files
- **Kallisto Index**: .idx file format

### Reference File Organization
```
/refs/
├── rrna.fa                    # rRNA reference
├── globin.fa                  # Globin RNA reference
├── GRCh38.fa                  # Human genome
├── GRCh38.fa.fai             # Genome index
├── GRCh38.dict               # Genome dictionary
├── gencode.v43.gtf           # Gene annotation
├── gencode.v43.collapse.gtf  # Collapsed annotation
├── star_151/                 # STAR index 151bp
├── star_150/                 # STAR index 150bp
├── star_120/                 # STAR index 120bp
├── star_100/                 # STAR index 100bp
├── star_75/                  # STAR index 75bp
├── star_69/                  # STAR index 69bp
├── kallisto.idx              # Kallisto index
├── irfinder_ref/             # IRFinder reference
├── xbp1.bed                  # XBP1 regions
└── mt.bed                    # Mitochondrial regions
```

## 🚨 Common Parameter Issues

### Missing References
```bash
# Error: Missing required parameter: --rrna_reference
# Solution: Provide rRNA reference path
--rrna_reference /path/to/rrna.fa
```

### Invalid File Paths
```bash
# Error: File not found
# Solution: Check file existence and permissions
ls -la /path/to/reference/file
```

### Docker Image Issues
```bash
# Error: Docker image not found
# Solution: Verify image availability
docker images | grep tool_name
```

## 📚 Related Documentation

- [README.md](README.md) - Complete pipeline documentation
- [QUICKSTART.md](QUICKSTART.md) - Quick start guide
- [nextflow.config](nextflow.config) - Configuration file reference
- [main.nf](main.nf) - Main workflow file

## 🆘 Need Help?

For parameter-related issues:
1. Check this parameter reference
2. Verify file paths and formats
3. Review error messages in logs
4. Contact: gcarvalhoneto@mednet.ucla.edu
