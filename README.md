# RNA-seq Workflow (nl-rna-seq_wf)

A comprehensive Nextflow-based RNA-seq analysis pipeline for processing and analyzing RNA sequencing data, featuring quality control, alignment, quantification, and advanced analysis capabilities including IRFinder, RNA-SeQC, and Kallisto quantification.

## 🚀 Features

- **Quality Control**: Fastp preprocessing with rRNA and globin RNA contamination detection
- **Alignment**: STAR alignment with multiple read length optimizations (69-151bp)
- **Quantification**: 
  - STAR gene-level counts
  - Kallisto transcript-level quantification
  - Subread featureCounts for gene expression
- **Quality Metrics**: RNA-SeQC comprehensive quality assessment
- **Advanced Analysis**: 
  - IRFinder for intron retention detection
  - Splice junction analysis with BAM2SJ
  - Coverage analysis with Mosdepth
- **Output Formats**: CRAM compression for efficient storage
- **Containerized**: Docker-based execution for reproducibility

## 🏗️ Architecture

```mermaid
graph TB
    subgraph Input
        A1[FASTQ Files] --> B1[Fastp QC]
        A2[Reference Data] --> B1
    end

    subgraph Contamination Check
        B1 --> C1[BWA rRNA Check]
        B1 --> C2[BWA Globin Check]
        C1 --> D1[Samtools Flagstat]
        C2 --> D2[Samtools Flagstat]
    end

    subgraph Alignment & Quantification
        B1 --> E1[STAR Alignment]
        B1 --> E2[Kallisto Quant]
        E1 --> F1[Samtools Index]
        F1 --> F2[Sambamba MarkDup]
        F2 --> G1[FeatureCounts]
        F2 --> G2[RNA-SeQC]
        F2 --> G3[IRFinder]
        F2 --> G4[Mosdepth Coverage]
        F2 --> G5[BAM2SJ]
    end

    subgraph Output
        G1 --> H1[Gene Counts]
        G2 --> H2[QC Reports]
        G3 --> H3[IRFinder Results]
        G4 --> H4[Coverage Files]
        G5 --> H5[Splice Junctions]
        F2 --> H6[CRAM Files]
    end

    style Input fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style Contamination Check fill:#f3e5f5,stroke:#4a148c,stroke-width:2px
    style Alignment fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px
    style Output fill:#fff3e0,stroke:#e65100,stroke-width:2px
```

## 📋 Prerequisites

### System Requirements
- **Operating System**: Linux-based (Ubuntu 18.04+ recommended)
- **CPU**: Minimum 48 cores for optimal performance
- **Memory**: Minimum 192GB RAM for STAR alignment
- **Storage**: Sufficient space for RNA-seq data processing
- **Docker**: Required for containerized execution

### Software Requirements

#### Docker
```bash
# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh

# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker

# Verify installation
docker --version
```

#### Nextflow
```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Verify installation
nextflow -version
```

## 🛠️ Installation

1. **Clone the repository**:
```bash
git clone https://github.com/uclanelsonlab/nl-rna-seq_wf.git
cd nl-rna-seq_wf/
```

2. **Set up Docker images**:
```bash
# Build or pull required Docker images
# Update nextflow.config with your Docker image paths
```

3. **Verify installation**:
```bash
nextflow run main.nf --help
```

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--fastq_r1` | Forward FASTQ file path | `/path/to/sample_R1.fastq.gz` |
| `--fastq_r2` | Reverse FASTQ file path | `/path/to/sample_R2.fastq.gz` |
| `--sample_name` | Sample identifier | `SAMPLE001` |
| `--rrna_reference` | rRNA reference genome | `/path/to/rrna_ref.fa` |
| `--globinrna_reference` | Globin RNA reference | `/path/to/globin_ref.fa` |
| `--gencode_gtf_path` | Gencode GTF annotation | `/path/to/gencode.v43.gtf` |
| `--human_fasta` | Human reference genome | `/path/to/GRCh38.fa` |
| `--star_index_151` | STAR index for 151bp reads | `/path/to/star_index_151` |
| `--kallisto_index` | Kallisto transcript index | `/path/to/kallisto.idx` |


## 📊 Output Structure

The pipeline generates organized outputs in the following structure:

```
results/
├── QC/                          # Quality control reports
│   ├── *.rrna.flagstat.txt     # rRNA contamination stats
│   ├── *.globinrna.flagstat.txt # Globin contamination stats
│   ├── *.html                  # Fastp QC reports
│   └── *.json                  # Fastp QC metrics
├── ALIGNMENT/                   # Alignment files
│   ├── *.ReadsPerGene.out.tab.gz # STAR gene counts
│   ├── *.Log.final.out         # STAR alignment stats
│   ├── *.SJ.out.tab.gz        # Splice junctions
│   ├── *.bam2SJ.out.tab.gz    # Reconstructed junctions
│   └── *.cram                  # Compressed alignment files
├── COUNTS/                      # Gene expression counts
│   ├── *.gene_id.exon.ct      # FeatureCounts output
│   └── *.gene_id.exon.ct.summary # Counts summary
├── QUANT/                       # Kallisto quantification
│   ├── *abundance.h5           # Transcript abundances h5
│   ├── *abundance.tsv           # Transcript abundances
│   └── *run_info.json           # Quantification info
├── IR/                          # IRFinder results
│   ├── *.txt.gz                # Intron retention data
│   └── *.bw                    # Coverage bigWig files
├── MOSDEPTH/                    # Coverage analysis
│   ├── *.mosdepth.global.dist.txt # Global coverage
│   └── *.regions.bed.gz        # Region-specific coverage
└── BAM2SJ/                      # Splice junction analysis
    ├── *_rare_junctions_all.tsv # All detected junctions
    └── *_rare_junctions_filtered.xlsx # Filtered junctions
```

## 🔧 Configuration

### Docker Images
Update `nextflow.config` with your Docker image paths:
```groovy
params {
    bwa_docker = "your-registry/bwa:latest"
    star_docker = "your-registry/star:latest"
    fastp_docker = "your-registry/fastp:latest"
    // ... other images
}
```

### Resource Allocation
Modify process resources in `nextflow.config`:
```groovy
process {
    withLabel: 'star_alignreads' {
        memory = 192.GB
        cpus = 48
    }
}
```

## 🧬 Analysis Components

### 1. Quality Control (Fastp)
- Adapter trimming
- Quality filtering
- Length filtering
- Duplicate detection

### 2. Contamination Detection
- rRNA contamination using BWA
- Globin RNA contamination analysis
- Statistical reporting

### 3. Alignment (STAR)
- Multi-read length optimization
- Splice-aware alignment
- Junction detection
- Gene-level quantification

### 4. Quantification
- **STAR**: Gene-level read counts
- **Kallisto**: Transcript-level abundance
- **FeatureCounts**: Exon-level counting

### 5. Quality Assessment (RNA-SeQC)
- Coverage metrics
- GC bias analysis
- Mapping quality
- Strand specificity

### 6. Advanced Analysis
- **IRFinder**: Intron retention detection
- **BAM2SJ**: Splice junction reconstruction
- **Mosdepth**: Coverage analysis

## 🐛 Troubleshooting

### Common Issues

1. **Memory Errors**:
   - Increase memory allocation in `nextflow.config`
   - Use `-resume` to restart from failed step

2. **Docker Issues**:
   - Verify Docker daemon is running
   - Check image availability and permissions

3. **Reference File Errors**:
   - Ensure all reference files are accessible
   - Verify file formats and integrity

## 📈 Performance

### Resource Recommendations
- **STAR Alignment**: 48 cores, 192GB RAM
- **IRFinder**: 16 cores, 32GB RAM
- **Other processes**: 8-16 cores, 16-32GB RAM

### Expected Runtime
- **Small dataset** (<50M reads): 2-4 hours
- **Medium dataset** (50-100M reads): 4-8 hours
- **Large dataset** (>100M reads): 8-16 hours

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👥 Contact

- **George Carvalho** - gcarvalhoneto@mednet.ucla.edu
- **UCLA Nelson Lab** - https://nelsonlab.ucla.edu/

## �� Acknowledgments

- Nextflow community for the workflow engine
- STAR developers for RNA-seq alignment
- IRFinder team for intron retention analysis
- All contributors to the bioinformatics tools used

---

## 📚 Additional Resources

### FRASER Analysis
For aberrant splicing analysis, see the FRASER configuration in `drop_files/config.yaml`.

### OUTRIDER Analysis
For aberrant expression analysis, use the R script in `script/run_outrider.R`.

### Batch Processing
Use the provided sample annotation files for processing multiple samples.

