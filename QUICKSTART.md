# Quick Start Guide

Get your RNA-seq analysis pipeline running in minutes! This guide provides the essential steps to start using the pipeline.

## 🚀 Prerequisites Check

Before starting, ensure you have:

- [ ] Linux-based operating system
- [ ] Docker installed and running
- [ ] Nextflow 22.04+ installed
- [ ] At least 48 CPU cores and 192GB RAM available
- [ ] Input FASTQ files ready
- [ ] Reference files prepared

## ⚡ Quick Installation

```bash
# 1. Clone the repository
git clone https://github.com/uclanelsonlab/nl-rna-seq_wf.git
cd nl-rna-seq_wf/

# 2. Verify Nextflow installation
nextflow -version

# 3. Verify Docker
docker --version
```

## 🔧 Quick Configuration

1. **Update Docker images** in `nextflow.config`:
```groovy
params {
    bwa_docker = "your-registry/bwa:latest"
    star_docker = "your-registry/star:latest"
    fastp_docker = "your-registry/fastp:latest"
    // ... update other images
}
```

2. **Prepare reference files**:
   - rRNA reference genome
   - Globin RNA reference
   - Human reference genome (FASTA + index)
   - Gencode GTF annotation
   - STAR indices for different read lengths
   - Kallisto transcript index

## 🏃‍♂️ Quick Run

### Minimal Command
```bash
nextflow run main.nf \
    --fastq_r1 /path/to/sample_R1.fastq.gz \
    --fastq_r2 /path/to/sample_R2.fastq.gz \
    --sample_name SAMPLE001 \
    --rrna_reference /path/to/rrna.fa \
    --globinrna_reference /path/to/globin.fa \
    --gencode_gtf_path /path/to/gencode.v43.gtf \
    --human_fasta /path/to/GRCh38.fa \
    --star_index_151 /path/to/star_151 \
    --kallisto_index /path/to/kallisto.idx
```

### Test Run
```bash
# Run with test data (if available)
nextflow run main.nf -profile test

# Or run with minimal parameters to check setup
nextflow run main.nf --help
```

## 📊 Expected Outputs

After successful completion, you'll find:

```
results/
├── QC/                    # Quality control reports
├── ALIGNMENT/            # STAR alignment files
├── COUNTS/               # Gene expression counts
├── QUANT/                # Kallisto quantification
├── IR/                   # IRFinder results
├── MOSDEPTH/             # Coverage analysis
└── BAM2SJ/              # Splice junction analysis
```

## 🐛 Quick Troubleshooting

### Common Issues & Solutions

1. **"Docker daemon not running"**
   ```bash
   sudo systemctl start docker
   sudo usermod -aG docker $USER
   newgrp docker
   ```

2. **"Out of memory"**
   - Increase memory in `nextflow.config`
   - Use `-resume` to restart from last successful step

3. **"Reference file not found"**
   - Check file paths and permissions
   - Verify file formats and integrity

4. **"Process failed"**
   ```bash
   # Check logs
   nextflow log
   
   # Resume from last successful step
   nextflow run main.nf -resume [session_id]
   ```

### Debug Mode
```bash
# Run with detailed logging
nextflow run main.nf -debug

# Check pipeline status
nextflow log
```

## 📈 Performance Tips

- **STAR alignment**: Use 48+ cores and 192GB+ RAM
- **Parallel processing**: Run multiple samples simultaneously
- **Storage**: Ensure sufficient disk space for intermediate files
- **Monitoring**: Use `htop` or `top` to monitor resource usage

## 🔄 Resume & Restart

```bash
# List available sessions
nextflow log

# Resume specific session
nextflow run main.nf -resume [session_id]

# Clean up and restart
rm -rf work/ results/
nextflow run main.nf [parameters]
```

## 📞 Need Help?

- **Documentation**: See full [README.md](README.md)
- **Issues**: Check existing issues or create new ones
- **Contact**: gcarvalhoneto@mednet.ucla.edu

## ✅ Success Checklist

- [ ] Pipeline runs without errors
- [ ] All output directories created
- [ ] Quality control reports generated
- [ ] Alignment files produced
- [ ] Quantification results available
- [ ] Analysis outputs complete

---

**🎯 You're all set!** The pipeline should now be processing your RNA-seq data. Monitor the progress and check outputs as they become available.
