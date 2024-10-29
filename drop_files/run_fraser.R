suppressPackageStartupMessages({
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(FRASER)
    library(optparse)
    library(DelayedArray)
    library(dplyr)
    library(readr)
})

option_list = list(
    make_option(c("-f", "--file"), type="character", default="sample_annotation_fraser_fibroblast_less.csv",
              help="Input file with BAM paths", metavar="character"),
    make_option(c("-t", "--tissue"), type="character", default="fibroblast",
              help="type of tissue [default= %default]", metavar="character"),
    make_option(c("-c", "--threads"), type="numeric", default=10,
              help="number of threads [default= %default]", metavar="numeric"),
    make_option(c("-d", "--directory"), type="character", default="FRASER_output",
              help="Output directory [default= %default]", metavar="character"),
    make_option(c("-n", "--name"), type="character", default="FRASER_analysis",
              help="Output name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tissue = opt$tissue
threads = opt$threads
directory = opt$directory
name = opt$name
file = opt$file
output = paste(name, tissue, sep = '_')
# hard coded parameters
assemblyVersion <- "hg38"
minExpressionInOneSample <- 20
quantileForFiltering <- 0.75
quantileMinExpression <- 10
minDeltaPsi <- 0.05
psiTypes <- c("jaccard")
implementation <- "PCA"
mp <- 3 # maxTestedDimensionProportion
padjCutoff <- 0.1
deltaPsiCutoff <- 0.1

if(.Platform$OS.type == "unix") { register(MulticoreParam(workers=min(35, multicoreWorkers())))
} else { register(SnowParam(workers=min(35, multicoreWorkers()))) }

message(date(), ": Loading sample annotation")
sampleTable <- fread(file)
settings <- FraserDataSet(colData=sampleTable, workingDir=directory)
strandSpecific(settings) <- sampleTable$strandSpecific
message(date(), ": Counting the reads")
fds <- countRNAData(settings)

# 1. PSIValues
fds <- calculatePSIValues(fds)
gc()
# 2. Filter
fds <- filterExpressionAndVariability(fds,
        minExpressionInOneSample = minExpressionInOneSample,
        quantile=quantileForFiltering,
        quantileMinExpression=quantileMinExpression,
        minDeltaPsi = minDeltaPsi,
        filterOnJaccard=TRUE,
        filter=FALSE)
gc()

# 3. Load PSI data
fitMetrics(fds) <- psiTypes
# Run hyper parameter optimization
message(date(), ": Fit the splicing model for each metric")
# fds <- FRASER(fds)
# Get range for latent space dimension
a <- 2
b <- min(ncol(fds), nrow(fds)) / mp   # N/mp
maxSteps <- 12
if(mp < 6){
  maxSteps <- 15
}
Nsteps <- min(maxSteps, b)
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique
fds <- optimHyperParams(fds, type=psiTypes,
                            implementation=implementation,
                            q_param=pars_q,
                            plot = FALSE)
gc()

# 4. Fit autoencoder
# run it for every type
message(date(), ": Running fit")
currentType(fds) <- psiTypes
q <- bestQ(fds, psiTypes) # 4
message(date(), ": Best Q: ", q)
verbose(fds) <- 3   # Add verbosity to the FRASER object
fds <- fit(fds, q=q, type=psiTypes, iterations=15, implementation=implementation)
gc()
fds <- saveFraserDataSet(fds)
# 5. Annotation
message(date(), ": Annotate introns with the HGNC symbols of the corresponding gene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
orgDb <- org.Hs.eg.db
gc()
# Annotate the fds with gene names and save it as a new object
fds <- annotateRangesWithTxDb(fds, txdb=txdb, orgDb=orgDb)
fds <- annotateIntronReferenceOverlap(fds, txdb)
gc()
fds <- saveFraserDataSet(fds)
# Calculate stats
# 7. Pvalues
message(date(), ": Calculate p-values")
fds <- calculatePvalues(fds, type=psiTypes)
gc()
fds <- saveFraserDataSet(fds)
# 8. Adjusted p-values
message(date(), ": Calculate p-adj")
fds <- calculatePadjValues(fds, type=psiTypes)
gc()
fds <- saveFraserDataSet(fds)

# 9. Results

message(date(), ": Saving project")
workingDir(fds) <- directory
name(fds) <- paste(name, tissue, sep = '_')
saveFraserDataSet(fds, dir=workingDir(fds), name=name(fds))

message(date(), ": Getting results")
# Extract results per junction
res_junc <- results(fds, psiType=psiTypes,
                    padjCutoff=padjCutoff,
                    deltaPsiCutoff=deltaPsiCutoff)
res_junc_dt   <- as.data.table(res_junc)
# number of samples per gene and variant
res_junc_dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
res_junc_dt[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end,strand"]
res_junc_dt <- annotatePotentialImpact(result=res_junc_dt, txdb=txdb, fds=fds)
res_junc_dt <- flagBlacklistRegions(result=res_junc_dt, assemblyVersion=assemblyVersion)
write_tsv(res_junc_dt, file=paste(output, "resultTableJunc.csv", sep = '_'))
print('Results per junction extracted')

# Extract full results by gene
res_gene <- results(fds, psiType=psiTypes,
                    aggregate=TRUE, collapse=FALSE,
                    all=TRUE)
res_genes_dt   <- as.data.table(res_gene)
print('Results per gene extracted')
write_tsv(res_genes_dt, file=paste(output, "resultTableGene_full.csv", sep = '_'))

# Subset gene results to aberrant
padj_cols <- grep("padjustGene", colnames(res_genes_dt), value=TRUE)
res_genes_dt <- res_genes_dt[do.call(pmin, c(res_genes_dt[,padj_cols, with=FALSE],
                                                list(na.rm = TRUE))) <= padjCutoff &
                                    abs(deltaPsi) >= deltaPsiCutoff &
                                    totalCounts >= 5,]
res_genes_dt <- annotatePotentialImpact(result=res_genes_dt, txdb=txdb, fds=fds)
res_genes_dt <- flagBlacklistRegions(result=res_genes_dt, assemblyVersion=assemblyVersion)
write_tsv(res_genes_dt, file=paste(output, "resultTableGene_aberrant.csv", sep = '_'))
