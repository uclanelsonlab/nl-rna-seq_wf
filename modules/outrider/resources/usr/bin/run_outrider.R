#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(magrittr)
    library(tools)
    library(optparse)
})
 
option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset with mutliple sample feature counts", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="outrider_results.tsv", 
              help="output file name [default= %default]", metavar="character"),
    make_option(c("-t", "--tissue"), type="character", default="fibroblast", 
              help="type of tissue [default= %default]", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
tissue = opt$tissue

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
message(date(), ": Creating OutriderDataSet object")
cts <- read.table(file = opt$file); # "featureCounts_fibroblast_24-07-18.tsv"
ods <- OutriderDataSet(countData = cts);
# Preprocessing
message(date(), ": Running filterExpression")
ods <- filterExpression(object = ods, minCounts = TRUE, filterGenes = TRUE);

message(date(), ": Running estimateSizeFactors")
ods <- estimateSizeFactors(ods)
## find optimal encoding dimension
# a <- 5 
# b <- min(ncol(ods), nrow(ods)) / 3   # N/3
# Nsteps <- min(20, b)   # Do at most 20 steps or N/3
# Do unique in case 2 were repeated
implementation = 'autoencoder'
# pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique
# message(date(), ": Running findEncodingDim")
# ods <- findEncodingDim(ods, params = pars_q, implementation = implementation)
# opt_q <- getBestQ(ods) # 22
if (tissue == 'blood') {
  opt_q <- 10;
} else if (tissue == 'fibroblast') {
  opt_q <- 22;
} else if (tissue == 'muscle') {
  opt_q <- 8;
} else if (tissue == 'fat') {
  opt_q <- 2;
} else if (tissue == 'liver') {
  opt_q <- 2;
}
## fit OUTRIDER
message(date(), ": Running SizeFactor estimation")
ods <- estimateSizeFactors(ods)
message(date(), ": Controlling for confounders")
implementation <- tolower(implementation)
ods <- controlForConfounders(ods, q=opt_q, implementation=implementation)
if(grepl("^(peer|pca)$", implementation)){
    message(date(), ": Fitting the data")
    ods <- fit(ods)
}
message("outrider fitting finished")

# P value calculation
message(date(), ": P-value calculation ...")
ods <- computePvalues(ods, alternative="two.sided", method="BY")
message(date(), ": Zscore calculation ...")
ods <- computeZscores(ods, peerResiduals=grepl('^peer$', "autoencoder"))
message(date(), ": Creating results table")
res <- results(ods, padjCutoff = 0.05, all = TRUE)
# Add fold change
res[, foldChange := round(2^l2fc, 2)]
# Subset to significant results
padj_cols <- grep("padjust", colnames(res), value=TRUE)
res <- res[do.call(pmin, c(res[,padj_cols, with=FALSE], list(na.rm = TRUE))) 
                <= 0.05]
# Save results
file <- paste0('ods_', opt$out, '.rds');
saveRDS(object = ods, file = file);
fwrite(res, opt$out, sep = "\t", quote = F)

message(date(), sessionInfo())