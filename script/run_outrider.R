library(OUTRIDER);
library(BoutrosLab.utilities);
library(optparse);

# Define command line options
option_list <- list(
  make_option(
    c("--path"),
    type = "character",
    default = NULL,
    help = "Path to the featureCounts directory",
    metavar = "character"
  ),
  make_option(
    c("--tissue"),
    type = "character",
    default = "fibroblast",
    help = "Tissue type [default= %default]",
    metavar = "character"
  )
);

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# Validate path
if (is.null(opt$path)) {
  print_help(opt_parser)
  stop("Please provide a path using --path argument")
}

path <- opt$path

# Validate path exists
if (!dir.exists(path)) {
  stop(sprintf("Directory '%s' does not exist", path))
}

files <- list.files(path = path);
sample.names <- substr(x = files, start = 0, stop = nchar(x = files) - 26);
df <- data.frame(file.name = files, sample.name = sample.names);
cts <- read.table(file = paste(path, df$file.name[1], sep = .Platform$file.sep), header = TRUE);
colnames(cts) <- c('Geneid', df$sample.name[1]);
for(i in 2:nrow(df)) {
  print(i);
  .cts <- read.table(file = paste(path, df$file.name[i], sep = .Platform$file.sep), header = TRUE);
  colnames(.cts) <- c('Geneid', df$sample.name[i]);
  cts <- merge(x = cts, y = .cts, by = 'Geneid');
}
rownames(cts) <- cts$Geneid
cts <- cts[,-1];

ods <- OutriderDataSet(countData = cts);
ods <- filterExpression(object = ods, minCounts = TRUE, filterGenes = TRUE);
ods <- estimateSizeFactors(object = ods);
# ods <- plotCountCorHeatmap(object = ods);

# # don't run this if you know the encoding dimensions
# ods <- findEncodingDim(ods)
# plotEncDimSearch(ods)

# assuming you know the encoding dimensions
# we can omit the BPPARAM parameter since parallelization works
# we can add a number of iterations to end earlier than the default of 15

encoding.dimensions <- 22;

ods <- controlForConfounders(
  ods = ods, 
  q = encoding.dimensions, 
  iterations = 5, 
  BPPARAM = SerialParam());

# ods <- plotCountCorHeatmap(object = ods, normalized = TRUE);
ods <- fit(object = ods);
ods <- computePvalues(object = ods, alternative = 'two.sided', method = 'BY');
ods <- computeZscores(ods = ods);

# file <- paste0('ods_', tissue.type, '.rds');
file <- generate.filename(project.stem = 'outrider', file.core = paste0(opt$tissue, '_hg38'), extension = 'rds');
saveRDS(object = ods, file = file);

res <- OUTRIDER::results(object = ods, padjCutoff = 0.05, zScoreCutoff = 1);
dim(x = res);
file <- generate.filename(project.stem = 'outrider', file.core = paste0(opt$tissue, '_hg38'), extension = 'txt');
write.table(
  x = res, 
  file = file, 
  quote = FALSE, 
  sep = '\t', 
  row.names = FALSE);
