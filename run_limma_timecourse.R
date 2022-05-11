#!/usr/bin/env Rscript
library(optparse)
library("limma")
library(foreach)
library(doParallel)


option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
        help="Input directory with counts", metavar="character"),
    make_option(c("-c", "--cores"), type="integer", default=16,
        help="Number of cores"),
    make_option(c("-o", "--output"), type="character", default=NULL,
        help="Output directory", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input)){
    print_help(opt_parser)
    stop("Input directory should be supplied", call.=FALSE)
}

if (is.null(opt$output)){
    print_help(opt_parser)
    stop("Output directory should be supplied", call.=FALSE)
}

registerDoParallel(opt$cores)
comp_ids <- read.csv(paste0(opt$input,"/","comp_ids.csv"), header = T)
ids <- as.vector(comp_ids$comp_ids)
foreach(i = 1:length(ids), .combine = "c") %dopar% 
{
	f <- ids[i]
	design_table <- read.table(paste0(opt$input,"/",f), header = F, sep = '\t')
	counts_table <- read.table(paste0(opt$input,"/",f,"_counts.csv"), header = T, sep = '\t', row.names = 1)
	counts_table_norm <- normalizeBetweenArrays(counts_table, method="quantile")
	treatment <- design_table$V2
	time <- design_table$V3
	design <- model.matrix(~0+treatment+treatment:time)
	colnames(design) <- c("treatmentcontrol", "treatmentcp", "treatmentcontrol_time", "treatmentcp_time")
	contrast <- makeContrasts(treatmentcontrol_time - treatmentcp_time, levels=colnames(design))
	fit <- lmFit(counts_table_norm, design)
	fit2 <- contrasts.fit(fit, contrast)
	fit <- eBayes(fit2)
	top_table <- topTable(fit, n = Inf)
	result <-top_table[which(top_table$adj.P.Val < 0.05),]
	write.table(result, paste0(opt$output,"/", f, "_limma_de_all_controls.csv"), sep = '\t',  quote = F,  row.names = T)
}
