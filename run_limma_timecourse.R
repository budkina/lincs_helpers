#!/usr/bin/env Rscript
library(optparse)
library("limma")
library(foreach)
library(doParallel)
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
        help="Input directory with counts", metavar="character"),
    make_option(c("-d", "--data"), type="character", default=NULL,
        help="Data directory with cmap data", metavar="character"),
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

if (is.null(opt$data)){
    print_help(opt_parser)
    stop("Input directory should be supplied", call.=FALSE)
}

if (is.null(opt$output)){
    print_help(opt_parser)
    stop("Output directory should be supplied", call.=FALSE)
}

registerDoParallel(opt$cores)
comp_ids <- read.csv(paste0(opt$input,"/","comp_ids.csv"), header = T, sep = '\t')
ids <- as.vector(comp_ids$comp_ids)
gene_info <- read.csv(paste0(opt$data,"/", "GSE92742_Broad_LINCS_gene_info.txt"), sep = '\t')
rownames(gene_info) <- gene_info$pr_gene_id 
foreach(i = 1:length(ids), .combine = "c") %dopar% 
{
	f <- ids[i]
	pert_id <- as.vector(comp_ids$pert_id)[i]
	pert_iname <- as.vector(comp_ids$pert_iname)[i]
	cell_id <- as.vector(comp_ids$cell_id)[i]
	pert_dose <- as.vector(comp_ids$pert_dose)[i]
	design_table <- read.table(paste0(opt$input,"/",f), header = F, sep = '\t')
	counts_table <- read.table(paste0(opt$input,"/",f,"_counts.csv"), header = T, sep = '\t', row.names = 1)
	sample_num <- unique(design_table$V5)
	for (j in sample_num)
	{
		design_table_sample <- design_table[design_table$V5 == j,]
		if (dim(design_table_sample)[1] == 0) {
			next
		}
		
		counts_table_sample <- counts_table[design_table_sample$V1]
		counts_table_norm <- normalizeBetweenArrays(counts_table_sample, method="quantile")
		treatment <- design_table_sample$V2
		time <- design_table_sample$V3
		batch<- design_table_sample$V4
		design <- model.matrix(~0+treatment+treatment:time)
		columns_num <- dim(design)[2]
		colnames(design)[columns_num] <- "treatmentcp_time"
		colnames(design)[columns_num-1] <- "treatmentcontrol_time"
		contrast <- makeContrasts(treatmentcontrol_time - treatmentcp_time, levels=colnames(design))
		fit <- lmFit(counts_table_norm, design)
		fit2 <- contrasts.fit(fit, contrast)
		fit <- eBayes(fit2)
		top_table <- topTable(fit, n = Inf)
		result <-top_table[which(top_table$adj.P.Val < 0.05),]
		result$gene_name <- gene_info[rownames(result),]$pr_gene_symbol
		signature_filename <- paste0(f, "_limma_de_sample_controls_",j,".csv")
		write.table(result, paste0(opt$output,"/", signature_filename), sep = '\t',  quote = F,  row.names = T)
	}
}

signature_filenames_column  <- c()
comp_ids_column <- c()
pert_inames_column <- c()
for(i in 1:length(ids))
{
	f <- ids[i]
	pert_iname <- as.vector(comp_ids$pert_iname)[i]
	design_table <- read.table(paste0("counts1","/",f), header = F, sep = '\t')
	sample_num <- unique(design_table$V5)
	for (j in sample_num)
	{
		signature_filename <- paste0(f, "_limma_de_sample_controls_",j,".csv")
		signature_filenames_column  <- c(signature_filenames_column,signature_filename)
		comp_ids_column <- c(comp_ids_column, f)
		pert_inames_column <- c(pert_inames_column, pert_iname)
	}
}
summary <- data.frame(signature_filenames_column, comp_ids_column, pert_inames_column)
colnames(summary) <- c("signature_filenames", "comp_ids", "pert_inames")
write.table(summary, paste0(opt$output,"/summary.csv"), sep = '\t',  quote = F,  row.names = T)