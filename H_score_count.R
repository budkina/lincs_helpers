#!/usr/bin/env Rscript
library(optparse)
library(foreach)
library(doParallel)
library(Hobotnica)
library(edgeR)

get_H<-function(signature_matrix, trt, sig_length)
{
    if (sig_length < 2)
    {
        return(0.5)
    }

    dist_matrix<-dist(t(signature_matrix))
    return(Hobotnica(dist_matrix, trt))
}

get_random_H <-function(H, sig, trt, data, sig_length, cores)
{
    if (sig_length < 2)
    {
        return(list(H=c(), pvalue=NA))
    }

    # P-value calculation
    registerDoParallel(cores)
    H_random<-foreach(i = 1:100000, .combine = "c") %dopar%
        {
            random_signature_matrix<-data[sample(nrow(data), sig_length), ]
            H_chunk<-get_H(random_signature_matrix, trt, sig_length)
            H_chunk
        }

    # pvalue calculation with pseudo count
    pvalue<-(sum(H_random > H)+1)/(length(H_random)+1)
    return(list(H=H_random, pvalue=pvalue))
}

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
        help="Input directory with signatures", metavar="character"),
    make_option(c("-t", "--compound"), type="character", default=NULL,
        help="Compound name for test", metavar="character"),
    make_option(c("-m", "--matrix"), type="character", default=NULL,
        help="Counts table filename", metavar="character"),
    make_option(c("-c", "--cores"), type="integer", default=16,
        help="Number of cores"),
    make_option(c("-d", "--design"), type="character", default=NULL,
        help="Design matrix filename", metavar="character"))
    make_option(c("-o", "--output"), type="character", default=NULL,
        help="Output prefix", metavar="character"))

summary <- read.table(paste0(opt$input,"/summary.csv"), header = T, sep = '\t')
compound_sigs <- summary[summary$pert_iname == opt$compound,]

# read counts table
data_matrix <- read.csv(opt$matrix, sep = ' ')

# read design matrix
design_matrix <- read.csv(opt$design, sep = ' ', header = F)
data_matrix = data_matrix[design_matrix$V1]

#### Get H-scores for compound signatures ####

signatures <- c()
H_scores <- c()
pvalues <- c()
for (file in compound_sigs$signature_filenames)
{
	sig_df <- read.csv(paste0(opt$input, '/', file), sep = '\t')
	sig <- sig_df$gene_name
	sig_length <- length(sig)
	if (sig_length < 2)
	{
		next
	}
	cm <- data_matrix[sig,]
	cm[is.na(cm)] <- 0
	if (any(data.frame(as.matrix(colSums(cm)))[,1] == 0))
	{
		next
	}
	submatrix_cpm <- cpm(cm)
	trt <- design_matrix$V2
	dist_matrix<-dist(t(submatrix_cpm))
	H<-Hobotnica(dist_matrix, trt)
	H_scores <- c(H_scores, H)
	result <- get_random_H(H, sig, trt, submatrix_cpm, sig_length, opt$cores)
	pvalues <- c(pvalues, result$pvalue)
	signatures <- c(signatures, file)
}

df_res <- data.frame(signatures, H_scores, pvalues)
write.table(df_res, paste0(opt$output,'_sig_H_scores.csv'), sep = '\t',  quote = F,  row.names = F)

#### Get H-scores for all signatures ####

signatures <- c()
H_scores <- c()
for (file in summary$signature_filenames)
{
	sig_df <- read.csv(paste0(opt$input, '/', file), sep = '\t')
	sig <- sig_df$gene_name
	sig_length <- length(sig)
	if (sig_length < 2)
	{
		next
	}
	cm <- data_matrix[sig,]
	cm[is.na(cm)] <- 0
	if (any(data.frame(as.matrix(colSums(cm)))[,1] == 0))
	{
		next
	}
	submatrix_cpm <- cpm(cm)
	trt <- design_matrix$V2
	dist_matrix<-dist(t(submatrix_cpm))
	H<-Hobotnica(dist_matrix, trt)
	H_scores <- c(H_scores, H)
	signatures <- c(signatures, file)
}

df_res <- data.frame(signatures, H_scores)
write.table(df_res, paste0(opt$output,'_all_H_scores.csv'), sep = '\t',  quote = F,  row.names = F)