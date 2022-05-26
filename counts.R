library(optparse)
library(Rsubread)

option_list = list(
    make_option(c("-w", "--working_dir"), type="character", default=NULL,
        help="Working directory", metavar="character"),
    make_option(c("-d", "--design"), type="character", default=NULL,
        help="Design dataframe: bam_filename, group_id columns", metavar="character"),
    make_option(c("-a", "--annotation"), type="character", default=NULL,
        help="Path to annotation", metavar="character"),
    make_option(c("-p", "--paired"), default=FALSE,
        help="Paired-end"),
    make_option(c("-c", "--cores"), type="integer", default=16,
        help="Number of cores")
    )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$design)){
    print_help(opt_parser)
    stop("Design matrix should be supplied", call.=FALSE)
}

if (is.null(opt$working_dir)){
    print_help(opt_parser)
    stop("Working directory should be supplied", call.=FALSE)
}

if (is.null(opt$annotation)){
    print_help(opt_parser)
    stop("Annotation file should be supplied", call.=FALSE)
}

# Set working directory
setwd(opt$working_dir)

design <- read.table(opt$design, header=TRUE, sep=' ')
bam_filenames <- design$bam_filename
counts_table <- list()
for (bam in bam_filenames)
{
    fc <- featureCounts(bam, annot.ext=opt$annotation, isGTFAnnotationFile=TRUE, nthreads=opt$cores, isPairedEnd = opt$paired,  GTF.attrType = 'gene_name')
    counts_table[[bam]]<-fc$counts
}

counts_table_df <- as.data.frame(counts_table)
write.table(counts_table_df, file="counts_table.csv", sep=" ", col.names=TRUE, row.names = TRUE, quote = FALSE)
