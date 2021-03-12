########################################
## Anke Busch
## a.busch@imb-mainz.de
##
## 2021
########################################

options(stringsAsFactors=FALSE)
library(rtracklayer)
library(optparse) 

# input arguments
option_list <- list(
   make_option(c("-c", "--countfile"), type="character", default=NULL, 
               help="count tables [default = %default]", metavar="character"),
   make_option(c("-o", "--outfile"), type="character", default=NULL, 
               help="output file [default = %default]", metavar="character"),
   make_option(c("-g", "--gtf"), type="character", default=NULL,
               help="gtf annotation file [default = %default]", metavar="character"),
   make_option(c("-f", "--feature"), type="character", default="exon",
               help="annotation feature to count reads in [default = %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# managing null arguments
if (is.null(opt$countfile) | is.null(opt$outfile) | is.null(opt$gtf)){
   print_help(opt_parser)
   stop("The count table (-c), the output file (-o) and the gtf annotation file (-g) must all be supplied.", call.=FALSE)
}

# check if gtf file exists
if(!file.exists(opt$gtf)) stop(paste("GTF file: ",opt$gtf," does NOT exist.\n"))

#example:
#opt$countfile <- "/PROJECTFOLDER/results/subread-count/sample.readcounts.tsv"
#opt$outfile   <- "/PROJECTFOLDER/results/TPMs/sample.tpm.tsv"
#opt$gtf       <- "/ANNOTATIONFOLDER/gencode.v31.annotation.gtf"
#opt$feature   <- "exon" 

# read count file
counts.df <- read.delim(opt$countfile,header=FALSE)
names(counts.df) <- c("gene_id","count")

# get the gene length (as sum of all features specified in opt$feature per gene)
gtf.feature <- import.gff(opt$gtf, format="gtf", feature.type=opt$feature)
gtf.flat.feature <- unlist(reduce(split(gtf.feature, elementMetadata(gtf.feature)$gene_id)))
gene.lengths.feature <- tapply(width(gtf.flat.feature), names(gtf.flat.feature), sum)

## add gene length to the count data frame
counts.df$gene.length <- gene.lengths.feature[counts.df$gene_id]

# define TPM function 
tpm <- function(counts, lengths) {
   rate <- counts / lengths
   rate / sum(rate) * 1e6
}

# get TPMs
tpm.df <- data.frame(gene_id = counts.df$gene_id,
                     TPM     = tpm(counts.df$count, counts.df$gene.length))

# write to file
write.table(tpm.df,file=opt$outfile,row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)



