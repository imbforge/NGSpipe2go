########################################
##
##  Merge STARR mRNA counts with endogenous mRNA counts for STARR 10X snRNA-seq
## --
## Who:  Martin Oti
## When: 22-March-2022
## --
## Reads in CellRanger count matrix with endogenous mRNA counts, and
## umicounts count matrix with STARR mRNA counts. Appends the sample
## number to the STARR cell barcodes, then merges them with the 10X
## counts and writes the output to a bioconductor SingleCellExperiment
## object.
## --
## Input:
##   <aggr=aggr/>
##   <umicounts=umicounts/>
##   <starrgtf=starr.gtf>
##   <sce=/path/to/output/sce.RDS>
## --
##
########################################
options(stringsAsFactors=F)
library(GenomicRanges)
library(SingleCellExperiment)
library(rtracklayer)

##
## Parse input parms
##
parseArgs <- function(args, string, default=NULL, convert="as.character") {

	if(length(i <- grep(string, args, fixed=T)) == 1) 
		return(do.call(convert, list(gsub(string, "", args[i]))))
    
	if(!is.null(default)) default else do.call(convert, list(NA))
}

args <- commandArgs(trailingOnly=T)
AGGR       <- parseArgs(args, "aggr=", "")          # CellRanger aggr output matrix directory
UMICOUNTS  <- parseArgs(args, "umicounts=", "")     # directory with umicount output files (.umicount.tsv.gz)
STARRGTF   <- parseArgs(args, "gtf=", "")           # gtf file with STARR screen sequences
OUTPUTSCE  <- parseArgs(args, "sce=" , "")    # output SingleCellExperiment RDS file

print(args)
if(length(args) == 0 | args[1] == "-h" | args[1] == "--help")
	stop("Rscript mergeSTARRCounts.R <aggr=aggr/> <umicounts=umicounts/> <starrgtf=genes.gtf> <sce=sce.RDS>")
if(!dir.exists(AGGR)) stop(paste("Directory", AGGR, "does NOT exist"))
if(!dir.exists(UMICOUNTS)) stop(paste("Directory", UMICOUNTS, "does NOT exist"))
if(!file.exists(STARRGTF)) stop(paste("File", STARRGTF, "does NOT exist"))



##
## Read in STARR screening constructs GTF and get gene names
##

starr_gtf <- import.gff(STARRGTF, format="gtf", feature.type="exon")
starr_names <- unique(as.data.frame(starr_gtf)[, c("gene_id", "gene_name")])


##
## Read in Cell Ranger counts matrix with endogenous gene counts
##

matrix_dir <- file.path(AGGR, "outs/count/filtered_feature_bc_matrix/")    # 10x sc/snRNA-seq
if (!dir.exists(matrix_dir)) {matrix_dir <- file.path(AGGR, "outs/filtered_feature_bc_matrix/")}    # multiome RNA component
barcode.path   <- file.path(matrix_dir, "barcodes.tsv.gz")
features.path  <- file.path(matrix_dir, "features.tsv.gz")
matrix.path    <- file.path(matrix_dir, "matrix.mtx.gz")
barcodes       <- read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
features       <- read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
cr_counts      <- readMM(file = matrix.path)
colnames(cr_counts) <- barcodes$V1
rownames(cr_counts) <- features$V1
cr_counts <- cr_counts[features$V3 == "Gene Expression",]


##
## Read in Umicounts counts matrix with STARR screening library mRNA counts
##

# Read in per-sample umicounts matrices
uc_files <- list.files(UMICOUNTS, pattern = ".umicount.tsv.gz", full.names = TRUE)
uc_mats <- lapply(uc_files, read.delim, header = TRUE)

# Add sample number as postfix to cell barcodes (to match with CellRanger's cell barcodes)
uc_mats <- lapply(1:length(uc_mats), function(i){
  colnames(uc_mats[i]) <- paste(colnames(uc_mats[i]), i, sep="-")
  uc_mats[i] })

# Ensure that all matrices have all mRNAs from the GTF (they normally only contain detected mRNAs)
uc_mats <- lapply(uc_mats, function(x){
  fullx <- matrix(0, ncol=ncol(x), nrow=length(starr_names))
  colnames(fullx) <- colnames(x)
  rownames(fullx) <- starr_names$gene_id
  fullx[rownames(x),] <- x
  })

# Merge per-sample matrices
uc_mat <- do.call("cbind", uc_mats)


##
## Merge Umicounts counts with CellRanger counts
##

# First create an empty umicounts matrix with same columns (cell barcodes) as 10X count matrix,
# and concatenate it underneath the 10x count matrix. Then fill in the concatenated empty umicounts
# matrix for those cell barcodes that are shared (detected cells)
detected_cells <- intersect(colnames(cr_counts), colnames(uc_mat))
empty_ucmat <- matrix(0, ncol=ncol(cr_counts), nrow=length(starr_names))
rownames(empty_ucmat) <- starr_names$gene_id
countsmat <- rbind(cr_counts, empty_ucmat)
countsmat[starr_names$gene_id, detected_cells] <- uc_mat[starr_names$gene_id, detected_cells]


##
## Create and save SingleCellExperiment object with counts
##

sce <- SingleCellExperiment(assays=list(counts=as.matrix(countsmat)))

# Add gene IDs & symbols from respective annotations to SingleCellExperiment object
rowData(sce)$gene_id <- rownames(rowData(sce))
rowData(sce)$gene_name  <- c(features$V2, starr_names$gene_name)

saveRDS(sce, file = OUTPUTSCE)

