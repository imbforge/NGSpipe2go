#####################################
##
## What: CapSTARRseq_FC.R
## Who : Martin Oti, based on DE_DESeq2.R by Sergi Sayols
## When: 01-03-2022
##
## Script to perform CapSTARR-seq fold change-based analysis, instead of differential expression analysis.
## Required when no STARR-seq plasmid/viral DNA was extracted from the samples, and RNA expression
## has to be compared with input DNA.
## Performs a one-sample t-test to check if the log2 fold changes are significantly different from 0.
## Therefore, it requires STARR-seq RNA replicates for every condition tested (DNA can be the same).
## (There is also an outlier test that does not require replicates, but this approach is much less
## reliable, so the script is not currently set up to use only this approach.)
##
## Args:
## -----
## targets=targets.txt      # file describing the targets 
##                          # must fit the format expected in ChIP-seq pipeline
## gtf=gene_model.gtf       # gene model in gtf format - for TPM calculation
## cwd=.                    # current working directory where the files .tsv files are located
## out=CapSTARRseq_FoldChange     # prefix filename for output
## pattern=","\\.readcounts.tsv"  # pattern for the count files
##
######################################
options(stringsAsFactors=FALSE)
library(qvalue)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(pheatmap)
library(viridis)

##
## get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
  
  if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  
  if(!is.null(default)) default else do.call(convert,list(NA))
}

args <- commandArgs(T)
ftargets     <- parseArgs(args,"targets=","targets.txt")     # file describing the targets
gene.model   <- parseArgs(args,"gtf=","")       # gtf gene model
cwd          <- parseArgs(args,"cwd=","./")     # current working directory
out          <- parseArgs(args,"out=","CapSTARRseq_FoldChange") # output filename
pattern      <- parseArgs(args,"pattern=","\\.readcounts.tsv") # output filename

runstr <- "Rscript CapSTARRseq_FC.R [targets=targets.txt] [gtf=] [cwd=.] [out=CapSTARRseq_FC] [pattern=RE]"
if(!file.exists(ftargets))   stop(paste("File",ftargets,"does NOT exist. Run with:\n",runstr))
if(!file.exists(cwd))        stop(paste("Dir",cwd,"does NOT exist. Run with:\n",runstr))
if(!file.exists(gene.model)) stop(paste("GTF File:", gene.model, " does NOT exist. Run with: \n", runstr))

if(!dir.exists(out)){dir.create(out, recursive = TRUE)}

##
## read in and process targets file
##
# load targets
targets <- read.delim(ftargets,head=T,colClasses="character",comment.char="#")
if(!all(c("IP", "IPname", "INPUT", "INPUTname", "group", "Replicate") %in% colnames(targets))) stop("targets file must have at least 6 columns and fit the format expected for the ChIP-seq pipeline (except the unnecessary 'PeakCaller' 7th column)")

# grep sample identifier in count file names
countfiles <- list.files(cwd)
countfiles <- countfiles[grep(pattern, countfiles)] # filter for valid count files

# remove file ending of target file names
targets_samples <- unique(c(targets$IP, targets$INPUT)) 

index_targetsamples <- sapply(paste0(targets_samples, "\\."), grep,  countfiles) # grep targets in countfiles

## check matching ambiguity
if(is(index_targetsamples, "list")) { # list means either zero or multiple matches
  if(any(x <- sapply(index_targetsamples, length)>1)) {
    stop(paste("\nA targets.txt entry matches multiple count file names"),
         "\ncount file names: ", paste(countfiles[unlist(index_targetsamples[x])], collapse=", "))
  }
  warning(paste("Entries in targets.txt which are not found in list of count files are removed from targets table:", 
                paste(targets$file[sapply(index_targetsamples, length)==0], collapse=", ")))
  targets_samples <- targets_samples[!(sapply(index_targetsamples, length)==0)] # remove target entries
  targets <- targets[targets$IP %in% targets_samples | targets$INPUT %in% targets_samples,]
  index_targetsamples <- sapply(targets_samples, grep, countfiles) # recreate index vector after removal of targets.txt entries
}

if (length(unique(index_targetsamples)) < length(countfiles)) { # check for count file names not included
  warning(paste("\nCount file names not included in targets.txt are ignored: ",  paste(countfiles[-index_targetsamples], collapse=", ")))
}

# replace IP & INPUT entries in targets data.frame by corresponding count file names
countfiles_targetsamples <- countfiles[index_targetsamples]
targets$IP <- countfiles_targetsamples[sapply(targets$IP, grep, countfiles_targetsamples)] 
targets$INPUT <- countfiles_targetsamples[sapply(targets$INPUT, grep, countfiles_targetsamples)] 

##
## calculate gene transcript lengths (union of all annotated exons) using rtracklayer & GenomicRanges
##
gtf <- import.gff(gene.model, format="gtf", feature.type="exon")
gtf.flat <- unlist(reduce(split(gtf, elementMetadata(gtf)$gene_id)))
gene.lengths <- tapply(width(gtf.flat), names(gtf.flat), sum)

# make sure we have a gene_name column (needed as output in the report later)
if(! "gene_name" %in% colnames(mcols(gtf))) {
  if("gene_id" %in% colnames(mcols(gtf))) {
    gtf$gene_name <- gtf$gene_id
  } else {
    gtf$gene_name <- NA
  }
}

# create a txdb object to collect the genes coordinates for later usage
txdb  <- makeTxDbFromGRanges(gtf)
genes <- as.data.frame(genes(txdb))


## Read in and preprocess counts data

exprDFs <- list()  # list to hold various expression matrices (counts, TPMs)

# First read in counts data
countslist_ip <- lapply(targets$IP, function(cf){
  df <- read.table(file.path(cwd, cf), header = FALSE, row.names = 1)
  df[, 1, drop = FALSE]
})
countslist_input <- lapply(targets$INPUT, function(cf){
  df <- read.table(file.path(cwd, cf), header = FALSE, row.names = 1)
  df[, 1, drop = FALSE]
})
exprDFs$Counts <- cbind(do.call("cbind", countslist_ip), do.call("cbind", countslist_input))
colnames(exprDFs$Counts) <- c(targets$IPname, targets$INPUTname)
rm(countslist_ip, countslist_input)

# Transcripts per million (TPMs)
# TPM function
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
exprDFs$TPM <- as.data.frame(apply(exprDFs$Counts, 2, tpm, gene.lengths[rownames(exprDFs$Counts)]))


## Log2 Fold Changes

exprDFs$Log2FC <- log2((exprDFs$TPM[,targets$IPname]+1) / (exprDFs$TPM[,targets$INPUTname]+1))


##
## Log2FC-based enhancer/repressor identification based on replicates: one-sample t-test
##

# Run t-tests
ttests <- lapply(unique(targets$group), function(grp){
  tt <- apply(exprDFs$Log2FC[, which(targets$group == grp)], 1, function(r){if(var(r)>0){t.test(r)}})
  means <- rowMeans(exprDFs$Log2FC[, which(targets$group == grp)], na.rm = FALSE)
  pvalues <- unlist(sapply(tt, "[[", "p.value"))
  qvals <- qvalue(pvalues)
  df <- data.frame("meanLog2FC" = means, "pvalue" = NA, "qvalue" = NA, "localFDR" = NA)
  df[names(pvalues), "pvalue"] <- pvalues
  df[names(qvals$qvalues), "qvalue"] <- qvals$qvalues
  df[names(qvals$lfdr), "localFDR"] <- qvals$lfdr
  df
})
names(ttests) <- unique(targets$group)

exprDFs$TTests <- do.call("cbind", ttests)


##
## Log2FC-based enhancer/repressor identification per sample: outlier test
## Does not require or use replicates, but is much less reliable than above t-test
##

# Functions for Median Absolute Deviation (MAD) calculations
# Function to calculate Relative Deviation (num MADs from median)
MAD_scores <- function(x) {
  medx <- stats::median(x)
  madx <- stats::mad(x)
  return((x-medx)/madx)
}
# Function to calculate outliers based on Median Absolute Deviation
MAD_outliers <- function(x, mads = 3) {
  medx <- stats::median(x)
  madx <- stats::mad(x)
  return(x > (medx+(mads*madx)) | x < (medx-(mads*madx)))
}

outliers <- list()
outliers$MADscores <- apply(exprDFs$Log2FC, 2, MAD_scores)
outliers$MADoutliers <- apply(exprDFs$Log2FC, 2, MAD_outliers)
outliers$MADoutliers[outliers$MADoutliers == FALSE] <- 0
outliers$MADoutliers[outliers$MADoutliers == TRUE] <- 1
for (mad in names(outliers)){
  colnames(outliers[[mad]]) <- paste(mad, colnames(outliers[[mad]]), sep = ".")
}

exprDFs$Outliers <- do.call("cbind", outliers)


## Annotate with genomic coordinates and write out data to Excel file

# extract the gene_name and genomic coordinates of each gene
locations.df <- data.frame(gene_name=gtf$gene_name[match(rownames(exprDFs$TPM), gtf$gene_id)], row.names=rownames(exprDFs$TPM))
i <- match(rownames(locations.df), genes$gene_id)
locations.df$chr    <- genes$seqnames[i]
locations.df$start  <- genes$start[i]
locations.df$end    <- genes$end[i]
locations.df$strand <- genes$strand[i]

# merge location and quantification for all expression-related data frames
exprDFs <- lapply(exprDFs, function(df){
  ann_df <- merge(locations.df, df, by=0)
  colnames(ann_df)[1]        <- "gene_id"
  ann_df
})


# write to Excel file
write.xlsx(exprDFs, file=paste0(out, "/CapSTARRseq_FoldChange_Results.xlsx"), row.names=F)


##########################################################################################
##
## TO DO: Plots (heatmaps, PCA plots, MA plots, Volcano plots)
## Use plotting code from DE_DESeq2.R as a guide
##
##########################################################################################


#save the sessionInformation and RData file
writeLines(capture.output(sessionInfo()),paste(out, "/CapSTARRseq_FoldChange_session_info.txt", sep=""))
save(exprDFs, gtf, file=paste0(out,"/CapSTARRseq_FoldChange.RData"))

