#####################################
##
## What: CRmotifCounts.R
## Who : Frank Rühle
## When: 01.06.2023
##
## Script to load motif counts and Peak-to-peak and peak-to-gene co-activity correlation data (feature linkage) from Cellranger results.
## To only include motif counts of cells passing QC, this module should be applied after QC filtering.
##
## Args:
## -----
## projectdir      # project directory
##
##
######################################

##
## get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
  
  if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  
  if(!is.null(default)) default else do.call(convert,list(NA))
}

args <- commandArgs(T)
projectdir    <- parseArgs(args,"project=") 
resultsdir    <- parseArgs(args,"res=")   
out           <- parseArgs(args,"outdir=") # output folder
cellranger_aggr_id <- parseArgs(args, "cellranger_aggr_id=") 

runstr <- "Rscript CRmotifCounts.R [projectdir=projectdir]"

# load R environment
renv::use(lockfile=file.path(projectdir, "NGSpipe2go/tools/sc_motifs/renv.lock"))
print(.libPaths())

library(tidyverse)
library(AnnotationDbi)
library(Biobase)
library(data.table)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(reshape2)
library(scater)
library(scran)
library(scuttle)
library(Seurat)
library(Signac)
library(uwot)

# set options
options(stringsAsFactors=FALSE)

# check parameter
print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))
print(paste("cellranger_aggr_id:", cellranger_aggr_id))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "ATAC"


### Read Peak motifs and Peak-to-peak and peak-to-gene co-activity correlation data from Cellranger (feature linkage).
# Note that the motif detection and feature linkage can also be performed using Signac. 
# This would be required if custom peak calling is done, rather than using Cellranger's peaks. 
# More information on this can be found in the Signac multiomic vignette (https://satijalab.org/signac/articles/pbmc_multiomic.html).

cellranger_dir <- file.path(resultsdir, cellranger_aggr_id, "outs")

peaks.path <- file.path(cellranger_dir, "atac_peaks.bed")
peak_annotation.path <- file.path(cellranger_dir, "atac_peak_annotation.tsv")
tf_barcodes.path <- file.path(cellranger_dir, "analysis/tf_analysis/filtered_tf_bc_matrix/barcodes.tsv.gz")
tf_motifs.path <- file.path(cellranger_dir, "analysis/tf_analysis/filtered_tf_bc_matrix/motifs.tsv")
tf_matrix.path <- file.path(cellranger_dir, "analysis/tf_analysis/filtered_tf_bc_matrix/matrix.mtx.gz")

## Cell-to-motif score matrix
motif_counts <- readMM(file = tf_matrix.path)
colnames(motif_counts) <- read.delim(tf_barcodes.path, header = FALSE, stringsAsFactors = FALSE)$V1
motifs <- read.delim(tf_motifs.path, header = FALSE, stringsAsFactors = FALSE)
rownames(motif_counts) <- motifs$V1

# Filter for cells in both (previously processed & filtered) RNA & ATAC
shared_barcodes <- intersect(colnames(motif_counts), colnames(sobj))
motif_counts <- motif_counts[,shared_barcodes]
write.table(as.data.frame(motif_counts), file = file.path(out, "Cellranger_motif_counts.txt"), quote=F, sep = "\t", row.names = T)


## Peak-to-motif mapping (10X Multiome only, even though it only uses ATAC data)
peaks2motifs.path <- file.path(cellranger_dir, "analysis/tf_analysis/peak_motif_mapping.bed")
peaks2motifs <- rtracklayer::import(peaks2motifs.path, format = "BED")
write.table(as.data.frame(peaks2motifs), file = file.path(out, "Cellranger_peaks2motifs.txt"), quote=F, sep = "\t", row.names = F)


## Feature linkage data
feature_linkage.path <- file.path(cellranger_dir, "analysis/feature_linkage/feature_linkage.bedpe")
feature_linkage <- rtracklayer::import(feature_linkage.path, format = "BEDPE")
colnames(mcols(feature_linkage))[3:5] <- c("significance", "distance", "linkage_type")
write.table(as.data.frame(feature_linkage), file = file.path(out, "Cellranger_feature_linkage.txt"), quote=F, sep = "\t", row.names = F)




#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/CRmotifCounts_session_info.txt"))
save(motif_counts, peaks2motifs, feature_linkage, file=paste0(out,"/CRmotifCounts.RData"))

