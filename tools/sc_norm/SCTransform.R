#####################################
##
## What: SCTransform.R
## Who : Frank Rühle
## When: 01.06.2023
##
## Script to normalize gene expression data in Seurat object.
##
## Args:
## -----
## targets=targets.txt      # file describing the targets 
##
##
######################################

renv::use(lockfile='NGSpipe2go/tools/sc_norm/renv.lock')
print(.libPaths())

options(stringsAsFactors=FALSE)
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
CORES <- 2


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

runstr <- "Rscript SCTransform.R [projectdir=projectdir]"

print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "RNA"

sobj <- SCTransform(sobj)
sobj <- RunPCA(sobj)
sobj <- RunUMAP(sobj, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')




#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/SCTransform_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(resultsdir, file=paste0(out,"/SCTransform.RData"))

