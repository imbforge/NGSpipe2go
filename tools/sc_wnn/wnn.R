#####################################
##
## What: wnn.R
## Who : Frank Rühle
## When: 02.06.2023
##
## Script to compute a joint neighbor graph that represent both the gene expression and DNA accessibility measurements using the weighted nearest neighbor (WNN) method.
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
knn           <- parseArgs(args,"knn=", convert="as.numeric")
knnRange      <- parseArgs(args,"knnRange=", convert="as.numeric")
clusterAlg    <- parseArgs(args,"clusterAlg=", convert="as.numeric")
clusterRes    <- parseArgs(args,"clusterRes=", convert="as.numeric")
skipFirstLSIcomp <- parseArgs(args,"skipFirstLSIcomp=", default=0, convert="as.numeric")
batchCorrection <- parseArgs(args,"batchCorrection=", default=FALSE, convert="as.logical")

runstr <- "Rscript wnn.R [projectdir=projectdir]"

# load R environment
renv::use(lockfile=file.path(projectdir, "NGSpipe2go/tools/sc_wnn/renv.lock"))
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
addTaskCallback(function(...) {set.seed(100);TRUE})

# check parameter
print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))
print(paste("knn:", knn))
print(paste("knnRange:", knnRange))
print(paste("clusterAlg:", clusterAlg))
print(paste("clusterRes:", clusterRes))
print(paste("skipFirstLSIcomp :", skipFirstLSIcomp ))
print(paste("batchCorrection : ", batchCorrection))

# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))


# We calculate a WNN graph (Weighted Nearest Neighbor), representing a weighted combination of RNA and ATAC-seq modalities. 
# We use this graph for UMAP visualization and clustering.
# These will be stored in the neighbors slot, 
# and can be accessed using sobj[['weighted.nn']]
# The WNN graph can be accessed at sobj[["wknn"]], 
# and the SNN graph used for clustering at sobj[["wsnn"]]
# Cell-specific modality weights can be accessed at sobj$RNA.weight

if(batchCorrection==TRUE) {
	assay2use = "integrated"
	reduction_list_to_use = list("pca.corrected", "lsi.corrected")
	pca_reduction = "pca.corrected"
} else {
	assay2use = "SCT"
	reduction_list_to_use = list("pca.sct", "lsi")
	pca_reduction = "pca.sct"
}

sobj <- FindMultiModalNeighbors(
  object = sobj,
  reduction.list = reduction_list_to_use, 
  dims.list = list(1:50, (1+skipFirstLSIcomp):50), 
  k.nn = knn, 
  knn.range= knnRange, 
  knn.graph.name = "wknn",
  snn.graph.name = "wsnn",
  weighted.nn.name = "weighted.nn",
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
sobj <- RunUMAP(
  object = sobj,
  nn.name = "weighted.nn",
  assay = assay2use,
  reduction = pca_reduction, 
  reduction.name = "umap.wnn",
  verbose = TRUE
)


# build a joint TSNE visualization
sobj <- RunTSNE(
  object = sobj,
  nn.name = "weighted.nn",
  assay = assay2use,
  reduction = pca_reduction, 
  reduction.name = "tsne.wnn",
  verbose = TRUE
)


# By default, cells are colored by their identity class (can be changed with the group.by parameter). 
umapSample <- DimPlot(sobj, label = TRUE, repel = TRUE, reduction = "umap.wnn", group.by = "sample") + ggtitle("WNN UMAP")
tsneSample <- DimPlot(sobj, label = TRUE, repel = TRUE, reduction = "tsne.wnn", group.by = "sample") + ggtitle("WNN TSNE")

ggsave(plot=umapSample, filename=file.path(out, "umap_anno_sample.pdf"))
ggsave(plot=tsneSample, filename=file.path(out, "tsne_anno_sample.pdf"))


### WNN clustering
sobj <- FindClusters(sobj, graph.name = "wsnn", algorithm = clusterAlg, resolution=clusterRes, verbose = FALSE)

# Save clusters in a different metadata column, because it gets overwritten every time FindClusters is called
sobj$clusters_wnn <- sobj$seurat_clusters
# Also annotate nuclei with per-sample clusters
sobj$Sample_WNNClusters <- paste(sobj$sample, sobj$clusters_wnn, sep = "_")

WNNumap <- DimPlot(sobj, reduction = "umap.wnn", group.by = "clusters_wnn", label = TRUE, repel = TRUE) + ggtitle("WNN clustering")
ggsave(plot=WNNumap, filename=file.path(out, "umap_anno_WNNcluster.pdf"))



#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/wnn_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(umapSample, tsneSample, WNNumap, file=paste0(out,"/wnn.RData"))

