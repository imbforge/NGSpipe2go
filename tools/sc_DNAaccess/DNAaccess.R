#####################################
##
## What: DNAaccess.R
## Who : Frank Rühle
## When: 02.06.2023
##
## Script to process the DNA accessibility assay by performing latent semantic indexing (LSI).
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
projectdir       <- parseArgs(args,"project=") 
resultsdir       <- parseArgs(args,"res=")   
out              <- parseArgs(args,"outdir=") # output folder
featureCutoff    <- parseArgs(args,"featureCutoff=")
skipFirstLSIcomp <- parseArgs(args,"skipFirstLSIcomp=", default=0, convert="as.numeric")

# featureCutoff for feature to be included in the VariableFeatures: either percentile specified as 'q' followed by the minimum percentile 
# or integer specifying the  for the feature to be included in the set of VariableFeatures.
if(!is.na(as.numeric(featureCutoff))) {featureCutoff <- as.numeric(featureCutoff)} 

runstr <- "Rscript DNAaccess.R [projectdir=projectdir]"

# load R environment
renv::use(lockfile=file.path(projectdir, "NGSpipe2go/tools/sc_DNAaccess/renv.lock"))
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
print(paste("featureCutoff:", featureCutoff))
print(paste("skipFirstLSIcomp :", skipFirstLSIcomp ))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "ATAC"

sobj <- RunTFIDF(sobj)                        # normalization
sobj <- FindTopFeatures(sobj, min.cutoff = featureCutoff) # feature selection (default 'q5')
sobj <- RunSVD(sobj)                          # dimensionality reduction (default n = 50)


depthCorPlot <- DepthCor(sobj) +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  theme(plot.subtitle=element_text(hjust=0.5)) +
  NoLegend()
ggsave(plot=depthCorPlot, filename = file.path(out, "atac_lsi_depth_corr_plot.pdf"))
ggsave(plot=depthCorPlot, filename = file.path(out, "atac_lsi_depth_corr_plot.png"))

sobj <- RunUMAP(sobj, reduction = 'lsi', dims = (1+skipFirstLSIcomp):50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
sobj <- RunTSNE(sobj, reduction = 'lsi', dims = (1+skipFirstLSIcomp):50, reduction.name = "tsne.atac", reduction.key = "atacTSNE_")



### ATAC data clustering
sobj <- FindNeighbors(object = sobj, reduction = 'lsi', dims = (1+skipFirstLSIcomp):50)
sobj <- FindClusters(object = sobj, verbose = FALSE, algorithm = 3)

# Save clusters in a different metadata column, because it gets overwritten every time FindClusters is called
sobj$clusters_atac <- sobj$seurat_clusters
# Also annotate nuclei with per-sample clusters
sobj$Sample_ATACClusters <- paste(sobj$sample, sobj$clusters_atac, sep = "_")
sobj$Group_ATACClusters <- paste(sobj$group, sobj$clusters_atac, sep = "_")


Idents(sobj) <- "clusters_atac"

atacUmapCluster <- DimPlot(object = sobj, label = TRUE, reduction = "umap.atac") +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("ATAC UMAP by cluster") 
ggsave(plot=atacUmapCluster, filename = file.path(out, "atac_umap_plot_clusters.pdf"))
ggsave(plot=atacUmapCluster, filename = file.path(out, "atac_umap_plot_clusters.png"))

Idents(sobj) <- "sample"

atacUmapSample <- DimPlot(object = sobj, label = FALSE, reduction = "umap.atac") +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("ATAC UMAP by sample")
ggsave(plot=atacUmapSample, filename = file.path(out, "atac_umap_plot_samples.pdf"))
ggsave(plot=atacUmapSample, filename = file.path(out, "atac_umap_plot_samples.png"))





#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/DNAaccess_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(depthCorPlot, atacUmapCluster, atacUmapSample, file=paste0(out,"/DNAaccess.RData"))

