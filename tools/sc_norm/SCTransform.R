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

runstr <- "Rscript SCTransform.R [projectdir=projectdir]"

# load R environment
renv::use(lockfile=file.path(projectdir, "NGSpipe2go/tools/sc_norm/renv.lock"))
print(.libPaths())

library(tidyverse)
library(Biobase)
library(data.table)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(reshape2)
library(Seurat)
library(Signac)
library(uwot)
library(glmGamPoi)

# set options
options(stringsAsFactors=FALSE)
addTaskCallback(function(...) {set.seed(100);TRUE})

# check parameter
print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "RNA"

# Apply sctransform normalization
# Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
# During normalization, we can also remove confounding sources of variation, for example, mitochondrial 
# mapping percentage by the using vars.to.regress parameter.
# Here we apply the updated version (‘v2’) of SCTransform. This update improves speed and memory consumption, 
# the stability of parameter estimates, the identification of variable features, and the the ability to perform 
# downstream differential expression analyses.
sobj <- SCTransform(sobj, assay = "RNA", new.assay.name = "SCT", vst.flavor = "v2") # 
sobj <- RunPCA(sobj, assay="SCT", dims = 1:50, reduction.name = 'pca.sct', reduction.key = 'PC_')
sobj <- RunUMAP(sobj, assay="SCT", dims = 1:50, reduction = "pca.sct", reduction.name = 'umap.sct', reduction.key = 'sctUMAP_')

# By default, cells are colored by their identity class (can be changed with the group.by parameter). 
pca_sct <- DimPlot(sobj, label = TRUE, repel = TRUE, reduction = "pca.sct", group.by = "sample") + ggtitle("SCT PCA")
umap_sct <- DimPlot(sobj, label = TRUE, repel = TRUE, reduction = "umap.sct", group.by = "sample") + ggtitle("SCT UMAP")

ggsave(plot=pca_sct, filename=file.path(out, "pca_sct_anno_sample.pdf"))
ggsave(plot=umap_sct, filename=file.path(out, "umap_sct_anno_sample.pdf"))


### SCT data clustering
sobj <- FindNeighbors(object = sobj, assay="SCT", reduction = 'pca.sct', dims = 1:50)
sobj <- FindClusters(object = sobj, verbose = FALSE, algorithm = 1)

# Save clusters in a different metadata column, because it gets overwritten every time FindClusters is called
sobj$clusters_sct <- sobj$seurat_clusters

Idents(sobj) <- "clusters_sct"
sctUmapCluster <- DimPlot(object = sobj, label = TRUE, reduction = "umap.sct") +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("SCT UMAP by cluster") 
ggsave(plot=sctUmapCluster, filename = file.path(out, "sct_umap_plot_clusters.pdf"))
ggsave(plot=sctUmapCluster, filename = file.path(out, "sct_umap_plot_clusters.png"))

Idents(sobj) <- "sample"
sctUmapSample <- DimPlot(object = sobj, label = FALSE, reduction = "umap.sct") +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("SCT UMAP by sample")
ggsave(plot=sctUmapSample, filename = file.path(out, "sct_umap_plot_samples.pdf"))
ggsave(plot=sctUmapSample, filename = file.path(out, "sct_umap_plot_samples.png"))



#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/SCTransform_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(pca_sct, umap_sct, sctUmapCluster, sctUmapSample, file=paste0(out,"/SCTransform.RData"))

