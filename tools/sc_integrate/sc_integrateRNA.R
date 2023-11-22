#####################################
##
## What: sc_integrateRNA.R
## Who : Sivarajan Karunanithi
## When: 10.11.2023
##
## Script to integrate multiple RNA datasets or batch-correct based on a batch variable.
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
resultsdir    <- parseArgs(args,"result=")   
batch         <- parseArgs(args,"batch=")   
n_features    <- parseArgs(args,"n_features=")   
out           <- parseArgs(args,"outdir=") # output folder
reduction_type<- parseArgs(args,"rdtype=") # output folder

runstr <- "Rscript sc_integrateRNA.R [projectdir=projectdir]"

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
options(future.globals.maxSize = 12e9)
addTaskCallback(function(...) {set.seed(100);TRUE})

# check parameter
print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("batch:", batch))
print(paste("n_features:", n_features))
print(paste("out:", out))
print(paste("reduction_type:", reduction_type))


# load sobj from previous module
sobj_noBC <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj_noBC) <- "RNA"

# start with a clean sobj_no_rep3 with out normalizations and reductions performed earlier
# But retains the post-filtered genes, cells, etc.
sobj_noBC_clean = DietSeurat(sobj_noBC,
                      assays = c("RNA","ATAC"),
                      counts = TRUE,
                      data = TRUE,
                      scale.data = FALSE,
                      features = NULL, # keeps all
                      dimreducs = NULL, # removes all
                      graphs = NULL, # removes all
                      misc = FALSE
                      )

rm(sobj_noBC)
gc()
print("creating list of objects")
# split the dataset into a list based on sample - as we want to correct for both condition and replicate
sobj.list <- SplitObject(sobj_noBC_clean, split.by = batch)


# Apply sctransform normalization
# Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
# During normalization, we can also remove confounding sources of variation, for example, mitochondrial 
# mapping percentage by the using vars.to.regress parameter.
# Here we apply the updated version (v2) of SCTransform. This update improves speed and memory consumption, 
# the stability of parameter estimates, the identification of variable features, and the the ability to perform 
# downstream differential expression analyses.

# Please note that, as there will be multiple SCT assays now created one for each sobject split based on the "batch" variable
# in several downstream steps one needs to take extra care to use the corresponding PrepSCT functions Eg. PrepSCTFindMarkers

print("Normalizing list of objects")
# normalize and identify variable features for each dataset independently
sobj.list <- lapply(X = sobj.list, FUN = SCTransform, vst.flavor = "v2", verbose = FALSE)

# select features that are repeatedly variable across datasets for integration
sobj.features <- SelectIntegrationFeatures(object.list = sobj.list, nfeatures = as.integer(as.character(n_features)), verbose = FALSE)

# Prepare the objects for integration
# Ensures that the sctransform residuals for the features specified to anchor.features are present in each object in the list. 
# This is necessary because the default behavior of SCTransform is to only store the residuals for the features determined to be variable. 
# Residuals are recomputed for missing features using the stored model parameters via the GetResidual function.
# Subsets the scale.data slot to only contain the residuals for anchor.features for efficiency in downstream processing.
sobj.list <- PrepSCTIntegration(object.list = sobj.list, anchor.features = sobj.features, verbose = FALSE)

# use the features identified above to identify anchors (cells) which will be used to integrate the data next.
# We should save the identified anchors for integrating the ATAC data later
sobj.anchors <- FindIntegrationAnchors(object.list = sobj.list, normalization.method = "SCT", anchor.features = sobj.features, reduction = reduction_type, dims=1:50, verbose = FALSE)

# integrate - essentially batch correct.
# please keep in mind that the results of batch-correction or integration is supposed to be used only to define the clusters. ALL downstream steps like DE analysis should use the RNA or SCT slots of the seurat object
sobj <- IntegrateData(anchorset = sobj.anchors, normalization.method = "SCT", verbose = FALSE)

# note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(sobj) <- "integrated"

sobj <- RunPCA(sobj, assay="integrated", npcs = 50, reduction.name = 'pca.corrected', reduction.key = 'rnaCorrectedPC_')
sobj <- RunUMAP(sobj, assay="integrated", dims = 1:50, reduction = "pca.corrected", reduction.name = 'umap.corrected.rna', reduction.key = 'rnaCorrectedUMAP_')

# By default, cells are colored by their identity class (can be changed with the group.by parameter). 
pca_correctedrna <- DimPlot(sobj, label = TRUE, repel = TRUE, reduction = "pca.corrected", group.by = "sample") + ggtitle("Corrected RNA - PCA")
umap_correctedrna <- DimPlot(sobj, label = TRUE, repel = TRUE, reduction = "umap.corrected.rna", group.by = "sample") + ggtitle("Corrected RNA - UMAP")

batch_pca_correctedrna <- DimPlot(sobj, label = TRUE, repel = TRUE, reduction = "pca.corrected", split.by = batch) + ggtitle("Corrected RNA - PCA")
batch_umap_correctedrna <- DimPlot(sobj, label = TRUE, repel = TRUE, reduction = "umap.corrected.rna", split.by = batch) + ggtitle("Corrected RNA - UMAP")

ggsave(plot=pca_correctedrna, filename=file.path(out, "pca_correctedrna_anno_sample.pdf"))
ggsave(plot=umap_correctedrna, filename=file.path(out, "umap_correctedrna_anno_sample.pdf"))

ggsave(plot=batch_pca_correctedrna, filename=file.path(out, "pca_correctedrna_anno_batch.pdf"))
ggsave(plot=batch_umap_correctedrna, filename=file.path(out, "umap_correctedrna_anno_batch.pdf"))

### Data clustering
sobj <- FindNeighbors(object = sobj, assay="integrated", reduction = 'pca.corrected', dims = 1:50)
sobj <- FindClusters(object = sobj, verbose = FALSE, algorithm = 1)

# Save clusters in a different metadata column, because it gets overwritten every time FindClusters is called
sobj$clusters_rna_corrected <- sobj$seurat_clusters

Idents(sobj) <- "clusters_rna_corrected"
rnacorrectedUmapCluster <- DimPlot(object = sobj, label = TRUE, reduction = "umap.corrected.rna") +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("Corrected RNA UMAP by cluster") 
ggsave(plot=rnacorrectedUmapCluster, filename = file.path(out, "correctedrna_umap_plot_clusters.pdf"))
ggsave(plot=rnacorrectedUmapCluster, filename = file.path(out, "correctedrna_umap_plot_clusters.png"))

Idents(sobj) <- "sample"
rnacorrectedUmapSample <- DimPlot(object = sobj, label = FALSE, reduction = "umap.corrected.rna") +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("Corrected RNA UMAP by sample")
ggsave(plot=rnacorrectedUmapSample, filename = file.path(out, "correctedrna_umap_plot_samples.pdf"))
ggsave(plot=rnacorrectedUmapSample, filename = file.path(out, "correctedrna_umap_plot_samples.png"))



#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/sc_integrateRNA_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(pca_correctedrna, umap_correctedrna, batch_pca_correctedrna, batch_umap_correctedrna, rnacorrectedUmapCluster, rnacorrectedUmapSample, file=paste0(out,"/sc_integrateRNA.RData"))
save(sobj.anchors, sobj.features, file=paste0(out,"/sc_integrateRNA_anchors.RData"))
