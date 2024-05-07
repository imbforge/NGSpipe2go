#####################################
##
## What: sc_integrateATAC.R
## Who : Sivarajan Karunanithi
## When: 13.11.2023
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
resultsdir    <- parseArgs(args,"res=")   
out           <- parseArgs(args,"outdir=") # output folder
featureCutoff    <- parseArgs(args,"featureCutoff=")
skipFirstLSIcomp <- parseArgs(args,"skipFirstLSIcomp=", default=0, convert="as.numeric")

runstr <- "Rscript sc_integrateATAC.R [projectdir=projectdir]"

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

# TODO: Have a check if integrateRNA was done - if not, run ATAC integration independently, which maybe useful for scATAC pipeline
# load the anchors we created in the RNAstep
if(!dir.exists(file.path(resultsdir, "sc_integrateRNA")))   stop(paste("Directory",file.path(resultsdir, "sc_integrateRNA"),"does NOT exist. Run with:\n",runstr))
load(file.path(resultsdir, "sc_integrateRNA" ,"sc_integrateRNA_anchors.RData"))

# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "ATAC"

sobj <- RunTFIDF(sobj)                        # normalization
sobj <- FindTopFeatures(sobj, min.cutoff = featureCutoff) # feature selection (default 'q5')
sobj <- RunSVD(sobj)                          # dimensionality reduction (default n = 50)

# we need a new object, otherwise the default integrated (if RNA was done before) will be overwritten
# integrate LSI embeddings
integrated_atac <- IntegrateEmbeddings(
  anchorset = sobj.anchors, # loaded from the integrateRNA step
  reductions = sobj[["lsi"]],
  new.reduction.name = "lsi.corrected",
  dims.to.integrate = 1:50
)

integrated_atac <- RunUMAP(integrated_atac, reduction = "lsi.corrected", dims = (1+skipFirstLSIcomp):50, reduction.name = 'umap.corrected.atac', reduction.key = "atacCorrectedUMAP_")
integrated_atac <- RunTSNE(integrated_atac, reduction = "lsi.corrected", dims = (1+skipFirstLSIcomp):50, reduction.name = "tsne.corrected.atac", reduction.key = "atacCorrectedTSNE_")

# copy the corrected atac embeddings to the common object
sobj@reductions$lsi.corrected = integrated_atac@reductions$lsi.corrected
sobj@reductions$umap.corrected.atac = integrated_atac@reductions$umap.corrected.atac
sobj@reductions$tsne.corrected.atac = integrated_atac@reductions$tsne.corrected.atac

### ATAC data clustering
sobj <- FindNeighbors(object = sobj, reduction = 'lsi.corrected', dims = (1+skipFirstLSIcomp):50)
sobj <- FindClusters(object = sobj, verbose = FALSE, algorithm = 3)

# Save clusters in a different metadata column, because it gets overwritten every time FindClusters is called
sobj$clusters_atac_corrected <- sobj$seurat_clusters
# Also annotate nuclei with per-sample clusters
sobj$Sample_ATACClusters_corrected <- paste(sobj$sample, sobj$clusters_atac_corrected, sep = "_")
sobj$Group_ATACClusters_corrected <- paste(sobj$group, sobj$clusters_atac_corrected, sep = "_")


Idents(sobj) <- "clusters_atac_corrected"

correctedatacUmapCluster <- DimPlot(object = sobj, label = TRUE, reduction = "umap.corrected.atac") +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("Corrected ATAC UMAP by cluster")
ggsave(plot=correctedatacUmapCluster, filename = file.path(out, "correctedatac_umap_plot_clusters.pdf"))
ggsave(plot=correctedatacUmapCluster, filename = file.path(out, "correctedatac_umap_plot_clusters.png"))

Idents(sobj) <- "sample"

correctedatacUmapSample <- DimPlot(object = sobj, label = FALSE, reduction = "umap.corrected.atac") +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("Corrected ATAC UMAP by sample")
ggsave(plot=correctedatacUmapSample, filename = file.path(out, "correctedatac_umap_plot_samples.pdf"))
ggsave(plot=correctedatacUmapSample, filename = file.path(out, "correctedatac_umap_plot_samples.png"))

#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/sc_integrateATAC_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(correctedatacUmapCluster, correctedatacUmapSample, file=paste0(out,"/sc_integrateATAC.RData"))
