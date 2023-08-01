#####################################
##
## What: CTannoSeurat.R
## Who : Frank Rühle
## When: 05.06.2023
##
## Script to annotate cell types by transferring cell labels from an existing reference dataset with Seurat.
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
run_custom_code <- function(x) {
  eval(parse(text=x))
}

args <- commandArgs(T)
projectdir    <- parseArgs(args,"project=") 
resultsdir    <- parseArgs(args,"res=")   
out           <- parseArgs(args,"outdir=") # output folder
pathRefDataset       <- parseArgs(args,"pathRefDataset=")
columnNameCelltypes  <- parseArgs(args,"columnNameCelltypes=")
assay2use     <- parseArgs(args,"assay=")   
norm_method   <- parseArgs(args,"norm_method=")  
dimReduction  <- parseArgs(args,"dimReduction=")
project_query <- parseArgs(args,"project_query=", default=F, convert="as.logical")
features2use  <- parseArgs(args,"features2use=", convert="run_custom_code")
l2_norm       <- parseArgs(args,"l2_norm=", default=F, convert="as.logical")
k_anchor      <- parseArgs(args,"k_anchor=", convert="as.numeric")
k_filter      <- parseArgs(args,"k_filter=", convert="as.numeric") 
k_score       <- parseArgs(args,"k_score=", convert="as.numeric")
if(length(features2use)==1 && is.na(features2use)) {features2use <- NULL}

runstr <- "Rscript CTannoSeurat.R [projectdir=projectdir]"

# load R environment
renv::use(lockfile=file.path(projectdir, "NGSpipe2go/tools/CTanno/renv.lock"))
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
print(paste("pathRefDataset:", pathRefDataset))
print(paste("columnNameCelltypes:", columnNameCelltypes))
print(paste("assay2use:", assay2use))
print(paste("norm_method:", norm_method))
print(paste("dimReduction:", dimReduction))
print(paste("project_query:", project_query))
print(paste("features2use:", paste(features2use, collapse=" ")))
print(paste("l2_norm:", l2_norm))
print(paste("k_anchor:", k_anchor))
print(paste("k_filter:", k_filter))
print(paste("k_score:", k_score))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- assay2use


# load reference data (must be SeuratObject)
cat("\nLoad reference dataset from", pathRefDataset, "with", columnNameCelltypes, "as cell type annotation column name.\n")
ref <- readRDS(file = file.path(pathRefDataset))
if(!"SCT" %in% names(ref) & assay2use == "SCT") {
  ref <- SCTransform(ref, assay = "RNA")
}
DefaultAssay(ref) <- assay2use


# Find a set of anchors between a reference and query object. 
# These anchors can later be used to transfer data from the reference to query object using the TransferData object. 
transfer_anchors <- FindTransferAnchors(
  reference = ref,
  query = sobj,
  normalization.method = norm_method, # "LogNormalize" or "SCT"
  reduction = dimReduction, # "pcaproject", "lsiproject", "rpca", "cca"
  reference.reduction = NULL,  # name of dim reduction to use from the ref if running pcaproject workflow. If NULL (default), use a PCA computed on ref object.
  project.query = project_query,
  features = features2use,
  l2.norm = l2_norm,
  dims = 1:30,
  k.anchor = k_anchor,
  k.filter = k_filter,
  k.score = k_score
  )

# Transfer categorical or continuous data across single-cell datasets.
predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = ref[[]][columnNameCelltypes][,],
  weight.reduction = dimReduction # "pcaproject", "lsiproject", "rpca", "cca" or custom DimReduc
)

sobj <- AddMetaData(
  object = sobj,
  metadata = predictions
) # cell annotation is stored in column 'predicted.id'


WNNumap_Seurat_CT <- DimPlot(sobj, reduction = "umap.wnn", label = TRUE, repel = TRUE, group.by = 'predicted.id') + ggtitle("Celltype annotation")
ggsave(plot=WNNumap_Seurat_CT, filename=file.path(out, "umap_anno_Seurat_celltypeanno.pdf"))
ggsave(plot=WNNumap_Seurat_CT, filename=file.path(out, "umap_anno_Seurat_celltypeanno.png"))

WNNumap_Seurat_CT_split_by_group <- DimPlot(sobj, reduction = "umap.wnn", label = TRUE, repel = TRUE, split.by="group", group.by = 'predicted.id') + 
  ggtitle("Celltype annotationg")
ggsave(plot=WNNumap_Seurat_CT_split_by_group, width = 10, height = 8, filename = file.path(out, "umap_anno_Seurat_celltypeanno_split_by_group.pdf"))
ggsave(plot=WNNumap_Seurat_CT_split_by_group, width = 10, height = 8, filename = file.path(out, "umap_anno_Seurat_celltypeanno_split_by_group.png"))


# write overview of cell annotations
metadat <- sobj[[]]
annotable <- as.data.frame.matrix(table(metadat$sample, metadat$predicted.id))
write.table(data.frame(sample=rownames(annotable), annotable), file=paste0(out, "/anno_seurat_overview.txt"), quote=F, row.names = F, sep="\t")
# store metadata
write.table(sobj[[]][, !sapply(sobj[[]], is.numeric)], file=paste0(out, "/anno_seurat_cells.csv"), sep=",", quote = F, row.names = F)


#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/CTannoSeurat_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(WNNumap_Seurat_CT, transfer_anchors, predictions, file=paste0(out,"/CTannoSeurat.RData"))

