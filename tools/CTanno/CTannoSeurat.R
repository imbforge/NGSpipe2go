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
## targets=targets.txt      # file describing the targets 
##
##
######################################

renv::use(lockfile='NGSpipe2go/tools/CTanno/renv.lock')
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
pathRefDataset           <- parseArgs(args,"pathRefDataset=")
columnNameCelltypes      <- parseArgs(args,"columnNameCelltypes=")

runstr <- "Rscript CTannoSeurat.R [projectdir=projectdir]"

print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))
print(paste("pathRefDataset:", pathRefDataset))
print(paste("columnNameCelltypes:", columnNameCelltypes))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "SCT"


# load reference data (must be SeuratObject)
cat("\nLoad reference dataset from", pathRefDataset, "with", columnNameCelltypes, "as cell type annotation column name.\n")
ref <- readRDS(file = file.path(pathRefDataset))
ref <- SCTransform(ref)
DefaultAssay(ref) <- "SCT"


# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = ref,
  query = sobj,
  normalization.method = "SCT",
  #reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:30
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = ref[[]][columnNameCelltypes][,],
  weight.reduction = sobj[['pca']],
  dims = 1:30
)

sobj <- AddMetaData(
  object = sobj,
  metadata = predictions
) # cell annotation is stored in column 'predicted.id'


WNNumap_Seurat_CT <- DimPlot(sobj, reduction = "umap.wnn", label = TRUE, repel = TRUE, group.by = 'predicted.id') + ggtitle("Celltype annotation")
ggsave(plot=WNNumap_Seurat_CT, filename=file.path(out, "umap_anno_Seurat_celltypeanno.pdf"))
ggsave(plot=WNNumap_Seurat_CT, filename=file.path(out, "umap_anno_Seurat_celltypeanno.png"))




#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/CTannoSeurat_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(WNNumap_Seurat_CT, transfer_anchors, predictions, file=paste0(out,"/CTannoSeurat.RData"))

