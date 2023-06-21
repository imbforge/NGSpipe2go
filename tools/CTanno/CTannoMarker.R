#####################################
##
## What: CTannoMarker.R
## Who : Sivarajan Karunanithi, Frank Rühle
## When: 05.06.2023
##
## Script to annotate cell types by marker genes with scType.
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
assay2use     <- parseArgs(args,"assay=")   
clusterVar    <- parseArgs(args,"clusterVar=")   
dbfile        <- parseArgs(args,"dbfile=")
tissue        <- parseArgs(args,"tissue=")

runstr <- "Rscript CTannoMarker.R [projectdir=projectdir]"

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
library(HGNChelper)

# set options
options(stringsAsFactors=FALSE)

# check parameter
print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))
print(paste("assay2use:", assay2use))
print(paste("clusterVar:", clusterVar))
print(paste("dbfile:", dbfile))
print(paste("tissue:", tissue))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- assay2use

set.seed(100)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file 
# TODO: This file can be provided in custom format, if we need to provide a different list of marker genes
# DB file should contain four columns (tissueType - tissue type, cellName - cell type, geneSymbolmore1 - positive marker genes, geneSymbolmore2 - marker genes not expected to be expressed by a cell type)

# TODO: Set this as a variable - an autodetect function is also available, perhaps it is safer to name it


# prepare gene sets
gs_list = gene_sets_prepare(dbfile, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = sobj[[assay2use]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(sobj@meta.data[,clusterVar]), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(sobj@meta.data[sobj@meta.data[,clusterVar]==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sobj@meta.data[,clusterVar]==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

sobj@meta.data$CTAnnotationSCType = "" # cell annotation is stored in column 'CTAnnotationSCType'
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sobj@meta.data$CTAnnotationSCType[sobj@meta.data[,clusterVar] == j] = as.character(cl_type$type[1])
}

WNNumap_scType_celltypeannot <- DimPlot(sobj, reduction = "umap.wnn", label = TRUE, repel = TRUE, group.by = 'CTAnnotationSCType') + ggtitle("Celltype anno based on WNN clustering")
ggsave(plot=WNNumap_scType_celltypeannot, filename=file.path(out, "umap_anno_scType_WNNcluster_celltypeanno.pdf"))
ggsave(plot=WNNumap_scType_celltypeannot, filename=file.path(out, "umap_anno_scType_WNNcluster_celltypeanno.png"))


WNNumap_scType_celltypeannot_split_by_group <- DimPlot(sobj, reduction = "umap.wnn", label = TRUE, repel = TRUE, split.by="group", group.by = 'CTAnnotationSCType') + 
  ggtitle("Celltype anno based on WNN clustering")
ggsave(plot=WNNumap_scType_celltypeannot_split_by_group, width = 10, height = 8, filename = file.path(out, "umap_anno_scType_WNNcluster_celltypeanno_split_by_group.pdf"))
ggsave(plot=WNNumap_scType_celltypeannot_split_by_group, width = 10, height = 8, filename = file.path(out, "umap_anno_scType_WNNcluster_celltypeanno_split_by_group.png"))



#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/CTannoMarker_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(WNNumap_scType_celltypeannot, sctype_scores, cL_resutls, file=paste0(out,"/CTannoMarker.RData"))

