#####################################
##
## What: diffPeaks.R
## Who : Frank Rühle
## When: 06.06.2023
##
## Script to identify differentially accessible peaks between all pairs of groups for each cluster and celltype.
##
## Args:
## -----
## targets=targets.txt      # file describing the targets 
##
##
######################################

renv::use(lockfile='NGSpipe2go/tools/sc_DNAaccess/renv.lock')
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
library(openxlsx)

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
minCells      <- parseArgs(args,"minCells=")
CTannoSelected <- parseArgs(args,"CTannoSelected=")

runstr <- "Rscript diffPeaks.R [projectdir=projectdir]"

print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))
print(paste("minCells:", minCells))
print(paste("CTannoSelected:", CTannoSelected))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "ATAC"


# differential accessibility between pairs of samples per calculated per cluster

sobj$label <- sobj$Group_ATACClusters
Idents(sobj) <- "label"

grp_counts <- table(sobj$Group_ATACClusters)

sample_pairs <- combn(sort(unique(sobj$group)), 2)
da_groups_atac <- list()
for (cl in sort(unique(sobj$clusters_atac))) {
  for (sp in 1:ncol(sample_pairs)) {
    grp1 <- paste(sample_pairs[1, sp], cl, sep="_")
    grp2 <- paste(sample_pairs[2, sp], cl, sep="_")
    if (grp1 %in% names(grp_counts) && grp_counts[grp1] >= minCells && 
        grp2 %in% names(grp_counts) && grp_counts[grp2] >= minCells) {
      group_pair <- paste(grp1, grp2, sep = ".")
      cat("Processing pair: ", group_pair)
      da_groups_atac[[group_pair]] <- FindMarkers(
        object = sobj,
        ident.1 = grp1,
        ident.2 = grp2,
        min.pct = 0.05,
        test.use = 'LR',
        latent.vars = 'nFeature_ATAC'
      )
      cat(paste0(" (contains ", nrow(da_groups_atac[[group_pair]][da_groups_atac[[group_pair]]$p_val_adj < 0.05, ]), " peaks with adj p-val < 0.05)\n"))
    }
  }
}

# Save results
if(length(da_groups_atac)>0) {
  openxlsx::write.xlsx(da_groups_atac, file = file.path(out, "peaks_diffacces_clusters.xlsx"), rowNames=T)
}

gc()



## differential accessibility between pairs of samples per calculated per celltype

# select the cell type annotation to be used in downstream analysis
if(!is.null(CTannoSelected) && !is.na(CTannoSelected)) {
  switch(CTannoSelected,
         Seurat={ celltype <- "predicted.id" },
         Marker={ celltype <- "CTAnnotationSCType" },
         stop(c("Don't find cell type annotation"))   
  )
  cat(paste0("We used the cell type annotation option '", CTannoSelected, "' for downstream analysis (i.e. annotation column '", celltype, "').\n"))  
  
} else {
  cat("no cell type annotation selected")
}

sobj$Group_Celltypes <- paste(sobj$group, sobj[[celltype]][,1], sep = "_")
sobj$label <- sobj$Group_Celltypes
Idents(sobj) <- "label"

grp_counts <- table(sobj$Group_Celltypes)

sample_pairs <- combn(sort(unique(sobj$group)), 2)
da_groups_atac_ct <- list()
for (cl in sort(unique(sobj[[celltype]][,1]))) {
  for (sp in 1:ncol(sample_pairs)) {
    grp1 <- paste(sample_pairs[1, sp], cl, sep="_")
    grp2 <- paste(sample_pairs[2, sp], cl, sep="_")
    if (grp1 %in% names(grp_counts) && grp_counts[grp1] >= minCells && 
        grp2 %in% names(grp_counts) && grp_counts[grp2] >= minCells) {
      group_pair <- paste(grp1, grp2, sep = ".")
      cat("Processing pair: ", group_pair)
      da_groups_atac_ct[[group_pair]] <- FindMarkers(
        object = sobj,
        ident.1 = grp1,
        ident.2 = grp2,
        min.pct = 0.05,
        test.use = 'LR',
        latent.vars = 'nFeature_ATAC'
      )
      cat(paste0(" (contains ", nrow(da_groups_atac_ct[[group_pair]][da_groups_atac_ct[[group_pair]]$p_val_adj < 0.05, ]), " peaks with adj p-val < 0.05)\n"))
    }
  }
}

# Save results
library(openxlsx)
if(length(da_groups_atac_ct)>0) {
  # sheetNames too long! Max length is 31 characters.
  names(da_groups_atac_ct) <- gsub("Proliferating cell", "Prolif", names(da_groups_atac_ct)) 
  names(da_groups_atac_ct) <- gsub("Oligodendrocyte", "ODC", names(da_groups_atac_ct))
  names(da_groups_atac_ct) <- gsub("Proliferating radial glia", "PRG", names(da_groups_atac_ct))
  names(da_groups_atac_ct) <- gsub("Radial glia", "Rad glia", names(da_groups_atac_ct))
  
  openxlsx::write.xlsx(da_groups_atac_ct, file = file.path(out, "peaks_diffacces_celltypes.xlsx"), rowNames=T)
}

gc()

#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/diffPeaks_info.txt"))
#saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(da_groups_atac, da_groups_atac_ct, file=paste0(out,"/diffPeaks.RData"))

