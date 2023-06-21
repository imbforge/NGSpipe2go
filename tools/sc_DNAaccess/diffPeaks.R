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
assay2use     <- parseArgs(args,"assay=")   
out           <- parseArgs(args,"outdir=") # output folder
minCells      <- parseArgs(args,"minCells=", convert="as.numeric")
clusterVar    <- parseArgs(args,"clusterVar=")
CTannoSelected <- parseArgs(args,"CTannoSelected=")
test2use      <- parseArgs(args,"test=")
latentVars    <- parseArgs(args,"latentVars=")
if(!test2use %in% c("LR", "negbinom", "poisson", "MAST")) {latentVars <- NULL}

runstr <- "Rscript diffPeaks.R [projectdir=projectdir]"

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
library(openxlsx)

# set options
options(stringsAsFactors=FALSE)

# check parameter
print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("assay2use:", assay2use))
print(paste("out:", out))
print(paste("minCells:", minCells))
print(paste("clusterVar:", clusterVar))
print(paste("CTannoSelected:", CTannoSelected))
print(paste("test2use:", test2use))
print(paste("latentVars:", latentVars))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- assay2use


# differential accessibility between pairs of samples per calculated per cluster
sobj$Group_Cluster <- paste(sobj$group, sobj[[clusterVar]][,1], sep = "_")
sobj$label <- sobj$Group_Cluster
Idents(sobj) <- "label"

grp_counts <- table(sobj$Group_Cluster)

sample_pairs <- combn(sort(unique(sobj$group)), 2)
da_groups_atac <- list()
for (cl in sort(unique(sobj[[clusterVar]][,1]))) {
  for (sp in 1:ncol(sample_pairs)) {
    grp1 <- paste(sample_pairs[1, sp], cl, sep="_")
    grp2 <- paste(sample_pairs[2, sp], cl, sep="_")
    if (grp1 %in% names(grp_counts) && grp_counts[grp1] >= minCells && 
        grp2 %in% names(grp_counts) && grp_counts[grp2] >= minCells) {
      group_pair <- paste(grp1, grp2, sep = ".")
      cat("Processing pair: ", group_pair)
      set.seed(100)
      da_groups_atac[[group_pair]] <- FindMarkers(
        object = sobj,
        assay = assay2use,
        ident.1 = grp1,
        ident.2 = grp2,
        min.pct = 0.05,
        test.use = test2use,
        latent.vars = latentVars
      )
      cat(paste0(" (contains ", nrow(da_groups_atac[[group_pair]][da_groups_atac[[group_pair]]$p_val_adj < 0.05, ]), " peaks with adj p-val < 0.05)\n"))
    }
  }
}

# Save results
if(length(da_groups_atac)>0) {
  names(da_groups_atac) <- substr(names(da_groups_atac), 1, apply(data.frame(nchar(names(da_groups_atac)),31), 1, min))
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
        set.seed(100)
        da_groups_atac_ct[[group_pair]] <- FindMarkers(
          object = sobj,
          assay = assay2use,
          ident.1 = grp1,
          ident.2 = grp2,
          min.pct = 0.05,
          test.use = test2use,
          latent.vars = latentVars
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
    names(da_groups_atac_ct) <- gsub("Proliferating radial glia", "PRG", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("GABAergic neurons", "GABA neur", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("Neural Progenitor cells", "n Prog", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("Neuroepithelial cells", "Neuroepith", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("Oligodendrocyte precursor cells", "ODC prec", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("Dopaminergic neurons", "Dop neur", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("Glutamatergic neurons", "Glutamat n", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("Mature neurons", "Mature n", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("Neuroblasts", "Neuroblast", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("Non myelinating Schwann cells", "NMSC", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("Radial glial cells", "Rad glia", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("Radial glia", "Rad glia", names(da_groups_atac_ct))
    names(da_groups_atac_ct) <- gsub("Oligodendrocyte", "ODC", names(da_groups_atac_ct))
    
    names(da_groups_atac_ct) <- substr(names(da_groups_atac_ct), 1, apply(data.frame(nchar(names(da_groups_atac_ct)),31), 1, min))
    openxlsx::write.xlsx(da_groups_atac_ct, file = file.path(out, "peaks_diffacces_celltypes.xlsx"), rowNames=T)
  }
  
  gc()

} else {
  cat("no cell type annotation selected")
}


#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/diffPeaks_info.txt"))
#saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(da_groups_atac, da_groups_atac_ct, file=paste0(out,"/diffPeaks.RData"))

