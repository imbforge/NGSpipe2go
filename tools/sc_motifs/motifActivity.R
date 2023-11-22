#####################################
##
## What: motifActivity.R
## Who : Frank Rühle
## When: 01.08.2023
##
## Script to determine differential activity scores between for each cluster or cell type, respectively, compared to all other cluster and celltypes
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
db            <- parseArgs(args,"db=")     
clusterVar    <- parseArgs(args,"clusterVar=")     
CTannoSelected <- parseArgs(args,"CTannoSelected=")     
motif2plot   <- parseArgs(args,"motif2plot=")     

runstr <- "Rscript motifActivity.R [projectdir=projectdir]"

# load R environment
renv::use(lockfile=file.path(projectdir, "NGSpipe2go/tools/sc_motifs/renv.lock"))
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
library(JASPAR2020)

# set options
options(stringsAsFactors=FALSE)
addTaskCallback(function(...) {set.seed(100);TRUE})

# check parameter
print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))
print(paste("db:", db))
print(paste("clusterVar:", clusterVar))
print(paste("CTannoSelected:", CTannoSelected))
print(paste("motif2plot:", motif2plot))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "ATAC"
load(file.path(resultsdir, "/gtf.RData"))


# load relevant BSgenome package (needed by Signac for motif analysis)
switch(db,
       hg38={ failed_BSgenome <- library("BSgenome.Hsapiens.UCSC.hg38") 
       BSgenome <- BSgenome.Hsapiens.UCSC.hg38
       JasparTaxGroup <- "vertebrates"
       },
       mm10={ failed_BSgenome <- library("BSgenome.Mmusculus.UCSC.mm10") 
       BSgenome <- BSgenome.Mmusculus.UCSC.mm10
       JasparTaxGroup <- "vertebrates"
       },
       stop(c("Don't find genome:", db))   
)

# subset BSgenome object for standard chromosomes
keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}
BSgenome_standard <-  keepBSgenomeSequences(BSgenome, standardChromosomes(BSgenome))


# sobj has data mapped to scaffolds named differently to the BSgenome. 
# The scaffolds aren't that useful for analysis, so you could solve this by removing peaks that are on scaffolds.
main.chroms <- standardChromosomes(BSgenome)
keep.peaks <- which(as.character(seqnames(granges(sobj))) %in% main.chroms)
sobj[["ATAC"]] <- subset(sobj[["ATAC"]], features = rownames(sobj[["ATAC"]])[keep.peaks])

# same is done for the expression assays
keep.genes <- which(rownames(sobj[["RNA"]]) %in% names(gtf))
sobj[["RNA"]] <- subset(sobj[["RNA"]], features = rownames(sobj[["RNA"]])[keep.genes])
keep.genes2 <- which(rownames(sobj[["SCT"]]) %in% names(gtf))
sobj[["SCT"]] <- subset(sobj[["SCT"]], features = rownames(sobj[["SCT"]])[keep.genes2])


# create motif object if not done before
if(is.null(attr(sobj[['ATAC']], "motifs")) || is.null(sobj[['ATAC']]@motifs)) {
  cat("\nCreate motifs object with JASPAR2020.\n")
  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- TFBSTools::getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = JasparTaxGroup, all_versions = FALSE)
  )
  
  # add motif information
  # Construct a Motif object containing DNA sequence motif information and add it to an existing Seurat object or ChromatinAssay. 
  # If running on a Seurat object, AddMotifs will also run RegionStats to compute the GC content of each peak and store the results in the feature metadata. 
  sobj <- AddMotifs(
    object = sobj,
    genome = BSgenome_standard,     # org-specific
    pfm = pfm,
    assay = "ATAC"
  )
} else {
  cat("\nuse existing motifs object.\n")
}



### differential motif activity analysis between groups of cells with ChromVar

# We can also compute a per-cell motif activity score by running chromVAR. This allows us to visualize motif activities per cell, 
# and also provides an alternative method of identifying differentially-active motifs between cell types. 
# ChromVAR identifies motifs associated with variability in chromatin accessibility between cells. 

sobj <- RunChromVAR(
  object = sobj,
  genome = BSgenome_standard
)

DefaultAssay(sobj) <- 'chromvar'

# plot activity of certain motif
if(!is.na(motif2plot)) {
  if(file.exists((motif2plot))) {
    m2plot <- read.delim(motif2plot, header = F, sep = "\t", comment.char = "#")[,1]
  } else {
    m2plot <- motif2plot
  }
  cat("\nmotif(s) to plot:", paste(m2plot, collapse=", "), "\n")
  
  for (m in m2plot) {
    motifActPlot <- FeaturePlot(
      object = sobj,
      reduction = "umap.wnn",
      features = m,
      min.cutoff = 'q5',
      max.cutoff = 'q95',
      pt.size = 0.1
    )
    ggsave(plot=motifActPlot, filename=file.path(out, paste0("Motif_activity_plot_", m, ".pdf")), width = 7, height = 7)
    ggsave(plot=motifActPlot, filename=file.path(out, paste0("Motif_activity_plot_", m, ".png")), width = 7, height = 7)
    }
}

# We can also directly test for differential activity scores between clusters or cell types. 
# This tends to give similar results as performing an enrichment test on differentially accessible peaks between the cell types (shown above).
# When performing differential testing on the chromVAR z-score, we can set mean.fxn=rowMeans and fc.name="avg_diff" 
# in the FindMarkers() function so that the fold-change calculation computes the average difference in z-score between the groups.

# extract motif object
motifobj <- SeuratObject::GetAssayData(object = sobj[['ATAC']], slot = "motifs")
#lgCMatrixObj <- GetMotifData(object = motifobj)
# convert ID to common name
ids <- colnames(motifobj)
names <- ConvertMotifID(object = motifobj, id = ids)
motifnames <- data.frame(id=ids, name=names)
write.table(motifnames, file =file.path(out, "motifnames_overview.txt"), quote = F, row.names = F,  sep="\t")


## differential activity per cluster
sobj$label <- sobj[[clusterVar]][,1]
Idents(sobj) <- "label"
diffAct_cluster <- list()
for (c in sort(unique(Idents(sobj)))) {
 if(nrow(sobj[[]][sobj[[]][,clusterVar]==c,])>=3) {
  diffAct_cluster[[c]] <- FindMarkers(
    object = sobj,
    ident.1 = c,
    ident.2 = NULL,
    only.pos = TRUE,
    mean.fxn = rowMeans,
    fc.name = "avg_diff"
  )
  diffAct_cluster[[c]]$motif_name <- ConvertMotifID(object = motifobj, id = rownames(diffAct_cluster[[c]]))
 }
}
openxlsx::write.xlsx(diffAct_cluster, file = file.path(out, "ChromVar_motif_enriched_cluster_vs_all.xlsx"), rowNames=T)


## differential activity per celltype
# select the cell type annotation to be used in downstream analysis
diffAct_ct <- NULL
if(!is.null(CTannoSelected) && !is.na(CTannoSelected)) {
  switch(CTannoSelected,
         Seurat={ celltype <- "predicted.id" },
         Marker={ celltype <- "CTAnnotationSCType" },
         stop(c("Don't find cell type annotation"))   
  )
  cat(paste0("We used the cell type annotation option '", CTannoSelected, "' for downstream analysis (i.e. annotation column '", celltype, "').\n"))  

  sobj$label <- sobj[[celltype]][,1]
  Idents(sobj) <- "label"
  diffAct_ct <- list()
  for (c in sort(unique(Idents(sobj)))) {
   if(nrow(sobj[[]][sobj[[]][,celltype]==c,])>=3) {
    diffAct_ct[[c]] <- FindMarkers(
      object = sobj,
      ident.1 = c,
      ident.2 = NULL,
      only.pos = TRUE,
      mean.fxn = rowMeans,
      fc.name = "avg_diff"
    )
    diffAct_ct[[c]]$motif_name <- ConvertMotifID(object = motifobj, id = rownames(diffAct_ct[[c]]))
   }
  }
  openxlsx::write.xlsx(diffAct_ct, file = file.path(out, "ChromVar_motif_enriched_celltype_vs_all.xlsx"), rowNames=T)
}


#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/motifActivity_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(diffAct_cluster, diffAct_ct, file=paste0(out,"/motifActivity.RData"))

