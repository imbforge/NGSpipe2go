#####################################
##
## What: motifFootprinting.R
## Who : Frank Rühle
## When: 06.06.2023
##
## Script to perform motif footprinting with enriched TF motifs.
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
  
  #if(!is.null(default)) default else do.call(convert,list(NA))
}

args <- commandArgs(T)
projectdir    <- parseArgs(args,"project=") 
resultsdir    <- parseArgs(args,"res=")   
out           <- parseArgs(args,"outdir=") # output folder
db            <- parseArgs(args,"db=")     
motifEnrich_dir <- parseArgs(args,"motifEnrich_dir=")
assay         <- parseArgs(args,"assay=")
motifsPerContrast <- parseArgs(args,"motifsPerContrast=", 0, convert="as.numeric")
motifsByName  <- parseArgs(args,"motifsByName=")
inPeaks       <- parseArgs(args,"inPeaks=", convert="as.logical")
upstream      <- parseArgs(args,"upstream=", 250, convert="as.numeric")
downstream    <- parseArgs(args,"downstream=", 250, convert="as.numeric")
groupPlotBy   <- parseArgs(args,"groupPlotBy=")   
splitPlotBy   <- parseArgs(args,"splitPlotBy=")   


runstr <- "Rscript motifFootprinting.R [projectdir=projectdir]"

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
library(motifmatchr)
library(TFBSTools)


# set options
options(stringsAsFactors=FALSE)
addTaskCallback(function(...) {set.seed(100);TRUE})

# check parameter
print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))
print(paste("db:", db))
print(paste("motifEnrich_dir:", motifEnrich_dir))
print(paste("assay:", assay))
print(paste("motifsPerContrast:", motifsPerContrast))
print(paste("motifsByName:", motifsByName))
print(paste("inPeaks:", inPeaks))
print(paste("upstream:", upstream))
print(paste("downstream:", downstream))


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


# For each enriched motif Now we can footprint any motif that we have positional information for. 
# By default, this includes every instance of the motif in the genome. We can instead use the 
# in.peaks = TRUE parameter to include only those motifs that fall inside a peak in the assay. 
# The Footprint() function gathers all the required data and stores it in the assay. We can then 
# plot the footprinted motifs using the PlotFootprint() function.

motifsPerCluster <- pCl <- NULL
motifsPerCelltype <- pCT <- NULL
motifsByNameLoaded <- pNL <- NULL

if(!is.null(motifEnrich_dir) && !is.na(motifEnrich_dir) && !is.null(motifsPerContrast) && motifsPerContrast>0) {
  cat(paste("\nMotif footprinting for top enriched motifs per cluster and celltype comparison\n"))
  
  load(file.path(motifEnrich_dir, "/motifEnrich.RData"))
  
  motifsPerCluster <- sapply(da_enriched_motifs, function(x) {c(x$motif.name[1:motifsPerContrast])})
  motifsPerCluster <- unique(as.vector(motifsPerCluster))
  
  motifsPerCelltype <- sapply(da_enriched_motifs_ct, function(x) {c(x$motif.name[1:motifsPerContrast])})
  motifsPerCelltype <- unique(as.vector(motifsPerCelltype))

}

if(!is.null(motifsByName) && !is.na(motifsByName) && file.exists(motifsByName)) {
  cat(paste("\nLoad motifs from file\n"))
  motifsByNameLoaded <- read.delim(motifsByName, header=F)
  motifsByNameLoaded <- as.vector(motifsByNameLoaded$V1)
}

motifs2footprint <- unique(c(motifsPerCluster, motifsPerCelltype, motifsByNameLoaded))
  
if (length(motifs2footprint)>0) {

  sobj <- Footprint(
    object = sobj,
    assay = assay,
    motif.name = motifs2footprint,
    genome = BSgenome,
    upstream = upstream,
    downstream = downstream,
    in.peaks = inPeaks
    )
  
  if(!is.null(motifsPerCluster)) {
    out_per_cluster <- file.path(out, "footprinting_motifs_enriched_per_cluster_comparison")
    if (!file.exists(file.path(out_per_cluster))) {dir.create(out_per_cluster)}
    for(i in 1:length(motifsPerCluster)) {
      cat(paste0("Print cluster motif ", motifsPerCluster[i], " (", i, " of ", length(motifsPerCluster), ")\n"))
      pCl <- PlotFootprint(sobj, assay=assay, features = motifsPerCluster[i], group.by = groupPlotBy, split.by = splitPlotBy)
      ggsave(plot=pCl, filename=file.path(out_per_cluster, paste0("motif_footprint_", motifsPerCluster[i], ".pdf"))) 
    }
  }
  
  if(!is.null(motifsPerCelltype)) {
    out_per_celltype <- file.path(out, "footprinting_motifs_enriched_per_celltype_comparison")
    if (!file.exists(file.path(out_per_celltype))) {dir.create(out_per_celltype)}
    for(i in 1:length(motifsPerCelltype)) {
      cat(paste0("Print celltype motif ", motifsPerCelltype[i], " (", i, " of ", length(motifsPerCelltype), ")\n"))
      pCT <- PlotFootprint(sobj, assay=assay, features = motifsPerCelltype[i], group.by = groupPlotBy, split.by = splitPlotBy)
      ggsave(plot=pCT, filename=file.path(out_per_celltype, paste0("motif_footprint_", motifsPerCelltype[i], ".pdf"))) 
    }
  }
  
  if(!is.null(motifsByNameLoaded)) {
    out_per_file <- file.path(out, "footprinting_motifs_loaded_by_file")
    if (!file.exists(file.path(out_per_file))) {dir.create(out_per_file)}
    for(i in 1:length(motifsByNameLoaded)) {
      cat(paste0("Print loaded motif ", motifsByNameLoaded[i], " (", i, " of ", length(motifsByNameLoaded), ")\n"))
      pNL <- PlotFootprint(sobj, assay=assay, features = motifsByNameLoaded[i], group.by = groupPlotBy, split.by = splitPlotBy)
      ggsave(plot=pNL, filename=file.path(out_per_file, paste0("motif_footprint_", motifsByNameLoaded[i], ".pdf"))) 
    }
  }

}



#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/motifFootprinting_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(pCl, pCT, pNL, file=paste0(out,"/motifFootprinting.RData"))

