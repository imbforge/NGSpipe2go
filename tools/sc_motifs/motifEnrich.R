#####################################
##
## What: motifEnrich.R
## Who : Frank Rühle
## When: 06.06.2023
##
## Script to identify enriched motifs in sets of differentially accessible peaks between all pairs of groups for each cluster and celltype.
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
diffPeaks_dir <- parseArgs(args,"diffPeaks_dir=")     
pval_thresh   <- parseArgs(args,"pval_thresh=", convert="as.numeric")
min_peaks     <- parseArgs(args,"min_peaks=", convert="as.numeric")

runstr <- "Rscript motifEnrich.R [projectdir=projectdir]"

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
print(paste("diffPeaks_dir:", diffPeaks_dir))
print(paste("pval_thresh:", pval_thresh))
print(paste("min_peaks:", min_peaks))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "ATAC"
load(file.path(resultsdir, "/gtf.RData"))
load(file.path(diffPeaks_dir, "/diffPeaks.RData"))


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

## Identify enriched motifs in differentially accessible peaks per cluster
# FindMotifs function: Find motifs over-represented in a given set of genomic features. 
# Computes the number of features containing the motif (observed) and compares this to the 
# total number of features containing the motif (background) using the hypergeometric test. 

cat(paste("\nLook for enriched motifs in differentially accessible peaks filtered for p_val_adj <", pval_thresh, "(motif enrichment is skipped if less than", min_peaks, "peaks)\t"))

da_enriched_motifs <- list()

for (gp in names(da_groups_atac)) {
  top_da_peaks <- rownames(da_groups_atac[[gp]][da_groups_atac[[gp]]$p_val_adj < pval_thresh, ]) # filter regions by significance
  cat("Finding motifs for:", gp, "(contains", length(top_da_peaks), "differentially accessible peaks after filtering)\n")
  if(length(top_da_peaks)>= min_peaks) {
    da_enriched_motifs[[gp]] <- FindMotifs(object = sobj, features = top_da_peaks)
    mplot <- MotifPlot(object = sobj, motifs = head(rownames(da_enriched_motifs[[gp]])))
    ggsave(plot=mplot, filename=file.path(out, paste0("Motif_position_weight_matrices_top_enriched_cluster_", gp, ".pdf")), width = 7, height = 7)
  }
}

# Save results
if(length(da_enriched_motifs)>0) {
  names(da_enriched_motifs) <- substr(names(da_enriched_motifs), 1, apply(data.frame(nchar(names(da_enriched_motifs)),31), 1, min))
  openxlsx::write.xlsx(da_enriched_motifs, file = file.path(out, "motif_enriched_clusters.xlsx"))
}
gc()


## Identify enriched motifs in differentially accessible peaks per celltype
# FindMotifs function: Find motifs over-represented in a given set of genomic features. 
# Computes the number of features containing the motif (observed) and compares this to the 
# total number of features containing the motif (background) using the hypergeometric test. 

cat(paste("\nLook for enriched motifs in differentially accessible peaks filtered for p_val_adj <", pval_thresh, "(motif enrichment is skipped if less than", min_peaks, "peaks)\t"))

da_enriched_motifs_ct <- list()

for (gp in names(da_groups_atac_ct)) {
  top_da_peaks_ct <- rownames(da_groups_atac_ct[[gp]][da_groups_atac_ct[[gp]]$p_val_adj < pval_thresh, ]) # filter regions by significance
  cat("Finding motifs for:", gp, "(contains", length(top_da_peaks_ct), "differentially accessible peaks after filtering)\n")
  if(length(top_da_peaks_ct)>= min_peaks) {
    da_enriched_motifs_ct[[gp]] <- FindMotifs(object = sobj, features = top_da_peaks_ct)
    mplot <- MotifPlot(object = sobj, motifs = head(rownames(da_enriched_motifs_ct[[gp]])))
    ggsave(plot=mplot, filename=file.path(out, paste0("Motif_position_weight_matrices_top_enriched_celltype_", gp, ".pdf")), width = 7, height = 7)
  }
}

# Save results
if(length(da_enriched_motifs_ct)>0) {
  names(da_enriched_motifs_ct) <- substr(names(da_enriched_motifs_ct), 1, apply(data.frame(nchar(names(da_enriched_motifs_ct)),31), 1, min))
  openxlsx::write.xlsx(da_enriched_motifs_ct, file = file.path(out, "motif_enriched_celltypes.xlsx"))
}
gc()




#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/motifEnrich_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(da_enriched_motifs, da_enriched_motifs_ct, file=paste0(out,"/motifEnrich.RData"))

