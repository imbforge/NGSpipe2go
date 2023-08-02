#####################################
##
## What: peaks2genes.R
## Who : Frank Rühle
## When: 05.06.2023
##
## Script to compute the correlation between gene expression and accessibility at nearby peaks.
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
projectdir     <- parseArgs(args,"project=") 
resultsdir     <- parseArgs(args,"res=")   
out            <- parseArgs(args,"outdir=") # output folder
db             <- parseArgs(args,"db=")     
genes2use      <- parseArgs(args,"genes2use=", convert="run_custom_code")
genes2plot     <- parseArgs(args,"genes2plot=", convert="as.character")
groupCellsInPlot <- parseArgs(args,"groupCellsInPlot=", convert="as.character")
plotUpstream   <- parseArgs(args,"plotUpstream=", default=0, convert="as.numeric")
plotDownstream <- parseArgs(args,"plotDownstream=", default=0, convert="as.numeric")
if(length(genes2use)==1 && is.na(genes2use)) {genes2use <- NULL}
if(length(genes2plot)==1 && is.na(genes2plot)) {genes2plot <- NULL}

runstr <- "Rscript peaks2genes.R [projectdir=projectdir]"

# load R environment
renv::use(lockfile=file.path(projectdir, "NGSpipe2go/tools/sc_peaks2genes/renv.lock"))
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
print(paste("genes2plot:", paste(genes2plot, collapse=" ")))
print(paste("genes2use:", paste(genes2use, collapse=" ")))
print(paste("groupCellsInPlot:", groupCellsInPlot))
print(paste("plotUpstream:", plotUpstream))
print(paste("plotDownstream:", plotDownstream))

# load relevant BSgenome package (needed by Signac for motif analysis)
switch(db,
       hg38={ failed_BSgenome <- library("BSgenome.Hsapiens.UCSC.hg38") 
       BSgenome <- BSgenome.Hsapiens.UCSC.hg38
       },
       mm10={ failed_BSgenome <- library("BSgenome.Mmusculus.UCSC.mm10") 
       BSgenome <- BSgenome.Mmusculus.UCSC.mm10
       },
       stop(c("Don't find genome:", db))   
)


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "ATAC"
load(file.path(resultsdir, "/gtf.RData"))


# first compute the GC content for each peak
sobj <- RegionStats(sobj, assay = "ATAC", genome = BSgenome)

sobj <- LinkPeaks(
  object = sobj,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = genes2use
)

# access with Links(sobj)
write.table(as.data.frame(Links(sobj)), file=file.path(out, "link_peaks2genes.txt"), sep="\t", quote = F, row.names = F)

sobj@assays$ATAC@annotation$tx_id <- sobj@assays$ATAC@annotation$transcript_id


# We can visualize these links using the CoveragePlot() function, or alternatively we could use the CoverageBrowser() function in an interactive analysis.

# Plot frequency of Tn5 insertion events (accessibility) for different groups of cells within given regions of the genome. 
# Tracks are normalized using a per-group scaling factor computed as the number of cells in the group multiplied by the mean 
# sequencing depth for that group of cells. This accounts for differences in number of cells and potential differences 
# in sequencing depth between groups. 

covplot <- NULL
if(!is.null(genes2plot)) {
  covplot <- CoveragePlot(
    object = sobj,
    assay = "ATAC",
    region = gtf[genes2plot],
    features = genes2plot,
    expression.assay = "SCT",
    split.assays = T,
    links = T,
    annotation = "gene", 
    group.by = groupCellsInPlot,
    #region.highlight =StringToGRanges(gtf["ACAP3"]) ,
    # expression.assay = "RNA",
    # expression.slot     = "data",
    #idents = idents.plot,
    extend.upstream = plotUpstream,
    extend.downstream = plotDownstream
  )

  ggsave(plot=covplot, filename=file.path(out, "peaks2genes_CoveragePlot.pdf"))
  ggsave(plot=covplot, filename=file.path(out, "peaks2genes_CoveragePlot.png"))
}

# Note that Cell Ranger ARC already does such peak-to-gene linking based on correlating peak accessibility to gene expression 
# (read in above into the `feature_linkage` object). This is the `Signac` version of the algorithm.



#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/peaks2genes_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(covplot, file=file.path(out,"peaks2genes.RData"))

