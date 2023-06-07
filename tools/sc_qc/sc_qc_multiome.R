#####################################
##
## What: sc_qc_multiome.R
## Who : Frank Rühle
## When: 23.05.2023
##
## Script to perform quality control metrics of 10X multiome assay data.
##
## Args:
## -----
## targets=targets.txt      # file describing the targets 
##
##
######################################

renv::use(lockfile='NGSpipe2go/tools/sc_qc/renv.lock')
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
db            <- parseArgs(args,"db=")     

runstr <- "Rscript sc_qc_multiome.R [projectdir=projectdir]"

print(paste("resultsdir:", resultsdir))
print(paste("out:", out))
print(paste("projectdir:", projectdir))
print(paste("db:", db))


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "ATAC"


# compute nucleosome signal score per cell
sobj <- NucleosomeSignal(sobj)
# compute TSS enrichment score per cell
sobj <- TSSEnrichment(sobj, fast = FALSE)

# count reads in peaks
sobj <- FRiP(sobj, assay = 'ATAC', total.fragments = 'fragments')

# count reads in blacklisted regions
# select relevant blacklist
blacklist <- switch(db,
       hg38={ blacklist_hg38_unified},
       hg19={ blacklist_hg19},
       mm10={ blacklist_mm10},
       dm6 ={ blacklist_dm6},
       dm3 ={ blacklist_dm3},
       ce11={ blacklist_ce11},
       ce10={ blacklist_ce10},
       warning(c("Don't find blacklist for genome:", db))   
      )

if(class(blacklist) == "GRanges") {
  sobj$blacklist_fraction <- FractionCountsInRegion(sobj, assay = 'ATAC', regions = blacklist)
}
# add blacklist ratio and fraction of reads in peaks
#sobj$pct_reads_in_peaks <- sobj$peak_region_fragments / sobj$passed_filters * 100
#sobj$blacklist_ratio <- sobj$blacklist_region_fragments / sobj$peak_region_fragments



# Nucleosome signal
fraghist <- FragmentHistogram(object = sobj, group.by = 'sample', region = 'chr1-1-10000000') +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("Fragment size histograms") +
  NoLegend()
ggsave(plot=fraghist, filename = file.path(out, "atac_qc_fragmentsizes_hist.pdf"))
ggsave(plot=fraghist, filename = file.path(out, "atac_qc_fragmentsizes_hist.png"))


# TSS enrichment
tss <- TSSPlot(sobj, group.by = 'sample') +
  theme(plot.title=element_text(face="bold", hjust=0.5)) +
  ggtitle("Fragment enrichment at TSS") +
  NoLegend()
ggsave(plot=tss, filename = file.path(out, "atac_qc_tss_enrich_profiles.pdf"))
ggsave(plot=tss, filename = file.path(out, "atac_qc_tss_enrich_profiles.png"))


# qc measures
features2plot <- c("nCount_RNA", "nCount_ATAC", "FRiP", "blacklist_fraction", "TSS.enrichment", "nucleosome_signal") 
features2plot <- features2plot[features2plot %in% colnames(sobj[[]])]
vln <- VlnPlot(
  object = sobj,
  features = features2plot,
  ncol = 3,
  #pt.size = 0,
  group.by="sample"
)
ggsave(plot=vln, filename = file.path(out, "qc_violin_plot.pdf"), height = 9, width = 9)
ggsave(plot=vln, filename = file.path(out, "qc_violin_plot.png"), height = 9, width = 9)




#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/sc_qc_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(fraghist, tss, vln, file=paste0(out,"/sc_qc.RData"))

