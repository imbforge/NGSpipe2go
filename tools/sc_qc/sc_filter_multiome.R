#####################################
##
## What: sc_filter_multiome.R
## Who : Frank Rühle
## When: 01.06.2023
##
## Script to apply quality control thresholds to 10X multiome assay data.
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
projectdir         <- parseArgs(args,"project=") 
resultsdir         <- parseArgs(args,"res=")   
out                <- parseArgs(args,"outdir=") # output folder
nCount_ATAC_min    <- parseArgs(args,"nCount_ATAC_min=", convert="as.numeric")     
nCount_ATAC_max    <- parseArgs(args,"nCount_ATAC_max=", convert="as.numeric")     
nCount_RNA_min     <- parseArgs(args,"nCount_RNA_min=", convert="as.numeric")     
nCount_RNA_max     <- parseArgs(args,"nCount_RNA_max=", convert="as.numeric")     
FRiPmin            <- parseArgs(args,"FRiPmin=", convert="as.numeric")     
FRiBLmax           <- parseArgs(args,"FRiBLmax=", convert="as.numeric")     
nucleosome_sig_max <- parseArgs(args,"nucleosome_sig_max=", convert="as.numeric")     
TSS_enrich_min     <- parseArgs(args,"TSS_enrich_min=", convert="as.numeric")     

runstr <- "Rscript sc_filter_multiome.R [projectdir=projectdir] "



print(paste("projectdir:", projectdir))
print(paste("resultsdir:", resultsdir))
print(paste("out:", out))
print(paste("nCount_ATAC_min:", nCount_ATAC_min))
print(paste("nCount_ATAC_max:", nCount_ATAC_max))
print(paste("nCount_RNA_min:", nCount_RNA_min))
print(paste("nCount_RNA_max:", nCount_RNA_max))
print(paste("FRiPmin:", FRiPmin))
print(paste("FRiBLmax:", FRiBLmax))
print(paste("nucleosome_sig_max:", nucleosome_sig_max))
print(paste("TSS_enrich_min:", TSS_enrich_min))


# overview table with thresholds:
qcfilt <- cbind(
  criterion=c(nCount_ATAC_min = paste("nCount ATAC >", nCount_ATAC_min),
              nCount_ATAC_max = paste("nCount ATAC <", nCount_ATAC_max),
              nCount_RNA_min  = paste("nCount RNA >", nCount_RNA_min),
              nCount_RNA_max  = paste("nCount RNA <", nCount_RNA_max),
              FRiPmin         = paste("Fraction of reads in peaks >", FRiPmin),
              FRiBLmax        = paste("Fraction of reads in blacklisted regions <", FRiBLmax),
              nucleosome_sig_max = paste("nucleosome signal <", nucleosome_sig_max),
              TSS_enrich_min = paste("TSS enrichment >", TSS_enrich_min)
  ))
write.table(qcfilt, file= file.path("/fsimb/groups/imb-bioinfocf/projects/cfb_internal/frank/multiome500/qc/sc_qc/", "qc_filtering.txt"), sep="\t", quote=F, row.names = F)              


# load sobj from previous module
sobj <- readRDS(file = file.path(resultsdir, "sobj.RDS"))
DefaultAssay(sobj) <- "ATAC"



## apply qc filters
sobj <- subset(
  sobj,
  subset = nCount_ATAC > nCount_ATAC_min & 
    nCount_ATAC < nCount_ATAC_max & 
    nCount_RNA > nCount_RNA_min & 
    nCount_RNA < nCount_RNA_max &
    FRiP > FRiPmin & 
    nucleosome_signal < nucleosome_sig_max & 
    TSS.enrichment > TSS_enrich_min
)

if("blacklist_fraction" %in% colnames(sobj[[]])) {
  sobj <- subset(
    sobj,
    subset = blacklist_fraction < FRiBLmax 
  )
} else {
  warning("blacklist_fraction was not found in sobj and therefore the FRiBLmax criterion could not be applied!")
}


# Plot qc metrics after filtering
features2plot <- c("nCount_RNA", "nCount_ATAC", "FRiP", "blacklist_fraction", "TSS.enrichment", "nucleosome_signal") 
features2plot <- features2plot[features2plot %in% colnames(sobj[[]])]
vlnFilt <- VlnPlot(
    object = sobj,
    features = features2plot,
    ncol = 3,
    #pt.size = 0,
    group.by="sample"
  )
ggsave(plot=vlnFilt, filename = file.path(out, "qc_filt_violin_plot.pdf"), height = 9, width = 9)
ggsave(plot=vlnFilt, filename = file.path(out, "qc_filt_violin_plot.png"), height = 9, width = 9)



#############################
# save the sessionInformation and R image
writeLines(capture.output(sessionInfo()),paste0(out, "/sc_filter_session_info.txt"))
saveRDS(sobj, file = file.path(resultsdir, "sobj.RDS"))
save(vlnFilt, file=paste0(out,"/sc_filter.RData"))

