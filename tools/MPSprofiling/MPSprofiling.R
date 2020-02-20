#####################################
##
## What: MPSprofiling.R
## Who : Frank Rühle
## When: 29-01-2020
##
## Processing barcode count data and calculation of protein stability indices (PSIs)
##
## Args:
## -----
## targets=targets.txt      # file describing the targets.
## prefix=RE                # prefix to remove from the sample name
## suffix=RE                # suffix to remove from the sample name (usually _readcounts.tsv)
## inputdir=.               # input directory where the files .tsv files are located
## out=DE.DESeq2            # prefix filename for output
##
######################################
options(stringsAsFactors=FALSE)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(scater)
library(reshape)
library(kableExtra)
library(plyr)
library(dplyr)
library(Biostrings)
library(mixtools)

library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(openxlsx)
library(rtracklayer)
library(pheatmap)

##
## get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {

    if(length(i <- grep(string,args,fixed=T)) == 1)
        return(do.call(convert,list(gsub(string,"",args[i]))))
    
    if(!is.null(default)) default else do.call(convert,list(NA))
}

args <- commandArgs(T)
ftargets     <- parseArgs(args,"targets=","targets.txt")     # file describing the targets
pre          <- parseArgs(args,"prefix=","")    # prefix to remove from the sample name
suf          <- parseArgs(args,"suffix=","")    # suffix to remove from the sample name
inputdir     <- parseArgs(args,"inputdir=","./")     # input files directory
out          <- parseArgs(args,"out=","MPSprofiling") # output directory

runstr <- "Rscript MPSprofiling.R [targets=targets.txt] [prefix=RE] [suffix=RE] [inputdir=.] [out=MPSprofiling]"
if(!file.exists(ftargets))   stop(paste("File",ftargets,"does NOT exist. Run with:\n",runstr))
if(!file.exists(inputdir))   stop(paste("Dir",inputdir,"does NOT exist. Run with:\n",runstr))

if(!dir.exists(out)) {dir.create(out, recursive = T)}

name_plotfolder <- file.path(out, "plots")
if(!dir.exists(name_plotfolder)) {dir.create(name_plotfolder, recursive = T)}


annoFactors <- "group"
qcmetrics <- c("sum", "detected")
NMADS =3 # number of absolute deviations from median. 

#######
# out <- "/fsimb/groups/imb-bioinfocf/projects/khmelinskii/imb_khmelinskii_2019_02_nieto_amplicon_test/results/MPSprofiling"
# inputdir <- "/fsimb/groups/imb-bioinfocf/projects/khmelinskii/imb_khmelinskii_2019_02_nieto_amplicon_test/results/barcode_count"
# ftargets <- "/fsimb/groups/imb-bioinfocf/projects/khmelinskii/imb_khmelinskii_2019_02_nieto_amplicon_test/targets.txt"
# pre <- ""
# suf <- ""
######




## load counts
f <- list.files(paste0(inputdir), pattern="\\.barcode_count\\.tsv$", full.names=TRUE) 
f <- f[!grepl("Undetermined", f)] # don't use count file with unmapped reads if present

counts <- mclapply(f, read.delim, header=T, row.names=1, comment.char = "#")
counts <- lapply(counts, "[", -1)

# merge count data objects into one dataset with union of barcodes
custom.merge <- function(x,y) { # custom merge function for Reduce
  z <- merge(x, y, all=T, by="row.names")
  rownames(z) <- z$Row.names
  colnames(z) <- 0:(ncol(z)-1)
  return(z[,-1])
}
counts <- Reduce(custom.merge, counts)

# rename colnames by samples
colnames(counts) <- basename(f)
if(pre != ""){gsub(pre, "", colnames(counts))}
if(suf != ""){gsub(suf, "", colnames(counts))}
colnames(counts) <- gsub(lcSuffix(colnames(counts)),"", colnames(counts))
colnames(counts) <- gsub(lcPrefix(colnames(counts)),"", colnames(counts))

# remove barcodes containg Ns
counts <- counts[!grepl("N", rownames(counts)), ]

# replace NA by zero
counts[is.na(counts)] <- 0

# load targets
targets <- read.delim(ftargets, header=T, comment.char="#", sep=",")
#if(!all(c("group", "file", "sample") %in% colnames(targets))) stop("targets file must have at least 3 columns and fit the format expected in DESeqDatcondition")
# what columns obligatory? group? experiment?
stopifnot(all(colnames(counts) %in% targets$unique_sample_id))
targets <- targets[match(colnames(counts), targets$unique_sample_id, nomatch=0), ]  # subset targets to match only samples with counts files


#### create SingleCellExperiment object
# to make use of sce visualisiation tools
sce_all <- SingleCellExperiment(assays=list(counts=as.matrix(counts)), colData = targets)
# one pool of cells (devided by fluorescence intensity) represents one cell in scRNAseq
# one experiment (8 pools and replicates) represents a group of cells in scRNAseq
# one cell barcode represents a gene in scRNAseq


######## quality control 
# prepare qc metrics
#sce_all <- calculateQCMetrics(sce_all, use_spikes=F) # deprecated.
sce_all <- scater::addPerCellQC(sce_all) # per column, creates sce_all$sum

# 1st filter step for extreme outlier
totalcount_outlier <- sce_all$total < 0.01* mean(sce_all$total)
cat("\nremove samples from dataset with < 1% of mean total counts:\n", colnames(sce_all)[totalcount_outlier], "\n")
sce_all <- sce_all[,!totalcount_outlier]

sce_all <- scater::addPerFeatureQC(sce_all) # per row

qc.frame <- colData(sce_all)
qc.frame$sum <- sce_all$sum/1e3 # library size in thousands
qc.frame <- as.data.frame(qc.frame)

# violinplots
qc.plots.violin <- lapply(qcmetrics, function(to.plot){ 
  if(to.plot=="sum"){
    ylabel <- "sum in thousands"
  }else{
    ylabel <- to.plot
  }
  
  plot.list <- lapply(annoFactors, function(separation){
    p <- ggplot(qc.frame, aes_string(separation,to.plot,color=separation))+
      geom_violin() +
       geom_quasirandom() +
      scale_fill_hue(l=40, c=40) +
      ylab(ylabel) +
      xlab(separation) +
      theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
    plot(p)
    return(p)
  })
  return(plot.list)
})

# QC filtering (not applied)
qc.drop <-  data.frame(row.names = colnames(sce_all))
qc.drop$libsize <- scater::isOutlier(sce_all$sum, nmads=NMADS, type="lower", log=TRUE)
qc.drop$feature_count <- scater::isOutlier(sce_all$detected, nmads=NMADS, type="lower", log=TRUE) 
qc.drop$pass <- !apply(qc.drop, 1, any)

qcfailed <- data.frame(criterion=c(paste("total counts <", NMADS, "MAD"),
                                   paste("total features <", NMADS, "MAD"),
                                   paste("remaining")), 
                       fractions=c(sum(qc.drop$libsize), sum(qc.drop$feature), sum(qc.drop$pass)))
kable(qcfailed, format="markdown", caption="QC filtering") %>% kable_styling()

# PCAplot mit QC metrics
# annoFactors <- c("group", "experiment")
# plotPCAfromQCmetrics(sce_all, metrics=qcmetrics, anno=annoFactors, qc.drop=qc.drop) +  
#   ggtitle(paste("PCA plot of QC metrics"))

# filter out cells failing QC (optional); avoid executing multiple times
# cat("\nremove indicated samples:\n",  rownames(qc.drop[!qc.drop$pass,]))
# if(length(qc.drop$pass) == ncol(sce_all)) { sce_all <- sce_all[, qc.drop$pass]} 
############### end qc


## pca plot of raw counts
set.seed(100)
sce_all <-scater::runPCA(sce_all, name = "PCA", ncomponents=2, ntop = nrow(sce_all), exprs_values="counts")
scater::plotReducedDim(sce_all, dimred = "PCA", colour_by=annoFactors[1],  point_size=2) +  
  ggtitle(paste("PCA plot of raw count data"))




####### new implementation

# initalise result lists for all experiments
raw_counts_all <- list()
bynuc_all <- list()
byaa_all <- list()


#### start analsis per experiment
for (e in unique(sce_all$experiment)) {
  
  sce <- sce_all[, sce_all$experiment==e]
  


# column "sub_experiment" used e.g. for strain names? 
  use_sub_experiment <- !( !"sub_experiment" %in% colnames(colData(sce)) || any(is.na(sce$sub_experiment)) || any(sce$sub_experiment=="" || all(sce$sub_experiment == sce$experiment)))
  if(use_sub_experiment) {helper_cat <- paste0(sce$experiment, sce$sub_experiment)} else {helper_cat <- sce$experiment}

# prepare transformed bin_index f per (sub)experiment
  #sce$totalbins <- sapply(helper_cat,  function(x) {sum(helper_cat==x)}) # total number of bins per (sub)experiment
  for (i in unique(helper_cat)) { # rank bins per (sub)experiment (in case bin index is not starting with 1)
    colData(sce)$bin_rank[helper_cat==i] <- rank(colData(sce)$bin[helper_cat==i])
    colData(sce)$totalbins[helper_cat==i] <- nlevels(factor(colData(sce)$bin[helper_cat==i]))
  }
  sce$f <- (sce$bin_rank -1) / (sce$totalbins -1) # transformed bin index [0,1]

# make MultiAssayExperiment for switching to long format
  raw_counts <- MultiAssayExperiment(ExperimentList(counts=assay(sce, "counts")), colData=DataFrame(colData(sce)))
  raw_counts <- longFormat(raw_counts, colDataCols=colnames(colData(raw_counts)))
  colnames(raw_counts)[colnames(raw_counts) %in% c("rowname", "value")] <- c("sequence", "counts")
  
  raw_counts <- as.data.frame(raw_counts)
  
  # in case there are sequences with zero counts, these entries are removed
  # (these entries may have emerged when the samples were merged together in the first place)
  logindex_nocounts <- raw_counts$counts==0
  if(sum(logindex_nocounts)>0) {
    cat("\n", sum(logindex_nocounts), "entries removed from raw count matrix due to zero counts:\n")
    print(raw_counts[logindex_nocounts, c("sample_name", "bin", "sequence")])
    cat("\n")
    raw_counts <- raw_counts[!logindex_nocounts,]
  }
  
  
  # helper column for background subtraction and plotting if column sub_experiment available
  if (use_sub_experiment) {raw_counts$helper_exp_col <- raw_counts$sub_experiment} else {
    raw_counts$helper_exp_col <- raw_counts$experiment
  }
  
### background subtraction
# and run downstream calculations separated for with and without background subtraction
  
   # To try to remove the noise, we will try to fit a Gaussian mixture model to each sample and subtract the mean of the background distribution. In addition,
  # we will threshold the raw counts to remove the background distribution. There are several cases to deal with: The case of well-separable distributions,
  # like wild type fraction E, the case of non-separable distributions like ubr1-mak3 fraction H, and the case where the background distribution completely
  # dominates, like ubr1-ufd4 fraction B. The threshold depends on the shape and position of the distributions. If the 97.5th percentile of the background
  # distribution is larger than the 2.5th percentile of the foreground distribution and the mean of the background distribution is lower than the mean of
  # the foreground distribution, we set the threshold at the 2.5th percentile of the foreground distribution to include all foreground reads, since this case
  # corresponds to non-separable distribution. Otherwise we threshold at the 97.5th percentile of the background distribution, which applies to well-separable
  # cases and cases with completely dominant background distribution.
  
  set.seed(100)
  raw_bgsubt <- raw_counts %>% 
    group_by(helper_exp_col, bin) %>%
    do({
      ret <- .
      mix <- normalmixEM(log10(ret$counts))
      ubounds <- 10^qnorm(0.975, mix$mu, mix$sigma)
      lbounds <- 10^qnorm(0.025, mix$mu, mix$sigma)
      bg_comp <- which.min(ubounds)
      if (ubounds[bg_comp] > lbounds[-bg_comp] && mix$mu[bg_comp] < mix$mu[-bg_comp]) {
        thresh <- lbounds[-bg_comp]
      } #else if (ubounds[bg_comp] > lbounds[-bg_comp] && mix$mu[bg_comp] < mix$mu[-bg_comp])
      else {
        thresh <- ubounds[bg_comp]
      }
      #ret$counts_bgsubt <- if_else(ret$counts > thresh, ret$counts, 0)
      ret$counts_bgsubt <- pmax(0, ret$counts - 10^mix$mu[bg_comp])
      ret$counts <- ret$counts_bgsubt
      ret$counts_bgsubt <- NULL
      ret
    })
   raw_bgsubt <- as.data.frame(raw_bgsubt)
  
   # in case there are sequences with zero counts after backgound subtraction, these entries are removed
   # (these entries may have emerged when background is subtracted from samm count numbers)
   logindex_nocounts_bgsubt <- raw_bgsubt$counts==0
   if(sum(logindex_nocounts_bgsubt)>0) {
     cat("\n", sum(logindex_nocounts_bgsubt), "entries removed from count matrix due to zero counts after background subtraction:\n")
     print(raw_bgsubt[logindex_nocounts_bgsubt, c("sample_name", "bin", "sequence")][1:min(20,sum(logindex_nocounts_bgsubt)),]) # print max 20 entries
     cat("\n")
     raw_bgsubt <- raw_bgsubt[!logindex_nocounts_bgsubt,]
   }
   
   
   
   
 ### run the following code separately for raw_counts and raw_bgsubt
  # raw_counts and raw_bgsubt have identical format and differ just in column "counts"
  
  for(bg in c("raw", "bgsubt")) {
  
    if(bg=="bgsubt") {raw_counts <- raw_bgsubt} 
   

    plot_histo <- raw_counts %>%
      filter(counts > 0) %>%
      ggplot(aes(counts)) +
      geom_histogram(bins=100) +
      facet_grid(helper_exp_col ~ bin) +
      scale_x_log10() +
      theme(strip.text.y = element_text(angle=0))
    
    ggsave(filename=file.path(name_plotfolder, paste0("histogram_",e , "_", bg, ".png")),
           plot=plot_histo, 
           width = 200, height = 150, 
           units = c("mm"),  dpi = 600, device="png")
   
      ## cell distribution per variant normalised by bin countsum and cell fraction:
      # divide raw counts by sum counts for all sequences from one strain (barcode_fw) and  bin (barcode_rev)
      # and multiply by cell fraction
      # This normalized count (aka Cp) is used by CombinatorialProfiler
        raw_counts$normalized_counts <- (raw_counts$counts / raw_counts$sum) * raw_counts$fraction
        
      # transformed bin index times normalised counts
       raw_counts$fxCp <- raw_counts$f * raw_counts$normalized_counts 
      
      # translate nucleotide to codon with Biostrings # confimed with CombinatorialProfiler
      # no.init.codon=T because otherwise the alternative initiation codons CTG and TTG are 
      # translated to M instead to L when they are located at the first position of the hexamer.
      # alternative from GENETIC_CODE_TABLE: getGeneticCode("Alternative Yeast Nuclear", full.search=T)
      # The code "Alternative Yeast Nuclear" differs from GENETIC_CODE only in alt_init_codons
        raw_counts$translation <- as.character(translate(DNAStringSet(as.character(raw_counts$sequence)),
                                              genetic.code=GENETIC_CODE, no.init.codon=T))
      
      # produce output table
        raw_counts <- data.frame(raw_counts[,c(colnames(targets), "bin_rank", "f", "fxCp", "sequence", "translation", "counts", "normalized_counts")])
      
      
      ## calculate PSI (aka dsi) per hexamer (bynuc)
          listcats4index <- list(sequence=raw_counts$sequence) # determine either 1 or 2 index categories for tapply (sequence, sub_experiment)
          if(use_sub_experiment) {
            listcats4index[["sub_experiment"]] <- raw_counts$sub_experiment # use sub_experiment only if given
            }
          
          bynuc <- tapply(raw_counts$fxCp, INDEX = listcats4index, sum) # calculate PSI numerator "fxCp_sum" 
          bynuc <- melt(bynuc) # for the case we have two INDEX variables
          colnames(bynuc) <- c(names(listcats4index), "fxCp_sum")
          bynuc$Cp_sum <- melt(tapply(raw_counts$normalized_counts, listcats4index, sum))[,"value"] # calculate PSI denominator "Cp_sum"
          bynuc$PSI <- bynuc$fxCp_sum / bynuc$Cp_sum # PSI
          
          # statistics
          statfunctions <- list("normalized_counts_min" = min, "normalized_counts_max" = max, "normalized_counts_mean" = mean, 
                             "normalized_counts_median" = median, "normalized_counts_std" = sd, "normalized_counts_sum" = sum)
          statfunctions_raw <- statfunctions
          names(statfunctions_raw) <- gsub("normalized_", "", names(statfunctions_raw))
          
          normcountstats <- sapply(statfunctions, function(x){
             melt(tapply(raw_counts$normalized_counts, listcats4index, x, na.rm = TRUE))[,"value"]} 
             )
          rawcountstats <- sapply(statfunctions_raw, function(x){
            melt(tapply(raw_counts$counts, listcats4index, x, na.rm = TRUE))[,"value"]} 
            )
          
          bynuc <- data.frame(experiment=e, 
                              translation=as.character(translate(DNAStringSet(as.character(bynuc$sequence)), genetic.code=GENETIC_CODE, no.init.codon=T)),
                              bynuc[,!colnames(bynuc) %in% c("fxCp_sum", "Cp_sum")], 
                              normcountstats,
                              rawcountstats,
                              nfractions= melt(tapply(factor(raw_counts$bin), listcats4index, nlevels))[,"value"])

      ## calculate PSI per di-residue (byaa)
          listcats4index_byaa <- list(translation=bynuc$translation)
          if(use_sub_experiment) {
            listcats4index_byaa[["sub_experiment"]] <- bynuc$sub_experiment # use sub_experiment only if given
          }
          
          byaa <- tapply(bynuc$PSI, listcats4index_byaa, median, na.rm=T)  
          byaa <- data.frame(rows=rownames(byaa), byaa, check.names=F, row.names = NULL)
          # transform to dataframe before melting because otherwise the di-residue "NA" is interpreted as NA
          byaa <- melt(byaa) 
          if(length(listcats4index_byaa)==1) {byaa <- byaa[,names(byaa) != "variable"]} # remove column variable=byaa
          colnames(byaa) <- c(names(listcats4index_byaa), "median_PSI")
          
          
        ## calculate pooled PSI per di-residue
          listcats4index_byaa_pooled <- list(translation=raw_counts$translation)
          if(use_sub_experiment) {
            listcats4index_byaa_pooled[["sub_experiment"]] <- raw_counts$sub_experiment # use sub_experiment only if given
          }
          
          byaa$pooled_fxCp_sum <- melt(tapply(raw_counts$fxCp, listcats4index_byaa_pooled, sum))[,"value"]
          #byaa$transl2 <- melt(tapply(raw_counts$fxCp, listcats4index_byaa_pooled, sum))[,"translation"] # just to check correct row order
          byaa$pooled_Cp_sum <- melt(tapply(raw_counts$normalized_counts, listcats4index_byaa_pooled, sum))[,"value"]
          byaa$pooled_PSI <- byaa$pooled_fxCp_sum / byaa$pooled_Cp_sum
          
          byaa_normcountstats <- sapply(statfunctions, function(x){
            melt(tapply(raw_counts$normalized_counts, listcats4index_byaa_pooled, x, na.rm = TRUE))[,"value"]} 
          )
          byaa_rawcountstats <- sapply(statfunctions_raw, function(x){
            melt(tapply(raw_counts$counts, listcats4index_byaa_pooled, x, na.rm = TRUE))[,"value"]} 
          )
          
          # create output table
          byaa <- data.frame(experiment=e, 
                             byaa[,!colnames(byaa) %in% c("pooled_fxCp_sum", "pooled_Cp_sum")], 
                             byaa_normcountstats,
                             byaa_rawcountstats,
                             nfractions= melt(tapply(factor(raw_counts$bin), listcats4index_byaa_pooled, nlevels))[,"value"],
                             nsequences= melt(tapply(bynuc$PSI, listcats4index_byaa, length))[,"value"]
                             )
          
          # store results in list
          raw_counts_all[[bg]][[e]] <- raw_counts
          bynuc_all[[bg]][[e]] <- bynuc
          byaa_all[[bg]][[e]] <- byaa
          
          write.table(raw_counts, file= file.path(out, paste0(e, "_", bg, "_counts.txt")), sep="\t", row.names = F, quote = F)
          write.table(bynuc, file= file.path(out, paste0(e, "_", bg, "_PSI_bynuc.txt")), sep="\t", row.names = F, quote = F)
          write.table(byaa, file= file.path(out, paste0(e, "_", bg, "_PSI_byaa.txt")), sep="\t", row.names = F, quote = F)
          
  } # end bg loop      
} # end e loop


for(bg in c("raw", "bgsubt")) {
  
  if(length(raw_counts_all[[bg]]) >1) {
  # rbind output tables
  raw_counts_all[[bg]][["all_exp"]] <- rbind.fill(raw_counts_all[[bg]])
  bynuc_all[[bg]][["all_exp"]] <- rbind.fill(bynuc_all[[bg]])
  byaa_all[[bg]][["all_exp"]] <- rbind.fill(byaa_all[[bg]])

  write.table(raw_counts_all[[bg]][["all_exp"]], file= file.path(out, paste0(bg, "_counts_all.txt")), sep="\t", row.names = F, quote = F)
  write.table(bynuc_all[[bg]][["all_exp"]], file= file.path(out, paste0(bg, "_PSI_bynuc_all.txt")), sep="\t", row.names = F, quote = F)
  write.table(byaa_all[[bg]][["all_exp"]], file= file.path(out, paste0(bg, "_PSI_byaa_all.txt")), sep="\t", row.names = F, quote = F)

  }
}
############################################################################



#save the sessionInformation
writeLines(capture.output(sessionInfo()),paste(out, "/MPSprofiling_session_info.txt", sep=""))
save(sce_all, raw_counts_all, bynuc_all, byaa_all, file=paste0(out,"/MPSprofiling.RData"))











