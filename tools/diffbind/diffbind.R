#####################################
##
## What: diffbind.R
## Who : Sergi Sayols
## When: 10-05-2016
##
## Script to perform differential binding analysis between 2 conditions (pairwise)
##
## Args:
## -----
## targets=targets.txt          # file describing the targets.
## contrasts=contrasts.txt      # file describing the contrasts
## cwd=./                       # current working directory where the files .tsv files are located
## bams=paste0(CWD, "/mapped")  # directory with the bam files
## peaks=paste0(CWD, "/results/macs2")  # directory with peak caller output
## out=paste0(CWD, "/results")  # prefix filename for output
## fragsize=200                 # average fragment size
## substractControl=TRUE        # substract input
## fullLibrarySize=TRUE         # use total number of reads in bam for normalization (FALSE=only peaks)
## tagwiseDispersion=TRUE       # calculate dispersion tagwise (use FALSE if no replicates)
## annotate=TRUE                # annotate after DB analysis?
## pe=FALSE                     # paired end experiment?
## tss=c(-3000,3000)            # region around the tss
## txdb=TxDb.Mmusculus.UCSC.mm9.knownGene   # Bioconductor transcript database, for annotation
## annodb=org.Mm.eg.db          # Bioconductor gene annotation database
##
######################################
options(stringsAsFactors=FALSE)
library(DiffBind)
library(openxlsx)
library(dplyr)
library(ggplot2)

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
FTARGETS   <- parseArgs(args,"targets=","targets.txt")     # file describing the targets
FCONTRASTS <- parseArgs(args,"contrasts=","contrasts_diffbind.txt") # file describing the contrasts
CWD        <- parseArgs(args,"cwd=","./")     # current working directory
BAMS       <- parseArgs(args,"bams=",paste0(CWD, "/mapped"))  # directory with the bam files
PEAKS      <- parseArgs(args,"peaks=",paste0(CWD, "/results/macs2"))  # directory with the peak files
OUT        <- parseArgs(args,"out=", paste0(CWD, "/results")) # directory where the output files will go
FRAGSIZE   <- parseArgs(args,"fragsize=", 200, "as.numeric")# fragment size
SUBSTRACTCONTROL  <- parseArgs(args,"substractControl=", TRUE, "as.logical")  # substract input
FULLLIBRARYSIZE   <- parseArgs(args,"fullLibrarySize=", TRUE, "as.logical")   # use total number of reads in bam for normalization (FALSE=only peaks)
TAGWISEDISPERSION <- parseArgs(args,"tagwiseDispersion=", TRUE, "as.logical") # calculate dispersion tagwise (use FALSE if no replicates)
ANNOTATE   <- parseArgs(args,"annotate=", TRUE, "as.logical") # annotate after DB analysis?
PE         <- parseArgs(args,"pe=", FALSE, "as.logical")      # paired end experiment?
TSS        <- parseArgs(args,"tss=", "c(-3000,3000)", "run_custom_code") # region around the tss
TXDB       <- parseArgs(args,"txdb=", "TxDb.Mmusculus.UCSC.mm9.knownGene") # Bioconductor transcript database, for annotation 
ANNODB     <- parseArgs(args,"annodb=", "org.Mm.eg.db") # Bioconductor gene annotation database

runstr <- "Rscript diffbind.R [targets=targets.txt] [contrasts=contrasts.txt] [cwd=./] [bams=./mapped] [peaks=./results/macs2] [out=./results] [fragsize=200] [annotate=TRUE] [pe=TRUE] [tss=c(-3000,3000)] [txdb=TxDb.Mmusculus.UCSC.mm9.knownGene] [annodb=org.Mm.eg.db]"
if(!file.exists(CWD))        stop("Dir",CWD,"does NOT exist. Run with:\n",runstr)
setwd(CWD)
if(!file.exists(FTARGETS))   stop("File",FTARGETS,"does NOT exist. Run with:\n",runstr)
if(!file.exists(FCONTRASTS)) stop("File",FCONTRASTS,"does NOT exist. Run with:\n",runstr)
if(!file.exists(BAMS))       stop("Dir",BAMS,"does NOT exist. Run with:\n",runstr)
if(!file.exists(PEAKS))      stop("Dir",PEAKS,"does NOT exist. Run with:\n",runstr)
if(!is.numeric(FRAGSIZE))    stop("Fragment size not numeric. Run with:\n",runstr)
if(!is.logical(ANNOTATE))    stop("Annotate not logical. Run with:\n",runstr)
if(!is.logical(PE))          stop("Paired end (pe) not logical. Run with:\n",runstr)
if(ANNOTATE & !is.numeric(TSS)) stop("Region around TSS not numeric. Run with:\n",runstr)
if(ANNOTATE & !require(TXDB, character.only=TRUE))   stop("Transcript DB", TXDB, "not installed\n")
if(ANNOTATE & !require(ANNODB, character.only=TRUE)) stop("Annotation DB", ANNODB, "not installed\n")


##
## make DB analysis
##
# load targets and make analysis
conts   <- read.delim(FCONTRASTS, head=F, comment.char="#")
targets <- read.delim(FTARGETS, head=T, colClasses="character", comment.char="#")

# determine file suffixes for targets 
donefiles <- list.files(PEAKS,pattern=".done$")
bam_suffix <- sub("^[^\\.]*\\.", "", gsub("_macs2.done$", "", donefiles[1]))
bam_suffix <- paste0(bam_suffix, ".bam")

peakfiles <- list.files(PEAKS,pattern=".xls")
isBlacklistFilt <- any(grepl("blacklist_filtered", peakfiles))
peak_suffix <- if(isBlacklistFilt) {"_macs2_blacklist_filtered_peaks.xls"} else {"_macs2_peaks.xls"}

# create modified targets file for diffbind
targets <- data.frame(
  SampleID= targets$IPname,
  Condition= targets$group,
  Replicate= targets$Replicate,
  bamReads= paste0(BAMS, "/", targets$IP, ".", bam_suffix),
  ControlID= targets$INPUTname,
  bamControl= paste0(BAMS, "/", targets$INPUT, ".", bam_suffix),
  Peaks= paste0(PEAKS, "/", targets$IPname, ".vs.", targets$INPUTname, peak_suffix),
  PeakCaller= targets$PeakCaller
)


# Construct DBA object
db <- dba(sampleSheet=targets, config=data.frame(fragmentSize=FRAGSIZE, bCorPlot=F, singleEnd=!PE))

# create DBA object containing consensus peaks per group
db2 <- dba.peakset(db, consensus=DBA_CONDITION)

# Heatmap using occupancy (peak caller score) data
png(paste0(OUT, "/heatmap_occupancy.png"), width = 150, height = 150, units = "mm", res=300)
dba.plotHeatmap(db, main="Correlation heatmap\noccupancy data")
dev.off()

# plot Overlap rate
png(paste0(OUT, "/Overlap_rate_plot.png"), width = 150, height = 150, units = "mm", res=300)
  olap.rate <- dba.overlap(db, mode=DBA_OLAP_RATE)
  plot(olap.rate, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets', main="Overlap rate plot")
dev.off()

# identify all overlapping peaks and derives a consensus peakset for the experiment. 
# Then count how many reads overlap each interval for each unique sample.
db <- dba.count(db, bUseSummarizeOverlaps=PE)  # bUseSummarizeOverlaps method slower and memory hungry, mandatory only for PE data

# apply the contrasts
 for (cont in conts[, 1]) {
  # parse formula in cont
  cont.form <- gsub("(.+)=(.+)", "\\2", cont)
  factors   <- unlist(strsplit(cont.form, "\\W"))
  factors   <- factors[factors != ""]
  c1 <- dba.mask(db, DBA_CONDITION, factors[1])
  c2 <- dba.mask(db, DBA_CONDITION, factors[2])
  db <- dba.contrast(db, group1=c1, group2=c2,  name1=factors[1], name2=factors[2], categories=DBA_CONDITION)
  }

  # run the diffbind analysis (DESeq2) for all contrasts
  db <- dba.analyze(db, bSubControl=SUBSTRACTCONTROL, bFullLibrarySize=FULLLIBRARYSIZE, bTagwise=TAGWISEDISPERSION)

# PCA plot
png(paste0(OUT, "/pca_plot_all_samples.png"), width = 150, height = 150, units = "mm", res=300)
  dba.plotPCA(db, DBA_CONDITION, label=DBA_CONDITION)
dev.off()

# Heatmap using affinity (read count) data
png(paste0(OUT, "/heatmap_affinity_scores.png"), width = 150, height = 150, units = "mm", res=300)
  dba.plotHeatmap(db, main="Correlation heatmap\naffinity scores")
dev.off()


# Boxplot of consensus peaks
  lpeaks <- db2$peaks # db2 contains Consensus peaks subsequent to samples 
  dfcolnames <- colnames(lpeaks[[1]]) # Consensus peak dfs don't have colnames, therefore name columns 
  lpeaks <- lapply(1:length(lpeaks), function(x) setNames(lpeaks[[x]], dfcolnames) )
  
  names(lpeaks) <- names(db2$masks$Consensus) # name list elements
  names(lpeaks)[db2$masks$Consensus] <- paste0("Consensus\n", names(lpeaks)[db2$masks$Consensus])
  
  lpeaks <- lapply(lpeaks, function(x) { # calculate peak sizes
    x$width <- abs(x$end - x$start)
    return(x)
  })
  
  for(x in unique(db2$samples$Condition)) { # concatenate peak sizes of samples from the same group
    lpeaks[[paste0("Replicates\n", x)]] <- dplyr::bind_rows(lpeaks[which(db2$samples$Condition==x)])
  }
  lpeaks <- lpeaks[grepl("(Consensus)|(Replicates)", names(lpeaks))] # samples not needed anymore
  
  plotdata <- dplyr::bind_rows(lpeaks, .id="group") # create data.frame for plotting
  
png(paste0(OUT, "/peak_width_boxplot.png"), width = 150, height = 150, units = "mm", res=300)
 try(ggplot(plotdata, aes(x=group, y=width, fill=group)) + geom_boxplot(show.legend = FALSE) + 
      theme(text = element_text(size=15), axis.text.x = element_text(angle=90, vjust=0.5)) + 
      labs(title="Size of consensus peaks") + theme(plot.title = element_text(hjust = 0.5)) )
dev.off()
  

# Venn diagram of all contrasts
png(paste0(OUT, "/venn_plot_all_contrasts.png"), width = 150, height = 150, units = "mm", res=300)
  try(dba.plotVenn(db, main="Binding Site Overlaps Per Contrast", contrast=1:nrow(dba.show(db, bContrast=T))))
dev.off()


# prepare results and plots for each contrast
  result <- lapply(1:nrow(conts), function(cont) {
    cont.name <- substr(gsub("(.+)=\\((.+)\\)", "\\2", conts[cont,1]), 1, 31)
    #cont.name <- gsub("(.+)=(.+)", "\\1", conts[cont,1])
    #cat(cont.name, fill=T)
    
    png(paste0(OUT, "/", cont.name, "_ma_plot.png"), width = 150, height = 150, units = "mm", res=300)
      try(dba.plotMA(db, contrast=cont))
    dev.off()
    png(paste0(OUT, "/", cont.name, "_boxplot_diff_sites.png"), width = 150, height = 150, units = "mm", res=300)
      try(dba.plotBox(db, contrast=cont))   # try, in case there are no significant peaks
    dev.off()
    png(paste0(OUT, "/", cont.name, "_volcano_plot.png"), width = 150, height = 150, units = "mm", res=300)
      rep <- try(dba.report(db, contrast=cont, th=1)) # scale dot size by Conc (max dot size set to 4)
      try(dba.plotVolcano(db, contrast=cont, dotSize=4*rep$Conc/max(rep$Conc)))
    dev.off()
    png(paste0(OUT, "/", cont.name, "_correlation_heatmap.png"), width = 150, height = 150, units = "mm", res=300)
      try(dba.plotHeatmap(db, contrast=cont))   # try, in case there are no significant peaks
    dev.off()
    png(paste0(OUT, "/", cont.name, "_pca_plot.png"), width = 150, height = 150, units = "mm", res=300)
      try(dba.plotPCA(db, contrast=cont, label=DBA_ID))   # try, in case there are no significant peaks
    dev.off()
    
    # plot consensus peaks
    vennmask <- dba.mask(db, DBA_CONDITION, factor(dba.show(db, bContrast=T)[cont,c("Group1", "Group2")]))
    if(sum(vennmask) <= 4) {   # if less than 2 replicates per group (or 4 replicates in total)
      png(paste0(OUT, "/", cont.name, "_venn_plot.png"), width = 150, height = 150, units = "mm", res=300)
        try(dba.plotVenn(db, mask=vennmask, main="Binding Site Overlaps Per Sample"))    # plot all relevant samples together
      dev.off()
    } else {  # generate a consensus peakset otherwise
       vennmask2 <- dba.mask(db2, DBA_CONDITION, factor(dba.show(db2, bContrast=T)[cont,c("Group1", "Group2")])) # db2 see above
       png(paste0(OUT, "/", cont.name, "_venn_plot.png"), width = 150, height = 150, units = "mm", res=300)
          try(dba.plotVenn(db2, main="Binding Site Overlaps Per Group",
                        mask=db2$masks$Consensus & names(db2$masks$Consensus) %in% factor(dba.show(db, bContrast=T)[cont,c("Group1", "Group2")]) ))
       dev.off()
    }

    tryCatch(dba.report(db, contrast=cont, bCalled=T), 
                    error=function(e) NULL) # dba.report crashes if there is exactly 1 significant hit to report
    
  })


##
## Annotate peaks
##
if(ANNOTATE) {
  library(ChIPseeker)
  txdb <- eval(parse(text=TXDB))
  result <- lapply(result, function(x) {
    tryCatch({
      x.ann <- annotatePeak(x, TxDb=txdb, annoDb=ANNODB, tssRegion=TSS, verbose=T)
      plotAnnoBar(x.ann)
      plotDistToTSS(x.ann)
      x.ann
    },
      error=function(e) NULL
    )
  })
}


writeLines(capture.output(sessionInfo()),paste(OUT, "/diffbind_session_info.txt", sep=""))
result <- lapply(result, as.data.frame)
names(result) <- substr(gsub("(.+)=\\((.+)\\)", "\\2", conts[,1]), 1, 31)
write.xlsx(result, file=paste0(OUT, "/diffbind.xlsx"))
saveRDS(result,  file=paste0(OUT, "/diffbind.rds"))
dba.save(db, dir=OUT, file='diffbind', pre="")

