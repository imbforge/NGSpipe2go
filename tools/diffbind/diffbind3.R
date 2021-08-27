#####################################
##
## What: diffbind3.R
## Who : Sergi Sayols, Frank Rühle
## When: 10-05-2016, updated 08-27-2021
##
## Script to perform differential binding analysis between 2 conditions (pairwise) using DiffBind v3
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
## librarySize=TRUE             # use total number of reads in bam for normalization (FALSE=only peaks)
## annotate=TRUE                # annotate after DB analysis?
## pe=FALSE                     # paired end experiment?
## tss=c(-3000,3000)            # region around the tss
## txdb=TxDb.Mmusculus.UCSC.mm9.knownGene   # Bioconductor transcript database, for annotation
## annodb=org.Mm.eg.db          # Bioconductor gene annotation database
##
######################################
options(stringsAsFactors=FALSE)
library(DiffBind) # load specific diffbind version #########################
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
DIFFBINDVERSION  <- parseArgs(args,"diffbindversion=", 3, "as.numeric") # selected DiffBind version
FTARGETS   <- parseArgs(args,"targets=","targets.txt")     # file describing the targets
FCONTRASTS <- parseArgs(args,"contrasts=","contrasts_diffbind.txt") # file describing the contrasts
CWD        <- parseArgs(args,"cwd=","./")     # current working directory
BAMS       <- parseArgs(args,"bams=",paste0(CWD, "/mapped"))  # directory with the bam files
PEAKS      <- parseArgs(args,"peaks=",paste0(CWD, "/results/macs2"))  # directory with the peak files
OUT        <- parseArgs(args,"out=", paste0(CWD, "/results/diffbind")) # directory where the output files will go
FRAGSIZE   <- parseArgs(args,"fragsize=", 200, "as.numeric")# fragment size
BLACKLIST  <- parseArgs(args,"blacklist=", FALSE, "as.logical") # if true, blacklist will be applied
GREYLIST   <- parseArgs(args,"greylist=", TRUE, "as.logical") # if true greylist will be generated for each Control
SUMMITS    <- parseArgs(args,"summits=", 200, "as.numeric") # summits for re-centering consensus peaks
FILTER     <- parseArgs(args,"filter=", 5, "as.numeric") # value to use for filtering intervals with low read counts
ANALYSISMETHOD  <- parseArgs(args,"analysisMethod=", "DESeq2", "as.character") # method for which to normalize (either "DESeq2" or "edgeRGLM")
LIBRARYSIZE     <- parseArgs(args,"librarySize=", "default", "as.character")   # use total number of reads in bam for normalization (FALSE=only peaks)
NORMALIZATION   <- parseArgs(args,"normalization=", "default", "as.character")   # use total number of reads in bam for normalization (FALSE=only peaks)
SUBSTRACTCONTROL<- parseArgs(args,"substractControl=", "default", "as.character")  # substract input
CONDITIONCOLUMN <- parseArgs(args,"conditionColumn=", "group", "as.character") # this targets column is interpreted as 'Condition' and is used as for defining the default design
FDR_TRESHOLD    <- parseArgs(args,"fdr_threshold=", 0.05, "as.numeric") # summits for re-centering consensus peaks
FOLD       <- parseArgs(args,"fold=", 0, "as.numeric") # summits for re-centering consensus peaks
ANNOTATE   <- parseArgs(args,"annotate=", TRUE, "as.logical") # annotate after DB analysis?
PE         <- parseArgs(args,"pe=", FALSE, "as.logical")      # paired end experiment?
TSS        <- parseArgs(args,"tss=", "c(-3000,3000)", "run_custom_code") # region around the tss
TXDB       <- parseArgs(args,"txdb=", "TxDb.Mmusculus.UCSC.mm9.knownGene") # Bioconductor transcript database, for annotation 
ANNODB     <- parseArgs(args,"annodb=", "org.Mm.eg.db") # Bioconductor gene annotation database
GENOMEDB   <- parseArgs(args,"genomedb=", "mm10") # Bioconductor gene annotation database

runstr <- "Rscript diffbind.R [targets=targets.txt] [contrasts=contrasts.txt] [cwd=./] [bams=./mapped] [peaks=./results/macs2] [out=./results] [fragsize=200] [summits=200] [annotate=TRUE] [pe=TRUE] [tss=c(-3000,3000)] [txdb=TxDb.Mmusculus.UCSC.mm9.knownGene] [annodb=org.Mm.eg.db]"
if(!file.exists(CWD))        stop("Dir",CWD,"does NOT exist. Run with:\n",runstr)
setwd(CWD)
if(!file.exists(FTARGETS))   stop("File",FTARGETS,"does NOT exist. Run with:\n",runstr)
if(!file.exists(FCONTRASTS)) stop("File",FCONTRASTS,"does NOT exist. Run with:\n",runstr)
if(!file.exists(BAMS))       stop("Dir",BAMS,"does NOT exist. Run with:\n",runstr)
if(!file.exists(PEAKS))      stop("Dir",PEAKS,"does NOT exist. Run with:\n",runstr)
if(!is.numeric(FRAGSIZE))    stop("Fragment size not numeric. Run with:\n",runstr)
if(!is.logical(BLACKLIST))   stop("blacklist not logical. Run with:\n",runstr)
if(!is.logical(GREYLIST))    stop("greylist not logical. Run with:\n",runstr)
if(!is.numeric(SUMMITS))     stop("Summits is not numeric. Run with:\n",runstr)
if(!is.numeric(FILTER))      stop("Filter threshold is not numeric. Run with:\n",runstr)
if(!is.numeric(FDR_TRESHOLD))      stop("FDR threshold is not numeric. Run with:\n",runstr)
if(!is.numeric(FOLD))        stop("Fold threshold is not numeric. Run with:\n",runstr)
if(!is.logical(ANNOTATE))    stop("Annotate not logical. Run with:\n",runstr)
if(!is.logical(PE))          stop("Paired end (pe) not logical. Run with:\n",runstr)
if(ANNOTATE & !is.numeric(TSS)) stop("Region around TSS not numeric. Run with:\n",runstr)
if(ANNOTATE & !require(TXDB, character.only=TRUE))   stop("Transcript DB", TXDB, "not installed\n")
if(ANNOTATE & !require(ANNODB, character.only=TRUE)) stop("Annotation DB", ANNODB, "not installed\n")

# check for DiffBind version
currentDiffbindVersion <- packageVersion('DiffBind')
DiffBindWarningText <- ""
if(currentDiffbindVersion < 3) {
  DiffBindWarningText <- "This script is made for diffbind version 3 and may not work properly for older versions."
  warning(paste0("You are using an outdated DiffBind version: ", currentDiffbindVersion, ". Beginning with version 3, DiffBind has included new functionalities and modified default settings. ", DiffBindWarningText))
  }


##
## make DB analysis
##

# load targets and make analysis
conts   <- read.delim(FCONTRASTS, head=T, stringsAsFactors = F, comment.char="#")
if(length(unique(conts$mmatrix))!=1) {stop("\nCan't use multiple design matrices in a single db object! Either unify the mmatrix entries in the contrast file or run diffbind in split mode.\n")}
targets_raw <- read.delim(FTARGETS, head=T, colClasses="character", comment.char="#")

# determine file suffixes for targets 
donefiles <- list.files(PEAKS,pattern=".done$")
bam_suffix <- sub("^[^\\.]*\\.", "", gsub("_macs2.done$", "", donefiles[1]))
bam_suffix <- paste0(bam_suffix, ".bam")

# check if blacklist filtering was applied
peakfiles <- list.files(PEAKS,pattern=".xls")
isBlacklistFilt <- any(grepl("blacklist_filtered", peakfiles)) 
peak_suffix <- if(isBlacklistFilt) {"_macs2_blacklist_filtered_peaks.xls"} else {"_macs2_peaks.xls"}

# check if any targets_raw$INPUT is indicated as "none". If so, Peak calling was done without Input samples.
isInputNone <- any(tolower(targets_raw$INPUT) == "none")

# create modified targets file for diffbind
targets <- data.frame(
  SampleID= targets_raw$IPname,
  Condition= targets_raw[,CONDITIONCOLUMN],
  Replicate= targets_raw$Replicate,
  bamReads= paste0(BAMS, "/", targets_raw$IP, ".", bam_suffix),
  ControlID= targets_raw$INPUTname,
  bamControl= paste0(BAMS, "/", targets_raw$INPUT, ".", bam_suffix),
  Peaks= paste0(PEAKS, "/", targets_raw$IPname, if(!isInputNone) {paste0(".vs.", targets_raw$INPUTname)}, peak_suffix),
  PeakCaller= targets_raw$PeakCaller)

  if(any(f<-colnames(targets_raw) %in% c("Tissue", "Factor", "Treatment"))) { # add additional factors allowed for modeling
    targets <- data.frame(targets, targets_raw[,f, drop=F])
  }
 
  if(isInputNone) {
  warning("You are running the DiffBind analysis without input control samples!\n")
  targets$ControlID <- NULL
  targets$bamControl <- NULL
  SUBSTRACTCONTROL <- FALSE
  }


# Construct DBA object
db <- dba(sampleSheet=targets, config=data.frame(AnalysisMethod=ANALYSISMETHOD, th=FDR_TRESHOLD, fragmentSize=FRAGSIZE,
                                                 doBlacklist=BLACKLIST, doGreylist=GREYLIST)) 

# create DBA object containing consensus peaks per group (needed later)
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


## process DBA object

  # apply black and grey lists and count reads
  tryGL <- try(dba.blacklist(db, blacklist=db$config$doBlacklist, greylist=db$config$doGreylist))  # crashes if greylist cannot be build
  if(class(tryGL)=="try-error") {
    db <- dba.blacklist(db, blacklist=db$config$doBlacklist, greylist=F) # apply blacklist only (doesn't crash if no blacklist found)
    db$config$doGreylist <- F
    cat("\nBuilding of greylist failed, analysis is performed without greylist.\n")
  } else {
    db <- tryGL
  }
  if(SUBSTRACTCONTROL=="default") {SUBSTRACTCONTROL_FINAL <- is.null(db$greylist)} else {SUBSTRACTCONTROL_FINAL <- as.logical(SUBSTRACTCONTROL)}
  blacklist_generated <- greylist_generated <- NULL
  if(BLACKLIST){blacklist_generated <- try(dba.blacklist(db, Retrieve=DBA_BLACKLIST))}
  if(GREYLIST) {greylist_generated  <- try(dba.blacklist(db, Retrieve=DBA_GREYLIST))}

    
  # identify all overlapping peaks and derives a consensus peakset for the experiment. 
  # Then count how many reads overlap each interval for each unique sample.
  db <- dba.count(db, minOverlap=2, score=DBA_SCORE_NORMALIZED, summits=SUMMITS, filter=FILTER, 
                  bScaleControl=TRUE, minCount=0, bUseSummarizeOverlaps=TRUE) 
  # score: which score to use in the binding affinity matrix. Note that all raw read counts are maintained for use by dba.analyze, 
  # regardless of how the score is set here. 
  # DBA_SCORE_NORMALIZED: normalized reads, as set by dba.normalize
 
  # retrieve dba object info
  infodb <- dba.show(db)
  infodb$Intervals <- NULL
  infodb$Caller <- NULL
  names(infodb)[names(infodb) == "Reads"] <- "fullLibSize"
  infodb$ReadPeaks <- round(infodb$fullLibSize*infodb$FRiP)

  
  # Normalizing the data
  db <- dba.normalize(db, method = db$config$AnalysisMethod, normalize = NORMALIZATION, library = LIBRARYSIZE, 
                      offsets = FALSE, bSubControl = SUBSTRACTCONTROL_FINAL, background = FALSE)
  # The major change in DiffBind version 3.0 is in how the data are modeled. In previous versions, a separate model was derived for each contrast, 
  # including data only for those samples present in the contrast. Model design options were implicit and limited to either a single factor, or 
  # a subset of two-factor "blocked" designs. Starting in version 3.0, the default mode is to include all the data in a single model, allowing 
  # for any allowable design formula and any set of allowable contrasts. 

  # retrieve normalization info with bRetrieve=TRUE
  infoNorm <- dba.normalize(db, method = db$config$AnalysisMethod, bRetrieve=TRUE)
  infodb$normFactors <- infoNorm$norm.factors
  infodb$normLibSize <- round(infodb$fullLibSize/infodb$normFactors)
  

  # apply the contrasts
   for (i in 1:nrow(conts)) {
   # parse formula in cont
    cont.name <- conts[i,1]
    cont.form <- conts[i,2]
    mmatrix   <- conts[i,3]
    mmatrix <- gsub(CONDITIONCOLUMN, "Condition", mmatrix)
    factors   <- gsub("(^\\s+|\\s+$)", "", unlist(strsplit(cont.form,"\\W")))
    factors   <- factors[factors != ""]
    contrast <- c("Condition", factors)
    if(length(factors) != 2) {
      warning(paste(conts[i,],"cannot deal with designs other than pairwise comparisons!"))
      return(NA)
    }
    db <- dba.contrast(db, design=mmatrix, contrast=contrast)
   }


# run the diffbind analysis (DESeq2) for all contrasts
  db <- dba.analyze(db) 


# PCA plot
png(paste0(OUT, "/pca_plot_all_samples.png"), width = 150, height = 150, units = "mm", res=300)
  dba.plotPCA(db, DBA_CONDITION, label=DBA_ID)
dev.off()

# Heatmap using affinity (read count) data
png(paste0(OUT, "/heatmap_affinity_scores.png"), width = 150, height = 150, units = "mm", res=300)
  dba.plotHeatmap(db, main="Correlation heatmap\naffinity scores")
dev.off()


# Boxplot of consensus peaks
  lpeaks <- db2$peaks # db2 contains group consensus peaksets additional to the sample peaksets
  dfcolnames <- colnames(lpeaks[[1]]) # Consensus peak dfs don't have colnames, therefore name columns 
  lpeaks <- lapply(1:length(lpeaks), function(x) setNames(lpeaks[[x]], dfcolnames) )
  
  names(lpeaks) <- names(db2$masks$Consensus) # name list elements
  names(lpeaks)[db2$masks$Consensus] <- paste0("Consensus\n", names(lpeaks)[db2$masks$Consensus])
  
  lpeaks <- lapply(lpeaks, function(x) { # calculate peak sizes
    x$width <- abs(x[,grep("end", dfcolnames, ignore.case = T)] - x[,grep("start", dfcolnames, ignore.case = T)]) # columns names may upper or lower case
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
    cont.name <- conts[cont,1]

    png(paste0(OUT, "/", cont.name, "_ma_plot.png"), width = 150, height = 150, units = "mm", res=300)
      try(dba.plotMA(db, contrast=cont, fold=FOLD))
    dev.off()
    png(paste0(OUT, "/", cont.name, "_boxplot_diff_sites.png"), width = 150, height = 150, units = "mm", res=300)
      try(dba.plotBox(db, contrast=cont))   # try, in case there are no significant peaks
    dev.off()
    png(paste0(OUT, "/", cont.name, "_volcano_plot.png"), width = 150, height = 150, units = "mm", res=300)
      rep <- try(dba.report(db, contrast=cont, th=1)) # scale dot size by Conc (max dot size set to 4)
      try(dba.plotVolcano(db, contrast=cont, fold=FOLD, dotSize=4*rep$Conc/max(rep$Conc)))
    dev.off()
    png(paste0(OUT, "/", cont.name, "_correlation_heatmap.png"), width = 150, height = 150, units = "mm", res=300)
      try(dba.plotHeatmap(db, contrast=cont))   # try, in case there are no significant peaks
    dev.off()
    png(paste0(OUT, "/", cont.name, "_pca_plot.png"), width = 150, height = 150, units = "mm", res=300)
      try(dba.plotPCA(db, contrast=cont, label=DBA_ID))   # try, in case there are no significant peaks
    dev.off()
    
    # plot consensus peaks
    # 1st group name between diffbind versions
    vennmask <- dba.mask(db, DBA_CONDITION, factor(dba.show(db, bContrast=T)[cont, grep("group", names(dba.show(db, bContrast=T)), ignore.case = T)])) 
    if(sum(vennmask) <= 4) {   # if less than 2 replicates per group (or 4 replicates in total)
      png(paste0(OUT, "/", cont.name, "_venn_plot.png"), width = 150, height = 150, units = "mm", res=300)
        try(dba.plotVenn(db, mask=vennmask, main="Binding Site Overlaps Per Sample"))    # plot all relevant samples together
      dev.off()
    } else {  # generate a consensus peakset otherwise
       vennmask2 <- dba.mask(db2, DBA_CONDITION, factor(dba.show(db2, bContrast=T)[cont, grep("group", names(dba.show(db2, bContrast=T)), ignore.case = T)])) # db2 see above
       png(paste0(OUT, "/", cont.name, "_venn_plot.png"), width = 150, height = 150, units = "mm", res=300)
          try(dba.plotVenn(db2, main="Binding Site Overlaps Per Group",
                        mask=db2$masks$Consensus & names(db2$masks$Consensus) %in% factor(dba.show(db, bContrast=T)[cont, grep("group", names(dba.show(db, bContrast=T)), ignore.case = T)]) ))
       dev.off()
    }

    tryCatch(dba.report(db, contrast=cont, bCalled=T, bUsePval=F, th=1, fold=0), # th=db$config$th, fold=FOLD # filtering is done later
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

  
  # create context dependent parameter table   
  explLibrarySize <- switch(infoNorm$lib.method,
                            full="Use the full library size for normalization",
                            RiP="Use the number of reads that overlap consensus peaks for normalization",
                            background="Use the total number of reads aligned to the chromosomes for which there is at least one peak (requires background bin calculation)"
                            )
  if(LIBRARYSIZE=="default") {explLibrarySize <- paste0("refers to method '", infoNorm$lib.method, "'. ", explLibrarySize)}
  
  explNormalization <- switch(infoNorm$norm.method,
                              RLE="RLE normalization (native to DESeq2)",
                              TMM="TMM normalization (native to EDGER)",
                              native="Use native method based on analysis method: RLE for DESeq2 or TMM for EDGER",
                              lib="Normalize by library size only",
                              offsets="Indicates that offsets have been specified using the offsets parameter, and they should be used without alteration",
                              "adjust offsets"="Indicates that offsets have been specified using the offsets parameter, and they should be adjusted for library size and mean centering before being used in a DESeq2 analysis"
                              )
  if(NORMALIZATION=="default") {explNormalization <- paste0("refers to method '", infoNorm$norm.method, "'. ", explNormalization)}
  
  # create overview table with diffbind settings
  diffbindSettings <- rbind(c(Parameter="DiffBind package version", Value=as.character(currentDiffbindVersion), Comment= paste(if(currentDiffbindVersion$major != floor(DIFFBINDVERSION)) {paste0("Major version not concordant with intended DiffBind v", DIFFBINDVERSION, ". Please check R module." )} else {""}, DiffBindWarningText)),
                            c("Apply auto-generated blacklist", BLACKLIST, if(BLACKLIST){if(class(blacklist_generated)=="try-error") {"Skipped because blacklist not available."} else {"Blacklist successfully generatedand applied."}} else {""}),
                            c("Apply auto-generated greylist", GREYLIST, if(GREYLIST){if(class(greylist_generated)=="try-error") {"Skipped because greylist not available."} else {"Greylist successfully generated and applied."}} else {""}),
                            c("Fragment size", FRAGSIZE, ""),
                            c("Summits", SUMMITS, if(SUMMITS==0) {"no re-centering of peaks."} else {paste0("re-center peaks around consensus summit with peak width 2x", SUMMITS, ".")}),
                            c("Filter threshold", FILTER, "threshold for filtering intervals with low read counts."),
                            c("Analysis method", ANALYSISMETHOD, ""),
                            c("Library size", LIBRARYSIZE, explLibrarySize),
                            c("Normalization method", NORMALIZATION, explNormalization),
                            c("Subtract control read counts", SUBSTRACTCONTROL, if(SUBSTRACTCONTROL=="default") {paste("set to", SUBSTRACTCONTROL_FINAL, "because greylist is", if(SUBSTRACTCONTROL_FINAL){"not"} else {""}, "available.")}),
                            c("Design", "contrast_diffbind.txt", paste(db$design, "(with", CONDITIONCOLUMN, "as Condition column)")),
                            c("FDR threshold", FDR_TRESHOLD, "significance threshold for differential binding analysis."),
                            c("Fold threshold", FOLD, "log Fold threshold for differential binding analysis.")
  )
  

writeLines(capture.output(sessionInfo()),paste(OUT, "/diffbind_session_info.txt", sep=""))
write.table(diffbindSettings, file=file.path(OUT, "diffbind_settings.txt"), row.names = F, quote = F, sep="\t")
write.table(infodb, file=file.path(OUT, "info_dba_object.txt"), row.names = F, quote = F, sep="\t")
result <- lapply(result, as.data.frame)
names(result) <- conts$contrast.name
write.xlsx(result, file=paste0(OUT, "/diffbind_all_sites.xlsx"))
result <- lapply(result, function(x) {x[x$FDR<=db$config$th & abs(x$Fold)>=FOLD, ]}) # filter result tables for significance
write.xlsx(result, file=paste0(OUT, "/diffbind.xlsx"))
saveRDS(result,  file=paste0(OUT, "/diffbind.rds"))
dba.save(db, dir=OUT, file='diffbind', pre="")

