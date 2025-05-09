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
library(parallel)
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
FTARGETS   <- parseArgs(args,"targets=","targets.txt")     # file describing the targets
FCONTRASTS <- parseArgs(args,"contrasts=","contrasts_diffbind.txt") # file describing the contrasts
CWD        <- parseArgs(args,"cwd=","./")     # current working directory
BAMS       <- parseArgs(args,"bams=",paste0(CWD, "/mapped"))  # directory with the bam files
PEAKS      <- parseArgs(args,"peaks=",paste0(CWD, "/results/macs2"))  # directory with the peak files
OUT        <- parseArgs(args,"out=", paste0(CWD, "/results/diffbind")) # directory where the output files will go
FRAGSIZE   <- parseArgs(args,"fragsize=", 200, "as.numeric") # fragment size
CORES      <- parseArgs(args,"cores=", 4, "as.numeric") # fragment size
BLACKLIST  <- parseArgs(args,"blacklist=", FALSE, "as.logical") # if true, blacklist will be applied
GREYLIST   <- parseArgs(args,"greylist=", FALSE, "as.logical") # if true greylist will be generated for each Control
SUMMITS    <- parseArgs(args,"summits=", 200, "as.numeric") # summits for re-centering consensus peaks
FILTER     <- parseArgs(args,"filter=", 1, "as.numeric") # value to use for filtering intervals with low read counts
MINOVERLAP <- parseArgs(args,"minOverlap=", 2, "as.numeric") # value to use for filtering intervals with low read counts
ANALYSISMETHOD  <- parseArgs(args,"analysisMethod=", "DESeq2", "as.character") # method for which to normalize (either "DESeq2" or "edgeRGLM")
LIBRARYSIZE     <- parseArgs(args,"librarySize=", "default", "as.character")   # use total number of reads in bam for normalization (FALSE=only peaks)
NORMALIZATION   <- parseArgs(args,"normalization=", "default", "as.character")   # use total number of reads in bam for normalization (FALSE=only peaks)
BACKGROUND      <- parseArgs(args,"background=", FALSE, "as.logical") # background may either be logical (default binsize 15000) or numeric with custom binsize.
if(is.na(BACKGROUND)) {BACKGROUND <- parseArgs(args,"background=", "15000", "as.numeric")}
SUBSTRACTCONTROL<- parseArgs(args,"substractControl=", "FALSE", "as.character")  # subtract input
CONDITIONCOLUMN <- parseArgs(args,"conditionColumn=", "group", "as.character") # this targets column is interpreted as 'Condition' and is used as for defining the default design
FDR_TRESHOLD    <- parseArgs(args,"fdr_threshold=", 0.05, "as.numeric") # summits for re-centering consensus peaks
LFC       <- parseArgs(args,"lfc=", 0, "as.numeric") # summits for re-centering consensus peaks
ANNOTATE   <- parseArgs(args,"annotate=", FALSE, "as.logical") # annotate after DB analysis?
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
if(!isTRUE(BLACKLIST))       {BLACKLIST <- FALSE} # if a specific blacklist is defined in essential.vars it is applied by the blacklist_filter module and not here. If set to TRUE in diffbind3.header diffbind applies a public blacklist if available. 
if(!is.logical(GREYLIST))    stop("greylist not logical. Run with:\n",runstr)
if(!is.numeric(SUMMITS))     stop("Summits is not numeric. Run with:\n",runstr)
if(!is.numeric(FILTER))      stop("Filter threshold is not numeric. Run with:\n",runstr)
if(!is.numeric(MINOVERLAP))  stop("minOverlap threshold is not numeric. Run with:\n",runstr)
if(!is.numeric(FDR_TRESHOLD))      stop("FDR threshold is not numeric. Run with:\n",runstr)
if(!is.numeric(LFC))         stop("Log2 fold change threshold is not numeric. Run with:\n",runstr)
if(!is.logical(ANNOTATE))    stop("Annotate not logical. Run with:\n",runstr)
if(!is.logical(PE))          stop("Paired end (pe) not logical. Run with:\n",runstr)
if(ANNOTATE & !is.numeric(TSS)) stop("Region around TSS not numeric. Run with:\n",runstr)
if(ANNOTATE & !require(ANNODB, character.only=TRUE)) stop("Annotation DB", ANNODB, "not installed\n")

if(ANNOTATE) {
  if(grepl("\\.gtf$", TXDB)){ # check the input format for the transcript annotation
    library(GenomicFeatures)
    txdb <- makeTxDbFromGFF(TXDB, format="gtf") # if the input format is gtf file, then this file will be used to create a TxDb object
  } else {
    if(!require(TXDB, character.only=TRUE)) stop("Transcript DB", TXDB, "not installed\n")
    txdb <- eval(parse(text=TXDB)) # if the input format is bioconductor, then the transcript annotation library will be used 
  }
}

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
contsAll   <- read.delim(FCONTRASTS, head=T, stringsAsFactors = F, comment.char="#", sep="\t")
if(length(unique(contsAll$mmatrix))!=1) {stop("\nCan't use multiple design matrices in a single db object! Either unify the mmatrix entries in the contrast file or run diffbind in split mode.\n")}
targets_raw <- read.delim(FTARGETS, head=T, colClasses="character", comment.char="#", sep="\t")

# determine file suffixes for targets 
donefiles <- list.files(PEAKS,pattern=".done$")
bam_suffix <- sub("^[^\\.]*\\.*", "", gsub("_macs2.done$", "", donefiles[1]))
bam_suffix <- ifelse(bam_suffix == "", paste0(bam_suffix, "bam"), paste0(bam_suffix, ".bam"))

# check if blacklist filtering was applied and if broad peak files are available
peakfiles <- list.files(PEAKS,pattern="(\\.xls$)|(Peak$)")
isBlacklistFilt <- any(grepl("blacklist_filtered", peakfiles)) 
isBroad <- any(grepl("broadPeak$", peakfiles)) 
peak_suffix <- sapply(targets_raw$PeakCaller, function(x) {switch(x,
                                                                  macs=if(isBlacklistFilt) {"_macs2_blacklist_filtered_peaks.xls"} else {"_macs2_peaks.xls"},
                                                                  bed= paste0("_macs2_peaks", if(isBlacklistFilt) {"_blacklist_filtered"}, if(isBroad) {".broadPeak"} else {".narrowPeak"})
)})

# check if any targets_raw$INPUT is indicated as "none". If so, Peak calling was done without Input samples.
isInputNone <- any(tolower(targets_raw$INPUT) == "none")

# create modified targets file for diffbind
targetsAll <- data.frame(
  SampleID= targets_raw$IPname,
  Condition= targets_raw[,CONDITIONCOLUMN],
  Replicate= targets_raw$Replicate,
  bamReads= paste0(BAMS, "/", targets_raw$IP, ".", bam_suffix),
  ControlID= targets_raw$INPUTname,
  bamControl= paste0(BAMS, "/", targets_raw$INPUT, ".", bam_suffix),
  Peaks= paste0(PEAKS, "/", targets_raw$IPname, if(!isInputNone) {paste0(".vs.", targets_raw$INPUTname)}, peak_suffix),
  PeakCaller= targets_raw$PeakCaller)

  if(any(f<-colnames(targets_raw) %in% c("Tissue", "Factor", "Treatment"))) { # add additional factors allowed for modeling
    targetsAll <- data.frame(targetsAll, targets_raw[,f, drop=F])
  }
 
if(isInputNone | (!is.na(as.logical(SUBSTRACTCONTROL)) & !as.logical(SUBSTRACTCONTROL))) {
  warning("You are running the DiffBind analysis without input control samples!\n")
  # if SUBSTRACTCONTROL=FALSE controls are set to NULL as well because the bSubControl parameter
  # in dba.count seems to be of no effect.
    targetsAll$ControlID <- NULL
    targetsAll$bamControl <- NULL
    SUBSTRACTCONTROL <- FALSE
  }



##
## diffbind crashes if peaksets with no peaks are included for generating consensus peakset
# remove targets which have no peaks
peaks <- sapply(targetsAll$Peaks, function(x) {
  tryCatch({
    nrow(read.delim(x, head=TRUE, comment="#"))
  }, error=function(e) 0)
})

if(!all(peaks > 0)) {
  warning("Sample(s) ", paste(targetsAll$SampleID[!(peaks > 0)], collapse=", "),
          " excluded from Diffbind because didn't have any peaks called")
  targetsAll <- targetsAll[peaks > 0, ] 
}

# and contrasts with "orphan" targets or groups containing only 1 replicate 
valid_conts <- sapply(strsplit(contsAll$contrast, "\\W"), function(factors) {
  factors <- gsub("(^\\s+|\\s+$)", "", factors)  # remove leading spaces
  factors <- factors[factors != ""]
  #all(factors %in% targetsAll$Condition)
  all(sapply(factors, function(x) {nrow(targetsAll[targetsAll$Condition==x,])>1}))
})

if(!all(valid_conts)) {
  warning("Contrast(s) ", paste(contsAll$contrast.name[!valid_conts], collapse=", "),
          " excluded from Diffbind because groups included which didn't have more than one valid peaksets")
  contsAll <- contsAll[valid_conts, ]
}


### start loop for sub_experiment 
result <- list()
minOverlapError <- c() # stores if dba.peakset crashes
if(any(is.null(contsAll$sub_experiment))){
  stop("Please specify sub_experiments in contrasts_diffbind.txt")
}
for (sub in unique(contsAll$sub_experiment)) {

  if(length(unique(contsAll$sub_experiment))>1) {
    print(paste("process sub_experiment", sub))
    subexpPrefix <- paste0("SubExp_", sub, "_")
  } else {
    subexpPrefix <- ""
  }
  
  conts <- contsAll[contsAll$sub_experiment==sub,]  
  groups2use   <- gsub("(^\\s+|\\s+$)", "", unlist(strsplit(conts$contrast,"\\W")))
  groups2use   <- groups2use[groups2use != ""]
  targets <- targetsAll[targetsAll$Condition %in% groups2use,]
  
  
  # Construct DBA object
  db <- dba(sampleSheet=targets, config=data.frame(AnalysisMethod=ANALYSISMETHOD, th=FDR_TRESHOLD, fragmentSize=FRAGSIZE,
                                                   RunParallel=TRUE, cores=CORES,
                                                   doBlacklist=BLACKLIST, doGreylist=GREYLIST)) 

  # create DBA object with added consensus peak sets PER GROUP (1 consensus peakset per group, this is needed later for plots)
  tryConsensus <- try(db2 <- dba.peakset(DBA=db, consensus=DBA_CONDITION, minOverlap=MINOVERLAP))  # crashes if minOverlap is too high for a certain group
  if(class(tryConsensus)=="try-error") {
    minOverlapError <- c(minOverlapError, sub)
    warning(paste("minOverlap=", MINOVERLAP, "in db2 object failed. Consider to reduce this value. This would effect the entire downstream analysis."))
    # db2 affects the Boxplot of consensus peaks and the consensus peakset venn diagram. 
    # But MINOVERLAP would also affect the dba.count call for db (when generating the consensus peakset for all samples).
  }
  
  # Heatmap using occupancy (peak caller score) data
  png(file.path(OUT, paste0(subexpPrefix, "heatmap_occupancy.png")), width = 150, height = 150, units = "mm", res=300)
  dba.plotHeatmap(db, main="Correlation heatmap\noccupancy data")
  dev.off()
  
  # plot Overlap rate
  png(file.path(OUT, paste0(subexpPrefix, "Overlap_rate_plot.png")), width = 150, height = 150, units = "mm", res=300)
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
    db <- dba.count(db, minOverlap=MINOVERLAP, score=DBA_SCORE_NORMALIZED, summits=SUMMITS, filter=FILTER, 
                    bScaleControl=TRUE, bSubControl = SUBSTRACTCONTROL_FINAL, minCount=0, bUseSummarizeOverlaps=TRUE) 
    # score: which score to use in the binding affinity matrix. Note that all raw read counts are maintained for use by dba.analyze, 
    # regardless of how the score is set here. 
    # DBA_SCORE_NORMALIZED: normalized reads, as set by dba.normalize
   
    # retrieve dba object info
    infodb <- dba.show(db)
    infodb$Intervals <- NULL
    infodb$Caller <- NULL
    names(infodb)[names(infodb) == "Reads"] <- "fullLibSize"
    if(all(c("fullLibSize", "FRiP") %in% names(infodb))) {infodb$ReadPeaks <- round(infodb$fullLibSize*infodb$FRiP)}
    
    # Normalizing the data
    db <- dba.normalize(db, method = db$config$AnalysisMethod, normalize = NORMALIZATION, library = LIBRARYSIZE, 
                        offsets = FALSE, background = BACKGROUND) 
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
  png(file.path(OUT, paste0(subexpPrefix, "pca_plot_all_samples.png")), width = 150, height = 150, units = "mm", res=300)
    dba.plotPCA(db, DBA_CONDITION, label=DBA_ID)
  dev.off()
  
  # Heatmap using affinity (read count) data
  png(file.path(OUT, paste0(subexpPrefix, "heatmap_affinity_scores.png")), width = 150, height = 150, units = "mm", res=300)
    dba.plotHeatmap(db, main="Correlation heatmap\naffinity scores")
  dev.off()
  
  
  # Boxplot of consensus peaks
  if(class(tryConsensus)=="DBA") { # check if db2 was successfully created
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
      
    p <- ggplot(plotdata, aes(x=group, y=width, fill=group)) + geom_boxplot(show.legend = FALSE) + 
          theme(text = element_text(size=15), axis.text.x = element_text(angle=90, vjust=0.5)) + 
          labs(title="Size of consensus peaks") + theme(plot.title = element_text(hjust = 0.5)) 
    ggsave(plot=p, filename=file.path(OUT, paste0(subexpPrefix, "peak_width_boxplot.png")))
  }
    
  
  # Venn diagram of all contrasts
  png(file.path(OUT, paste0(subexpPrefix, "venn_plot_all_contrasts.png")), width = 150, height = 150, units = "mm", res=300)
    try(dba.plotVenn(db, main="Binding Site Overlaps Per Contrast", contrast=1:nrow(dba.show(db, bContrast=T))))
  dev.off()
  
  
  # prepare results and plots for each contrast
    result[[sub]] <- lapply(1:nrow(conts), function(cont) {
      cont.name <- conts[cont,1]
  
      png(file.path(OUT, paste0(subexpPrefix, cont.name, "_ma_plot.png")), width = 150, height = 150, units = "mm", res=300)
        try(dba.plotMA(db, contrast=cont, fold=LFC))
      dev.off()
      png(file.path(OUT, paste0(subexpPrefix, cont.name, "_boxplot_diff_sites.png")), width = 150, height = 150, units = "mm", res=300)
        try(dba.plotBox(db, contrast=cont))   # try, in case there are no significant peaks
      dev.off()
      png(file.path(OUT, paste0(subexpPrefix, cont.name, "_volcano_plot.png")), width = 150, height = 150, units = "mm", res=300)
        rep <- try(dba.report(db, contrast=cont, th=1)) # scale dot size by Conc (max dot size set to 4)
        try(dba.plotVolcano(db, contrast=cont, fold=LFC, dotSize=4*rep$Conc/max(rep$Conc)))
      dev.off()
      png(file.path(OUT, paste0(subexpPrefix, cont.name, "_correlation_heatmap.png")), width = 150, height = 150, units = "mm", res=300)
        try(dba.plotHeatmap(db, contrast=cont))   # try, in case there are no significant peaks
      dev.off()
      png(file.path(OUT, paste0(subexpPrefix, cont.name, "_pca_plot.png")), width = 150, height = 150, units = "mm", res=300)
        try(dba.plotPCA(db, contrast=cont, label=DBA_ID))   # try, in case there are no significant peaks
      dev.off()
      
      # plot consensus peaks
      # 1st group name between diffbind versions
      vennmask <- dba.mask(db, DBA_CONDITION, factor(as.character(dba.show(db, bContrast=T)[cont, grep("group", names(dba.show(db, bContrast=T)), ignore.case = T)])))
      if(sum(vennmask) <= 4) {   # if less than 2 replicates per group (or 4 replicates in total)
        png(file.path(OUT, paste0(subexpPrefix, cont.name, "_venn_plot.png")), width = 150, height = 150, units = "mm", res=300)
          try(dba.plotVenn(db, mask=vennmask, main="Binding Site Overlaps Per Sample"))    # plot all relevant samples together
        dev.off()
      } else {  # generate a consensus peakset otherwise
        if(class(tryConsensus)=="DBA") { # check if db2 was successfully created
          vennmask2 <- dba.mask(db2, DBA_CONDITION, factor(as.character(dba.show(db2, bContrast=T)[cont, grep("group", names(dba.show(db2, bContrast=T)), ignore.case = T)]))) # db2 see above
           png(file.path(OUT, paste0(subexpPrefix, cont.name, "_venn_plot.png")), width = 150, height = 150, units = "mm", res=300)
              try(dba.plotVenn(db2, main="Binding Site Overlaps Per Group",
                            mask=db2$masks$Consensus & names(db2$masks$Consensus) %in% factor(dba.show(db, bContrast=T)[cont, grep("group", names(dba.show(db, bContrast=T)), ignore.case = T)]) ))
           dev.off()
        }
      }
  
      tryCatch({
        x=dba.report(db, contrast=cont, bCalled=T, bUsePval=F, th=1, fold=0) # th=db$config$th, fold=LFC # filtering is done later
        colnames(mcols(x)) <- plyr::revalue(colnames(mcols(x)), c("Fold" = paste0("lfc_", cont.name), "FDR" = "BHadj_pvalue")) # rename column names
        x <- x[,!colnames(mcols(x)) %in% c("p-value")] # remove column with unadjusted p-value
        return(x)
      }, error=function(e) GRanges()) # dba.report crashes if there is exactly 1 significant hit to report
      
    })
  
    names(result[[sub]]) <- paste0(subexpPrefix, conts$contrast.name)
    result[[sub]] <- result[[sub]][sapply(result[[sub]], length) !=0] # remove empty GRanges if present
    
  ##
  ## Annotate peaks
  ##
  if(ANNOTATE) {
    library(ChIPseeker)
    result[[sub]] <- lapply(result[[sub]], function(x) {
      tryCatch({
        x.ann <- ChIPseeker::annotatePeak(x, TxDb=txdb, annoDb=ANNODB, tssRegion=TSS, verbose=T)
      },
        error=function(e) NULL
      )
    })
  }

dba.save(db, dir=OUT, file='diffbind', pre=subexpPrefix)
write.table(infodb, file=file.path(OUT, paste0(subexpPrefix, "info_dba_object.txt")), row.names = F, quote = F, sep="\t")

} # end sub loop
  
  # create context dependent parameter table   
  explLibrarySize <- switch(infoNorm$lib.method,
                            full="Use the full library size for normalization.",
                            RiP="Use the number of reads that overlap consensus peaks for normalization.",
                            background="Use the total number of reads aligned to the chromosomes for which there is at least one peak (requires background bin calculation)."
                            )
  if(LIBRARYSIZE=="default") {explLibrarySize <- paste0("refers to method '", infoNorm$lib.method, "'. ", explLibrarySize)}
  
  explNormalization <- switch(infoNorm$norm.method,
                              RLE="RLE normalization (native to DESeq2).",
                              TMM="TMM normalization (native to EDGER).",
                              native="Use native method based on analysis method: RLE for DESeq2 or TMM for EDGER.",
                              lib="Normalize by library size only.",
                              offsets="Indicates that offsets have been specified using the offsets parameter, and they should be used without alteration.",
                              "adjust offsets"="Indicates that offsets have been specified using the offsets parameter, and they should be adjusted for library size and mean centering before being used in a DESeq2 analysis."
                              )
  if(NORMALIZATION=="default") {explNormalization <- paste0("refers to method '", infoNorm$norm.method, "'. ", explNormalization)}

  explBackground <- if(isTRUE(infoNorm$background)) {
                          if(is.numeric(BACKGROUND)) {
                            paste("background bins used for normalization coumputed with bin size", BACKGROUND, "bp.")} else {
                            "background bins used for normalization coumputed with default bin size 15000 bp."}
                        } else {
                            if(is.numeric(infoNorm$background)) { # if BACKGROUND is numeric, this may be stored here in future DiffBind versions
                            paste("background bins used for normalization coumputed with bin size", infoNorm$background, "bp.")
                            } else {"no background bins computed."}
                        }
  

  # create overview table with diffbind settings
  diffbindSettings <- rbind(c(Parameter="DiffBind package version", Value=as.character(currentDiffbindVersion), Comment=DiffBindWarningText),
                            c("external blacklist applied", isBlacklistFilt, if(isBlacklistFilt) {"MACS2 peak files have already been blacklist or greylist filtered"} else {""}),
                            c("Apply auto-generated blacklist", BLACKLIST, if(BLACKLIST){if(class(blacklist_generated)=="try-error") {"Skipped because blacklist not available."} else {"Blacklist successfully generatedand applied."}} else {""}),
                            c("Apply auto-generated greylist", GREYLIST, if(GREYLIST){if(class(greylist_generated)=="try-error") {"Skipped because greylist not available."} else {"Greylist successfully generated and applied."}} else {""}),
                            c("Broad peaks used", isBroad, paste("loaded peak files with suffix:", unique(peak_suffix), collapse = ", ")),
                            c("Fragment size", FRAGSIZE, ""),
                            c("Summits", SUMMITS, if(SUMMITS==0) {"no re-centering of peaks."} else {paste0("re-center peaks around consensus summit with peak width 2x", SUMMITS, ".")}),
                            c("Filter threshold", FILTER, "threshold for filtering intervals with low read counts."),
                            c("min Overlap", MINOVERLAP, if(length(minOverlapError) >0) {paste("Per-group consensus peaksets not available for subexperiment", paste(minOverlapError, collapse=", "), "(only relevant for some of the plots).")} else {paste("include peak in consensus if present in at least", MINOVERLAP, "peaksets.")} ),
                            c("Analysis method", ANALYSISMETHOD, ""),
                            c("Library size", LIBRARYSIZE, explLibrarySize),
                            c("Normalization method", NORMALIZATION, explNormalization),
                            c("Background", BACKGROUND, explBackground),
                            c("Subtract control read counts", SUBSTRACTCONTROL, if(SUBSTRACTCONTROL=="default") {paste("set to", SUBSTRACTCONTROL_FINAL, "because greylist is", if(SUBSTRACTCONTROL_FINAL){"not"} else {""}, "available.")} else {""}),
                            c("Design", "contrast_diffbind.txt", paste(db$design, "(with", CONDITIONCOLUMN, "as Condition column).")),
                            c("FDR threshold", FDR_TRESHOLD, "significance threshold for differential binding analysis."),
                            c("LFC threshold", LFC, "log2 fold change threshold for differential binding analysis.")
  )
  

writeLines(capture.output(sessionInfo()),paste(OUT, "/diffbind_session_info.txt", sep=""))
write.table(diffbindSettings, file=file.path(OUT, "diffbind_settings.txt"), row.names = F, quote = F, sep="\t")
result <- unlist(result, recursive = F) # flatten the nested list
result <- lapply(result, as.data.frame)
write.xlsx(result, file=paste0(OUT, "/diffbind_all_sites.xlsx"))
result <- lapply(result, function(x) {
  tryCatch({
    x[x$BHadj_pvalue<=db$config$th & abs(x[,grep("^lfc_", colnames(x))])>=LFC, ]   # filter result tables for significance
  }, error=function(e) x)
})
write.xlsx(result, file=paste0(OUT, "/diffbind.xlsx"))
saveRDS(result,  file=paste0(OUT, "/diffbind.rds"))
print(warnings())
