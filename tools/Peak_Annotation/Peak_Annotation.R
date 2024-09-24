#############################################################################
##
## What: Peak_Annotation.R
## Who: 
## When: 
##
## Script to annotate ChIPseq peaks using result from MACS2
##
## Args: 
## -----
## peakData=                                        # path to the xls file result from MACS2
## transcriptType=Bioconductor                      # define the transcript annotation type for the analysis
## transcriptDb=TxDb.Hsapiens.UCSC.hg19.knownGene   # either GTF file or transcript annotation database as transcript annotation database
## orgDb=org.Hs.eg.db                               # optionally genome wide annotation 
## regionTSS=3000                                   # TSS region parameter from -3000 to +3000
## targets=targets.txt                              # tab separated file describing the samples to be processed
## orderby=group                                    # column in targets to order samples in plots
## out=Peak_Annotation                              # output directory
############################################################################
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')
library(ChIPseeker)
library(clusterProfiler)
library(GenomicFeatures)
library(openxlsx)

##
# get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
 
   if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  if(!is.null(default)) default else do.call(convert,list(NA))
  
}

args           <- commandArgs(T)
peakData       <- parseArgs(args,"peakData=")                        # .xls result from MACS2
transcriptType <- parseArgs(args, "transcriptType=", "Bioconductor") # transcript annotation type
transcriptDb   <- parseArgs(args, "transcriptDb=", "")               # transcript annotation database
orgDb          <- parseArgs(args, "orgDb=", "")                      # genome wide annotation
regionTSS      <- parseArgs(args, "regionTSS=", 3000, "as.numeric")  # TSS region parameter 
ftargets       <- parseArgs(args, "targets=", "")                    # tab separated file describing the samples to be processed
orderby        <- parseArgs(args, "orderby=", "")                    # column in targets to order samples in plots
out            <- parseArgs(args,"out=", "Peak_Annotation") # output directory

runstr <- paste0("Called with: Rscript Peak_Annotation.R",
                   " peakData="      , peakData      ,
                   " transcriptType=", transcriptType,
                   " transcriptDb="  , transcriptDb  ,
                   " orgDb="         , orgDb         ,
                   " regionTSS="     , regionTSS     ,
                   " targets="       , ftargets      ,
                   " orderby="       , orderby       ,
                   " out="           , out           )
cat(runstr, fill=TRUE)
if (!is.numeric(regionTSS)) stop("regionTSS not numeric.")
if(ftargets != "" && !file.exists(ftargets)) stop("targets file", ftargets, "doesn't exist.")
targets <- read.delim(ftargets)
if(ftargets != "" && orderby != "" && !(orderby %in% colnames(targets))) stop("targets file doesn't have a column named", orderby)

if(length(list.files(peakData,pattern="_macs2_blacklist_filtered_peaks.xls")) > 0) {
      peakFiles <-list.files(peakData,pattern="_macs2_blacklist_filtered_peaks.xls", full.names = TRUE)
      filename <- strsplit(basename(peakFiles), "_macs2_blacklist_filtered_peaks.xls") # take the filenames and put it as names for the plots	
} else {
      peakFiles <-list.files(peakData,pattern="_macs2_peaks.xls", full.names = TRUE)
      filename <- strsplit(basename(peakFiles), "_macs2_peaks.xls") # take the filenames and put it as names for the plots     
}

# remove targets which have no peaks
peakcount <- sapply(peakFiles, function(x) {
  tryCatch({
    nrow(read.delim(x, head=TRUE, comment="#"))
  }, error=function(e) 0)
})
if(!all(peakcount > 0)) {
  warning("Sample(s) ", paste(basename(peakFiles)[!(peakcount > 0)], collapse=", "),
          " excluded because didn't have any peaks called")
  peakFiles <- peakFiles[peakcount > 0 ] 
  filename <- filename[peakcount > 0 ]
}


peaks <- lapply(peakFiles, readPeakFile) # read all the xls files using 'readPeakFile' function
# bug in ChIPseeker: MACS xls files (1-based) are read as 0-based. Modify condition if fixed in future version:
if(packageVersion('ChIPseeker')>0) { # bug in ChIPseeker: MACS xls files (1-based) are read as 0-based. Modify condition if fixed in future version.
  peaks <- lapply(peaks, function(x) {
    BiocGenerics::start(x) <- BiocGenerics::start(x)-1
    return(x)})
} 


if(transcriptType!="Bioconductor"){ # check the input format for the transcript annotation
   txdb <- makeTxDbFromGFF(transcriptDb, format="gtf") # if the input format is gtf file, then this file will be used to create a TxDb object
} else {
   library(transcriptDb, character.only = TRUE) # if the input format is bioconductor, then the transcript annotation library will be used 
   txdb <- eval(parse(text=transcriptDb))
  
}

if(orgDb!=""){ # check if genome wide annotation should be used
  peakAnno <- lapply(peaks, annotatePeak, TxDb=txdb, tssRegion=c(-regionTSS,regionTSS), annoDb=orgDb)
} else {
  peakAnno <- lapply(peaks, annotatePeak, TxDb=txdb, tssRegion=c(-regionTSS,regionTSS))  
}

names(peakAnno) <- filename

# sort samples if requested
if(ftargets != "" && orderby != "") {
  # built macs output name `NAMEoutput="\${IPname}.vs.\${INPUTname}"` as in NGSpipe2go/modules/ChIPseq/macs2.groovy
  targets$NAMEoutput <- paste0(targets$IPname, ".vs.", targets$INPUTname)
  targets$NAMEoutput <- sub("\\.vs\\.none", "", targets$NAMEoutput)   # fix macs2 outputs with no input

  # calculate order, first by column `orderby`, then by the original contrast name
  i <- order(targets[, orderby][match(names(peakAnno), targets$NAMEoutput)], names(peakAnno))
  peakAnno <- peakAnno[i]
}

# create barplot showing the feature distribution
png(file=paste0(out, "/ChIPseq_Feature_Distribution_Barplot.png"), width = 700, height = 500)
plot(plotAnnoBar(peakAnno))
dev.off()

# create barplot showing the feature distribution related to TSS
png(file=paste0(out, "/ChIPseq_Feature_Distribution_Related_to_TSS_Barplot.png"), width = 700, height = 500)
plot(plotDistToTSS(peakAnno))
dev.off()

# create upsetplot 
for(i in 1:length(peakAnno)){
  png(file=paste0(out, "/", names(peakAnno)[i], "_ChIPseq_UpSetplot.png"), width = 700, height = 500)
  print(upsetplot(peakAnno[[i]]))
  dev.off()
}

# create ChIP peaks coverage plot
for(i in 1:length(peakAnno)){
  png(file=paste0(out, "/", names(peakAnno)[i], "_ChIPseq_Peaks_Coverageplot.png"), width = 700, height = 500)
  plot(covplot(as.GRanges(peakAnno[[i]]),weightCol="X.log10.pvalue."))
  dev.off()
}

# create xls output contains the peak annotation
outputData <- lapply(peakAnno, as.data.frame)
names(outputData) <- make.names(substr(names(outputData), 1, 30), unique=TRUE)
write.xlsx(outputData, file=paste0(out, "/Peak_Annotation.xlsx"))

# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/ChIPseq_Peak_Annotation_session_info.txt", sep=""))
save(peakFiles,peaks,peakAnno,filename,outputData, file=paste0(out,"/Peak_Annotation.RData"))
