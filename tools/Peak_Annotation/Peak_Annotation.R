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
## out=                                             # output directory
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

args <- commandArgs(T)
peakData <- parseArgs(args,"peakData=") # .xls result from MACS2
transcriptType <- parseArgs(args, "transcriptType=", "Bioconductor") # transcript annotation type
transcriptDb <- parseArgs(args, "transcriptDb=", "TxDb.Hsapiens.UCSC.hg19.knownGene") # transcript annotation database
orgDb <- parseArgs(args, "orgDb=", "org.Hs.eg.db") # genome wide annotation
regionTSS <- parseArgs(args, "regionTSS=", 3000, "as.numeric") # TSS region parameter 
out <- parseArgs(args,"out=", "Peak_Annotation") # output directory

runstr <- paste0("Call with: Rscript Peak_Annotation.R [peakData=",peakData,"] [transcriptType=",transcriptType,"] [transcriptDb=",transcriptDb,"] [orgDb=",orgDb,"] [regionTSS=",regionTSS,"][out=",out,"]")
cat(runstr)
if (!is.numeric(regionTSS)) stop("regionTSS not numeric. Run with:\n",runstr)

peakFiles <-list.files(peakData,pattern=".xls", full.names = TRUE) # list of the full path of the .xls file 

if(length(list.files(peakData,pattern="_macs2_blacklist_filtered_peaks.xls")) > 0) {
      peakFiles <-list.files(peakData,pattern="_macs2_blacklist_filtered_peaks.xls", full.names = TRUE)
      filename <- strsplit(basename(peakFiles), "_macs2_blacklist_filtered_peaks.xls") # take the filenames and put it as names for the plots	
} else {
      peakFiles <-list.files(peakData,pattern="_macs2_peaks.xls", full.names = TRUE)
      filename <- strsplit(basename(peakFiles), "_macs2_peaks.xls") # take the filenames and put it as names for the plots     
}

peaks <- lapply(peakFiles, readPeakFile) # read all the xls files using 'readPeakFile' function

if(transcriptType!="Bioconductor"){ # check the input format for the transcript annotation
   txdb <- makeTxDbFromGFF(transcriptDb, format="gtf") # if the input format is gtf file, then this file will be used to create a TxDb object
} else {
   library(transcriptDb, character.only = TRUE) # if the input format is bioconductor, then the transcript annoation library will be used 
   txdb <- eval(parse(text=transcriptDb))
  
}

if(orgDb!=""){ # check if genome wide annotation should be used
  peakAnno <- lapply(peakFiles, annotatePeak, TxDb=txdb, tssRegion=c(-regionTSS,regionTSS), annoDb= orgDb)
} else {
  peakAnno <- lapply(peakFiles, annotatePeak, TxDb=txdb, tssRegion=c(-regionTSS,regionTSS))  
}

names(peakAnno) <- filename

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
  png(file=paste0(out, "/", filename[[i]], "_ChIPseq_UpSetplot.png"), width = 700, height = 500)
  print(upsetplot(peakAnno[[i]]))
  dev.off()
}

# create ChIP peaks coverage plot
for(i in 1:length(peakAnno)){
  png(file=paste0(out, "/", filename[[i]], "_ChIPseq_Peaks_Coverageplot.png"), width = 700, height = 500)
  plot(covplot(peaks[[i]],weightCol="X.log10.pvalue."))
  dev.off()
}

# create xls output contains the peak annotation
outputData <- lapply(peakAnno, as.data.frame)
write.xlsx(outputData, file=paste0(out, "/Peak_Annotation.xlsx"))

# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/ChIPseq_Peak_Annotation_session_info.txt", sep=""))
save(peakFiles,peaks,peakAnno,filename,outputData, file=paste0(out,"/Peak_Annotation.RData"))
