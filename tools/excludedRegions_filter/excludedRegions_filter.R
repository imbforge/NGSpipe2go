## What: excludedRegions_filter.R
## Who: Giuseppe Petrosino
## When: 2018-02-08
##
## Script to filter out peaks overlapping genomic regions to be excluded
##
## Args: 
## -----
## peakData=                       # path to the xls file result from MACS2
############################################################################
options(stringsAsFactors=FALSE)
library(ChIPseeker)
library(rtracklayer)
library(tidyr)


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
excludedRegionsBED <- parseArgs(args, "excludedRegionsBED=")
out <- parseArgs(args,"out=", "macs2")

runstr <- paste0("Call with: Rscript excludedRegions_filter.R [peakData=",peakData,"] [excludedRegionsBED=",excludedRegionsBED,"] [out=",out,"]")
cat(runstr)


if(file.exists(excludedRegionsBED)) {


peakFiles <-list.files(peakData,pattern=".xls", full.names = TRUE) # list of the full path of the .xls file 

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
}

peaks <- lapply(peakFiles, readPeakFile) # read all the xls files using 'readPeakFile' function
if(packageVersion('ChIPseeker')>0) { # bug in ChIPseeker: MACS xls files (1-based) are read as 0-based. Modify condition if fixed in future version.
        peaks <- lapply(peaks, function(x) {
          BiocGenerics::start(x) <- BiocGenerics::start(x)-1
          return(x)})
        } 

exclRegions <- import(excludedRegionsBED)

# remove peaks overlapping regions to exclude
peaks.wo.exclReg <- lapply(peaks, function(x) {
	m <- x[!x %over% exclRegions]
})


# create a summary table
filename <- strsplit(basename(peakFiles), "_macs2_peaks.xls") # take the filenames 
peaks.number <- sapply(peaks, function(x) { mm <-length(x) })
peaks.wo.exclReg.number <- sapply(peaks.wo.exclReg, function(x) { mm <-length(x) })
peaks.in.exclReg.percentage <- 100-(100/(peaks.number)*peaks.wo.exclReg.number)
summary.table <- cbind(filename,peaks.number,peaks.wo.exclReg.number,peaks.in.exclReg.percentage)
names(peaks.wo.exclReg) <- paste0(filename, "_macs2_excludedRegions_filtered_peaks")


# create xls output containing the peaks wo excluded regions
columnNames2replace <- c(seqnames="chr", X.log10.pvalue.="-log10(pvalue)", X.log10.FDR="-log10(qvalue)", X.log10.qvalue.="-log10(qvalue)")
outputData <- lapply(peaks.wo.exclReg, function(x) {
  x <- as.data.frame(x)
  x$strand <- x$width <- NULL # the 'strand' column may cause error when read in again with ChIPseeker::readPeakFile
  colnames(x) <- dplyr::recode(colnames(x), !!!columnNames2replace) # rename column names for conformity with unfiltered peak file
  return(x)
  })

sapply(names(outputData), 
 function (x) write.table(outputData[[x]], file=paste0(out, "/", x, ".xls"), row.names=F, sep="\t"))

write.csv(summary.table, file=paste0(out, "/peaks_detected_table.csv"), row.names=F)


## filter bed files
bedfile_suffixes <- c(".narrowPeak", ".broadPeak")
bedOutputData <- list()

for (i in bedfile_suffixes) {
  
  bedFiles <-list.files(peakData,pattern=paste0(i, "$"), full.names = TRUE) # list of the full path of the bed file
  if(length(bedFiles)==0) {next}
  beds <- lapply(bedFiles, import) # read all the bed files using 'import' function
  
  # remove peaks overlapping regions to exclude
  bed.wo.exclReg <- lapply(beds, function(x) {
    m <- x[!x %over% exclRegions]
  })
  names(bed.wo.exclReg) <- gsub(i, paste0("_excludedRegions_filtered",i), basename(bedFiles))
  
  bed.df <- lapply(bed.wo.exclReg, as.data.frame, stringsAsFactors=F)
  sapply(names(bed.wo.exclReg), function (x) export.bed(bed.wo.exclReg[[x]], con=paste0(out, "/", x)))
  
  bedOutputData[[paste0("excludedRegions_filtered",i)]] <- bed.df
}


# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/ChIPseq_excludedRegions_filter_session_info.txt", sep=""))
save(peakFiles,peaks,exclRegions,filename,outputData, bedOutputData, file=paste0(out,"/excludedRegions_filter.RData"))

} else {
  cat("\nExcludedRegions filtering skipped because no valid BED file provided\n")
  save(excludedRegionsBED, file=paste0(out,"/excludedRegions_filter.RData"))
}

