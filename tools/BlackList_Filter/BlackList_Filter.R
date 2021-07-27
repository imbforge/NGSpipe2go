## What: BlackList_Filter.R
## Who: Giuseppe Petrosino
## When: 2018-02-08
##
## Script to filter out peaks overlapping blacklisted genomic regions
##
## Args: 
## -----
## peakData=                       # path to the xls file result from MACS2
############################################################################
options(stringsAsFactors=FALSE)
library(ChIPseeker)
library(rtracklayer)


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
blacklistRegions <- parseArgs(args, "blacklistRegions=")
out <- parseArgs(args,"out=", "macs2")

runstr <- paste0("Call with: Rscript BlackList_Filter.R [peakData=",peakData,"] [blacklistRegions=",blacklistRegions,"] [out=",out,"]")
cat(runstr)


if(file.exists(blacklistRegions)) {


peakFiles <-list.files(peakData,pattern=".xls", full.names = TRUE) # list of the full path of the .xls file 
peaks <- lapply(peakFiles, readPeakFile) # read all the xls files using 'readPeakFile' function
blacklist <- import(blacklistRegions)

# remove peaks overlapping blacklist regions
peaks.wo.blacklst <- lapply(peaks, function(x) {
	m <- x[!x %over% blacklist]
})


# create a summary table
filename <- strsplit(basename(peakFiles), "_macs2_peaks.xls") # take the filenames 
peaks.number <- sapply(peaks, function(x) { mm <-length(x) })
peaks.wo.blacklst.number <- sapply(peaks.wo.blacklst, function(x) { mm <-length(x) })
peaks.in.blacklst.percentage <- 100-(100/(peaks.number)*peaks.wo.blacklst.number)
summary.table <- cbind(filename,peaks.number,peaks.wo.blacklst.number,peaks.in.blacklst.percentage)
names(peaks.wo.blacklst) <- paste0(filename, "_macs2_blacklist_filtered_peaks")


# create xls output contains the peaks wo blacklist regions
outputData <- lapply(peaks.wo.blacklst, as.data.frame)
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
  
  # remove peaks overlapping blacklist regions
  bed.wo.blacklst <- lapply(beds, function(x) {
    m <- x[!x %over% blacklist]
  })
  names(bed.wo.blacklst) <- gsub(i, paste0("_blacklist_filtered",i), basename(bedFiles))
  
  bed.df <- lapply(bed.wo.blacklst, as.data.frame, stringsAsFactors=F)
  sapply(names(bed.wo.blacklst), function (x) export.bed(bed.wo.blacklst[[x]], con=paste0(out, "/", x)))
  
  bedOutputData[[paste0("blacklist_filtered",i)]] <- bed.df
}


# save the sessionInformation
writeLines(capture.output(sessionInfo()), paste(out, "/ChIPseq_BlackList_Filter_session_info.txt", sep=""))
save(peakFiles,peaks,blacklist,filename,outputData, bedOutputData, file=paste0(out,"/BlackList_Filter.RData"))

} else {
  cat("\nBlacklist filtering skipped because no valid blacklist provided\n")
  save(blacklistRegions, file=paste0(out,"/BlackList_Filter.RData"))
}

