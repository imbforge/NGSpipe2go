###############################################
##
## dupRadar duplication rate analysis in RNAseq experiments
## 
## Who:  Sergi Sayols
## When: 2-Oct-2014
##
###############################################
library(dupRadar)
library(Cairo)

##
## Parse input parms
##
parseArgs <- function(args, string, default=NULL, convert="as.character") {

	if(length(i <- grep(string, args, fixed=T)) == 1) 
		return(do.call(convert, list(gsub(string, "", args[i]))))
    
	if(!is.null(default)) default else do.call(convert, list(NA))
}

##
## get name patterns from command line
##
args     <- commandArgs(T)
bam      <- parseArgs(args, "bam=")	# the bam file to analyse
gtf      <- parseArgs(args, "gtf=")	# usually, same GTF file as used in htseq-count
stranded <- parseArgs(args, "stranded=", "no")# no|yes|reverse
paired   <- parseArgs(args, "paired=", "no")	# is a paired end experiment"
outdir   <- parseArgs(args, "outdir=", "./")  # output directory
threads  <- parseArgs(args, "threads=", 1, "as.integer") # number of threads to be used

runstr <- "Call with: Rscript dupRadar.R bam=<file.bam> gtf=<genes.gtf> stranded=[no|yes|reverse] paired=[no|yes] outdir=./ threads=1\n"
if(length(args) == 0) { cat(runstr); quit(save="no") }
if(any(grepl("^-h|^--help", args))) { cat(runstr); quit(save="no") }
if(is.na(bam)) { cat(runstr); quit(save="no") }
if(is.na(gtf)) { cat(runstr); quit(save="no") }
if(is.na(stranded) | !(grepl("no|yes|reverse", stranded))) stop("Stranded has to be no|yes|reverse")
if(is.na(paired)   | !(grepl("no|yes", paired))) stop("Paired has to be no|yes")
if(is.na(threads)) stop("Threads has to be an integer number")
if(!file.exists(bam)) stop(paste("File", bam, "does NOT exist"))
if(!file.exists(gtf)) stop(paste("File", gtf, "does NOT exist"))
if(!file.exists(outdir)) stop(paste("Dir", outdir, "does NOT exist"))

##
## analyze duprates and plot
##
stranded <- if(stranded == "no") 0 else if(stranded == "yes") 1 else 2

dm <- analyzeDuprates(bam, gtf, stranded, (paired == "yes"), threads, autosort=FALSE)
CairoPNG(paste0(outdir, "/", gsub("\\.bam", "", basename(bam)), "_dupRadar.png"))
duprateExpDensPlot(DupMat=dm)
write.table(dm, sep="\t", quote=FALSE, row.names=FALSE,
            file=paste0(outdir, "/", gsub("(\\.dupmarked|)\\.bam", "", basename(bam)), "_dupRadar.tsv")) # also remove the .duprm prefix
dev.off()
