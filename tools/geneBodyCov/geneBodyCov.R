########################################
##
##  Reimplementation of geneBodyCoverage
## --
## Who:  Sergi Sayols
## When: 25-aug-2016
## --
## Input:
##   <bam=x.bam>
##   <gtf=genes.gtf>
##   <paired=yes|no>
##   <stranded=yes|no|reverse>
##   <multimappers=yes|no>
##   <outdir=./>
##   <threads=1>
## --
## Todo:
##   -Multiple input bam files
##   -Extend some bp to the 5' and 3' ends
##
########################################
options(stringsAsFactors=F)
library(ShortRead)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(Cairo)

##
## Parse input parms
##
parseArgs <- function(args, string, default=NULL, convert="as.character") {

	if(length(i <- grep(string, args, fixed=T)) == 1) 
		return(do.call(convert, list(gsub(string, "", args[i]))))
    
	if(!is.null(default)) default else do.call(convert, list(NA))
}

args <- commandArgs(trailingOnly=T)
BAM      <- parseArgs(args, "bam=", "")	# the bam file
GENESGTF <- parseArgs(args, "gtf=", "")   # the gtf file
PAIRED   <- parseArgs(args, "paired=", "no") # is the experiment strand specific?
STRANDED <- parseArgs(args, "stranded=", "no") # is the experiment strand specific?
MMAPPERS <- parseArgs(args, "multimappers=", "no") # include multimapping reads?
OUTDIR   <- parseArgs(args, "outdir=" , "./") # output directory
THREADS  <- parseArgs(args, "threads=", 1, "as.numeric") # number of threads to be used

print(args)
if(length(args) == 0 | args[1] == "-h" | args[1] == "--help")
	stop("Rscript geneBodyCov.R <bam=x.bam> <gtf=genes.gtf> <stranded=no> <outdir=./> <threads=1>")
if(!file.exists(BAM)) stop(paste("File", BAM, "does NOT exist"))
if(!file.exists(GENESGTF)) stop(paste("File", GENESGTF, "does NOT exist"))
if(is.na(PAIRED)   | !(grepl("no|yes", PAIRED))) stop("Paired has to be no|yes")
if(is.na(STRANDED) | !(grepl("no|yes|reverse", STRANDED))) stop("Stranded has to be no|yes|reverse")
if(is.na(MMAPPERS) | !(grepl("no|yes", MMAPPERS))) stop("Multimappers has to be no|yes")
if(is.na(THREADS))  stop("Threads has to be a number")

##
## read input bam file and calculate the coverage
##
cvg <-
    if(PAIRED == "no") {
        aln <- if(MMAPPERS == "yes") readGAlignments(BAM) else readGAlignments(BAM, param=ScanBamParam(tagFilter=list("NH"=1)))
        coverage(switch(STRANDED,
                        no=unstrand(aln),
                        yes=aln,
                        reverse=invertStrand(aln)))
    } else {
        aln <- if(MMAPPERS == "yes") readGAlignmentPairs(BAM) else readGAlignmentPairs(BAM, param=ScanBamParam(tagFilter=list("NH"=1)))
        strandMode(aln) <- switch(STRANDED,    # only works in GAlignmentPairs-class objects
                                  no=0,
                                  yes=1,
                                  reverse=2)
        coverage(aln)
    }
rm(aln); gc()   # free memory

##
## Read and flatten the gtf file
##
gtf <- import.gff(GENESGTF, format="gtf", feature.type="exon")
gtf <- reduce(split(gtf, elementMetadata(gtf)$gene_id))
gtf <- gtf[sapply(gtf, function(x) sum(width(x))) > 100]  # kick out genes shorter than 100bp
gtf <- keepSeqlevels(gtf, intersect(seqlevels(gtf), seqlevels(cvg)), pruning.mode="coarse")   # drop genes in chromosomes which we don't have coverage

##
## subset from the coverage only the gene regions, and calculate the binned coverage
##
rangeCov <- mclapply(gtf, function(gene){
  tryCatch({
    # get the absolute covarage for the current gene
    # strandedness of library is already taken into account in the coverage analysis (cvg <- ...)

    # minus strand genes need to be flipped
    if(unique(decode(strand(gene))) == "-") {
      x <- rev(unlist(cvg[gene], use.names=FALSE))
    } else {
      x <- unlist(cvg[gene], use.names=FALSE) 
    }
    
    # split the gene in 100 bins and calculate the avg coverage per bin
    bins <- cut(1:length(x), 100)
    x <- tapply(decode(x), bins, mean)

    # calculate the percentage of coverage per position
    if(max(x) > 0) x / max(x) else x
  }, error=function(e) numeric(100))  # gene outside the boundaries of the chromosome, returning a vector of 100 zeros
}, mc.cores=THREADS)

##
## calculate the per bin average across all genes
##
filename <- (paste0(OUTDIR, "/", gsub(".bam$", "_geneBodyCov.png", basename(BAM))))
png(filename, type="cairo")
try({
  rangeCov <- do.call(rbind, rangeCov)  # flatten the list
  rangeCov <- rangeCov[apply(rangeCov, 1, function(x) any(x > 0)), ]  # suppress not expressed genes
  avg <- apply(rangeCov, 2, mean) # and calculate the average per bin
  avg <- avg / max(avg) # which is then normalized again, as it seems to be in geneBodyCoverage.py from RSeQC
  names(avg) <- 1:100
  write.csv(avg, gsub("png", "csv", filename))
  # and plot
  plot(1:100, avg, type="l", ylim=c(0, 1), main=basename(BAM), xlab="Gene body percentile 5'->3'", ylab="Average normalized coverage")
  lines(lowess(1:100, avg, f=1/4), col="red", lwd=2)
})
dev.off()
