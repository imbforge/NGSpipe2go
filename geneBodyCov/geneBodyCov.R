########################################
##
## Calcultate the gene body coverage:
## reimplementation of geneBodyCoverage from RSeQC 'cause it sucks for anything
## else than their precompiled gene models (human, mouse, fly, fish)
##
########################################
options(stringsAsFactors=F)
library(ShortRead)
library(GenomicRanges)
library(rtracklayer)
library(parallel)

##
## Parse input parms
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {

	if(length(i <- grep(string,args,fixed=T)) == 1) 
		return(do.call(convert,list(gsub(string,"",args[i]))))
    
	if(!is.null(default)) default else do.call(convert,list(NA))
}

args <- commandArgs(trailingOnly=T)
BAM      <- parseArgs(args,"bam=","")	# the bam file
GENESGTF <- parseArgs(args,"gtf=","")   # the gtf file
PAIRED   <- parseArgs(args,"paired=","no") # is the experiment strand specific?
STRANDED <- parseArgs(args,"stranded=","no") # is the experiment strand specific?
OUTDIR   <- parseArgs(args,"outdir=" ,"./") # output directory
THREADS  <- parseArgs(args,"threads=",1,"as.numeric") # number of threads to be used

print(args)
if(length(args) == 0 | args[1] == "-h" | args[1] == "--help")
	stop("Rscript geneBodyCov.R <bam=x.bam> <gtf=genes.gtf> <stranded=no> <outdir=./> <threads=1>")
if(!file.exists(BAM)) stop(paste("File",BAM,"does NOT exist"))
if(!file.exists(GENESGTF)) stop(paste("File",GENESGTF,"does NOT exist"))
if(is.na(PAIRED)   | !(grepl("no|yes",PAIRED))) stop("Paired has to be no|yes")
if(is.na(STRANDED) | !(grepl("no|yes|reverse",STRANDED))) stop("Stranded has to be no|yes|reverse")
if(is.na(THREADS))  stop("Threads has to be a number")

##
## Read and flatten the gtf file
##
# our own implementation of the tile function, to include also the gene name
# check the original code with showMethod("tile") and selectMethod("tile","GRanges")
xtile <- function (x, n, width, ...) { 
	sn <- seqnames(x) 
	strand <- strand(x) 
	x <- ranges(x)
	genes <- names(x)
	tiles <- IRanges::tile(x,n,width,...)
	gr <- GRanges(rep(sn, elementLengths(tiles)),		# seqnames
				  unlist(tiles),						# ranges
				  rep(strand,elementLengths(tiles)),	# strand
				  rep(genes,elementLengths(tiles)))		# gene names
#	relist(gr, tiles)	# do not provide as a GRangesList. GRanges more suitable for countOverlaps
}

# read and strip input gtf file
gtf <- import.gff(GENESGTF,format="gtf",asRangedData=F,feature.type="exon")
gtf <- unlist(reduce(split(gtf,elementMetadata(gtf)$gene_id)))
gtf <- unlist(xtile(gtf,width=1))

##
## read input bam file count overlaps between bam and gtf
##
aln <- switch(PAIRED,
			  no =readGAlignments    (BAM,param=ScanBamParam(tag="NH")),
			  yes=readGAlignmentPairs(BAM,param=ScanBamParam(tag="NH")))
# discard multihits
aln  <- switch(PAIRED,
			   no =aln[elementMetadata(aln)$NH == 1,],
			   yes=aln[elementMetadata(first(aln))$NH == 1,])	# both mates have the same NH
# change strand if reverse
if(STRANDED == "reverse") strand(aln) <- ifelse(strand(aln) == "+","-","+")
# and count reads on features
counts <- countOverlaps(gtf,aln,ignore.strand=(STRANDED == "no"))

# Normalize gene lengths to 0-100 bins
ncounts <- tapply(counts,elementMetadata(gtf)[,1],function(x) {

	# kick out gene if it was shorter than 100bp
	if(length(x) < 100) return(rep(0,100))

	# calculate the number of positions included in each bin
	w <- length(x) / 100	# has to be corrected for decimals
	r <- rep(floor(w),100)	# number of bp without correcting for decimals
	p <- sample(1:length(r),length(x) - sum(r),replace=F)
	r[p] <- floor(w) + 1	# add 1 bp to some selected bins

	# assign every bp to a bin and calculate the average
	xnorm <- tapply(x,rep(1:100,r),mean)
})

#ncounts <- matrix(unlist(ncounts),nrow=length(ncounts),ncol=100,byrow=T)
ncounts <- do.call(rbind,ncounts)

##
## plot the gene body coverage
##
# get rid of non expressed genes
ncounts <- ncounts[rowSums(ncounts) > 0,]

# make per position (per bin) average coverage across all genes
ncounts.avg <- apply(ncounts,2,mean)

# and plot
pdf(paste0(OUTDIR,"/",gsub(".bam","_geneBodyCov.pdf",basename(BAM))))
plot(1:100,ncounts.avg,type='l',yaxt="n",main=basename(BAM),xlab="gene length percentage",ylab="average coverage")
lines(lowess(1:100,ncounts.avg,f=1/4),col='red',lwd=2)
dev.off()
