####################################
##
## Plot IPstrength in the way CHANCE does
## Ref: Diaz, A et al. Normalization, bias correction, and peak calling for ChIP-seq.
##      Stat Appl Genet Mol Biol. 2012 March 31
## --
## Who:  Sergi Sayols
## When: 13-jan-2014
## --
## The process follows these steps:
##   1-Read bam files and count reads per bin
##   2-Sort bins in IP, and sort input in the same order as IP
##   3-Calculate for IP and input the cumulative sum of reads
##   4-Calculate for IP and input the cumulative percentages
##   5-Calculate the alpha scaling factor (max diff between he cumulative percentages)
##   6-Plot IPstrength
## --
## Input:
##   arg1: IP BAM file
##   arg2: IP sample name
##   arg3: Input BAM file
##   arg4: Input sample name
##   arg5: output file name
## --
## REMEMBER: Change at chunk1 the reference organism!! set to human (hg19) by default
##
####################################
library(ShortRead)
library(Rsamtools)
library(parallel)

##
## Static parms
##
BIN_SIZE = 1000			# in human genome, ~1500000 bins
STRAND_SPECIFIC = "*"	# for strand-blind sample prep protocol

##
## Get IP and input from command line
##
args  <- commandArgs(T)
IP         <- args[1]
IP.name    <- args[2]
input      <- args[3]
input.name <- args[4]
output     <- args[5]
bsgenome   <- args[6]

if(length(args) != 6)   stop("Rscript IPstrength.R <IP> <IP.name> <input> <input.name> <output> <bsgenome>")
if(!file.exists(IP))    stop(paste("File",IP   ,"does NOT exist"))
if(!file.exists(input)) stop(paste("File",input,"does NOT exist"))
if(!any(grepl(bsgenome,list.files(.libPaths())))) warning(paste(bsgenome,"not available in",.libPaths()))
org <- unlist(strsplit(bsgenome,"\\."))[2]

cat("Program called with args:",args,fill=T)

##
## 1-Read bam files and count reads per bin
##
if(!require(bsgenome,character.only=T)) {
    cat("Tiling genome from",IP,"...\n")
    aln  <- BamFile(IP)
    bins <- tileGenome(seqinfo(aln),tilewidth=BIN_SIZE,cut.last.tile.in.chrom=TRUE)
} else {
    bins <- tileGenome(seqinfo(get(org)),tilewidth=BIN_SIZE,cut.last.tile.in.chrom=TRUE)
}

counter <- function(fl,bins)
{
	cat("Reading",fl,"...\n")
	aln <- readGAlignments(fl)
	strand(aln) <- STRAND_SPECIFIC
	hits   <- countOverlaps(aln,bins)
	counts <- countOverlaps(bins,aln[hits==1])	# discard reads hitting more than one gene
	names(counts) <- names(bins)
	counts
}
counts <- mclapply(c(IP,input),counter,bins,mc.cores=2)
names(counts) <- c("IP","input")

##
##   2-Sort bins in IP, and sort input in the same order as IP
##
o <- order(counts[[1]])
counts[[2]] <- counts[[2]][o]	# sort input
counts[[1]] <- counts[[1]][o]	# sort IP

##
##   3-Calculate for IP and input the cumulative sum of reads
##
cs <- mclapply(counts,cumsum,mc.cores=2)

##
##   4-Calculate for IP and input the cumulative percentages
##
csp <- mclapply(cs,function(x) x / x[length(x)],mc.cores=2)

##
##   5-Calculate the alpha scaling factor (max diff between he cumulative percentages)
##
k <- which.max(abs(csp[[2]] - csp[[1]]))
alpha      <- cs[[1]][k] / cs[[2]][k]	# scaling factor: ratio between the cumulative sums at k
enrichment <- 1 - (k / length(bins))

##
##   6-Plot IPstrength
##
x <- 1:length(bins)
x <- x / x[length(x)]

#df <- data.frame(x=1:length(bins),y=csp[[1]],class="Input")
#df <- rbind(df,data.frame(x=1:length(bins),y=csp[[2]],class="IP"))
#p <- ggplot(df,aes(x=x,y=y,color=class)) + geom_point()
#print(p)

png(paste0(output,".png"))
plot  (x=x,y=csp[[1]],col="blue",xlab="% of bins",ylab="% of tags",type='l',lwd=6,
	   main=paste0("alpha=",round(alpha,4),", enrichment=",round(enrichment*100,2),"%"))	# IP
lines(x=x,y=csp[[2]],col="red" ,lwd=6)	# input
abline(v=x[k],col="green",lty=2,lwd=2)			# location of the multiplicative scaling factor alpha
legend("topleft",legend=c(input.name,IP.name),fill=c("red","blue"))
dev.off()

cat(IP.name,input.name,alpha,enrichment,sep=",",file=paste0(output,".csv"),fill=T,append=T)
