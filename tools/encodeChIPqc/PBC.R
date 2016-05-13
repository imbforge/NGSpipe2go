###################################
##
## PCR Bottleneck Coefficient (PBC):
## A measure of library complexity, i.e. how skewed the distribution of read counts per location
## is towards 1 read per location.
## 
## Who : Sergi Sayols
## When: 26 Sept 2014
##
## Defined in the ENCODE guidelines (https://genome.ucsc.edu/ENCODE/qualityMetrics.html) as:
##
## PBC = N1/Nd
## 
## (where N1= number of genomic locations to which EXACTLY  one unique mapping read maps, and
##        Nd= number of genomic locations to which AT LEAST one unique mapping read maps, i.e.
##  the number of non-redundant, unique mapping reads).
## 
## PBC is further described on the ENCODE Software Tools page. Provisionally, 0-0.5 is severe
## bottlenecking, 0.5-0.8 is moderate bottlenecking, 0.8-0.9 is mild bottlenecking, while 0.9-1.0
## is no bottlenecking. Very low values can indicate a technical problem, such as PCR bias, or a
## biological finding, such as a very rare genomic feature. Nuclease-based assays (DNase, MNase)
## detecting features with base-pair resolution (transcription factor footprints, positioned 
## nucleosomes) are expected to recover the same read multiple times, resulting in a lower PBC
## score for these assays. Note that the most complex library, random DNA, would approach 1.0, 
## thus the very highest values can indicate technical problems with libraries. It is the practice
## for some labs outside of ENCODE to remove redundant reads; after this has been done, the value
## for this metric is 1.0, and this metric is not meaningful. 82% of TF ChIP, 89% of His ChIP, 77%
## of DNase, 98% of FAIRE, and 97% of control ENCODE datasets have no or mild bottlenecking.
library(GenomicAlignments)

##
## Get IP and input from command line
##
args    <- commandArgs(T)
IP      <- args[1]

if(length(args) != 1)   stop("Rscript PBC.R <bam file>")
if(!file.exists(IP))    stop(paste("File",IP,"does NOT exist"))

cat("Program called with args:",args,fill=T)

##
## 1-Read bam files and count reads per position
##
system.time( {
	cat("Reading ",IP,"...\n")
	aln <- as.data.frame(readGAlignments(IP))
	cat("Summarizing tags on genomic positions\n")
	x <- table(ifelse(aln$strand == '+',
					  paste(aln$seqnames,aln$start),
					  paste(aln$seqnames,aln$end)))
	PBC <- sum(x == 1) / length(x)
})

write.csv(data.frame(IP,PBC),file=gsub("\\.bam$","_PBC.csv",IP),row.names=F)
