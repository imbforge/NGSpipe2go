###################################
##
## Shortened version of the spp.R QC to report only the phantompeak profile
## Based on the SPP tutorial from http://compbio.med.harvard.edu/Supplements/ChIP-seq/tutorial.html
##
## Who : Sergi Sayols
## When: 11 nov 2014
##
## Types of analysis:
##   *Calculate binding characteristics
##   *Report QC metrics RSC and NSC
##
##################################

# Get cmd arguments
argv        <- commandArgs(trailingOnly=T)
SAMPLE_FILE <- argv[1]
SAMPLE      <- argv[2]
MIN.SHIFT   <- as.integer(argv[3])
MAX.SHIFT   <- as.integer(argv[4])
BIN.SIZE    <- as.integer(argv[5])
READ.LEN    <- as.integer(argv[6])
CORES       <- if(as.integer(argv[7]) > 16) 16 else as.integer(argv[7])
print(argv)

# load the library
library(spp)
#library(snow)
library(parallel)
cluster <- makeCluster(CORES)

##
## Loading tag data, selecting choosing alignment quality, removing anomalies
##
chip.data  <- read.bam.tags(SAMPLE_FILE)

# get binding info from cross-correlation profile
# srange gives the possible range for the size of the protected region;
# srange should be higher than tag length; making the upper boundary too high will increase calculation time
#
# bin - bin tags within the specified number of basepairs to speed up calculation;
# increasing bin size decreases the accuracy of the determined parameters
#
# IMPORTANT: add remove.tag.anomalies=F if duplicates were already removed
#
binding.characteristics <- get.binding.characteristics(chip.data,srange=c(MIN.SHIFT,MAX.SHIFT),bin=BIN.SIZE,remove.tag.anomalies=F,cluster=cluster)

stopCluster(cluster)

# select informative tags based on the binding characteristics
chip.data  <- select.informative.tags(chip.data,binding.characteristics)

# restrict or remove singular positions with very high tag counts
chip.data  <- remove.local.tag.anomalies(chip.data)

##
## Report QC metrics RSC and NSC
## The Normalized Strand Coeficient, NSC, is the normalized ratio between the fragment-length 
## cross-correlation peak and the background cross-correlation.
## The Relative Strand Correlation, RSC, is the ratio between the fragment-length peak and the 
## read-length peak.
## ENCODE cutoff: NSC values > 1.05 and RSC values > 0.8
##

## Detect the peak at the fragment length
cc    <- binding.characteristics$cross.correlation

# detect all the candidate peaks by comparing cc$y[i] to cc$y[i+/-bw]
bw    <- ceiling(2/BIN.SIZE)	# adjust if bin size was smaller than 2
slope <- as.numeric(cc$y [(1+bw):length(cc$y )] - cc$y [1:(length(cc$y) -bw)] < 0) # 0:pos, 1:neg
# 1: a peak (from left to right, shifting from pos to neg slope, -1: a valley shifting from neg to pos slope
shift <- as.numeric(slope[(1+bw):length(slope)] - slope[1:(length(slope)-bw)]) 
peaks <- which(shift == 1) + bw	# phantompeak tools busquen el peak amb -1!!

# Remove fake peaks from 10:READ.LEN+10 bp
# Correct the peak detected by get.binding.characteristics putting to one found discarding the area around READ.LEN
# If the highest peak is within the discarded area, it means a problematic IP, and the max peak doesnt correspond to the fragment size
peaks <- peaks[(cc$x[peaks] < 10) | cc$x[peaks] > READ.LEN+10]
binding.characteristics$peak$x <- cc$x[peaks[which.max(cc$y[peaks])]]
binding.characteristics$peak$y <- cc$y[peaks[which.max(cc$y[peaks])]]

## Detect the peak at the read length (phantom peak)
peaks <- which((cc$x >= (READ.LEN - round(2*BIN.SIZE))) &
			   (cc$x <= (READ.LEN + round(2*BIN.SIZE))))
binding.characteristics$phantompeak$x <- cc$x[peaks[which.max(cc$y[peaks])]]
binding.characteristics$phantompeak$y <- cc$y[peaks[which.max(cc$y[peaks])]]

## Detect the minim correlation within the window, which will be the background cross-correlation
binding.characteristics$back.cc$x <- cc$x[which.min(cc$y)]
binding.characteristics$back.cc$y <- cc$y[which.min(cc$y)]

## plot cross-correlation profile
NSC <- (binding.characteristics$peak$y / binding.characteristics$back.cc$y)
RSC <- (binding.characteristics$peak$y        - binding.characteristics$back.cc$y) /
	   (binding.characteristics$phantompeak$y - binding.characteristics$back.cc$y)

#pdf(file=paste0(SAMPLE,".crosscorrelation.pdf"),width=5,height=5)
png(paste0(SAMPLE,"_phantompeak.png"))
plot(binding.characteristics$cross.correlation$x,runmean(binding.characteristics$cross.correlation$y,BIN.SIZE),
	 type='l',xlab="strand shift",ylab="cross-correlation",
	 main=paste(SAMPLE,"\n","NSC",round(NSC,2),"RSC",round(RSC,2),
				"Read",binding.characteristics$phantompeak$x,"bp",
				"Fragment",binding.characteristics$peak$x,"bp"))
lines(x=c(binding.characteristics$phantompeak$x,binding.characteristics$phantompeak$x,min(cc$x)),
	  y=c(min(cc$y),binding.characteristics$phantompeak$y,binding.characteristics$phantompeak$y),
	  lty=2,col="green")
lines(x=c(binding.characteristics$peak$x,binding.characteristics$peak$x,min(cc$x)),
	  y=c(min(cc$y),binding.characteristics$peak$y,binding.characteristics$peak$y),
	  lty=2,col="blue")
lines(x=c(binding.characteristics$back.cc$x,binding.characteristics$back.cc$x,min(cc$x)),
	  y=c(min(cc$y),binding.characteristics$back.cc$y,binding.characteristics$back.cc$y),
	  lty=2,col="red")
legend("topright",legend=c("phantom peak at read length","peak at fragment length","background cross-correlation"),fill=c("green","blue","red"))
dev.off()
