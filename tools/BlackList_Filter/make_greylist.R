#####################################
##
## What: make_greylist.R
## Who : Frank Rühle
## When: 05-19-2022
##
## Script to create grey lists from control bam files to filter genomic regions
##
## Args:
## -----
## FTARGETS    # file describing the targets.
## BAMS        # directory with the bam files
## PEAKS       # directory with peak caller output
## OUT         # prefix filename for output
## KAR         # file containing chromosome sizes for the reference genome of interest, one per line, as "chromName chromLength" pairs
## MAXGAP      # If the distance between neighbouring grey regions is less than or equal to maxGap, the regions will be merged into one big region
##
######################################
options(stringsAsFactors=FALSE)
library(GreyListChIP) 
library(GenomicRanges)


##
## get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
  
  if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  
  if(!is.null(default)) default else do.call(convert,list(NA))
}
run_custom_code <- function(x) {
  eval(parse(text=x))
}

args <- commandArgs(T)
FTARGETS   <- parseArgs(args,"targets=","targets.txt")
PEAKS      <- parseArgs(args,"peaks=",paste0(CWD, "/results/macs2"))  # directory with the peak files
BAMS       <- parseArgs(args,"bams=",paste0(CWD, "/mapped"))  # directory with the bam files
OUT        <- parseArgs(args,"out=", paste0(CWD, "/results")) # directory where the output files will go
KAR        <- parseArgs(args,"kar=","genome.fa.fai")
MAXGAP     <- parseArgs(args,"maxgap=","10000", convert="as.numeric")


#### pre-calculate greylist ###
# see DiffBind vignette section 6 Blacklists and Greylists:
# https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
# and GreyListChIP vignette to prepare a greylist
# https://bioconductor.org/packages/release/bioc/vignettes/GreyListChIP/inst/doc/GreyList-demo.pdf

  # GreyListChIP vignette: The greylist object is a list with two elements. The first, greylist$master,
  # contains the full greylist to be applied. The other element, greylist$controls,
  # is a GRangesList containing the individual greylists computed for each of the control tracks  

  
  targets <- read.delim(FTARGETS, head=T, colClasses="character", comment.char="#")
  # determine file suffixes for targets 
  donefiles <- list.files(PEAKS,pattern=".done$")
  bam_suffix <- sub("^[^\\.]*\\.*", "", gsub("_macs2.done$", "", donefiles[1]))
  bam_suffix <- ifelse(bam_suffix == "", paste0(bam_suffix, "bam"), paste0(bam_suffix, ".bam"))
  targets$bamControl= paste0(BAMS, "/", targets$INPUT, ".", bam_suffix)

  control_files <- unique(targets$bamControl)
  gl_controls <- GRangesList()
  for(i in control_files) { 
    cname <- gsub("\\..*$", "", basename(i))
    print(cname)
    gl <- new("GreyList",karyoFile=KAR) # generating a tiling of the genome
    gl <- countReads(gl,bamFile=i) # counting reads from a BAM file for the tiling
    gl <- calcThreshold(gl,reps=100,sampleSize=30000,p=0.99,cores=1) # sampling from the counts and fitting the samples to the negative binomial distribution to calculate the read count threshold
    gl <- makeGreyList(gl,maxGap=MAXGAP) # filtering the tiling to identify regions of high signal and create the greylist
    gl_controls[[cname]] <- gl@regions # extract the region GRanges from the greylist object
  }
  
  gl_master <- Reduce(GenomicRanges::union, gl_controls) # the master greylist is just the combination of control greylists
  export(gl_master,con=file.path(OUT, "masterGreyList.bed"))
  greylistMnC <- list(master=gl_master, controls=gl_controls)
  
  saveRDS(greylistMnC,  file=paste0(OUT, "/greylist.rds"))
