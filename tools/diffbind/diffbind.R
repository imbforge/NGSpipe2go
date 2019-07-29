#####################################
##
## What: diffbind.R
## Who : Sergi Sayols
## When: 10-05-2016
##
## Script to perform differential binding analysis between 2 conditions (pairwise)
##
## Args:
## -----
## targets=targets.txt          # file describing the targets.
## contrasts=contrasts.txt      # file describing the contrasts
## cwd=./                       # current working directory where the files .tsv files are located
## bams=paste0(CWD, "/mapped")  # directory with the bam files
## out=paste0(CWD, "/results")  # prefix filename for output
## fragsize=200                 # average fragment size
## substractControl=TRUE        # substract input
## fullLibrarySize=TRUE         # use total number of reads in bam for normalization (FALSE=only peaks)
## tagwiseDispersion=TRUE       # calculate dispersion tagwise (use FALSE if no replicates)
## annotate=TRUE                # annotate after DB analysis?
## pe=FALSE                     # paired end experiment?
## tss=c(-3000,3000)            # region around the tss
## txdb=TxDb.Mmusculus.UCSC.mm9.knownGene   # Bioconductor transcript database, for annotation
## annodb=org.Mm.eg.db          # Bioconductor gene annotation database
##
######################################
options(stringsAsFactors=FALSE)
library(DiffBind)
library(openxlsx)

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
FTARGETS   <- parseArgs(args,"targets=","targets.txt")     # file describing the targets
FCONTRASTS <- parseArgs(args,"contrasts=","contrasts.txt") # file describing the contrasts
CWD        <- parseArgs(args,"cwd=","./")     # current working directory
BAMS       <- parseArgs(args,"bams=",paste0(CWD, "/mapped"))  # directory with the bam files
OUT        <- parseArgs(args,"out=", paste0(CWD, "/results")) # directory where the output files will go
FRAGSIZE   <- parseArgs(args,"fragsize=", 200, "as.numeric")# fragment size
SUBSTRACTCONTROL  <- parseArgs(args,"substractControl=", "TRUE", "as.logical")  # substract input
FULLLIBRARYSIZE   <- parseArgs(args,"fullLibrarySize=", "TRUE", "as.logical")   # use total number of reads in bam for normalization (FALSE=only peaks)
TAGWISEDISPERSION <- parseArgs(args,"tagwiseDispersion=", "TRUE", "as.logical") # calculate dispersion tagwise (use FALSE if no replicates)
ANNOTATE   <- parseArgs(args,"annotate=", TRUE, "as.logical") # annotate after DB analysis?
PE         <- parseArgs(args,"pe=", FALSE, "as.logical")      # paired end experiment?
TSS        <- parseArgs(args,"tss=", c(-3000,3000), "run_custom_code") # region around the tss
TXDB       <- parseArgs(args,"txdb=", "TxDb.Mmusculus.UCSC.mm9.knownGene") # Bioconductor transcript database, for annotation 
ANNODB     <- parseArgs(args,"annodb=", "org.Mm.eg.db") # Bioconductor gene annotation database

runstr <- "Rscript diffbind.R [targets=targets.txt] [contrasts=contrasts.txt] [cwd=./] [bams=./mapped] [out=./results] [fragsize=200] [annotate=TRUE] [pe=TRUE] [tss=c(-3000,3000)] [txdb=TxDb.Mmusculus.UCSC.mm9.knownGene] [annodb=org.Mm.eg.db]"
if(!file.exists(CWD))        stop("Dir",CWD,"does NOT exist. Run with:\n",runstr)
setwd(CWD)
if(!file.exists(FTARGETS))   stop("File",FTARGETS,"does NOT exist. Run with:\n",runstr)
if(!file.exists(FCONTRASTS)) stop("File",FCONTRASTS,"does NOT exist. Run with:\n",runstr)
if(!file.exists(BAMS))       stop("Dir",BAMS,"does NOT exist. Run with:\n",runstr)
if(!is.numeric(FRAGSIZE))    stop("Fragment size not numeric. Run with:\n",runstr)
if(!is.logical(ANNOTATE))    stop("Annotate not logical. Run with:\n",runstr)
if(!is.logical(PE))          stop("Paired end (pe) not logical. Run with:\n",runstr)
if(ANNOTATE & !is.numeric(TSS)) stop("Region around TSS not numeric. Run with:\n",runstr)
if(ANNOTATE & !require(TXDB, character.only=TRUE))   stop("Transcript DB", TXDB, "not installed\n")
if(ANNOTATE & !require(ANNODB, character.only=TRUE)) stop("Annotation DB", ANNODB, "not installed\n")

pdf(paste0(OUT, "/diffbind.pdf"))

##
## make DB analysis
##
# load targets and make analysis
conts   <- read.delim(FCONTRASTS, head=F, comment.char="#")
targets <- read.delim(FTARGETS, head=T, colClasses="character", comment.char="#", )
db <- dba(sampleSheet=targets, config=data.frame(fragmentSize=FRAGSIZE, bCorPlot=F, singleEnd=!PE))
db <- dba.count(db, bUseSummarizeOverlaps=PE)  # bUseSummarizeOverlaps method slower and memory hungry, mandatory only for PE data
dba.plotPCA(db, DBA_CONDITION, label=DBA_CONDITION)

# parse the formula in cont and do the analysis
result <- lapply(conts[, 1], function(cont) {
  cat(cont, fill=T)
  cont.name <- gsub("(.+)=(.+)", "\\1", cont)
  cont.form <- gsub("(.+)=(.+)", "\\2", cont)
  factors   <- unlist(strsplit(cont.form, "\\W"))
  factors   <- factors[factors != ""]

  c1 <- dba.mask(db, DBA_CONDITION, factors[1])
  c2 <- dba.mask(db, DBA_CONDITION, factors[2])
  db <- dba.contrast(db, group1=c1, group2=c2,  name1=factors[1], name2=factors[2], categories=DBA_CONDITION)
  db <- dba.analyze(db, bSubControl=SUBSTRACTCONTROL, bFullLibrarySize=FULLLIBRARYSIZE, bTagwise=TAGWISEDISPERSION)
  dba.plotMA(db)
  dba.plotBox(db)
  dba.plotVenn(db, c1 | c2)

  dba.report(db, bCalled=T)
})

##
## Annotate peaks
##
if(ANNOTATE) {
  library(ChIPseeker)
  txdb <- eval(parse(text=TXDB))
  result <- lapply(result, annotatePeak, TxDb=txdb, annoDb=ANNODB, tssRegion=TSS, verbose=T)
  lapply(result, plotAnnoBar)
  lapply(result, plotDistToTSS)
}

dev.off()

writeLines(capture.output(sessionInfo()),paste(OUT, "/diffbind_session_info.txt", sep=""))
result <- lapply(result, as.data.frame)
names(result) <- substr(gsub("(.+)=\\((.+)\\)", "\\2", conts[,1]), 1, 31)
write.xlsx(result, file=paste0(OUT, "/diffbind.xlsx"))
saveRDS(result,  file=paste0(OUT, "/diffbind.rds"))

