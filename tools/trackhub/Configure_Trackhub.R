############################################################################
##
## What: Configure_Trackhub.R
## Who: Martin Oti
## When: 2018-07-17
##
## Script to create UCSC Track Hub configuration file
## with tracks from specific project.
## Configuration file is in YAML format with two sections:
## - Global: General configuration data (genome assembly, file locations, etc)
## - Tracks: Per-track configuration data (for 'trackDb.txt' file)
##
## Args: 
## -----
## ESSENTIAL_PROJECT=       # main project folder
## ESSENTIAL_DB=            # UCSC assembly (e.g. mm9, hg19)
## ESSENTIAL_CHROMSIZES=    # chromosome sizes file
## TRACKHUB_TRACKSDIR=      # project tracks folder
## TRACKHUB_CONFIG=         # trackhub configuration file
## TRACKHUB_TARGETS=        # file with sample info
## TRACKHUB_FTPBASE=        # Public FTP root folder
## TRACKHUB_FTPURLBASE=     # Public FTP root URL
## TRACKHUB_UCSCCFG=        # UCSC configuration files folder (with narrowPeak.as/broadPeak.as)
## ESSENTIAL_STRANDED=      # set to "yes" for stranded BigWigs
## TRACKHUB_PEAKSDIR=       # ChIP-seq peaks dir
##
############################################################################
options(stringsAsFactors=FALSE)

library(glue)        # for string interpolation
library(configr)     # for (YAML) config file writing


##
# get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
   if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  if(!is.null(default)) default else do.call(convert,list(NA))
}


args <- commandArgs(T)
ESSENTIAL_PROJECT   <- parseArgs(args,"ESSENTIAL_PROJECT=","")                           # main project folder
ESSENTIAL_DB        <- parseArgs(args,"ESSENTIAL_DB=","")                                # UCSC assembly (e.g. mm9, hg19)
ESSENTIAL_CHROMSIZES   <- parseArgs(args,"ESSENTIAL_CHROMSIZES=","")                     # chromosome sizes file
TRACKHUB_TRACKSDIR  <- parseArgs(args,"TRACKHUB_TRACKSDIR=",paste0(ESSENTIAL_PROJECT,"/tracks"))  # project tracks folder
TRACKHUB_CONFIG     <- parseArgs(args,"TRACKHUB_CONFIG=",paste0(ESSENTIAL_PROJECT,"/trackhub.yaml"))  # trackhub configuration file
TRACKHUB_TARGETS    <- parseArgs(args,"TRACKHUB_TARGETS=",paste0(ESSENTIAL_PROJECT,"/targets.txt"))  # file with sample info
TRACKHUB_FTPBASE    <- parseArgs(args,"TRACKHUB_FTPBASE=","")                            # Public FTP root folder
TRACKHUB_FTPURLBASE <- parseArgs(args,"TRACKHUB_FTPURLBASE=","")                         # Public FTP root URL
TRACKHUB_UCSCCFG    <- parseArgs(args,"TRACKHUB_UCSCCFG=","")                            # UCSC configuration files folder
ESSENTIAL_STRANDED  <- parseArgs(args,"ESSENTIAL_STRANDED=","no")                        # set to "yes" for stranded RNA-seq
TRACKHUB_PEAKSDIR   <- parseArgs(args,"TRACKHUB_PEAKSDIR=","")                           # ChIP-seq peaks dir


runstr <- paste0("Call with: Rscript Configure_Trackhub.R ESSENTIAL_PROJECT=",ESSENTIAL_PROJECT,
                 " ESSENTIAL_DB=",ESSENTIAL_DB,
                 " ESSENTIAL_CHROMSIZES=",ESSENTIAL_CHROMSIZES,
                 " [TRACKHUB_TRACKSDIR=",TRACKHUB_TRACKSDIR,
                 "] [TRACKHUB_CONFIG=",TRACKHUB_CONFIG,
                 "] [TRACKHUB_TARGETS=",TRACKHUB_TARGETS,
                 "] [TRACKHUB_FTPBASE=",TRACKHUB_FTPBASE,
                 "] [TRACKHUB_FTPURLBASE=",TRACKHUB_FTPURLBASE,
                 "] [TRACKHUB_UCSCCFG=",TRACKHUB_UCSCCFG,
                 "] [ESSENTIAL_STRANDED=",ESSENTIAL_STRANDED,
                 "] [TRACKHUB_PEAKSDIR=",TRACKHUB_PEAKSDIR,
                 "]")

if(ESSENTIAL_PROJECT == "")     stop(paste("ESSENTIAL_PROJECT is required. Run with:\n",runstr,"\n"))
if(ESSENTIAL_DB == "")          stop(paste("ESSENTIAL_DB (genome assembly) is required. Run with:\n",runstr,"\n"))

# Echo run command for logging purposes
cat(runstr)


# Hierarchical list of configuration variables
# Will be saved as configuration file at the end
cfg <- list()

# Global variables
cfg$Global <- list()
cfg$Global$TRACKHUB_FTPURLBASE <- TRACKHUB_FTPURLBASE
cfg$Global$TRACKHUB_FTPBASE <- TRACKHUB_FTPBASE
cfg$Global$TRACKHUB_UCSCCFG <- TRACKHUB_UCSCCFG
cfg$Global$TRACKHUB_ASSEMBLY <- ESSENTIAL_DB
cfg$Global$TRACKHUB_TRACKSDIR <- TRACKHUB_TRACKSDIR
if ( TRACKHUB_PEAKSDIR != "" ) {
  if ( ESSENTIAL_CHROMSIZES != "" ) {
    cfg$Global$TRACKHUB_PEAKSDIR <- TRACKHUB_PEAKSDIR
    cfg$Global$TRACKHUB_CHROMSIZES <- ESSENTIAL_CHROMSIZES
  } else {
    stop("ESSENTIAL_CHROMSIZES is also required if TRACKHUB_PEAKSDIR is given.\n")
  }
}

# Set track FTP location-related variables (for FTP folder creation & FTP URL construction)
cfg$Global$TRACKHUB_PROJECT <- basename(ESSENTIAL_PROJECT)
cfg$Global$TRACKHUB_GROUP <- basename(dirname(ESSENTIAL_PROJECT))


# Tracks
cfg$Tracks <- list()


## Create the trackhub files


# Read in 'targets.txt' file and determine if this is a ChIP-seq or DNA-/RNA-seq project based on its format
# ChIP-seq files should have the INPUT- & IP-related fields
# DNA-/RNA-seq files should have at least fields named "sample" (sample name) & "file" (bigwig file name)
if (file.exists(TRACKHUB_TARGETS)) {
  targets <- read.delim(TRACKHUB_TARGETS)
  if ( length(setdiff(c("IP","IPname","INPUT","INPUTname"), colnames(targets))) > 0 ) {
    TRACKHUB_EXPTYPE <- "xnaseq"        # targets file should have "sample" & "file" fields
  } else {
    TRACKHUB_EXPTYPE <- "chipseq"       # targets file should have all IP & INPUT fields above
  }
} else {
  TRACKHUB_EXPTYPE <- "unknown"
}

# Generate per-track configuration data for the ChIP-seq/RNA-seq tracks

if(TRACKHUB_EXPTYPE == "xnaseq") {
  if (ESSENTIAL_STRANDED != "no") {   # stranded DNA-/RNA-seq (group strands in multiWig)
    apply(targets, 1, function(x) {
      # multiWig container track entry
      containerName <- unname(x["sample"])
      cfg$Tracks[[containerName]] <<- list()
      cfg$Tracks[[containerName]]$container <<- "multiWig"
      cfg$Tracks[[containerName]]$shortLabel <<- substr(containerName,1,17)
      cfg$Tracks[[containerName]]$longLabel <<- substr(containerName,1,80)
      
      # Strand-specific bigWig tracks
      for (strand in c("fwd", "rev")) {
        trackName <- paste(unname(x["sample"]), strand, sep = ".")
        fileName <- unname(x["file"])
        if (grepl(".readcounts.tsv", fileName)) {        # if RNA-seq, base bigwig name on readcounts name
          fileName <- sub(".readcounts.tsv", glue(".{strand}.bw"), fileName)
        } else {                                         # else "file" should hold bigwig name excl extension
          fileName <- sub("$", glue(".{strand}.bw"), fileName)
        }
        cfg$Tracks[[trackName]] <<- list()
        cfg$Tracks[[trackName]]$type <<- "bigWig"
        cfg$Tracks[[trackName]]$file <<- fileName
        cfg$Tracks[[trackName]]$parent <<- containerName
        cfg$Tracks[[trackName]]$shortLabel <<- trackName
        cfg$Tracks[[trackName]]$longLabel <<- trackName
        if (length(trackName) > 17) {
          cfg$Tracks[[trackName]]$shortLabel <<- substr(trackName,1,13)
          cfg$Tracks[[trackName]]$shortLabel <<- paste(cfg$Tracks[[trackName]]$shortLabel, strand, sep = ".")
        }
        if (length(trackName) > 80) {
          cfg$Tracks[[trackName]]$longLabel <<- substr(trackName,1,76)
          cfg$Tracks[[trackName]]$longLabel <<- paste(cfg$Tracks[[trackName]]$longLabel, strand, sep = ".")
        }
        if (strand == "fwd") {
          cfg$Tracks[[trackName]]$color <<- "32,128,32"
        } else {
          cfg$Tracks[[trackName]]$color <<- "128,32,32"
        }
      }
    })
  } else {                            # unstranded DNA-/RNA-seq (separate tracks)
    apply(targets, 1, function(x) {
      trackName <- unname(x["sample"])
      fileName <- unname(x["file"])
      if (grepl(".readcounts.tsv", x["file"])) {            # if RNA-seq, base bigwig name on readcounts name
        fileName <- sub(".readcounts.tsv", ".bw", fileName)
      } else {                                              # else "file" should hold bigwig name excl extension
        fileName <- sub("$", ".bw", fileName)
      }
      cfg$Tracks[[trackName]] <<- list()
      cfg$Tracks[[trackName]]$type <<- "bigWig"
      cfg$Tracks[[trackName]]$file <<- fileName
      cfg$Tracks[[trackName]]$shortLabel <<- substr(trackName,1,17)
      cfg$Tracks[[trackName]]$longLabel <<- substr(trackName,1,80)
      cfg$Tracks[[trackName]]$color <<- "32,32,128"
    })
  }
} else {
  if(TRACKHUB_EXPTYPE == "chipseq") {   # ChIP-seq (separate tracks, no input duplication)
    apply(targets, 1, function(x) {
      # Extract relevant fields from line, convert to unnamed strings
      INPUTname <- unname(x["INPUTname"])
      IPname <- unname(x["IPname"])
      INPUT <- unname(x["INPUT"])
      IP <- unname(x["IP"])
      
      # Create input signal track if it does not already exist
      if (! glue("{INPUTname}_signal") %in% names(cfg$Tracks)) {
        INPUTbw <- glue("{INPUTname}_signal")
        cfg$Tracks[[INPUTbw]] <<- list()
        cfg$Tracks[[INPUTbw]]$type <<- "bigWig"
        cfg$Tracks[[INPUTbw]]$file <<- sub(".bam",".bw",INPUT)
        cfg$Tracks[[INPUTbw]]$shortLabel <<- substr(INPUTbw,1,17)
        cfg$Tracks[[INPUTbw]]$longLabel <<- substr(INPUTbw,1,80)
        cfg$Tracks[[INPUTbw]]$color <<- "128,128,128"
      }
      # Create IP signal track if not already created
      if (! glue("{IPname}_signal") %in% names(cfg$Tracks)) {
        IPbw <- glue("{IPname}_signal")
        cfg$Tracks[[IPbw]] <<- list()
        cfg$Tracks[[IPbw]]$type <<- "bigWig"
        cfg$Tracks[[IPbw]]$file <<- sub(".bam",".bw",IP)
        cfg$Tracks[[IPbw]]$shortLabel <<- substr(IPbw,1,17)
        cfg$Tracks[[IPbw]]$longLabel <<- substr(IPbw,1,80)
        cfg$Tracks[[IPbw]]$color <<- "128,128,255"
      }
      # Create IP-versus-input peaks track if not already created
      if (! glue("{IPname}.vs.{INPUTname}_peaks") %in% names(cfg$Tracks)) {
        IPbb <- glue("{IPname}.vs.{INPUTname}_peaks")
        cfg$Tracks[[IPbb]] <<- list()
        cfg$Tracks[[IPbb]]$type <<- "bigBed"
        cfg$Tracks[[IPbb]]$file <<- glue("{IPname}.vs.{INPUTname}_macs2_peaks.bb")
        cfg$Tracks[[IPbb]]$shortLabel <<- substr(IPbb,1,17)
        cfg$Tracks[[IPbb]]$longLabel <<- substr(IPbb,1,80)
        cfg$Tracks[[IPbb]]$color <<- "128,128,255"
      }
    })
  } else {
    cat("No targets file found, all tracks will be ungrouped uniform tracks.\n")
  }
}

# Add all remaining bigWig/bigBed tracks to trackhub as uniform grey ungrouped tracks

trackhub_files <- unlist(lapply(cfg$Tracks, "[[", "file"))
trackdir_files <- list.files(path = TRACKHUB_TRACKSDIR, 
                              pattern = "\\.bw$|\\.bb$", 
                              full.names = FALSE)

lapply(setdiff(trackdir_files, trackhub_files), function(track_file) {
  
  track_name <- gsub(".", "_", track_file, fixed = TRUE)
  
  if (length(grep("\\.bw$", track_file)) > 0) track_type <- "bigWig"
  else if (length(grep("\\.bb$", track_file)) > 0) track_type <- "bigBed"
  else track_type <- "unknown"    # should never occur; fix "list.files()" above if so
  
  cfg$Tracks[[track_name]] <<- list()
  cfg$Tracks[[track_name]]$type <<- track_type
  cfg$Tracks[[track_name]]$file <<- track_file
  cfg$Tracks[[track_name]]$shortLabel <<- substr(track_name,1,17)
  cfg$Tracks[[track_name]]$longLabel <<- substr(track_name,1,80)
  cfg$Tracks[[track_name]]$color <<- "128,128,128"
})


# Write out configuration file
write.config(cfg, file.path = TRACKHUB_CONFIG, write.type = "yaml")





