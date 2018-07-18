############################################################################
##
## What: Make_Trackhub.R
## Who: Martin Oti
## When: 2018-07-17
##
## Script to create UCSC Track Hub with tracks from specific project.
## Reads configuration data from a pre-generated configuration file.
## Configuration file should be in YAML format with sections:
## - Global: All non-track-specific configuration data
## - Tracks: Per-track configuration data (for 'trackDb.txt' file)
##
## Args: 
## -----
## ESSENTIAL_PROJECT=   # main project folder
## TRACKHUB_CONFIG=     # trackhub generator configuration file
##
## ----------------------------------------------------------------
## Trackhub configuration file must contain these global variables:
## ----------------------------------------------------------------
##  TRACKHUB_FTPURLBASE: # Main FTP site URL
##  TRACKHUB_FTPBASE:    # File system location of main FTP site
##  TRACKHUB_UCSCCFG:    # File system location of UCSC track type specification files
##  TRACKHUB_ASSEMBLY:   # UCSC genome assembly (e.g. hg38)
##  TRACKHUB_CHROMSIZES: # Chromosome sizes file for the relevant UCSC genome assembly
##  TRACKHUB_TRACKSDIR:  # Directory containing bigWig/bigBed track files
##  TRACKHUB_PEAKSDIR:   # Directory containing narrowPeak/broadPeak files (optional)
##
############################################################################
options(stringsAsFactors=FALSE)

library(glue)        # for string interpolation
library(configr)     # for (YAML) config file parsing



##
# get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
   if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  if(!is.null(default)) default else do.call(convert,list(NA))
}


args <- commandArgs(T)
ESSENTIAL_PROJECT  <- parseArgs(args,"ESSENTIAL_PROJECT=","")                                        # main project folder
TRACKHUB_CONFIG    <- parseArgs(args,"TRACKHUB_CONFIG=",paste0(ESSENTIAL_PROJECT,"/trackhub.yaml"))  # trackhub configuration file


runstr <- glue("Call with: Rscript Make_Trackhub.R ESSENTIAL_PROJECT={ESSENTIAL_PROJECT} TRACKHUB_CONFIG={TRACKHUB_CONFIG}")

if(ESSENTIAL_PROJECT == "")       stop(paste("ESSENTIAL_PROJECT is required. Run with:\n",runstr))
if(!file.exists(TRACKHUB_CONFIG)) stop(paste("TRACKHUB_CONFIG file",TRACKHUB_CONFIG,"does NOT exist. Run with:\n",runstr))

# Echo run command for logging purposes
cat(runstr)


setwd(ESSENTIAL_PROJECT)


# Read in the trackhub configuration data
cfg <- read.config(TRACKHUB_CONFIG)

# Ensure all required configuration variables are present in config file
required_vars <- c("TRACKHUB_FTPURLBASE", 
                   "TRACKHUB_FTPBASE", 
                   "TRACKHUB_UCSCCFG", 
                   "TRACKHUB_ASSEMBLY", 
                   "TRACKHUB_CHROMSIZES", 
                   "TRACKHUB_TRACKSDIR")
if ( !("Global" %in% names(cfg)) ) {
  stop(glue("Config file does not have a 'Global' section with the global trackhub variables.\n"))
}
missing_vars <- setdiff(required_vars, names(cfg$Global))
if ( length(missing_vars) > 0 ) {
  stop(glue("Missing required trackhub config file variable: {missing_vars}; \n"))
}



# Create relevant derived variables from pre-existing variables
# Assumes ESSENTIAL_PROJECT path ends with GROUP and then PROJECT names
# FTP public folder location and URL are created based on these
PROJECT <- basename(ESSENTIAL_PROJECT)
GROUP <- basename(dirname(ESSENTIAL_PROJECT))
PROJECT_FTPDIR <- file.path(cfg$Global$TRACKHUB_FTPBASE, GROUP, PROJECT)
PROJECT_FTPURL <- paste(cfg$Global$TRACKHUB_FTPURLBASE, GROUP, PROJECT, sep='/')

if (!dir.exists(file.path(PROJECT_FTPDIR, "tracks/ucsc_track_hub", cfg$Global$TRACKHUB_ASSEMBLY))) {
  dir.create(file.path(PROJECT_FTPDIR, "tracks/ucsc_track_hub", cfg$Global$TRACKHUB_ASSEMBLY), recursive=TRUE)
}

# If peak files are available, convert them to bigBed format
# Peaks should be in narrowPeak or broadPeak format (MACS2 generates these)
# MACS2 puts values > 1000 in score field; cap these at 1000 to satisfy bedToBigBed requirements
# Note: 'bedToBigBed' requires narrowPeak/broadPeak field description files for bigBed conversion
#       They can be obtained from:
#         https://raw.githubusercontent.com/ENCODE-DCC/encValData/master/as/narrowPeak.as
#         https://raw.githubusercontent.com/ENCODE-DCC/encValData/master/as/broadPeak.as
if (dir.exists(cfg$Global$TRACKHUB_PEAKSDIR)) {
  # .narrowPeak files if any
  lapply(list.files(cfg$Global$TRACKHUB_PEAKSDIR, pattern = ".narrowPeak$", full.names = TRUE),
         function(PEAKFILE) {
           BBFILE <- sub(".narrowPeak$", ".bb", file.path(cfg$Global$TRACKHUB_TRACKSDIR, basename(PEAKFILE)))
           if ( !file.exists(BBFILE) || (file.mtime(BBFILE) < file.mtime(PEAKFILE)) ) {
             system(sprintf("cat %s | awk -v OFS='\t' '{if($5>1000)$5=1000; print $0}' > %s.fixed", PEAKFILE, PEAKFILE))
             system(sprintf("bedToBigBed -type=bed6+4  -as=%s/narrowPeak.as  %s.fixed  %s  %s",
                            cfg$Global$TRACKHUB_UCSCCFG, PEAKFILE, cfg$Global$TRACKHUB_CHROMSIZES, BBFILE))
             system(sprintf("rm %s.fixed", PEAKFILE))
           }
         })
  # .broadPeak files if any
  lapply(list.files(cfg$Global$TRACKHUB_PEAKSDIR, pattern = ".broadPeak$", full.names = TRUE),
         function(PEAKFILE) {
           BBFILE <- sub(".broadPeak$", ".bb", file.path(cfg$Global$TRACKHUB_TRACKSDIR, basename(PEAKFILE)))
           if ( !file.exists(BBFILE) || (file.mtime(BBFILE) < file.mtime(PEAKFILE)) ) {
             system(sprintf("cat %s | awk -v OFS='\t' '{if($5>1000)$5=1000; print $0}' > %s.fixed", PEAKFILE, PEAKFILE))
             system(sprintf("bedToBigBed -type=bed6+3  -as=%s/broadPeak.as  %s.fixed  %s  %s",
                            cfg$Global$TRACKHUB_UCSCCFG, PEAKFILE, cfg$Global$TRACKHUB_CHROMSIZES, BBFILE))
             system(sprintf("rm %s.fixed", PEAKFILE))
           }
         })
}

# Copy the tracks to the FTP location
# Only copy the files given in the trackhub configuration file
trackfiles <- unlist(lapply(cfg$Tracks, "[[", "file"))
file.copy(file.path(cfg$Global$TRACKHUB_TRACKSDIR, trackfiles), file.path(PROJECT_FTPDIR, "tracks"))




####
## Create the trackhub files
####


# "genomes.txt"
GENOMESTXT_FILE <- file.path(PROJECT_FTPDIR, "tracks/ucsc_track_hub/genomes.txt")
GENOMESTXT_STR <- glue(
  'genome {cfg$Global$TRACKHUB_ASSEMBLY}\n',
  'trackDb {cfg$Global$TRACKHUB_ASSEMBLY}/trackDb.txt')
cat(GENOMESTXT_STR, file = GENOMESTXT_FILE, sep = "")

# "hub.txt" (truncate project name to fit in shortLabel & longLabel)
shortLabel <- unlist(strsplit(PROJECT, "_"))
shortLabel <- paste(shortLabel[2:length(shortLabel)], collapse = "_")
shortLabel <- substr(shortLabel, 1, 17)
longLabel <- substr(PROJECT, 1, 80)
HUBTXT_FILE <- file.path(PROJECT_FTPDIR, "tracks/ucsc_track_hub/hub.txt")
HUBTXT_STR <- glue(
  'hub {PROJECT}\n',
  'shortLabel {shortLabel}\n',
  'longLabel {longLabel}\n',
  'genomesFile genomes.txt\n',
  'email imb-cf-bioinformatics@lists.uni-mainz.de\n')
cat(HUBTXT_STR, file = HUBTXT_FILE, sep = "")

# "trackDb.txt"
TRACKDBTXT_FILE <- file.path(PROJECT_FTPDIR, "tracks/ucsc_track_hub", cfg$Global$TRACKHUB_ASSEMBLY, "trackDb.txt")

# Remove "trackDb.txt" file if it exists, as commands below all append to it rather than overwriting,
# which can lead to duplication of entries that break the trackhub
if (file.exists(TRACKDBTXT_FILE)) {
  file.remove(TRACKDBTXT_FILE)
}



##
# Helper functions for trackDb.txt file writing
##

# multiWig
# Required params:
#   params$shortLabel
#   params$longLabel
write_track_multiWig <- function(trackdbfile, trackname, params) {
  defaultparams <- list(
    'container' = 'multiWig',
    'type' = 'bigWig',
    'viewLimits' = '-10:10',
    'visibility' = 'full',
    'smoothingWindow' = '5',
    'maxHeightPixels' = '150:30:11',
    'aggregate' = 'transparentOverlay',
    'showSubtrackColorOnUi' = 'off',
    'windowingFunction' = 'mean',
    'priority' = '1.4',
    'configurable' = 'on',
    'autoScale' = 'on'
  )
  trackparams <- config.list.merge(defaultparams, params)
  tracklines <- paste(names(trackparams), trackparams, sep = " ")
  cat(glue("track {trackname}"), file = trackdbfile, sep = "\n", append = TRUE)
  cat(tracklines, "", file = trackdbfile, sep = "\n", append = TRUE)
}

# bigWig
# Required params:
#   params$file
#   params$shortLabel
#   params$longLabel
# Frequently required optional params:
#   params$parent
#   params$color
write_track_bigWig <- function(trackdbfile, trackname, params) {
  defaultparams <- list(
    'type' = 'bigWig',
    'color' = '64,64,128',
    'viewLimits' = '-10:10',
    'visibility' = 'full',
    'smoothingWindow' = '5',
    'maxHeightPixels' = '150:30:11',
    'windowingFunction' = 'mean',
    'priority' = '1.4',
    'configurable' = 'on',
    'autoScale' = 'on'
  )
  trackparams <- config.list.merge(defaultparams, params)
  tracklines <- paste(names(trackparams), trackparams, sep = " ")
  cat(glue("track {trackname}"), file = trackdbfile, sep = "\n", append = TRUE)
  cat(tracklines, "", file = trackdbfile, sep = "\n", append = TRUE)
}

# bigBed
# Required params:
#   params$file
#   params$shortLabel
#   params$longLabel
# Frequently required optional params:
#   params$color
write_track_bigBed <- function(trackdbfile, trackname, params) {
  defaultparams <- list(
    'type' = 'bigBed',
    'color' = '64,64,128',
    'visibility' = 'pack',
    'maxHeightPixels' = '150:30:11',
    'configurable' = 'on',
    'autoScale' = 'on'
  )
  trackparams <- config.list.merge(defaultparams, params)
  tracklines <- paste(names(trackparams), trackparams, sep = " ")
  cat(glue("track {trackname}"), file = trackdbfile, sep = "\n", append = TRUE)
  cat(tracklines, "", file = trackdbfile, sep = "\n", append = TRUE)
}


##
# Write tracks out to 'trackDb.txt'
# Strategy: for each track entry in configuration file, call the relevant
# writer function (multiWig/bigWig/bigBed) to write out track entry
##


lapply(seq_along(cfg$Tracks), function(i) {
  
  track_params <- cfg$Tracks[[i]]
  trackname <- names(cfg$Tracks)[[i]]
  
  # Convert track file name to 'bigDataUrl' track field
  if ("file" %in% names(track_params)) {
    track_params$bigDataUrl <- glue("{PROJECT_FTPURL}/tracks/{track_params$file}")
    track_params$file <- NULL
  }
  
  # Call appropriate track writer function based on track 'type' field
  # (or 'container' field if it is a container track)
  if ("container" %in% names(track_params)) {
    track_writer_func <- glue("write_track_{track_params$container}")
  } else {
    track_writer_func <- glue("write_track_{track_params$type}")
  }
  if (exists(track_writer_func)) {
    get(track_writer_func)(TRACKDBTXT_FILE, trackname, track_params)
  } else {
    print(glue("Unrecognized track type: {track_params$type}"))
  }
})


####
## Round off
####

# Make everything in the FTP location world-readable
system(glue("find {PROJECT_FTPDIR} -type d -exec chmod a+rx {{}} \\; &&
             find {PROJECT_FTPDIR} -type f -exec chmod a+r {{}} \\;"))

# Finally mark as finished
cat("Trackhub successfully created. URL for UCSC GB:",
    file = file.path(ESSENTIAL_PROJECT,"trackhub.done"), sep = "\n")
# Also add trackhub URL for UCSC GB to the 'trackhub.done' file
cat(glue("{PROJECT_FTPURL}/tracks/ucsc_track_hub/hub.txt"), 
    file = file.path(ESSENTIAL_PROJECT,"trackhub.done"), sep = "\n", append = TRUE)


