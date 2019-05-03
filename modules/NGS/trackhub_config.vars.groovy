//vars for task trackhub, version 1
TRACKHUB_FTPURLBASE="https://hpc1.imb.uni-mainz.de/public"  // public FTP base URL
TRACKHUB_FTPBASE="/fsimb/services/ftp/public/"  // public FTP root folder
TRACKHUB_UCSCCFG="/fsimb/common/tools/ucsc/config/"  // folder with UCSC tools configuration files (e.g. narrowPeak.as/broadPeak.as)
TRACKHUB_TARGETS=ESSENTIAL_PROJECT + "/targets.txt"  // targets file describing the samples for ChIP-seq
TRACKHUB_PEAKSDIR=RESULTS + "/macs2"  // location of peak files for ChIP-seq (comment out if no peak files)
TRACKHUB_TRACKSDIR=TRACKS  // location of track files for putting into trackhub
TRACKHUB_CONFIG=ESSENTIAL_PROJECT + "/trackhub.yaml"  // trackhub configuration file

// Variables located in "essential.vars.groovy" that are relevant to this module:
// ESSENTIAL_DB  (UCSC genome assembly, e.g. "hg19")
// ESSENTIAL_CHROMSIZES  (chromosome sizes file)
// ESSENTIAL_STRANDED  (stranded sequencing or not, for strans-specific bigwig creation)
// TRACKS  (full path of project subdirectory containing tracks)
// RESULTS  (full path of project subdirectory containing ChIP-seq peak files)

