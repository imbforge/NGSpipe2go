// breaktag ESSENTIAL VARIABLES

// Define essential variables here.
// Further module-specific variables can be adjusted in the corresponding ".header" files for each module.
//

// General parameters
ESSENTIAL_PROJECT="your/project/folder"
ESSENTIAL_SAMPLE_PREFIX="" 
ESSENTIAL_THREADS=16

// Mapping parameters
ESSENTIAL_BWA_REF="your/genome/references/human/GRCh38.canonical_assembly.genome.fa"
ESSENTIAL_PAIRED="yes"        // paired end design
ESSENTIAL_QUALITY=60          // min mapping quality of reads to be kept. Defaults to 60

// further optional pipeline stages to include
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes")

// project folders
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
RAWDATA=PROJECT + "/rawdata"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
TARGETS=PROJECT + "/targets.txt"

