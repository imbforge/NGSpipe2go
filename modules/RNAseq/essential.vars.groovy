ESSENTIAL_PROJECT="/project/"
ESSENTIAL_STAR_REF="/annotation/mm9/star_genome"
ESSENTIAL_GENESGTF="/annotation/mm9/gencode.vM1.annotation.gtf"
ESSENTIAL_GENESBED="/annotation/mm9/mm9_UCSC_knownGene.bed"
ESSENTIAL_FEATURETYPE="gene_type" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_SAMPLE_PREFIX=""
ESSENTIAL_PAIRED="no"           // paired end design
ESSENTIAL_STRANDED="reverse"    // strandness: no|yes|reverse
ESSENTIAL_ORG="mouse"           // UCSC organism
ESSENTIAL_DB="mm9"              // UCSC assembly version
ESSENTIAL_READLENGTH=50         // added for STAR version > 2.4.1a
ESSENTIAL_THREADS=4             // number of threads for parallel tasks

//global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
