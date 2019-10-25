//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!
//
ESSENTIAL_PROJECT="./"
ESSENTIAL_BOWTIE_REF="./ref/mmu"
ESSENTIAL_GENOME_REF="./ref/mmu.fa"

ESSENTIAL_GENESGTF="./ref/mmu.gtf"
ESSENTIAL_RRNA_BOWTIE_REF="./ref/rrna"

ESSENTIAL_SPECIES="Mouse"       // necessary for miRDeep2, used to refer to UCSC
ESSENTIAL_SAMPLE_PREFIX=""
ESSENTIAL_FEATURETYPE="gene_biotype" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_PAIRED="no"           // paired end design
ESSENTIAL_STRANDED="yes"        // strandness: no|yes|reverse
ESSENTIAL_THREADS=4             // number of threads for parallel tasks

ESSENTIAL_READLENGTH=51         // actual read length in original raw data (incl. insert, UMIs, adapter)
ESSENTIAL_MINADAPTEROVERLAP=5   // minimal overlap with adapter
ESSENTIAL_MINREADLENGTH=15      // remaining read length plus UMIs (2x4) 
ESSENTIAL_UMI_LENGTH=4          // (2x4bp)
ESSENTIAL_ADAPTER_SEQUENCE="TGGAATTCTCGGGTGCCAAGG" // needed for cutadapt adapter trimming

// vars for mirDeep2
ESSENTIAL_MATURE_MIRNA="./ref/mature.fa"
ESSENTIAL_HAIRPIN_MIRNA="./ref/hairpin.fa"

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
PROCESSED=PROJECT + "/rawdata_processed"
MAPPED=PROJECT + "/mapped"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"

// optional pipeline stages to include

