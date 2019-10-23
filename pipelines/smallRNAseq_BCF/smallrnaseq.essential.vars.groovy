//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!
//
ESSENTIAL_PROJECT="/project/"
ESSENTIAL_BOWTIE_REF="/data/mus_musculus/gencode/release-M10_GRCm38.p4/full/index/bowtie/1.1.2/GRCm38.p4"
ESSENTIAL_GENOME_REF="/data/mus_musculus/gencode/release-M10_GRCm38.p4/full/genome/GRCm38.p4.genome_chr.fa" // necessary for miRDeep2

ESSENTIAL_GENESGTF="/data/mus_musculus/gencode/release-M10_GRCm38.p4/full/annotation/gencode.vM10.chr_patch_hapl_scaff.annotation.gtf"
ESSENTIAL_RRNA_BOWTIE_REF="/data/rrna/mouse/mouse_all_rRNA" // necessary for fastqscreen

ESSENTIAL_SPECIES="Mouse"      // necessary for miRDeep2 and fastqscreen
ESSENTIAL_SAMPLE_PREFIX="Sample_"
ESSENTIAL_FEATURETYPE="gene_type" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_PAIRED="no"          // paired end design
ESSENTIAL_STRANDED="yes"       // strandness: no|yes|reverse
ESSENTIAL_THREADS=4            // number of threads for parallel tasks

ESSENTIAL_READLENGTH=51        // actual read length in original raw data (incl. insert, UMIs, adapter)
ESSENTIAL_MINADAPTEROVERLAP=5  // minimal overlap with adapter
ESSENTIAL_MINREADLENGTH=26     // remaining read length plus UMIs (2x4) 
ESSENTIAL_UMI_LENGTH=8         // (2x4bp)
ESSENTIAL_ADAPTER_SEQUENCE="TGGAATTCTCGGGTGCCAAGG" // needed for cutadapt adapter trimming

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

