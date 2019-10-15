ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/tmp/ngspipe2go_scRNAseq_smartseq2"
ESSENTIAL_STAR_REF="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/tmp/ngspipe2go_scRNAseq_smartseq2/ref/"
ESSENTIAL_GENESGTF="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/tmp/ngspipe2go_scRNAseq_smartseq2/ref/genes_chr19.gtf"
ESSENTIAL_GENESBED="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/tmp/ngspipe2go_scRNAseq_smartseq2/ref/genes_chr19.bed"
ESSENTIAL_CHROMSIZES=""  // chromosome sizes file
ESSENTIAL_FEATURETYPE="gene_biotype" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_SAMPLE_PREFIX=""
ESSENTIAL_PAIRED="no"           // paired end design
ESSENTIAL_STRANDED="yes"    // strandness: no|yes|reverse
ESSENTIAL_ORG="human"           // UCSC organism
ESSENTIAL_DB="hg19"              // UCSC assembly version
ESSENTIAL_READLENGTH=50         // added for STAR version > 2.4.1a
ESSENTIAL_THREADS=4             // number of threads for parallel tasks
ESSENTIAL_FRAGMENT_USAGE="no"   //should fragments be reconstituted? should always be no for rnaseq
ESSENTIAL_FILTER_CHR=""         //chromosomes to include in post-mapping analysis.
ESSENTIAL_BAMCOVERAGE="--binSize 1 --skipNonCoveredRegions --normalizeUsing CPM" // NO smoothing should be done for RNAseq
//Barcodes
ESSENTIAL_WHITELIST="single_cell_barcodes_whitelist.csv" //the barcode list should be a list of valid barcodes separated by newline
ESSENTIAL_BCPATTERN="CCCCCCCNNNNNNNN" //barcode pattern as it is present in MARS-Seq data. 
                                     //C stands for Cellbarcode and N is part of the random barcode
                                     //Please be aware that the addumibarcode module which uses this variable will only 
                                     //work for MARS-Seq where the barcode and the umi are present in the second reads otherwise it has to be modified! 
//Adapter trimming
ESSENTIAL_ADAPTER_SEQUENCE="TruSeqLTHT=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" //standard sequence to trim illumina reads
ESSENTIAL_MINADAPTEROVERLAP=5
ESSENTIAL_MINREADLENGTH=30

//global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
TRIMMED=PROJECT + "/trimmed"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
FUSION=PROJECT + "/fusion"

// optional pipeline stages to include
RUN_TRACKHUB=false
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes")
