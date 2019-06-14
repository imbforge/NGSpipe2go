ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/scRNA_MARSseq_test/"
ESSENTIAL_STAR_REF="/fsimb/common/genomes/mus_musculus/ensembl/grcm38/canonical/index/star_ERCC/2.5.2b/"
ESSENTIAL_GENESGTF="/fsimb/common/genomes/mus_musculus/ensembl/grcm38/canonical/annotation/Mus_musculus.GRCm38.88+ERCC92.gtf"
ESSENTIAL_GENESBED="/fsimb/common/genomes/mus_musculus/ensembl/grcm38/canonical/annotation/Mus_musculus.GRCm38.89.bed"
ESSENTIAL_CHROMSIZES=""  // chromosome sizes file
ESSENTIAL_FEATURETYPE="gene_biotype" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_SAMPLE_PREFIX="test_"
ESSENTIAL_PAIRED="no"           // paired end design
ESSENTIAL_STRANDED="yes"    // strandness: no|yes|reverse
ESSENTIAL_ORG="mouse"           // UCSC organism
ESSENTIAL_DB="mm10"              // UCSC assembly version
ESSENTIAL_READLENGTH=70         // added for STAR version > 2.4.1a
ESSENTIAL_THREADS=4             // number of threads for parallel tasks
ESSENTIAL_FRAGMENT_USAGE="no"   //should fragments be reconstituted? should always be no for rnaseq
ESSENTIAL_FILTER_CHR=""         //chromosomes to include in post-mapping analysis.
//Maybe add
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
