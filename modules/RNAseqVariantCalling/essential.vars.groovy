//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!
//
ESSENTIAL_PROJECT="/local/scratch1/imb-kettinggr/adomingues/projects/snps-splicing"
ESSENTIAL_STAR_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/StarIndex2_4_1d_modified/"
ESSENTIAL_GENOME_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/chr_sequences/chr.clean.fa"
ESSENTIAL_READLENGTH=101

// ESSENTIAL_GENESGTF="/fsimb/groups/imb-bioinfocf/common-data/igenomes_reference/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"        // "/fsimb/groups/imb-bioinfocf/common-data/igenomes_reference/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf" // used for gene mapping (STAR), counting (HTSEQ), duprate analysis
// ESSENTIAL_GENESGTF2="/fsimb/groups/imb-bioinfocf/common-data/annotation/mm9/gencode.vM1.annotation.gtf.gz"   // "/fsimb/groups/imb-bioinfocf/common-data/annotation/mm9/gencode.vM1.annotation.gtf.gz"  // to use with RNAtypes
// ESSENTIAL_GENESBED="/fsimb/groups/imb-bioinfocf/common-data/annotation/mm9/mm9_UCSC_knownGene.bed"      // "/fsimb/groups/imb-bioinfocf/common-data/annotation/mm9/mm9_UCSC_knownGene.bed"
// ESSENTIAL_SAMPLE_PREFIX="Sample_imb_richly_2014_05_"     // "Sample_imb_richly_2014_05_"
// ESSENTIAL_PAIRED="yes"                   // paired end design
// ESSENTIAL_STRANDED="reverse"    // strandness: no|yes|reverse
// ESSENTIAL_ORG="mouse"           //"mouse"               //UCSC organism
// ESSENTIAL_DB="mm9"              //"mm9"                 //UCSC assembly version
// ESSENTIAL_READLENGTH="50" // added for STAR version > 2.4.1a

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/results/mapped"
QC=PROJECT + "data/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=MAPPED + "/tracks"
