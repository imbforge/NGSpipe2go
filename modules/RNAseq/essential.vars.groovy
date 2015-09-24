//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!
//
ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/rnaseq_pipe_modular"          // "/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imb_cf_2013_04_sayols_infrastructure_pipelines/test"
ESSENTIAL_STAR_REF="/fsimb/groups/imb-bioinfocf/common-data/star_genomes/Mus_musculus/UCSC/mm9/star2.4.2a_noGTF/"        // "/fsimb/groups/imb-bioinfocf/common-data/star_genomes/mm9/"
ESSENTIAL_GENESGTF="/fsimb/groups/imb-bioinfocf/common-data/igenomes_reference/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"        // "/fsimb/groups/imb-bioinfocf/common-data/igenomes_reference/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf" // used for gene mapping (STAR), counting (HTSEQ), duprate analysis
ESSENTIAL_GENESGTF2="/fsimb/groups/imb-bioinfocf/common-data/annotation/mm9/gencode.vM1.annotation.gtf.gz"   // "/fsimb/groups/imb-bioinfocf/common-data/annotation/mm9/gencode.vM1.annotation.gtf.gz"  // to use with RNAtypes
ESSENTIAL_GENESBED="/fsimb/groups/imb-bioinfocf/common-data/annotation/mm9/mm9_UCSC_knownGene.bed"      // "/fsimb/groups/imb-bioinfocf/common-data/annotation/mm9/mm9_UCSC_knownGene.bed"
ESSENTIAL_SAMPLE_PREFIX="Sample_imb_richly_2014_05_"     // "Sample_imb_richly_2014_05_"
ESSENTIAL_PAIRED="no"                   // paired end design
ESSENTIAL_STRANDED="reverse"    // strandness: no|yes|reverse
ESSENTIAL_ORG="mouse"           //"mouse"               //UCSC organism
ESSENTIAL_DB="mm9"              //"mm9"                 //UCSC assembly version
ESSENTIAL_READLENGTH="50" // added for STAR version > 2.4.1a

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
