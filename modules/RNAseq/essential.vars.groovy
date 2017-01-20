ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge/rnaseq_go_test"
ESSENTIAL_STAR_REF="/fsimb/groups/imb-bioinfocf/common-data/star_genomes/Mus_musculus/UCSC/mm9/star2.4.2a_noGTF"
ESSENTIAL_GENESGTF="/fsimb/common/genomes/mus_musculus/ucsc/mm9/canonical/annotation/gencode.vM1.annotation.gtf"
ESSENTIAL_GENESBED="/fsimb/groups/imb-bioinfocf/common-data/annotation/mm9/mm9_UCSC_knownGene.bed"
ESSENTIAL_FEATURETYPE="gene_type" //gencode uses gene_type; ensemble uses gene_biotype
ESSENTIAL_SAMPLE_PREFIX="Sample_imb_richly_2014_05_"
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
