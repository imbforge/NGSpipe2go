ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/tmp/ngspipe2go_dnaseq_test"
ESSENTIAL_BWA_REF="/fsimb/common/genomes/saccharomyces_cerevisiae/ensembl/R64/canonical/index/bwa/0.7.15/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa"
ESSENTIAL_CALL_REGION=null
ESSENTIAL_PAIRED="yes"
ESSENTIAL_SAMPLE_PREFIX=""
ESSENTIAL_KNOWN_VARIANTS="/fsimb/groups/imb-bioinfocf/common-data/GATK_resources/yeast/ensembl_98/saccharomyces_cerevisiae.vcf"
ESSENTIAL_THREADS=4

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"

// optional pipeline stages to include
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes")

