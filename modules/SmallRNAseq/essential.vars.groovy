//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!
//
ESSENTIAL_PROJECT="/local/scratch1/imb-kettinggr/adomingues/projects/bpipe_small_rna"          // "/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imb_cf_2013_04_sayols_infrastructure_pipelines/test"
ESSENTIAL_BOWTIE_PATH="~/bin/bowtie-0.12.8/bowtie"
ESSENTIAL_BOWTIE_REF="/home/adomingu/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv9/Sequence/BowtieIndexChr/chr"

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
QC=PROJECT + "data/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
MAPPED=RESULTS + "/mapped"
TMP=PROJECT + "/tmp"
TRACKS=MAPPED + "/tracks"

// binaries
TOOL_SAMTOOLS="/home/adomingu/bin/samtools"
