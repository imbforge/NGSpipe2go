//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!
//
ESSENTIAL_PROJECT="/local/scratch1/imb-kettinggr/adomingues/projects/bpipe_small_rna"
ESSENTIAL_BOWTIE_PATH="~/bin/bowtie-0.12.8/bowtie"
ESSENTIAL_BOWTIE_REF="/home/adomingu/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv9/Sequence/BowtieIndexChr/chr"
ESSENTIAL_GENOME_REF="/home/adomingu/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv9/Sequence/chr_sequences/chr.clean.fa"
ESSENTIAL_FEATURES="/local/scratch1/imb-kettinggr/adomingues/projects/bpipe_small_rna/data/annotations/features.bed"

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
QC=PROJECT + "/data/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
MAPPED=RESULTS + "/mapped"
TMP=PROJECT + "/tmp"
TRACKS=MAPPED + "/tracks"

// binaries
TOOL_SAMTOOLS="/home/adomingu/bin/samtools"
