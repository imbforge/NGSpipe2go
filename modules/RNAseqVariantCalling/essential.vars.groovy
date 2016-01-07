//Pipeline GATK RNA-seq variant calling
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!
//
ESSENTIAL_PROJECT="/local/scratch1/imb-kettinggr/adomingues/projects/snps-splicing"
ESSENTIAL_STAR_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/StarIndex2_4_1d_modified/"
ESSENTIAL_GENOME_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/chr_sequences/chr.clean.fa"
ESSENTIAL_READLENGTH=101


//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/results/mapped"
QC=PROJECT + "data/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=MAPPED + "/tracks"
