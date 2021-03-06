//Pipeline GATK RNA-seq variant calling
ESSENTIAL_PROJECT="/local/scratch1/imb-kettinggr/adomingues/projects/snps-splicing"
ESSENTIAL_STAR_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/StarIndex2_4_1d_modified/"
ESSENTIAL_GENOME_REF="/fsimb/groups/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Sequence/chr_sequences/chr.clean.fa"
ESSENTIAL_VCF_REF="/home/adomingu/imb-kettinggr/genomes/Danio_rerio/Ensembl/Zv10/Annotation/variation/Danio_rerio.vcf.gz"
ESSENTIAL_READLENGTH=101
ESSENTIAL_THREADS=4

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "data/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=MAPPED + "/tracks"

// optional pipeline stages to include
