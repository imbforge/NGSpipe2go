//Pipeline GATK RNA-seq variant calling
ESSENTIAL_PROJECT="/tmp/ngspipe2go_rnaseqvariantcalling_test"
ESSENTIAL_STAR_REF="/tmp/ngspipe2go_rnaseqvariantcalling_test/ref/"
ESSENTIAL_GENOME_REF="/tmp/ngspipe2go_rnaseqvariantcalling_test/ref/ref.fa"
ESSENTIAL_VCF_REF="/tmp/ngspipe2go_rnaseqvariantcalling_test/knowVariants.vcf"
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
