//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!
//
ESSENTIAL_PROJECT="/project"
ESSENTIAL_ORG="mouse"
ESSENTIAL_DB="mm9"
ESSENTIAL_BWA_REF="/data/igenomes_reference/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/genome.fa"
ESSENTIAL_GATK_REF=""
ESSENTIAL_CALL_REGION="/project/dnaseq/call.region.bed" // OPTIONALLY limit variant calling to a limited region, e.g. exonic sequences or chr3 only. Will save a lot of compute time.
ESSENTIAL_PAIRED="yes"                   // paired end design
ESSENTIAL_SAMPLE_PREFIX="NA12877_"     // "Sample_imb_richly_2014_05_"
ESSENTIAL_KNOWN_VARIANTS="/data/GATK_resources/2.8/hg19/dbsnp_138.hg19.vcf" // KnownVariants could dbSNP from GATK resource bundle. This is crucial for BaseQualityRecalibration step! vcf file needs to be unzipped
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
