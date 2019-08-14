ESSENTIAL_PROJECT="/project"
ESSENTIAL_BWA_REF="/data/igenomes_reference/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/genome.fa"
ESSENTIAL_CALL_REGION="/project/dnaseq/call.region.bed" // OPTIONALLY limit variant calling to a limited region, e.g. exonic sequences or chr3 only. Will save a lot of compute time.
ESSENTIAL_PAIRED="yes"                   // paired end design
ESSENTIAL_SAMPLE_PREFIX="Sample_" 
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

// optional pipeline stages to include

