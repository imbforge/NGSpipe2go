//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!
//
ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge/test/dnaseq"          // "/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imb_cf_2013_04_sayols_infrastructure_pipelines/test"
ESSENTIAL_BWA_REF="/fsimb/groups/imb-bioinfocf/common-data/igenomes_reference/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa"        // "/fsimb/groups/imb-bioinfocf/common-data/star_genomes/mm9/"
ESSENTIAL_GATK_REF=""
ESSENTIAL_PAIRED="yes"                   // paired end design
ESSENTIAL_SAMPLE_PREFIX="NA12877_"     // "Sample_imb_richly_2014_05_"
ESSENTIAL_KNOWN_VARIANTS="/fsimb/groups/imb-bioinfocf/common-data/GATK_resources/2.8/hg19/dbsnp_138.hg19.vcf" // KnownVariants could dbSNP from GATK resource bundle. This is crucial for BaseQualityRecalibration step! vcf file needs to be unzipped

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
