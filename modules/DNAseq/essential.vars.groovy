//Pipeline generated with command line: ./imb-pip.pl --tasks-pip=1 --force
//By: ssayolsp At: Fr 17 Okt 2014 17:12:41 CEST
//
// REMEMBER TO CHANGE THESE ESSENTIAL VARS!!
//
ESSENTIAL_PROJECT="project-folder"
ESSENTIAL_BWA_REF=""
ESSENTIAL_CALL_REGION="" // OPTIONAL! limit variant calling to a limited region, e.g. exonic sequences or chr3 only. Will save a lot of compute time.
ESSENTIAL_SAMPLE_PREFIX=""  //
ESSENTIAL_KNOWN_VARIANTS="" // KnownVariants could dbSNP from GATK resource bundle. This is crucial for BaseQualityRecalibration step! vcf file needs to be unzipped

//global vars that will be reused in some global vars
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
