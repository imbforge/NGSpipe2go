ESSENTIAL_PROJECT="/project"
ESSENTIAL_BWA_REF="/data/igenomes_reference/Mus_musculus/UCSC/mm9/Sequence/BWAIndex/genome.fa"
ESSENTIAL_CALL_REGION="/project/dnaseq/call.region.bed" // OPTIONALLY limit variant calling to a limited region, e.g. exonic sequences or chr3 only. Will save a lot of compute time.
ESSENTIAL_PAIRED="yes"                   // paired end design
ESSENTIAL_SAMPLE_PREFIX="Sample_" 
ESSENTIAL_KNOWN_VARIANTS="/data/GATK_resources/2.8/hg19/dbsnp_138.hg19.vcf" // KnownVariants could dbSNP from GATK resource bundle. This is crucial for BaseQualityRecalibration step! vcf file needs to be unzipped
ESSENTIAL_HAPMAP_VARIANTS="/data/GATK_resources/2.8/hg19/hapmap_3.3.hg19.sites.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
ESSENTIAL_OMNI_VARIANTS="/data/GATK_resources/2.8/hg19/1000G_omni2.5.hg19.sites.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
ESSENTIAL_MILLS_VARIANTS="/data/GATK_resources/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
ESSENTIAL_THOUSAND_GENOMES_VARIANTS="/data/GATK_resources/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
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

