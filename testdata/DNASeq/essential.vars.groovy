ESSENTIAL_PROJECT="/fsimb/groups/imb-bioinfocf/projects/cfb_internal/tmp/ngspipe2go_dnaseq_test"
ESSENTIAL_BWA_REF="/fsimb/groups/imb-bioinfocf/common-data/GATK_resources/gatk_bundle_hg38_v0/Homo_sapiens_assembly38.fasta"
ESSENTIAL_CALL_REGION=null
ESSENTIAL_PAIRED="yes"
ESSENTIAL_SAMPLE_PREFIX=""
ESSENTIAL_KNOWN_VARIANTS="/fsimb/groups/imb-bioinfocf/common-data/GATK_resources/gatk_bundle_hg38_v0/dbsnp_146.hg38.vcf.gz" // KnownVariants could dbSNP from GATK resource bundle. This is crucial for BaseQualityRecalibration step! vcf file needs to be unzipped
ESSENTIAL_HAPMAP_VARIANTS="/fsimb/groups/imb-bioinfocf/common-data/GATK_resources/gatk_bundle_hg38_v0/hapmap_3.3.hg38.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
ESSENTIAL_OMNI_VARIANTS="/fsimb/groups/imb-bioinfocf/common-data/GATK_resources/gatk_bundle_hg38_v0/1000G_omni2.5.hg38.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
ESSENTIAL_MILLS_VARIANTS="/fsimb/groups/imb-bioinfocf/common-data/GATK_resources/gatk_bundle_hg38_v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
ESSENTIAL_THOUSAND_GENOMES_VARIANTS="/fsimb/groups/imb-bioinfocf/common-data/GATK_resources/gatk_bundle_hg38_v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz" // varirants provided by the GATK bundle. Essential for Variant Score Recalibration
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

